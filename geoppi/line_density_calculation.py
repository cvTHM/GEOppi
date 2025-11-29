# %%

# -*- coding: utf-8 -*-

import pandas as pd
import geopandas as gp
import numpy as np

from geoppi.auxFunctions import (create_circumferential_points, create_point_splitter_npoints, closest_objects_to_points, assign_attr_by_max_intersection_area, )

def closest_lines_to_polygons(
        polygons:gp.GeoDataFrame,
        lines:gp.GeoDataFrame,
        maxDistances:np.array
        )->gp.GeoDataFrame:
    
    """
    Function that matches index of closest line object to each polygon provided in **polygons** based on spatial indexing with pbox.

    :param polygons: GeoDataFrame containing polygon objects.\n
    :param lines: GeoDataFrame containting line objects in which to searhc for closest objects.\n
    :param maxDistances: array-like or float/integer defining the polygon-wise or constant maximum value for the distance to apply matching. If no matching object within the maximum distance is found, np.nan is entered. \n

    :return: polygons as GeoDataFrame, complemented by additional column for closest line object index and distance
    """

    # Plausibility checks
    if type(maxDistances) in (int, float):
        maxDistances = np.ones(polygons.shape[0]) * maxDistances    

    # Create unique matching id (temporary) for polygons
    matchID_temp = 'matchID_temp'
    polygons[matchID_temp] = np.arange(len(polygons))

    # Create circumferential points around polygons
    points = create_circumferential_points(polygons = polygons, npoints = 10, lmin = 5)
    points[matchID_temp] = polygons[matchID_temp]

    ## Match corresponding max distance to points
    points['maxDistance'] = maxDistances
    points = points.explode(index_parts = True).reset_index(drop = True)

    maxDistances = np.array(points['maxDistance'])

    points.drop(columns = ['maxDistance'], inplace = True)

    print('\n### All circumferential points around polygons are generated ###')

    # Plausibility checks
    if len(points.geom_type.unique()) == 1:
        if not points.geom_type.unique()[0] == 'Point':
            print('\n### Please provide point objects only. ###')
            return points
    else:
        print('\n### Please provide consistent geometry objects with only geometry type "Point". ###')
        return points    

    # Find closest object to each point (indicating index of closest line segment)
    idxs, distances = closest_objects_to_points(points = points, geomObjects = lines, maxDist = maxDistances)

    seriesIdx = pd.Series(idxs)
    distances = pd.Series(distances)      
        
    # Assign nearest street segment to each point at buildings' boundaries
    points['nearestline']                            = seriesIdx.values
    points['distance']                               = distances.values

    # Find index of points which feature the shortes distance to adjacent line segment within group of points for each polygon
    idx_ps_min_distances = points[~points['nearestline'].isnull()].groupby(by=matchID_temp)['distance'].idxmin().values
    points_ed = points.loc[idx_ps_min_distances]

    # Assign ID of closest line segment to each polygon by matching informatin from respective points
    polygons = pd.merge(polygons, points_ed[[matchID_temp, 'nearestline', 'distance']], left_on=matchID_temp, right_on=matchID_temp, how = 'left')
    
    polygons.drop(columns = [matchID_temp], inplace = True)

    if polygons['nearestline'].isnull().any():
        import warnings
        warnings.warn('Some objects were not attached to the network. '
                      'Set larger min_size. {} affected elements'.format(sum(polygons['nearestline'].isnull())))

    return polygons

def summed_attributes_on_lines(
        lines:gp.GeoDataFrame,
        polygons:gp.GeoDataFrame,
        spatial_distribution:gp.GeoDataFrame = None,
        spatial_connection_ratio:str = None,
        target_attr:str = None,
        additional_attr:list = None,
        rand_sampling:bool = False,
        nSamples:int = 1,
        target_connection_ratio:float = 1,
        nPolygons_min:int = 1,
        func_max_distance = None
    )->gp.GeoDataFrame:

    ### Initializations
    cs = lines.crs

    polygons.to_crs(cs, inplace = True)

    if spatial_distribution is not None:
        spatial_distribution.to_crs(cs, inplace = True)

    # Check for availability of defined attributes
    if target_attr not in polygons.columns:
        print(f'\n### Attribute {target_attr} not found in polygons attributes. Aborting...')
        return
       
    summed_target_attrs = set([col for col in additional_attr if col in polygons.columns] + [target_attr])

    # Define unique IDs for identification and matching
    polygons_uniqueID = 'unique_ID_polys'
    lines_uniqueID = 'unique_ID_lines'
    spatial_uniqueID = 'unique_ID_spatial'

    polygons[polygons_uniqueID] = np.arange(len(polygons)).astype(int)
    lines[lines_uniqueID] = np.arange(len(lines)).astype(int)
    
    if spatial_distribution is not None:
        spatial_distribution[spatial_uniqueID] = np.arange(len(spatial_distribution)).astype(int)

        if (spatial_connection_ratio is None) | (spatial_connection_ratio not in spatial_distribution.columns):
            print(f'\n Attribute {spatial_connection_ratio} is set to None or not found in DataFrame columns. Please provide proper attribute name in DataFrame spatial_distribution.')

        polygons = assign_attr_by_max_intersection_area(gp1 = polygons, gp_source=spatial_distribution, gp1_id = polygons_uniqueID, attr=spatial_uniqueID)

    # Define maximum searching distance by which to assign closets line object to polygon (function-wise)
    if func_max_distance is not None:
        maxDistancesArray = np.array(list(map(lambda x: func_max_distance(x), np.array(polygons[target_attr]))))

    else:
        maxDistancesArray = 500 # Default value

    ### Start calculation    
    # Assign unique ID of closest line object to polygons
    polygons = closest_lines_to_polygons(polygons = polygons, lines = lines, maxDistances = maxDistancesArray)

    # Transfer unique ID of lines to polygons instead of line index
    dictIdxIDLines = dict(zip(list(lines.index), list(lines[lines_uniqueID])))
    polygons['nearestline'] = polygons['nearestline'].apply(lambda x: dictIdxIDLines[x] if x in dictIdxIDLines.keys() else x)
    
    # Create copies
    lines_out = lines.copy()        

    # Select buildings according to current filter
    selected_polygons = gp.GeoDataFrame()

    if rand_sampling:
        print(f'\n### A randomized sampling for each poylgon in spatial_distribution is performed to define the buildings which are included in the calculation of the line density and its spatial distribution.\n The target ratio of connection is set to {int(target_connection_ratio*100)} %.')

        # Spatial distribution of different connection ratios is desired
        if spatial_distribution is not None:

            # Initialisations
            spatial_distribution_out = spatial_distribution.copy()        
            spatial_distribution_out = pd.merge(spatial_distribution_out, polygons.groupby(by = spatial_uniqueID)[spatial_uniqueID].count().rename('nPolygons'), left_on = spatial_uniqueID, right_index = True)

            # Temporarily change matching with spatial distribution to 'NA' for polygons without match
            polygons[spatial_uniqueID] = polygons[spatial_uniqueID].fillna('NA')

            # for n, HEX in spatial_distribution_out.iterrows():
            for (spatial_ID, spatial_CR) in zip(list(spatial_distribution_out[spatial_uniqueID]) + ['NA'],list(spatial_distribution_out[spatial_connection_ratio]) + [target_connection_ratio]):

                # Initializations
                sampleList = list()
                res_arr_target_attr = np.zeros(nSamples)

                # Detect spatial distribution-individual connection ratio
                # CR = HEX[spatial_connection_ratio] if (spatial_connection_ratio is not None and spatial_connection_ratio in HEX.index) else target_connection_ratio
                CR = spatial_CR

                # Create filter for buildings which shall be INCLUDED and EXCLUDED from sampling
                # idxs_include = polygons[(polygons[spatial_uniqueID] == HEX[spatial_uniqueID])].index
                idxs_include = polygons[(polygons[spatial_uniqueID] == spatial_ID)].index

                if len(idxs_include) == 0:
                    continue

                if len(idxs_include) == 0:
                    fraction = 0

                else:
                    fraction = 1 if len(idxs_include) <= nPolygons_min else min(1, max(0, len(idxs_include) * CR  / len(idxs_include)))


                for jj in range(nSamples):
                    selection = polygons.loc[idxs_include, :].sample(frac = fraction)            

                    sampleList.append(selection)

                    res_arr_target_attr[jj]             = np.round(np.nansum(selection[target_attr]), 1)

                # Get randomized selection with median heat demand
                # Index of median (no. of sample featuring median heat demand)

                med                                                    = np.sort(res_arr_target_attr)[len(res_arr_target_attr)//2-1]
                idxmedian                                              = np.argmin(abs(res_arr_target_attr - med))

                # Use selected sample of buildings in network to assign heat demand to single pipe sections
                selected_polygons = pd.concat((selected_polygons, sampleList[idxmedian]))

        # No spatial distribution desired - target_connection ratio applies to all polygons        
        else:
            print(
                '\n### Line density and its spatial distribution is calculated using every building found in the examined area and with respect to the desired connection ratio. ###')

            # Initializations
            sampleList = list()
            res_arr_target_attr = np.zeros(nSamples)

            # Create filter for buildings which shall be INCLUDED and EXCLUDED from sampling
            idxs_include = polygons.index

            if len(idxs_include) == 0:
                fraction = 0

            else:
                fraction = 1 if len(idxs_include) <= nPolygons_min else min(1, max(0, (len(idxs_include)) * target_connection_ratio ) / len(idxs_include))


            for jj in range(nSamples):
                selection = polygons.loc[idxs_include, :].sample(frac = fraction)            

                sampleList.append(selection)

                res_arr_target_attr[jj]             = np.round(np.nansum(selection[target_attr]), 1)

            # Get randomized selection with median heat demand
            # Index of median (no. of sample featuring median heat demand)

            med                                                    = np.sort(res_arr_target_attr)[len(res_arr_target_attr)//2-1]
            idxmedian                                              = np.argmin(abs(res_arr_target_attr - med))

            # Use selected sample of buildings in network to assign heat demand to single pipe sections
            selected_polygons = pd.concat((selected_polygons, sampleList[idxmedian]))

    # Change back matching with spatial distribution from NA to None
    polygons.loc[polygons[spatial_uniqueID] == 'NA', spatial_uniqueID] = None
    selected_polygons.loc[selected_polygons[spatial_uniqueID] == 'NA', spatial_uniqueID] = None

    ### Transfer results for summed attributes to line sections
    polygons['usage_ld_calc'] = False # Initialisation
    polygons.loc[selected_polygons.index, 'usage_ld_calc'] = True

    selected_polygons_transfer = polygons[polygons['usage_ld_calc'] == True].copy()

    # Transfer results to line objects
    lines_out['nPolygons'] = selected_polygons_transfer.groupby('nearestline')[target_attr].count()
    lines_out['nPolygons'] = lines_out['nPolygons'].fillna(0)

    for sa in summed_target_attrs:
        # Line density (MWh/m)
        lines_out[f'summed_{sa}'] = selected_polygons_transfer.groupby('nearestline')[sa].sum()            

        # Remove nans
        lines_out[f'summed_{sa}'] = lines_out[f'summed_{sa}'].fillna(0)
    

    return lines_out, selected_polygons

# %%

if __name__ == '__main__':

    import os
    from pathlib import Path

    cs = 'EPSG:25832'

    ## Load example data
    flp = os.getcwd() / Path(r'examples/data/exampleNetwork/')
    lines = gp.read_file(flp / Path(r'streetsRaw.gpkg')).to_crs(cs)

    buildings = gp.read_file(flp / Path(r'buildings.gpkg')).to_crs(cs).drop(columns = 'nearestline')

    hex = gp.read_file(flp / Path(r'hex_div.gpkg')).to_crs(cs)

    # Define whether a randomized sampling in each hexagon with aim connection ratio shall take place
    rand_sampling = True

    # Number of randomized samples for selection of buildings to attain aimAG (sample with median heat demand is chosen)
    nSamples = 1

    # Define minimum number of buildings which shall be found within single hexagon below which the aimAG is ignored and ALL buildings are considered for calculation
    nBuildings_min = 1

    # Define aim of connection ratio within each subdivision of the regarded area in hex_div
    aimAG = 1

    # Define hexagon-specific connection ratio as attribute name of its value in hex
    attr_hex_CR = 'connectionRatio'

    # Define attribute from buildings for calculation of line density
    target_attr = 'demand_use_th'

    # Define additional attributes from the buildings layer which shall be summed on closest line objects
    summed_target_attrs = ['demand_use_th_2045_san']

    # Define maximum distance between buildings and line objects to include them into calculation of line density
    # -> Example uses function-like distance calculation
    def calc_distance_from_heat_demand(
        heat_demand, 
        min_line_density = 1000, # At least line density of 1 MWh/m at house connection line to assign closest line to polygon
        minLimit = 10,
        maxLimit = 100):

        return(min(maxLimit, max(minLimit, heat_demand/min_line_density)) )

    maxDistancesArray = np.array(list(map(lambda x: calc_distance_from_heat_demand(x, min_line_density = 1.7, minLimit = 10), np.array(buildings[target_attr]))))

    ### Start calculation

    lines_out, polys_out = summed_attributes_on_lines(
        polygons = buildings,
        lines = lines,
        spatial_distribution = hex,
        rand_sampling = rand_sampling,
        nSamples = nSamples,
        nPolygons_min = nBuildings_min,
        target_connection_ratio = aimAG,
        spatial_connection_ratio = attr_hex_CR,
        target_attr = 'demand_use_th',
        additional_attr = summed_target_attrs,
        func_max_distance = calc_distance_from_heat_demand
    )

lines_out.to_file(flp / Path(f'results/streetsRaw_lineDensity_calc_final.gpkg'), driver  ='GPKG')

polys_out.to_file(flp / Path(f'results/buildings_lineDensity_calc_final.gpkg'), driver = 'GPKG')



#     import os
#     from pathlib import Path

#     cs = 'EPSG:25832'

#     ## Load example data
#     flp = os.getcwd() / Path(r'examples/data/exampleNetwork/')
#     lines = gp.read_file(flp / Path(r'streetsRaw.gpkg')).to_crs(cs)

#     buildings = gp.read_file(flp / Path(r'buildings.gpkg')).to_crs(cs).drop(columns = 'nearestline')

#     hex = gp.read_file(flp / Path(r'hex_div.gpkg')).to_crs(cs)

#     # Define whether a randomized sampling in each hexagon with aim connection ratio shall take place
#     rand_sampling = True

#     # Number of randomized samples for selection of buildings to attain aimAG (sample with median heat demand is chosen)
#     nSamples = 10

#     # Define minimum number of buildings which shall be found within single hexagon below which the aimAG is ignored and ALL buildings are considered for calculation
#     nBuildings_min = 1

#     # Define aim of connection ratio within each subdivision of the regarded area in hex_div
#     aimAG = 1

#     # Define hexagon-specific connection ratio as attribute name of its value in hex
#     attr_hex_CR = 'connectionRatio'

#     # Define attribute from buildings for calculation of line density
#     target_attr = 'demand_use_th'

#     # Define additional attributes from the buildings layer which shall be summed on closest line objects
#     summed_target_attrs = []

#     # Define maximum distance between buildings and line objects to include them into calculation of line density
#     # -> Example uses function-like distance calculation
#     def calc_distance_from_heat_demand(
#         heat_demand, 
#         min_line_density = 1.5*1e3, 
#         minLimit = 10,
#         maxLimit = 100):

#         return(min(maxLimit, max(minLimit, heat_demand/min_line_density)) )

#     maxDistancesArray = np.array(list(map(lambda x: calc_distance_from_heat_demand(x, min_line_density = 1.7, minLimit = 10), np.array(buildings[target_attr]))))

#     ### Start calculation

#     # Initialisations
#     buildings_uniqueID = 'unique_ID'
#     lines_uniqueID = 'unique_ID_lines'
#     hex_uniqueID = 'unique_ID_hex'

#     buildings['unique_ID'] = np.arange(len(buildings)).astype(int)
#     lines[lines_uniqueID] = np.arange(len(lines)).astype(int)

#     summed_target_attrs = set(summed_target_attrs + [target_attr])
    
#     if hex is not None:
#         hex[hex_uniqueID] = np.arange(len(hex)).astype(int)

#         buildings = assign_attr_by_max_intersection_area(gp1 = buildings, gp_source=hex, gp1_id = buildings_uniqueID, attr=hex_uniqueID)

    
#     # Assign unique ID of closest line bject to buildings
#     buildings = closest_lines_to_polygons(polygons = buildings, lines = lines, maxDistances = maxDistancesArray)

#     # Transfer unique ID of lines to buildings instead of line index
#     dictIdxIDLines = dict(zip(list(lines.index), list(lines[lines_uniqueID])))
#     buildings['nearestline'] = buildings['nearestline'].apply(lambda x: dictIdxIDLines[x] if x in dictIdxIDLines.keys() else x)
    
#     # Create copies
#     lines_out = lines.copy()        

#     # Select buildings according to current filter
#     selected_builds = gp.GeoDataFrame()

#     if rand_sampling:
#         print(f'\n### A randomized sampling for each poylgon in spatial_distribution is performed to define the buildings which are included in the calculation of the line density and its spatial distribution.\n The target ratio of connection is set to {int(aimAG*100)} %.')

#         # if hex is not None:

#         # Initialisations
#         hex_out = hex.copy()        
#         hex_out = pd.merge(hex_out, buildings.groupby(by = hex_uniqueID)[hex_uniqueID].count().rename('nBuildings'), left_on = hex_uniqueID, right_index = True)

#         for n, HEX in hex_out.iterrows():

#             # Initializations
#             sampleList = list()
#             res_arr_heat_demand = np.zeros(nSamples)

#             # Detect hex-individual connection ratio
#             hex_CR = HEX[attr_hex_CR] if (attr_hex_CR is not None and attr_hex_CR in HEX.index) else aimAG

#             # Create filter for buildings which shall be INCLUDED and EXCLUDED from sampling
#             idxs_include = buildings[(buildings[hex_uniqueID] == HEX[hex_uniqueID])].index

#             if len(idxs_include) == 0:
#                 continue

#             if len(idxs_include) == 0:
#                 fraction = 0

#             else:
#                 fraction = 1 if len(idxs_include) <= nBuildings_min else min(1, max(0, len(idxs_include) * hex_CR  / len(idxs_include)))


#             for jj in range(nSamples):
#                 selection = buildings.loc[idxs_include, :].sample(frac = fraction)            

#                 sampleList.append(selection)

#                 res_arr_heat_demand[jj]             = np.round(np.nansum(selection[target_attr]), 1)

#             # Get randomized selection with median heat demand
#             # Index of median (no. of sample featuring median heat demand)

#             med                                                    = np.sort(res_arr_heat_demand)[len(res_arr_heat_demand)//2-1]
#             idxmedian                                              = np.argmin(abs(res_arr_heat_demand - med))

#             # Use selected sample of buildings in network to assign heat demand to single pipe sections
#             selected_builds = pd.concat((selected_builds, sampleList[idxmedian]))

#     # If rand_sampling == False then all buildings which meet the specified conditions are taken for the calculation of line density
#     else:
#         print(
#             '\n### Line density and its spatial distribution is calculated using every building found in the examined area and with respect to the desired connection ratio. ###')

#         # Initializations
#         sampleList = list()
#         res_arr_heat_demand = np.zeros(nSamples)

#         # Create filter for buildings which shall be INCLUDED and EXCLUDED from sampling
#         idxs_include = buildings.index

#         if len(idxs_include) == 0:
#             fraction = 0

#         else:
#             fraction = 1 if len(idxs_include) <= nBuildings_min else min(1, max(0, (len(idxs_include)) * aimAG ) / len(idxs_include))


#         for jj in range(nSamples):
#             selection = buildings.loc[idxs_include, :].sample(frac = fraction)            

#             sampleList.append(selection)

#             res_arr_heat_demand[jj]             = np.round(np.nansum(selection[target_attr]), 1)

#         # Get randomized selection with median heat demand
#         # Index of median (no. of sample featuring median heat demand)

#         med                                                    = np.sort(res_arr_heat_demand)[len(res_arr_heat_demand)//2-1]
#         idxmedian                                              = np.argmin(abs(res_arr_heat_demand - med))

#         # Use selected sample of buildings in network to assign heat demand to single pipe sections
#         selected_builds = pd.concat((selected_builds, sampleList[idxmedian]))

#     ### Transfer results for summed attributes to line sections
#     buildings['usage_ld_calc'] = False # Initialisation
#     buildings.loc[selected_builds.index, 'usage_ld_calc'] = True

#     selected_builds_transfer = buildings[buildings['usage_ld_calc'] == True].copy()

#     # Transfer results to line objects
#     lines_out['nBuilds'] = selected_builds_transfer.groupby('nearestline')[target_attr].count()
#     lines_out['nBuilds'] = lines_out['nBuilds'].fillna(0)

#     for sa in summed_target_attrs:
#         # Line density (MWh/m)
#         lines_out[f'summed_{sa}'] = selected_builds_transfer.groupby('nearestline')[sa].sum()            

#         # Remove nans
#         lines_out[f'summed_{sa}'] = lines_out[f'summed_{sa}'].fillna(0)

#     lines_out[f'ld_{target_attr}_per_m'] = lines_out[f'summed_{target_attr}']/lines_out.geometry.length

#     # lines_out.to_file(flp / Path(f'results/streetsRaw_lineDensity_calc.gpkg'), driver  ='GPKG')
#     # buildings.to_file(flp / Path(f'results/buildings_lineDensity_calc.gpkg'), driver = 'GPKG')




