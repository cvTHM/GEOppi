# -*- coding: utf-8 -*-

import pandas as pd
import geopandas as gp
import numpy as np

from geoppi.auxFunctions import (create_circumferential_points, closest_objects_to_points, assign_attr_by_max_intersection_area, )




def closest_lines_to_polygons(
        polygons:gp.GeoDataFrame,
        lines:gp.GeoDataFrame,
        maxDistances:np.array
        )->gp.GeoDataFrame:
    
    """
    Function that matches index of closest line object to each polygon provided in **polygons** based on spatial indexing with pbox and points along polygon boundaries.

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

def sum_attributes_on_lines(
        lines:gp.GeoDataFrame,
        polygons:gp.GeoDataFrame,
        spatial_distribution:gp.GeoDataFrame = None,
        spatial_connection_ratio:str = None,
        target_attr:str = None,
        agg_func:str = 'median',
        additional_attr:list = None,
        rand_sampling:bool = False,
        nSamples:int = 1,
        target_connection_ratio:float = 1,
        nPolygons_min:int = 1,
        func_max_distance = None
    )->gp.GeoDataFrame:

    """
    Function that sums up specified attributes of polygons on closest line segments. May be used to calculate line densities of heating or arbeitrary other demands.\n
    Allows for randomized sampling of chosen polygons selected for the final summary of attributes on line segments if desired connection ratio is not equal to 1.\n
    Allows for spatial distribution of different connection ratios in provided polygon objects *spatial_distribution*.\n

    :param lines: GeoDataFrame of line objects (e.g. street segments) on which attributes summary shall be trasnferred.\n
    :param polygons: GeoDataFrame of polygon objects from which attributes shall be transferred to closest line segments.\n
    :param spatial_distribution: GeoDataFrame of polygon objects containing information on spatially varying target connection ratios. Polygons are chose nbased on their matching object in *spatial_distribution* to determine if they are chosen for the calculation of summed attributes.\n
    :param spatial_connection_ratio: str denoting the attribute name in *spatial_distribution* in which infromation of target connection ratio in the objects is stored.\n
    :param target_attr: str denoting the attribute from *polygons* which shall be used to perform randomized sampling if target connection ratio is less than 1. NOTE: The sample with the median, min or max summed value is chosen based on the provided *agg_func*.\n
    :param agg_func: str denoting the final selection of samples, defaults to 'median'.\n
    :param additional_attr: list of str denoting additional attributes apart from *target_attr* that shall be summed and transferred to line segments.\n
    :param nSamples: int denoting the number of samples to perform.\n
    :param target_connection_ratio: float denoting the desired connection ratio of a) polygons which do not intersect with provided polygon objects from *spatial_distribution* or b) if *spatial_distribution* is not provided.\n
    :param nPolygons_min: int denoting a minimum number of polygons to perform randomized sampling. If fewer than *nPolygons_min* polygons are found to match with an object from *spatial_distribution* all matching polygons are chosen.\n
    :param func_max_distance: function-like that can be provided mapping minimum line densities of *target_attr* in the connection lines from *lines* to *polygons* to consider the polygons. If the line density is not sufficient in their connection lines (based on straight distance), the polygons are not taken intop account for summary of any of the attributes.\n

    :returns: lines_out denoting GeoDataFrame of edited line objects with summary of specified attributes AND selected_polygons as GeoDataFrame containting polygon objects used for the attribute summary.
    """

    ### Initializations
    cs = lines.crs

    polygons.to_crs(cs, inplace = True)

    if spatial_distribution is not None:
        spatial_distribution.to_crs(cs, inplace = True)

    # Check for availability of defined attributes
    if target_attr not in polygons.columns:
        print(f'\n### Attribute {target_attr} not found in polygons attributes. Aborting...')
        return
    
    if 'nearestline' in polygons.columns:
        nearestLineAttr = 'nearestline_1'
        print(rf'\n... Attention! Matching of closest line to polygons is done with attribute {nearestLineAttr}')
        
    else:
        nearestLineAttr = 'nearestline'
       
    summed_target_attrs = set([col for col in additional_attr if (col in polygons.columns) and (additional_attr is not None)] + [target_attr])

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
    polygons[nearestLineAttr] = polygons[nearestLineAttr].apply(lambda x: dictIdxIDLines[x] if x in dictIdxIDLines.keys() else x)
    
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
                if agg_func == 'median':
                    med                                                    = np.sort(res_arr_target_attr)[len(res_arr_target_attr)//2-1]                
                    idx_final                                              = np.argmin(abs(res_arr_target_attr - med))

                elif agg_func == 'max':
                    idx_final                                              = np.argmax(abs(res_arr_target_attr))

                elif agg_func == 'min':
                    idx_final                                              = np.argmin(abs(res_arr_target_attr))

                # Use selected sample of buildings in network to assign heat demand to single pipe sections
                selected_polygons = pd.concat((selected_polygons, sampleList[idx_final]))

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

            if agg_func == 'median':
                med                                                    = np.sort(res_arr_target_attr)[len(res_arr_target_attr)//2-1]                
                idx_final                                              = np.argmin(abs(res_arr_target_attr - med))

            elif agg_func == 'max':
                idx_final                                              = np.argmax(abs(res_arr_target_attr))

            elif agg_func == 'min':
                idx_final                                              = np.argmin(abs(res_arr_target_attr))

            # Use selected sample of buildings in network to assign heat demand to single pipe sections
            selected_polygons = pd.concat((selected_polygons, sampleList[idx_final]))


    else: # No randomized sampling
        print(f'\n... No randomized sampling is desired. A target conection ratio of 100% is applied.')
        selected_polygons = polygons.copy()

    # Change back matching with spatial distribution from NA to None
    if spatial_uniqueID in polygons.columns:    
        polygons.loc[polygons[spatial_uniqueID] == 'NA', spatial_uniqueID] = None
        selected_polygons.loc[selected_polygons[spatial_uniqueID] == 'NA', spatial_uniqueID] = None

    ### Transfer results for summed attributes to line sections
    polygons['usage_ld_calc'] = False # Initialisation
    polygons.loc[selected_polygons.index, 'usage_ld_calc'] = True

    selected_polygons_transfer = polygons[polygons['usage_ld_calc'] == True].copy()

    # Transfer results to line objects
    lines_out['nPolygons'] = selected_polygons_transfer.groupby(nearestLineAttr)[target_attr].count()
    lines_out['nPolygons'] = lines_out['nPolygons'].fillna(0)

    for sa in summed_target_attrs:
        # Line density (MWh/m)
        lines_out[f'summed_{sa}'] = selected_polygons_transfer.groupby(nearestLineAttr)[sa].sum()            

        # Remove nans
        lines_out[f'summed_{sa}'] = lines_out[f'summed_{sa}'].fillna(0)
    

    return lines_out, selected_polygons


