# %%
# -*- coding: utf-8 -*-


### Imports
import numpy as np
import pandas as pd
from pathlib import Path
import geopandas as gp
from auxFunctions import (assign_attr_by_max_intersection_area, closest_objects_to_points,)
from line_density_calculation import *
from create_network_topology import *
import logging


# %%


cs = "EPSG:25832"

### Load data
flp_in = Path(r"C:\Git\geoppi\geoppi\examples\data\example_gas_network")

streets = gp.read_file(flp_in / Path(r"ATKIS_DLM50_Straßenachse.gpkg")).to_crs(cs)
# streets_name = gp.read_file(r"Z:\41_Wärmeplanung_Hessen\06_GIS\ATKIS\Digitales Landschaftsmodell\ATKIS Basis DLM Hessen (shape)\dlm_42002_ax_strasse.shp", encoding ='utf-8').to_crs(cs)

buildings = gp.read_file(flp_in / Path(r"buildings_final.gpkg")).to_crs(cs)
parcels = gp.read_file(flp_in / Path(r"parcels_utilisation.gpkg")).to_crs(cs)

gas_network = gp.read_file(flp_in / Path(r"gas_network.gpkg")).to_crs(cs)
gas_network = gas_network.drop(columns = [n for n in gas_network.columns if not any([k in n.lower() for k in ("nennweite1", "druckstufe", "leitungsfu", "geometry")])]).reset_index(drop = True)
gas_network["id"] = gas_network.index.copy()

dh_expectation = gp.read_file(flp_in / Path(r"dh_expectation.gpkg")).to_crs(cs)


# %%
## Create convex hull and cut layers to its extent
hull = gas_network.unary_union.convex_hull
hull = gp.GeoDataFrame(geometry = [hull], crs = gas_network.crs)

streets = gp.sjoin(streets, hull, how = "inner", predicate  = "intersects")
streets = streets[~streets["index_right"].isna()].drop(columns = ["index_right"])

# streets_name = gp.sjoin(streets_name, hull, how = "inner", predicate  = "intersects")
# streets_name = streets_name[~streets_name["index_right"].isna()].drop(columns = ["index_right"])

buildings = gp.sjoin(buildings, hull, how = "inner", predicate  = "intersects")
buildings = buildings[~buildings["index_right"].isna()].drop(columns = ["index_right"])

parcels = gp.sjoin(parcels, hull, how = "inner", predicate  = "intersects")
parcels = parcels[~parcels["index_right"].isna()].drop(columns = ["index_right"])

gas_network = gp.sjoin(gas_network, hull, how = "inner", predicate  = "intersects")
gas_network = gas_network[~gas_network["index_right"].isna()].drop(columns = ["index_right"])


# %% Transfer of attributes to gas network lines

## Filter gas network to distribution lines
mask = [False if any([n in x.lower() for n in ("anschluss",)]) else True for x in gas_network["LEITUNGSFU"]]
gas_network_distribution = gas_network[mask]


## Cut gas network lines with itself
_, gas_network_distribution = createUniqueJunctions(gas_network_distribution.explode(index_parts = False), splitLines = True)



## Transfer street names to gas network sections
# streets_temp = streets_name.copy()
# streets_temp.geometry = streets_name.geometry.buffer(5)

# gas_network_distribution = assign_attr_by_max_intersection_area(gp1 = gas_network_distribution, gp_source = streets_temp, attr = ["nam_name", "wdm_widmun"], gp1_id = "id")


# %% Edit buildings layer before processing
import re

buildings["address"] = buildings["address"].replace("null", None)
buildings["street_name"] = buildings["address"].apply(lambda x: re.sub(r'^\s+', '', re.search(r'^([^ ]+)\s\d', x).group(1)) if isinstance(x, str) and re.search(r'^([^ ]+)\s\d', x) else x)


# %% Test Anliegerstraßen
def get_adjacent_lines(
    polygons:gp.GeoDataFrame,
    lines:gp.GeoDataFrame,
    dist:float = 100,
    nLines:int = 10,
    resolution:float = 2,
    blocking_polygons:gp.GeoDataFrame = None,
    return_connection_lines:bool = False,
    unique_ID_polys:str = 'unID_polys',
    unique_ID_lines:str = 'unID_lines',
    col_adjacentLineIDs:str = 'adjacentLine_IDs',
    col_adjacentDists:str = 'adjacentDists',
    col_adjacentConnectionLines:str = 'adjacentConnectionLines'
    ):

    from auxFunctions import nearest_points

    """
    Function to determine adjacent lines to provided polygons in such way that either the *nLines* nearest lines within *dist* are returned or only those where the direct connection from polygons to lines does not intersect any object from *parcels*.\n

    :param polygons: GeoDataFrame of polygons to determine adjacent lines for.\n
    :param lines: GeoDataFrame of lines to determine adjacent lines from.\n
    :param dist: Distance within which to search for adjacent lines.\n
    :param nLines: Maximum number of adjacent lines to return.\n
    :param resolution: Resolution of points along plolygon boundaries.\n
    :param blocking_polygons: GeoDataFrame of polygons to identify non-adjacent lines (blocking the way).\n
    :param return_connection_lines: Boolean to return connection lines.\n
    :param unique_ID_polys: Name of unique ID column in *polygons*.\n
    :param unique_ID_lines: Name of unique ID column in *lines*.\n
    :param col_adjacentLineIDs: Name of column to store IDs of adjacent lines.\n
    :param col_adjacentDists: Name of column to store distances to adjacent lines.\n
    :return: GeoDataFrame of polygons with adjacent lines as attribute.
    """

    ### Plausibility checks
    cs = polygons.crs
    if cs != lines.crs:
        lines = lines.copy().to_crs(cs)
        logging.warning(f'\nAttention! Polygons and lines do not share same crs!')

    if unique_ID_lines not in lines.columns:
        lines[unique_ID_lines] = np.arange(len(lines))
    if unique_ID_polys not in polygons.columns:
        polygons[unique_ID_polys] = np.arange(len(polygons))

    # Create temporary copies and reset indices
    lines = lines.copy()
    lines.reset_index(inplace = True, drop = True)
    dict_idx_idLines = dict(zip(lines.index, lines[unique_ID_lines]))

    ### Create circumferential points around polygons
    circumferentialPoints = create_circumferential_points(polygons = polygons, distance = resolution)
    circumferentialPoints[unique_ID_polys] = polygons[unique_ID_polys].copy()
    circumferentialPoints = circumferentialPoints.explode(index_parts = True).reset_index(drop = True)

    ### Find closest object to each point (indicating index of closest line segment)
    idxs, distances = closest_objects_to_points(points = circumferentialPoints, geomObjects = lines, nObjects = nLines, maxDist = np.ones(len(circumferentialPoints))*dist)
        
    ### Assign nearest line segments to each point at polys' boundaries
    circumferentialPoints['nearestLines']                            = idxs
    circumferentialPoints['distanceToLines']                         = distances
    
    connectionLines = [[LineString(nearest_points(geom1 = row.geometry, geom2 = lines.iloc[int(k)].geometry)) if k is not None and not (isinstance(k, float) and np.isnan(k)) else None for k in row["nearestLines"]] for _, row in circumferentialPoints.iterrows()]

    circumferentialPoints["connectionLines"] = connectionLines # Geometries of connection lines
    
    # Check lines for intersection with other polygons
    # If no intersection is found, empty list is returned []
    if blocking_polygons is not None:        

        circumferentialPoints["blocking_polygons"] = None
        circumferentialPoints["unblocked_line"] = [[True]*len(j) for j in circumferentialPoints["connectionLines"]]
        polys_sindex = polygons.sindex

        for idx, cL in circumferentialPoints.iterrows():
            blocking_polys_ids = []

            own_poly_id = cL[unique_ID_polys]

            valid_gs = cL["unblocked_line"].copy()

            for kk, g in enumerate(cL["connectionLines"]):
                blocking_polys_ids_single = []
                

                if g is None or g != g:
                    continue
                
                blocking_polys_candidates_idx = list(polys_sindex.query(g, predicate = "intersects"))

                for ps in blocking_polys_candidates_idx:
                    poly_row = polygons.iloc[ps]

                    if poly_row[unique_ID_polys] == own_poly_id:
                        continue

                    intersection = g.intersection(poly_row.geometry)

                    if intersection.length > 0:
                        blocking_polys_ids_single.append(poly_row[unique_ID_polys])
                        valid_gs[kk] = False

                blocking_polys_ids.append(blocking_polys_ids_single)

            circumferentialPoints.at[idx, "blocking_polygons"] = blocking_polys_ids
            circumferentialPoints.at[idx, "unblocked_line"] = valid_gs
            mask_unblocked = np.array(circumferentialPoints.at[idx, "unblocked_line"])
            circumferentialPoints.at[idx, "nearestLines"] = list(np.array(circumferentialPoints.at[idx, "nearestLines"])[mask_unblocked].flatten())
            circumferentialPoints.at[idx, "distanceToLines"] = list(np.array(circumferentialPoints.at[idx, "distanceToLines"])[mask_unblocked].flatten())
            circumferentialPoints.at[idx, "connectionLines"] = list(np.array(circumferentialPoints.at[idx, "connectionLines"])[mask_unblocked].flatten())

    
    groupedLineIDs = circumferentialPoints.groupby(by = unique_ID_polys)["nearestLines"].sum()
    groupedDists = circumferentialPoints.groupby(by = unique_ID_polys)["distanceToLines"].sum()    
    groupedConnectionLines = circumferentialPoints.groupby(by = unique_ID_polys)["connectionLines"].sum()

    ### Get unique lines belonging to each polygon
    adjacentLines, adjacentDists, adjacentConnectionLines = [], [], []


    breakpoint()
    
    for nn, ll in enumerate(groupedLineIDs.values):
        mins = {}
        lines_res = {}
        
        idxnoNaN = ~np.isnan(ll)
        ll = np.array(ll)[idxnoNaN] # Filter out nan line ids
        dd = np.array(groupedDists.values[nn])[idxnoNaN]
        cc = np.array(groupedConnectionLines.values[nn])[idxnoNaN] if return_connection_lines else np.full(len(ll), None)[idxnoNaN]

        for idx, (i, d) in enumerate(zip(ll, dd)):
            if d < mins.get(i, float("inf")):
                mins[i] = d
                lines_res[i] = cc[idx]

        if len(mins.keys()) == 0:
            uniqueLines, shortestDists = tuple(), tuple()
            uniqueLines_geoms = tuple()
        else:            
            uniqueLines, shortestDists = zip(*mins.items())
            _, uniqueLines_geoms = zip(*lines_res.items())

        adjacentLines.append(uniqueLines)
        adjacentDists.append(shortestDists)
        adjacentConnectionLines.append(uniqueLines_geoms)

   
    # Store in pandas DataFrame
    if not return_connection_lines and blocking_polygons is None:
        df_adjacentLines = pd.DataFrame(index = groupedLineIDs.index, data = {col_adjacentLineIDs: adjacentLines, col_adjacentDists: adjacentDists})
    else:
        df_adjacentLines = pd.DataFrame(index = groupedLineIDs.index, data = {col_adjacentLineIDs: adjacentLines, col_adjacentDists: adjacentDists, col_adjacentConnectionLines: adjacentConnectionLines})

        df_connectionLines = groupedConnectionLines.copy()
        df_connectionLines[:] = adjacentConnectionLines # Assign only nObjects shortest geometries of single LineStrings for each building
        df_connectionLines = df_connectionLines.explode()
        df_connectionLines = gp.GeoDataFrame(data = {unique_ID_polys:df_connectionLines.index}, geometry = df_connectionLines.values, crs = cs)
    

    # df_connectionLines["valid_adjacent_line_geoms"] = df_connectionLines["blocking_polygons"].apply(lambda x: isinstance(x, list) and len(x) == 0)
    # grouped2 = df_connectionLines[df_connectionLines["valid_adjacent_line_geoms"] == True].groupby(by = unique_ID_polys)["geometry"].agg(list)


    # df_adjacentLines[col_adjacentConnectionLines] = grouped2

   
    ### Assign results to polygons
    polygons = pd.merge(polygons, df_adjacentLines, how = 'left', left_on = unique_ID_polys, right_index = True)

    if return_connection_lines:
        return polygons, df_connectionLines
    else:
        return polygons


# %%


polys_out, connLines =  get_adjacent_lines(
    polygons = parcels.head(50),
    unique_ID_polys = "oid",
    lines = streets,
    return_connection_lines = True,
    blocking_polygons = parcels,
    resolution = 5,
    dist = 50,
    nLines = 5
)


# polys_out['adjacentLine_IDs'] = polys_out['adjacentLine_IDs'].astype(str)
# polys_out['adjacentDists'] = polys_out['adjacentDists'].astype(str)
# connLines['blocking_polygons'] = connLines['blocking_polygons'].astype(str)

polys_out.to_file(flp_in / Path(r"polys_test.gpkg"), driver = "GPKG")
connLines.to_file(flp_in / Path(r"connectionLines_test.gpkg"), driver = "GPKG")
# lines_out.to_file(flp_in / Path(r"lines_test.gpkg"), driver = "GPKG")


# %% Matching parcel IDs and adjacent line IDs to buildings

buildings_out = pd.merge(buildings, polys_out[[]])




# %%

gas_network_distribution, buildings_ed = sum_attributes_on_lines(
    lines = gas_network,
    polygons = buildings,
    target_attr = "demand_use_th_gas",
    rand_sampling = False,
    func_max_distance = lambda x: 200
)





