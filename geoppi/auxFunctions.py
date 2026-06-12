# -*- coding: utf-8 -*-
import geopandas as gp
import pandas as pd
import numpy as np
import math
import rasterio as rio
import networkx as nx
import libpysal
import warnings
from tqdm import tqdm
import shapely
import logging
from shapely.geometry import Point, MultiPoint, LineString, MultiLineString, MultiPolygon
from shapely.ops import transform
from shapely.strtree import STRtree
from shapely.prepared import prep

# %%

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

    circumferentialPoints["nearestLines_list"] = circumferentialPoints["nearestLines"].apply(lambda x: list(x) if not isinstance(x, list) else x)
    circumferentialPoints["distanceToLines_list"] = circumferentialPoints["distanceToLines"].apply(lambda x: list(x) if not isinstance(x, list) else x)
    
    connectionLines = [[LineString(nearest_points(geom1 = row.geometry, geom2 = lines.iloc[int(k)].geometry)) if k is not None and not (isinstance(k, float) and np.isnan(k)) else None for k in row["nearestLines_list"]] for _, row in circumferentialPoints.iterrows()]

    # connectionLines = [[LineString(nearest_points(geom1 = row.geometry, geom2 = lines.iloc[int(row["nearestLines"])].geometry)) if row["nearestLines"] is not None and not (np.isnan(row["nearestLines"])) else None] for _, row in circumferentialPoints.iterrows()]

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
            circumferentialPoints.at[idx, "nearestLines_list"] = list(np.array(circumferentialPoints.at[idx, "nearestLines_list"])[mask_unblocked].flatten())
            circumferentialPoints.at[idx, "distanceToLines_list"] = list(np.array(circumferentialPoints.at[idx, "distanceToLines_list"])[mask_unblocked].flatten())
            circumferentialPoints.at[idx, "connectionLines"] = list(np.array(circumferentialPoints.at[idx, "connectionLines"])[mask_unblocked].flatten())

    
    groupedLineIDs = circumferentialPoints.groupby(by = unique_ID_polys)["nearestLines_list"].sum() # Append all entries of line IDs from each circumferential point, grouped by polygon
    groupedDists = circumferentialPoints.groupby(by = unique_ID_polys)["distanceToLines_list"].sum()    
    groupedConnectionLines = circumferentialPoints.groupby(by = unique_ID_polys)["connectionLines"].sum()

    ### Get unique lines belonging to each polygon
    adjacentLines, adjacentDists, adjacentConnectionLines = [], [], []

    for nn, ll in enumerate(groupedLineIDs.values):
        mins = {}
        lines_res = {}
        
        idxnoNaN = ~np.isnan(ll)
        ll = np.array(ll)[idxnoNaN].astype(int) # Filter out nan line ids
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


        # Sort unique line IDs, geometries and distances with ascending distance values
        sort_idx = np.argsort(shortestDists)

        shortestDists = list(np.array(shortestDists)[sort_idx])
        uniqueLines = list(np.array(uniqueLines)[sort_idx])
        uniqueLines_geoms = list(np.array(uniqueLines_geoms)[sort_idx])

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
      
    ### Assign results to polygons
    polygons = pd.merge(polygons, df_adjacentLines, how = 'left', left_on = unique_ID_polys, right_index = True)

    if return_connection_lines:
        return polygons, df_connectionLines
    else:
        return polygons

def assign_attr_by_max_intersection_area(gp1:gp.GeoDataFrame, gp_source:gp.GeoDataFrame, attr:str, gp1_id:str='id'):

    '''
    Function to assign attr from GeoDataFrame gp_source to input-GeoDataFrame gp1 by maximum intersection area between each object of gp1 and gp_source.

    :param gp1: GeoDataFrame containing objects to which attributes from gp_source shall be matched.\n
    :param gp_source: GeoDataFrame containing attributes which shall be transferred to objects in gp1 with max. intersection.\n
    :param attr: str or list of strings denoting the attributes that shall be transferred.\n
    :param gp1_id: str denoting a unique identifier column name for objects in gp1.\n        
    '''

    import geopandas as gp

    # Mapping from ID to index in gp1
    dictemp = dict(zip(gp1[gp1_id], gp1.index))

    # Ensure attr is a list and restrict to columns that actually exist in gp_source
    if not isinstance(attr, list):
        if isinstance(attr, str):
            attr = [attr]

    attr = [at for at in attr if at in gp_source.columns]

    # Remove attributes which already exist in gp1
    ex_attr = [at for at in attr if at in gp1.columns]
    if len(ex_attr) > 0:
        print(f'\nAttributes {ex_attr} already contained in gp1. Columns are dropped and overwritten.')

    gp1 = gp1.copy()
    gp1.drop(columns=ex_attr, inplace=True, errors="ignore")

    # ---------- spatial pre-selection with spatial index ----------
    # Build spatial index from gp_source
    sidx = gp_source.sindex

    # For each object in gp1, get the potentially intersecting objects from gp_source
    # using bounding box intersection (fast)
    candidate_lists = gp1.geometry.bounds.apply(
        lambda b: list(sidx.intersection(b)), axis=1
    )

    # Reduce gp_source to the union of all candidate indices
    all_candidate_idx = sorted(set(i for lst in candidate_lists for i in lst))
    if len(all_candidate_idx) == 0:
        # No potential intersections: just append empty columns and return
        for at in attr:
            gp1[at] = None
        return gp1

    gp_source_reduced = gp_source.iloc[all_candidate_idx]

    # ---------- actual overlay only with preselected subset ----------
    gtemp = gp.overlay(
        gp1,
        gp_source_reduced[attr + ['geometry']],
        how='intersection',
        keep_geom_type=False
    )

    # Define metric for "maximum intersection"
    gtemp['temporary'] = gtemp.length

    if (
        (all(gp1.geom_type.isin(['LineString', 'MultiLineString'])) and
         all(gp_source.geom_type.isin(['Polygon', 'MultiPolygon'])))
        or
        (all(gp_source.geom_type.isin(['LineString', 'MultiLineString'])) and
         all(gp1.geom_type.isin(['Polygon', 'MultiPolygon'])))
    ):
        # Line–polygon combinations: use intersection length
        gtemp['temporary'] = gtemp.length
    elif (
        all(gp1.geom_type.isin(['Polygon', 'MultiPolygon'])) and
        all(gp_source.geom_type.isin(['Polygon', 'MultiPolygon']))
    ):
        # Polygon–polygon combinations: use intersection area
        gtemp['temporary'] = gtemp.area

    if gtemp.empty:
        # No real intersections -> return with empty attribute columns
        for at in attr:
            gp1[at] = None
        return gp1

    # For each gp1 object, keep the record with the largest intersection area/length
    idxmax = gtemp.groupby(by=gp1_id)['temporary'].idxmax()
    ls = gtemp.loc[idxmax.values, [gp1_id] + attr].set_index(gp1_id, drop=True)
    ls.index = [dictemp[p] for p in ls.index]

    # Write new attributes back to gp1 (aligned by index)
    gp1.loc[ls.index, attr] = ls[attr]

    return gp1

def assign_attr_by_max_intersection_area_agg(
    gp1:gp.GeoDataFrame, # Polygons
    gp_source:gp.GeoDataFrame,
    attr:str | list,
    gp1_id:str='id',
    gp_source_id:str | None = None,
    min_share:float=0.5,
    agg_func : str = "sum"
)->gp.GeoDataFrame:
    
    """
    Function to assign and aggregate *attr* from GeoDataFrame gp_source to input-GeoDataFrame gp1 by min. relative intersection area of X % between each object of gp1 and gp_source.

    :param gp1: GeoDataFrame containing objects to which attributes from gp_source shall be matched.\n
    :param gp_source: GeoDataFrame containing attributes which shall be transferred to objects in gp1 with max. intersection.\n
    :param attr: str or list of strings denoting the attributes that shall be transferred.\n
    :param gp1_id: str denoting a unique identifier column name for objects in gp1. If not provided, it is created with "id" as default.\n  
    :param gp_source_id: str denoting a unique identifier column name for objects in gp1. If not provided, it is created with "__source_id__" as default.\n   
    :param min_share: flaot denoting the minimum share of intersection area between metching polygons from gp_source to polygon in gp1 to transfer and aggregate *attr* on gp1.\n
    :param agg_func: str denoting the desired aggregation function (either 'sum' or 'mean').\n
    """

    # Plausbibility checks
    if agg_func not in ("sum", "mean", None):
        raise ValueError(f"\n... argument agg_func must be either 'sum' or 'mean' or None.")
    
    if not all(gp1.geom_type.isin(['Polygon', 'MultiPolygon'])):
        raise ValueError(f"\n... All objects from gp1 must be of type Polygon or MultiPolygon.")

    # --- Validate gp1_id ---
    if gp1_id not in gp1.columns:
        gp1[gp1_id] = np.arange(len(gp1))
        gp1[gp1_id] = gp1[gp1_id].astype(int)

    # --- Validate / create gp_source_id ---
    if gp_source_id is not None:
        if gp_source_id not in gp_source.columns:
            gp_source[gp_source_id] = np.arange(len(gp_source))
            gp_source[gp_source_id] = gp_source[gp_source_id].astype(int)
    else:
        gp_source = gp_source.copy()
        gp_source_id = "__source_id__"
        gp_source[gp_source_id] = np.arange(len(gp_source))

    # --- Ensure attr is a list and exists in gp_source ---
    if isinstance(attr, str):
        attr = [attr]
    attr = [a for a in attr if a in gp_source.columns]

    if not attr:
        raise ValueError("No valid attributes found in gp_source")

    # --- Remove existing attributes in gp1 ---
    ex_attr = [a for a in attr if a in gp1.columns]
    if ex_attr:
        print(f'Attributes {ex_attr} already in gp1 → overwritten')
        gp1 = gp1.drop(columns=ex_attr)

    # --- Precompute source areas (important for performance!) ---        
    if all(gp_source.geom_type.isin(['LineString', 'MultiLineString'])):  
        source_areas = gp_source.geometry.length
        all_sources_lines = True
    elif all(gp_source.geom_type.isin(['Polygon', 'MultiPolygon'])):
        source_areas = gp_source.geometry.area
        all_sources_lines = False

    # --- Build spatial index for gp1 ---
    sindex = gp1.sindex

    matches = []

    # --- Loop over gp_source geometries ---
    print(f"\n... Processing objects from source GeoDataFrame")
    for idx, geom in tqdm(zip(gp_source.index, gp_source.geometry), total = len(gp_source.index)):

        # Step 1: bounding box preselection
        possible_idx = list(sindex.intersection(geom.bounds))
        if not possible_idx:
            continue

        candidates = gp1.iloc[possible_idx]

        # Step 2: compute intersections
        inter = candidates.geometry.intersection(geom)

        # Remove empty intersections
        mask = ~inter.is_empty
        if not mask.any():
            continue

        inter = inter[mask]
        candidates = candidates.loc[mask]

        # Step 3: compute intersection areas
        areas = inter.area if not all_sources_lines else inter.length

        # --- NEW: filter by minimum share of source polygon ---
        share = areas / source_areas.loc[idx]

        valid_mask = share >= min_share
        if not valid_mask.any():
            continue

        areas = areas[valid_mask]
        candidates = candidates.loc[valid_mask]

        # Step 4: select gp1 polygon with maximum intersection
        max_idx = areas.idxmax()
        gp1_match_id = candidates.loc[max_idx, gp1_id]

        # Step 5: store mapping + attributes from gp_source
        row = {gp1_id: gp1_match_id}
        for a in attr:
            row[a] = gp_source.loc[idx, a]

        matches.append(row)

    # --- Convert matches to DataFrame ---
    df_matches = pd.DataFrame(matches)

    if df_matches.empty:
        print("No intersections found after applying share threshold.")
        return gp1

    # --- Aggregate attributes by gp1 polygon ---
    if agg_func == "sum":
        df_agg = df_matches.groupby(gp1_id)[attr].sum()
    elif agg_func == "mean":
        df_agg = df_matches.groupby(gp1_id)[attr].mean()
    elif agg_func == None:
        df_agg = df_matches.groupby(gp1_id)[attr].first()

    # --- Assign aggregated values back to gp1 ---
    gp1 = gp1.set_index(gp1_id)
    gp1[attr] = df_agg
    gp1 = gp1.reset_index()

    return gp1



def calc_intersection_area(
    gp_source:gp.GeoDataFrame,
    gp1:gp.GeoDataFrame,
    new_col: str = "intersection_area"
):
    
    """
    Fucntion that takes two GeoDataFrames *gp1* and *gp_source* and calculates the sum of intersection areas between overlaping elements of these GeoDataFrames.\n

    :param gp_source: GeoDataFrame of Polygon OR Line objects shall be checked for intersections with *gp1*.\n
    :param gp1: GeoDataFrame in which tio transfer intersection areas of each element with elemnets from *gp_source.\n
    :param new_col: str denoting the target attribute name in gp1 where to store values of intersection areas.\n

    :return: GeoDatFrame
    """

    if all([k in (['LineString', 'MultiLineString']) for k in gp_source.geom_type]):  
        all_sources_lines = True
    else:
        all_sources_lines = False

    # Build spatial index on SOURCE (more efficient this way)
    sindex = gp_source.sindex

    # Prepare result array
    result = np.zeros(len(gp1))

    print("\n... Processing gp1 geometries")

    for i, (_, geom) in enumerate(tqdm(zip(gp1.index, gp1.geometry), total=len(gp1))):

        # Get candidate geometries from gp_source using bounding box
        possible_idx = list(sindex.intersection(geom.bounds))
        if not possible_idx:
            continue

        candidates = gp_source.iloc[possible_idx]

        # Compute actual intersections
        inter = candidates.geometry.intersection(geom)

        # Remove empty geometries
        inter = inter[~inter.is_empty]
        if len(inter) == 0:
            continue

        # Sum up intersection areas
        if all_sources_lines:
            result[i] = inter.length.sum()
        else:
            result[i] = inter.area.sum()

    # Attach result as new column
    gp1[new_col] = result

    return gp1



def nearest_points(geom1, geom2):
    """
    Function to determine the two closest points of two input geometries..\n
    
    :param geom1: first input geometry.\n
    :param geom2: second input geometry.\n
    :returns: tuple containing the calculated nearest points in the input geometries.\n
    The points are returned in the same order as the input geometries.
    """
    seq = shapely.shortest_line(geom1, geom2)
    if seq is None:
        if geom1.is_empty:
            raise ValueError("The first input geometry is empty")
        else:
            raise ValueError("The second input geometry is empty")

    p1 = shapely.get_point(seq, 0)
    p2 = shapely.get_point(seq, 1)
    return (p1, p2)



def round_coords_2D(geom, ndigits:int=3):
    """
    Function that rounds entire geometries of 2D LineString objects to desired number of decimal digits.\n
    
    :param geom: LineString geometry.\n
    :param ndigits: number of decimals after rounding operation.\n
    :return: roundeed geometry.\n
    """
    def _round(x, y, z=None):
        if z is None:
            return (round(x, ndigits), round(y, ndigits))
        else:
            return (round(x, ndigits), round(y, ndigits), round(z, ndigits))
    
    return transform(_round, geom)

def create_circumferential_points(polygons:gp.GeoDataFrame, distance = None, npoints:int = 10, lmin:float = 5)->gp.GeoDataFrame:

    """
    Function creating MultiPoint object with circumferential points around boundary of each polygon object in **polygons**.

    :param polygons: GeoDataFrame of polygon objects.\n
    :param npoints: integer defining the desired number of points to create per polygon.\n
    :param lmin: float denoting the minimum distance between points to preserve.\n

    :return: GeoDataFrame containting MultiPoint object for each polygon.
    """

    # Imports
    from tqdm import tqdm

    # Create boundary line
    bounds = polygons[['geometry']].copy()
    bounds['geometry'] = bounds.geometry.boundary
   
    print('\n### Generating circumferential points ... ###\n')
    bounds['circumferentialPoints'] = bounds.geometry.apply(lambda x: create_point_splitter_npoints(x, npoints = npoints, distance = distance, lmin = lmin))  

    points = gp.GeoDataFrame(bounds[['circumferentialPoints']].rename(columns = {'circumferentialPoints':'geometry'}), geometry = 'geometry')

    # points = gp.GeoDataFrame(geometry = [create_points_along_boundary(p, dist = distance, return_Point_object = True) for p in tqdm(polygons.geometry)], crs = polygons.crs)

    return points


def create_points_along_boundary(polygon, return_Point_object:bool = False, dist:float = 1):

        """
        Function creating circumferential points along boundary of a polygon object.\n

        :param polygon: shapely.Polygon object.\n
        :param dist: float denoting the desired distance between circumferential points along boundary to create.\n

        :return: list of points along boundary of the polygon
        """

        boundary = polygon.boundary
        if boundary.length == 0:
            return []
        
        numPoints = int(boundary.length // dist)
        if return_Point_object:
            return MultiPoint([Point(np.round(boundary.interpolate(ii * dist).x, 2), np.round(boundary.interpolate(ii * dist).y, 2)) for ii in range(numPoints+1)])
        else:
            return [(np.round(boundary.interpolate(ii * dist).x, 2), np.round(boundary.interpolate(ii * dist).y, 2)) for ii in range(numPoints+1)]

def create_point_splitter_npoints(line, npoints:int = 8, distance:float = None, lmin:float = 5):

    '''
    Create MultiPoint layer with points along idx_ps_min_distances given LineString geometry at idx_ps_min_distances predefined number of points (respecting minimum values for single line segments between points).

    :param line: shaely.goemetry line object.\n
    :param npoints: int denoting the number of points that shall be created along line.\n
    :param lmin: float denoting the minimum length of line segments between consecutive points along the line.\n
    
    :returns: shapely.MultiPoint object for splitting points on the provided line object.
    '''
    import numpy as np
    from shapely.ops import unary_union

    if distance is None:
        distances = np.round(np.linspace(0, line.length, npoints+1), 1)
    else:
        distances = np.round(np.arange(0, line.length, distance), 1)
    
    # If last element is shorter than lmin, skip last splitting point
    if line.length - distances[-1] < lmin:
        distances = distances[:-1]


    points = [line.interpolate(d) for d in distances]
    result = unary_union(points)

    return(result)


def closest_objects_to_points(points:gp.GeoDataFrame, geomObjects:gp.GeoDataFrame, nObjects:int = 1, **kwargs)->tuple[np.array, np.array]:

    """
    Function that returns the index of the closest geometry object to each point in **points**.\n
    A kwarg can be "maxDist", either as array for each point in **points** or as a constant float/integer applied to all points.\n
    If no matching object within the maximum distance is found, np.nan is entered.\n
    If nObjects > 1, found objects within the provided threshold distance are returned as a list with ascending distance.\n

    :param points: GeoDataFrame of point objects for which to find closest geometry object in **geomObjects**\n
    :param geomObjects: GeoDataFrame containting geometry objects of arbitrary type in whcih to look for closest objects using *shapely.distance*.\n
    :param nObjects: int (optional), denoting the max. number of geometric objects to return for each point.\n
    :kwargs: if 'maxDist' is contained in kwargs.keys(), the maximum searching distance for the closest geometry object is restricted to this value. Can either be provided point-wise ias an array or as a constant value (float, integer) for all points. If no matching object within the maximum distance is found, np.nan is entered.\n

    :return: Two numpy.arrays with (i) indicating the matching index of the closest object found in **geomObjects** and (ii) the corresponding distances.
    """

    # Checking kwargs
    defaultDist = 100

    if 'maxDist' in kwargs.keys():
        maxDist = kwargs['maxDist']

        if type(maxDist) == float | type(maxDist) == int:
            maxDist = np.ones(points.shape[0]) * maxDist

    else:
        maxDist = np.ones(points.shape[0]) * defaultDist
        print('\n### No array for maximum searching distance is defined. Using default value of {} m. ###'.format(defaultDist))

    # Initialisation of output
    result = []
    dists = []

    # Initialise parameters for distance search
    INFTY = 1000000000000

    print('\n### Generating rtree... ###')
    idx = geomObjects.sindex

    for ix, r in tqdm(points.iterrows(), total = points.shape[0], desc = 'Snapping in pBox...'):

        p = r.geometry

        MIN_SIZE = maxDist[ix]
        max_dist = MIN_SIZE

        pbox = (p.x - MIN_SIZE, p.y - MIN_SIZE, p.x + MIN_SIZE, p.y + MIN_SIZE)

        hits = list(idx.intersection(pbox))
        d = INFTY
        nid = None

        if nObjects == 1:
            for h in hits:
                new_d = p.distance(geomObjects.geometry.loc[h])

                if (d >= new_d) & (new_d <= max_dist):
                    d = new_d
                    nid = geomObjects.index[h]

            if nid is None:
                result.append(np.nan)
                dists.append(d)

            else:
                result.append(nid)
                dists.append(d)


        elif nObjects > 1: # Multiple closest objects shall be returned
            nid_list = []
            d_list = []

            for h in hits:
                new_d = p.distance(geomObjects.geometry.loc[h])

                if (new_d <= max_dist):
                    d = new_d
                    nid = geomObjects.index[h]

                    nid_list.append(nid)
                    d_list.append(d)

            if len(nid_list) == 0:
                result.append([np.nan])
                dists.append([np.nan])

            else:
                sort_idx = np.argsort(d_list)
                nid_list = np.array(nid_list)[sort_idx]
                d_list = np.array(d_list)[sort_idx]

                result.append(list(nid_list[:nObjects]))
                dists.append(list(d_list[:nObjects]))       

    return result, dists

def calc_thermalLoss_pipe(
        net,
        inplace:bool = False
    ):

    """
    Function that calculates pipe-specific thermal loss power in pandapipes network model.\n

    :param net: pandapipes network model with existing res_pipe DataFrame (available after thermal pipeflow).\n
    :param inplace: boolean if results shall be directly transferred to results DataFrame. Then, the net object is returned. If False, the array of thermal loss for all pipes is returned.\n
    :return: network instance
    """

    # Initializations
    if hasattr(net, 'res_pipe'):
        mdot = net.res_pipe['mdot_from_kg_per_s'].values

        reverseFlow = np.where(mdot < 0)
        mask = np.ones(mdot.shape[0], dtype=bool)
        mask[reverseFlow] = False

        tFrom = net.res_pipe['t_from_k'].values
        tFrom[~mask] = net.res_pipe['t_to_k'].values[~mask]

        tTo = net.res_pipe['t_outlet_k'].values

        cp = net.fluid.get_heat_capacity((tFrom + tTo)/2)

        qloss = abs(mdot) * cp * (tTo - tFrom)

        if inplace:
            net.res_pipe['Pthermal_W'] = qloss
            return net
        else:
            return qloss

    else:
        print(f"\nnet.res_pipe is not found...")
        if inplace:
            return net
        else:
            return None



def extractPointsFromLines(
        lines:gp.GeoDataFrame,
        onlyIntermediate:bool = False
        ):
    
    """
    Function that extracts points from given GeoDataFrame of line objects.

    :param lines: GeoDataFrame containing LineString objects
    :param onlyIntermediate: Boolean if only intermediate points (no start- and endpoints) shall be extracted

    :return: tuple of lists with (shapely.Point objects, coordinates of points)
    """

    ### Imports
    from shapely.geometry import Point

    # Gather all points in list
    pointsList = []
    coordsList = []

    # Only intermediat points are extracted
    if onlyIntermediate:
        for row, pi in lines.iterrows():
            
            if len(pi.geometry.coords) > 2: # Intermediate points only exist when more then one start- and endpoint are existent in line geometry
                coords_i = [Point(pcoords) for pcoords in list(pi.geometry.coords)[1:-1]]
                coords_j = [pcoords for pcoords in list(pi.geometry.coords)[1:-1]]
        
                pointsList += coords_i
                coordsList += coords_j

    # All points are extracted
    else:
        for row, pi in lines.iterrows():
            
            coords_i = [Point(pcoords) for pcoords in list(pi.geometry.coords)[:]]
            coords_j = [pcoords for pcoords in list(pi.geometry.coords)[:]]
    
            pointsList += coords_i
            coordsList += coords_j

    return pointsList, coordsList


def measure_street_width(
    axis_lines: gp.GeoDataFrame,
    street_parcels: gp.GeoDataFrame,
    step_size: float = 1.0,
    max_range: float = 50.0,
    max_valid_width: float = 15.0,
    fan_angles: list = None,
    overlap_threshold: int = 0,
    min_boundary_angle: float = 0.0,
):
    """
    Function that measures street width by casting rays from sample points along given axis lines against the boundary of the merged street parcel polygon.\n

    :param axis_lines: GeoDataFrame containing street axis LineString objects.\n
    :param street_parcels: GeoDataFrame containing street parcel polygons (merged internally via union_all()).\n
    :param step_size: float denoting the distance between sample points along each axis line.\n
    :param max_range: float denoting the maximum ray length in meters.\n
    :param max_valid_width: float denoting the maximum valid total width; measurements >= this value are discarded.\n
    :param fan_angles: list of angles in degrees relative to the perpendicular; per point the angle yielding the smallest total width wins. Defaults to [0.0].\n
    :param overlap_threshold: integer; if > 0, a measurement line is dropped when it intersects with more than this many other measurement lines (filters chaotic measurements).\n
    :param min_boundary_angle: float; if > 0, a measurement line is dropped if its endpoint hits the street boundary outside [min_boundary_angle, 180 - min_boundary_angle] degrees.\n

    :return: tuple (axis_lines, profile_gdf, lines_gdf) of GeoDataFrames; axis_lines holds original line gdf, supplemented with attribute of narrowest passage found in its course, profile_gdf holds sample points with width attributes, lines_gdf holds the measurement lines.
    """
    if fan_angles is None:
        fan_angles = [0.0]
    if axis_lines.crs != street_parcels.crs:
        street_parcels = street_parcels.to_crs(axis_lines.crs)

    # Prepare merged street polygon (one geometry instead of many parcels)
    road          = street_parcels.geometry.union_all()
    road_prep     = prep(road)
    road_boundary = road.boundary

    # Create temp copy
    axis_lines = axis_lines.copy()

    # ---- nested helpers ----

    def _flatten_lines(geoms):
        """Flatten an iterable of LineString/MultiLineString geometries into a list of plain LineStrings."""
        out = []
        for g in geoms:
            if g is None or g.is_empty:
                continue
            if g.geom_type == 'MultiLineString':
                out.extend(g.geoms)
            else:
                out.append(g)
        return out

    def _tangent(line, dist, delta=0.1):
        """Unit tangent vector of line at position dist (None if direction is undefined)."""
        a = line.interpolate(max(0, dist - delta))
        b = line.interpolate(min(line.length, dist + delta))
        dx, dy = b.x - a.x, b.y - a.y
        L = math.hypot(dx, dy)
        if L == 0:
            return None
        return dx / L, dy / L

    def _rotate(vec, angle_deg):
        """Rotate a 2D vector by angle_deg degrees (counter-clockwise)."""
        a = math.radians(angle_deg)
        c, s = math.cos(a), math.sin(a)
        return vec[0] * c - vec[1] * s, vec[0] * s + vec[1] * c

    def _cast_ray(start, vec):
        """Distance from start in direction vec until the boundary is hit (or max_range)."""
        end = (start.x + vec[0] * max_range, start.y + vec[1] * max_range)
        ray = LineString([(start.x, start.y), end])
        hits = ray.intersection(road_boundary)
        if hits.is_empty:
            return max_range
        if hits.geom_type == 'Point':
            points = [hits]
        elif hits.geom_type == 'MultiPoint':
            points = list(hits.geoms)
        else:
            return max_range
        distances = [start.distance(p) for p in points if start.distance(p) > 1e-6]
        return min(distances) if distances else max_range

    def _ray_crosses_other_axis(start, vec, length):
        """True if the ray crosses another axis line further than 5 cm from its start."""
        if length < 0.01 or axis_tree is None:
            return False
        end = (start.x + vec[0] * length, start.y + vec[1] * length)
        ray = LineString([(start.x, start.y), end])
        for j in axis_tree.query(ray, predicate='intersects'):
            inter = ray.intersection(axis_geoms[int(j)])
            if not inter.is_empty and inter.distance(start) > 0.05:
                return True
        return False

    def _filter_overlap(profile, lines):
        """Drop measurement lines that intersect with >= overlap_threshold other measurement lines."""
        geoms = list(lines.geometry)
        tree = STRtree(geoms)
        keep = []
        for i, g in enumerate(geoms):
            n = sum(1 for j in tree.query(g, predicate='intersects') if int(j) != i)
            keep.append(n < overlap_threshold)
        return profile.loc[keep].reset_index(drop=True), lines.loc[keep].reset_index(drop=True)

    def _filter_boundary_angle(profile, lines, search_radius=0.5):
        """Drop measurement lines whose endpoint hits the street boundary at an angle outside [min_boundary_angle, 180 - min_boundary_angle]."""
        max_angle = 180.0 - min_boundary_angle
        polys = list(road.geoms) if isinstance(road, MultiPolygon) else [road]
        starts_list, vecs_list = [], []
        for poly in polys:
            for ring in [poly.exterior, *poly.interiors]:
                coords = np.asarray(ring.coords)
                starts_list.append(coords[:-1])
                vecs_list.append(np.diff(coords, axis=0))
        starts = np.vstack(starts_list)
        vecs   = np.vstack(vecs_list)
        seg_geoms = [LineString([s, s + v]) for s, v in zip(starts, vecs)]
        seg_tree  = STRtree(seg_geoms)

        keep = []
        for line in lines.geometry:
            c = np.asarray(line.coords)
            m_vec  = c[1] - c[0]
            m_norm = np.linalg.norm(m_vec)
            ok = True
            if m_norm > 1e-10:
                for endpoint in c:
                    ep = Point(endpoint)
                    cand = seg_tree.query(ep.buffer(search_radius), predicate='intersects')
                    if len(cand) == 0:
                        continue
                    dists = [seg_geoms[int(j)].distance(ep) for j in cand]
                    nearest = int(cand[int(np.argmin(dists))])
                    b_vec  = vecs[nearest]
                    b_norm = np.linalg.norm(b_vec)
                    if b_norm < 1e-10:
                        continue
                    cos_a = np.clip(np.dot(m_vec, b_vec) / (m_norm * b_norm), -1, 1)
                    angle = math.degrees(math.acos(cos_a))
                    if angle < min_boundary_angle or angle > max_angle:
                        ok = False
                        break
            keep.append(ok)
        return profile.loc[keep].reset_index(drop=True), lines.loc[keep].reset_index(drop=True)

    # ---- setup ----

    axis_geoms = _flatten_lines(axis_lines.geometry)
    axis_tree  = STRtree(axis_geoms) if axis_geoms else None

    print(f"\n... Starting measurement with {len(fan_angles)} ray(s) per point")
    print(f"... Fan angles in use: {np.round(fan_angles, 1)}")

    profile_rows, line_rows = [], []

    # ---- main loop ----

    for road_id, axis in tqdm(zip(axis_lines.index, axis_lines.geometry),
                              total=len(axis_lines), desc="Measuring streets"):
        if axis is None or axis.length == 0:
            continue

        for dist in np.arange(0, axis.length, step_size):
            point = axis.interpolate(dist)
            if not road_prep.intersects(point):
                continue

            tangent = _tangent(axis, dist)
            if tangent is None:
                continue
            normal = (-tangent[1], tangent[0])

            # Pick the SMALLEST valid total width across all fan angles
            best = None
            for angle in fan_angles:
                v_left  = _rotate(normal, angle)
                v_right = (-v_left[0], -v_left[1])

                d_left  = _cast_ray(point, v_left)
                d_right = _cast_ray(point, v_right)

                if (_ray_crosses_other_axis(point, v_left,  d_left) or
                    _ray_crosses_other_axis(point, v_right, d_right)):
                    continue

                total = d_left + d_right
                if best is None or total < best[0]:
                    best = (total, d_left, d_right, v_left)

            if best is None or best[0] >= max_valid_width:
                continue

            total, d_left, d_right, v_left = best
            p_left  = (point.x + v_left[0] * d_left,  point.y + v_left[1] * d_left)
            p_right = (point.x - v_left[0] * d_right, point.y - v_left[1] * d_right)

            profile_rows.append({
                'geometry':      point,
                'strasse_id':    road_id,
                'station':       round(dist, 2),
                'breite_links':  round(d_left, 2),
                'breite_rechts': round(d_right, 2),
                'breite_m':      round(total, 2),
            })
            line_rows.append({
                'geometry':   LineString([p_left, p_right]),
                'strasse_id': road_id,
                'station':    round(dist, 2),
                'breite_m':   round(total, 2),
            })

    profile_gdf = gp.GeoDataFrame(profile_rows, crs=axis_lines.crs)
    lines_gdf   = gp.GeoDataFrame(line_rows,    crs=axis_lines.crs)

    # ---- optional post-filters ----

    if overlap_threshold > 0 and len(lines_gdf):
        print("... Applying overlap filter")
        before = len(lines_gdf)
        profile_gdf, lines_gdf = _filter_overlap(profile_gdf, lines_gdf)
        print(f"    removed: {before - len(lines_gdf)}")

    if min_boundary_angle > 0 and len(lines_gdf):
        print(f"... Applying boundary-angle filter (>= {min_boundary_angle} deg)")
        before = len(lines_gdf)
        profile_gdf, lines_gdf = _filter_boundary_angle(profile_gdf, lines_gdf)
        print(f"    removed: {before - len(lines_gdf)}")

    # Transfer results with narrowest passage to original line DataFrame
    axis_lines["narrowPassage_m"] = profile_gdf.groupby(by = "strasse_id").agg({"breite_m":"min"})

    return axis_lines, profile_gdf, lines_gdf

def detect_lines_in_narrow_passages(
    lines:gp.GeoDataFrame,
    polygons:gp.GeoDataFrame,
    merge_touching_polygons:bool = False,
    threshDistance:float = 10,
    distPointsCircumference:float = 1,
    nNeighbours:int = 10,
    col:str = 'narrowPassage_m' 
    )->gp.GeoDataFrame:

    """
    Function that marks line objects if they are positioned between polygons with smaller distance than defined thresh value.\n

    :param lines: geopandas.GeoDataFrame with line objects (e.g. streets, ...).\n
    :param polygons: geopandas.GeoDataFrame with polygon objects (e.g. buildings, ...).\n
    :param merge_touching_polygons: Boolean selection if touching polygons shall be merged before creation of closest connecting lines between different polygons within threshDistance.\n
    :param threshDistance: float denoting the defined minimum distance between polygon outer boundaries between which lines are marked.\n
    :param distPointsCircumference: float denoting the distance of points that are created along the boundaries of polygons (!Affects quality of the result!).\n
    :param nNeighbours: integer denoting the number of neighbouring points to search within the defined thresh distance (!Affects quality of the result!).\n
    :param col: string denoting the column name in which the distance within the narrowest passage is entered for original line objects.\n

    :return: geopandas.GeoDataFrame of line objects with additional column *col*
    """

    # Imports
    import libpysal
    import numpy as np
    from shapely import LineString

    # Plausibility checks
    nNeighbours = max(1, nNeighbours)

    cs = lines.crs

    # Definitions
    def qweights_wrapper(df, col : str = None):

        if col is None:
            col = 'qweights'

        # testing: This does broadcast
        # queen_weights = libpysal.weights.Queen.from_dataframe(df)
        # return (queen_weights.component_labels + df.index[0]).astype(int)

        # ... maybe this does?
        queen_weights = libpysal.weights.Queen.from_dataframe(df)
        df[col] = (queen_weights.component_labels + df.index[0]).astype(int)
        return df


    def create_points_along_boundary(polygon, dist:float = 1):

        """
        Function creating circumferential points along boundary of a polygon object.\n

        :param polygon: shapely.Polygon object.\n
        :param dist: float denoting the desired distance between circumferential points along boundary to create.\n

        :return: list of points along boundary of the polygon
        """

        boundary = polygon.boundary
        if boundary is None or boundary.is_empty or boundary.length == 0:
            return []
        
        numPoints = int(boundary.length // dist)
        return [(np.round(boundary.interpolate(ii * dist).x, 2), np.round(boundary.interpolate(ii * dist).y, 2)) for ii in range(numPoints+1)]
           
    ## Data preparation
    # Merge all touching/overlapping polygons to reduce number of processed polygons
    if merge_touching_polygons:
        print(f'\n... Merging touching/overlapping polygons\n')
        polygons = qweights_wrapper(polygons, col = 'qWeight')
        polygons = polygons.dissolve(by = 'qWeight', aggfunc = 'first')

    # Create circumferenital points along boundary of polygons
    all_points = [create_points_along_boundary(poly, dist = distPointsCircumference) for poly in polygons.geometry]

    # Initialize container for shortest line objects between neighbouring points
    shortestLines = []

    ## Start looping through all points
    for n, pts in enumerate(all_points):

        # Initialize container for lines for each polygon
        singleLines = []

        # Exclude points belonging to the same polygon from the search for neighbours
        otherpts = [pt for j in range(len(all_points)) if j != n for pt in all_points[j]]
        lOthers = len(otherpts)

        if lOthers > 0 and len(pts) > 0:
            idx, _ = nnearest(A = np.array(pts), B = np.array(otherpts), distance = threshDistance, n = nNeighbours)

            for nn, p in enumerate(idx):
                singleLines.append( [(pts[nn], otherpts[ix]) for ix in p if ix < lOthers] )
                singleLines = [b for b in singleLines if b]

            shortestLines.append(singleLines)

    # Flatten list of lines
    shortestLines = [x for xs in shortestLines for x in xs]
    shortestLineObjects = gp.GeoDataFrame(geometry = [LineString(coords) for ls in shortestLines for coords in ls])
    shortestLineObjects.set_crs(cs, inplace = True)
    
    # Order all geometries in canonical form and drop duplicates (Connecting lines feature only one start- and one endpoint)
    shortestLineObjects['geometry'] = shortestLineObjects.normalize()
    shortestLineObjects.drop_duplicates()

    shortestLineObjects['length_m'] = shortestLineObjects.geometry.length

    # Transfer results to output DF of lines
    lines = lines.copy()

    joined = gp.sjoin(lines, shortestLineObjects[['geometry', 'length_m']], how = 'left', predicate = 'intersects')
    joined['index_left'] = joined.index
    joined.reset_index(drop = True, inplace = True)
    valid = joined.dropna(subset=['length_m'])
    idxs_min = valid.groupby('index_left')['length_m'].idxmin().dropna()
    lines.loc[idxs_min.index, col] = joined.loc[idxs_min.values, 'length_m'].values
    shortestLines = shortestLineObjects.loc[joined.loc[idxs_min.values, 'index_right'].values]

    return lines, shortestLines


def sort_with_permutation(lst:list, key, **kwargs):

    """
    Function that applies sorting to a list and returns sorted list and permutation order of this sorted list.\n

    :param lst: list to be sorted.\n
    :param key: function-like definition fopr which key to sort the list.\n
    :param kwargs: kwargs passed to sorted(...)., e.g. "reverse = True"\n

    :return: permutation order, sorted list
    """

    indexed = list(enumerate(lst))
    sorted_indexed = sorted(indexed, key = key, **kwargs)
    perm = [i for i, _ in sorted_indexed]
    sorted_lst = [val for _, val in sorted_indexed]
    
    return perm, sorted_lst


def nnearest(A:np.array, B:np.array, distance:float=5, n:int=2)->np.array:
    
    """
    Function that gives the indices of nearest coordinate pair in **B** to each coordinate pair in **A**.

    :param A: numpy.array containing x- and y- coordinates of points for which to search for nearest points in **B**
    :param B: numpy.array containgin x- and y- coordinates of points in which to search for the nearest points.
    :param distance: float denoting the search distance around each coordinate pair in **A**
    :param n: integer denoting the maximum number of nearest coordinate pairs in **B** within **distance** for each pair in **A**

    :return idx: numpy.array with indices of nearest matches in **B** to each pair in **A** at corresponding array index of **idx**
    """
    
    
    from scipy.spatial import cKDTree

    btree = cKDTree(B)
    dist, idx = btree.query(A, k=n, distance_upper_bound=distance)

    return (idx, dist)

def split_lines_at_points(
        line, 
        points
    ):

    """
    Split LineString at given MultiPoint objects
    Information from line DataFrame are currently dropped and only geometry is preserved!

    :param line: shapely.LineString object
    :param points: shapely.MultiPoint object
    :return: shapely.LineString object of split line object
    """

    # Internal function definitions
    def split_line_by_point(line, point, tolerance: float = 1.0e-1):

        # Imports
        from shapely.ops import split, snap
        return split(snap(line, point, tolerance), point)

    result = split_line_by_point(line, points)
    return(result)

def extend_line(
        line:LineString, 
        offset:float
        )->LineString:
    
    """
    Function that extends a shapely LineString object by a fixed distance *offset* in both directions at either end.\n
    
    :param line: shapely LineString object.\n
    :param offset: float denoting the fixed distance by which to extend line.\n
    :return: Extended line as shapely LineString object.
    """

    coords = list(line.coords)
    # Get the direction vector at the start of the LineString
    start_vec = np.array(coords[1]) - np.array(coords[0])
    start_dir = start_vec / np.linalg.norm(start_vec)
    # Calculate the new start point by moving backwards along the direction
    new_start = np.array(coords[0]) - offset * start_dir

    # Get the direction vector at the end of the LineString
    end_vec = np.array(coords[-1]) - np.array(coords[-2])
    end_dir = end_vec / np.linalg.norm(end_vec)
    # Calculate the new end point by moving forwards along the direction
    new_end = np.array(coords[-1]) + offset * end_dir

    # Build a new coordinate list with the extended endpoints
    new_coords = [tuple(new_start)] + coords[1:-1] + [tuple(new_end)]
    
    return LineString(new_coords)

def split_lines(
        lines:gp.GeoDataFrame, # DataFrame of shapely.LineString objects
        points:gp.GeoDataFrame # DataFrame of shapely.Point objects
        ):

    """
    Function that splits all line objects in **lines** at the provided point objects in **points**

    :param lines: GeoDataFrame containing line objects
    :param points: GeoDataFrame containing point objects at which to split lines

    :return: GeoDataFrame of split ine objects
    """  
    
    
    # Imports
    import pandas as pd   


    # Initialize output data
    lines['original_index'] = lines.index

    lines_out = lines.head(0).copy()
    objs = []

    for row, pip in lines.iterrows():

        res = split_lines_at_points(pip['geometry'], points)
        objs += [n for n in res.geoms]

        nr_of_objs = len(res.geoms)

        lines_out = pd.concat((lines_out, pd.concat([pd.DataFrame(pip).T] * nr_of_objs, ignore_index = True)), ignore_index = True)

    lines_out['geometry'] = objs

    # Make sure to maintain crs
    lines_out.set_crs(lines.crs, inplace = True)

    return (lines_out)


def split_lines_at_length(
        lines:gp.GeoDataFrame,
        distance:float,
        min_distance_last_segment:float,
        geom_col:str = 'geometry',
        keep_cols:bool = True,
        keep_original_line_idx:bool = False,
        return_splittingPoints:bool = False
        ) -> gp.GeoDataFrame:
    
    '''
    Function to split shapely.LineString objects in geopandas GeoDataFrame at given lengths along their path.\n

    :param lines: geopandas.GeoDataFrame with shapely.LineString objects.\n
    :param distance: float defining the desired equidistant lengths at which lkines shall be split.\n
    :param min_distance: flaot denoting the minimum distance the last segment of splitted lines shall maintain.\n
    :param geom_col: str denoting the nemae of the geometry columns. Defaults to "geometry"\n
    :param keep_cols: Boolean iof all attributes from GeoDataFrame lines shall be kept.\n
    :param keep_original_line_idx: Boolean if original row index of lines shall be transferred to output splittingPoints.\n
    :param return_splittingPointS: Bool if GeoDataFrame with MultiPoint objects shall be returned as well.\n
    :return: GeoDataFrame of splitted line objects (and optionally GeoDataFrame of splitting points)
    '''

    import shapely

    cs = lines.crs
    
    # Create temporary Multipoint object at which to split lines
    lines['splitter']                   = lines[geom_col].apply(lambda x: create_point_splitter_npoints(x, distance = distance, lmin = min_distance_last_segment))

    # Initialize output data
    lines['original_index'] = lines.index

    lines_split = lines.head(0).copy()
    objs = []

    for row, pip in lines.iterrows():

        if not isinstance(pip['splitter'], (shapely.geometry.multipoint.MultiPoint, shapely.geometry.point.Point)):
            splitpoints = None
        else:
            splitpoints = pip['splitter']

        try:
            res = split_lines_at_points(pip['geometry'], splitpoints) if splitpoints is not None else pip['geometry']
        
        except:
            print('Stop')

        objs += [n for n in res.geoms]

        nr_of_objs = len(res.geoms)

        lines_split = pd.concat((lines_split, pd.concat([pd.DataFrame(pip).T] * nr_of_objs, ignore_index = True)), ignore_index = True)

    lines_split['geometry'] = objs
    lines_split.set_crs(cs, inplace = True)    

    if return_splittingPoints:
        splittingPoints = gp.GeoDataFrame(data = {'original_line_idx':lines.index} if keep_original_line_idx else None, geometry = lines['splitter'].to_list())
        #lines_split['splitter'].copy().rename({'splitter':'geometry'})
        lines_split.drop(columns = ['splitter'], inplace = True)

        return lines_split, gp.GeoDataFrame(splittingPoints, geometry = 'geometry', crs = lines.crs)

    else:
        lines_split.drop(columns = ['splitter'], inplace = True)
        return lines_split

def extractRasterValsAtPoints(
        j:gp.GeoDataFrame, 
        raster:rio.DatasetReader,
        statsType:str = 'min', 
        buff_distance:float = 2,
        attrCol:str = 'rasterVal')->gp.GeoDataFrame:

    """
    Function that extracts heights at specified points **j** from digital surface model **dgm**.
    Within a radius of **buff_distance** around each junction the **statsType** value (allowed arguments: ['min', 'max', 'mean', 'median', 'majority']) is searched in **dgm** and assigned to the junction.
    
    :param j: GeoDataFrame of shapely.Point objects
    :param raster: raster dataset
    :param statsType: str denoting the desired type of zonal statistics applied to raster stats at each point
    :param buff_distance: float denoting the desired buffer around each point in **j** from which to extract data
    :param attrCol: str denoting the desired attribute name in **j** in which to store raster data values
    :return: GeoDataFrame of points **j**
    """ 


    # Import
    from rasterstats import zonal_stats

    # Convert to array
    arr = raster.read(1)

    # Create affine transformation
    affine = raster.transform

    # Read resolution from affine transformation
    (x, y) = (affine[0], -affine[4])

    if buff_distance <= min(x,y):
        buff_distance_edited = 1.5 * min(x,y)
        print('\n### Buffer distance is extended from %f.2 to %f.2 because min. raster side length is %f.2 ###' % (buff_distance, buff_distance_edited, min(x,y)))
    else:
        buff_distance_edited = buff_distance

    # Data preparation
    j_buff = j.copy()
    j_buff['geometry'] = j_buff.buffer(distance = buff_distance_edited)

    # Find minimum value (e.g. outside house ground surface) within buffer distance around point, ignoring nodata/missing data (np.nan in input array)
    zs = zonal_stats(j_buff, arr, affine = affine, stats=statsType, nodata = -999)
    j[attrCol] = [c[statsType] for c in zs]

    return j

def match_polygons_to_points_by_intersection(
    polygons:gp.GeoDataFrame,
    points:gp.GeoDataFrame,
    attributes:list = [])->gp.GeoDataFrame:

    '''
    Function that creates matching of polygons to points by intersection. If multiple points within the provided entity of points intersects with a polygon, the first match is taken.
    Search for intersection is based on geopdnas operation.

    :param polygons: GeoDataFrame of polygons (shapely.Polygon)
    :param points: GeoDataFrame of points (shapely.Point)
    :param attributes: ist of attributes within **points** that shall be transferred to the matching output GeoDataFrame of polygons.
    
    :return: GeoDataFrame **polygons** with added attribute for connected junction
    '''

    # Plausbility checks       
    if 'index_right' in polygons.columns:
        polygons.drop(columns = 'index_right', inplace = True)
        print('\n### Attention! Existing column "index_right" found in polygons DataFrame. It is overwritten. ###')

    if 'geometry' not in attributes:
        attributes.append('geometry')

    # Code
    points = points.copy()
    polygons = polygons.copy()

    polygons.reset_index(drop = True, inplace = True)
    polygons['indexLeftTemp'] = polygons.index

    poly = gp.sjoin(left_df = polygons, right_df = points[attributes], predicate = 'intersects', how = 'left')

    attrs = attributes.copy()
    attrs.remove('geometry')
    attrs.append('index_right')

    s = poly.groupby(by = 'indexLeftTemp').nth(0)[attrs]

    for attr in attrs:
        polygons[attr] = s[attr]

    polygons.drop(columns = 'indexLeftTemp', inplace = True)

    print('\n### Function "match_polygons_to_points_by_intersection": %.0f of %.0f polygons are found with matching points. ###' % (len(polygons[~polygons['index_right'].isnull()]), len(polygons)))

    # Find polygons with multiple matches for points
    multipleMatches = poly.groupby(by = 'indexLeftTemp').transform('count')
    multipleMatches.rename(columns = {'index_right':'matches'}, inplace = True)

    if len(multipleMatches[multipleMatches['matches'] > 1]) > 0:
        print('\n### Polygons with multiple results for matching points by intersection are: ', multipleMatches[multipleMatches['matches'] > 1],' ###')

    return polygons

def geodata_from_geometry(df:gp.GeoDataFrame, typ:str):

    """
    Function that converts geometry column of shapely objects into pandapipes-compatible geodata column.

    :param df: GeoDataFrame containing geometry-column
    :param typ: str denoting the type of geometry object
    :return: GeoDataFrame with edited geodata column
    """
    
    if (typ == 'line'):
        df['geodata'] = df['geometry'].astype(str)
        df['geodata'] = df['geodata'].str.replace('LINESTRING', '')
        df['geodata'] = df['geodata'].str.lstrip()
        df['geodata'] = df['geodata'].str.replace(', ', '),(')
        df['geodata'] = df['geodata'].str.replace(' ', ',')
        df['geodata'] = '[' + df['geodata'] + ']'

    elif (typ == 'point'):
        df['geodata'] = df['geometry'].astype(str)
        df['geodata'] = df['geodata'].str.replace('POINT', '')
        df['geodata'] = df['geodata'].str.lstrip()
        df['geodata'] = df['geodata'].str.replace(' ', ', ')

    elif (typ == "polygon"):
        coords_list = []
        for poly in df["geometry"]:            
            coords = [[x, y] for x, y in poly.exterior.coords]
            coords_list.append(coords)
        df["geodata"] = coords_list

    def lit(x):
        import ast
        try:
            return (ast.literal_eval(x))
        except:
            return (x)

    df['temp'] = df['geodata'].apply(lambda x: lit(x))
    df['geodata'] = df['temp']
    df.drop(columns='temp', inplace=True)

    return (df)

def get_closest_value(
        input_value:float, 
        input_list:list
        )->float:
    
    """
    Function that returns the closest matching value from an input list to a given value

    :param input_value: float for which closest matching value in **input_list** shall be searched.
    :param input_list: list-like in which to search for closoest value to **input_value**
    :return: closest value from **input_list**
    """

    arr = np.asarray(input_list)

    i = (np.abs(arr - input_value)).argmin()

    return arr[i]

def get_next_higher_value(
        input_value:float, 
        input_list:list
    )->float:
    """
    Function that returns the next higher value from an input list compared to a given value.
    If no higher value exists, returns the maximum value from the list.

    :param input_value: float for which next higher value in **input_list** shall be searched.
    :param input_list: list-like in which to search for next higher value to **input_value**
    :return: next higher value from **input_list** or maximum value if none exists
    """

    arr = np.asarray(input_list)
    higher_values = arr[arr > input_value]
    
    if len(higher_values) > 0:
        return higher_values.min()
    else:
        return arr.max()

def cycleEdges(
        G:nx.graph,
        cycleNodes:list, 
        weights:str
        )->tuple:
    
    """
    Function that takes a networkx graph G and a list of a list of nodes for each cycle found in G.

    :param G: networkx graph
    :param cycleNodes: list-like with lists of nodes for each cycle found in G, e.g. [[0, 1, 2], [2, 4, 3]]
    :param weights: string denoting the edge attributes which shall be summed up for each cycle
    :return: list of edges within each cycle, list of summed weights within each cycle
    """

    # Initialisations
    cycleEdges_list = []
    cycleEdgesWeights = []

    for nodes in cycleNodes:
        summedWeight = 0
        edgesInCycle = []
        for u,v,a in G.edges(data = True):
            if u in nodes and v in nodes:
                summedWeight += a[weights]
                edgesInCycle.append((u,v))
                
        cycleEdges_list.append(edgesInCycle)
        cycleEdgesWeights.append(summedWeight)

    return cycleEdges_list, cycleEdgesWeights

def checkConnectivity(lines:gp.GeoDataFrame, source:str = 'from_junction', target:str = 'to_junction'):
    
    """
    Function that takes line object network **lines** as input argument and checks if unconnected lines are present

    :param lines: GeoDataFrame of line objects.\n
    :returns: Boolean if all lines are connected.\n
    """

    # Imports
    import networkx as nx
    from shapely import set_precision

    # Create temporary copies and set precision
    lines = lines.copy()
    lines['geometry'] = set_precision(lines['geometry'], 0.01)

    G = gdf_to_nx(lines, approach="primal")
    # G = nx.from_pandas_edgelist(lines, source = source, target = target, edge_attr = True)

    # Create subgraphs of connected components in G
    if not nx.is_connected(G):
        subs = [G.subgraph(c).copy() for c in nx.connected_components(G)]
        isConnected = False
    else:
        subs = [G]
        isConnected = True

    # Convert graphs back to GeoDataFrames
    gdf_out = gp.GeoDataFrame()
    for n, sub in enumerate(subs):
        edges = nx_to_gdf(sub, lines=True, points = False)
        edges = edges[[col for col in edges.columns if col in lines.columns]]

        edges['subgraph'] = n

        gdf_out = pd.concat((gdf_out, edges), axis = 0)

    return isConnected, gdf_out


def closest_point(points:np.array, target:tuple, threshDistance:float = np.inf):
    """
    Find the point in the array that is closest to the target point.

    :param points: A 2D array where each row represents a point (x, y).
    :param target: A tuple representing the target point (x, y).
    :param threshDistance: Distance threshold, if target is farther away than threshold, no matching closest point is returns (e.g. empty tuple)

    :returns:
    np.ndarray: The point from the array that is closest to the target.
    """
    # Calculate the Euclidean distances from the target to each point
    distances = np.linalg.norm(points - target, axis=1)
    
    # Find the index of the minimum distance
    closest_index = np.argmin(distances)
    
    # Return the closest point
    if distances[closest_index] <= threshDistance:
        return points[closest_index]
    else: # Return empty tuple
        return tuple()

def get_dict_from_aggregated_groups(
        df:pd.DataFrame,
        groupCol:str,
        val:str,
        func:str = 'max'
    ):

    """
    Function that returns dictionaries from grouped attributes in pandas DataFrame for following matchings:
    - dictionary with row indices as keys and aggregated value in column **val** by **func** (e.g. 'max', 'mean') from affiliated **groupCol** as value
    - dictionary with aggregated value in column **val** from affiliated **groupCol** as keys and list of row indices as values

    :param df: pandas.DataFrame
    :param groupCol: str denoting column in which to group
    :param val: str denoting column from which to take values
    :param func: function which shall be applied to the grouped values (implemented: 'max', 'mean', 'min')
    :returns: two dictionaries
    """

    if func not in ('max', 'mean', 'min'):
        print('\n### Please choose proper function for aggregation of grouped values. ###')
        return

    df['newCol'] = df.groupby(by = groupCol)[val].transform(func)
    idx = df.index.to_series()
    
    dict1 = df['newCol'].to_dict()
    dict2 = idx.groupby(df['newCol']).agg(list).to_dict()

    return dict1, dict2

def drop_elements(net, element_type:str, element_index):
    """
    Drops element, result and group entries from the pandapipes net.

    See Also
    --------
    drop_elements : providing more generic usage (inter-element connections considered)
    """

    ### Imports
    import pandapipes as ppi

    ### Definitions
    def ensure_iterability(var, len_=None):
        """
        Ensures iterability of a variable (and also the length if given).

        Examples
        --------
        >>> ensure_iterability([1, 2])
        [1, 2]
        >>> ensure_iterability(1)
        [1]
        >>> ensure_iterability("Hi")
        ["Hi"]
        >>> ensure_iterability([1, 2], len_=2)
        [1, 2]
        >>> ensure_iterability([1, 2], len_=3)
        ValueError("Length of variable differs from 3.")
        """
        if hasattr(var, "__iter__") and not isinstance(var, str):
            if isinstance(len_, int) and len(var) != len_:
                raise ValueError("Length of variable differs from %i." % len_)
        else:
            len_ = len_ or 1
            var = [var] * len_
        return var

    element_index = ensure_iterability(element_index)
    net[element_type] = net[element_type].drop(element_index)

    # res_element
    res_element_type = "res_" + element_type
    if res_element_type in net.keys() and isinstance(net[res_element_type], pd.DataFrame):
        drop_res_idx = net[res_element_type].index.intersection(element_index)
        net[res_element_type] = net[res_element_type].drop(drop_res_idx)

def translate_element_geometries(
    net,
    element_type:str = 'junction',
    elementIdx:list = None,
    trans:list = [0,0]
    ):


    """
    Function that applies a 2D translation to elements **junction** or **pipe** of pandapipes network.\n
    If available, a geometry column in the component Dataframe is updated.
    In any case, the DataFrame in net.element_type_geodata is tried to be updated.\n
    
    :param net: pandapipes network object\n
    :param element_type: either "junction" or "pipe"\n
    :param elementIdx: iterable for speecific element row indices at which geometries shall be updated. If None, all elements will be translated\n
    :param trans: list-like of values for translation in x- and y-direction [trans_x, trans_y]

    :returns: Inplace updates
    """

    # Plausibility checks
    if not hasattr(net, element_type):
        print(f'\n### Specified network element type {element_type} not included in network. Aborting... ###')
    else:

        if elementIdx is None:
            elementIdx = net[element_type].index

        # Update geometry
        if 'geometry' not in net[element_type]:
            print(f'\n### Column "geometry" with shapely geometry information not included in network element {element_type}. ###')
            updatedGeometry = False

        else:
            net[element_type].loc[elementIdx, 'geometry'] = gp.GeoSeries(net[element_type].loc[elementIdx, 'geometry']).translate(trans[0], trans[1], 0)
            print(f'\n### shapely.geometry information in network component {element_type} updated/trasnlated. ###\n')

        # Update geodata
        if hasattr(net, element_type+'_geodata'):

            if element_type == 'junction':                
                net[element_type+'_geodata'].loc[elementIdx, 'x'] += trans[0]
                net[element_type+'_geodata'].loc[elementIdx, 'y'] += trans[1]

            if element_type == 'pipe':
                net[element_type+'_geodata']['coords'] = net[element_type+'_geodata'].apply(lambda x: [tuple(np.array(x['coords'][k]) + np.array(trans)) for k in range(len(x['coords']))] if x.name in elementIdx else x['coords'], axis = 1)
            
            print(f'\n### Columns in {element_type}_geodata with geometry information for pandapipes element {element_type} are updated. ###')

### *--- Functions for pressure loss calculations ---*
def func_nikuradse(
        d:float, 
        mdot:float=25, 
        dp_per_l:float=100, 
        rho:float=1000, 
        nu:float=0.413e-6, 
        k:float=0.0469
        ):
    
    """
    Implicit formulation of specific pressure loss according to Nikuradse formulae.

    Params all in SI units if not specified otherwise...

    :param d: float denoting diameter of the pipe.\n
    :param mdot: float denoting the mass flow rate in the pipe.\n
    :param dp_per_l: float denoting specific pressure loss as set value.\n
    :param rho: density of the medium inside the pipe.\n
    :param nu: kinematic viscosity of the medium inside the pipe.\n
    :param k: float denoting the wall roughness at the inner pipe wall (mm).\n

    :return: float difference denoting the differnce between set value of specific pressure loss and calculated value by other arguments.\n
    """

    difference = dp_per_l - (8 * mdot ** 2 / (rho * d**5 * np.pi**2) * (
                64 * rho * nu * d * np.pi / (4 * mdot) + 1 / (-2 * np.log10(k/1e3 / (3.71 * d)))** 2))
    return difference

def func_swameejain(
    d:float, 
    mdot:float = 25, 
    dp_per_l:float = 100, 
    rho:float = 1000, 
    nu:float = 0.413e-06, 
    k:float = 0.0469
    ):

    """
    Implicit formulation of specific pressure loss according to Swamee-Jain formulae.

    Params all in SI units if not specified otherwise...

    :param d: float denoting diameter of the pipe.\n
    :param mdot: float denoting the mass flow rate in the pipe.\n
    :param dp_per_l: float denoting specific pressure loss as set value.\n
    :param rho: density of the medium inside the pipe.\n
    :param nu: kinematic viscosity of the medium inside the pipe.\n
    :param k: float denoting the wall roughness at the inner pipe wall (mm).\n

    :return: float difference denoting the differnce between set value of specific pressure loss and calculated value by other arguments.\n
    """

    difference = dp_per_l - rho / (2*d) * (mdot * 4 / (rho*d**2*np.pi))**2 * 0.25 / (np.log10(k/1e3/(3.71*d) + (5.74/(4*mdot/(d*np.pi*nu*rho))**(0.9))))**2

    return difference

### Function adopted from momepy module

def _primal_to_gdf(net, points, lines, spatial_weights, node_id):
    """Generate gdf(s) from a primal network. Helper for ``nx_to_gdf``."""
    """Function adopted from momepy: Fleischmann, M. (2019) ‘momepy: Urban Morphology Measuring Toolkit’, Journal of Open Source Software, 4(43), p. 1807. doi: 10.21105/joss.01807."""
    if points is True:
        gdf_nodes = _points_to_gdf(net)

        if spatial_weights is True:
            weights = libpysal.weights.W.from_networkx(net)
            weights.transform = "b"

    if lines is True:
        gdf_edges = _lines_to_gdf(net, points, node_id)

    if points is True and lines is True:
        if spatial_weights is True:
            return gdf_nodes, gdf_edges, weights
        return gdf_nodes, gdf_edges
    if points is True and lines is False:
        if spatial_weights is True:
            return gdf_nodes, weights
        return gdf_nodes
    return gdf_edges

def _points_to_gdf(net):
    """Generate a point gdf from nodes. Helper for ``nx_to_gdf``."""
    """Function adopted from momepy: Fleischmann, M. (2019) ‘momepy: Urban Morphology Measuring Toolkit’, Journal of Open Source Software, 4(43), p. 1807. doi: 10.21105/joss.01807."""
    node_xy, node_data = zip(*net.nodes(data=True), strict=True)
    if isinstance(node_xy[0], int) and "x" in node_data[0]:
        geometry = [Point(data["x"], data["y"]) for data in node_data]  # osmnx graph
    else:
        geometry = [Point(*p) for p in node_xy]
    gdf_nodes = gp.GeoDataFrame(list(node_data), geometry=geometry)
    if "crs" in net.graph:
        gdf_nodes.crs = net.graph["crs"]
    return gdf_nodes

def _lines_to_gdf(net, points, node_id):
    """Generate a linestring gdf from edges. Helper for ``nx_to_gdf``."""
    """Function adopted from momepy: Fleischmann, M. (2019) ‘momepy: Urban Morphology Measuring Toolkit’, Journal of Open Source Software, 4(43), p. 1807. doi: 10.21105/joss.01807."""
    starts, ends, edge_data = zip(*net.edges(data=True), strict=True)
    gdf_edges = gp.GeoDataFrame(list(edge_data))

    if points is True:
        gdf_edges["node_start"] = [net.nodes[s][node_id] for s in starts]
        gdf_edges["node_end"] = [net.nodes[e][node_id] for e in ends]

    if "crs" in net.graph:
        gdf_edges.crs = net.graph["crs"]

    return gdf_edges

def gdf_to_nx(
    gdf_network,
    approach="primal",
    length="mm_len",
    multigraph=True,
    directed=False,
    angles=True,
    angle="angle",
    oneway_column=None,
):
    """
    Convert a LineString GeoDataFrame to a ``networkx.MultiGraph`` or other
    Graph as per specification. Columns are preserved  as edge or node
    attributes (depending on the ``approach``). Index is not preserved.

    See the User Guide page :doc:`../../user_guide/graph/convert` for details.

    Parameters
    ----------
    gdf_network : GeoDataFrame
        A GeoDataFrame containing objects to convert.
    approach : str, default 'primal'
        Allowed options are ``'primal'`` or ``'dual'``. Primal graphs represent
        endpoints as nodes and LineStrings as edges. Dual graphs represent
        LineStrings as nodes and their topological relation as edges. In such a
        case, it can encode an angle between LineStrings as an edge attribute.
    length : str, default 'mm_len'
        The attribute name of segment length (geographical)
        which will be saved to the graph.
    multigraph : bool, default True
        Create a ``MultiGraph`` of ``Graph`` (potentially directed).
        ``MutliGraph`` allows multiple edges between any pair of nodes,
        which is a common case in street networks.
    directed : bool, default False
        Create a directed graph (``DiGraph`` or ``MultiDiGraph``).
        Directionality follows the order of LineString coordinates.
    angles : bool, default True
        Capture the angles between LineStrings as an attribute of a dual graph.
        Ignored if ``approach='primal'``.
    angle : str, default 'angle'
        The attribute name of the angle between LineStrings which will
        be saved to the graph. Ignored if ``approach='primal'``.
    oneway_column : str, default None
        Create an additional edge for each LineString which allows bidirectional
        path traversal by specifying the boolean column in the GeoDataFrame. Note,
        that the reverse conversion ``nx_to_gdf(gdf_to_nx(gdf, directed=True,
        oneway_column="oneway"))`` will contain additional duplicated geometries.

    Returns
    -------
    net : networkx.Graph, networkx.MultiGraph, networkx.DiGraph, networkx.MultiDiGraph
        Graph as per specification.

    See also
    --------
    nx_to_gdf

    Examples
    --------
    >>> import geopandas as gp
    >>> df.head(5)
                                                geometry
    0  LINESTRING (1603585.640 6464428.774, 1603413.2...
    1  LINESTRING (1603268.502 6464060.781, 1603296.8...
    2  LINESTRING (1603607.303 6464181.853, 1603592.8...
    3  LINESTRING (1603678.970 6464477.215, 1603675.6...
    4  LINESTRING (1603537.194 6464558.112, 1603557.6...

    Primal graph:

    >>> G = gdf_to_nx(df)
    >>> G
    <networkx.classes.multigraph.MultiGraph object at 0x7f8cf90fad50>

    >>> G_directed = gdf_to_nx(df, directed=True)
    >>> G_directed
    <networkx.classes.multidigraph.MultiDiGraph object at 0x7f8cf90f56d0>

    >>> G_digraph = gdf_to_nx(df, multigraph=False, directed=True)
    >>> G_digraph
    <networkx.classes.digraph.DiGraph object at 0x7f8cf9150c10>

    >>> G_graph = gdf_to_nx(df, multigraph=False, directed=False)
    >>> G_graph
    <networkx.classes.graph.Graph object at 0x7f8cf90facd0>

    Dual graph:

    >>> G_dual = gdf_to_nx(df, approach="dual")
    >>> G_dual
    <networkx.classes.multigraph.MultiGraph object at 0x7f8cf9150fd0>

    """
    """Function adopted from momepy: Fleischmann, M. (2019) ‘momepy: Urban Morphology Measuring Toolkit’, Journal of Open Source Software, 4(43), p. 1807. doi: 10.21105/joss.01807."""
    gdf_network = gdf_network.copy()
    if "key" in gdf_network.columns:
        gdf_network.rename(columns={"key": "__key"}, inplace=True)

    if multigraph and directed:
        net = nx.MultiDiGraph()
    elif multigraph and not directed:
        net = nx.MultiGraph()
    elif not multigraph and directed:
        net = nx.DiGraph()
    else:
        net = nx.Graph()

    net.graph["crs"] = gdf_network.crs
    gdf_network[length] = gdf_network.geometry.length
    fields = list(gdf_network.columns)

    if approach == "primal":
        if oneway_column and not directed:
            raise ValueError(
                "Bidirectional lines are only supported for directed graphs."
            )

        _generate_primal(net, gdf_network, fields, multigraph, oneway_column)

    elif approach == "dual":
        if directed:
            raise ValueError("Directed graphs are not supported in dual approach.")

        _generate_dual(
            net, gdf_network, fields, angles=angles, multigraph=multigraph, angle=angle
        )

    else:
        raise ValueError(
            f"Approach '{approach}' is not supported. Use 'primal' or 'dual'."
        )

    return net

def nx_to_gdf(
    net, points=True, lines=True, spatial_weights=False, nodeID="nodeID"  # noqa
):
    """
    Convert a ``networkx.Graph`` to a LineString GeoDataFrame and Point GeoDataFrame.

    Automatically detects an ``approach`` of the graph and assigns
    edges and nodes to relevant geometry type.

    See the User Guide page :doc:`../../user_guide/graph/convert` for details.

    Parameters
    ----------
    net : networkx.Graph
        A ``networkx.Graph`` object.
    points : bool (default is ``True``)
        Export point-based gdf representing intersections.
    lines : bool (default is ``True``)
        Export line-based gdf representing streets.
    spatial_weights : bool (default is ``False``)
        Set to ``True`` to export a libpysal spatial weights
        for nodes (only for primal graphs).
    nodeID : str
        The name of the node ID column to be generated.

    Returns
    -------
    GeoDataFrame
       The  Selected gdf or tuple of both gdfs or tuple of gdfs and weights.

    See also
    --------
    gdf_to_nx

    Examples
    --------
    >>> import geopandas as gpd
    >>> df = gp.read_file(datasets.get_path('bubenec'), layer='streets')
    >>> df.head(2)
                                                geometry
    0  LINESTRING (1603585.640 6464428.774, 1603413.2...
    1  LINESTRING (1603268.502 6464060.781, 1603296.8...
    >>> G = gdf_to_nx(df)

    Converting the primal Graph to points as intersections and lines as street segments:

    >>> points, lines = nx_to_gdf(graph)
    >>> points.head(2)
       nodeID                         geometry
    0       1  POINT (1603585.640 6464428.774)
    1       2  POINT (1603413.206 6464228.730)
    >>> lines.head(2)
                         geometry      mm_len  node_start  node_end
    0  LINESTRING (1603585.640...  264.103950           1         2
    1  LINESTRING (1603561.740...   70.020202           1         9

    Storing the relationship between points/nodes as a libpysal W object:

    >>> points, lines, W = nx_to_gdf(graph, spatial_weights=True)
    >>> W
    <libpysal.weights.weights.W object at 0x7f8d01837210>

    Converting the dual Graph to lines. The dual Graph does not export edges to GDF:

    >>> G = gdf_to_nx(df, approach="dual")
    >>> lines = nx_to_gdf(graph)
    >>> lines.head(2)
                                                geometry      mm_len
    0  LINESTRING (1603585.640 6464428.774, 1603413.2...  264.103950
    1  LINESTRING (1603607.303 6464181.853, 1603592.8...  199.746503
    """
    """Function adopted from momepy: Fleischmann, M. (2019) ‘momepy: Urban Morphology Measuring Toolkit’, Journal of Open Source Software, 4(43), p. 1807. doi: 10.21105/joss.01807."""
    # generate nodes and edges geodataframes from graph
    primal = None
    if "approach" in net.graph:
        if net.graph["approach"] == "primal":
            primal = True
        elif net.graph["approach"] == "dual":
            return _dual_to_gdf(net)
        else:
            raise ValueError(
                f"Approach '{net.graph['approach']}' is not supported. "
                "Use 'primal' or 'dual'."
            )

    if not primal:
        warnings.warn(
            message="Approach is not set. Defaulting to 'primal'.",
            category=UserWarning,
            stacklevel=2,
        )

    for nid, n in enumerate(net):
        net.nodes[n][nodeID] = nid

    return _primal_to_gdf(
        net,
        points=points,
        lines=lines,
        spatial_weights=spatial_weights,
        node_id=nodeID,
    )


def _generate_primal(graph, gdf_network, fields, multigraph, oneway_column=None):
    """Generate a primal graph. Helper for ``gdf_to_nx``."""
    """Function adopted from momepy: Fleischmann, M. (2019) ‘momepy: Urban Morphology Measuring Toolkit’, Journal of Open Source Software, 4(43), p. 1807. doi: 10.21105/joss.01807."""
    graph.graph["approach"] = "primal"

    msg = (
        " This can lead to unexpected behaviour. "
        "The intended usage of the conversion function "
        "is with networks made of LineStrings only."
    )

    if "LineString" not in gdf_network.geom_type.unique():
        warnings.warn(
            message="The given network does not contain any LineString." + msg,
            category=RuntimeWarning,
            stacklevel=3,
        )

    if len(gdf_network.geom_type.unique()) > 1:
        warnings.warn(
            message="The given network consists of multiple geometry types." + msg,
            category=RuntimeWarning,
            stacklevel=3,
        )

    key = 0
    for row in gdf_network.itertuples():
        first = row.geometry.coords[0]
        last = row.geometry.coords[-1]

        data = list(row)[1:]
        attributes = dict(zip(fields, data, strict=True))
        if multigraph:
            graph.add_edge(first, last, key=key, **attributes)
            key += 1

            if oneway_column:
                oneway = bool(getattr(row, oneway_column))
                if not oneway:
                    graph.add_edge(last, first, key=key, **attributes)
                    key += 1
        else:
            graph.add_edge(first, last, **attributes)

def _generate_dual(graph, gdf_network, fields, angles, multigraph, angle):
    """Generate a dual graph. Helper for ``gdf_to_nx``."""
    """Function adopted from momepy: Fleischmann, M. (2019) ‘momepy: Urban Morphology Measuring Toolkit’, Journal of Open Source Software, 4(43), p. 1807. doi: 10.21105/joss.01807."""
    graph.graph["approach"] = "dual"
    key = 0

    sw = libpysal.weights.Queen.from_dataframe(gdf_network, silence_warnings=True)
    cent = gdf_network.geometry.centroid
    gdf_network["temp_x_coords"] = cent.x
    gdf_network["temp_y_coords"] = cent.y

    for i, row in enumerate(gdf_network.itertuples()):
        centroid = (row.temp_x_coords, row.temp_y_coords)
        data = list(row)[1:-2]
        attributes = dict(zip(fields, data, strict=True))
        graph.add_node(centroid, **attributes)

        if sw.cardinalities[i] > 0:
            for n in sw.neighbors[i]:
                start = centroid
                end = (
                    gdf_network["temp_x_coords"].iloc[n],
                    gdf_network["temp_y_coords"].iloc[n],
                )
                p0 = row.geometry.coords[0]
                p1 = row.geometry.coords[-1]
                geom = gdf_network.geometry.iloc[n]
                p2 = geom.coords[0]
                p3 = geom.coords[-1]
                points = [p0, p1, p2, p3]
                shared = [x for x in points if points.count(x) > 1]
                if shared:  # fix for non-planar graph
                    remaining = [e for e in points if e not in [shared[0]]]
                    if len(remaining) == 2:
                        if angles:
                            angle_value = _angle(remaining[0], shared[0], remaining[1])
                            if multigraph:
                                graph.add_edge(
                                    start, end, key=0, **{angle: angle_value}
                                )
                                key += 1
                            else:
                                graph.add_edge(start, end, **{angle: angle_value})
                        else:
                            if multigraph:
                                graph.add_edge(start, end, key=0)
                                key += 1
                            else:
                                graph.add_edge(start, end)

def _angle(a, b, c):
    """
    Measure the angle between a-b, b-c (in degrees). Helper for ``gdf_to_nx``.
    Adapted from cityseer's implementation.
    """
    """Function adopted from momepy: Fleischmann, M. (2019) ‘momepy: Urban Morphology Measuring Toolkit’, Journal of Open Source Software, 4(43), p. 1807. doi: 10.21105/joss.01807."""
    a1 = math.degrees(math.atan2(b[1] - a[1], b[0] - a[0]))
    a2 = math.degrees(math.atan2(c[1] - b[1], c[0] - b[0]))
    return abs((a2 - a1 + 180) % 360 - 180)

def _dual_to_gdf(net):
    """Generate a linestring gdf from a dual network. Helper for ``nx_to_gdf``."""
    """Function adopted from momepy: Fleischmann, M. (2019) ‘momepy: Urban Morphology Measuring Toolkit’, Journal of Open Source Software, 4(43), p. 1807. doi: 10.21105/joss.01807."""
    starts, edge_data = zip(*net.nodes(data=True), strict=True)
    gdf_edges = gp.GeoDataFrame(list(edge_data))
    gdf_edges.crs = net.graph["crs"]
    return gdf_edges
