
# -*- coding: utf-8 -*-

### Imports

import geopandas as gp
import pandas as pd
import pandapipes.toolbox as ppitlbx
import numpy as np
from shapely import (set_precision, STRtree)
from shapely.ops import (unary_union,)
from shapely.geometry import (LineString, Point, MultiPoint)
from collections import Counter

from geoppi.auxFunctions import (extractPointsFromLines, nnearest, split_lines, split_lines_at_points, create_circumferential_points, closest_objects_to_points,             match_polygons_to_points_by_intersection, geodata_from_geometry, extractRasterValsAtPoints, checkConnectivity, create_point_splitter_npoints, extend_line, assign_attr_by_max_intersection_area, round_coords_2D)

def create_connection_lines(
        lines:gp.GeoDataFrame,
        polys:gp.GeoDataFrame,
        unique_ID_polys:str = 'unID_polys',
        unique_ID_lines:str = 'unID_lines',
        dist_connection_points:float = 2.5,
        max_dist_connection_lines:float = 100
    )->gp.GeoDataFrame:

    """
    Function to automatically create connection lines from polygons to lines.\n
    Connection lines are crated along straight paths (shotrest distance) from polygons boundaries to outer lines provided as GeoDataFrame. An overlap of the start/end point of connection lines and the polygons is ensured. Lines are split at crossing junctions after creating the connection lines.\n

    :param lines: GeoDataFrame of line objects.\n
    :param polys: GeoDataFrame of polygon objects.\n
    :param unique_ID_polys: str denoting the desired attribute name for temporary unique ID of polygons.\n
    :param unique_ID_lines: str denoting temporary unique ID for line objects.\n
    :param dist_connection_points: float denoting the desired distance between possbiel connection points along the line objects. Shorter distances increase the computation time.\n
    :param max_dist_connection_lines: float denoting the maximum distance from polygons to the lines to consider when creating connection lines.\n
    :return: GeoDataFrame of all line objects (original (split at crossing junctions) + connection lines).
    """
    
    ## Initializations    
    cs = lines.crs
    if cs != polys.crs:
        polys = polys.copy().to_crs(cs)
        print(f'\nAttention! Polygons and lines do not share same crs!')

    minOffset_connectionLines = 0.001

    # Create unique IDs
    lines[unique_ID_lines] = np.arange(len(lines))
    polys[unique_ID_polys] = np.arange(len(polys))

    ## Creation of connection lines from polygons to closets points on provided lines
    # Negative buffer of polygons to ensure overlap with newly created connection lines
    polys_buff = polys.copy()
    polys_buff['geometry'] = polys_buff['geometry'].buffer(-0.1)

    # Create circumferential points around polygon boundaries
    circumferentialPoints = create_circumferential_points(polygons = polys_buff, npoints = 10, lmin = 2)
    circumferentialPoints[unique_ID_polys] = polys[unique_ID_polys].copy()
    circumferentialPoints = circumferentialPoints.explode(index_parts = True).reset_index(drop = True)

    # Create splitting points on provided lines
    allJunctionPoints = lines['geometry'].apply(lambda x: create_point_splitter_npoints(x, distance = dist_connection_points, lmin = dist_connection_points*1.5))
    allJunctionPoints = allJunctionPoints[~(allJunctionPoints.is_empty)].reset_index(drop = True)
    allJunctionPoints = allJunctionPoints.explode(index_parts = False).reset_index(drop = True)

    # Find closest object to each point (indicating index of closest line segment)
    idxs, distances = closest_objects_to_points(points = circumferentialPoints, geomObjects = allJunctionPoints, maxDist = np.ones(len(circumferentialPoints))*max_dist_connection_lines)
        
    # Assign nearest line segment to each point at polys' boundaries
    circumferentialPoints['nearestPoint']                            = idxs
    circumferentialPoints['distanceToPoint']                         = distances

    # Find index of points which feature the shortes distance to adjacent line segment within group of points for each polygon
    idx_ps_min_distances = circumferentialPoints[~circumferentialPoints['nearestPoint'].isnull()].groupby(by=unique_ID_polys)['distanceToPoint'].idxmin().values
    circumferentialPoints_ed = circumferentialPoints.loc[idx_ps_min_distances]
    circumferentialPoints_ed['nearestPoint_geometry'] = allJunctionPoints.iloc[circumferentialPoints_ed['nearestPoint'].values].geometry.values

    # Create valid lines between polys circumferential points and splitting points
    # Overlap of 0.001m is used to ensure proper splitting at connecting points
    connectionLines = gp.GeoDataFrame(geometry = [extend_line(LineString([p1, p2]), offset = minOffset_connectionLines) for p1, p2 in zip(circumferentialPoints_ed['geometry'], circumferentialPoints_ed['nearestPoint_geometry'])]).set_crs(cs)

    # Separate lines at intersections
    allLines_temp = unary_union(pd.concat((lines.geometry, connectionLines.geometry), axis = 0))
    allLines = gp.GeoDataFrame(geometry = list(allLines_temp.geoms)).set_crs(cs)

    # Filter only those connection lines which feature a single common start- or endpoint and have a short length
    shortConnectionLines = allLines[allLines.length < minOffset_connectionLines + 0.0001]

    # Extract endpoints of short lines
    endpoints_short = []
    for line in shortConnectionLines.geometry:
        coords = np.array(line.coords)
        endpoints_short.append((Point(coords[0]), Point(coords[-1])))

    start_pts_short = pd.Series([pt[0] for pt in endpoints_short], index=shortConnectionLines.index)
    end_pts_short = pd.Series([pt[1] for pt in endpoints_short], index=shortConnectionLines.index)

    # Gather all start- and endpoints and all coordinates of intermediate points in lines
    all_points = []
    for line in allLines.geometry:
        if not line.is_empty:
            coords = list(line.coords)
            all_points.extend([Point(xy) for xy in coords])

    # Create Counter object
    point_counter = Counter(all_points)

    # Frequency of occurrences for all start- and endpoints
    start_freqs = np.array([point_counter.get(pt, 0) for pt in start_pts_short])
    end_freqs = np.array([point_counter.get(pt, 0) for pt in end_pts_short])

    # Only short connection lines whose start- or endpoints occur max. 1 times in the entire amount of points in the lines shall be selected
    condition = (start_freqs <= 1) | (end_freqs <= 1)
    mask = shortConnectionLines[condition]

    # Negation of the conditions is the sum of all lines excluding the short connection lines. Possibly created lines with zero length are removed
    allLines = allLines[(~allLines.index.isin(mask.index)) & (allLines.geometry.length >= 1e-04)]

    ## Transfer all attributes from original lines to newly created lines by max. intersection (shapely.unary_union only returns geometries)
    allLines[f'{unique_ID_lines}_2'] = np.arange(len(allLines))

    lines_buff = lines.copy()
    lines_buff.geometry = lines_buff.geometry.buffer(0.05)

    allLines_buff = allLines.copy()
    allLines_buff.geometry = allLines_buff.geometry.buffer(0.05)

    allLines_buff = assign_attr_by_max_intersection_area(gp1 = allLines_buff, gp_source = lines_buff, gp1_id = f'{unique_ID_lines}_2', attr = list(lines_buff.columns.drop('geometry')))

    allLines_buff.geometry = allLines.geometry
    allLines_buff.drop(columns = [f'{unique_ID_lines}_2'], inplace = True)

    return allLines_buff

def assign_default_values_ppi(
        gdf:gp.GeoDataFrame, 
        type:str = 'junction', 
        drop_old_cols:bool = False,
        tidy_up:bool = False)->gp.GeoDataFrame:
    
    """
    Function that assigns default values to lines, junctions and valves as components for pandapipe networks.
    The function takes GeoDataFrames as inputs (with geometry column 'geometry') and returns them with the assigned default values.
    If any of the default columns for the pandapipes structure is found in the columns of the provided GeoDataFrame, the data are not overwritten!
    **drop_old_cols** drops all columns (except 'geometry') from DF before assignment of new columns!
    **tidy_up** drops all columns which are not included in the default set of columns for use in pandapipes datastructure after assignment.
    
    Example:
    gdf of type 'pipe' contains columns ['from_junction', 'to_junction', 'ID']
    drop_old_cols = False, tidy_up = True

    -> gdf is returned with 
        * old values (from input state of gdf) in columns ['from_junction', 'to_junction', 'ID', 'geometry']
        * default values in all other columns from dict_attr_val.keys()
        * no other additional columns
    
    """
    
    # Imports and function definitions
    import numpy as np

    # Check for correct type
    types_implemented = ['junction', 'pipe', 'sink', 'valve', 'heat_consumer']

    if type not in types_implemented:
        print(f'### Attention! Please provide proper type from {types_implemented}. Aborting ...###')
        return()

    if drop_old_cols:
        cols = list(gdf.columns)
        cols.remove('geometry')
        gdf = gdf.drop(columns = cols)

    if type == 'junction':
        dict_attr_val = {
            'ID':np.arange(len(gdf)).astype(int),
            'height_m':[0],
            'name':['KNO' + str(n) for n in np.arange(len(gdf)).astype(int)],
            'type':[type],
            'in_service':[True],
            'pn_bar':[1],
            'tfluid_k':[273.15+10]
        }

    if type == 'valve':
        dict_attr_val = {
            'ID':np.arange(len(gdf)).astype(int),
            'name':['valve' + str(n) for n in np.arange(len(gdf)).astype(int)],
            'type':[type],
            'opened':[True],
            'from_junction':[0],
            'to_junction':[0],
            'diameter_m':[0.05],
            'loss_coefficient':[0]
        }

    if type == 'pipe':

        # Assign default value for line length if no geometry is provided
        if 'geometry' not in gdf.columns:
            print('\n### Column "geometry" could not be found in gdf.columns. Length values are set to default value and do not represent\nactual geometric properties of lines. ###')
            gdf['length_m'] = 1 
        else:
            gdf['length_m'] = gdf['geometry'].length

        dict_attr_val = {
            'ID':np.arange(len(gdf)).astype(int),
            'from_junction':[0],
            'to_junction':[0],
            'name':['rlaV' + str(n) for n in np.arange(len(gdf)).astype(int)],
            'type':[type],
            'loss_coefficient':[0],
            'sections':[1],
            'text_k':[273.15+10],
            'qext_w':[0],
            'in_service':[True],
            'nominalWidth':[50],
            'length_km':gdf['length_m'] / (1e03),
            'diameter_m':[0.05],
            'k_mm':[0.0469],
            'alpha_w_per_m2k':[1.433506]
        }

    if type == 'sink':
        dict_attr_val = {
            'ID':np.arange(len(gdf)).astype(int),
            'junction':[0],
            'name':['sink' + str(n) for n in np.arange(len(gdf)).astype(int)],
            'type':[type],
            'in_service':[True],
            'mdot_kg_per_s':[0.02],
            'scaling':[1]
        }    

    if type == 'heat_consumer':
        dict_attr_val = {
            'ID':np.arange(len(gdf)).astype(int),
            'from_junction':[0],
            'to_junction':[0],
            'controlled_mdot_kg_per_s':[0.05],
            'deltat_k':[None],
            'diameter_m':[0.05],
            'treturn_k':[None],
            'qext_w':[1000],
            'name':['heat_consumer'+str(int(n)) for n in range(len(gdf))]
        }


    existing_cols = []

    for attr in dict_attr_val:
        if attr not in gdf.columns:
            gdf[attr] = dict_attr_val[attr] if len(dict_attr_val[attr]) == len(gdf.index) else dict_attr_val[attr]*len(gdf.index)
        else:
            existing_cols += [attr]

    if len(existing_cols) > 0:
        print(f'\n### The following columns were found in the provided gdf {type} and are NOT REPLACED with default values:\n{existing_cols} ###\n')

    if tidy_up:
        gdf = gdf[[n for n in dict_attr_val.keys()] + ['geometry']]

    return (gdf)

def createUniqueJunctions(
        lines:gp.GeoDataFrame, 
        splitLines:bool = True
    ):

    """
    Function that creates a set of unique junctions at all start- and end-coordinates of the given GeoDataFrame of line objects **lines**

    :param lines: GeoDataFrame containing lines
    :param splitLines: Boolean if lines shall be splitted at identified crossings of lines (if existent).
    :return: GeoDataFrame of all **junctions** and GeoDataFrame of (separated) lines.
    """

    # Imports
    from shapely.geometry import Point
    import pandas as pd

    # Extract all intermediate points from lines
    intermediatePoints, intermediateCoords = extractPointsFromLines(lines = lines, onlyIntermediate = True)

    points_list_a               = [[Point(pp.coords[0]), Point(pp.coords[-1])] for pp in lines['geometry']]  # List of (start/end) points of newly added lines in observed area
    points_list_b               = [val for sublist in points_list_a for val in sublist]
    points_list_b               += intermediatePoints
    points_list                 = np.array(points_list_b)

    # A_new                       = np.array([(k.x, k.y) for k in points_list])
    A_new                       = np.array([(Point(k).x, Point(k).y) for k in points_list])

    neighbours, _                = nnearest(A_new, A_new, distance=0.1, n=10)

    a                           = [tuple(sorted(list(l))) for l in neighbours]
    s                           = set(a) # -> Remove duplicate points
    idxs                        = [a.index(i) for i in s]
    points_new_list             = [Point(points_list[i]) for i in idxs] # -> All newly created points, including start-, end- and intermediate points of lines. Duplicates are already removed!


    # Identify all start- and endpoints BEFORE SPLITTING
    all_startpoints         = [(Point(pp.coords[0]).x, Point(pp.coords[0]).y) for pp in lines['geometry']]
    all_endpoints           = [(Point(pp.coords[-1]).x, Point(pp.coords[-1]).y) for pp in lines['geometry']]

    startendpoints = list(set(all_startpoints + all_endpoints))

    startendjunctions = gp.GeoDataFrame(geometry = [Point(n) for n in startendpoints])

    # Intermediate points
    intermediate_points = list(set(intermediateCoords))
    intermediatejunctions = gp.GeoDataFrame(geometry = [Point(n) for n in intermediate_points])

    # Identify crossing junctions (later used for separation of lines at these junctions) -> Crossing junctions can only exist within intermediate points
    alljunctions                     = gp.GeoDataFrame(geometry = points_new_list)
    alljunctions['splitting_points'] = alljunctions.apply(lambda p: 
                                        True if ((p.geometry.x, p.geometry.y) in intermediateCoords and (p.geometry.x, p.geometry.y) in startendpoints) or (intermediateCoords.count((p.geometry.x, p.geometry.y)) > 1) else False, axis = 1)


    # Creating multipoint object of splitting points/crossings and split lines
    splitting_points = alljunctions[alljunctions['splitting_points'] == True].dissolve()

    # Split lines at splitting_points if desired
    if (len(splitting_points) > 0) & (splitLines):    
        lines = split_lines(lines, splitting_points.loc[0, 'geometry'])

    # Put all start points' geometries of pipes into an array
    startpoints = [(np.round(Point(pp.coords[0]).x,2), np.round(Point(pp.coords[0]).y,2)) for pp in lines['geometry']]
    endpoints = [(np.round(Point(pp.coords[-1]).x,2), np.round(Point(pp.coords[-1]).y,2)) for pp in lines['geometry']]

    startendpoints = list(set(startpoints + endpoints))

    junctions = gp.GeoDataFrame(geometry = [Point(n) for n in startendpoints])


    return junctions, lines

def assignJunctionsToLines(
        lines:gp.GeoDataFrame,
        junctions:gp.GeoDataFrame,
        FromToAttributes:list = ('from_junction', 'to_junction')
    ):

    """
    Function that assigns the index of corresponding junctions at start- and end-coordinates at each line in **lines**.
    Point objects are given in GeoDataFrame **junctions**

    :param lines: GeoDataFrame with line objects
    :param junctions: GeoDataFrame with point objects
    :param FromToAttributes: tuple of attribute names at which to insert junction indices to **lines**
    :return: Updated GeoDataFrame of **lines**
    """

    # Imports
    from shapely.geometry import Point

    # Put all junctions into an array
    all_junctions = np.array([(j.x, j.y) for j in junctions['geometry']])

    startpoints = [(np.round(Point(pp.coords[0]).x,4), np.round(Point(pp.coords[0]).y,4)) for pp in lines['geometry']]
    endpoints = [(np.round(Point(pp.coords[-1]).x,4), np.round(Point(pp.coords[-1]).y,4)) for pp in lines['geometry']]

    # Find closest junctions to start- and endpoints of lines within a distance of 0.02m
    from_junction_idxs, _ = nnearest(startpoints, all_junctions, distance = 0.02, n=1)
    to_junctions_idxs, _ = nnearest(endpoints, all_junctions, distance = 0.02, n=1)

    lines[FromToAttributes[0]] = from_junction_idxs
    lines[FromToAttributes[1]] = to_junctions_idxs

    # Check if all junctions are connected to a pipe
    if len(set(list(np.concatenate((lines[FromToAttributes[0]].values,lines[FromToAttributes[1]].values))))) - len(junctions) != 0:
        print('### Warning! There are unconnected junctions! ###')

    else:
        print('\n### The set of start- and end-junctions of all pipes coincides with the total sum of created junctions. Continuing... ###\n')

    
    return lines

def assignJunctionsToBuildings(
        junctions:gp.GeoDataFrame,
        buildings:gp.GeoDataFrame,
        junctionAttribute:str = 'junction'
    ):

    """
    Function matching junction row indices from GeoDataFrame **junctions** to buildings,r epresented by polygons in **buildings**.

    :param junctions: GeoDataFrame of junctions as shapely.Point objects
    :param buildings: GeoDataFrame of polygons as shapely.Polygon objects

    :return: GeoDataFrame of building polygons
    """

    buildings_out = match_polygons_to_points_by_intersection(polygons = buildings, points = junctions)

    buildings_out.rename(columns = {'index_right':junctionAttribute}, inplace = True)

    return buildings_out

def create_basic_network_topology(
        lines:gp.GeoDataFrame,
        split_lines:bool = True,
        ID_col:str = 'ID',
        cs = 'EPSG:25832'
    ) -> gp.GeoDataFrame:

    """
    Function that takes simple geometric information of line objects in **lines** and converts them into a network topology of lines and junctions.\n
    Global precision is limtied to three decimal digits for simplification of junction snapping and line splitting.\n

    :param lines: GeoDataFrame of line objects
    :param split_lines: Boolean if lines shall be separated at identified crossings with other lines
    :param ID_col: str denoting attribute name for unique ID which is added to the lines
    :param cs: str denoting coordinate system which shall be applied to the output DataFrames
    :return: Final network topology of lines and junctions and their relative connection
    """

    ### Data preparation ###

    # Remove empty geometries in lines
    lines                           = lines[(~lines['geometry'].isnull())]
    lines                           = lines.explode(index_parts = False)
    lines.reset_index(drop = True, inplace = True)
    lines = lines.set_crs(cs)

    # Round to 3 decimals to reduce global precision and simplify junction snapping
    lines.geometry = lines.geometry.apply(lambda x: round_coords_2D(geom = x, ndigits = 3))
    print(f'\n... Attention! Global precision of line coordinates is limited to 3 decimal digits!')
    
    ### Geometric operations ###

    # Create set of junctions from all sart- and endpoints in lines and remove duplicates
    junctions, lines = createUniqueJunctions(lines, splitLines = split_lines)

    junctions = junctions.set_crs(cs)

    ### Assign junction IDs to lines ###
    lines = assignJunctionsToLines(lines = lines, junctions = junctions, FromToAttributes=('from_junction', 'to_junction'))

    lines[ID_col] = np.arange(len(lines))
    lines[ID_col] = lines[ID_col].astype(int)

    lines.reset_index(drop = True, inplace = True)

    print('######### Assignment of junctions to pipes DONE #########') 
    
    
    return(lines, junctions)

def insertValvesInLines(
        valves:gp.GeoDataFrame, 
        lines:gp.GeoDataFrame, 
        junctions:gp.GeoDataFrame, 
        FromToAttributes:list = ('from_junction', 'to_junction'),
        cs:str = 'EPSG:25832'
        ):

    """
    Function that inserts valves into lines at specified positions.

    :param valves: GeoDataFrame of valves as shapely.Point objects. To insert these into the lines, they have to\ncoincide with intermediate points of **lines**.
    :param lines: GeoDataFrame containing lines.
    :param junctions: GeoDataFrame containing junctions as  shapely.Point objects
    :param FromToAttributes: tuple containing attribute names in **lines** indicating the connection to\nrow indices for start- and end-junctions from **junctions**
    :param cs: coordinate system specified for the output GeoDataFrames
    :return: GeoDtaaFrames of lines, junctions and valves with newly created elements after splitting and inserting.
    """   

    
    # Imports
    import numpy as np
    import pandas as pd
    import geopandas as gp


    # Initialize output data
    lines = lines.copy()

    junctions_out = junctions.copy()
    lines_out = lines.copy()
    valves_out = valves.copy().set_crs(cs)

    # Check for unplausible vale geometries
    valves_out = valves_out[~valves_out.geometry.isnull()]

    # Initialisation
    reRun = True

    while reRun:

        # List of pipe indices which have intersections with multiple valves -> Only one valve intersection per pipe can be split up in this function
        original_line_idx_matches           = []
        original_idx_line_multiple_valves   = []


        for idx, valve in valves_out.iterrows():
            point = valve.geometry

            ### Check if valve is located on intermediate point within pipe
            line = lines[lines['geometry'].distance(point) < 0.01]
            junction = junctions[junctions['geometry'].distance(point) < 0.01]

            if (len(line) == 1) & (len(junction) == 0):

                # Check if pipe has already been edited and plit up with valves
                original_line_idx = line.index

                if original_line_idx[0] in original_line_idx_matches:
                    original_idx_line_multiple_valves.append(original_line_idx[0])
                    continue

                original_line_idx_matches.append(original_line_idx[0])

                ### Get information from original network components pipe and junctions
                from_junction_idx = line[FromToAttributes[0]].values[0]
                to_junction_idx = line[FromToAttributes[1]].values[0]
                
                # Create two new junctions which are added to the output DF of junctions
                max_ID_junctions = max(junctions_out.index)
                new_junctions = pd.DataFrame(columns = junctions.columns)

                # 2 new junctions per valve are created
                new_junctions = pd.concat((new_junctions, pd.DataFrame(valve).transpose()[['geometry']]))
                new_junctions = pd.concat((new_junctions, pd.DataFrame(valve).transpose()[['geometry']]))

                new_junctions.index = np.arange(max_ID_junctions + 1, len(new_junctions) + max_ID_junctions + 1)

                # Create DF for new pipes -> One pipe becomes two pipes after splitting
                max_ID_lines = max(lines_out.index)
                new_lines = line.copy()
                new_lines = pd.concat((new_lines, line)).reset_index(drop = True)

                # Split the pipe which contains the valve location
                res = split_lines_at_points(line['geometry'].values[0], point)
                objs = [n for n in res.geoms]

                new_lines['geometry'] = objs

                # Assign corresponding new junctions to newly added pipes
                new_lines.loc[0, FromToAttributes[0]] = from_junction_idx
                new_lines.loc[0, FromToAttributes[1]] = new_junctions.index[0]

                new_lines.loc[1, FromToAttributes[0]] = new_junctions.index[1]
                new_lines.loc[1, FromToAttributes[1]] = to_junction_idx

                new_lines.set_index(np.arange(max_ID_lines +1, len(new_lines) + max_ID_lines + 1), inplace = True)
                new_lines = new_lines.set_crs(cs)

                # Add newly added pipes to lines_out, drop original pipe
                lines_out = pd.concat((lines_out, new_lines))
                lines_out.drop(original_line_idx, inplace = True)

                # Add newly added junctions to junctions_out
                junctions_out = pd.concat((junctions_out, gp.GeoDataFrame(new_junctions, crs= cs)), ignore_index = False)   

                # Connect valve element to newly created junctions
                valves_out.loc[idx, FromToAttributes[0]] = new_junctions.index[0]
                valves_out.loc[idx, FromToAttributes[1]] = new_junctions.index[1]

        if len(original_idx_line_multiple_valves) > 0:
            print('\n### Attention! Multiple valves have been found to intersect with lines with the following original indices:\n',
                'Line indices: ' + str(original_idx_line_multiple_valves) + '\n',
                'The function is re-run in order to correctly separate all lines. ###\n')

            reRun = True

        else:
            reRun = False

        # Reset index of pipes DataFrame
        lines_out.reset_index(drop = True, inplace = True)


    return (junctions_out, lines_out, valves_out)

def identify_connection_lines(
        lines:gp.GeoDataFrame,
        buildings:gp.GeoDataFrame,
        connectionTypeAttr:str = 'connectionType',
        connAttrLines:list = ['from_junction', 'to_junction'],
        connAttrBuildings:str = 'junction',
        naming:list = ('distribution', 'houseConnection')
        )->gp.GeoDataFrame:

    """
    Function to assign attributes to **lines** based on their topological connection to junctions referred to in **buildings**' attribute **connAttrBuildings**.
    If identical values found in **buildings**' attribute **junction** and any of the attribtues given in **lines**' attributes **connAttrLines**, the line is marked as connectionType 'houseConnection', 'distribution' else

    :param  lines: GeoDataFrame of lines
    :param buildings: GeoDataFrame of buildings, defining the discrimination between distribution lines and house connection lines
    :param connectionTypeAttr: str denoting the attribute name which is added to **lines** for discriminating between distribution lines and house connection lines
    :param connAttrBuildings: str denoting attribute name in buildings indicating the row index of connected junction
    :param connAttrLines: str denoting the row indices of connected junctions in **lines** - correspondent attribute to **connAttrBuildings** in **buildings**
    :param naming: list of strings denoting the type distribution line (first entry) or house connection line (second entry) which is added to the output line DataFrame.

    :return: GeoDataFrame of lines 
    """

    # Initialise column
    lines[connectionTypeAttr] = naming[0] # Default value is naming[0]

    # Create temporary copy
    ls = lines.copy()
    ls.reset_index(drop = True, inplace = True)

    # Create dicitonary translating indices
    dictTemp = dict(zip(list(ls.index), list(lines.index)))

    lens = len(ls)

    BJs = np.array(buildings[connAttrBuildings])
    allJs = np.array(ls[connAttrLines])

    # Index, from- and to- junctions of all lines which share a from- or to-junction with the buildings
    idxBs = np.array([ [idx, allJs[idx, 0], allJs[idx, 1]] for idx in range(lens) if any([allJs[idx, k] in BJs for k in [0, 1]]) ])

    # Create container
    conns = []

    # Iterate through all lines which have intersection of from- or to-junction with buildings
    for jj in idxBs:
        otherJs = set(np.delete(allJs, jj[0], axis = 0).flatten())

        # Check if any from- or to-junction icoincides with buildings and is NOT contained in set of all remaining junctions of other lines
        if ((jj[1] not in otherJs) & (jj[1] in BJs)) | ((jj[2] not in otherJs) & (jj[2] in BJs)):
            conns.append(naming[1])
        else:
            conns.append(naming[0])

    trans = [dictTemp[k] for k in idxBs[:, 0]]
    
    # Transfer to original DataFrame
    lines.loc[trans, connectionTypeAttr] = conns

    return lines

def create_ppi_network_components(
    network_name:str = 'net', 
    fluid:str = 'water',
    **components):

    """
    Function that creates pandapipes network object from provided pandas.DataFrames of desired network components.

    :param network_name: str denoting the network name.
    :param fluid: str denoting the fluid employed in the network.
    :param components: kwargs for network components which shall be added to the network. These have to be passed as pandas.DataFrame
        possible kwargs for network components:
            - 'junctions', 'pipes', 'valves', 'sinks', 'heat_consumers'
            - 'dict_attributes': dictionary of key, value-pairs for additional attributes from the components' DataFrames which shall be added to the pandapipes network component's DataFrame.
    """


    # Imports
    import pandapipes as ppi
    
    
    net = ppi.create_empty_network(name=network_name, fluid=fluid)

    if 'dict_attributes' in components.keys():
        dic = components['dict_attributes']
        del components['dict_attributes']
    else:
        dic = {
            'defaultKey':'defaultValue'
        }

    # Example structure of dict_attributes
    # dict_attributes = {
    #     'junctions':column_name,
    #     'pipes':column_name2
    # }

    for cc in components.keys():
        
        if (components[cc] is not None):
            if (not components[cc].empty):

                print('Modelling of ' + cc)

                ### Junctions ###
                if cc == 'junctions':

                    # Add geodata from geometry if not existent
                    if ('geometry' in components[cc].columns) & ('geodata' not in components[cc].columns):
                        components[cc] = geodata_from_geometry(components[cc], typ = 'point')

                    ppi.create_junctions(
                        net=net,
                        nr_junctions=len(components[cc]['ID']),
                        pn_bar=list(components[cc]['pn_bar']),
                        tfluid_k=list(components[cc]['tfluid_k']),
                        geodata=list(components[cc]['geodata']),
                        # geometry=list(components[cc]['geometry']),
                        height_m=list(components[cc]['height_m']),
                        name=list(components[cc]['name']),
                        ID=list(components[cc]['ID'])
                    )

                    if cc in list(dic.keys()):
                        for val in dic[cc]:
                            
                            if val not in components[cc].columns:
                                print('\n### Attention! Attribute ' + str(val) + ' not found in '+ str(cc) + ' columns.\n',
                                        'Skipping this attribute. ###\n')                            
                                continue

                            else:
                                dtyp = components[cc][val].dtype
                                net.junction[val] = list(components[cc][val])
                                net.junction[val] = net.junction[val].astype(dtyp)
                
                ### Pipes ###
                if cc == 'pipes':

                    # Add geodata from geometry if not existent
                    if ('geometry' in components[cc].columns) & ('geodata' not in components[cc].columns):
                        components[cc] = geodata_from_geometry(components[cc], typ = 'line')

                    ppi.create_pipes_from_parameters(
                        net,
                        from_junctions=list(components[cc]['from_junction']),
                        to_junctions=list(components[cc]['to_junction']),
                        length_km=list(components[cc]['length_km']),
                        diameter_m=list(components[cc]['diameter_m']),
                        k_mm=list(components[cc]['k_mm']),
                        alpha_w_per_m2k=list(components[cc]['alpha_w_per_m2k']),
                        loss_coefficient=list(components[cc]['loss_coefficient']),
                        sections=list(components[cc]['sections']),
                        geodata=list(components[cc]['geodata']),
                        # geometry=list(components[cc]['geometry']),
                        in_service=list(components[cc]['in_service']),
                        name=list(components[cc]['name']),
                        ID=list(components[cc]['ID'])
                    )

                    if cc in list(dic.keys()):
                        for val in dic[cc]:
                            if val not in components[cc].columns:
                                print('\n### Attention! Attribute ' + str(val) + ' not found in '+ str(cc) + ' columns.\n',
                                        'Skipping this attribute. ###\n')                            
                                continue

                            else:
                                dtyp = components[cc][val].dtype
                                net.pipe[val] = list(components[cc][val])
                                net.pipe[val] = net.pipe[val].astype(dtyp)

                ### Valves ###
                if cc == 'valves':
                    ppi.create_valves(
                        net=net,
                        from_junctions=list(components[cc]['from_junction']),
                        to_junctions=list(components[cc]['to_junction']),
                        diameter_m=list(components[cc]['diameter_m']),
                        opened=list(components[cc]['opened']),
                        loss_coefficient=list(components[cc]['loss_coefficient']),
                        name=list(components[cc]['name']),
                        # geometry=list(components[cc]['geometry']),
                        ID=list(components[cc]['ID'])
                    )

                    if cc in list(dic.keys()):
                        for val in dic[cc]:
                            if val not in components[cc].columns:
                                print('### Attention! Attribute ' + str(val) + ' not found in '+ str(cc) + ' columns.\n',
                                        'Skipping this attribute. ###')                            
                                continue

                            else:
                                dtyp = components[cc][val].dtype
                                net.valve[val] = list(components[cc][val])
                                net.valve[val] = net.valve[val].astype(dtyp)

                ### Mass flow controllers ###
                if cc == 'flow_controllers':
                    ppi.create_flow_controls(
                        net = net,
                        from_junctions = list(components[cc]['from_junction']),
                        to_junctions=list(components[cc]['to_junction']),
                        diameter_m=list(components[cc]['diameter_m']),
                        controlled_mdot_kg_per_s = list(components[cc]['mdot_kg_per_s']),
                        control_active = list(components[cc]['control_active']),
                        ID=list(component['ID']),
                        name=list(component['name']),
                        in_service=list(component['in_service'])
                    )
                    
                    if cc in list(dic.keys()):
                        for val in dic[cc]:
                            if val not in components[cc].columns:
                                print('\n### Attention! Attribute ' + str(val) + ' not found in '+ str(cc) + ' columns.\n',
                                        'Skipping this attribute. ###\n')                            
                                continue

                            else:

                                dtyp = components[cc][val].dtype
                                net.flow_control[val] = list(components[cc][val])
                                net.flow_control[val] = net.flow_control[val].astype(dtyp)



                ### Heat exchangers ###
                if cc == 'heat_exchangers':
                    for idx, heat in components[cc].iterrows():
                        ppi.create_heat_exchanger(
                            net=net,
                            from_junction=heat['from_junction'],
                            to_junction=heat['to_junction'],
                            diameter_m=heat['diameter_m'],
                            qext_w=heat['qext_w'],
                            loss_coefficient=heat['loss_coefficient'],
                            name=heat['name'],
                            ID=heat['ID'],
                            in_service=heat['in_service']
                        )


                ### Circ pump pressure ###
                if cc == 'circ_pump_pressure':
                    for idx, producer in components[cc].iterrows():
                        ppi.create_circ_pump_const_pressure(
                            net=net,
                            from_junction=producer['to_junction'],
                            to_junction=producer['from_junction'],
                            p_bar=producer['p_bar'],
                            plift_bar=producer['plift_bar'],
                            t_k=producer['t_k'],
                            in_service=producer['in_service'],
                            type=producer['type'],
                            ID=producer['ID'],
                            name=producer['name']
                        )


                ### Circ pump mass flow ###
                if cc == 'circ_pump_mass_flow':
                    for idx, producer in components[cc].iterrows():
                        ppi.create_circ_pump_const_mass_flow(
                            net=net,
                            from_junction=producer['to_junction'],
                            to_junction=producer['from_junction'],
                            p_bar=producer['p_bar'],
                            mdot_kg_per_s=producer['mdot_kg_per_s'],
                            t_k=producer['t_k'],
                            in_service=producer['in_service'],
                            type=producer['type'],
                            ID=producer['ID'],
                            name=producer['name']
                        )


                ### Pumps ###
                if cc == 'pumps':
                    for idx, pump in components[cc].iterrows():
                        ppi.create_pump(
                            net=net,
                            from_junction=pump['from_junction'],
                            to_junction=pump['to_junction'],
                            std_type=pump['std_type'],
                            in_service=pump['in_service'],
                            type=pump['type'],
                            ID=pump['ID'],
                            name=pump['name']
                        )


                ### Sources ###
                if cc == 'sources':
                    ppi.create_sources(
                        net=net,
                        junctions=list(components[cc]['junction']),
                        mdot_kg_per_s=list(components[cc]['mdot_kg_per_s']),
                        scaling=list(components[cc]['scaling']),
                        ID=list(components[cc]['ID']),
                        name=list(components[cc]['name']),
                        in_service=list(components[cc]['in_service']),
                        type='source'
                    )


                ### Sinks ###
                if cc == 'sinks':

                    component = components[cc][~components[cc]['junction'].isna()].reset_index(drop = True)

                    ppi.create_sinks(
                        net=net,
                        junctions=list(component['junction']),
                        mdot_kg_per_s=list(component['mdot_kg_per_s']),
                        scaling=list(component['scaling']),
                        ID=list(component['ID']),
                        name=list(component['name']),
                        in_service=list(component['in_service']),
                        type='sink'
                    )

                    if cc in list(dic.keys()):
                        for val in dic[cc]:
                            
                            if val not in components[cc].columns:
                                print('\n### Attention! Attribute ' + str(val) + ' not found in '+ str(cc) + ' columns.\n',
                                        'Skipping this attribute. ###\n')                            
                                continue

                            else:
                                dtyp = components[cc][val].dtype
                                net.sink[val] = list(component[val])
                                net.sink[val] = net.sink[val].astype(dtyp)


                ### Ext grids ###
                if cc == 'ext_grid':
                    for idx, ext_grid in components[cc].iterrows():
                        ppi.create_ext_grid(
                            net=net,
                            junction=ext_grid['junction'],
                            p_bar=ext_grid['p_bar'],
                            t_k=ext_grid['t_k'],
                            name=ext_grid['name'],
                            in_service=ext_grid['in_service'],
                            type=ext_grid['type']
                        )


                ### Heat consumers (pandapipes version 0.10.0 onwards)
                if cc == 'heat_consumers':
                    ppi.create_heat_consumers(
                        net = net,
                        from_junctions = list(components[cc]['from_junction']),
                        to_junctions = list(components[cc]['to_junction']),
                        diameter_m = list(components[cc]['diameter_m']),
                        controlled_mdot_kg_per_s = list(components[cc]['controlled_mdot_kg_per_s']),
                        qext_w = list(components[cc]['qext_w'])              
                    )

                    if cc in list(dic.keys()):
                        for val in dic[cc]:
                            
                            if val not in components[cc].columns:
                                print('\n### Attention! Attribute ' + str(val) + ' not found in '+ str(cc) + ' columns.\n',
                                        'Skipping this attribute. ###\n')                            
                                continue

                            else:
                                dtyp = components[cc][val].dtype
                                net.heat_consumer[val] = list(components[cc][val])
                                net.heat_consumer[val] = net.heat_consumer[val].astype(dtyp)

    return(net)

def create_ppi_network(
        network_name:str = 'net', 
        fluid:str = 'water',
        modelling_type:str = 'feed_line',
        consumerFeedLineJunctionIdx:list = [],
        **components):
    
    
    """
    Function that automatically creates pandapipes networks from GeoDataFrames with topology information (pipes-junction-structure) as inputs.
    Modelling is based on symmetrical feed line and reflux line within the network!

    :param network_name: str denoting the desired name of the pandapipes network.
    :param fluid: str denoting the desired fluid employed in the pandapipes network.
    :param modelling_type: str denoting the moedlling type, either "feed_line", "feed_and_reflux_line" or "feed_and_reflux_line_circPumps".
    :param consumerFeedLineJunctionIdx: When modelling feed and reflux line based on a nexisting network of feed line topology, the indexes of heat consumers' connected junctions can be passed here.
    :param components: kwargs for network components which shall be added to the network. These have to be passed as pandas.DataFrame
        possible kwargs for network components:
            - 'junctions', 'pipes', 'valves', 'sinks', 'heat_consumers'
            - 'dict_attributes': dictionary of key, value-pairs for additional attributes from the components' DataFrames which shall be added to the pandapipes network component's DataFrame.
    
    **modelling_type** can either be "feed_line" only or "feed_and_reflux_line".
    "components" additionally takes "dict_attributes" as input argument. This adds all contained attributes (values) for the component (key) to the ppi network structure.
        - "feed_line"
            * consumers/connected buildings are modelled as mass flow sinks
            * only basic hydraulic calculations may possible with this modelling approach
        - "feed_and_reflux_line"
            * consumers/connected buildings are modelled as mass flow controllers in combination with heat exchangers (heat_consumers)
            * production sites are modelled as circulation pumpus with constant pressure and set temperature (peak load producer)
            * further production sites (base load producers) are modelled as flow controllers and heat exchangers with supply to the feed_line
            * For further information consult https://github.com/e2nIEE/pandapipes/issues/309
        - "feed_and_reflux_line_circPumps"
            * consumers/connected buildings are modelled as mass flow controllers in combination with heat exchangers (heat_consumers)
            * production sites are modelled as circulation pumps with constant pressure and set temperature (peak load producer)
            * further production sites (base load producers) are modelled as circulation pumps with constant mass flow and downstream flow controllers.
            * For further information consult https://github.com/e2nIEE/pandapipes/issues/309
    
    For further information on default and pre-installed fluids and data structure, consult https://pandapipes.readthedocs.io/en/latest/
    """

    ### Imports
    import pandapipes as ppi

    if 'dict_attributes' in components.keys():
        dic = components['dict_attributes']
    else:
        dic = {
            'defaultKey':'defaultValue'
        }

    
    ### Return network or extend it with reflux line according to chosen modelling type
    if modelling_type in ('feed_line'):     

        # Create feed line network
        net = create_ppi_network_components(
                network_name=network_name,
                fluid=fluid,
                **components
                )    

        print('\n### Pandapipes network for one-directional flow (feed line only) is created.\n',
              'Junctions, pipes, valves are created.\n',
              'Consumers/buildings are modelled as mass flow sinks.\n',
              'Production sites are modelled as external grids (peak load producer) or as mass flow sources. ###\n')

        return (net)
    

    elif modelling_type in ('feed_and_reflux_line', 'feed_and_reflux_line_circPumps'):

        # Temporarily remove components which can only be handled in networks with feed- and reflux line
        if 'heat_consumers' in components:
            components_feedLine = components.copy()
            del components_feedLine['heat_consumers']

        # Create feed line network
        net = create_ppi_network_components(
                network_name=network_name,
                fluid=fluid,
                **components_feedLine
                )    

        ## a) Copy feed line network and reindex all components
        net2 = net.deepcopy()
        ppitlbx.create_continuous_junction_index(net2, start = len(net.junction.index))

        ## b) Create network containing both feed line and reflux line
        # Store additional information in pandapipes structure
        net.junction['Layer'] = 'feedLine'
        net2.junction['Layer'] = 'refluxLine'

        net.pipe['Layer'] = 'feedLine'
        net2.pipe['Layer'] = 'refluxLine'

        # Add feed line and reflux line components to components_cpl
        components_cpl = dict()

        junctions_cpl = pd.concat((net.junction, net2.junction))
        pipes_cpl = pd.concat((net.pipe, net2.pipe))

        components_cpl['junctions'] = junctions_cpl
        components_cpl['pipes'] = pipes_cpl

        # Add valves if existent
        if hasattr(net, 'valve'):
            net.valve['Layer'] = 'feedLine'
            net2.valve['Layer'] = 'refluxLine'

            valves_cpl = pd.concat((net.valve, net2.valve))
            components_cpl['valves'] = valves_cpl


        # Create dictionaries for matching junctions between feed line and reflux line (symmetric topology) and update them as 'to_junction' attribute in heat_consumers
        if len(consumerFeedLineJunctionIdx) == 0:
            print('\n### No iterable for feed line junction indices of consumers is provided. ###')

            if hasattr(net, 'sink'):
                consFeedIdx = net.sink['junction']
                print('\n### Information about connected feed line junctions is taken from provided DataFrame of sinks in feed line network. ###')

            else:
                print('\n### No sinks provided for pandapipes network from which to draw information on feed line junction indices. Aborting ...###')
                return net

        else:
            consFeedIdx = consumerFeedLineJunctionIdx

        dict_consumer_junctions_VLRL = {
            idx:idx+len(net.junction.index) for idx in consFeedIdx
        }      

        if 'heat_consumers' in components.keys():
            # Update 'to_junction' attribute in heat_consumers
            components_cpl['heat_consumers'] = components['heat_consumers']
            components_cpl['heat_consumers']['from_junction'] = list(dict_consumer_junctions_VLRL.keys())
            components_cpl['heat_consumers']['to_junction'] = list(dict_consumer_junctions_VLRL.values())


        ## Create basic network topology with feed line and reflux line
        net_cpl = create_ppi_network_components(
                network_name=network_name,
                fluid=fluid,
                dict_attributes = dic,
                **components_cpl
            )  
            

        ## Printing
        if modelling_type == 'feed_line':
            textstr = '\n### Pandapipes network with feed line only is created.\n \
              Junctions, pipes, valves are created.\n \
              Consumers/buildings are modelled as mass flow sinks (sinks).\n \
              Production sites are modelled as external grids (peak load producer) or as mass flow sources (base load producer). ###\n'

        elif modelling_type == 'feed_and_reflux_line':
            textstr = '\n### Pandapipes network with feed and reflux line is created.\n \
              Junctions, pipes, valves are created.\n \
              Consumers/buildings are modelled as heat exchangers + mass flow controllers (heat_consumers).\n \
              Production sites are modelled as circ_pump_pressure (peak load producer) or as inverse heat consumer (base load producer). ###\n'

        elif modelling_type == 'feed_and_reflux_line_circPumps':
            textstr = '\n### Pandapipes network with feed and reflux line is created.\n \
              Junctions, pipes, valves are created.\n \
              Consumers/buildings are modelled as heat exchangers + mass flow controllers (heat_consumers).\n \
              Production sites are modelled as circ_pump_pressure (peak load producer) or as circ_pump_mass + flow controllers (base load producer). ###\n'

        print(textstr)
        return(net_cpl)
    
    else:
        print('\n### Attention! Please choose proper modelling type.\n',
              '"feed_line" and "feed_and_reflux_line" are valid entries.\n',
              'Aborting... ###\n')

def create_ppi_producers_from_GIS_data(
        net,        # Pandapipes network
        production_sites:gp.GeoDataFrame,
        modelling_type:str  = 'feed_line',
        dict_col_layer:dict = {'Layer':{'feed':'V', 'reflux':'R'}},
        dict_col_type:dict = {'Type':{'peak':'peak', 'base':'base'}},  
        namecol:str = 'name',      
        add_cols:list = ['Pth_kW'],
        cs = 'EPSG:25832',
        reindex_junctions:bool = False
    ):
    
    """
    Function that uses junction component table from pandapipes network **net** and gp.GeoDataFrame **production_sites** to match production sites to corresponding junction indices in the network.
    
    :param net: pandapipes network in which producer components shall be implemented
    :param modelling_type: either 'feed_line' or 'feed_and_reflux_line' or 'feed_and_reflux_line_circPumps' to control different modelling approach for peak and base load producers.    
    :param dict_col_layer: specifies the column (key of 1st dict) where to find information on feed line and reflux line junctions (keys of 2nd dict)
    :param dict_col_type: specifies the column (key of 1st dict) where to find type information on peak load and base load producers (keys of 2nd dict)
    :param add_cols: specifies additional/optional columns which shall be transferred from GeoDataFrame "production_sites" to the network components.
    :param reindex_junctions: optionally cretes continuous junctions index via pandapipes.toolbox.create_continuous_junction_index(net, start = 0)
    
    For further information, please consult https://github.com/e2nIEE/pandapipes/issues/309
    """
    
    ### Imports
    import geopandas as gp
    import pandas as pd
    import pandapipes as ppi
    import pandapipes.toolbox as ppitlbx

    ### Function code
    # Layer columns
    lay_col             = list(dict_col_layer.keys())[0]
    lay_col_attr_feed   = dict_col_layer[lay_col]['feed']
    lay_col_attr_reflux = dict_col_layer[lay_col]['reflux']

    # Type columns
    typ_col             = list(dict_col_type.keys())[0]
    typ_col_attr_peak   = dict_col_type[typ_col]['peak']
    typ_col_attr_base   = dict_col_type[typ_col]['base']

    # Create counter
    counter_baseload = 0

    # Name columns
    namelist = list(production_sites[namecol]) if (len(namecol) > 0) & (namecol in production_sites.columns) else ['prod' + str(nn) for nn in np.arange(len(production_sites))]

    # Identify junctions in feed line and reflux line connected to production sites
    production_sites.reset_index(drop = True, inplace = True)
    

    def _check_producer_junction_count(prods, junctions_gdf):
        # Print a specific error if a producer polygon intersects more than one junction (only the first is used).
        counts = gp.sjoin(prods[['geometry']], junctions_gdf[['geometry']], predicate = 'intersects', how = 'inner').groupby(level = 0).size()
        for idx, cnt in counts.items():
            if cnt > 1:
                print(f'\n### Error! Producer at index {idx} intersects {cnt} junctions. A producer must intersect exactly one junction; only the first is used. ###\n')
    ### Differ between modeling types: "feed_line" and "feed_and_reflux_line"
    if modelling_type == 'feed_line':
        print('\n### Peak load producers are modelled as external grids.\n \
              Base load producers are modelled as mass flow sources. ###')
        
        production_sites = assignJunctionsToBuildings(
            buildings = production_sites.to_crs(cs),
            junctions = gp.GeoDataFrame(net.junction, geometry = 'geometry').set_crs(cs),
            junctionAttribute = 'flow_junction'
            )

        _check_producer_junction_count(production_sites, gp.GeoDataFrame(net.junction, geometry = 'geometry').set_crs(cs))

        ### Create external grid at production sites with peak load producer
        for n, p in production_sites.iterrows():
            
            if p[typ_col] == typ_col_attr_peak:

                ppi.create_ext_grid(
                    net = net,
                    junction = p['flow_junction'],
                    p_bar = 1,
                    t_k = 273.15 + 70,
                    type = 'auto',
                    in_service = True,
                    name = namelist[n]
                )

                for val in add_cols:
                    if val not in production_sites.columns:
                        print(f'\n### Attention! Attribute {val} not found in columns.\n \
                                Skipping this attribute while creating producer with index {p.index}. ###\n')                           
                        continue

                    else:
                        dtyp = production_sites[val].dtype
                        net.ext_grid.loc[net.ext_grid.index[-1], val] = p[val]

            elif p[typ_col] == typ_col_attr_base:

                ppi.create_source(
                    net = net,
                    junction = p['flow_junction'],
                    mdot_kg_per_s = 1,
                    in_service = True,
                    ID = n,
                    name = namelist[n]
                )

                counter_baseload += 1

                for val in add_cols:
                    if val not in production_sites.columns:
                        print(f'\n### Attention! Attribute {val} not found in columns.\n \
                                Skipping this attribute while creating producer with index {p.index}. ###\n')                            
                        continue

                    else:
                        dtyp = production_sites[val].dtype
                        net.source.loc[net.source.index[-1], val] = p[val]

            else:
                print('\n### Attention! Producer at index ', str(n), ' does not contain information on peak or base load characteristics. Skipping this producer. ###')




        ### Printing
        if hasattr(net, 'ext_grid'):
            print('\n### Production sites in the network are created.\n',
                'Peak load producers are created as follows (' + str(int(len(net.ext_grid))) + ' in total):\n\n')
            print(net.ext_grid)
            print('\n')
        
        if hasattr(net, 'source'):
            print('\nBase load producers are created as follows (' + str(int(counter_baseload)) + ' in total):\n\n')
            print(net.source.iloc[-counter_baseload:])



    elif modelling_type == 'feed_and_reflux_line':
        print('\n### Peak load producers are modelled as circ_pump_pressure.\n \
                Base load producers are modelled as heat_exchangers + flow_controller. ###')
        
        production_sites = assignJunctionsToBuildings(
            buildings = production_sites.to_crs(cs),
            junctions = gp.GeoDataFrame(net.junction[net.junction[lay_col] == lay_col_attr_feed], geometry = 'geometry').set_crs(cs),
            junctionAttribute = 'flow_junction'
            )
        
        production_sites = assignJunctionsToBuildings(
            buildings = production_sites.to_crs(cs),
            junctions = gp.GeoDataFrame(net.junction[net.junction[lay_col] == lay_col_attr_reflux], geometry = 'geometry').set_crs(cs),
            junctionAttribute = 'reflux_junction'
            )

        _check_producer_junction_count(production_sites, gp.GeoDataFrame(net.junction[net.junction[lay_col] == lay_col_attr_feed], geometry = 'geometry').set_crs(cs))

        # Create counter
        counter_baseload = 0

        ### Create circulation pumps with constant pressure at peak load production sites
        for n, p in production_sites.iterrows():
            
            if p[typ_col] == typ_col_attr_peak:

                ppi.create_circ_pump_const_pressure(
                    net = net,
                    return_junction = p['reflux_junction'],
                    flow_junction = p['flow_junction'],
                    p_flow_bar = 10,
                    plift_bar = 5,
                    t_flow_k = 273.15 + 70,
                    type = 'auto',
                    in_service = True,
                    ID = n,
                    name = namelist[n]
                )

                
                for val in add_cols:
                    if val not in production_sites.columns:
                        print(f'\n### Attention! Attribute {val} not found in columns.\n \
                                Skipping this attribute while creating producer with index {p.index}. ###\n')                           
                        continue

                    else:
                        dtyp = production_sites[val].dtype
                        net.circ_pump_pressure.loc[net.circ_pump_pressure.index[-1], val] = p[val]

            ### Create heat exchangers and flow controller at base load production sites
            elif p[typ_col] == typ_col_attr_base:                

                ppi.create_heat_consumer(
                    net = net,
                    from_junction = p['reflux_junction'],
                    to_junction = p['flow_junction'],
                    diameter_m = 0.05,
                    qext_w = -1,
                    controlled_mdot_kg_per_s = 1,
                    name = namelist[n],
                    ID = n,
                    in_service = True
                )

                counter_baseload += 1


                for val in add_cols:
                    if val not in production_sites.columns:
                        print(f'\n### Attention! Attribute {val} not found in columns.\n \
                                Skipping this attribute while creating producer with index {p.index}. ###\n')                       
                        continue

                    else:
                        dtyp = production_sites[val].dtype
                        # net.circ_pump_mass.loc[net.circ_pump_mass.index[-1], val] = p[val] -> Modelling approach with circ_pump_mass... currently not implemented/not working
                        net.heat_consumer.loc[net.heat_consumer.index[-1], val] = p[val]

            else:
                print('\n### Attention! Producer at index ', str(n), ' does not contain information on peak or base load characteristics. Skipping this producer. ###')

        ### Printing
        if hasattr(net, 'circ_pump_pressure'):
            print('\n### Production sites in the network are created.\n',
                'Peak load producers are created as follows (' + str(int(len(net.circ_pump_pressure))) + ' in total):\n\n')
            print(net.circ_pump_pressure)
            print('\n')
        
        if hasattr(net, 'heat_consumer'):
            print('\nBase load producers are created as follows (' + str(int(counter_baseload)) + ' in total):\n\n')
            if counter_baseload != 0:
                print(net.heat_consumer.iloc[-counter_baseload:])

    elif modelling_type == 'feed_and_reflux_line_circPumps':
        print('\n### Peak load producers are modelled as circ_pump_pressure.\n \
                Base load producers are modelled as circ_pump_mass + downstream flow_controller. ###')
        
        production_sites = assignJunctionsToBuildings(
            buildings = production_sites.to_crs(cs),
            junctions = gp.GeoDataFrame(net.junction[net.junction[lay_col] == lay_col_attr_feed], geometry = 'geometry').set_crs(cs),
            junctionAttribute = 'flow_junction'
            )
        
        production_sites = assignJunctionsToBuildings(
            buildings = production_sites.to_crs(cs),
            junctions = gp.GeoDataFrame(net.junction[net.junction[lay_col] == lay_col_attr_reflux], geometry = 'geometry').set_crs(cs),
            junctionAttribute = 'reflux_junction'
            )
        
        _check_producer_junction_count(production_sites, gp.GeoDataFrame(net.junction[net.junction[lay_col] == lay_col_attr_feed], geometry = 'geometry').set_crs(cs))

        # Create counter for base load producers
        counter_baseload = 0
        
        ### Create circulation pumps with constant pressure at peak load production sites
        for n, p in production_sites.iterrows():
            
            if p[typ_col] == typ_col_attr_peak:

                ppi.create_circ_pump_const_pressure(
                    net = net,
                    return_junction = p['reflux_junction'],
                    flow_junction = p['flow_junction'],
                    p_flow_bar = 10,
                    plift_bar = 5,
                    t_flow_k = 273.15 + 70,
                    type = 'auto',
                    in_service = True,
                    ID = n,
                    name = namelist[n]
                )
                
                for val in add_cols:
                    if val not in production_sites.columns:
                        print(f'\n### Attention! Attribute {val} not found in columns.\n \
                                Skipping this attribute while creating producer with index {p.index}. ###\n')                            
                        continue

                    else:
                        dtyp = production_sites[val].dtype
                        net.circ_pump_pressure.loc[net.circ_pump_pressure.index[-1], val] = p[val]

            ### Create heat exchangers and flow controller at base load production sites
            elif p[typ_col] == typ_col_attr_base:

                # Identify junction connected to base load producer reflux
                reflux_js = net.junction[net.junction.index == int(p['reflux_junction'])]

                if len(reflux_js) > 1: # Ambiguous number of reflux junctions found at this producer
                    reflux_j = reflux_j.iloc[[0]]
                    print(f'\n### There are multiple possible reflux junctions found at base load producer {p.index} \
                          . Junction {reflux_j.index} is used. ###\n')
                    
                elif len(reflux_js) == 0:
                    print(f'\n### No possible reflux junctions found at base load producer {p.index} \
                          . Please check for proper intersection with producer geometry. ###\n')
                    continue

                elif len(reflux_js) == 1:
                    reflux_j = assign_default_values_ppi(gdf = reflux_js, type = 'junction', drop_old_cols = False)

                # Create additional (intermediate) junction at base load producers
                ppi.create_junction(
                    net = net,
                    pn_bar = reflux_j['pn_bar'].values[0],
                    tfluid_k = reflux_j['tfluid_k'].values[0],
                    height_m = reflux_j['height_m'].values[0],
                    name = 'KNO_intermediate_' + namelist[n],
                    type = 'junction',
                    in_service = True,
                    geodata = (net['junction_geodata'].loc[reflux_j.index, 'x'].values[0], net['junction_geodata'].loc[reflux_j.index, 'y'].values[0]) if hasattr(net, 'junction_geodata') else (p['geometry'].centroid.x, p['geometry'].centroid.y),
                    geometry = reflux_j['geometry'].values[0] if 'geometry' in reflux_j.columns else p['geometry'].centroid
                )

                intermediate_idx = net.junction.index[-1]

                # Create circulation pump with specified mass flow and downstream flow_controller for base load producer
                ppi.create_circ_pump_const_mass_flow(
                    net = net,
                    return_junction = p['reflux_junction'],
                    flow_junction = intermediate_idx,
                    p_flow_bar = 10,
                    mdot_flow_kg_per_s = 5,
                    t_flow_k = 273.15 + 70,
                    type = 'auto',
                    name = namelist[n],
                    in_service = True,
                    ID = n
                )

                for val in add_cols:
                    if val not in production_sites.columns:
                        print(f'\n### Attention! Attribute {val} not found in columns.\n \
                                Skipping this attribute while creating producer with index {p.index}. ###\n')                          
                        continue

                    else:
                        dtyp = production_sites[val].dtype
                        net.circ_pump_mass.loc[net.circ_pump_mass.index[-1], val] = p[val]

                # Create mass flow controller downstream of the circulation pump
                ppi.create_flow_control(
                    net = net,
                    from_junction = intermediate_idx,
                    to_junction = p['flow_junction'],
                    controlled_mdot_kg_per_s = 5,
                    control_active = True,
                    name = namelist[n] + '_flow_control',
                    in_service = True,
                    ID = n
                )

                counter_baseload += 1

                for val in add_cols:
                    if val not in production_sites.columns:
                        print('\n### Attention! Attribute ' + str(val) + ' not found in columns.\n',
                                'Skipping this attribute. ###\n')                            
                        continue

                    else:
                        dtyp = production_sites[val].dtype
                        net.flow_control.loc[net.flow_control.index[-1], val] = p[val]


                
                




        
        ## Create additional junction at base load producer site




    else:
        print('\n### Attention! Please choose proper modelling type.\n',
              '"feed_line" and "feed_and_reflux_line" are valid entries.\n',
              'Aborting... ###\n')
        
        
    ### Optional: Reindex junctions if desired
    if reindex_junctions:        
        ppitlbx.create_continuous_junction_index(net, start = 0)

    


    return(net)

def create_ppi_network_from_gdf(
        pipes:gp.GeoDataFrame = None,
        buildings:gp.GeoDataFrame = None,
        valves:gp.GeoDataFrame = None,
        producers:gp.GeoDataFrame = None,
        rasterHeight = None, # rasterio.Dataset
        existingNetworkFeedLine = None,
        modelling_type:str = 'feed_line',
        buildings_uniqueID:str = 'build_ID',
        heatingDemandAttr:str = 'demand_use_th',
        producerTypeDict:dict = {'Type':{'peak':'peak', 'base':'base'}},
        producerAddAttr:list = [],
        networkName:str = 'net1',
        networkFluid:str = 'water',
        checkUnconnectedComponents:bool = False
    ):

    """
    Function that creates complete pandapipes network object from geopandas.GeoDataFrame of lines, buildings, valves, production sites.

    :param pipes: GeoDataFrame of pipes geometry (shapely.Linestring)
    :param buildings: GeoDataFrame of buildings geometries representing heat consumers (shapely.Polygon)
    :param valves: GeoDataFrame of valve geometries (shapely.Point)
    :params producers: GeoDataFrame of production site geometries (shapely.Polygon)
    :param rasterHeight: rasterio.Dataset of raster values containing height values which shall be written to junctions in the network.
    :param existingNetworkFeedLine: Optional existing pandapipes network object of feed line topology from which network object with feed and reflux line shall be generated.
    :param modelling_type: str denoting the desired modelling_type ("feed_line", "feed_and_reflux_line", "feed_and_reflux_line_circPumps")
    :param buildings_uniqueID: str denoting attribute name of unique identifier for buildings.
    :param heatingDemandAttr: str denoting attribute name of heating demand in buildings which shall be transferred to network components.
    :param networkName: str denoting the desired network name.
    :param networkFluid: str denoting the fluid employed in the network.
    :param producerTypeDict: dictionary for matching of producer type to attribute name in DataFrame. Values 'peak' and 'base' are differentiated in column with name of the dictionary key.
    :param producerAddAttr: list of additional attribute names which shall be added to the pandapipes network components of producers.
    
    :returns: pandapipes network object
    """

    ### Case a): No existing feed line network is provided for creating symmetric reflux line ###
    if existingNetworkFeedLine is None:
        print('\n### Pandapipes network (feed line only) is generated from input GeoDataFrames of network components. ###')

        # Plausibility check for modelling_type
        if modelling_type not in ('feed_line'):
            print(f'\n### Attention! desired modelling_type {modelling_type} requires provided feed line network. Aborting... ###\n')
            return
        
        # Create basic network topology of pipes and junctions
        pipes, junctions = create_basic_network_topology(lines = pipes, split_lines = True)

        # Apply precision to line and point objects
        for obj in (pipes, valves, junctions):
            if obj is not None:
                obj = obj.copy()
                obj['geometry'] = set_precision(obj['geometry'], 0.01)

        # Insert valves if existent
        if valves is not None:
            junctions, pipes, valves = insertValvesInLines(lines = pipes, junctions = junctions, valves = valves, FromToAttributes=['from_junction', 'to_junction'])

        # Extract height values if raster data are passed
        if rasterHeight is not None:
            junctions = extractRasterValsAtPoints(j = junctions, raster = rasterHeight, statsType = 'min', buff_distance = 2, attrCol = 'height_m')

            if any(junctions['height_m'].isna()):
                print(f'\n### Attention!: Some junctions feature None values for height attribute "height_m". Check if provided raster dataset has correct CRS. ###\n')

        # Connect buildings to network junctions
        if buildings is not None:
            buildings_out = assignJunctionsToBuildings(buildings = buildings, junctions = junctions, junctionAttribute = 'junction')
            
            if any(buildings_out['junction'].isna()):
                print('\n### Attention! Buildings with following unique IDs could not be matched with junction: ', str(list(buildings_out.loc[buildings_out['junction'].isna(), buildings_uniqueID])), ' ###')
            
            sinks = buildings_out[[heatingDemandAttr, buildings_uniqueID, 'junction']].copy()
            sinks = sinks[~sinks['junction'].isna()].reset_index(drop = True)

            pipes = identify_connection_lines(
                lines = pipes, 
                buildings = buildings_out, 
                connectionTypeAttr = 'connectionType', 
                connAttrLines = ['from_junction', 'to_junction'], 
                connAttrBuildings = 'junction',
                naming = ('distribution', 'houseConnection')
                )
            
            print(f'\n### Pandapipes network pipes receive attribute "connectionType" to discriminate between "distribution" and "connection" pipes. ###')
            
        else:
            buildings_out = None
            sinks = None
        
        # Assign default values to network components
        junctions   = assign_default_values_ppi(gdf = junctions, type = 'junction', drop_old_cols = False)
        pipes       = assign_default_values_ppi(gdf = pipes, type = 'pipe', drop_old_cols = False, tidy_up = False)
        valves      = assign_default_values_ppi(gdf = valves, type = 'valve', drop_old_cols = False, tidy_up = False) if valves is not None else valves
        sinks       = assign_default_values_ppi(gdf = sinks, type = 'sink', drop_old_cols = False, tidy_up = False) if sinks is not None else sinks

        # Insert geodata information (geometry representation for pandapipes plotting)
        pipes       = geodata_from_geometry(pipes, typ = 'line')
        junctions   = geodata_from_geometry(junctions, typ = 'point')


        # Check for unconnected pipes
        if checkUnconnectedComponents:
            isConnected, subgdfs = checkConnectivity(pd.concat((pipes, valves), axis = 0) if valves is not None else pipes, source = 'from_junction', target = 'to_junction')

            if not isConnected:
                pipes['subgraph'] = list(subgdfs['subgraph'])

                print(f'\n### {len(set(subgdfs["subgraph"]))} unconnected networks are found in the provided line objects. Information are stored in column "subgraph" of line DataFrame. ###')


        ## Create pandapipes network for feed line
        # Attribute names in dict_attributes.values() are optional attrbutes which shall be added to the ppi network component table if they
        # exist in the provided component DataFrame.
        netFeedLine = create_ppi_network(
            network_name = networkName,
            fluid = networkFluid,
            modelling_type = modelling_type,
            junctions = junctions,
            pipes = pipes,
            valves = valves,
            sinks = sinks,
            dict_attributes = {
                'pipes':['connectionType', 'nominalWidth', 'Layer', 'geometry', 'subgraph'],
                'junctions':['Layer', 'ID', 'geometry'],
                'sinks':['Layer', 'geometry', heatingDemandAttr, buildings_uniqueID],
                'valves':['Layer', 'geometry']
            }
        )

        # Add production sites to network
        if producers is not None:
            netFeedLine = create_ppi_producers_from_GIS_data(
                net = netFeedLine,
                production_sites = producers,
                modelling_type = modelling_type,
                dict_col_type = producerTypeDict,
                add_cols = producerAddAttr
                )
        
        net = netFeedLine.deepcopy()

    ### Case b): Existing feed line network is provided for creating symmetric reflux line ###    
    else:
        print('\n### Pandapipes network (feed and reflux line) is generated from input network of feed line. ###')

        # Plausibility check for modelling_type
        if modelling_type not in ('feed_and_reflux_line', 'feed_and_reflux_line_circPumps'):
            print(f'\n### Attention! Desired modelling_type {modelling_type} should only be used for creating feed line network. The modelling_type is changed to \
                  feed_and_reflux_line_circPumps by default ... ###\n')
            modelling_type = 'feed_and_reflux_line_circPumps'

        netFeedLine = existingNetworkFeedLine.deepcopy()

        # Assign attribute "Layer" to junctions and pipes and valves to identify feed line and reflux line components
        if 'Layer' not in netFeedLine.junction.columns:
            netFeedLine.junction['Layer'] = 'feedLine'

        if 'Layer' not in netFeedLine.pipe.columns:
            netFeedLine.pipe['Layer'] = 'feedLine'

        if hasattr(netFeedLine, 'valve'):
            if 'Layer' not in netFeedLine.valve.columns:
                netFeedLine.valve['Layer'] = 'feedLine'

        # Assign trail_ID to pipes to match parallel pipes in feed line and reflux line 
        netFeedLine.pipe['trail_ID'] = np.arange(len(netFeedLine.pipe))
        netFeedLine.pipe['trail_ID'] = netFeedLine.pipe['trail_ID'].astype(int)

        # Re-set all pipes in service
        netFeedLine.pipe['in_service'] = True

        # Connect buildings to network junctions
        if buildings is not None:
            buildings_out = assignJunctionsToBuildings(buildings = buildings, junctions = gp.GeoDataFrame(netFeedLine.junction).set_crs(buildings.crs), junctionAttribute = 'junction')
            heat_consumers = buildings_out.copy()
            heat_consumers = heat_consumers[~heat_consumers['junction'].isna()].reset_index(drop = True)
            heat_consumers.rename(columns = {'junction':'from_junction'}, inplace = True)            

            heat_consumers = assign_default_values_ppi(
                gdf = heat_consumers,
                type = 'heat_consumer',
                tidy_up = False
            )

        else:
            buildings_out = None
            heat_consumers = None

        # Create pandapipes network for combined feed and reflux line
        # Attribute names in dict_attributes.values() are optional attrbutes whcih shall be added to the ppi network component table if they
        # exist in the provided component DataFrame.
        netFeedRefluxLine = create_ppi_network(
            network_name = networkName,
            fluid = networkFluid,
            modelling_type = modelling_type,
            consumerFeedLineJunctionIdx = list(heat_consumers['from_junction']),
            junctions=netFeedLine.junction.copy(),
            pipes=netFeedLine.pipe.copy(),
            valves=netFeedLine.valve.copy() if hasattr(netFeedLine, 'valve') else valves,
            heat_consumers = heat_consumers,
            dict_attributes = {
                'pipes':['connectionType', 'nominalWidth', 'Layer', 'ID', 'trail_ID', 'geometry'],
                'junctions':['Layer', 'ID', 'geometry'],
                'heat_consumers':[buildings_uniqueID, heatingDemandAttr, 'geometry', 'name'],
                'valves':['Layer', 'geometry']
            }
        )

        if producers is not None:
            netFeedRefluxLine = create_ppi_producers_from_GIS_data(
                net = netFeedRefluxLine,
                production_sites = producers,
                modelling_type = modelling_type,
                dict_col_layer = {'Layer':{'feed':'feedLine', 'reflux':'refluxLine'}},
                dict_col_type = producerTypeDict,
                add_cols = producerAddAttr
            )

        net = netFeedRefluxLine.deepcopy()

    
    return net