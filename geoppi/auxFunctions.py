# -*- coding: utf-8 -*-
import geopandas as gp
import pandas as pd
import numpy as np
import math
import rasterio as rio
import networkx as nx
import libpysal
import warnings
from shapely.geometry import Point

def calc_thermalLoss_pipe(
        net
    ):

    """
    Function that calculates pipe-specific thermal loss power in pandapipes network model.\n

    :param net: pandapipes network model with existing res_pipe DataFrame (available after thermal pipeflow).\n
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

        net.res_pipe['Pthermal_W'] = qloss

    else:
        pass

    return net

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


def detect_lines_in_narrow_passages(
    lines:gp.GeoDataFrame,
    polygons:gp.GeoDataFrame,
    threshDistance:float = 10,
    distPointsCircumference:float = 1,
    nNeighbours:int = 10,
    col:str = 'narrowPassage' 
    )->gp.GeoDataFrame:

    """
    Function that marks line objects if they are positioned between polygons with smaller distance than defined thresh value.\n

    :param lines: geopandas.GeoDataFrame with line objects (e.g. streets, ...).\n
    :param polygons: geopandas.GeoDataFrame with polygon objects (e.g. buildings, ...).\n
    :param threshDistance: float denoting the defined minimum distance between polygon outer boundaries between which lines are marked.\n
    :param distPointsCircumference: float denoting the distance of points that are created along the boundaries of polygons (!Affects quality of the result!).\n
    :param nNeighbours: integer denoting the number of neighbouring points to search within the defined thresh distance (!Affects quality of the result!).\n
    :param col: string denoting the column name in which either False (line not within minimum diastance in passage) or True is entered for original line objects.\n

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
        if boundary.length == 0:
            return []
        
        numPoints = int(boundary.length // dist)
        return [(np.round(boundary.interpolate(ii * dist).x, 2), np.round(boundary.interpolate(ii * dist).y, 2)) for ii in range(numPoints+1)]
           
    ## Data preparation
    # Merge all touching/overlapping polygons to reduce number of processed polygons
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

        if lOthers > 0:
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

    # Transfer results to output DF of lines
    lines = lines.copy()
    lines = lines.sjoin(shortestLineObjects, how = 'left', predicate = 'intersects')
    idxs = lines[~lines['index_right'].isna()].drop_duplicates(subset = 'geometry').index

    lines[col] = False
    lines.loc[idxs, col] = True


    return lines


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
    :param points: shapley.MultiPoint object
    :return: shapely.LineString object of split line object
    """

    # Internal function definitions
    def split_line_by_point(line, point, tolerance: float = 1.0e-1):

        # Imports
        from shapely.ops import split, snap
        return split(snap(line, point, tolerance), point)

    result = split_line_by_point(line, points)
    return(result)

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
    Function that creates matching of polygons to points by intersection. If multiple points within the provided entity of points interect with a polygon, the first match is taken.
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
