# -*- coding: utf-8 -*-
import geopandas as gp
import pandas as pd
import numpy as np
import rasterio as rio
import networkx as nx

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

    return (idx)

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

    # Find minimum value (e.g. outside house ground surface) within buffer distance around point
    zs = zonal_stats(j_buff, arr, affine = affine, stats=statsType)
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
    import momepy
    from shapely import set_precision

    # Create temporary copies and set precision
    lines = lines.copy()
    lines['geometry'] = set_precision(lines['geometry'], 0.01)

    G = momepy.gdf_to_nx(lines, approach="primal")
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
        edges = momepy.nx_to_gdf(sub, lines=True, points = False)
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


    # difference = dp_per_l * d**5 * rho * np.pi**2 / (mdot**2 * 8) - (0.25 / (np.log10(k/1e3/(3.71*d) + (5.74/(4*mdot/(d*np.pi*nu*rho))**(0.9))))**2)

    difference = dp_per_l - rho / (2*d) * (mdot * 4 / (rho*d**2*np.pi))**2 * 0.25 / (np.log10(k/1e3/(3.71*d) + (5.74/(4*mdot/(d*np.pi*nu*rho))**(0.9))))**2

    return difference

### *--- Functions not meant for publication ---*
def plot_ppi_network_v3(
        net, # pandapipes network
        plab:bool = False,
        jlab:bool = False,
        lengthlab:bool = False,
        figsize:tuple = (16,10),
        returnAxes:bool = False,
        **kwargs
    ):

    """
    Function that controls simple plot as inherent plotting function of pandapipes.

    Exemplary keyword arguments are:
    respect_valves=False, 
    respect_in_service=True, 
    pipe_width=2.0,
    junction_size=1.0, 
    ext_grid_size=1.0, 
    plot_sinks=False, 
    plot_sources=False,
    sink_size=1.0, 
    source_size=1.0, 
    valve_size=1.0, 
    pump_size=1.0,
    heat_exchanger_size=1.0, 
    pressure_control_size=1.0, 
    compressor_size=1.0, 
    flow_control_size=1.0,
    scale_size=True, 
    junction_color="r", 
    pipe_color='silver', 
    ext_grid_color='orange',
    valve_color='silver', 
    pump_color='silver', 
    heat_exchanger_color='silver',
    pressure_control_color='silver', 
    compressor_color='silver', 
    flow_control_color='silver',
    library="igraph", 
    show_plot=True

    """

    ### Imports
    import pandapipes.plotting as ppplot
    import matplotlib.pyplot as plt
    import numpy as np

    ### Default kwargs
    defaultValues = {
        'respect_valves':False, 
        'respect_in_service':True, 
        'pipe_width':2.0,
        'junction_size':0.25, 
        'ext_grid_size':0.5, 
        'plot_sinks':False, 
        'plot_sources':False,
        'sink_size':0.5, 
        'source_size':0.5, 
        'valve_size':0.25, 
        'pump_size':0.25,
        'plot_heat_consumer':False,
        'heat_consumer_size':10,
        'heat_exchanger_size':0.5, 
        'pressure_control_size':0.5, 
        'compressor_size':0.5, 
        'flow_control_size':0.5,
        'scale_size':True, 
        'junction_color':"r", 
        'pipe_color':'blue', 
        'ext_grid_color':'orange',
        'valve_color':'silver', 
        'pump_color':'silver', 
        'heat_exchanger_color':'black',
        'pressure_control_color':'silver', 
        'compressor_color':'silver', 
        'flow_control_color':'black',
        'heat_consumer_color':'black',
        'library':"igraph", 
        'show_plot':True
    }

    for k, v in defaultValues.items():
        if k not in kwargs.keys():
            kwargs[k] = v


    ### Create figure and axes environment
    fig, ax = plt.subplots(figsize = figsize)

    ### Create simple plot
    ax = ppplot.simple_plot(
        net = net, 
        ax = ax,
        respect_valves=kwargs['respect_valves'], 
        respect_in_service=kwargs['respect_in_service'], 
        pipe_width=kwargs['pipe_width'],
        junction_size=kwargs['junction_size'], 
        ext_grid_size=kwargs['ext_grid_size'], 
        plot_sinks=kwargs['plot_sinks'], 
        plot_sources=kwargs['plot_sources'],
        sink_size=kwargs['sink_size'], 
        source_size=kwargs['source_size'], 
        valve_size=kwargs['valve_size'], 
        pump_size=kwargs['pump_size'],
        plot_heat_consumer = kwargs['plot_heat_consumer'],
        heat_exchanger_size=kwargs['heat_exchanger_size'], 
        pressure_control_size=kwargs['pressure_control_size'], 
        compressor_size=kwargs['compressor_size'], 
        flow_control_size=kwargs['flow_control_size'],
        scale_size=kwargs['scale_size'], 
        junction_color=kwargs['junction_color'], 
        pipe_color=kwargs['pipe_color'], 
        ext_grid_color=kwargs['ext_grid_color'],
        valve_color=kwargs['valve_color'], 
        pump_color=kwargs['pump_color'], 
        heat_exchanger_color=kwargs['heat_exchanger_color'],
        pressure_control_color=kwargs['pressure_control_color'], 
        compressor_color=kwargs['compressor_color'], 
        flow_control_color=kwargs['flow_control_color'],
        library=kwargs['library'], 
        show_plot=False
        )

    
    ### Control plotting of heat_consumer
    if kwargs['plot_heat_consumer']:
        if hasattr(net, 'heat_consumer'):
            js = net.junction[net.junction.index.isin(np.concatenate((net.heat_consumer['from_junction'],net.heat_consumer['to_junction'])))].index
            x = [net.junction_geodata['x'][id] for id in js]
            y = [net.junction_geodata['y'][id] for id in js]
            ax.scatter(
                x = x,
                y = y,
                marker  = 'o',
                c = kwargs['heat_consumer_color'],
                s = kwargs['heat_consumer_size'],
                zorder = 10
            )


    ### Control behaviour of additional annotations to junction and pipes
    if jlab:    
        offset = 0.05
    
        for n, id in enumerate(net.junction_geodata.index):        
            ax.text(net.junction_geodata.x[id] + 2 * offset,
                            net.junction_geodata.y[id] + 2 * offset,
                            "J " +
                            str(net.junction_geodata.index[id])
                            
                            )
        
    if plab:
        offset = 0.05

        for n, id in enumerate(net.pipe.index):
            x_mid = 0.5*(net.junction_geodata.x[net.pipe.from_junction[id]]+net.junction_geodata.x[net.pipe.to_junction[id]])
            y_mid = 0.5*(net.junction_geodata.y[net.pipe.from_junction[id]]+net.junction_geodata.y[net.pipe.to_junction[id]])
            ax.text(x_mid + 2 * offset,
                    y_mid + 2 * offset,
                    "P " + 
                    str(net.pipe_geodata.index[n])
                    )

    if lengthlab:
        offset = 0.05

        for n, id in enumerate(net.pipe.index):
            x_mid = 0.5*(net.junction_geodata.x[net.pipe.from_junction[id]]+net.junction_geodata.x[net.pipe.to_junction[id]])
            y_mid = 0.5*(net.junction_geodata.y[net.pipe.from_junction[id]]+net.junction_geodata.y[net.pipe.to_junction[id]])
            ax.text(x_mid + 1 * offset,
                    y_mid + 1 * offset,
                    str(np.round(net.pipe.length_km[id], 3))
                    )
    
    if kwargs['show_plot']:
        plt.show()


    if returnAxes:
        return ax
