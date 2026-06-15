# %%
# -*- coding: utf-8 -*-

import geopandas as gp
import pandas as pd
import numpy as np
import networkx as nx

from geoppi.auxFunctions import (closest_point, sort_with_permutation, gdf_to_nx, nx_to_gdf, )

### Functions for analysis of shortest path lengths between specified sets of start- and end nodes and summed weights
def shortest_paths_pathweights(
        G:nx.graph,
        startNodes:list = None,
        endNodes:list = None,
        weight:str = None
    ):
    """
    Function to return shortest paths as lists of visited nodes and corresponding path lengths in networkx MultiGraph objects.\n

    :param G: networkx graph object (possible: MultiGraph)\n
    :param startNodes: list-like denoting indicators of start nodes in graph for path definition, defaults to None (all nodesi ngraph)\n
    :param endNodes: list-like denoting spcified end nodes in graph for path definition, defaults to None (all nodes in graph)\n
    :param weight: str denoting edge attribute for weights of path lengths, defaults to None\n
    :return: two dictionaries with indicators of end nodes as keys and 1) paths (list of visited nodes) as values and 2) path lengths
    """

    # Plausbility checks
    if startNodes is None:
        startNodes = list(G.nodes())

    if endNodes is None:
        endNodes = list(G.nodes())

    # Initializations
    paths = []
    lengths = []

    # Extract all paths from start nodes to specified end nodes
    ## Store dictionaries for all start nodes in list
    for n in startNodes:
        path = dict(nx.single_source_all_shortest_paths(G = G, source = n, weight = weight))
        length = dict(nx.shortest_path_length(G = G, source = n, weight = weight))

        path_new = {k:v for k,v in path.items() if k in endNodes}
        length_new = {k:v for k,v in length.items() if k in endNodes}

        # Cut paths to those ending at the end nodes
        paths.append(path_new)
        lengths.append(length_new)

    # Initialize dictionaries for shotrets paths and lengths
    alleEndNodes = set().union(*(d.keys() for d in lengths))
    shortestLength = {}
    shortestPath = {}

    for key in alleEndNodes:
        min_len = float('inf')
        min_idx = -1

        for i, d in enumerate(lengths):
            if key in d and d[key] < min_len:
                min_len = d[key]
                min_idx = i
            if min_idx >= 0:
                shortestLength[key] = min_len
                shortestPath[key] = paths[min_idx][key]

    return shortestPath, shortestLength


def sum_attrs_along_shortest_paths(
        G:nx.graph,
        startNodes:list = None,
        endNodes:list = None,
        weight:str = None,
        output_attr:str = None,
        dictAttrsEndNodes:dict = None
    ):
    """
    Function that sums up attributes specified at end nodes on all edges belonging to the shoretst path between end nodes and specified start nodes. The attributes at the end nodes are stored in the dictionary **dictAttrsEndNodes** with node indicators as keys and attributes as values.\n

    :param G: networkx graph object (possible: MultiGraph).\n
    :param startNodes: list-like denoting indicators of start nodes in graph for path definition, defaults to None (all nodesi ngraph)\n
    :param endNodes: list-like denoting spcified end nodes in graph for path definition, defaults to None (all nodes in graph)\n
    :param weight: str denoting edge attribute for weights of path lengths, defaults to None\n
    :param output_attr: str denoting the desired output attribute name on graph edges with summed attributes, defaults to None\n
    :param dictAttrsEndNodes: dict containing end node indicators as keys and attributes as values. These attr8ibutes are summed along the edges of the shortest paths between the end nodes and the start nodes., defaults to None\n
    :return: networkx graph object G with newly created edge attributes and dict containing edge attributes (values of dictionary) for edges with unique key (keys of dictionary).\n
    """

    # Initializations
    outAttr_weight = 'sum_' if output_attr is None else output_attr

    # Check if G is multigraph object and may contain edge keys
    multi = True if isinstance(G,  (nx.MultiGraph, nx.MultiDiGraph)) else False

    # Calculate shortest path between each end node and the closest starting node and respective shortest path weights
    shortestPath, _ = shortest_paths_pathweights(G = G, startNodes = startNodes, endNodes = endNodes, weight = weight)

    if multi:  # mulitgraph object, edge keys for parallel edges may be present
        for k, path in shortestPath.items():
            for u, v in zip(path[0][:-1], path[0][1:]):

                edges = G[u][v] # -> All edges betwen nodes u and v

                # Address possible parallel edges
                min_key, _ = min(edges.items(), key=lambda item: item[1].get(weight, float('inf')))

                G[u][v][min_key][outAttr_weight] += dictAttrsEndNodes[k]

        # Initialize final output dictionary for matching of line ID and summed Attribute (summed weight)
        edge_attr_dict = {}
        for u, v, k, attr in G.edges(keys=True, data=True):
            edge_attr_dict[k] = attr.get(outAttr_weight)

    else:
        for k, path in shortestPath.items():
            for u, v in zip(path[0][:-1], path[0][1:]):

                edges = G[u][v] # -> All edges betwen nodes u and v

                G[u][v][outAttr_weight] += dictAttrsEndNodes[k]

        # Initialize final output dictionary for matching of line ID and summed Attribute (summed weight)
        edge_attr_dict = {}
        for u, v, attr in G.edges(data=True):
            edge_attr_dict[(u,v)] = attr.get(outAttr_weight)

    return G, edge_attr_dict


def sum_heat_demands_to_closest_supplier(
        edgelist:pd.DataFrame,
        sources:str = 'from_junction',
        targets:str = 'to_junction',
        startNodes:list = None,
        endNodes:list = None,
        edge_key:str = 'unique_ID',
        weight:str = None,
        outAttr_weight:str = None,
        dictAttrsEndNodes:dict = None
    ):

    """
    Function for summing up attributes (e.g. heat demands) at specified nodes found in **dictAttrsEndNodes** along the shortest path between all end nodes and all start nodes.\n

    :param edgelist: pandas.DataFrame denoting edgelist for conversion into networkx graph object.\n
    :param sources: str denoting the column containing indices of source nodes for each edge.\n
    :param targets: str denoting the column containing indices of source nodes for each edge.\n
    :param startNodes: list containing indices of all start nodes from which shortest paths shall be searched.\n
    :param endNodes: list containing indices of all end nodes for which shortest paths shall be searched.\n
    :param edge_key: str denoting column name with uniqwie identifier for each edge in edgelist.\n
    :param weight: str denoting column name and edge attribute defining the edges' weights.\n
    :param outAttr: str denoting the output column name of the summed attributes along the shortest paths.\n
    :param dictAttrsEndNodes: dict containing end node indices as keys and attributes at these nodes as values.\n

    :returns: pandas.DataFrame with modified edgelist and summed attributes along shortest paths.\n
    """

    # Plausibility checks
    if edge_key not in edgelist.columns:
        edgelist[edge_key] = np.arange(len(edgelist))
        print(f'\n### Parameter edge_key as unique identifier for edges is not found. {edge_key} is added to the edgelist. ###')

    elif edgelist[edge_key].nunique() != len(edgelist):
        edgelist[edge_key] = np.arange(len(edgelist))
        print(f'\n### Provided parameter edge_key as unique identifier for edges is not unique!. {edge_key} is added to the edgelist. ###')

    # Create Multigraph from provided edgelist
    G = nx.from_pandas_edgelist(df = edgelist, source = sources, target = targets, edge_attr = weight, edge_key = edge_key, create_using = nx.MultiGraph())

    # Initialize summed weight attribute
    outAttr_weight = 'sum_' if outAttr_weight is None else outAttr_weight
    nx.set_edge_attributes(G, values = 0, name = outAttr_weight)

    # Sum weights along the sortest path between each end node and provided start nodes
    ## edge_attr_dict contains unique edge keys as keys and summed attributes on edges as values
    G, edge_attr_dict = sum_attrs_along_shortest_paths(G = G, startNodes = startNodes, endNodes = endNodes, weight = weight, output_attr = outAttr_weight, dictAttrsEndNodes = dictAttrsEndNodes)

    edgelist_out = edgelist.copy()

    edgelist_out[outAttr_weight] = edgelist_out.apply(lambda x: edge_attr_dict[x[edge_key]] if x[edge_key] in edge_attr_dict.keys() else 0, axis = 1)

    return edgelist_out


def count_downstream_buildings(
        lines:gp.GeoDataFrame,
        buildings:gp.GeoDataFrame,
        producers:gp.GeoDataFrame,
        output_attr:str = 'n_buildings_downstream',
        weight:str = 'length',
        line_id:str = 'unique_ID_lines',
        max_distance:float = 50.0,
    ):
    """
    Calculate the number of buildings downstream of each line segment in a radial network.

    The function uses the existing shortest-path utilities in GEOppi to create a producer-to-consumer
    orientation from the supplied line network. Producer locations are treated as sources and buildings
    as terminal consumers. Every line segment receives a cumulative count of buildings in its downstream
    subtree, including the buildings directly attached to that segment itself.

    :param lines: GeoDataFrame of line objects representing the distribution network.
    :param buildings: GeoDataFrame of building polygons to be counted downstream of each line.
    :param producers: GeoDataFrame of producer/starting-point geometries used as graph sources.
    :param output_attr: str denoting the output column name for downstream building counts.
    :param weight: str denoting the edge attribute used for shortest-path weighting.
    :param line_id: str denoting the unique line identifier in *lines*.
    :param max_distance: float denoting the search radius for attaching buildings to the nearest line.

    :returns: GeoDataFrame of line objects with the added downstream-building count attribute.
    """

    if lines.crs != buildings.crs:
        buildings = buildings.to_crs(lines.crs)

    if lines.crs != producers.crs:
        producers = producers.to_crs(lines.crs)

    if line_id not in lines.columns:
        lines = lines.copy()
        lines[line_id] = np.arange(len(lines), dtype=int)

    # Attach each building to its nearest line segment.
    buildings_attached = gp.sjoin_nearest(
        buildings[["geometry"]].copy(),
        lines[[line_id, "geometry"]].copy(),
        how="left",
        distance_col="distance_to_line",
    )

    buildings_attached = buildings_attached.rename(columns={line_id: line_id})
    buildings_attached = buildings_attached[buildings_attached[line_id].notna()].copy()

    buildings_per_line = (
        buildings_attached.groupby(line_id, dropna=False)
        .size()
        .rename("n_buildings_attached")
        .astype(int)
    )

    if buildings_per_line.empty:
        lines = lines.copy()
        lines[output_attr] = 0
        return lines

    # Build a simple primal graph from the line GeoDataFrame.
    G = gdf_to_nx(lines, approach="primal", multigraph=False, directed=False)

    # Use the existing shortest-path machinery to orient the network from producers to consumers.
    producer_coords = np.array([tuple(point) for point in zip(producers.geometry.x, producers.geometry.y)])
    graph_nodes = np.asarray(G.nodes)

    producer_nodes = []
    for point in producer_coords:
        nearest = closest_point(points=graph_nodes, target=point, threshDistance=max_distance)
        if len(nearest) > 0:
            producer_nodes.extend([tuple(node) for node in nearest])

    producer_nodes = list(dict.fromkeys(producer_nodes))

    if not producer_nodes:
        lines = lines.copy()
        lines[output_attr] = 0
        return lines

    path_map = nx.multi_source_dijkstra_path(G, sources=producer_nodes, weight=weight)

    # Create a directed production-to-consumption tree.
    D = nx.DiGraph()
    edge_line_map = {}

    for u, v, data in G.edges(data=True):
        edge_line_map[(u, v)] = data.get(line_id, None)

    for target, path in path_map.items():
        if len(path) < 2:
            continue
        parent = path[-2]
        child = path[-1]
        D.add_edge(parent, child)
        if (parent, child) in edge_line_map:
            D[parent][child][line_id] = edge_line_map[(parent, child)]

    # Propagate consumer counts from leaves towards the producer side.
    downstream_counts = {node: 0 for node in D.nodes}

    for node in reversed(list(nx.topological_sort(D))):
        edge_count = 0
        for child in D.successors(node):
            edge_count += downstream_counts[child]

        line_id_here = None
        for parent, child, data in D.in_edges(node, data=True):
            line_id_here = data.get(line_id, None)
            break

        count_here = int(buildings_per_line.get(line_id_here, 0)) if line_id_here is not None else 0
        downstream_counts[node] = count_here + edge_count

    lines_out = lines.copy()
    lines_out[output_attr] = 0

    for parent, child, data in D.edges(data=True):
        line_id_here = data.get(line_id)
        if line_id_here is not None:
            lines_out.loc[lines_out[line_id] == line_id_here, output_attr] = int(downstream_counts.get(child, 0))

    # Ensure every line has a sensible default value.
    lines_out[output_attr] = lines_out[output_attr].fillna(0).astype(int)

    return lines_out

### BFS ###
def network_span_bfs(
        lines:gp.GeoDataFrame,
        startpoints:gp.GeoDataFrame = None,
        val_Attr:str = 'ld',
        att_pth:str = 'Pth',
        length_Attr:str = 'length',
        att_nB:str = None,
        min_nB:int = np.inf,
        att_heatDemand:str = 'demand_use_th',
        sortByAttr:str = None,
        sortMethod:str = 'ascending',
        considerPowerBudget:bool = True,
        energyBudget:list = [np.inf],
        powerBudget:list = [np.inf],
        adaptThermalLoss = None,
        adaptSimultaneityFactor = None,
        adaptThermalPowerLoss = None,
        minvalAttr_maxLength:dict = {-1:1000},
        createMST:bool = False,
        weightMST:str = None,
        subsetNodesMST:gp.GeoDataFrame = None,
        removeDuplicateGraphs:bool = True,
        removeSmallNetworks:bool = True,
        mergeTouchingGraphs:bool = True,
        returnAddResults:bool = False,
        reduceInputGraphMultipleStartpoints:bool = False
    ):


    """
    Function to create spanned tree/network based on a defined energy budget that can be distributed from a chosen starting point.
    Breadth-first-search-based algorithm.

    :param lines: geopandas.GeoDataFrame of LineString objects representing the potential network routing.\n
    :param startpoints: geopandas.GeoDataFrame of Point objects used as starting points for the network routing algorithm\n
    :param val_Attr: string, edge attribute containing length-specific information (e.g. line density).\n
    :param att_pth: string, edge attribute containting information on summed heating power on line segments (MW).\n
    :param length_Attr: string, edge attribute containing information on line length (m)\n
    :param att_nB: string, denotes the edge attribute indicating th number of buildings connected to the line.\n
    :param min_nB: int, denoting the threshhold value for minimum number of buildings which shall be connected to a network.\n
    :param energyBudget: List of floats defining the energy budget (MWh) which may not be exceeded by the sum of (valAttr*lengthAttr) over all selected lines in the graph. List-like for iteration over each provided start point.\n
    :param powerBudget: List of floats, defining the thermal power budget (MW) which may not be exceeded over all selected lines in the graph. List-like for iteration over each provided start point.\n
    :param adaptThermalLoss: Function-like, Switch if thermal losses relative to annual heat demand in the form of thermalLoss/heatDemand shall be adapted with line density. A function can be provided. If None, no thermal losses are considered.\n
    :param adaptSimultaneityFactor: Function-like, Switch if varying simultaneity factor for summed thermal power demands shall be applied. A function of number of connected buildings in attribute **att_nB** can be provided. If None, the factor is set to 1.\n
    :param att_heatDemand: string, edge attribute containing information on summed heat demand in this line (MWh). By default it is None and the heat demand is calculated as a product of valAttr and lengthAttr.\n
    :param sortByAttr: string, denoting the attribute by which all edges (separately for each adjacent node) shall be sorted before starting the algorithm with sorting method **sortMethod**.\n
    :param sortMethod: string, either 'ascending' or 'descending' based on desired sorting method for the attribute **sortByAttr**.\n
    :param minvalAttr_maxLength: dict, Dictionary mapping the maximum allowed length of a line (m) to an upper threshold of **valAttr**\n
    :param createMST: Boolean if minimum spanning tree of all resulting graphs shall be calculated.\n
    :param weightMST: attribute name for calculation of minimum spanning tree.\n
    :param subsetNodesMST: geopandas.GeoDataFrame with subset of nodes (point objects) which shall be considered for caclulation of MST. If None or empty, MST for complete graph is calculated.\n
    :param removeDuplicateGraphs: Boolean switch if postprocessing function shall be applied which removes duplicate graphs from result list.\n
    :param removeSmallNetworks: boolean switch if postprocessing function shall be applied which removes graphs that are too small (refering to attribute att_nB and threshold min_nB).\n
    :param mergeTouchingGraphs: boolean switch if postprocessing function shall be applied which merges touchign graphs and removes duplicate edges.\n
    :param returnAddResults: boolean switch if additional results from the BFS algorithm shall be returned (SFList, usedEnergyBudgetList, usedPowerBudgetList, currentThermalLossFactorList, summedLengthList, avlineDensityList).\n
    :param reduceInputGraphMultipleStartpoints: boolean switch to control if subgraphs of already chosen edges shall be removed from the original graph before it is passed to the next starting location (for multiple starting locations only).\n
    
    :returns: edges of graph resulting from breadth-first-search algorithm
    
    Reminder for structure of networkx graphs edges:
    node = (x1,y1)
    G[(x1,y1)] yields adjacency view for node (x1,y1) in G to adjacent nodes (x2,y2), (x3,y3) ...
        * Adjacency view = {
            (x1,y1):{K1:{attributes:values}},
            (x2,y2):{K2:{attributes:values}},
            ...
            }

        with K1, K2, K3... being keys for the respective connection between the point. These are unique even in self-loops, allowing for correct identification of parallel lines.

    :return: _description_
    """

    # Initialize column for identification of lines
    enum_ID = 'enumeration_ID'
    lines[enum_ID] = np.arange(len(lines)) # Create unique ID for enumeration and identification of lines

    # Plausibility checks
    if any([expr not in lines.columns for expr in (val_Attr, att_pth, length_Attr, att_nB, att_heatDemand)]):
        print('\n### Provide all necessary columns for processing!. ###')
        return    
    
    # Create networkx graph from line representation
    if sortByAttr is not None:
        if sortMethod.lower() not in ('ascending', 'descending'):
            sortMethod = 'ascending'
            print('### No proper sort method for attribute ', sortByAttr, ' is chosen. Ascending is used as default method. ###')
        elif sortMethod.lower() == 'ascending':
            sortMethod = 'ascending'
        else:
            sortMethod = 'descending'

        G = gdf_to_nx(lines.sort_values(by = sortByAttr, ascending = True if sortMethod == 'ascending' else False), approach="primal")
    else:
        G = gdf_to_nx(lines, approach="primal") 


    # Defining startpoints
    if startpoints is None:
        print('### Please provide GeoDataFrame of startpoints for network extension search. Aborting... ###')
        return

    else:
        arrA = np.array([np.round([x,y],3) for x,y in zip(startpoints.geometry.x, startpoints.geometry.y)])
        arrB = np.array(G.nodes)

        startLocs = list(set((tuple(closest_point(points = arrB, target = poin, threshDistance = 100)),0) for poin in arrA))

    # Initialization of optional additional outputs
    list_G = []
    orderDictList = []
    SFList = []
    usedEnergyBudgetList = []
    usedPowerBudgetList = []
    currentThermalLossFactorList = []
    summedLengthList = []
    avlineDensityList = []

    # Initialization of input set of edges/original graph for BFS algo
    G_start = G.copy() 

    for nn, startLoc in enumerate(startLocs):
        print(f'\n### Network extension algorithm, startpoint {nn+1} of {len(startLocs)}.\n')
        # res contains tuples of form (edge), ID, order for each edge added to the final graph ((edge), enum_ID, order, currentSF, usedEnergyBudget, usedPowerBudget, currentThermalLossFactor, summedLength, avAtt)
        res  = list(
            bfs_network_extension(
                G = G_start,
                start = startLoc[0],
                valAttr = val_Attr,
                att_pth = att_pth,
                lengthAttr = length_Attr,
                att_nB = att_nB,
                att_heatDemand = att_heatDemand,
                considerPowerBudget = considerPowerBudget,
                energyBudget = energyBudget[nn],
                powerBudget = powerBudget[nn],
                adaptThermalLoss = adaptThermalLoss,
                adaptThermalPowerLoss = adaptThermalPowerLoss,
                adaptSimultaneityFactor = adaptSimultaneityFactor,
                minvalAttr_maxLength = minvalAttr_maxLength,
                enumeration_ID = enum_ID
                )
        )

        # Cut original graph to chosen list of edges to obtain network graph object
        G_add = cut_graph_to_edgelist(G = G, edgelist = [ed[0] for ed in res], keySensitive = True)

        # Create dictionary for order of addition (values) to graph for each line segment ID (keys)
        orderDictList.append({r[1]:r[2] for r in reversed(res)})

        # Gather additional results in lists
        SFList.append([r[3] for r in res])
        usedEnergyBudgetList.append([r[4] for r in res])
        usedPowerBudgetList.append([r[5] for r in res])
        currentThermalLossFactorList.append([r[6] for r in res])
        summedLengthList.append([r[7] for r in res])
        avlineDensityList.append([r[8] for r in res])

        # Add graph to output list of graphs
        list_G.append(G_add)

        # Cut G_start to remaining, unchosen edges from original graph G
        if reduceInputGraphMultipleStartpoints:
            G_start = cut_graph_to_edgelist(G = G_start, edgelist = [ed[0] for ed in res], keySensitive = True, reverseSelection = True)

    # list_G_unedited = list_G.copy()
    # print(f'no. of edges in unedited graph = {len(list(list_G_unedited[0].edges))}')
    # print(f'no. of edges in graph = {len(list(list_G[0].edges))}')

    # Deletion of graphs which are identical or subgraphs of each other
    # list_G is sorted in descending order from the biggest graph. If any subsequent graph is a complete subgraph (a duplicate), it will become empty. The permutation order of sorting is returned as well to keep consistency to the input order of starting locations.
    if removeDuplicateGraphs:
        list_G = remove_duplicate_graphs(graphList = list_G)

    # Merge touching graphs
    if mergeTouchingGraphs:
        list_G = merge_touching_graphs(graphList = list_G, G = G)

    # Remove too small networks
    if removeSmallNetworks:
        list_G = remove_small_networks(graphList = list_G, valAttrvalThresh = {att_nB:min_nB})

    # Calculate MST for graphs
    if createMST:
        arrA = np.array([np.round([x,y],3) for x,y in zip(subsetNodesMST.geometry.x, subsetNodesMST.geometry.y)])
        
        for n, G_b in enumerate(list_G):            
            arrB = np.array(G_b.nodes)
            MSTnodes = list(filter(None, set(tuple(closest_point(points = arrB, target = poin, threshDistance = 10)) for poin in arrA)))
            list_G[n] = MST_graph_subset(g = G_b, weight = weightMST, subsetNodes = MSTnodes)

    ### Convert graphs into GeoDataFrames and return them
    for n, G_b in enumerate(list_G):

        nodes, edges, sw = nx_to_gdf(G_b, points=True, lines=True, spatial_weights=True)
        edges = edges[[col for col in edges.columns if col in list(lines.columns) + ['order']]]

        # Assign additional columns
        edges['network_ID'] = n

        # Correct enumeration_ID for lines added in postprocessing
        orderDict = orderDictList[n]
        edges['order'] = [orderDict[x] if x in orderDict.keys() else 999999 for x in edges[enum_ID]]
        edges.drop(columns = enum_ID, inplace = True)

        if n == 0:
            edges_collector_gdf = edges.copy()
        else:
            edges_collector_gdf = pd.concat((edges_collector_gdf, edges))
    
    if returnAddResults:
        return edges_collector_gdf, list_G, SFList, usedEnergyBudgetList, usedPowerBudgetList, currentThermalLossFactorList, summedLengthList, avlineDensityList
    
    else:
        return edges_collector_gdf, list_G

def bfs_network_extension(
        G:nx.graph,
        start:tuple = None,
        valAttr:str = 'ld_demand_use_th',
        att_pth:str = None,
        lengthAttr:str = 'length',
        att_nB:str = None,
        att_heatDemand:str = None,
        energyBudget:float = 1000.0,
        powerBudget:float = np.inf,
        adaptThermalLoss = None,
        adaptSimultaneityFactor = None,
        adaptThermalPowerLoss:bool = None,
        considerPowerBudget:bool = True,
        minvalAttr_maxLength:dict = {1:100},
        enumeration_ID:str = None
    ):


    """
    Function to create spanned tree/network based on a defined energy budget that can be distributed from a chosen starting point.
    Breadth-first-search-based algorithm.

    :param G: networkx.graph object containing all edges and nodes that can possibly be part of the output graph and shall be searched.\n
    :param start: Tuple of (x,y) point values for starting point of search.\n
    :param val_Attr: Edge attribute containing length-specific information (e.g. line density).\n
    :param att_pth: Edge attribute containting information on summed heating power on line segments (MW).\n
    :param length_Attr: Edge attribute containing inforation on line length (m)\n
    :param att_nB: denotes the edge attribute indicating th number of buildings connected to the line
    :param energyBudget: Float defining the energy budget (MWh) which may not be exceeded by the sum of (valAttr*lengthAttr) over all selected lines in the graph\n
    :param powerBudget: Float defining the thermal power budget (MW) which may not be exceeded over all selected lines in the graph.\n
    :param adaptThermalLoss: Function-like or None. Function takes line density of entire current network as input and returns proportional value of summed annual heat demand in the current network. Switch if **thermal losses relative to annual heat demand** in the form of thermalLoss/heatDemand shall be adapted with line density. This factor is also applied to thermal power demand. A function can be provided. If None, no thermal losses are considered. 
    :param adaptSimultaneityFactor: Function-like. Switch if varying simultaneity factor for summed thermal power demands shall be applied. A function of number of connected buildings in attribute *att_nB* can be provided. If None, the factor is set to 1.
    :param adaptThermalPowerLoss: function-like. Function takes power density (heat demand power / trail length) (MW/m) and returns relatiove factor of thermal power loss.\n
    :param att_heatDemand: Edge attribute containing information on summed heat demand in this line (MWh). By default it is None and the heat demand is calculated as a product of valAttr and lengthAttr.\n
    :param minvalAttr_maxLength: Dictionary mapping the maximum allowed length of a line (m) to an upper threshold of **valAttr**\n
    :param unitEnergy: By default, energy budget is expressed as MWh, line density is expressed as MWh/m\n
    :param unitPower: By default, thermal power budget is expressed as MW\n

    :returns: edges of graph resulting from breadth-first-search algorithm
    
    Reminder for structure of networkx graphs edges:
    node = (x1,y1)
    G[(x1,y1)] yields adjacency view for node (x1,y1) in G to adjacent nodes (x2,y2), (x3,y3) ...
        * Adjacency view = {
            (x1,y1):{K1:{attributes:values}},
            (x2,y2):{K2:{attributes:values}},
            ...
            }

        with K1, K2, K3... being keys for the respective connection between the point. These are unique even in self-loops, allowing for correct identification of parallel lines.

        
    Author: C. Völzel, 2024-08
    """

    from collections import deque

    nodes = list(G.nbunch_iter(start)) # Create list of nodes, ordered with start node at the beginning
    if not nodes:
        return

    directed = G.is_directed()
    orientation = None
    kwds = {"data": False}

    if G.is_multigraph() is True:
        kwds["keys"] = True

    # Set up edge lookup
    def edges_from(node):
        return iter(G.edges(node, **kwds))

    if directed:
        def edge_id(edge):
            # remove direction indicator
            return edge[:-1] if orientation is not None else edge

    else:
        def edge_id(edge):
            return (frozenset(edge[:2]),) + edge[2:] # Extract unique id for each edge by accessing all additional attributes (data dictionary) besides (x,y) coordinates of start and end points
    
    # Validity checks
    if considerPowerBudget and powerBudget == np.inf:
        print('\n### considerPowerBudget is set to True, but provided powerBudget is np.inf. Please check again. Aborting... ###')
        return

    # Initialisations
    eBudget = energyBudget
    usedEnergyBudget = 0

    pBudget = powerBudget
    usedPowerBudget = 0

    order = 0


    # start BFS
    visited_nodes = set(nodes)

    visited_edges = set()
    queue = deque([(n, edges_from(n)) for n in nodes]) # Create list of collections with all nodes and adjacent edges in graph G

    potentialEdges = []

    # Plausibility checks
    if att_pth is None:
        print(f'"att_pth" is not provided in graph representation. It is calculated from {att_heatDemand}/1700. ###')



    while queue and (usedEnergyBudget <= eBudget) and (usedPowerBudget <= pBudget or not considerPowerBudget):

        parent, children_edges = queue.popleft() # Get parent node and children edges by queue.popleft()-command

        for edge in children_edges:
            child = edge[1]
           
            if child not in visited_nodes:
                visited_nodes.add(child)
                queue.append((child, edges_from(child)))

            edgeid = edge_id(edge)

            # Check if edge has not been visited yet (key-sensitive)
            # and check if current edge shares at least one common node with the set of potential edges (output graph) to ensure coherence of the output graph
            if (edgeid not in visited_edges) and (not any(potentialEdges) or any(set(edge[:-1]) & set(b for e in potentialEdges for b in e[:-1]))):

                ### Loop over each keys of child -> Necessary to identify self-loops            
                for key in list(G[parent][child].keys()):
                    
                    if key == edgeid[1]:
                        edge_att_ld        = G[parent][child].get(key, {}).get(valAttr)
                        edge_att_length    = G[parent][child].get(key, {}).get(lengthAttr)
                        edge_att_hd        = G[parent][child].get(key, {}).get(att_heatDemand) if att_heatDemand is not None else edge_att_ld*edge_att_length
                        edge_att_pth       = G[parent][child].get(key, {}).get(att_pth) if considerPowerBudget else 0
                        enum_ID            = G[parent][child].get(key, {}).get(enumeration_ID)

                        ### Calculate statistics in current partial graph
                        G_new               = cut_graph_to_edgelist(G = G, edgelist = potentialEdges, keySensitive = True)
                        summedLength        = G_new.size(weight = lengthAttr) # summed network length
                        summedHeatDemand    = G_new.size(weight = att_heatDemand) # summed heat demand (without losses)
                        summedB             = G_new.size(weight = att_nB) if att_nB is not None else G_new.size() # Number of edges if att_nB not existent
                        avAtt               = summedHeatDemand/max(0.001, summedLength) # Average line density (MWh/m)

                        if att_pth is None:
                            summedPth      = G_new.size(weight = att_heatDemand / 1700)                            
                        else:                            
                            summedPth      = G_new.size(weight = att_pth) # Summed thermal power demand (no losses)

                        currentThermalLossFactor = adaptThermalLoss(max(0, avAtt)) if adaptThermalLoss is not None else 0
                        summedHeatDemand *= (1+currentThermalLossFactor) # summed heat demand (incl. losses)

                        currentThermalPowerLossFactor = adaptThermalPowerLoss(summedPth / max(1, summedLength)) if adaptThermalPowerLoss is not None else 0
                        currentSF = adaptSimultaneityFactor(summedB) if adaptSimultaneityFactor is not None else 1

                        summedPth *= (currentSF * (1+currentThermalPowerLossFactor)) # summed thermal power demand (incl. losses)

                        # print('Current thermal loss factor = '+str(currentThermalLossFactor))
                        # print('Current simultaneity factor = '+str(currentSF))
                        # print('PTh = '+str(summedPth))

                        # Check if line density attribute is high enough for adding edge to graph
                        if edge_att_ld > max(minvalAttr_maxLength.keys()):

                            # Check if energy budget and power budget are not exhausted
                            if (eBudget - usedEnergyBudget) >= edge_att_hd * (1+currentThermalLossFactor) and ((pBudget - usedPowerBudget) >= currentSF * (1+currentThermalPowerLossFactor) * edge_att_pth or not considerPowerBudget): 

                                visited_edges.add(edgeid)
                                usedEnergyBudget = summedHeatDemand
                                usedPowerBudget = summedPth

                                potentialEdges.append(edge)

                                order += 1
                                yield edge, enum_ID, order, currentSF, usedEnergyBudget, usedPowerBudget, currentThermalLossFactor, summedLength, avAtt

                        elif edge_att_ld <= max(minvalAttr_maxLength.keys()):
                            if (eBudget - usedEnergyBudget) >= edge_att_hd * (1+currentThermalLossFactor) and ((pBudget - usedPowerBudget) >= currentSF * (1+currentThermalPowerLossFactor) * edge_att_pth or not considerPowerBudget): 

                                # Find first value in minvalAttr_maxLength which is higher than the line's line density
                                nextid = list(map(lambda i: i > edge_att_ld, list(minvalAttr_maxLength.keys()))).index(True)
                                selectedKey = list(minvalAttr_maxLength.keys())[nextid]

                                if edge_att_length <= minvalAttr_maxLength[selectedKey]:
                                    visited_edges.add(edgeid)
                                    usedEnergyBudget = summedHeatDemand
                                    usedPowerBudget = summedPth

                                    potentialEdges.append(edge)

                                    order += 1
                                    yield edge, enum_ID, order, currentSF, usedEnergyBudget, usedPowerBudget, currentThermalLossFactor, summedLength, avAtt
                        
                        else:
                            continue


    ### Calculate statistics in final graph
    G_new               = cut_graph_to_edgelist(G = G, edgelist = potentialEdges, keySensitive = True)
    summedLength        = G_new.size(weight = lengthAttr)
    summedHeatDemand    = G_new.size(weight = att_heatDemand)
    summedB             = G_new.size(weight = att_nB) if att_nB is not None else G_new.size() # Number of edges if att_nB not existent
    summedPth           = G_new.size(weight = att_heatDemand / 1700) if att_pth is None else G_new.size(weight = att_pth) 

    currentSF = adaptSimultaneityFactor(summedB) if adaptSimultaneityFactor is not None else 1
    currentThermalPowerLossFactor = adaptThermalPowerLoss(summedPth / max(1, summedLength)) if adaptThermalPowerLoss is not None else 0
    summedPth *= (currentSF * (1+currentThermalPowerLossFactor))

    currentThermalLossFactor = adaptThermalLoss(summedHeatDemand/max(1,summedLength)) if adaptThermalLoss is not None else 0
    summedHeatDemand *= (1+currentThermalLossFactor)


    print(f'\n### Remaining energy budget = {(eBudget - summedHeatDemand):.2f} MWh, \nNumber of visited nodes = {len(visited_nodes)}')
    print(f'\n### Remaining thermal power budget = {pBudget - summedPth:.2f} MW')
    print(f'\n### Total number of connected buildings = {summedB}')
    print(f'\n### Total network length = {summedLength:.2f} m')
    print(f'\n### Average line density in network = {summedHeatDemand/max(1e-10,summedLength):.2f} MWh/m')
    print(f'\n### Potential lines are:\n {potentialEdges}\nNumber of potential lines is {len(potentialEdges)}')
    print(f'\n### Number of visited edges is {len(visited_edges)}')

    return

def get_biggest_coherent_subgraph(
        graphList:list
        ):
    
    """
    Function that searches for the biggest subgraph for each graph in a list of graphs.
    Searches by greatest number of **connected components** and returns this biggest coherent subgraph.

    :param graphsList: List of networkx.graph objects
    :returns: (modified)

    Author: C. Völzel, 2024-08
    """

    graphs_out = graphList.copy()

    if not isinstance(graphs_out, list):
        graphs_out = [graphs_out]

    for n, gr in enumerate(graphs_out):
        graphs_out[n] = gr.subgraph(max(nx.connected_components(gr), key=len))

    return graphs_out


### DFS ###
def network_span_dfs_level_search(
        lines:gp.GeoDataFrame,
        startpoints:gp.GeoDataFrame = None,
        sortByAttr:str = None,
        sortMethod:str = 'descending',
        valAttr:str = 'ld',
        lengthAttr:str = None,
        nminAttr:dict = {'nBuildings':10},
        nLevels:int = 3,
        createMST:bool = False,
        weightMST:str = None,
        subsetNodesMST:gp.GeoDataFrame = None,
        removeDuplicateGraphs:bool = True,
        removeSmallNetworks:bool = True,
        complementSingleStreets:bool = True,
        mergeTouchingGraphs:bool = True,
        **parameters
    ):
    
    """
    Function that creates coherent GeoDataFrame of lines from defined inputs by depth-first-search.

    :param lines: geopandas.GeoDataFrame of LineString objects
    :param startpoints: geopandas.GeoDataFrame of Point objects used as starting points (optional)
    :param sortByAttr: str denoting the column name by which the DataFrame shall be sorted in order *sortMethod* before conversion to networkx graph (sorting of adjacent nodes per node).\n
    :param sortMethod: str denoting method for sorting, either *ascending* or *descending*.\n
    :param valAttr: Attribute in **lines** as length-specific values. This is the primary attribute for controlling the dfs!
    :param lengthAttr: Attribute in **lines** denoting the length of the lines (optional: replaced by lines.geometry.length if not provided)
    :param nminAttr: dictionary providing additional attribute names and threshold values. Eventually, final network graphs' attributes are summed up. If the sum is below the threshold, the network graph is marked as not useable.
    :param nLevels: int defining the allowed search level depth (1-3)
    :param createMST: Boolean if minimum spanning tree of all resulting graphs shall be calculated.
    :param weightMST: attribute name for calculation of minimum spanning tree.
    :param subsetNodesMST: geopandas.GeoDataFrame with subset of nodes (point objects) which shall be considered for caclulation of MST. If None or empty, MST for complete graph is calculated.
    :param parameters: kwargs for additional arguments controlling the behaviour of the search algorithm. Possible kwrds are:\n
        - valStart          :float, denotes min. value of **valAttr** for determining starting lines if **startpoints** is not provided\n
        - valFirstlvl       :float, denotes value for 1st level search for next suitable line in network extension\n
        - valMeanSecondlvl  :float, denotes mean value of lines in 1st and 2nd level search for next suitable line in network extension \n
        - valMeanThirdlvl   :float, denotes mean value of lines in 1st, 2nd and 3rd level search for next suitable line in network extension \n
        - valRemove         :float, denotes minimum value of **valAttr** below which lines are sorted in ascending order of their length. Lines with bigger length and low values are deleted first (as long as the network is not separated).\n
        - minLength         :float, denotes minimum length of lines for search in 2nd and 3rd level.\n
        - minValJump        :float, denotes minimum value in **valAttr** to jump across the lines to a next level\n
        - maxLengthJump     :float, denotes the maximum length of lines to jump to higher level search\n
    
    :returns: tuple of GeoDataFrames of resulting edges in graphs and list of networkx graphs
    """
    # Create temporary copy
    lines = lines.copy()

    # Extracting default values
    if 'valStart' not in parameters.keys():
        valStart = 4
    else:
        valStart = parameters['valStart']

    if 'valFirstlvl' not in parameters.keys():
        valFirstlvl = 2.5
    else:
        valFirstlvl = parameters['valFirstlvl']

    if 'valMeanSecondlvl' not in parameters.keys():
        valMeanSecondlvl = 2
    else:
        valMeanSecondlvl = parameters['valMeanSecondlvl']

    if 'valMeanThirdlvl' not in parameters.keys():
        valMeanThirdlvl = 1.5
    else:
        valMeanThirdlvl = parameters['valMeanThirdlvl']

    if 'valRemove' not in parameters.keys():
        valRemove = 1
    else:
        valRemove = parameters['valRemove']

    if 'minLength' not in parameters.keys():
        minLength = 50
    else:
        minLength = parameters['minLength']

    if 'minValJump' not in parameters.keys():
        minValJump = 0
    else:
        minValJump = parameters['minValJump']

    if 'maxLengthJump' not in parameters.keys():
        maxLengthJump = 1000
    else:
        maxLengthJump = parameters['maxLengthJump']

    # Initialize column for identification of lines
    enum_ID = 'enumeration_ID'
    lines[enum_ID] = np.arange(len(lines)) # Create unique ID for enumeration and identification of lines

    # Creating networkx graph from input of lines GeoDataFrame  
    if lengthAttr is None:
        lines['length_dfs']     = lines.geometry.length
        print('\n### Length attribute ', lengthAttr, ' is not included in provided DataFrame. It is now calculated from the lines geometry. ###')

    # Create networkx graph from line representation
    if sortByAttr is not None:
        if sortMethod.lower() not in ('ascending', 'descending'):
            sortMethod = 'ascending'
            print('### No proper sort method for attribute ', sortByAttr, ' is chosen. Ascending is used as default method. ###')
        elif sortMethod.lower() == 'ascending':
            sortMethod = 'ascending'
        else:
            sortMethod = 'descending'

        G = gdf_to_nx(lines.sort_values(by = sortByAttr, ascending = True if sortMethod == 'ascending' else False), approach="primal")
    else:
        G = gdf_to_nx(lines, approach="primal") 
    
    # Defining startpoints
    if startpoints is None:
        startLocs = [(u, v, k) for u, v, k, data in G.edges(keys=True, data=True) if data.get(valAttr, 0) >= valStart]

    else:
        arrA = np.array([np.round([x,y],3) for x,y in zip(startpoints.geometry.x, startpoints.geometry.y)])
        arrB = np.array(G.nodes)

        startLocs = list(set((tuple(closest_point(points = arrB, target = poin, threshDistance = 100)),0) for poin in arrA))

    if len(startLocs) == 0:
        print('\n### No suitable starting point could be found by the parameters provided. Aborting... ###')
        return lines
    
    # Start algorithm: graphs are created for each starting point
    # All graphs are stored in a list **list_G**

    list_G = []
    orderDictList = []

    for j, startLoc in enumerate(startLocs):
        # res contains tuples of form ((parent), (child), key, ID, order, level) for each edge added to the final multigraph
        res = list(
            dfs_level_search(
                G=G, 
                start=startLoc[0], 
                valAttr=valAttr, 
                lengthAttr=lengthAttr, 
                valFirstlvl=valFirstlvl,
                valMeanSecondlvl=valMeanSecondlvl,
                valMeanThirdlvl= valMeanThirdlvl,
                minLength=minLength, 
                minValJump=minValJump,
                maxLengthJump = maxLengthJump,
                nLevels = nLevels,
                enumeration_ID = enum_ID
                )
                )
        
        G_add = cut_graph_to_edgelist(G = G, edgelist = [ed[:-3] for ed in res])

        # Create dictionary for order of addition to graph for each line segment
        orderDictList.append({r[3]:(r[4], r[5]) for r in reversed(res)})

        list_G.append(G_add)

    list_G_original = list_G.copy()   
    
    # Remove duplicate graphs from the list of graphs
    if removeDuplicateGraphs:  
        list_G = remove_duplicate_graphs(graphList = list_G)

    # unerkannte Einzelstraßen hinzufügen, die 
    if complementSingleStreets:
        list_G = [add_adjacent_edges_from_G(G = G, subG = G_b, valAttrvalThresh = {valAttr:valFirstlvl}) for G_b in list_G]

    # Merge all graphs that touch / share common edge or node
    if mergeTouchingGraphs:
        list_G = merge_touching_graphs(graphList = list_G, G = G)

    # Remove edges with valAttr below threshold
    list_G = [remove_edges(
        G = G_b, 
        maxVal = valRemove, 
        minValRemove = valRemove,
        valAttr = valAttr, 
        sortAttr = lengthAttr,
        sortReverse = True) for G_b in list_G]

    # Remove small networks
    if removeSmallNetworks:
        list_G = remove_small_networks(graphList = list_G, valAttrvalThresh = nminAttr)
    

    # Calculate MST for graphs
    if createMST:
        if subsetNodesMST is not None:
            arrA = np.array([np.round([x,y],3) for x,y in zip(subsetNodesMST.geometry.x, subsetNodesMST.geometry.y)])
        
        for n, G_b in enumerate(list_G):
            if subsetNodesMST is None:
                arrA = np.array(G_b.nodes)        
            arrB = np.array(G_b.nodes)
            MSTnodes = list(filter(None, set(tuple(closest_point(points = arrB, target = poin, threshDistance = 1)) for poin in arrA)))
            list_G[n] = MST_graph_subset(g = G_b, weight = weightMST, subsetNodes = MSTnodes)

    # Convert graphs to geopandas DataFrame and return them
    for n, G_b in enumerate(list_G):

        nodes, edges, sw = nx_to_gdf(G_b, points=True, lines=True, spatial_weights=True)
        edges = edges[[col for col in edges.columns if col in list(lines.columns) + ['order', 'level']]]
        
        # Assign additional columns
        edges['network_ID'] = n

        # Correct enumeration_ID for lines added in postprocessing
        orderDict = orderDictList[n]
        edges[['order', 'level']] = [[orderDict[x][0], orderDict[x][1]] if x in orderDict.keys() else [999999, 999999] for x in edges[enum_ID]]
        edges.drop(columns = enum_ID, inplace = True)

        if n == 0:
            edges_collector_gdf = edges.copy()
        else:
            edges_collector_gdf = pd.concat((edges_collector_gdf, edges), axis = 0)


    return edges_collector_gdf, list_G


def dfs_level_search(
        G:nx.graph, 
        start:tuple=None, 
        valAttr:str = 'ld_demand_use_th', 
        lengthAttr:str = 'length',
        valFirstlvl:float = 0, 
        valMeanSecondlvl:float = 0, 
        valMeanThirdlvl:float = 0,
        minLength:float = 0, 
        minValJump:float = 0, 
        maxLengthJump:float = 9999, 
        nLevels:int = 3,
        enumeration_ID:str = None
        ):

    """
    Function that performs depth-first-search through provided networkx graph G based on parameter inputs (controls)

    :param G: networkx graph object
    :param start: tuple defining start point coordinates (or edge)
    :param valAttr: Attribute in **G** as length-specific values. This is the primary attribute for controlling the dfs!
    :param lengthAttr: Attribute in **G** denoting the length of the lines (optional: replaced by lines.geometry.length if not provided)
    :param nLevels: int defining the allowed search level depth (1-3)
    :param valFirstlvl       :float, denotes value for 1st level search for next suitable line in network extension
    :param valMeanSecondlvl  :float, denotes mean value of lines in 1st and 2nd level search for next suitable line in network extension 
    :param valMeanThirdlvl   :float, denotes mean value of lines in 1st, 2nd and 3rd level search for next suitable line in network extension 
    :param minLength         :float, denotes minimum length of lines for search in 2nd and 3rd level.
    :param minValJump        :float, denotes minimum value in **valAttr** to jump across the lines to a next level
    :param maxLengthJump     :float, denotes the maximum length of lines to jump to higher level search
    
    :returns: list of edges within the spanned graph from G
    """
        
    # Initialisation of start node if none is provided
    start = random.choice(list(G.edges(keys=True)))[0] if start is None else start

    # Initialization of visited as empty set to track visited egdes
    visited = set()
    
    # Initialize empty list of potential egdes for 2nd level search (one step ahead)
    potentials_one_step = []
    
    # Initialize empty list of potential egdes for 3rd level search (two steps ahead)
    potentials_two_steps = []
    
    # Plausibility check for max. 3 levels for search
    if nLevels not in [1,2,3]:
        print(f'\n### Error: nLevels must not exceed 3. Aborting... ###')
        return
    
    # Initialization of counter for addition to graph
    order = 0
       

    """
    1st level
    """
         
    # Initialize stack with tuple (start, iter(G[start])). iter(G[start]) returns all neighbouring nodes of start.
    stack = [(start, iter(G[start]))]

    level = 1
    
    # While loop as long as stack is not empty
    while stack:
        
        # parent == current node
        # child == node accessible by parent
        parent, children = stack[-1]
        
        # Search for all children adjacent to the parent
        for child in children:

            # Initialize boolean switch if child can be added to the stack as new parent node in next iterations
            child_append_stack = False

            # Loop through parent - child - keys (multigraph)
            for key in list(G[parent][child].keys()):
 
                # Check if edge is already visited
                if (parent, child, key) not in visited:     

                    # -> key is always analysed as well in addition to parent and child node to properly capture rings in the network!
                    # Reason is that edges routing between identical nodes can only be discriminated from each other by their keys (e.g. indices in the edgelist)
  
                    # Check if valAttr fulfills first criterion
                    if G[parent][child].get(key, {}).get(valAttr) > valFirstlvl:
                          
                        # Return edge (defined by nodes)                 
                        order += 1
                        enumID = G[parent][child].get(key, {}).get(enumeration_ID)
                        yield parent, child, key, enumID, order, level
                        
                        # Add parent, child, key to visited
                        visited.add((parent, child, key))

                        child_append_stack = True
                         
                    # Check if edge can be added to potential edges for next level search 
                    elif G[parent][child].get(key, {}).get(lengthAttr) <= maxLengthJump\
                     and G[parent][child].get(key, {}).get(valAttr) >= minValJump\
                     and len(list(iter(G[child]))) > 0:

                        # Add egde to list of potential edges for next level search
                        potentials_one_step.append((parent, child, key))
                        
            # Add child to stack if boolean switch is True     
            if child_append_stack:

                stack.append((child, iter(G[child])))
                break
        
        # Remove top element from stack  
        else:
            stack.pop()
 
    """
    2nd level
    """
    
    # Check for existing potential edges for 2nd level search
    if potentials_one_step and nLevels > 1:

        level = 2

        potentials_one_step = list(set(potentials_one_step))
    
        # Loop through edges contained in list of potential egdes
        for potential in potentials_one_step:
            
            # Extract parameters of potential edge
            # Naming parent_org because current child can also become parent (multigraph, self-loops)
            grandparent, parent_org, gp_p_key = potential
            
            # Initialize stack with tuple (start, iter(G[parent_org])). iter(G[parent_org]) returns all neighbouring nodes of parent_org.
            stack = [(parent_org, iter(G[parent_org]))]
             
            # While loop as long as stack is not empty
            while stack:
                
                # parent == current node
                # child == node accessible by parent
                parent, children = stack[-1]
                
                # Search for all children adjacent to the parent
                for child in children:
                    
                    # Initialize boolean switch if child can be added to the stack as new parent node in next iterations
                    child_append_stack = False
                
                    # Loop through parent - child - keys (multigraph)
                    for key in list(G[parent][child].keys()):

                        # Check if edge is already visited   
                        if (parent,child,key) not in visited or parent_org == parent:     

                            # Check if valAttr fulfills first criterion
                            if G[parent][child].get(key, {}).get(valAttr) > valFirstlvl:  

                                # Initialize boolean switch if egdge may be added
                                add_edge_okay = False

                                # Additional check for 2nd level criteria
                                if parent == parent_org:

                                    summedLength = G[grandparent][parent].get(gp_p_key, {}).get(lengthAttr) + G[parent][child].get(key, {}).get(lengthAttr)

                                    meanVal = (G[grandparent][parent].get(gp_p_key, {}).get(valAttr) * G[grandparent][parent].get(gp_p_key, {}).get(lengthAttr) + \
                                    G[parent][child].get(key, {}).get(valAttr) * G[parent][child].get(key, {}).get(lengthAttr)) / summedLength                           

                                    if meanVal > valMeanSecondlvl\
                                        and G[parent][child].get(key, {}).get(lengthAttr) > minLength:
                                            
                                            # Edges fulfill criteria
                                            add_edge_okay = True

                                else:
                                    if ((grandparent, parent, gp_p_key) not in visited) & (parent in G[grandparent]):
                                        add_edge_okay = True
                                        
                                # If edge is suited
                                if add_edge_okay:
                                                                           
                                    # Return edge gp - p (defined by nodes)    
                                    order += 1
                                    enumID = G[grandparent][parent].get(gp_p_key, {}).get(enumeration_ID)
                                    yield grandparent, parent, gp_p_key, enumID, order, level
                                    
                                    # Add grandparent, parent, key to visited
                                    visited.add((grandparent, parent, gp_p_key))
                                    
                                    # Return edge p - c (defined by nodes)
                                    order += 1
                                    enumID = G[parent][child].get(key, {}).get(enumeration_ID)
                                    yield parent, child, key, enumID, order, level
                                    
                                    # Add parent, child, key to visited
                                    visited.add((parent, child, key))
                                    
                                    child_append_stack = True
                                
                            else:
                                if parent_org == parent:
                                    summedLength = G[grandparent][parent].get(gp_p_key, {}).get(lengthAttr) + G[parent][child].get(key, {}).get(lengthAttr)
                                    meanVal = (G[grandparent][parent].get(gp_p_key, {}).get(valAttr) * G[grandparent][parent].get(gp_p_key, {}).get(lengthAttr) + \
                                    G[parent][child].get(key, {}).get(valAttr) * G[parent][child].get(key, {}).get(lengthAttr)) / summedLength

                                    if summedLength <= maxLengthJump and meanVal >= minValJump and len(list(iter(G[child]))) > 0:
                                 
                                        # Add connecting edge to list of potential edges for search in 3rd level
                                        potentials_two_steps.append((grandparent, parent, child, gp_p_key, key))
                                    
                    # Add child to stack if boolean switch is True      
                    if child_append_stack:
                        
                        # child zum Stack hinzufügen
                        stack.append((child, iter(G[child])))
                        break
                    
                # Remove top element from stack    
                else:
                    stack.pop()
                                
    """
    3rd level
    """
    
    # Check for existing potential edges for 3rd level search
    if potentials_two_steps  and nLevels > 2:
        level = 3
            
        # Loop through edges contained in list of potential egdes
        for potential in potentials_two_steps:
            
            # Extract parameters of potential edge
            # Naming parent_org because current child can also become parent (multigraph, self-loops)
            grandgrandparent, grandparent, parent_org, ggp_gp_key , gp_p_key = potential
            
            # Initialize stack with tuple (start, iter(G[parent_org])). iter(G[parent_org]) returns all neighbouring nodes of parent_org.
            stack = [(parent_org, iter(G[parent_org]))]
            
            # While loop as long as stack is not empty
            while stack:
                
                # parent == current node
                # child == node accessible by parent
                parent, children = stack[-1]
                
                # Search for all children adjacent to the parent
                for child in children:
                    
                    # Initialize boolean switch if child can be added to the stack as new parent node in next iterations
                    child_append_stack = False

                    # Loop through parent - child - keys (multigraph)
                    for key in list(G[parent][child].keys()):
                
                        # Check if edge is already visited   
                        if (parent,child, key) not in visited or parent == parent_org:

                            # Check if edge fulfills criterion
                            if G[parent][child].get(key, {}).get(valAttr) > valFirstlvl:

                                # Initialize boolean switch if egdge may be added
                                add_edge_okay = False
                                
                                # Additional check for 3rd level criteria
                                if parent == parent_org:

                                    summedLength = G[grandgrandparent][grandparent].get(ggp_gp_key, {}).get(lengthAttr) + G[grandparent][parent].get(gp_p_key, {}).get(lengthAttr) +\
                                        G[parent][child].get(key, {}).get(lengthAttr)
                                    
                                    meanVal = (G[grandgrandparent][grandparent].get(ggp_gp_key, {}).get(valAttr) * G[grandgrandparent][grandparent].get(ggp_gp_key, {}).get(lengthAttr) + \
                                                G[grandparent][parent].get(gp_p_key, {}).get(valAttr) * G[grandparent][parent].get(gp_p_key, {}).get(lengthAttr) + \
                                                G[parent][child].get(key, {}).get(valAttr) * G[parent][child].get(key, {}).get(lengthAttr)) / summedLength                          

                                    if meanVal > valMeanThirdlvl\
                                        and G[parent][child].get(key, {}).get(lengthAttr) > minLength:
                                    
                                        # Edges fulfill criteria
                                        add_edge_okay = True
                                      
                                else:
                                    if ((grandgrandparent,grandparent, ggp_gp_key) not in visited) and ((grandparent, parent, gp_p_key) not in visited) and (grandparent in G[grandgrandparent]) and (parent in G[grandparent]):
                                        add_edge_okay = True
                                    
                                # If edges iare suited   
                                if add_edge_okay:
                                        
                                    # Return edge ggp - gp (defined by nodes)
                                    order += 1
                                    enumID = G[grandgrandparent][grandparent].get(ggp_gp_key, {}).get(enumeration_ID)
                                    yield grandgrandparent, grandparent, ggp_gp_key, enumID, order, level
                                    
                                    # Add grandgrandparent, grandparent, key to visited
                                    visited.add((grandgrandparent, grandparent, ggp_gp_key))                                    
                                    
                                    # Return edge gp - p (defined by nodes)
                                    order += 1
                                    enumID = G[grandparent][parent].get(gp_p_key, {}).get(enumeration_ID)
                                    yield grandparent, parent, gp_p_key, enumID, order, level
                                    
                                    # Add grandparent, parent, key to visited
                                    visited.add((grandparent, parent, gp_p_key))
                                    
                                    # Return edge p - c (defined by nodes)
                                    order += 1                       
                                    enumID = G[parent][child].get(key, {}).get(enumeration_ID)     
                                    yield parent, child, key, enumID, order, level
                                    
                                    # Add parent, child, key to visited
                                    visited.add((parent, child, key))
                                    
                                    child_append_stack = True
                                
                    # Add child to stack if boolean switch is True       
                    if child_append_stack:

                        stack.append((child, iter(G[child])))
                        break
                    
                # Remove top element from stack    
                else:
                    stack.pop()


def cut_graph_to_edgelist(G:nx.graph,edgelist:list,keySensitive:bool = True, reverseSelection:bool = False):
    
    """
    Function that takes **G** as networkx graph and creates subgraph from **G** containing only edges provided in **edgelist**

    :param G: networkx graph object\n
    :param edgelist: list of edges from which subgraph of G shall be created\n
    :param keySensitive: boolean switch if check for identical keys shall be performed or only adjacent nodes are relevant.\n
    :param reverseSelection: boolean defining whether all edges not contained in edgelist shall be removed from original graph (reverseSelection == False) or if edges contained in edgelist shall be removed (reverseSelection == True)\n
    :returns: G_copy as subgraph of G\n

    """
    
    # 1) create copy of G
    # 2) all edges not contained in edgelist are removed from G 
    # 3) all nodes without a connection are removed

    # Create copy
    G_copy = G.copy()

    if not reverseSelection:
        # remove edges not contained in edgelist
        if keySensitive: # Used for MultiGraphs
            for edge in list(G_copy.edges(keys=True)):
                u, v, key = edge  # Extract nodes and keys

                if (u, v, key) not in edgelist and (v, u, key) not in edgelist:
                    G_copy.remove_edge(u, v, key)

        else:
            for edge in list(G_copy.edges(keys=False)):
                u, v = edge  # Extract nodes

                if (u, v) not in edgelist and (v, u) not in edgelist:
                    G_copy.remove_edge(u, v)


    else:
        # remove edges contained in edgelist
        if keySensitive: # Used for MultiGraphs
            for edge in list(G_copy.edges(keys=True)):
                u, v, key = edge  # Extract nodes and keys

                if (u, v, key) in edgelist:
                    G_copy.remove_edge(u, v, key)
                if (v, u, key) in edgelist:
                    G_copy.remove_edge(v, u, key)

        else:
            for edge in list(G_copy.edges(keys=False)):
                u, v = edge  # Extract nodes

                if (u, v) in edgelist:
                    G_copy.remove_edge(u, v)
                if (v, u) in edgelist:
                    G_copy.remove_edge(v, u)
             
    # Remove nodes without a connection
    isolated_nodes = [node for node, degree in dict(G_copy.degree()).items() if degree == 0]
    G_copy.remove_nodes_from(isolated_nodes)
    
    return G_copy

def remove_edges(
    G:nx.graph, 
    maxVal:float = 1, 
    minValRemove:float = 0, 
    valAttr:str = 'ld_demand_use_th',
    sortAttr:str = 'length',
    sortReverse:bool = True
    ):

    """
    Function that removes edges from provided networkx graph object based on defined parameters for attributes of the egdes.
    Edges are not removed if this separated the graph.

    On 1st level, edges are prioritised for removal if their value of **valAttr** is below **minValRemove**.
    On 2nd level, edges with values of **valAttr** between **minValRemove** and **maxVal** are checked for removal.

    :param G: networkx graph object
    :param maxVal: float defining the upper limit of values in **valAttr**. If values are higher than this value, the edge is not removed.
    :param minValRemove: float defining the lower limit of values in **valAttr** below which the edge may be removed.
    :param valAttr: attribute name by which to decide if edge shall be removed. Contains values that are compared to the input thresholds **maxVal** and **minValRemove**
    :param sortAttr: attribute name by which edges are sorted.
    :param sortReverse: boolean if edges with highest value in **sortAttr** should be prioritised for removal (True) or those with lowest value (False)

    :return: (modified) graph
    """
    
    # Plausibiity checks
    if minValRemove >= maxVal: 
        minValRemove = maxVal-0.00000000000000001
    
    # Create copy of original graph
    G_cop = G.copy()

    edge_to_removed_first = set()
    edge_to_removed_second = set()
    
    # Loop through all edges in graph
    for edge in sorted(G_cop.edges(keys=True), key=lambda edge: G_cop.get_edge_data(*edge).get(valAttr)):
        
        # Search for edges affiliated with current parent-child combination
        pa_ch_keys = list(iter(G_cop.get_edge_data(edge[0],edge[1]).keys()))
          
        # Loop through keys
        for ekey in pa_ch_keys:

            # Extract value of valAttr
            currVal = G_cop.get_edge_data(edge[0],edge[1]).get(ekey, {}).get(valAttr)
                 
            # Check thresholds
            if minValRemove <= currVal < maxVal:
                            
                # edge mit len unter Grenzwert zur Liste hinzufügen
                edge_to_removed_second.add((edge[0], edge[1], ekey, G_cop[edge[0]][edge[1]].get(ekey, {}).get(sortAttr)))

            elif currVal < minValRemove:

                # edge mit len unter Grenzwert zur Liste hinzufügen
                edge_to_removed_first.add((edge[0], edge[1], ekey, G_cop[edge[0]][edge[1]].get(ekey, {}).get(sortAttr)))
                
          
    # Sorting for fourth attribute (sortAttr)
    edge_to_removed_first_sort = sorted(list(edge_to_removed_first), key=lambda x: x[3], reverse = True if sortReverse else False)
     
    # Loop through edge in list of edges prioritised for removal; Remove edges in this order
    for edge in edge_to_removed_first_sort:
        
        # Remove edge and (if necessary) unconnected nodes
        G_cop.remove_edge(edge[0],edge[1], key=edge[2])
        G_cop.remove_nodes_from([node for node, degree in dict(G_cop.degree()).items() if degree == 0])
        
        # Check if graph is still coherent
        if not nx.is_connected(G_cop):
            
            # Edge data
            edge_data = G.get_edge_data(edge[0],edge[1], key=edge[2])
            
            # If graph is separated by removel, re-add the edge
            G_cop.add_edge(edge[0],edge[1], key=edge[2], **edge_data)
        
    # Loop through edge in list of edges prioritised for removal on 2nd level; Remove edges in this order
    for edge in edge_to_removed_second:
        
        # Remove edge and (if necessary) unconnected nodes
        G_cop.remove_edge(edge[0],edge[1], key=edge[2])
        G_cop.remove_nodes_from([node for node, degree in dict(G_cop.degree()).items() if degree == 0])
        
        # Check if graph is still coherent
        if not nx.is_connected(G_cop):#
        
            # Edge data
            edge_data = G.get_edge_data(edge[0],edge[1], key=edge[2])
            
            # If graph is separated by removel, re-add the edge
            G_cop.add_edge(edge[0],edge[1], key=edge[2], **edge_data)
            
    return G_cop

def remove_duplicate_graphs(graphList:list):

    """
    Function that removes all duplicate graphs from a list of networkx graphs **graphList**

    :param graphList: list-like of networkx graph objects
    :return: list of unique graphs
    """
    
    # Initialize list of unique graphs 
    unique_graphs = []
    unique_graph_indices = []
    
    # Initialize stack of sorted graphs from list
    # stack = sorted(graphList, key=lambda graph: len(graph.edges(keys=True, data = True)), reverse=True)
    perm, stack = sort_with_permutation(lst = graphList, key = lambda x: len(x[1].edges(keys=True, data = True)), reverse = True)
  
    


    while stack:
        
        # Current graph
        curr_graph = stack[0]
        original_idx = perm[len(unique_graphs)]        

        # Append to list of unique graphs
        unique_graphs.append(curr_graph)
        unique_graph_indices.append(original_idx)
        
        # other graphs to compare to
        list_comp_graphs = stack[1:] # Reverse list
        comp_indices = perm[len(unique_graphs):]
       
        # empty stack
        stack = []
        perm_rest = []
        
        # Loop through all other graphs to compare to
        for n, comp_graph in enumerate(list_comp_graphs):
            
            # Get edges of current graph and graph to compare to
            edges_cur_set = set(curr_graph.edges(keys=True))
            edges_com_set = set(comp_graph.edges(keys=True))
            
            merged_set = edges_cur_set | edges_com_set
            
            # If set of edges is equal, then identical graphs or subgraph
            if not len(merged_set) == len(edges_cur_set):
                stack.append(comp_graph)
                perm_rest.append(comp_indices[n])

        perm = perm[:len(unique_graphs)] + perm_rest

        
    restore_order = sorted(zip(unique_graph_indices, unique_graphs), key=lambda x: x[0])    
                
    return [g for i, g in restore_order]

def merge_touching_graphs(graphList:list, G:nx.graph):

    """
    Function that merges all networkx graph objects contained in the provided list of graphs **graphList** to a single graph.

    :param graphList: list-like of networkx graph objects.
    :param G: networkx graph defining outer boundaries / max-. extension for all contained graphs

    :return: list of merged touchign graph objects
    """

    A_0 = graph_connections(graphList)
    while len(A_0) > 1:
        match = False
        for k in range(len(A_0)):
            if A_0[k] != []:
                c = A_0[k][0]
                match = True
                break
        if match:
            res = list(set(set(graphList[k].edges(keys=True)) | set(graphList[c].edges(keys=True))))
            if k > c: k, c = c, k
            del graphList[k]
            del graphList[c - 1]
            graphList.append(cut_graph_to_edgelist(G, res))
        A_0_new = A_0.copy()
        A_0 = graph_connections(graphList)
        if A_0_new == A_0:
            break
    return graphList

def remove_small_networks(
    graphList:list,
    valAttrvalThresh:dict = None
    ):

    """
    Function that removes graphs from a list of graphs which do not meet specified conditions-
    The sum of attributes of all edges, specified in **valAttrvalThresh.keys()** have to be higher than the corresponding threshold values in **valAttrvalThres.values()**

    :param graphList: list-like of networkx graph objects
    :param valAttrvalThresh: dictionary defining pairs of attribute names and corresponding minimum sum values across the graph

    :return: (modified) list of graph objects
    """

    list_G_final = []
    for G_b in graphList:

        # Initialization
        keepGraph = []

        for k, v in valAttrvalThresh.items():
            sumVals = sum(G_b.get_edge_data(edge[0], edge[1]).get(edge[2], {}).get(k, 0) for edge in G_b.edges(keys=True))
            
            if sumVals >= v:
                keepGraph.append(True)
            else:
                keepGraph.append(False)
        
        if all(keepGraph):
            list_G_final.append(G_b)

    return list_G_final

def add_adjacent_edges_from_G(
        G:nx.graph, 
        subG:nx.graph, 
        valAttrvalThresh:dict = None
        ):
    
    """
    Function that adds edges from superior graph G to subgraph subG under specified conditions:
    - the edge has to touch (share common node) with subG
    - the edge must not be part of subG already
    - the edge attributes **valAttr.keys()** have to exceed a provided threshhold (optional)

    :param G: networkx graph object
    :parm subG: networkx graph object
    :param valAttrvalThresh: dictionary matching attribute names to their minimum threshold value for decision to add edge to **subG**

    :return: (modified) subgraph
    """

    # Loop across all edges in G
    for edge in G.edges(keys=True):

        # Loop across all edge attributes of the edge
        for key in list(G[edge[0]][edge[1]].keys()):

            # Initialise boolean for adding edge to subG
            addEdge = False

            # Check if edge shares common node with subG, but is not yet contained in subG
            if (subG.has_node(edge[0]) or subG.has_node(edge[1])) and not subG.has_edge(edge[0], edge[1], key):

                if valAttrvalThresh is not None:

                    # Initialise counter
                    measure = []            

                    # Check if attribute values exceed thresholds for addition
                    for dictKey, dictVal in valAttrvalThresh.items():
                        if G[edge[0]][edge[1]][key].get(dictKey, 0) > dictVal:   
                            measure.append(True)
                        else:
                            measure.append(False)
                
                    if all(measure):
                        addEdge = True

                else:
                    addEdge = True

                if addEdge:                
                    # edge data
                    edge_data = G.get_edge_data(edge[0], edge[1], key=key)
                    # add edge
                    subG.add_edge(edge[0], edge[1], key=key, **edge_data)

    return subG

def graph_connections(graphs:list):
    """
    Function that searches for connected graphs within a list/iterable containing networkx.graph objects

    :param graphs: list-like of networkx.graph objects
    :return: nested list-like object indicating connections between graphs in provided list of graphs
    """
    connections = {}
    for i, graph1 in enumerate(graphs):
        connections[i] = []
        for j, graph2 in enumerate(graphs):
            if i != j and are_connected(graph1, graph2):
                connections[i].append(j)
    return connections 

def dfs(graph:nx.graph, start:tuple, visited):
    """
    Function performing depth-first-search in networkx.graph object graph 

    :param graph: networkx.graph object
    :param start: tuple of coordinates of start node (x,y)
    :param visited: list of visited nodes [(x,y), (...)]
    """
    visited.add(start)
    for neighbor in graph[start]:
        if neighbor not in visited:
            dfs(graph, neighbor, visited)

def are_connected(graph1, graph2):
    """
    Function that checks if two provided networkx.graph objects are connected, e.g. share common node(s)

    :param graph1: networkx.graph object
    :param graph2: networkx.graph object
    :return: Boolean indicating whether provided graph ibjects are connected.
    """
    if not graph1 or not graph2:  # Check if any graph is empty.
        return False    
    visited = set()
    start_node = next(iter(graph1)) if graph1 else None  # Choose starting node from graph1
    if start_node is not None:
        dfs(graph1, start_node, visited)
    
    con = False
    # Check if node from graph2 is included in visited nodes of graph1
    for node in graph2:
        if node in visited:
            con = True
    return con


def MST_graph_subset(g:nx.graph, weight:str = None, subsetNodes:set = None):

    """
    Function that creates a minimum spanning tree from the input graph **g** with the edge weight **weight**.
    Optional input of subset of nodes **subsetNodes** which shall be considered in the MSt.

    :param g: networkx graph
    :param weight: string denoting the edge attribute by which to assess MST weights
    :param subsetNode: set of obligatory nodes from which the MST shall be extracted.

    :returns: minimum spanning tree of g with the defined inputs
    """

    # Plausibility checks and set boolean for usage of subsetNodes
    if type(subsetNodes) == set | type(subsetNodes) == list:
        if len(subsetNodes) <= 1:
            calcMSTEntireGraph = True
            print('\n### Provided subset of nodes in graph g contains one or less nodes. MST for entire set of nodes in graph is calculated. ###')

    if subsetNodes is None:
        calcMSTEntireGraph = True

    # Initialisations
    calcMSTEntireGraph = False
    edgeset = set()

    if calcMSTEntireGraph:
        gout = nx.minimum_spanning_tree(g, weight = weight)

    else:
        for target in subsetNodes:
            if subsetNodes[0] == target:
                continue

            # Calculate shortest path from first node in subsetNodes to all other nodes contained in the subset
            path = nx.shortest_path(g, source = subsetNodes[0], target=target, weight=weight)

            # Extract edges from paths from subetNode[0] to all other contained nodes in subsetNodes
            edges = zip(path[:-1], path[1:])
            edgeset.update((u, v) for u, v in edges)

        # Create temporary new graph from edges in calculated paths
        g1 = nx.from_edgelist(edgeset)

        # Calculate MST from temporary graph
        gout = nx.minimum_spanning_tree(g.subgraph(g1.nodes), weight = weight)

    return gout

def MST_gdf_subset(lines:gp.GeoDataFrame, weightMST:str = None, subsetNodesMST:gp.GeoDataFrame = None):

    """
    Function that creates minimum spanning tree from input geopandas.GeoDataFrame of line objects.

    :param lines: geopandas.GeoDataFrame of line objects with attribute **weightMST** in columns.\n
    :param weightMST: column name/attribute for weighting MST.\n
    :param subsetNodesMST: geopandas.GeoDataFrame of point objects containing all obligatory points which shall be included in the MST representation.\n
    :return: edited GeoDataFrame of line objectsas MST of input GeoDataFrame
    """

    # Create temporary copy
    gdf = lines.copy()

    # Create graph representation of line GeoDataFrame
    G = gdf_to_nx(gdf, approach="primal")

    if subsetNodesMST is not None:
        arrA = np.array([np.round([x,y],3) for x,y in zip(subsetNodesMST.geometry.x, subsetNodesMST.geometry.y)])              
        arrB = np.array(G.nodes)
        MSTnodes = list(filter(None, set(tuple(closest_point(points = arrB, target = poin, threshDistance = 10)) for poin in arrA)))

    else:
        MSTnodes = None

    Gout = MST_graph_subset(g = G, weight = weightMST, subsetNodes = MSTnodes)

    # Convert networkx graph to GeoDataFrame
    _, lines_out, _ = nx_to_gdf(Gout, points=True, lines=True, spatial_weights=True)
    lines_out = lines_out[[col for col in lines_out.columns if col in lines.columns]]

    return lines_out


# %% *--- Minimal example for networkRouting algorithm ---*

if __name__ == '__main__':

    # Attribute for line density as main evaluation attribute on line objects
    att_ld = 'ld_demand_use_th'#total_exp_ren_1.5_2045_climate_demo'

    # Attributes for length of street segment
    lengthAttr = 'length'

    # Attribute for number of buildings connected to street segment
    att_nB = 'nBuildsdemand_use_th'#total_exp_ren_1.5_2045_climate_demo'

    # Minimum value for automatic selection of start edges
    dem_val_start = 4

    # Minimum value for 1st level search
    valFirstlvl = 2.1

    # Mean line density for 2nd level search (second Level)
    valMeanSecondlvl = 2.5

    # Mean line density for 3rd level search (third Level)
    valMeanThirdlvl = 2

    # Minimum length of street segments considered in 2nd and 3rd level
    minLength = 0.1

    # Minimum length of street segments omitted
    minValJump = 0

    # Maximum length of street segments omitted
    maxLengthJump = 1000

    # Threshold for line density attribute below which algorithm removes street segments sorted with descending length
    dem_val_remove_edge = 1

    # Minimum number of buildings connected to entire networks
    n_build = 17

    print(f'\n###\nExecuting minimal example for network routing algorithm with the following parameters:\n\
          valFirstlvl = {valFirstlvl}\n \
          valMeanSecondlvl = {valMeanSecondlvl}\n \
          valMeanThirdlvl = {valMeanThirdlvl}\n \
          dem_val_remove_edge = {dem_val_remove_edge}\n \
          ###')

    from shapely.geometry import Point, LineString
    import matplotlib.pyplot as plt

    linegeos = [
        LineString([[0, 0], [0, 1]]),
        LineString([[0, 1], [0, 2]]),
        LineString([[0, 2], [0, 3]]),
        LineString([[0, 2], [1, 2]]),
        LineString([[2, 2], [2, 3]]),
        LineString([[0, 3], [1, 3]]),
        LineString([[0, 2], [-1, 2], [-1, 3], [0, 3]]),
        LineString([[1, 2], [2, 2]]),
        LineString([[1, 3], [2, 3]])
    ]

    lines = gp.GeoDataFrame(geometry = linegeos)
    lines['ld_demand_use_th'] = [4, 3, 2, 2, 3, 2, 3, 2, 2]
    lines['nBuildsdemand_use_th'] = 20

    lines['length'] = lines.geometry.length

    startpoints = gp.GeoDataFrame(geometry = [Point(0, 0)])

    print(f'\n### Before algorithm ...')
    fig, ax2 = plt.subplots()
    lines.plot(
        ax = ax2, 
        column = 'ld_demand_use_th',
        legend = True,
        cmap = 'hot',
        vmin = 0,
        vmax = 5
        )
    
    startpoints.plot(ax = ax2, color = 'red', zorder = 10)

    gdf_graphs, graphList = network_span_dfs_level_search(
        lines = lines,
        startpoints = startpoints,
        sortByAttr = att_ld,
        sortMethod = 'descending',
        valAttr = att_ld,
        lengthAttr = lengthAttr,
        nminAttr = {att_nB:n_build},
        nLevels = 3,
        createMST = False,
        valFirstlvl = valFirstlvl,
        valMeanSecondlvl = valMeanSecondlvl,
        valMeanThirdlvl = valMeanThirdlvl,
        minValJump = minValJump,
        maxLengthJump = maxLengthJump,
        minLength = minLength,
        valRemove = dem_val_remove_edge,
        valStart = dem_val_start
    )

    print(f'\n### After algorithm ...')

    fig, ax = plt.subplots()
    gdf_graphs.plot(
        ax = ax, 
        column = 'ld_demand_use_th',
        legend = True,
        cmap = 'hot',
        vmin = 0,
        vmax = 5
        )

    startpoints.plot(ax = ax, color = 'red', zorder = 10)




