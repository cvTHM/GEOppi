# %%
import geopandas as gp
import pandapipes as ppi
from pandapipes.control import run_control

from geoppi.internal_auxFunctions import implement_controllers

from geoppi.internal_auxFunctions import (extract_FluidProperties_ppi, transfer_LoadPoint_ppi)

cs = 'EPSG:25832'


tFeed = 273.15 + 10
tReflux = 273.15 + 2


# %% Implement control strategy for substation with HP and defined COP at heat consumers

#### Preparations ###

# netFeedReflux = ppi.from_pickle(str(flp_out / name) + '_netFeedReflux_dimensioned_manually_new.p')
netFeedReflux = ppi.from_pickle(r'examples/data/exampleNetwork/results/netFeedReflux_dimensioned_manually_final_controlled.p')

cp_fluid, rho_fluid, nu_fluid, g = extract_FluidProperties_ppi(net = netFeedReflux, t = 0.5 * (tFeed + tReflux))

# ppi.drop_elements_at_junctions(netFeedReflux, netFeedReflux.flow_control[['from_junction', 'to_junction']].values.flatten())

ppi.set_user_pf_options(net = netFeedReflux, mode="sequential", friction_model = 'swamee-jain', quit_on_inconsistency_connectivity=True, reset = True)


# %% Scratch function for summing up heat demands along path to closest heat supply site
import numpy as np
import networkx as nx
import pandas as pd

# Length attribute definition
weightAttr = 'length_km'

colsToKeep = ['from_junction', 'to_junction', weightAttr, 'geometry']

# Prepare dataframe to convert
## Keep opened valves

if hasattr(netFeedReflux, 'valve'):

    DF_edges = pd.concat((
        netFeedReflux.pipe[netFeedReflux.pipe['Layer'] == 'feedLine'][[n for n in colsToKeep if n in netFeedReflux.pipe.columns]], 
        netFeedReflux.valve[(netFeedReflux.valve['Layer'] == 'feedLine') & (netFeedReflux.valve['opened'] == True)][[n for n in colsToKeep if n in netFeedReflux.valve.columns]]
        ), axis = 0)
    
else:
    DF_edges = netFeedReflux.pipe[netFeedReflux.pipe['Layer'] == 'feedLine'][[n for n in colsToKeep if n in netFeedReflux.pipe.columns]]

DF_edges.loc[DF_edges[weightAttr].isna(), weightAttr] = 0


# %%



nCons = netFeedReflux.junction[netFeedReflux.junction.index.isin(netFeedReflux.heat_consumer['from_junction'])]

nPros = netFeedReflux.junction[
    netFeedReflux.junction.index.isin(
        netFeedReflux.circ_pump_pressure['flow_junction'].to_list() + 
        [netFeedReflux.flow_control.loc[0, 'to_junction']]
        )
    ]

startLocs = list(nPros.index)
endLocs = list(nCons.index)


# %% Vorgehen mit pandas edgelist

DF_edges['ID'] = np.arange(len(DF_edges))




summedAttrs = dict(zip(netFeedReflux.heat_consumer['from_junction'], netFeedReflux.heat_consumer['demand_use_th'].fillna(0)))

## Find shortest paths from starting locations to end locations and their respective path lengths




# %%


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



# %%

DF_out = sum_heat_demands_to_closest_supplier(
    edgelist = DF_edges, 
    sources = 'from_junction', 
    targets = 'to_junction', 
    startNodes = list(nPros.index), 
    endNodes = list(nCons.index), 
    edge_key = 'ID', 
    weight = 'length_km', 
    outAttr_weight = 'summed_demand_use_th', 
    dictAttrsEndNodes=dict(zip(netFeedReflux.heat_consumer['from_junction'], netFeedReflux.heat_consumer['demand_use_th'].fillna(0)))
    )


DF_out


# %%


gp.GeoDataFrame(DF_edges, geometry = 'geometry').set_crs(cs).to_file(r'examples/data/exampleNetwork/results/testoutput.gpkg', driver = 'GPKG')

# %%
import matplotlib.pyplot as plt


G = nx.MultiGraph()

# Vier Knoten
G.add_nodes_from([1, 2, 3, 4])

# Masche: zwei parallele Kanten zwischen 1 und 2 (Schleife/Masche)
G.add_edge(1, 2, key=1, len=1)
G.add_edge(3, 4, key=4, len=1)
G.add_edge(2, 3, key = 2,len=2)
G.add_edge(2, 3, key = 3,len=1)

# Initlaization
nx.set_edge_attributes(G, values = 0, name = 'sum_weight')


# Zeichnen (optional)
pos = {1: (0, 0), 2: (1, 0), 3: (2, 0), 4: (3, 0)}
nx.draw(G, pos, with_labels=True, connectionstyle="arc3,rad=0.2")
plt.show()

startLocs = [4]
endLocs = [2, 1]


# %%

## Find shortest paths from starting locations to end locations and their respective path lengths

paths = []
lengths = []

for n in startLocs:

    path = dict(nx.single_source_all_shortest_paths(G = G, source = n, weight = 'len'))
    length = dict(nx.shortest_path_length(G = G, source = n, weight = 'len'))

    path_new = {k:v for k,v in path.items() if k in endLocs}
    length_new = {k:v for k,v in length.items() if k in endLocs}

    paths.append(path_new)
    lengths.append(length_new)

allKeys = set().union(*(d.keys() for d in lengths))
shortestLength = {}
shortestPath = {}

summedAttrs = {k:np.random.randint(1, 10)*100 for k in allKeys}


for key in allKeys:
    min_len = float('inf')
    min_idx = -1

    for i, d in enumerate(lengths):
        if key in d and d[key] < min_len:
            min_len = d[key]
            min_idx = i
        if min_idx >= 0:
            shortestLength[key] = min_len
            shortestPath[key] = paths[min_idx][key]



for k, path in shortestPath.items():
    for u, v in zip(path[0][:-1], path[0][1:]):

        edges = G[u][v] # -> All edges betwen nodes u and v

        # Address possible parallel edges
        min_key, min_attr = min(edges.items(), key=lambda item: item[1].get('len', float('inf')))

        G[u][v][min_key]['sum_weight'] += summedAttrs[k]

# Initialize final output dictionary for matching of line ID and summed Attribute (sumemd weight)
edge_attr_dict = {}
for u, v, k, attr in G.edges(keys=True, data=True):
    edge_attr_dict[k] = attr.get('sum_weight')




# %%


