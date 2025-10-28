# -*- coding: utf-8 -*-


# %% Imports


import geopandas as gp
from shapely import set_precision
from pathlib import Path
import time

from suitable_network_routing import (network_span_dfs_level_search)
from internal_auxFunctions import (animate_networkGeneration_gdf,)

# %% Load data

flp = Path('C:/GitLab/paper_generic_network_creation/data/lineDensity')
flp_out = flp

# Coordinate system
cs = 'EPSG:25832'

lines               = gp.read_file(flp / Path(r'Streets_LineDensity_bestand_renov_v1.gpkg')).explode(index_parts = False).to_crs(cs)
lines['length']     = lines.geometry.length
lines['geometry']   = set_precision(lines.geometry, 0.01)
lines['inverse_ld_demand_use_th'] = 1/lines['ld_demand_use_th']

startpoints         = gp.read_file(flp / Path(r'startpointsProducers.gpkg')).to_crs(cs)
 
consumerConnections = gp.read_file(flp / Path(r'Points_consumerIntersection.gpkg')).to_crs(cs)


# %%  Limits & attributes

# Attribute for line density as main evaluation attribute on line objects
att_ld = 'ld_demand_use_th'#total_exp_ren_1.5_2045_climate_demo'

# Attributes for length of street segment
lengthAttr = 'length'

# Attribute for number of buildings connected to street segment
att_nB = 'nBuildsdemand_use_th'#total_exp_ren_1.5_2045_climate_demo'

# Minimum value for automatic selection of start edges
dem_val_start = 4

# Minimum value for 1st level search
valFirstlvl = 1.499

# Mean line density for 2nd level search (second Level)
valMeanSecondlvl = 1.499

# Mean line density for 3rd level search (third Level)
valMeanThirdlvl = 2

# Minimum length of street segments considered in 2nd and 3rd level
minLength = 0.1

# Minimum length of street segments omitted
minValJump = 0

# Maximum length of street segments omitted
maxLengthJump = 250

# Threshold for line density attribute below which algorithm removes street segments sorted with descending length
dem_val_remove_edge = 0.9

# Minimum number of buildings connected to entire networks
n_build = 17


# %% Perform algorithm

starttime = time.time()

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
    valStart = dem_val_start,
    removeDuplicateGraphs = True,
    mergeTouchingGraphs = True,
    complementSingleStreets = True,
    removeSmallNetworks = True
)

endtime = time.time()
timeElapsed = endtime - starttime
print(f'\n### Time elapsed for algorithm is {timeElapsed:.2} seconds ###')

import matplotlib.pyplot as plt

fig, ax = plt.subplots()
gdf_graphs.plot(ax = ax, column = att_ld)

startpoints.plot(ax = ax, color = 'red', zorder = 10)


# %% *--- Template: Create animated visualization of development for network routing ---*

lines_final = gdf_graphs.copy()
lines_final.sort_values(by = 'order', inplace = True, ascending = True)
lines_final.reset_index(drop = True, inplace = True)

lines_final['length_cum'] = lines_final['length'].cumsum()
lines_final['heat_Demand'] = lines_final['ld_demand_use_th'] * lines_final['length']
lines_final['heat_Demand_cum'] = lines_final['heat_Demand'].cumsum()

lines_final['ld_demand_use_th_mean_ordered'] = lines_final['heat_Demand_cum'] / lines_final['length_cum']


animate_networkGeneration_gdf(
    gdf = lines_final,
    valAttr = 'ld_demand_use_th',
    cbarLabel = 'line density (MWh/m)',
    cmapMinMax = (0, 5),
    plotOriginalgdf = True,
    legendAttr = {'order ': 'order', 'level ':'level', 'mean ld (MWh/m) = ':'ld_demand_use_th_mean_ordered', 'current ld (MWh/m) = ':'ld_demand_use_th', 'summed length (m) = ':'length_cum'},
    savePath=Path(r'C:\Users\Voelzel\Desktop\videos'),
    filename = 'video.gif'
)


# %% *--- Save output data ---*

gdf_graphs['year'] = 2025

# Save edges
output_path = flp_out / Path(f'graphs_val1{str(valFirstlvl).replace(".", "p")}_val2{str(valMeanSecondlvl).replace(".", "p")}_val3{str(valMeanThirdlvl).replace(".", "p")}_valRemove{str(dem_val_remove_edge).replace(".", "p")}_2025_test20250512.gpkg')

gdf_graphs.to_file(output_path, driver="GPKG")

