
# %% Imports
import os
import sys

parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if parentdir not in sys.path:
    sys.path.insert(0, parentdir)

import geopandas as gp
import time
from pathlib import Path

from geoppi.suitable_network_routing import (network_span_bfs)
from geoppi.internal_auxFunctions import (simultaneity_DH, relThermalLoss_DH, relThermalLossPower_DH,)

parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if parentdir not in sys.path:
    sys.path.insert(0, parentdir)

# %% Load data

flp = parentdir / Path(r'examples/data/exampleNetwork')
flp_out = flp / Path(r'results')


# Coordinate system
cs = 'EPSG:25832'

lines                   = gp.read_file(flp / Path(r'networkRaw.gpkg')).explode(index_parts = False).to_crs(cs)
lines['length']         = lines.geometry.length
lines['demand_use_th']  = lines['ld_demand_use_th'] * lines.geometry.length
lines['p_use_th']       = lines['demand_use_th'] / 1700

startpoints             = gp.read_file(flp / Path(r'startpoints_producers2.gpkg')).to_crs(cs)

backgroundLines         = lines.copy()

### Limits & attributes

# Attribute for line density
att_ld = 'ld_demand_use_th'

# Attribute for heat demand
att_heatDemand = 'demand_use_th'

# Attribute for thermal power on line (MW)
att_pth = 'p_use_th'

# Attribute for length of line
att_len = 'length'

# Attribute für die Anzahl an Gebäuden pro Straßenzug
att_nB = 'nBuildsdemand_use_th'

# Energy budget for production at defined starting points (MWh)
energyBudget = [300]

# Thermal power budget for production at defined starting points (MW)
powerBudget = [0.5]

# mindest Anzahl an Gebäuden pro ein Netz
n_build = 5

# Sort edges by edge attribute in descending order (sorting per node adjacent edges)
sortByAttr = att_ld


# %% Start algorithm

    
starttime = time.time()

gdf_graphs, graphList = network_span_bfs(
    lines = lines,
    startpoints = startpoints.loc[[0]],
    val_Attr = att_ld,
    att_pth = att_pth,
    length_Attr = att_len,
    att_nB = att_nB,
    min_nB = n_build,
    att_heatDemand = att_heatDemand,
    sortByAttr = sortByAttr,
    sortMethod = 'descending',
    considerPowerBudget = True,
    energyBudget = energyBudget,
    powerBudget = powerBudget,
    adaptThermalLoss = relThermalLoss_DH,
    adaptSimultaneityFactor = simultaneity_DH,
    adaptThermalPowerLoss = relThermalLossPower_DH,
    minvalAttr_maxLength = {-0.1:25},
    createMST = False,
    returnAddResults = False
)

endtime = time.time()
timeElapsed = endtime - starttime
print(f'\n### Time elapsed for algorithm is {timeElapsed:.2} seconds ###')


## Plotting
import matplotlib.pyplot as plt

fig, ax = plt.subplots()
gdf_graphs.plot(ax = ax, column = att_ld, vmin = 0, vmax = 5, legend = True)

startpoints.plot(ax = ax, color = 'red', zorder = 10)

# Export
if not os.path.exists(flp_out):
    os.makedirs(flp_out)
    print(f'### New directory for storing results is created! ###')

# Save edges
output_path = flp_out / Path(f'exmapleNetwork_extension.gpkg')

# gdf_graphs.to_file(output_path, layer = f'MW_Pthprod2_{powerBudget[0]:.1f}', driver="GPKG")

