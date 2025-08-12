# %% Imports


import geopandas as gp
import pandas as pd
import os
import time
from pathlib import Path

from suitable_network_routing import (network_span_bfs)
from internal_auxFunctions import (simultaneity_DH, relThermalLoss_DH, relThermalLossPower_DH, animate_networkGeneration_gdf,)


# %% Load data

flp = Path(r'C:\GitLab\paper_generic_network_creation\data\lineDensity')
flp_out = Path(r'C:\Users\Voelzel\Desktop\tempRes')



# Coordinate system
cs = 'EPSG:25832'

lines               = gp.read_file(flp / Path(r'Streets_LineDensity_bestand_renov_v1.gpkg')).explode(index_parts = False).to_crs(cs)
lines['length']     = lines.geometry.length
lines['demand_use_th'] = lines['ld_demand_use_th'] * lines.geometry.length
lines['p_use_th'] = lines['demand_use_th'] / 1700

mask                = gp.read_file(flp / Path('mask_Networks.gpkg')).to_crs(cs)
lines               = gp.sjoin(left_df = lines, right_df = mask[['geometry']], how = 'left', predicate = 'intersects')
lines               = lines[~lines['index_right'].isna()].reset_index(drop = True).drop(columns = 'index_right')

startpoints         = gp.read_file(flp / Path(r'startpointsProducers.gpkg')).to_crs(cs)

backgroundLines = gp.read_file(flp / Path(r'Streets_LineDensity_bestand_renov_v1.gpkg')).to_crs(cs)

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
energyBudget = [[10e04]]#, 10e04]]#[[8000, 10e04, 10000], [16000, 10e04, 20000], [22000, 10e04, 20000]]

# Thermal power budget for production at defined starting points (MW)
powerBudget = [[12.5]]#, 12.5]]#, 12.5, 8]]#[[4, 12.5, 4], [8, 12.5, 8], [18, 12.5, 8]]

# mindest Anzahl an Gebäuden pro ein Netz
n_build = 17

# Sort edges by edge attribute in descending order (sorting per node adjacent edges)
sortByAttr = att_ld



# %% Scratch

redGraphs = gp.read_file(Path(r'C:\Users\Voelzel\Desktop\tempRes\graphs_extension_2025_v2_reducedInputGraph.gpkg'))

# %% Scratch

avLD = redGraphs['demand_use_th'].sum() / redGraphs.geometry.length.sum()

Qloss = relThermalLoss_DH(avLD)*redGraphs['demand_use_th'].sum()
Qdem = (Qloss + redGraphs['demand_use_th'].sum())/(1e03)

Pdem = redGraphs['p_use_th'].sum() * 0.7

print(f'\n Pdem = {Pdem} MW')
print(f'\n Qdem = {Qdem} MWh')

print(f'\n Qsupply = {10e04 + 22000} MWh')
print(f'\n Psupply = {18+12.5} MW')


# %% Start algorithm
# def relThermalLossPower_DH(pd, referToInput:bool = False):

    # return 0.05

# def simultaneity_DH(n:int, bottomLim:float = 0.7):

    # return 1

# def relThermalLoss_DH(ld, referToInput:bool = False):

    # return 0.15

for r, run in enumerate(powerBudget):
    
    starttime = time.time()

    gdf_graphs, graphList, SFList1, usedEnergyBudgetList1, usedPowerBudgetList1, currentThermalLossFactorList1, summedLengthList1, avlineDensityList1 = network_span_bfs(
        lines = lines,
        startpoints = startpoints.loc[[1]],
        val_Attr = att_ld,
        att_pth = att_pth,
        length_Attr = att_len,
        att_nB = att_nB,
        min_nB = n_build,
        att_heatDemand = att_heatDemand,
        sortByAttr = sortByAttr,
        sortMethod = 'descending',
        considerPowerBudget = True,
        energyBudget = energyBudget[r],
        powerBudget = powerBudget[r],
        adaptThermalLoss = relThermalLoss_DH,
        adaptSimultaneityFactor = simultaneity_DH,
        adaptThermalPowerLoss = relThermalLossPower_DH,
        minvalAttr_maxLength = {0.5:25},
        createMST = False,
        returnAddResults = True,
        reduceInputGraphMultipleStartpoints = True
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
        print('### New directory for storing results is created! ###')

    # Save edges
    # output_path = flp_out / Path(f'graphs_extension_2025_v3_relLoss0p05_relLoss0p15_simultaneity1.gpkg')
    output_path = flp_out / Path(f'graphs_extension_2025_v2_OnlyProd2.gpkg')

    # gdf_graphs.to_file(output_path, layer = f'Pthprod1_{powerBudget[r][0]:.1f}MW_Pthprod2_{powerBudget[r][1]:.1f}', driver="GPKG")

    gdf_graphs.to_file(output_path, layer = f'Pthprod1_{powerBudget[r][0]:.1f}', driver="GPKG")

# %% *--- Template: Create plot over number of street segments in network extension analysis ---*

# Plot appearance
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import numpy as np



xfont = {'fontname': 'Times New Roman', 'fontsize':9}
legendfont = {'family': 'Times New Roman', 'size': 9}

fig, ax = plt.subplots(1, 1, sharex=True)
fig.set_size_inches(3.22, 2, forward=True)
fig.set_dpi(300)

# Data
x1 = summedLengthList1[0]
x2 = summedLengthList2[0]

ax2 = ax.twinx()

# Plot SF
ax.plot(x1, SFList1[0], color = '0.05', label = r'simultaneity factor sf', linestyle = '-')
ax.plot(x2, SFList2[0], color = '0.05', label = r'simultaneity factor sf', linestyle = ':')

# Plot used energy budget
ax.plot(x1, np.array(usedEnergyBudgetList1[0])/22000, color = '0.3', label = r'proportion of used energy budget $\frac{Q_{dem}+Q_{loss}}{Q_{supply,potential}}$', linestyle = '-')
ax.plot(x2, np.array(usedEnergyBudgetList2[0])/22000, color = '0.3', label = r'proportion of used energy budget $\frac{Q_{dem}+Q_{loss}}{Q_{supply,potential}}$', linestyle = ':')

# Plot thermal energy losses
ax.plot(x1, currentThermalLossFactorList1[0], color = '0.6', label = r'rel. thermal loss f$_{rel}$', linestyle = '-')
ax.plot(x2, currentThermalLossFactorList2[0], color = '0.6', label = r'rel. thermal loss f$_{rel}$', linestyle = ':')

# Plot av. line density
ax2.plot(x1[1:], avlineDensityList1[0][1:], color = '0.85', label = r'average line density ld', linestyle = '-')
ax2.plot(x2[1:], avlineDensityList2[0][1:], color = '0.85', label = r'average line density ld', linestyle = ':')


# Set axes limits

ax.set_yticks(ax.get_yticks(), labels = [str(tick) for tick in ax.get_yticks()], **legendfont)
ax2.set_yticks(np.array([0, 1.5, 3, 4.5, 6]), labels = [str(tick) for tick in np.array([0, 1.5, 3, 4.5, 6])], **legendfont)

ax.set_xlabel(r'Total trail length (km)', **xfont)
axXticks = np.append(np.unique(np.arange(10)*1e03).round(0).astype(int), int(max(max(x1), max(x2))))
ax.set_xticks(axXticks, labels = [str(np.round(tick/1e03, 1)) for tick in axXticks], **legendfont)

ax.set_ylabel('sf;' + '   ' + 'f$_{rel}$;' + '   ' + r'$\frac{Q_{dem}+Q_{loss}}{Q_{supply,potential}}$', **xfont)
ax2.set_ylabel(r'Line density ($\frac{MWh}{m}$)', **xfont)

ax.set_xlim([0, max(max(x1), max(x2))])
ax.set_ylim([0, 1.005])
ax2.set_ylim([0, 6])

ax.grid('minor')


# Create legend
lines_labels = ax.get_legend_handles_labels()
lines_labels2 = ax2.get_legend_handles_labels()

handles = lines_labels[0][0::2] + lines_labels2[0][0::2]
labels = lines_labels[-1][0::2] + lines_labels2[-1][0::2]

ax2.legend(handles, labels, loc='lower center', bbox_to_anchor = (0.5, -0.9), ncol = 1, prop = legendfont)

# fig.savefig(flp_out / Path('fig_networkExtension_variation.png'), dpi = 300, bbox_inches = 'tight')


# %% *--- Template: Create animated visualization of development for network routing ---*
import numpy as np

lines_final = gdf_graphs.copy()
lines_final.sort_values(by = 'order', inplace = True, ascending = True)
lines_final.reset_index(drop = True, inplace = True)

lines_final['length_cum'] = lines_final['length'].cumsum()
lines_final['heat_Demand'] = lines_final['ld_demand_use_th'] * lines_final['length']
lines_final['heat_Demand_cum'] = lines_final['heat_Demand'].cumsum()

lines_final['ld_demand_use_th_mean_ordered'] = lines_final['heat_Demand_cum'] / lines_final['length_cum']

lines_final['prop_usedEnergyBudget'] = np.array(usedEnergyBudgetList1[0][:74]) / 15000
lines_final['prop_usedPowerBudget'] = np.array(usedPowerBudgetList1[0][:74]) / 7
lines_final['SF'] = SFList1[0][:74]
lines_final['Qloss_factor'] = currentThermalLossFactorList1[0][:74]


animate_networkGeneration_gdf(
    gdf = lines_final,
    valAttr = 'ld_demand_use_th',
    underlyinggdf = backgroundLines,
    cbarLabel = 'line density (MWh/m)',
    cmapMinMax = (0, 3),
    legendAttr = {'order =': 'order', 'av. ld (MWh/m) =':'ld_demand_use_th_mean_ordered', 'current ld (MWh/m) =':'ld_demand_use_th', 'network length (m) =':'length_cum', 'prop. used energy budget =':'prop_usedEnergyBudget', 'prop. used power budget =':'prop_usedPowerBudget','SF =':'SF'},
    interval = 500,
    savePath=Path(r'C:\Users\const\Desktop\videos'),
    filename = 'video_bfs.gif'
)

