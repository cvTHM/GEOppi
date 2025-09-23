# %% Imports
import os
import sys

parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if parentdir not in sys.path:
    sys.path.insert(0, parentdir)

import rasterio
from pathlib import Path
import geopandas as gp
import pandapipes as ppi
import numpy as np
from geoppi.auxFunctions import (get_dict_from_aggregated_groups)
from geoppi.create_network_topology import (create_ppi_network_from_gdf,)
from geoppi.suitable_network_routing import (MST_gdf_subset, sum_heat_demands_to_closest_supplier)
from geoppi.dimension_network_pipes import (update_ppi_results, assign_insulation_type, hydraulicDimensioningNetwork_singleLoadPoint, assign_nominal_widths_manually)
from geoppi.internal_auxFunctions import (extract_FluidProperties_ppi, transfer_LoadPoint_ppi, implement_controllers)
from pandapipes.control import run_control

from geoppi.pipe_characteristics import load_hydraulic_pipe_characteristics


# %% Load data

# Define coordinate system
cs = 'EPSG:25832'

# Define network name
name = 'exampleNetwork'

# Select whather network shall be saved as pickle file
save_network = True

### Filepaths
# Input path
flp = parentdir / Path(r'examples/data/exampleNetwork')

# output path
flp_out = flp / Path(r'results')

if not os.path.exists(flp_out):
    os.makedirs(flp_out)
    print('### New directory for saving pandapipes network is created! ###')
            
# output name
fln_out = 'network_' + name

# Load pipes
pipes = gp.read_file(flp / Path('networkRaw.gpkg')).to_crs(cs)

# Load Buildings layer
buildings = gp.read_file(flp / Path('buildings.gpkg')).to_crs(cs)

# Load valve Layer
valves = gp.read_file(flp / ('valves.gpkg')).to_crs(cs)

# Load production site layer
producers = gp.read_file(flp / Path('producers2.gpkg')).to_crs(cs)

# Load elevation raster layer
dgm = rasterio.open(flp / Path('dgm.tif'))

# Load optional point layer for creating MST
# points_MST = gp.read_file(flp / ('points_ConsumersNetworkCalc.gpkg'))
# points_MST2 = producers.geometry.centroid

# import pandas as pd
# points_MST = pd.concat((points_MST, points_MST2), axis = 0)
# pipes = MST_gdf_subset(lines = pipes, subsetNodesMST = points_MST)

print('######### Load data DONE #########')

# %%

# *--- Parameters for network creation ---*

# Define unique ID for identification of buildings
buildings_uniqueID = 'build_ID'

if buildings_uniqueID not in buildings.columns:
    buildings[buildings_uniqueID] = np.arange(len(buildings))

# Define attribute containing heating demand information
heatingDemandAttr = 'demand_use_th'

# * --- Parameters for network pipe dimensioning ---*
### Network type (either "KMR" or "5GDHC")
networkType = 'KMR'

### Set temperatures in the network (°C)
tFeed = 65
tReflux = 40

# Set network pressure at peak load producers (bar)
pNetwork = 8

### Valves
# Selection of valves to be closed ("-1": all closed)
valves_closed = [-1]

# Full load hours for calculation of heating power from consumers' heating demand
flh = (1700)

# Define attribute name containing information about thermal power (kW) of producers in producer models
thermalPowerAttr = 'Pth_kW'

# Define attribute name containing infotmation about distribution ('VE') and house connection pipes ('HA')
connectionTypeAttr = 'connectionType'

# Define naming convention for base load producers in heat_consumer component DataFrame
baseLoadProdNaming = 'prod'

# Target specific pressure loss in pipes (Pa/m, usually 60-100 Pa/m)
dp_spec = 100

# Selection of subset for nominal widths for pipes
subsetDistribution = [50, 65, 80, 100, 125, 150, 200, 225, 250, 300, 350]
subsetConnection = [32, 40, 50, 65, 80, 100, 125, 150, 200, 225, 300]

# Desired insulation for pipes
insulationType = 'std'

# Medium in the network
fluid = 'water'

# Factor for elongation of pipes to consider influence of armatures etc. on pressure loss
factor_l = 1.3

# Pressure loss calculation model for pandapipes
friction_model = 'swamee-jain'

# Selection if real heights at junctions shall be considered
respect_height = False

# Access network type-dependent data for pipes
dictNominalInnerDiameter, dictWallRoughness = load_hydraulic_pipe_characteristics(networkType)

# Global setting to suppress numpy warnings
np.seterr(divide='ignore')


# %% Create pandapipes network for feed line

netFeed = create_ppi_network_from_gdf(
    pipes = pipes,
    buildings = buildings,
    valves = valves,
    producers = producers,
    rasterHeight = dgm,
    buildings_uniqueID = 'build_ID',
    heatingDemandAttr = heatingDemandAttr,
    networkName = 'exampleNet',
    networkFluid = fluid,
    producerTypeDict = {'Type':{'peak':'peak','base':'base'}},
    producerAddAttr = ['Pth_kW', 'dp_bar_internal'],
    modelling_type = 'feed_line'
)


# %% Set parameters for described load point

netFeed = transfer_LoadPoint_ppi(
    net = netFeed,
    tFeed = tFeed + 273.15,
    tReflux = tReflux + 273.15,
    pNetwork = pNetwork,
    flh = flh,
    heatingDemandAttr = heatingDemandAttr,
    thermalPowerAttr = thermalPowerAttr,
    baseLoadProdNaming = baseLoadProdNaming,
    check_for_mass_flow_exceeding = True
)

# %% Start algorithm for iterative network dimensioning

cp_fluid, rho_fluid, nu_fluid, g = extract_FluidProperties_ppi(net = netFeed, t = 0.5 * (tFeed + tReflux))

# Adapt user defined options of the network to hydraulic and thermal calculation mode
ppi.set_user_pf_options(net = netFeed, mode="hydraulics", friction_model = friction_model, quit_on_inconsistency_connectivity=True, reset = True)

netFeed = hydraulicDimensioningNetwork_singleLoadPoint(
    net = netFeed,
    networkType = networkType,
    dpSpec = dp_spec,
    dictNominalInnerDiameter = {'distributionPipes':dictNominalInnerDiameter, 'connectionPipes':dictNominalInnerDiameter},
    dictWallRoughness = dictWallRoughness,
    subsetDistribution = subsetDistribution,
    subsetConnection = subsetConnection,
    connectionTypeAttr = connectionTypeAttr,
    elongationFactorPipes = factor_l,
    nominalWidthAttr = 'nominalWidth',
    closedValvesIdx = valves_closed,
    frictionModel = friction_model,
    respectHeight = False,
    fluidProperties = {'rho':rho_fluid, 'cp':cp_fluid, 'nu':nu_fluid},
    greaterWidthFinalStep={'distributionPipes':0, 'connectionPipes':0}
)


### Assign final heat transfer coefficient
netFeed.pipe = assign_insulation_type(
    df_pipes = netFeed.pipe, 
    colNameNominalWidth = 'nominalWidth', 
    network_type = networkType, 
    insulation_type = insulationType if networkType == 'KMR' else None
)

netFeed = update_ppi_results(
    net = netFeed,
    user_pf_options={'mode':'hydraulics', 'friction_model':'swamee-jain'},
    respectHeight = False,
    fluidProperties = {'rho':rho_fluid, 'cp':cp_fluid, 'nu':nu_fluid},
    elongationFactorPipes = 1.3,
    resultAttributes={ 'pipe': ['Pa_per_m', 'mdot_from_kg_per_s'],'junction': ['t_k', 'p_bar'] }
)

# %% *--- Optional: Intermediate plots ---*
import matplotlib.pyplot as plt

gdf = gp.GeoDataFrame(netFeed.pipe, geometry = 'geometry').set_crs(cs)
gdf['mdot_kg_per_s'] = abs(gdf['mdot_from_kg_per_s'])

fig, ax1 = plt.subplots(3, figsize = (30,20))

styles = {'linewidth':5}
xfont = {'fontsize': 20, 'weight':'normal'}  # Font for axes
tfont = {'fontsize': 20}  # Font for title

ax1[0].set_title('Nominal widths', **tfont)
gdf[(gdf['in_service'] ==True) & (gdf[connectionTypeAttr] == 'distribution')].plot(
    ax = ax1[0], 
    column = 'nominalWidth', 
    cmap=plt.get_cmap('Paired'),
    categorical = True,
    missing_kwds = dict(color = 'grey', label = '-'),
    legend = True,
    legend_kwds = {'title':'nominal widths.', 
                   'loc':'lower right',
                   'markerscale':1, 
                    'title_fontsize':'xx-large', 
                    'fontsize':'xx-large'},
    **styles)
leg = ax1[0].get_legend()
leg.set_bbox_to_anchor((1.1,0,0.2,0.5))

ax1[1].set_title('Affiliation with loop', **tfont)
gdf[(gdf['in_service'] ==True) & (gdf[connectionTypeAttr] == 'distribution')].plot(
    ax = ax1[1], 
    column = 'loop', 
    cmap=plt.get_cmap('Paired'),
    legend = True,
    categorical = True,
    missing_kwds = dict(color = 'grey', label = '-'),
    legend_kwds = {'title':'Loop no.', 
                   'loc':'lower right',
                   'markerscale':1, 
                    'title_fontsize':'xx-large', 
                    'fontsize':'xx-large'},
    **styles)
leg = ax1[1].get_legend()
leg.set_bbox_to_anchor((1.1,0,0.2,0.5))

ax1[2].set_title('Spec. pressure loss', **tfont)
gdf[(gdf['in_service'] ==True) & (gdf[connectionTypeAttr] == 'distribution')].plot(
    ax = ax1[2], 
    column = 'Pa_per_m', 
    cmap=plt.get_cmap('magma'),
    vmin = 0,
    vmax = 200, 
    legend = True,
    **styles)



# %% *--- Manually adapt nominal widths ---*
_, dictNWs = get_dict_from_aggregated_groups(
    df = netFeed.pipe,
    groupCol = 'loop',
    val = 'nominalWidth',
    func = 'max'
)

netFeed_dimensioned = assign_nominal_widths_manually(
    net = netFeed,
    dictNW = dictNWs,
    dictNominalInnerDiameter = dictNominalInnerDiameter,
    dictWallRoughness = dictWallRoughness,
    nominalWidthAttr = 'nominalWidth'
)

### Assign final heat transfer coefficient
netFeed_dimensioned.pipe = assign_insulation_type(
    df_pipes = netFeed_dimensioned.pipe, 
    colNameNominalWidth = 'nominalWidth', 
    network_type = networkType, 
    insulation_type = insulationType if networkType == 'KMR' else None
)



# %% *--- Create complete Pandapipes network from pre-dimensioned feed line network ---*

netFeedReflux = create_ppi_network_from_gdf(
    existingNetworkFeedLine = netFeed_dimensioned,
    buildings = buildings,
    heatingDemandAttr = heatingDemandAttr,
    producers = producers,
    networkName = 'exampleNet',
    networkFluid = 'water',
    producerTypeDict = {'Type':{'peak':'peak','base':'base'}},
    producerAddAttr = ['Pth_kW', 'dp_bar_internal'], # Necessary thermal power attribute to write into pandapipes tables
    modelling_type = 'feed_and_reflux_line_circPumps'
)

# Set user-defined pipeflow option for sequential calculation!
ppi.set_user_pf_options(net = netFeedReflux, reset = True, mode='sequential', friction_model = friction_model, quit_on_inconsistency_connectivity=True)

netFeedReflux = transfer_LoadPoint_ppi(
    net = netFeedReflux,
    tFeed = tFeed + 273.15,
    tReflux = tReflux + 273.15,
    text_pipes = 280,
    pNetwork = pNetwork,
    flh = flh,
    heatingDemandAttr = heatingDemandAttr,
    thermalPowerAttr = thermalPowerAttr,
    baseLoadProdNaming = baseLoadProdNaming
)

netFeedReflux = update_ppi_results(
    net = netFeedReflux,
    respectHeight=False,
    fluidProperties={'cp':cp_fluid, 'rho':rho_fluid, 'nu':nu_fluid},
    elongationFactorPipes=factor_l,
    resultAttributes = {'pipe':['Pa_per_m', 'mdot_from_kg_per_s'], 'junction':['t_k', 'p_bar'], 'heat_consumer':['t_from_k', 't_to_k']}
)


# %% Analysis of summed heat demands along paths from heat consumers to closest heat supplier
import pandas as pd

# Prepare dataframe to convert
## Keep opened valves
weightAttr = 'length_km'
colsToKeep = ['from_junction', 'to_junction', weightAttr, 'geometry']

if hasattr(netFeedReflux, 'valve'):

    DF_edges = pd.concat((
        netFeedReflux.pipe[netFeedReflux.pipe['Layer'] == 'feedLine'][[n for n in colsToKeep if n in netFeedReflux.pipe.columns]], 
        netFeedReflux.valve[(netFeedReflux.valve['Layer'] == 'feedLine') & (netFeedReflux.valve['opened'] == True)][[n for n in colsToKeep if n in netFeedReflux.valve.columns]]
        ), axis = 0)
    
else:
    DF_edges = netFeedReflux.pipe[netFeedReflux.pipe['Layer'] == 'feedLine'][[n for n in colsToKeep if n in netFeedReflux.pipe.columns]]

DF_edges.loc[DF_edges[weightAttr].isna(), weightAttr] = 0

DF_out = sum_heat_demands_to_closest_supplier(
    edgelist = DF_edges,
    sources = 'from_junction',
    targets = 'to_junction',
    startNodes = netFeedReflux.circ_pump_pressure['flow_junction'].to_list() + [netFeedReflux.flow_control.loc[0, 'to_junction']],
    endNodes = list(netFeedReflux.heat_consumer['from_junction']),
    weight = weightAttr,
    outAttr_weight = 'summed_demand_use_th',
    dictAttrsEndNodes = dict(zip(netFeedReflux.heat_consumer['from_junction'], netFeedReflux.heat_consumer['demand_use_th'].fillna(0)))
)


gdf = gp.GeoDataFrame(DF_out, geometry = 'geometry').set_crs(cs)
fig, ax1 = plt.subplots(1, figsize = (30,20))

styles = {'linewidth':5}
xfont = {'fontsize': 20, 'weight':'normal'}  # Font for axes
tfont = {'fontsize': 20}  # Font for title

ax1.set_title('Summed annual heat demand (kWh)\nfrom closest heat supplier to heat consumers', **tfont)
gdf.plot(
    ax = ax1, 
    column = 'summed_demand_use_th', 
    cmap=plt.get_cmap('magma'),
    vmin = 0,
    vmax = gdf['summed_demand_use_th'].max()*0.75, 
    legend = True,
    **styles)


# %% Implement controllers
# Define control targets
circ_pump_pressure_idx = 0

dpmin_target = 1.2 # bar
pmin_target = 1.5 # bar

# Create controllers 
netFeedReflux = implement_controllers(
    net = netFeedReflux,
    drop_all = True,
    pminCtrlDict = {
        'create':True, 
        'pmin_target':pmin_target, 
        'circ_pump_pressure_idx':circ_pump_pressure_idx,
        'abs_tol':0.1,
        'order':1,
        'level':5
        },
    dpminCtrlDict = { # Lowest control order in level 5 = highest priority in level 5
        'create':True, 
        'dpmin_target':dpmin_target, 
        'circ_pump_pressure_idx':circ_pump_pressure_idx,
        'index':1,
        'abs_tol':0.1,
        'order':0,
        'level':5
        },
    TRefluxHeatConsumerCtrlDict = { # Lowest control level = highest priority
        'create':True,
        'T_target':tReflux + 273.15,
        'min_mdot':0.015,
        'min_dT':3,
        'abs_tol':1,
        'order':0,
        'level':-1
        },
    PthLimitedCtrlDict = {
        'create':True,
        'circ_pump_mass_idxs':netFeedReflux.circ_pump_mass.index, 
        'circ_pump_pressure_index':circ_pump_pressure_idx, 
        'flow_controller_idxs':netFeedReflux.flow_control.index,
        'Pth_target_kW':netFeedReflux.circ_pump_mass['Pth_kW'].values, 
        'priority_list':[0],
        'abs_tol':5000,
        'order':5,
        'level':5
        }
    )

print(netFeedReflux.controller)

# Run pipeflow with controls
run_control(net = netFeedReflux,  max_iter = 50)


# %%
# Print exemplary results
netFeedReflux.res_heat_consumer['dp_bar'] = netFeedReflux.res_heat_consumer['p_from_bar'] - netFeedReflux.res_heat_consumer['p_to_bar']
dpmin_reached, idx_dpmin = min(list(zip(netFeedReflux.res_heat_consumer['dp_bar'], netFeedReflux.res_heat_consumer.index)), key = lambda x: x[0])

buildIDD = netFeedReflux.heat_consumer['build_ID'].at[idx_dpmin]

print(f'\n### dpmin reached = {dpmin_reached} bar. Target value is {dpmin_target} bar; at heat_consumer with index {idx_dpmin} and build_ID = {buildIDD} ###\n')

###
pmin_reached = min(netFeedReflux.res_junction['p_bar'])
print(f'\n### pmin reached is {pmin_reached} bar. Target value is {pmin_target} bar.')

###
import pandas as pd
pd.cut(netFeedReflux.res_heat_consumer['t_to_k'], np.arange(tReflux + 273.15 - 2, tReflux + 273.15 + 2, 0.25)).value_counts(normalize = True)

# netFeedReflux.pipe['mdot_from_kg_per_s'] = netFeedReflux.res_pipe['mdot_from_kg_per_s']

# gp.GeoDataFrame(netFeedReflux.pipe, geometry = 'geometry').set_crs(cs).to_file(flp_out / Path('netFeedReflux_dimensioned_manually_final_controlled_pipes.gpkg'), driver = 'GPKG')

# netFeedReflux.heat_consumer['t_from_k'] = netFeedReflux.res_heat_consumer['t_from_k']
# netFeedReflux.heat_consumer['t_to_k'] = netFeedReflux.res_heat_consumer['t_to_k']

# netFeedReflux.heat_consumer['dp_bar'] = netFeedReflux.res_heat_consumer['p_from_bar'] - netFeedReflux.res_heat_consumer['p_to_bar']

# gp.GeoDataFrame(netFeedReflux.heat_consumer, geometry = 'geometry').set_crs(cs).to_file(flp_out / Path('netFeedReflux_dimensioned_manually_final_controlled_heat_consumers.gpkg'), driver = 'GPKG')

# netFeedReflux.junction['t_k'] = netFeedReflux.res_junction['t_k']
# netFeedReflux.junction['p_bar'] = netFeedReflux.res_junction['p_bar']

# gp.GeoDataFrame(netFeedReflux.junction, geometry = 'geometry').set_crs(cs).to_file(flp_out / Path('netFeedReflux_dimensioned_manually_final_controlled_junctions.gpkg'), driver = 'GPKG')