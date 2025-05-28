# %% Imports


import geopandas as gp
import pandapipes as ppi
import numpy as np
import os

import rasterio
from pathlib import Path

from auxFunctions import (get_dict_from_aggregated_groups)
from create_network_topology import (create_ppi_network_from_gdf, MST_gdf_subset)
from dimension_network_pipes import (update_ppi_results, assign_insulation_type, hydraulicDimensioningNetwork_singleLoadPoint, assign_nominal_widths_manually)
from internal_auxFunctions import (extract_FluidProperties_ppi, transfer_LoadPoint_ppi, implement_controllers)
from pandapipes.control import run_control

from pipe_characteristics import load_hydraulic_pipe_characteristics


# %% Load data

# Define coordinate system
cs = 'EPSG:25832'

# Define network name
name = 'Herborn'

# Select whather network shall be saved as pickle file
save_network = True

### Filepaths
# Input path
flp = Path(r'C:\Users\const\Desktop\Herborn_test')

# output path
flp_out = Path(r'C:\Users\const\Desktop\Herborn_test\results')

if not os.path.exists(flp_out):
    os.makedirs(flp_out)
    print('### New directory for saving pandapipes network is created! ###')
            
# output name
fln_out = 'network_' + name

# Load pipes
pipes = gp.read_file(flp / Path('networkRaw.gpkg')).to_crs(cs)

# Load Buildings layer
buildings = gp.read_file(flp / Path('buildings_networkCalc.gpkg')).to_crs(cs)
buildings = buildings[buildings['demand_use_th'] > 1000].reset_index(drop = True)

# Load valve Layer
valves = gp.read_file(flp / ('valves.gpkg')).to_crs(cs)

# Load production site layer
producers = gp.read_file(flp / Path('producers.gpkg')).to_crs(cs)

# Load elevation raster layer
# dgm = rasterio.open('C:/Users/const/Desktop/Beispielnetz_1/exampleNetwork/dgm_merged_Buseck.tif')
dgm = None#rasterio.open(flp / Path('dgm.tif'))

# Load optional point layer for creating MST
# points_MST = gp.read_file(flp / ('points_ConsumersNetworkCalc.gpkg'))
# points_MST2 = producers.geometry.centroid

# import pandas as pd
# points_MST = pd.concat((points_MST, points_MST2), axis = 0)
# pipes = MST_gdf_subset(lines = pipes, weightMST = 'demand_use_total_exp_ren_0.75_2030_climate_demo', subsetNodesMST = points_MST)

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
flh = (1700/0.78)

# Define attribute name containing information about thermal power (kW) of producers in producer models
thermalPowerAttr = 'Pth_kW'

# Define attribute name containing infotmation about distribution ('VE') and house connection pipes ('HA')
connectionTypeAttr = 'connectionType'

# Define naming convention for base load producers in heat_consumer component DataFrame
baseLoadProdNaming = 'prod'

# Target specific pressure loss in pipes (Pa/m, usually 60-100 Pa/m)
dp_spec = 100

# Selection of subset for nominal widths for pipes
subsetDistribution = [65, 80, 100, 125, 150, 200, 225, 250, 300, 350]
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


### Update results with new insulation types
netFeed = update_ppi_results(
    net = netFeed,
    user_pf_options = {'mode':'hydraulics', 'friction_model':friction_model,'quit_on_inconsistency_connectivity':True},
    respectHeight=False,
    fluidProperties={'cp':cp_fluid, 'rho':rho_fluid, 'nu':nu_fluid},
    elongationFactorPipes=factor_l,
    resultAttributes = {'pipe':['Pa_per_m', 'mdot_from_kg_per_s'], 'junction':['t_k', 'p_bar']}
)

# %%

ppi.to_pickle(netFeed, filename = str(flp_out / name) + '_netFeed_dimensioned.p')
netFeed = ppi.from_pickle(filename = str(flp_out / name) + '_netFeed_dimensioned.p')
gp.GeoDataFrame(netFeed.pipe).set_crs(cs).to_file(str(flp_out / name) + '_netFeed_pipes_dimensioned.gpkg', driver = 'GPKG')

# netFeed_dimensioned = netFeed.deepcopy()
# netFeed = ppi.from_pickle(str(flp_out / name) + '_netFeed_dimensioned.p')

# %%
### *--- Optional: Assign nominal widths by manual selection of groups of pipes. Can be done before creating coupled feed- and reflux line network. ---*

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

# pipesNW100 = gp.read_file(r'C:\GitLab\paper_generic_network_creation\data\heatingNetwork_Herborn\results_v2\pipes_NW_100_feedLine.gpkg')
# idxsNW100 = netFeed_dimensioned.pipe[netFeed_dimensioned.pipe['ID'].isin(pipesNW100['ID'])].index

# netFeed_dimensioned = assign_nominal_widths_manually(
#     net = netFeed_dimensioned,
#     dictNW = {100:list(idxsNW100)},
#     dictNominalInnerDiameter = dictNominalInnerDiameter,
#     dictWallRoughness = dictWallRoughness,
#     nominalWidthAttr = 'nominalWidth'
# )


### Assign final heat transfer coefficient
netFeed_dimensioned.pipe = assign_insulation_type(
    df_pipes = netFeed_dimensioned.pipe, 
    colNameNominalWidth = 'nominalWidth', 
    network_type = networkType, 
    insulation_type = insulationType if networkType == 'KMR' else None
)

netFeed_dimensioned = update_ppi_results(
    net = netFeed_dimensioned,
    user_pf_options = {'mode':'hydraulics', 'friction_model':friction_model,'quit_on_inconsistency_connectivity':True},
    respectHeight=False,
    fluidProperties={'cp':cp_fluid, 'rho':rho_fluid, 'nu':nu_fluid},
    elongationFactorPipes=factor_l,
    resultAttributes = {'pipe':['Pa_per_m', 'mdot_from_kg_per_s'], 'junction':['t_k', 'p_bar']}
)

ppi.to_pickle(netFeed_dimensioned, filename = str(flp_out / name) + '_netFeed_dimensioned_manually_new.p')
gp.GeoDataFrame(netFeed.pipe).set_crs(cs).to_file(str(flp_out / name) + '_netFeed_pipes_dimensioned_manually_new.gpkg', driver = 'GPKG')


# %% Plotting

import matplotlib.pyplot as plt

fig, ax = plt.subplots(2, figsize = (20,15))

netFeed.pipe['nominalWidth'].hist(ax = ax[0])
netFeed_dimensioned.pipe['nominalWidth'].hist(ax = ax[1])


gdf = gp.GeoDataFrame(netFeed_dimensioned.pipe, geometry = 'geometry').set_crs(cs)
gdf['mdot_kg_per_s'] = abs(gdf['mdot_from_kg_per_s'])

fig, ax1 = plt.subplots(3, figsize = (20,20))

styles = {'linewidth':5}
xfont = {'fontsize': 20, 'weight':'normal'}  # Font for axes
tfont = {'fontsize': 20}  # Font for title

ax1[0].set_title('Nominal widths', **tfont)
gdf[(gdf['in_service'] ==True) & (gdf[connectionTypeAttr] == 'distribution')].plot(
    ax = ax1[0], 
    column = 'nominalWidth', 
    cmap=plt.get_cmap('Paired'),
    categorical = True,
    legend_kwds = {'title':'Nominal width', 
                   'loc':'lower right',
                   'markerscale':1, 
                    'title_fontsize':'xx-large', 
                    'fontsize':'xx-large'},
    legend = True,
    **styles)
leg = ax1[0].get_legend()
leg.set_bbox_to_anchor((1.1,0,0.2,0.5))

ax1[1].set_title('Spec. pressure loss', **tfont)
gdf[(gdf['in_service'] ==True) & (gdf[connectionTypeAttr] == 'distribution')].plot(
    ax = ax1[1], 
    column = 'Pa_per_m', 
    cmap=plt.get_cmap('magma'),
    vmin = 0,
    vmax = 200, 
    legend = True,
    **styles)


ax1[2].set_title('Affiliation with loop', **tfont)
gdf[(gdf['in_service'] ==True) & (gdf[connectionTypeAttr] == 'distribution')].plot(
    ax = ax1[2], 
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
leg = ax1[2].get_legend()
leg.set_bbox_to_anchor((1.1,0,0.2,0.5))




# %% Create complete Pandapipes network from pre-dimensioned feed line network

# netFeed_dimensioned = ppi.from_pickle(str(flp_out / name) + '_netFeed_dimensioned_manually.p')

# cp_fluid, rho_fluid, nu_fluid, g = extract_FluidProperties_ppi(net = netFeed_dimensioned, t = 0.5 * (tFeed + tReflux))

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
ppi.set_user_pf_options(net = netFeedReflux, mode="sequential", friction_model = friction_model, quit_on_inconsistency_connectivity=True, reset = True)

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

ppi.to_pickle(netFeedReflux, filename = str(flp_out / name) + '_netFeedReflux_dimensioned_manually_new.p')

gp.GeoDataFrame(netFeedReflux.pipe).set_crs(cs).to_file(flp_out / 'networkFeedReflux_dimensioned_manually_pipes_new.gpkg', driver = 'GPKG')
gp.GeoDataFrame(netFeedReflux.junction).set_crs(cs).to_file(flp_out / 'networkFeedReflux_dimensioned_manually_junctions_new.gpkg', driver = 'GPKG')

# netFeedReflux = ppi.from_pickle(str(flp_out) + '/exampleNetwork_netFeedReflux_circPumps.p')

# %% *--- Controller creation and exemplified application ---*

netFeedReflux = ppi.from_pickle(str(flp_out / name) + '_netFeedReflux_dimensioned_manually.p')
cp_fluid, rho_fluid, nu_fluid, g = extract_FluidProperties_ppi(net = netFeedReflux, t = 0.5 * (tFeed + tReflux))

# Initialize
netFeedReflux.heat_consumer['Pth_kW'] = netFeedReflux.heat_consumer['demand_use_th'] / 1700 * 0.1

# Output DataFrames
import pandas as pd
results = pd.DataFrame(
    columns = ['partialLoad', 'Pth_demand_kW', 'Pth_loss_kW', 'Pth_prod1_kW', 'Pth_prod2_kW', 'Pth_prod3_kW', 'P_pump_prod1_kW', 'P_pump_prod2_kW', 'P_pump_prod3_kW', 'dp_pump_prod1_bar', 'dp_pump_prod2_bar', 'dp_pump_prod3_bar', 'p_supply_prod1_bar', 'p_supply_prod2_bar', 'p_supply_prod3_bar', 'buildID_dpmin']
)

results['partialLoad'] = results['partialLoad'].astype(float)
results['Pth_demand_kW'] = results['Pth_demand_kW'].astype(float)
results['Pth_loss_kW'] = results['Pth_loss_kW'].astype(float)
results['Pth_prod1_kW'] = results['Pth_prod1_kW'].astype(float)
results['Pth_prod2_kW'] = results['Pth_prod2_kW'].astype(float)
results['Pth_prod3_kW'] = results['Pth_prod3_kW'].astype(float)
results['P_pump_prod1_kW'] = results['P_pump_prod1_kW'].astype(float)
results['P_pump_prod2_kW'] = results['P_pump_prod2_kW'].astype(float)
results['P_pump_prod3_kW'] = results['P_pump_prod3_kW'].astype(float)
results['dp_pump_prod1_bar'] = results['dp_pump_prod1_bar'].astype(float)
results['dp_pump_prod2_bar'] = results['dp_pump_prod2_bar'].astype(float)
results['dp_pump_prod3_bar'] = results['dp_pump_prod3_bar'].astype(float)
results['p_supply_prod1_bar'] = results['p_supply_prod1_bar'].astype(float)
results['p_supply_prod2_bar'] = results['p_supply_prod2_bar'].astype(float)
results['p_supply_prod3_bar'] = results['p_supply_prod3_bar'].astype(float)
results['buildID_dpmin'] = results['buildID_dpmin'].astype(int)


# Initialize output array for saving heat consumer results
# Rows: heat consumers
# 2nd dim: array elements from 19 variation runs
# 3rd dim: contains dp values, Tsupply values and Treflux values
resultsConsumers = np.zeros((len(netFeedReflux.heat_consumer), 10, 3))

try:
    if len(netFeedReflux.controller) > 0:
        netFeedReflux.controller['in_service'] = False
        print(f'\n### Controllers in network.controller are set out of service. ###')
except:
    pass

# Define control targets
circ_pump_pressure_idx = 0

dpmin_target = 1.2 # bar
pmin_target = 1.5 # bar

# Create controllers 
netFeedReflux = implement_controllers(
    net = netFeedReflux,
    pminCtrlDict = {
        'create':True, 
        'pmin_target':pmin_target, 
        'circ_pump_pressure_idx':circ_pump_pressure_idx,
        'abs_tol':0.1,
        'order':1,
        'level':5
        },
    dpminCtrlDict = {
        'create':True, 
        'dpmin_target':dpmin_target, 
        'circ_pump_pressure_idx':circ_pump_pressure_idx,
        'index':1,
        'abs_tol':0.1,
        'order':0,
        'level':5
        },
    TRefluxHeatConsumerCtrlDict = {
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
        'priority_list':[0, 1],
        'abs_tol':50000,
        'order':5,
        'level':5
        }
    )

print(netFeedReflux.controller)

print(f'\nLänge der Rohre (gesamt, geometry) = {np.nansum(gp.GeoDataFrame(netFeedReflux.pipe).geometry.length):.4} m')
print(f'\nLänge der Rohre (gesamt, length) = {np.nansum(netFeedReflux.pipe["length_km"]*1000):.4} m')

# %%

sf = 0.78

# Loop for variation of thermal power demand
for nn, factorize in enumerate(np.linspace(0.1, 1, 10)):#enumerate([0.1, 0.5, 1]): # enumerate(np.linspace(0.1, 1, 19)):

    ## Initializations
    # Storing intermediate values
    resultRows = []

    original_mdots = netFeedReflux.heat_consumer['controlled_mdot_kg_per_s'].copy()
    original_mdots_circ_pump_mass = netFeedReflux.circ_pump_mass['mdot_flow_kg_per_s'].copy()
    Pth_target = netFeedReflux.circ_pump_mass['Pth_kW'].copy()

    netFeedReflux.heat_consumer['Pth_kW'] = netFeedReflux.heat_consumer['demand_use_th'] / 1700 * sf * factorize

    netFeedReflux = transfer_LoadPoint_ppi(
        net = netFeedReflux,
        tFeed = tFeed + 273.15,
        tReflux = tReflux + 273.15,
        text_pipes = 280,
        pNetwork = pNetwork,
        flh = flh,
        heatingDemandAttr = heatingDemandAttr,
        thermalPowerAttr = thermalPowerAttr,
        thermalPowerAttrConsumers = 'Pth_kW',
        baseLoadProdNaming = baseLoadProdNaming
    )

    # *--- Running calculations with controllers implemented ---*

    # Run pipeflow with consideration of controllers
    run_control(net = netFeedReflux,  max_iter = 50)

    # Start another pipeflow to compare results (identical)
    ppi.pipeflow(net = netFeedReflux, mode = 'sequential', friction_model = friction_model)

    #
    # *--- Printing ---*
    ###
    meanT_base = (netFeedReflux.res_circ_pump_mass['t_to_k'].values + netFeedReflux.res_circ_pump_mass['t_from_k'].values)/2
    Pth_reached = netFeedReflux.res_circ_pump_mass['mdot_from_kg_per_s'].values * netFeedReflux.fluid.get_heat_capacity(meanT_base) * (netFeedReflux.res_circ_pump_mass['t_to_k'].values - netFeedReflux.res_circ_pump_mass['t_from_k'].values)

    print(f'Thermal power output in res DF is {Pth_reached/(1e03)} kW. Entered target value is\n {Pth_target} kW')

    meanT_peak = (netFeedReflux.res_circ_pump_pressure['t_to_k'].values + netFeedReflux.res_circ_pump_pressure['t_from_k'].values) / 2
    Pth_peakLoad = netFeedReflux.res_circ_pump_pressure.loc[0, 'mdot_from_kg_per_s'] * netFeedReflux.fluid.get_heat_capacity(meanT_peak/2+273.15) * (netFeedReflux.res_circ_pump_pressure.loc[0, 't_to_k'] - netFeedReflux.res_circ_pump_pressure.loc[0, 't_from_k'])

    print(f'\n### Thermal output power of peak load producer is\n {Pth_peakLoad/(1e03)} kW. ###\n')

    difference_mdots_circ_pump_mass = original_mdots_circ_pump_mass - netFeedReflux.res_circ_pump_mass['mdot_from_kg_per_s']

    print(f'\n### Changed mass flows at circ pump mass compared to initial values are\n{difference_mdots_circ_pump_mass} ###\n')


    ###
    # Calculate intermediate results
    Pth_demand_kW = np.nansum(netFeedReflux.res_heat_consumer['qext_w']) / 1e03
    Pth_loss_kW = Pth_demand_kW - (np.nansum(Pth_peakLoad) + np.nansum(Pth_reached)) / 1e03

    P_pump_prod1 = (netFeedReflux.res_circ_pump_pressure.loc[0,'mdot_from_kg_per_s'] / rho_fluid * (netFeedReflux.res_circ_pump_pressure.loc[0, 'p_to_bar'] - netFeedReflux.res_circ_pump_pressure.loc[0, 'p_from_bar']) * 1e05) / 1e03

    P_pump_prod2 = (netFeedReflux.res_circ_pump_mass.loc[0, 'mdot_from_kg_per_s'] / rho_fluid * (netFeedReflux.res_flow_control.loc[0, 'p_to_bar'] - netFeedReflux.res_circ_pump_mass.loc[0, 'p_from_bar']) * 1e05) / 1e03

    P_pump_prod3 = (netFeedReflux.res_circ_pump_mass.loc[1, 'mdot_from_kg_per_s'] / rho_fluid * (netFeedReflux.res_flow_control.loc[1, 'p_to_bar'] - netFeedReflux.res_circ_pump_mass.loc[1, 'p_from_bar']) * 1e05) / 1e03

    dp_pump_prods = [
        netFeedReflux.res_circ_pump_pressure.loc[0, 'p_to_bar'] - netFeedReflux.res_circ_pump_pressure.loc[0, 'p_from_bar'],
        netFeedReflux.res_flow_control.loc[0, 'p_to_bar'] - netFeedReflux.res_circ_pump_mass.loc[0, 'p_from_bar'],
        netFeedReflux.res_flow_control.loc[1, 'p_to_bar'] - netFeedReflux.res_circ_pump_mass.loc[1, 'p_from_bar']
        ]
    
    p_supply_prods = [
        netFeedReflux.res_circ_pump_pressure.loc[0, 'p_to_bar'],
        netFeedReflux.res_flow_control.loc[0, 'p_to_bar'],
        netFeedReflux.res_flow_control.loc[1, 'p_to_bar']
    ]

    netFeedReflux.res_heat_consumer['dp_bar'] = netFeedReflux.res_heat_consumer['p_from_bar'] - netFeedReflux.res_heat_consumer['p_to_bar']
    dpmin_reached, idx_dpmin = min(list(zip(netFeedReflux.res_heat_consumer['dp_bar'], netFeedReflux.res_heat_consumer.index)), key = lambda x: x[0])

    buildIDD = netFeedReflux.heat_consumer['build_ID'].at[idx_dpmin]

    print(f'\n### dpmin reached = {dpmin_reached} bar. Target value is {dpmin_target} bar; at heat_consumer with index {idx_dpmin} and build_ID = {buildIDD} ###\n')


    ###
    pmin_reached = min(netFeedReflux.res_junction['p_bar'])
    print(f'\n### pmin reached is {pmin_reached} bar. Target value is {pmin_target} bar.')


    ###
    difference_mdot = original_mdots - netFeedReflux.res_heat_consumer['mdot_from_kg_per_s']
    difference_mdot_inactive_consumers = difference_mdot[netFeedReflux.heat_consumer[netFeedReflux.heat_consumer['qext_w'].isna()].index]


    print(f'\n### Changed mass flows at heat consumers compared to initial values  are\n{difference_mdot}\n')

    print(f'\n### Changed mass flows at inactive heat consumers compared to initial values  are\n{difference_mdot_inactive_consumers}')

    nIterations = netFeedReflux.controller.loc[3, 'object'].iterations
    print(f'\n### Controller for reflux temperatures needed {nIterations} iteration. ###\n')

    # Write result data to output DataFrame  
    resultRows.append(
        {
            'partialLoad':factorize,
            'Pth_demand_kW':Pth_demand_kW,
            'Pth_loss_kW':Pth_loss_kW,
            'Pth_prod1_kW':Pth_peakLoad[0]/1e03,
            'Pth_prod2_kW':Pth_reached[0]/1e03,
            'Pth_prod3_kW':Pth_reached[1]/1e03,
            'P_pump_prod1_kW':P_pump_prod1,
            'P_pump_prod2_kW':P_pump_prod2,
            'P_pump_prod3_kW':P_pump_prod3,
            'dp_pump_prod1_bar':dp_pump_prods[0],
            'dp_pump_prod2_bar':dp_pump_prods[1],
            'dp_pump_prod3_bar':dp_pump_prods[2],
            'p_supply_prod1_bar':p_supply_prods[0],
            'p_supply_prod2_bar':p_supply_prods[1],
            'p_supply_prod3_bar':p_supply_prods[2],
            'buildID_dpmin':buildIDD          
        }
    )

    results = pd.concat((results, pd.DataFrame(resultRows)), ignore_index = True)

    # Write data for heat consumers to numpy array
    resultsConsumers[:, nn, 0] = np.array(netFeedReflux.res_heat_consumer['p_from_bar']-netFeedReflux.res_heat_consumer['p_to_bar'])
    resultsConsumers[:, nn, 1] = np.array(netFeedReflux.res_heat_consumer['t_from_k'])
    resultsConsumers[:, nn, 2] = np.array(netFeedReflux.res_heat_consumer['t_to_k'])

    # Saving data
    # ppi.to_pickle(netFeedReflux, filename = str(flp_out / name) + '_netFeedReflux_dimensioned_manually_controlled.p')

    gp.GeoDataFrame(netFeedReflux.pipe).set_crs(cs).to_file(flp_out / 'pipes_heatDemand_variation.gpkg',layer = f'thermalPowerfactor{factorize:.1f}', driver = 'GPKG')
    gp.GeoDataFrame(netFeedReflux.junction).set_crs(cs).to_file(flp_out / 'networkFeedReflux_dimensioned_manually_controlled_junctions.gpkg', driver = 'GPKG')

    gp.GeoDataFrame(netFeedReflux.heat_consumer).set_crs(cs).to_file(flp_out / 'heat_consumers_heatDemand_variation.gpkg', layer = f'thermalPowerfactor{factorize:.1f}', driver = 'GPKG')


# results.to_excel(flp_out / Path(f'results_heatDemand_variation_sf{sf:.2}.xlsx'))
# np.save(flp_out / Path(f'results_heatDemand_variation_heatConsumers_sf{sf:.2}.npy'), resultsConsumers)

print(f'\nLänge der Rohre (gesamt, geometry) = {np.nansum(gp.GeoDataFrame(netFeedReflux.pipe).geometry.length):.4} m')
print(f'\nLänge der Rohre (gesamt, length) = {np.nansum(netFeedReflux.pipe["length_km"]*1000):.4} m')

# %%
import pandas as pd
import numpy as np

# results = pd.read_excel(flp_out / Path(f'results_heatDemand_variation_sf0.78.xlsx'))
# resultsConsumers = np.load(flp_out / Path(f'results_heatDemand_variation_heatConsumers_sf0.78.npy'))

# Plot appearance
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

xfont = {'fontname': 'Times New Roman', 'fontsize':9}
legendfont = {'family': 'Times New Roman', 'size': 9}

# Data
data = results.copy()
data[['Pth_prod2_kW', 'Pth_prod3_kW', 'Pth_prod1_kW']] *= 1/1e03 # MW
data['rel_th_loss'] = abs(data['Pth_loss_kW'] / data['Pth_demand_kW'])
data['rel_pumping_power'] = data[['P_pump_prod1_kW', 'P_pump_prod2_kW', 'P_pump_prod3_kW']].sum(axis = 1) / data['Pth_demand_kW']
data['rel_pumping_power_el_incldpProducers'] = \
    (
        data['P_pump_prod1_kW'] / data['dp_pump_prod1_bar'] * (data['dp_pump_prod1_bar'] + 0.8) + \
        data['P_pump_prod2_kW'] / data['dp_pump_prod2_bar'] * (data['dp_pump_prod2_bar'] + 0.8) + \
        data['P_pump_prod3_kW'] / data['dp_pump_prod3_bar'] * (data['dp_pump_prod3_bar'] + 0.8)
    ) / 0.7 / data['Pth_demand_kW']

# data[['P_pump_prod1_kW', 'P_pump_prod2_kW', 'P_pump_prod3_kW']].sum(axis = 1) / data['Pth_demand_kW']

# Create plot

# PLot size for single column layout is (3.22, 2)

fig, (ax, ax_bottom) = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [0.8, 1]})
fig.set_size_inches(6.44, 4, forward=True)
fig.set_dpi(300)

ax.stackplot(
    data['partialLoad'], 
    data['Pth_prod2_kW'], 
    data['Pth_prod3_kW'], 
    data['Pth_prod1_kW'],
    labels = ['$P_{th,bl1}$', '$P_{th,bl2}$', '$P_{th,pl}$'],
    colors = [(255/255, 0/255, 0/255), (0/255, 0/255, 50/255), (255/255, 255/255, 0/255)],
    alpha = 0.4
    )

axXticks = np.unique(data['partialLoad'].values.round(1))
ax.set_ylim([0, 30])
ax.tick_params(axis='both', which='major', labelsize=8)
ax.set_ylabel(r'Thermal power supply (MW)', **xfont)
ax.set_xticks(axXticks, labels = [str(tick) for tick in axXticks], **legendfont)
ax.set_yticks(ax.get_yticks(), labels = [str(tick) for tick in ax.get_yticks()], **legendfont)

ax.grid('minor')

# Ax2
ax2 = ax.twinx()
line_c, = ax2.plot(data['partialLoad'], data['rel_th_loss'], color='red', label=r'$\frac{P_{th,loss}}{P_{th,dem}}$')
line_d, = ax2.plot(data['partialLoad'], data['rel_pumping_power_el_incldpProducers'], color='red', linestyle = '-.', label=r'$\frac{P_{el,pump}}{P_{th,dem}}$')

ax2.set_ylim([0, 0.15])
ax2.set_yticks(ax2.get_yticks(), labels = [str(tick) for tick in ax2.get_yticks()], **legendfont)
ax2.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

ax2.set_ylabel('Rel. thermal loss,\nRel. el. pumping power', **xfont)

## Bottom axes
# Plot pressure differences at heat consumers
scat_dp = ax_bottom.scatter(
    x = np.array(data['partialLoad']).repeat(np.shape(resultsConsumers)[0]),
    y = np.array(resultsConsumers[:, :, 0]).flatten('F'),
    marker = '_',
    facecolors = 'green',
    s = 120,
    linewidths = 0.2,
    alpha = 1,
    label = '$\Delta p$',
    zorder = 11
)


# Plot supply and return temperatures at heat consumers
ax_bottom2 = ax_bottom.twinx()

# Supply temperature
scat_supply = ax_bottom2.scatter(
    x = np.array(data['partialLoad']).repeat(np.shape(resultsConsumers)[0]),
    y = np.array(resultsConsumers[:, :, 1]).flatten('F') - 273.15,
    marker = 'o',
    facecolors = 'none',
    edgecolors = 'black',
    linewidths = 0.25,
    alpha = 1,
    label = '$T_{supply}$',
    zorder = 10
)

## Return temperature
scat_return = ax_bottom2.scatter(
    x = np.array(data['partialLoad']).repeat(np.shape(resultsConsumers)[0]),
    y = np.array(resultsConsumers[:, :, 2]).flatten('F') - 273.15,
    marker = '^',
    edgecolors = 'orange',
    facecolors = 'none',
    linewidths = 0.5,
    alpha = 1,
    label = '$T_{return}$',
    zorder = 9
)

ax_bottom.set_xlabel(r'Rel. thermal demand load $\frac{P_{th,dem}}{P_{th,dem,max}}$', **xfont)
ax_bottom.set_xlim([0.05, 1.05])
ax_bottom.set_ylim([0.2, 6])
ax_bottomYticks = np.arange(7)

ax_bottom2.set_ylim([20, 80])

ax_bottom.set_ylabel(r'Pressure difference $\Delta p$ (bar)', **xfont)
ax_bottom.set_xticks(ax_bottom.get_xticks(), labels = [str(tick) for tick in ax_bottom.get_xticks()], **legendfont)
ax_bottom.set_yticks(ax_bottomYticks, labels = [str(tick) for tick in ax_bottomYticks], **legendfont)
ax_bottom.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax_bottom.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

ax_bottom2.set_yticks(ax_bottom2.get_yticks(), labels = [str(tick) for tick in ax_bottom2.get_yticks()], **legendfont)
ax_bottom2.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

ax_bottom2.set_ylabel(r'Supply / Return temperature (°C)', **xfont)

ax_bottom.grid('minor')

# Legende kombinieren
lines_labels = ax.get_legend_handles_labels()
lines_labels2 = ax2.get_legend_handles_labels()
scat_labels2 = ax_bottom2.get_legend_handles_labels()
scat_labels = ax_bottom.get_legend_handles_labels()

handles = lines_labels[0] + [line_c, line_d] + [scat_supply, scat_return] + [scat_dp]
labels = lines_labels[-1] + lines_labels2[-1] + scat_labels2[-1] + scat_labels[-1]

ax_bottom2.legend(handles, labels, loc='upper left', bbox_to_anchor = (1.125, 1.5), ncol = 1, prop = legendfont)

# fig.savefig(flp_out / Path('fig_heatDemandLoad_variation.png'), dpi = 300, bbox_inches = 'tight')


# %%

# netFeedReflux.res_heat_consumer['t_to_k'].value_counts(bins = [317, 320, 325, 330, 335, 340, 345, 350], normalize = True)