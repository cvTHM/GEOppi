# %% Imports


import geopandas as gp
import pandapipes as ppi
import numpy as np
import os

import rasterio
from pathlib import Path

from auxFunctions import (get_dict_from_aggregated_groups)
from create_network_topology import (create_ppi_network_from_gdf,)
from suitable_network_routing import (MST_gdf_subset,)
from dimension_network_pipes import (update_ppi_results, assign_insulation_type, hydraulicDimensioningNetwork_singleLoadPoint, assign_nominal_widths_manually)
from internal_auxFunctions import (extract_FluidProperties_ppi, transfer_LoadPoint_ppi, implement_controllers)
from pandapipes.control import run_control

from pipe_characteristics import load_hydraulic_pipe_characteristics


# %% Load data

# Define coordinate system
cs = 'EPSG:25832'

# Define network name
name = 'exampleNetwork'

# Select whather network shall be saved as pickle file
save_network = True

### Filepaths
# Input path
flp = Path(r'exampleNetwork')

# output path
flp_out = Path(r'exampleNetwork\results')

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
ppi.set_user_pf_options(net = netFeed, mode="hydraulics", friction_model = friction_model, quit_on_inconsistency_connectivity=True)

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
    user_pf_options = {'mode':'sequential', 'friction_model':friction_model,'quit_on_inconsistency_connectivity':True},
    respectHeight=False,
    fluidProperties={'cp':cp_fluid, 'rho':rho_fluid, 'nu':nu_fluid},
    elongationFactorPipes=factor_l,
    resultAttributes = {'pipe':['Pa_per_m', 'mdot_from_kg_per_s'], 'junction':['t_k', 'p_bar'], 'heat_consumer':['t_from_k', 't_to_k']}
)


ppi.to_pickle(net = netFeedReflux, filename = str(flp_out / Path('netFeedReflux_dimensioned_manually_final.p')))
gp.GeoDataFrame(netFeedReflux.pipe, geometry = 'geometry').to_file(flp_out / Path('netFeedReflux_dimensioned_manually_final_pipes.gpkg'), driver = 'GPKG')