# %%
import geopandas as gp
import pandapipes as ppi
from pandapipes.control import run_control

from internal_auxFunctions import implement_controllers

from internal_auxFunctions import (extract_FluidProperties_ppi, transfer_LoadPoint_ppi)



tFeed = 273.15 + 10
tReflux = 273.15 + 5


# %% Implement control strategy for substation with HP and defined COP at heat consumers

# netFeedReflux = ppi.from_pickle(str(flp_out / name) + '_netFeedReflux_dimensioned_manually_new.p')
netFeedReflux = ppi.from_pickle(r'exampleNetwork/results/exampleNetwork_netFeedReflux_dimensioned_manually.p')

cp_fluid, rho_fluid, nu_fluid, g = extract_FluidProperties_ppi(net = netFeedReflux, t = 0.5 * (tFeed + tReflux))

ppi.drop_elements_at_junctions(netFeedReflux, netFeedReflux.flow_control[['from_junction', 'to_junction']].values.flatten())



# %%


# Initialize
netFeedReflux.heat_consumer['Pth_kW'] = netFeedReflux.heat_consumer['demand_use_th'] / 1700


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
    dpminCtrlDict = {
        'create':True, 
        'dpmin_target':dpmin_target, 
        'circ_pump_pressure_idx':circ_pump_pressure_idx,
        'index':1,
        'abs_tol':0.05,
        'order':0,
        'level':5
        },
    COPConversionQdemConsumerCtrlDict = {
        'create':True,
        'abs_tol':1,
        'order':0,
        'level':4,
        'Tsink':273.15+60,
        'efficiency':0.5
    }
    )

print(netFeedReflux.controller)

# %%

netFeedReflux.heat_consumer['Pth_kW'] = netFeedReflux.heat_consumer['demand_use_th'] / 1700

netFeedReflux = transfer_LoadPoint_ppi(
        net = netFeedReflux,
        tFeed = tFeed + 273.15,
        tReflux = tReflux + 273.15,
        text_pipes = 280,
        pNetwork = 2,
        flh = 1700,
        heatingDemandAttr = 'demand_use_th',
        thermalPowerAttr = 'Pth_kW',
        thermalPowerAttrConsumers = 'Pth_kW',
        baseLoadProdNaming = 'prod'
    )

# *--- Running calculations with controllers implemented ---*

# Run pipeflow with consideration of controllers
run_control(net = netFeedReflux,  max_iter = 50)