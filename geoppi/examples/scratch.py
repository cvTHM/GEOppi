# %%
import pandapipes as ppi
import geoppi


net = ppi.from_pickle(r'Z:\02_Mitarbeiter\Constantin_Voelzel\THESA_THM\98_Konferenzen\NEIS_2025_Hamburg\data\heatingNetwork_Herborn\results_v2\Herborn_netFeedReflux_dimensioned_manually.p')

# %%

len(net.heat_consumer)
# %%

net.pipe.loc[net.pipe['connectionType']=='distribution', 'length_km'].sum()

# %%

net.controller
circ_pump_pressure_idx = 0

# Create controllers 
net2 = geoppi.implement_controllers(
    net = net,
    drop_all = True,
    pminCtrlDict = {
        'create':True, 
        'pmin_target':1.5, 
        'circ_pump_pressure_idx':circ_pump_pressure_idx,
        'abs_tol':0.1,
        'order':1,
        'level':5
        },
    dpminCtrlDict = { # Lowest control order in level 5 = highest priority in level 5
        'create':True, 
        'dpmin_target':0.8, 
        'circ_pump_pressure_idx':circ_pump_pressure_idx,
        'index':1,
        'abs_tol':0.1,
        'order':0,
        'level':5
        },
    TRefluxHeatConsumerCtrlDict = { # Lowest control level = highest priority
        'create':True,
        'T_target':45 + 273.15,
        'min_mdot':0.015,
        'min_dT':3,
        'abs_tol':0.1,
        'order':0,
        'level':-1
        },
    PthLimitedCtrlDict = {
        'create':True,
        'circ_pump_mass_idxs':net.circ_pump_mass.index, 
        'circ_pump_pressure_index':circ_pump_pressure_idx, 
        'flow_controller_idxs':net.flow_control.index,
        'Pth_target_kW':net.circ_pump_mass['Pth_kW'].values, 
        'priority_list':[0, 1],
        'abs_tol':5000,
        'order':5,
        'level':5
        }
    )

print(net2.controller)

# %%
from pandapipes.control import run_control

run_control(net = net2, mode = 'bidirectional')

# %%

net2.res_circ_pump_mass