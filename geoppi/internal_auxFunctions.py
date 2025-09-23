####################### Internal auxiliary functions ######################
################## Not meant for publication ##############################

import pandapipes as ppi
import geopandas as gp
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import os
import matplotlib.animation
from matplotlib.collections import LineCollection, PathCollection
from pathlib import Path

<<<<<<< HEAD:internal_auxFunctions.py
from customControllers import (ppi_CircPumpMassPthermalCtrl, ppi_CircPumpMassPthermalCtrl_limited, ppi_HeatConsumerSetTempCtrl, ppi_HeatConsumersMinDiffPressureCtrl, ppi_JunctionsMinAbsolutePressureCtrl, ppi_HeatConsumerCOPConversionQdem)
=======
from geoppi.customControllers import (ppi_CircPumpMassPthermalCtrl, ppi_CircPumpMassPthermalCtrl_limited, ppi_HeatConsumerSetTempCtrl, ppi_HeatConsumersMinDiffPressureCtrl, ppi_JunctionsMinAbsolutePressureCtrl, ppi_HeatConsumerCOPConversionQdem)
>>>>>>> develop_maintainer:geoppi/internal_auxFunctions.py

def implement_controllers(
        net,
        drop_all:bool = False,
        pminCtrlDict:dict = {'create':False, 'pmin_target':1, 'abs_tol':0.1, 'circ_pump_pressure_idx':0,
                              'order':1, 'level':5, 'index':None},
        dpminCtrlDict:dict = {'create':False,'dpmin_target':0.7, 'abs_tol':0.1, 'circ_pump_pressure_idx':0, 'order':0, 'level':5, 'index':None},
        PthCtrlDict:dict = {'create':False,'circ_pump_mass_idxs':None, 'flow_controller_idxs':None, 'Pth_target_kW':1000,'abs_tol':100,  'order':2, 'level':5, 'index':None},
        PthLimitedCtrlDict:dict = {'create':False,'circ_pump_mass_idxs':None, 'circ_pump_pressure_index':None, 'flow_controller_idxs':None, 'Pth_target_kW':1000, 'abs_tol':100, 'order':2, 'level':5, 'priority_list':None, 'index':None},
        TRefluxHeatConsumerCtrlDict:dict = {'create':False,'T_target':273.15+40, 'abs_tol':1, 'order':0, 'level':-1, 'min_dT':None, 'min_mdot':None, 'index':None},
        COPConversionQdemConsumerCtrlDict:dict = {'create':False,'efficiency':0.5, 'Tsink':273.15+60, 'deltaTsource':4, 'abs_tol':1, 'order':0, 'level':-1, 'index':None}
        ):
    

    """
    Function to implement common controllers in pandapipes network for district heating simulations. \n
    The aim is to ensure a network operation with (sorted with descending priority, e.g. in the course of their implementation):
        - a controlled return temperature at the heat consumers, thus varying their mass flow rate to match a specified return temperature based on their thermal extraction power
        - a controlled minimum pressure difference at the "worst supplied" heat consumer (identified by the lowest pressure difference between feed and reflux line)
        - a controlled minimum absolute pressure level in the network
        - a controlled mass flow rate for base load producers to match their specified thermal output power.\n

    Prerequisites for applying this function:\n
    - closed-loop pandapipes network (feed line and reflux line) with modelling approaches:
    1) heat_consumer models for heat consumer
    2) ONE circ_pump_const_pressure component as peak load producer defining the overall network pressure level and supplying thermal power that is not covered by all other base load producers
    3) arbitrary number of base load producers modelled as circ_pump_const_mass_flow
    
    :param net: pandapipes network object.\n
    :param drop_all: Boolöean denoting if all exiting controllers ahll be dropped first.\n
    :param pminCtrlDict: dictionary containing keys to parameterize the controller for minimum absolute pressure. Defaults to {'pmin_target':1, 'circ_pump_pressure_idx':0, 'order':0, 'level':2}\n
    :param dpminCtrlDict: dictionary containing keys to parameterize the controller for minimum pressure difference at heat consumers. Defaults to {'dpmin_target':0.7, 'circ_pump_pressure_idx':0, 'order':0, 'level':1}.\n
    :param PthCtrlDict: dictionary containing keys to parameterize the controller for defined thermal power at base load producers (circ_pump_mass components). Defaults to {'circ_pump_mass_idxs':None, 'flow_controller_idxs':None, 'Pth_target_kW':1000, 'order':0, 'level':-5}.\n
    :param TRefluxHeatConsumerCtrlDict: dictionary containing keys to parameterize the controller for defined return temperature at heat consumers. Defaults to {'T_target':273.15+40, 'order':0, 'level':0}.\n
    :return: pandapipes network
    """
    
    # Create temporary copy of network
    net_out = net.deepcopy()

    if drop_all:
        net_out.controller.drop(index = net_out.controller.index, inplace = True)

    # Assigning default values for order and level of control strategy
    ### Creating controllers

    ## Control of COP conversion at heat consumer source side thermal power
    if 'order' not in COPConversionQdemConsumerCtrlDict.keys():
        COPConversionQdemConsumerCtrlDict['order'] = 0
    if 'level' not in COPConversionQdemConsumerCtrlDict.keys():
        COPConversionQdemConsumerCtrlDict['level'] = 0
    if 'abs_tol' not in COPConversionQdemConsumerCtrlDict.keys():
        COPConversionQdemConsumerCtrlDict['abs_tol'] = 1
    if 'Tsink' not in COPConversionQdemConsumerCtrlDict.keys():
        COPConversionQdemConsumerCtrlDict['Tsink'] = 273.15+60
    if 'efficiency' not in COPConversionQdemConsumerCtrlDict.keys():
        COPConversionQdemConsumerCtrlDict['efficiency'] = 0.5
    if 'create' not in COPConversionQdemConsumerCtrlDict.keys():
        COPConversionQdemConsumerCtrlDict['create'] = True
    if 'index' not in COPConversionQdemConsumerCtrlDict.keys():
        COPConversionQdemConsumerCtrlDict['index'] = None
    if 'deltaTsource' not in COPConversionQdemConsumerCtrlDict.keys():
        COPConversionQdemConsumerCtrlDict['deltaTsource'] = 4

    if COPConversionQdemConsumerCtrlDict['create']:
        controller_COPConversion_Qdem = ppi_HeatConsumerCOPConversionQdem(
            net = net_out, 
            efficiency = COPConversionQdemConsumerCtrlDict['efficiency'], 
            proportional_gain = 0.5, 
            abs_tol = COPConversionQdemConsumerCtrlDict['abs_tol'], 
            order = COPConversionQdemConsumerCtrlDict['order'], 
            level = COPConversionQdemConsumerCtrlDict['level'], 
            index = COPConversionQdemConsumerCtrlDict['index'],
            Tsink_heatingSystem = COPConversionQdemConsumerCtrlDict['Tsink'],
            deltaTsource = COPConversionQdemConsumerCtrlDict['deltaTsource']
            ) 

    ### Control of TReflux in network heat consumers
    if 'order' not in TRefluxHeatConsumerCtrlDict.keys():
        TRefluxHeatConsumerCtrlDict['order'] = 0
    if 'level' not in TRefluxHeatConsumerCtrlDict.keys():
        TRefluxHeatConsumerCtrlDict['level'] = 0
    if 'abs_tol' not in TRefluxHeatConsumerCtrlDict.keys():
        TRefluxHeatConsumerCtrlDict['abs_tol'] = 1
    if 'min_dT' not in TRefluxHeatConsumerCtrlDict.keys():
        TRefluxHeatConsumerCtrlDict['min_dT'] = None
    if 'min_mdot' not in TRefluxHeatConsumerCtrlDict.keys():
        TRefluxHeatConsumerCtrlDict['min_mdot'] = None
    if 'create' not in TRefluxHeatConsumerCtrlDict.keys():
        TRefluxHeatConsumerCtrlDict['create'] = True
    if 'index' not in TRefluxHeatConsumerCtrlDict.keys():
        TRefluxHeatConsumerCtrlDict['index'] = None


    if TRefluxHeatConsumerCtrlDict['create']:
        controller_Treflux = ppi_HeatConsumerSetTempCtrl(
            net = net_out, 
            target_T = TRefluxHeatConsumerCtrlDict['T_target'], 
            junction_type = 'to_junction', 
            proportional_gain = 0.5, 
            abs_tol = TRefluxHeatConsumerCtrlDict['abs_tol'], 
            order = TRefluxHeatConsumerCtrlDict['order'], 
            level = TRefluxHeatConsumerCtrlDict['level'], 
            index = TRefluxHeatConsumerCtrlDict['index'],
            min_dT = TRefluxHeatConsumerCtrlDict['min_dT'],
            min_mdot = TRefluxHeatConsumerCtrlDict['min_mdot']
            )    

    ### Control of dpmin in network heat consumers
    if 'order' not in dpminCtrlDict.keys():
        dpminCtrlDict['order'] = 0
    if 'level' not in dpminCtrlDict.keys():
        dpminCtrlDict['level'] = 1
    if 'abs_tol' not in dpminCtrlDict.keys():
        dpminCtrlDict['abs_tol'] = 0.05
    if 'create' not in dpminCtrlDict.keys():
        dpminCtrlDict['create'] = True
    if 'index' not in dpminCtrlDict.keys():
        dpminCtrlDict['index'] = None


    if dpminCtrlDict['create']:

        # Check for single circ pump pressure component in network
        if hasattr(net_out, 'circ_pump_pressure'):
            if (len(net_out.circ_pump_pressure) > 1) and ('circ_pump_pressure_idx' not in dpminCtrlDict.keys()):
                print(f'\n### More than one circ_pump_pressure_component was found in the network, but no index to define the controlled component. Please provide proper index in dictionary dpminCtrlDict with key "circ_pump_pressure_idx". ###\n')
                return net
            
            elif (len(net_out.circ_pump_pressure) == 1) and ('circ_pump_pressure_idx' not in dpminCtrlDict.keys()):
                circ_pump_pressure_idx = net_out.circ_pump_pressure.index[0]

            else:
                circ_pump_pressure_idx = dpminCtrlDict['circ_pump_pressure_idx']

        controller_dpmin = ppi_HeatConsumersMinDiffPressureCtrl(
            net = net_out, 
            circPumpPressure_idx = circ_pump_pressure_idx, 
            heatConsumer_idxs = None, 
            proportional_gain = 0.5,
            target_dpmin_bar = dpminCtrlDict['dpmin_target'], 
            abs_tol = dpminCtrlDict['abs_tol'], 
            order = dpminCtrlDict['order'], 
            level = dpminCtrlDict['level'], 
            index = dpminCtrlDict['index']
            )


    # Control of pmin in network
    if 'order' not in pminCtrlDict.keys():
        pminCtrlDict['order'] = 0
    if 'level' not in pminCtrlDict.keys():
        pminCtrlDict['level'] = 2
    if 'abs_tol' not in pminCtrlDict.keys():
        pminCtrlDict['abs_tol'] = 0.1
    if 'create' not in pminCtrlDict.keys():
        pminCtrlDict['create'] = True
    if 'index' not in pminCtrlDict.keys():
        pminCtrlDict['index'] = None


    if pminCtrlDict['create']:
        # Check for single circ pump pressure component in network
        if hasattr(net_out, 'circ_pump_pressure'):
            if (len(net_out.circ_pump_pressure) > 1) and ('circ_pump_pressure_idx' not in pminCtrlDict.keys()):
                print(f'\n### More than one circ_pump_pressure_component was found in the network, but no index to define the controlled component. Please provide proper index in dictionary pminCtrlDict with key "circ_pump_pressure_idx". ###\n')
                return net
            
            elif (len(net_out.circ_pump_pressure) == 1) and ('circ_pump_pressure_idx' not in pminCtrlDict.keys()):
                circ_pump_pressure_idx = net_out.circ_pump_pressure.index[0]

            else:
                circ_pump_pressure_idx = pminCtrlDict['circ_pump_pressure_idx']

        controller_pmin = ppi_JunctionsMinAbsolutePressureCtrl(
            net = net_out, 
            circPumpPressure_idx = circ_pump_pressure_idx, 
            junction_idxs = None, 
            target_pmin_bar = pminCtrlDict['pmin_target'], 
            proportional_gain = 0.5,
            abs_tol = pminCtrlDict['abs_tol'], 
            order = pminCtrlDict['order'], 
            level = pminCtrlDict['level'], 
            index = pminCtrlDict['index']
            )
    

    # Control of thermal power at base load producers in network
    if 'order' not in PthCtrlDict.keys():
        PthCtrlDict['order'] = 0
    if 'level' not in PthCtrlDict.keys():
        PthCtrlDict['level'] = -5
    if 'abs_tol' not in PthCtrlDict.keys():
        PthCtrlDict['abs_tol'] = 100
    if 'create' not in PthCtrlDict.keys():
        PthCtrlDict['create'] = True
    if 'index' not in PthCtrlDict.keys():
        PthCtrlDict['index'] = None

    if 'circ_pump_mass_idxs' not in PthCtrlDict.keys():
        circ_pump_mass_idxs = None
    else:
        circ_pump_mass_idxs = PthCtrlDict['circ_pump_mass_idxs']

    if 'flow_controller_idxs' not in PthCtrlDict.keys():
        flow_controller_idxs = None
    else:
        flow_controller_idxs = PthCtrlDict['flow_controller_idxs']


    if PthCtrlDict['create']:
        controller_Pth = ppi_CircPumpMassPthermalCtrl(
            net = net_out, 
            circPumpMass_idx = circ_pump_mass_idxs, 
            flow_controller_idx = flow_controller_idxs, 
            target_Pth = PthCtrlDict['Pth_target_kW']*1e03, 
            abs_tol = PthCtrlDict['abs_tol'], 
            order = PthCtrlDict['order'], 
            level = PthCtrlDict['level'], 
            index = PthCtrlDict['index']
            )
        
    # Control of limited thermal power at base load producers in network
    if 'order' not in PthLimitedCtrlDict.keys():
        PthLimitedCtrlDict['order'] = 0
    if 'level' not in PthLimitedCtrlDict.keys():
        PthLimitedCtrlDict['level'] = -5
    if 'abs_tol' not in PthLimitedCtrlDict.keys():
        PthLimitedCtrlDict['abs_tol'] = 100
    if 'create' not in PthLimitedCtrlDict.keys():
        PthLimitedCtrlDict['create'] = True 
    if 'index' not in PthLimitedCtrlDict.keys():
        PthLimitedCtrlDict['index'] = None
    

    if 'priority_list' not in PthLimitedCtrlDict.keys():
        PthLimitedCtrlDict['priority_list'] = list(net.circ_pump_mass.index)
        print(f'\n### Attentaion: no priority list is provided in dictionary PthLimitedCtrlDict. Alternatively, the index of the circ_pump_mass components is used as priorisation. ###\n')

    if 'circ_pump_mass_idxs' not in PthLimitedCtrlDict.keys():
        circ_pump_mass_idxs = net.circ_pump_mass.index
    else:
        circ_pump_mass_idxs = PthLimitedCtrlDict['circ_pump_mass_idxs']

    if 'flow_controller_idxs' not in PthLimitedCtrlDict.keys():
        flow_controller_idxs = net.flow_control.index
    else:
        flow_controller_idxs = PthLimitedCtrlDict['flow_controller_idxs']

    # Create controller
    if PthLimitedCtrlDict['create']:
        # Check for single circ pump pressure component in network
        if hasattr(net_out, 'circ_pump_pressure'):
            if (len(net_out.circ_pump_pressure) > 1) and ('circ_pump_pressure_idx' not in PthLimitedCtrlDict.keys()):
                print(f'\n### More than one circ_pump_pressure_component was found in the network, but no index to define the controlled component. Please provide proper index in dictionary PthLimitedCtrlDict with key "circ_pump_pressure_idx". ###\n')
                return net
            
            elif (len(net_out.circ_pump_pressure) == 1) and ('circ_pump_pressure_idx' not in PthLimitedCtrlDict.keys()):
                circ_pump_pressure_idx = net_out.circ_pump_pressure.index[0]

            else:
                circ_pump_pressure_idx = PthLimitedCtrlDict['circ_pump_pressure_idx']
        
        # Iterate through all existing circ_pump_mass components and create controllers for each of them
        for n, idx in enumerate(circ_pump_mass_idxs):

            index_counter = int(PthLimitedCtrlDict['index'] + PthLimitedCtrlDict['priority_list'][n] - min(PthLimitedCtrlDict['priority_list'])) if PthLimitedCtrlDict['index'] is not None else PthLimitedCtrlDict['index']

            flow_controller_idx = flow_controller_idxs[n]

            controller_Pth_limited = ppi_CircPumpMassPthermalCtrl_limited(
                net = net_out,
                target_Pth = PthLimitedCtrlDict['Pth_target_kW'][n]*1e03,
                circPumpMass_idx = idx,
                flow_controller_idx = flow_controller_idx,
                circPumpPressure_idx = circ_pump_pressure_idx,
                proportional_gain = 0.9,
                priority_list = PthLimitedCtrlDict['priority_list'],
                priority = PthLimitedCtrlDict['priority_list'][n],
                abs_tol = PthLimitedCtrlDict['abs_tol'], 
                order = int(PthLimitedCtrlDict['order'] + PthLimitedCtrlDict['priority_list'][n]),
                level = PthLimitedCtrlDict['level'],
                index = index_counter
                )

    return net_out



def transfer_LoadPoint_ppi(
    net, # Pandapipes network
    fluid:str = None,
    tFeed:float = 273.15 + 70, # K
    tReflux:float = 273.15 + 40, # K
    pNetwork:float = 10,
    text_pipes:float = 283.15,
    flh:float = 1700,
    heatingDemandAttr:str = 'demand_use_th',
    thermalPowerAttr = 'Pth_kW',
    thermalPowerAttrConsumers:str = None,
    baseLoadProdNaming:str = 'prod',
    check_for_mass_flow_exceeding:bool = False
    ):

    ### Define residual constant values
    small = 1e-05

    ### Set parameters for network operation
    if fluid is not None:
        ppi.create_fluid_from_lib(
            net = net, 
            name = fluid, 
            overwrite = True)
    else:
        fluid = net.fluid.name

    # Extract fluid properties
    cp_fluid, rho_fluid, nu_fluid, g = extract_FluidProperties_ppi(net, t = (tFeed + tReflux)/2)

    # Discretisation of network pipes
    net.pipe['sections'] = net.pipe.apply(lambda x: 1 if x['length_km'] <= 0.004 else 5, axis = 1)

    # Assign external temeprature boundary condition to pipes
    net.pipe['text_k'] = text_pipes


    ### Set parameters for heat consumers
    if hasattr(net, 'sink'):
        if thermalPowerAttrConsumers is not None:
            net.sink[thermalPowerAttr]      = net.sink[thermalPowerAttrConsumers]
        else:
            net.sink[thermalPowerAttr]          = net.sink[heatingDemandAttr] / flh

        net.sink['mdot_kg_per_s']           = net.sink[thermalPowerAttr] * 1e3 / (cp_fluid * (tFeed - tReflux))
        net.sink['mdot_kg_per_s']           = net.sink['mdot_kg_per_s'].fillna(0)
        net.sink.loc[net.sink['mdot_kg_per_s'] == 0, 'in_service'] = False

        sumSinks_mdot = np.nansum(net.sink['mdot_kg_per_s'])

        print('\n### Total heating demand of consumers connected to network: %.1f MWh ###' % (net.sink[heatingDemandAttr].sum()/(1e3)))
        print('\n### Total thermal power of consumers connected to network at load point: %.1f kW ###' % (net.sink[thermalPowerAttr].sum()))


    if hasattr(net, 'heat_consumer'):
        tempMask = net.heat_consumer[~net.heat_consumer['name'].str.contains(baseLoadProdNaming)].index
        
        if thermalPowerAttrConsumers is not None:
            net.heat_consumer.loc[tempMask, thermalPowerAttr]      = net.heat_consumer.loc[tempMask, thermalPowerAttrConsumers]
        else:
            net.heat_consumer.loc[tempMask, thermalPowerAttr]                     = \
                net.heat_consumer.loc[tempMask, heatingDemandAttr] / flh
        
        net.heat_consumer.loc[tempMask, 'qext_w']                             = \
            net.heat_consumer.loc[tempMask, thermalPowerAttr] * 1e3            

        net.heat_consumer.loc[tempMask, 'controlled_mdot_kg_per_s']           = \
            net.heat_consumer.loc[tempMask, 'qext_w'] / (cp_fluid * (tFeed - tReflux))      

        sumHeatConsumers_mdot = np.nansum(net.heat_consumer.loc[tempMask, 'controlled_mdot_kg_per_s'])

        net.heat_consumer['controlled_mdot_kg_per_s']           = net.heat_consumer['controlled_mdot_kg_per_s'].fillna(0)

        print('\n### Total heating demand of consumers connected to network: %.1f MWh ###' % (net.heat_consumer.loc[tempMask, heatingDemandAttr].sum()/(1e3)))
        print('\n### Total thermal power of consumers connected to network at load point: %.1f kW ###' % (net.heat_consumer.loc[tempMask, thermalPowerAttr].sum()))

    ### Set parameters for producers
    # Peak load producers
    if hasattr(net, 'ext_grid'):
        tempMask = net.ext_grid[net.ext_grid['name'].str.startswith(baseLoadProdNaming)].index
        net.ext_grid.loc[tempMask, 't_k']    = tFeed
        net.ext_grid.loc[tempMask,'p_bar']   = pNetwork


    if hasattr(net, 'circ_pump_pressure'):
        tempMask = net.circ_pump_pressure[net.circ_pump_pressure['name'].str.startswith(baseLoadProdNaming)].index
        net.circ_pump_pressure.loc[tempMask, 't_flow_k'] = tFeed
        net.circ_pump_pressure.loc[tempMask, 'p_flow_bar'] = pNetwork


    ## Base load producers
    if hasattr(net, 'source'):
        tempMask = net.source[net.source['name'].str.contains(baseLoadProdNaming)].index
        
        net.source.loc[tempMask, 'mdot_kg_per_s'] = \
            net.source.loc[tempMask, thermalPowerAttr] * 1e3 / (cp_fluid * (tFeed - tReflux))
        

        if check_for_mass_flow_exceeding:
            if np.nansum(net.source.loc[tempMask, 'mdot_kg_per_s']) >= sumSinks_mdot:
                # Generate relative values in order to scale down mass flows if no exiting flows out of the external grid are desired
                # A residual flow into the network from the external grid is ensured with the constant **small**
                ratiosSource_mdot = net.source.loc[tempMask, 'mdot_kg_per_s'] / np.nansum(net.source.loc[tempMask, 'mdot_kg_per_s'])
                net.source.loc[tempMask, 'mdot_kg_per_s'] = (sumSinks_mdot - small) * ratiosSource_mdot

                print(f'\n### The specified load point resulted in higher cumulated thermal power and mass flows of base load producers in net.source than necessary for all aggregated demands of sinks. The mass flows of source components are adapted so that no reverse flow occurs in the peak load producer ext_grid leaving the network. ###\n')
        
        net.source['mdot_kg_per_s'] = net.source['mdot_kg_per_s'].fillna(0)
        
        net.source.loc[net.source['mdot_kg_per_s'] == 0, 'in_service'] = False

    # Case modelling_type == 'feed_and_reflux'
    if (hasattr(net, 'heat_consumer')):
        if any(net.heat_consumer['name'].str.contains(baseLoadProdNaming)):
            tempMask = net.heat_consumer[net.heat_consumer['name'].str.contains(baseLoadProdNaming)].index

            net.heat_consumer.loc[tempMask, 'qext_w'] = \
                net.heat_consumer.loc[tempMask, thermalPowerAttr] * (-1e3)
            
            net.heat_consumer.loc[tempMask, 'controlled_mdot_kg_per_s'] = \
                abs(net.heat_consumer.loc[tempMask, 'qext_w']) / (cp_fluid * (tFeed - tReflux))
            
            if check_for_mass_flow_exceeding:
                if np.nansum(net.heat_consumer.loc[tempMask, 'controlled_mdot_kg_per_s']) >= sumHeatConsumers_mdot:
                    # Generate relative values in order to scale down mass flows if no reverse flows in circ_pump_pressure component is desired
                    # A residual flow into the feed line of the network from the circ_pump_pressure is ensured with the constant **small**
                    ratiosHeatConsumers_mdot = net.heat_consumer.loc[tempMask, 'controlled_mdot_kg_per_s'] / np.nansum(net.heat_consumer.loc[tempMask, 'controlled_mdot_kg_per_s'])
                    net.heat_consumer.loc[tempMask, 'controlled_mdot_kg_per_s'] = (sumHeatConsumers_mdot - small) * ratiosHeatConsumers_mdot

                    print(f'\n### The specified load point resulted in higher cumulated thermal power and mass flows of base load producers in net.heat_consumer than necessary for all aggregated demands of real heat consumers. The mass flows of base load producer components are adapted so that no reverse flow occurs in the peak load producer circ_pump_pressure as this might lead to convergence problems. ###\n')
                        
            net.heat_consumer['controlled_mdot_kg_per_s']           = net.heat_consumer['controlled_mdot_kg_per_s'].fillna(0)



    # Case modelling_type == 'feed_and_reflux_circPumps'
    if hasattr(net, 'circ_pump_mass'):
        if any(net.circ_pump_mass['name'].str.contains(baseLoadProdNaming)):
            tempMask = net.circ_pump_mass.loc[net.circ_pump_mass['name'].str.contains(baseLoadProdNaming)].index
            
            net.circ_pump_mass.loc[tempMask, 't_flow_k'] = tFeed
            net.circ_pump_mass.loc[tempMask, 'p_flow_bar'] = pNetwork
            net.circ_pump_mass.loc[tempMask, 'mdot_flow_kg_per_s'] = \
                (net.circ_pump_mass.loc[tempMask, thermalPowerAttr] * (1e3)) / (cp_fluid * (tFeed - tReflux))
            
            if check_for_mass_flow_exceeding:
                if np.nansum(net.circ_pump_mass.loc[tempMask, 'mdot_flow_kg_per_s']) >= sumHeatConsumers_mdot:
                    # Generate relative values in order to scale down mass flows if no reverse flows in circ_pump_pressure component is desired
                    # A residual flow into the feed line of the network from the circ_pump_pressure is ensured with the constant **small**
                    ratiosCircPumpMass_mdot = net.circ_pump_mass.loc[tempMask, 'mdot_flow_kg_per_s'] / np.nansum(net.circ_pump_mass.loc[tempMask, 'mdot_flow_kg_per_s'])
                    net.circ_pump_mass.loc[tempMask, 'mdot_flow_kg_per_s'] = (sumHeatConsumers_mdot - small) * ratiosCircPumpMass_mdot
                    print(f'\n### The specified load point resulted in higher cumulated thermal power and mass flows of base loas producers in net.circ_pump_mass than necessary for all aggregated demands of heat consumers. The mass flows of circ_pump_mass components and corresponding flow_controls are adapted so that no reverse flow occurs in the peak load producer circ_pump_pressure as this might lead to convergence problems. ###\n')
            
            net.circ_pump_mass['mdot_flow_kg_per_s'] = net.circ_pump_mass['mdot_flow_kg_per_s'].fillna(0)
            
        if any(net.flow_control['name'].str.contains(baseLoadProdNaming)):
            tempMask = net.flow_control.loc[net.flow_control['name'].str.contains(baseLoadProdNaming)].index

            net.flow_control.loc[tempMask, 'controlled_mdot_kg_per_s'] = \
                (net.flow_control.loc[tempMask, thermalPowerAttr] * (1e3)) / (cp_fluid * (tFeed - tReflux))
            
            if check_for_mass_flow_exceeding:
                if np.nansum(net.flow_control.loc[tempMask, 'controlled_mdot_kg_per_s']) >= sumHeatConsumers_mdot:
                    # Generate relative values in order to scale down mass flows if no reverse flows in circ_pump_pressure component is desired
                    # A residual flow into the feed line of the network from the circ_pump_pressure is ensured with the constant **small**
                    ratiosflowControl_mdot = net.flow_control.loc[tempMask, 'controlled_mdot_kg_per_s'] / np.nansum(net.flow_control.loc[tempMask, 'controlled_mdot_kg_per_s'])
                    net.flow_control.loc[tempMask, 'controlled_mdot_kg_per_s'] = (sumHeatConsumers_mdot - small) * ratiosflowControl_mdot
            
            net.flow_control['controlled_mdot_kg_per_s']           = net.flow_control['controlled_mdot_kg_per_s'].fillna(0)
            

    return net


def extract_FluidProperties_ppi(    
    net, # Pandapipes network
    t:float
    ):

    """
    Function that returns fluid properties at temperature **t** from pandapipes network **net**

    :param net: pandapipes network object
    :param temperature: temperature (K) at which fluid properties shall be extracted
    :returns: tuple of values for (spec. heat capacity, density, kinematic viscosity, gravitaional constant)
    """

    fluid = net.fluid.name

    cp_fluid          = ppi.call_lib(fluid).get_heat_capacity(t)
    rho_fluid         = ppi.call_lib(fluid).get_density(t)
    nu_fluid          = ppi.call_lib(fluid).get_viscosity(t) / rho_fluid
    g                 = ppi.constants.GRAVITATION_CONSTANT

    return(cp_fluid, rho_fluid, nu_fluid, g)


def simultaneity_DH(
        n:int,
        bottomLim:float = 0.7
    )->float:

    """
    Function approximating simultaneous thermal power demands in heating network with *n* connected buildings (DHW and space heating).
    
    Sources:
    Based on equation by Winter et al. (Untersuchungen der Gleichzeitigkeit in kleinen und mittleren Nahwärmenetzen, Winter, W. and Haslauer, T. and Obernberger, I., In: Euroheat & Power 09&10, 2001)
    Also used in Fröhler, R. A. (Nahwärmenetze - Entwicklung eines heuristischen Verfahrens zur Prognose von Nahwärmenetzverlusten mit statistischen Methoden, Diss. TU München, 2022)
    
    :param n: int, number of connected buildings/consumer
    :param bottomLim: float, denotes a bottom limit for the minimum allowed simultaneity factor
    :returns: simultaneity factor of sum of thermal powers in heating network, dependent on the number of connected consumers

    Author: C. Völzel
    """

    a = 0.449677646267461
    b = 0.551234688
    c = 53.8438239
    d = 1.762743268

    GZF = a + b/(1+(n/c)**d)

    if bottomLim is None:
        return GZF
    else:
        return max(bottomLim, GZF)

# Optional:_ Create plot of simultaneity factor
if __name__ == '__main__':

    # Plot appearance
    import matplotlib.pyplot as plt
    xfont = {'fontname': 'Times New Roman', 'fontsize':9}

    # Data
    n = np.linspace(0, 250, 251)
    gzf = np.array(list(map(lambda x: simultaneity_DH(x, bottomLim = None), n)))

    # Create plot
    fig, ax = plt.subplots()
    fig.set_size_inches(3.22, 2, forward=True)
    fig.set_dpi(300)

    ax.plot(n, gzf, color = 'gray')
    
    ax.set_xlabel(r'No. of buildings $n$', **xfont)
    ax.set_xlim([0, 250])
    ax.tick_params(axis='both', which='major', labelsize=8)

    ax.set_ylabel(r'Simultaneity factor $sf$', **xfont)

    ax.grid('major')

    # Create annotation with equation
    ax.annotate(
        r'$sf = 0.4497\cdot\frac{0.5512}{1+(\frac{n}{53.84})^{1.763}}$', 
        xy = (100,0.85),
        xytext=(100, 0.85), 
        fontsize = 9,
        bbox=dict(facecolor='white', alpha=1, edgecolor='black')
        )
    
    # fig.savefig('simultaneityFactor.png', dpi = 'figure', bbox_inches = 'tight')

    

def relThermalLossPower_DH(
        pd:float,
        referToInput:bool = False
    )->float:
    
    """
    Function describing the instantaneous thermal power loss in a DH network (relative to output/input energy) in heating networks.

    :param pd: float denoting the power line density of the network (MW/m).\n
    :param referToInput: boolean selecting if relative factor of tzhermal power loss shall be related to input or output thermal power.\n

    :return: float of relative thermal power loss in DH network
    """

    factor = 0.05 if pd else 0
    pThermalLoss = factor if referToInput else factor / (1-factor)
    
    return pThermalLoss

def relThermalLoss_DH(
        ld:float,
        referToInput:bool = False,
        minVal:float = 0.01
        )->float:
    
    """
    Function describing the partition of thermal losses (relative to output/input thermal energy) in heating networks.
    Data are based on various studies examined in Fröhler, R. A. (Nahwärmenetze - Entwicklung eines heuristischen Verfahrens zur Prognose von Nahwärmenetzverlusten mit statistischen Methoden, Diss. TU München, 2022, p. 94)
    
    :param ld: float, Line density of network (MWh/m)
    :param referToInput: boolean, Boolean if losses shall be provided relative to input thermal energy (True) or relative to demanded thermal energy (False)
    :returns: float of relative value of thermal losses in district heating network
    
    Author: C. Völzel
    """
    
    a = -0.093815
    b = 0.169281

    factor = max(minVal, a * np.log(ld) + b if ld else 0)
    pThermalLoss = factor if referToInput else factor / (1-factor)    
    
    return pThermalLoss

# Optional:_ Create plot of simultaneity factor
if __name__ == '__main__':

    # Plot appearance
    import matplotlib.pyplot as plt
    xfont = {'fontname': 'Times New Roman', 'fontsize':9}

    # Data
    ld = np.linspace(0.1, 5, 50)
    frel = np.array(list(map(lambda x: relThermalLoss_DH(x, referToInput = True), ld)))

    # Create plot
    fig, ax = plt.subplots()
    fig.set_size_inches(3.22, 2, forward=True)
    fig.set_dpi(300)

    ax.plot(ld, frel, color = 'gray')
    
    ax.set_xlabel(r'Line density $ld$ ($MWh/m$)', **xfont)
    ax.set_xlim([0, 5])
    ax.tick_params(axis='both', which='major', labelsize=8)

    ax.set_ylabel(r'Rel. thermal loss $f_{rel,loss}$', **xfont)

    ax.grid('major')

    # Create annotation with equation
    ax.annotate(
        r'$f_{rel,loss} = -0.0938\cdot ln(ld)+0.1693$',
        xy = (0.8, 0.35),
        xytext=(0.8, 0.3), 
        fontsize = 8.5,
        bbox=dict(facecolor='white', alpha=1, edgecolor='black')
        )
    
    fig.savefig('relThermalLoss.png', dpi = 'figure', bbox_inches = 'tight')


def plot_Gs(G:nx.graph, 
            G_list:list, 
            startNodes:list = [],
            plot_start_nodes:bool = False,
            title:str = 'Titel',
            figsize:tuple = (15,8)):
    
    """
    Plotting function to plot generated networkx graphs from list G_copy_list above an (original) graph G.

    :param G: networkx graph object of original graph
    :param startNodes: List of starting nodes to plot
    :param plot_start_nodes: Boolean if start nodes shall be plottes. A list of start nodes has to bve provided!
    """

    # General plot appearance
    plt.figure(figsize=figsize)
    edgeColor_G = 'black'
    edgeColor_Gcopy = 'red'
    nodeColor = 'red'
    nodeSize = 120

    tfont = {'fontsize': 20}  # Font for title

    # Define positions for plotting nodes
    positions = {n: [n[0], n[1]] for n in list(G.nodes)}
    
    # Plotte den Graphen G
    nx.draw(G, positions, node_size=5, edge_color = edgeColor_G)
    
    for G_copy in G_list:        
    
        # Plotte den Graphen G_copy über den Graphen G in Rot
        nx.draw_networkx_edges(G_copy, positions, edge_color=edgeColor_Gcopy)
        nx.draw_networkx_nodes(G_copy, positions, node_size=5, node_color='red')
    
    if plot_start_nodes and len(startNodes)>0:
        # Hervorhebe den random_node in Rot
        nx.draw_networkx_nodes(G_copy, positions, nodelist=startNodes, node_size=nodeSize, node_color=nodeColor)
    elif plot_start_nodes and len(startNodes)==0:
        print('\n### If start nodes shall be plotted, provide proper list of start nodes as iterable of (x,y) pairs. ###')
    else:
        pass
    
    # Stelle die Achsenbeschriftung aus
    plt.axis("off")
    
    # Setze den Titel für den Plot
    plt.title(title, **tfont)
        
    # Zeige den Plot an
    plt.show()


def animate_networkGeneration_gdf(
        gdf:gp.GeoDataFrame,
        legendAttr:dict = {'level':'level'},
        valAttr:str = 'demand_use_th',
        underlyinggdf:gp.GeoDataFrame = None,
        cbarLabel:str = 'heat demand (MWh)',
        legendfont:dict = None,
        cbarfont:dict = None,
        cmap = plt.cm.RdYlGn_r,
        cmapMinMax:tuple = None,
        interval:float = 400,
        savePath:Path = None,
        filename:str = 'video.gif',
        dpi:int = 300,
        figsize:tuple = (1920/300, 1080/300)
    ):


    """
    Function that creates animated, consecutive plots of data in geopandas DataFrames as gif.\n
    Plotting order is defined by the nmumerical index of the DataFrame rows (from first row to last row).\n

    :param gdf: geopandas DataFrame containing data and geometry objects.\n
    :param orderAttr: str denoting atribute name by which to define the plotting order.\n
    :param legendattr: dictionary of strings with keys being the strings shown in the legend and values being the attribute names where values in **gdf** are stored.\n
    :param valAttr: str denoting atribute name by which to define colorbar and to which to assign colormap.\n
    :param underlyinggdf: gdf of data which shall be plotted as underlying data.\n
    :param cbarLabel: str denoting the colorbar label.\n
    :param cmap: matplotlib colormap object. Defaults to plt.cm.RdYlGn_r.\n
    :param cmapMinMax: tuple of floats defining the colorbar range. Refers to **valAttr**. Default to found (min, max) values in **gdf**\n
    :param interval: float defining the duration of a single frame (ms). Defaults to 400 ms.\n
    :param savePath: Path object defining the absolute path where to save the gif-animation.

    :returns: No returns
    """

    # Initializations
    if cmapMinMax is None:
        cmapMinMax = (np.nanmin(gdf[valAttr]), np.nanmax(gdf[valAttr]))
    
    # Initialize plotting figure
    fig, ax = plt.subplots(figsize = figsize, dpi = dpi)
    ax.set_axis_off()

    if legendfont is None:
        legendfont = {'family': 'Times New Roman', 'size': 10}

    if cbarfont is None:
        cbarfont = {'labelfontfamily': 'Times New Roman', 'labelsize': 10}

    # Create container for single frames
    frames = []
    frameArtists = []
    edgeCollections = []
    nodeCollections = []

    # Plot original (complete) graph
    if underlyinggdf is not None:
        edgeArtist_orig = underlyinggdf.plot(
                ax = ax,
                column = valAttr, 
                vmin = cmapMinMax[0], 
                vmax = cmapMinMax[1],
                linewidth = 1,
                alpha = 0.5,
                cmap = cmap,
                zorder = 0
                )
    
    # Plot gdf data in order of the specified attribute **orderAttr**
    for nn, row in gdf.iterrows():   
       
        # Plot graph up to current edge
        edgeArtist = gdf.loc[:nn].plot(
            ax = ax,
            column = valAttr, 
            vmin = cmapMinMax[0], 
            vmax = cmapMinMax[1],
            linewidth = 2.5,
            alpha = 1,
            cmap = cmap,
            zorder = 10
            )
        
        # Create colorbar
        if nn == 1:                
            cbar = fig.colorbar(edgeArtist.collections[0], ax=ax, location='right')
            cbar.set_label(cbarLabel, **legendfont)
            cbar.ax.tick_params(**cbarfont)

        werte = [row[v] if type(row[v]) == float else str(row[v]) for k,v in legendAttr.items()]
        max_len = max([len(f"{w:.3f}" if isinstance(w, float) else str(w)) for w in werte])
        
        legString = []
        for k,v in legendAttr.items():
            legString.append(f'{k} {row[v]:>{max_len}.3}') if type(row[v]) == float else legString.append(f'{k} {str(row[v]):>{max_len}}')

        legString_final = '\n'.join(legString)

        an1 = ax.annotate(
            legString_final,
            xy=(0.01, 0.99),
            xycoords="figure fraction",
            va="top", 
            ha="left",
            bbox=dict(boxstyle="round", fc="w"),
            **legendfont
            )
        
        # Gather all artiosts within axes
        line_collections = [artist for artist in ax.get_children() if isinstance(artist, LineCollection)]
        path_collections = [artist for artist in ax.get_children() if isinstance(artist, PathCollection)]                

        edgeCollections += line_collections
        nodeCollections += path_collections

        frameArtists = edgeCollections.copy() + nodeCollections.copy() + [an1]
        frames.append(frameArtists)

    # Saving animation
    if savePath is not None:
    
        if not os.path.exists(savePath):
            os.makedirs(savePath)
            print(f'### New directory for saving visualization videos is created in {savePath} ###')

        ani = matplotlib.animation.ArtistAnimation(fig=fig, artists=frames, interval=interval)
        ani.save(filename=savePath / Path(filename), writer="pillow", dpi = dpi)


    return