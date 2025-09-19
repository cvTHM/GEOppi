# -*- coding: utf-8 -*-
import numpy as np

from pandapower.control.basic_controller import (BasicCtrl)
from pandapower.control.util.auxiliary import (drop_same_type_existing_controllers, log_same_type_existing_controllers)
from pandapower.auxiliary import (get_free_id)

### Definition of plausibility checks for user-defined pf options ###
def checkUserPFOptionsThermal(net, ctrlName:str):
        
    """
    Function that checks for calculation modes entered in the pandapipes network **net**. If a defined controller requires a calculation mode with heat transfer calculation to work properly, this function can print a warning message.\n

    :param net: pandapipes network object
    :param ctrlName: str denoting the name of the controller class
    :returns: No return
    """
        
    # Plausibility check for user-defined pipeflow options
    if hasattr(net, 'user_pf_options'):
        pf_options = net.user_pf_options
        if 'mode' in pf_options.keys():
            if pf_options['mode'] not in ('sequential', 'bidirectional', 'all', 'heat'): # Check if heat transfer calculation is included in user-defined pipeflow options
                modee = pf_options['mode']
                print(f'\n### Attention: Defined mode in user_pf_options of the network {modee} does not include heat transfer calculations. Make sure to set mode to any of ("sequential", "bidirectional", "all", "heat") as this is necessary for the defined controller type {ctrlName}. ###\n')
        else:
            print(f'\n### Attention: No user-defined option "mode" in user_pf_options found. Make sure to include mode from any of ("sequential", "bidirectional", "all", "heat") in user_pf_options as this is necessary for the defined controller type {ctrlName}. ###\n')
    else:
        print(f'\n### Attention: No user-defined options in user_pf_options found. Make sure to include mode from any of ("sequential", "bidirectional", "all", "heat") in user_pf_options as this is necessary for the defined controller type {ctrlName}. ###\n')


#####################################################################


### Definition of custom controllers ###

class ppi_HeatConsumerCOPConversionQdem(BasicCtrl):

    def __init__(self, net, Tsink_heatingSystem:list = [273.15 + 60], efficiency:list = [0.5], deltaTsource:list = [4], heatConsumer_idxs:list = None, abs_tol:float = 0.1, proportional_gain:float = 0.5, order:int = 1, level:int = 1, index:int = None, **kwargs):

        super(ppi_HeatConsumerCOPConversionQdem, self).__init__(net, **kwargs)

        self.Tsink_heatingSystem = Tsink_heatingSystem
        self.efficiency = efficiency
        self.initialCOP = 4
        self.deltaTsource = deltaTsource
        self.heatConsumer_idxs = heatConsumer_idxs if heatConsumer_idxs is not None else net.heat_consumer.index
        self.abs_tol = abs_tol
        self.proportional_gain = proportional_gain
        self.iterations = 0
        self.convergence = False

        if index is None and "controller" in net.keys():
            index = get_free_id(net.controller)

        self.index = self.add_controller_to_net(net = net, in_service = True, initial_run = True, index = index, order = order, level = level, recycle = False, overwrite = True, drop_same_existing_ctrl = True, **kwargs)

        # Assign Tsink_heatingSystem and efficiencies to heat consumers individually (if desired)
        # Default values for missing entries:
        ## Tsink_heatingSystem = 273.15+60 K
        ## efficiency = 0.5
        ## deltaTsource = 4 K

        if isinstance(Tsink_heatingSystem, list):
            if len(Tsink_heatingSystem) > len(self.heatConsumer_idxs):
                self.Tsink_heatingSystem = np.array(Tsink_heatingSystem[:len(self.heatConsumer_idxs)])
            elif len(Tsink_heatingSystem) < len(self.heatConsumer_idxs):
                self.Tsink_heatingSystem = np.array(Tsink_heatingSystem.extend(list(np.ones(len(self.heatConsumer_idxs)-len(Tsink_heatingSystem))*(273.15+60))))
        else:
            self.Tsink_heatingSystem = np.array([Tsink_heatingSystem] * len(self.heatConsumer_idxs))

        if isinstance(efficiency, list):
            if len(efficiency) > len(self.heatConsumer_idxs):
                self.efficiency = np.array(efficiency[:len(self.heatConsumer_idxs)])
            elif len(efficiency) < len(self.heatConsumer_idxs):
                self.efficiency = np.array(efficiency.extend(list(np.ones(len(self.heatConsumer_idxs)-len(efficiency))*(0.5))))
        else:
            self.efficiency = np.array([efficiency] * len(self.heatConsumer_idxs))   

        if isinstance(deltaTsource, list):
            if len(deltaTsource) > len(self.heatConsumer_idxs):
                self.deltaTsource = np.array(deltaTsource[:len(self.heatConsumer_idxs)])
            elif len(deltaTsource) < len(self.heatConsumer_idxs):
                self.deltaTsource = np.array(deltaTsource.extend(list(np.ones(len(self.heatConsumer_idxs)-len(deltaTsource))*(4))))
        else:
            self.deltaTsource = np.array([deltaTsource] * len(self.heatConsumer_idxs))

        # Write 

        # Save initial conditions of network
        ## During control step, ccontrolled_mdoit_kg_per_s and qext_w in heat_consumer components are changed
        net.heat_consumer['controlled_mdot_kg_per_s_init'] = net.heat_consumer['controlled_mdot_kg_per_s'].copy() 
        net.heat_consumer['qext_w_init'] = net.heat_consumer['qext_w'].copy()       


    def add_controller_to_net(self, net, in_service, initial_run, order, level, index, recycle,
                              drop_same_existing_ctrl, overwrite, **kwargs):
        """
        adds the controller to net['controller'] dataframe.

        INPUT:
            **in_service** (bool) - in service status

            **order** (int) - order

            **index** (int) - index

            **recycle** (bool) - if controller needs a new bbm (ppc, Ybus...) or if it can be used \
                                 with prestored values. This is mostly needed for time series \
                                 calculations

        """
        if drop_same_existing_ctrl:
            drop_same_type_existing_controllers(net, type(self), index=index, **kwargs)
        else:
            log_same_type_existing_controllers(net, type(self), index=index, **kwargs)

        # use base class method to raise an error if the object is in DF and overwrite = False
        # if the index is None, the base class is in charge of obtaining the next free index in the data frame
        fill_dict = {"in_service": in_service, "initial_run": initial_run, "recycle": recycle,
                     "order": order, "level": level}
        added_index = super().add_to_net(net=net, element='controller', index=index, overwrite=overwrite,
                           fill_dict=fill_dict, preserve_dtypes=True)
        return added_index

    def initialize_control(self, net):
        # Extract indices of heat consumers which are active or inactive -> Definition is based on specified thermal demand
        self.heatConsumer_active_idxs = list(net.heat_consumer.loc[(net.heat_consumer['qext_w'] > 0) & (~net.heat_consumer['controlled_mdot_kg_per_s'].isna())].index)
        self.heatConsumer_inactive_idxs = list(set(net.heat_consumer.index).difference(set(self.heatConsumer_active_idxs)))

        # Extract total demanded heating power (sink side) at active consumers
        ## Assumption: Current thermal demand load on sink side is always written in attribute qext_w for each time step
        self.heatConsumer_active_Qdem = net.heat_consumer.loc[self.heatConsumer_active_idxs, 'qext_w'].values
       
        # Initialization of target source side thermal power extraction
        self.heatConsumer_active_Qdem_source_target = self.heatConsumer_active_Qdem / self.initialCOP * (self.initialCOP-1)

    
    def get_heatConsumer_states(self, net, idxs:list):

        T_from = net.res_heat_consumer.loc[idxs, 't_from_k'].values
        qext_source = net.res_heat_consumer.loc[idxs, 'qext_w'].values

        return T_from, qext_source
    

    def control_step(self, net):

        self.iterations += 1

        ### Control step is performed for all currently active heat consumers

        # Get current states at active heat consumers
        T_from_current, qext_source_current = self.get_heatConsumer_states(net = net, idxs = self.heatConsumer_active_idxs)

        # Update current heat capacity of fluid
        cp = net.fluid.get_heat_capacity(273.15+10)

        # Determine current set source-side extracted power, dependent on current temperature T_from
        COP_set = (self.Tsink_heatingSystem[self.heatConsumer_active_idxs]) / (self.Tsink_heatingSystem[self.heatConsumer_active_idxs] - T_from_current) * self.efficiency[self.heatConsumer_active_idxs]
        qext_source_set = self.heatConsumer_active_Qdem / COP_set * (COP_set-1)
        
        qext_error = qext_source_set - qext_source_current        

        # Set new values for control step
        self.heatConsumer_active_Qdem_source_target = qext_source_set

        qext_source_new = qext_source_current + qext_error * self.proportional_gain
        mdot_source_new = qext_source_new / (cp * self.deltaTsource[self.heatConsumer_active_idxs])

        # Transfer new set values to heat consumer models
        net.heat_consumer.loc[self.heatConsumer_active_idxs, 'qext_w'] = qext_source_new
        net.heat_consumer.loc[self.heatConsumer_active_idxs, 'controlled_mdot_kg_per_s'] = mdot_source_new

        return super(ppi_HeatConsumerCOPConversionQdem, self).control_step(net)
    

    def is_converged(self, net):

        convergence = False

        # Extract results at active consumers
        _, qext_source_current = self.get_heatConsumer_states(net = net, idxs = self.heatConsumer_active_idxs)

        qext_source_error = qext_source_current - self.heatConsumer_active_Qdem_source_target

        if (all(abs(qext_source_error) <= self.abs_tol)):
            convergence = True
            self.convergence = True

        return convergence
    
    def finalize_control(self, net):
        if self.convergence:
            T_from_current, qext_source_current = self.get_heatConsumer_states(net = net, idxs = self.heatConsumer_active_idxs)

            # Temporary outputs for test purposes
            net.res_heat_consumer.loc[self.heatConsumer_active_idxs, 'efficiency'] = self.efficiency[self.heatConsumer_active_idxs]
            net.res_heat_consumer.loc[self.heatConsumer_active_idxs, 'COP'] = self.Tsink_heatingSystem[self.heatConsumer_active_idxs] / (self.Tsink_heatingSystem[self.heatConsumer_active_idxs] - T_from_current) * self.efficiency[self.heatConsumer_active_idxs]




class ppi_HeatConsumerSetTempCtrl(BasicCtrl):
    """
    Class of controller in pandapipes network for...
    - controlling fixed temperature values at junctions of components "heat_consumers"\n
    - varying the mass flow at heat consumer components until the target temperatures are reached.\n
    - A minimum value for the controlled mass flow can be specified with min_mdot.\n
        - Convergence is reached either if all heat consumers show a deviation of the set target temperature at to_junction or from_junction < tolerance OR if heat consumers do not meet the set target temperature due to the minimum mass flow specified (e.g. in low load points).\n

    The junction at which to control the temperature can be specified, so it can be either the "from_junction" or the "to_junction".\n
    Indexes of heat consumer components at which to control the temperature can be specified as "heatConsumer_idxs". If None, the controller will be applied to all heat consumer components.\n
    Accessing results from heat transfer calculations in pandapipes network.\n
    Inherits from class **BasicCtrl** in pandapower.control.basic_controller.\n

    Prerequisits for the use of this controller:
    - heat consumers feature positive mass flow from feed line (from_junction) to return line (to_junction)
    - heat consumers feature positive values of **qext_w**, resulting in heat extraction from the network (heating application)
    """

    def __init__(self, net, target_T:float = 300,  min_dT:float = None, min_mdot:float = None, junction_type:str = 'to_junction', proportional_gain:float = 0.5, abs_tol:float = 0.5, heatConsumer_idxs:list = None, index:int = None, order:int = 1, level:int = 1, **kwargs):
        """
        Init function of the controller **ppi_HeatConsumerSetTempCtrl**

        :param net: pandapipes network object.\n
        :param target_T: float denoting the target value of the controlled temperature at junctions **junction_type**, defaults to 300\n
        :param min_dT: float denoting minimum temperature difference between **target_T** and corresponding complementary temperature at the other connected node of the heat consumer. If None, this control is not applied.\n
        :param min_mdot: float denoting the minimum mass flow to maintain for each heat consumer. If None, this defaults to 0 to ensure positive mass flow from **from_junction** to **to_junction**.\n
        :param junction_type: str denoting which junction temperature shall be controlled (either "from_junction" or "to_junction"); default "to_junction"\n
        :param proportional_gain: float denoting a proportional gain for the control step, defaults to 0.5\n
        :param abs_tol: float denoting the absolute allowed tolerance at each heat consumer, defaults to 0.5\n
        :param heatConsumer_idxs: list-like, indexes of the heat consumers at which to adapt the mass flow, defaults to None\n
        :param index: int defaults to None\n
        :param order: int defaults to 1\n
        :param level: int defaults to 1\n
        """

        super(ppi_HeatConsumerSetTempCtrl, self).__init__(net, **kwargs)

        self.target_T = target_T
        self.min_dT = min_dT if min_dT is not None else 5
        self.min_mdot = min_mdot if min_mdot is not None else 0
        self.heatConsumer_idxs = heatConsumer_idxs if heatConsumer_idxs is not None else net.heat_consumer.index
        self.junction_type = junction_type
        self.proportional_gain = proportional_gain
        self.iterations = 0
        self.abs_tol = abs_tol

        if index is None and "controller" in net.keys():
            index = get_free_id(net.controller)

        self.index = self.add_controller_to_net(net = net, in_service = True, initial_run = True, index = index, order = order, level = level, recycle = False, overwrite = True, drop_same_existing_ctrl = True, **kwargs)

        # Plausibility check for user-defined pipeflow options
        checkUserPFOptionsThermal(net = net, ctrlName = "ppi_HeatConsumerSetTempCtrl")

        # Save initial conditions of network
        ## During control step, ccontrolled_mdoit_kg_per_s in heat_cosnumer components is changed
        net.heat_consumer['controlled_mdot_kg_per_s_init'] = net.heat_consumer['controlled_mdot_kg_per_s'].copy()


    def add_controller_to_net(self, net, in_service, initial_run, order, level, index, recycle,
                              drop_same_existing_ctrl, overwrite, **kwargs):
        """
        adds the controller to net['controller'] dataframe.

        INPUT:
            **in_service** (bool) - in service status

            **order** (int) - order

            **index** (int) - index

            **recycle** (bool) - if controller needs a new bbm (ppc, Ybus...) or if it can be used \
                                 with prestored values. This is mostly needed for time series \
                                 calculations

        """
        if drop_same_existing_ctrl:
            drop_same_type_existing_controllers(net, type(self), index=index, **kwargs)
        else:
            log_same_type_existing_controllers(net, type(self), index=index, **kwargs)

        # use base class method to raise an error if the object is in DF and overwrite = False
        # if the index is None, the base class is in charge of obtaining the next free index in the data frame
        fill_dict = {"in_service": in_service, "initial_run": initial_run, "recycle": recycle,
                     "order": order, "level": level}
        added_index = super().add_to_net(net=net, element='controller', index=index, overwrite=overwrite,
                           fill_dict=fill_dict, preserve_dtypes=True)
        return added_index
    
    def initialize_control(self, net):
        # Extract indices of heat consumers which are active or inactive -> Definition is based on specified thermal demand
        self.heatConsumer_active_idxs = list(net.heat_consumer.loc[(net.heat_consumer['qext_w'] > 0) & (~net.heat_consumer['controlled_mdot_kg_per_s'].isna())].index)
        self.heatConsumer_inactive_idxs = list(set(net.heat_consumer.index).difference(set(self.heatConsumer_active_idxs)))
        
        self.cp = net.fluid.get_heat_capacity(self.target_T)
        self.junction_idxs = net.heat_consumer['from_junction'].values if self.junction_type == 'from_junction' else net.heat_consumer['to_junction'].values

        # Initialize target temperature array for each heat consumer (necessary if individual target has to be modified, e.g. to meet min. temperature difference instead of set temperature)      
        self.target_T_arr = self.target_T * np.ones(len(self.heatConsumer_active_idxs))

        # Initialize target_T_arr in heat_consumer DF
        net.heat_consumer['target_T'] = np.nan

    def is_converged(self, net):

        convergence = False

        # Extract current results
        _, mdot_from_current, T_from_current, T_to_current = self.get_heatConsumer_states(net = net, idxs = self.heatConsumer_active_idxs)
        
        # Control temperature at "from_junction"
        if self.junction_type == 'from_junction':
            T_error = self.target_T_arr - T_from_current

        # Control temperature at "to_junction"
        else:
            T_error = self.target_T_arr - T_to_current

        if (all(abs(T_error) <= self.abs_tol)) | (all(mdot_from_current[np.where(abs(T_error) > self.abs_tol)[0]] == self.min_mdot)):
            convergence = True

        return convergence
    
    def get_heatConsumer_states(self, net, idxs:list):        
        qext_w = net.res_heat_consumer.loc[idxs, 'qext_w'].values
        mdot = net.res_heat_consumer.loc[idxs, 'mdot_from_kg_per_s'].values
        T_from = net.res_heat_consumer.loc[idxs, 't_from_k'].values
        T_to = net.res_heat_consumer.loc[idxs, 't_to_k'].values

        return qext_w, mdot, T_from, T_to

    def control_step(self, net):

        self.iterations += 1

        ### Control step is performed for all currently active heat consumers

        # Get current states at heat consumers
        qext_current, mdot_current, T_from_current, T_to_current = self.get_heatConsumer_states(net = net, idxs = self.heatConsumer_active_idxs)

        # Control temperature at "from_junction"
        if self.junction_type == 'from_junction':            
            if (self.min_dT is not None) | (self.target_T_arr - T_to_current < 0):
                idxs_temp = np.where(self.target_T_arr - T_to_current < self.min_dT)[0]
                self.target_T_arr[idxs_temp] = T_from_current[idxs_temp] + self.min_dT

            mdot_error = (qext_current/(self.cp * (self.target_T_arr - T_to_current)) - qext_current/(self.cp*(T_from_current - T_to_current)))

        # Control temperature at "to_junction"
        else:
            if any(T_from_current - self.target_T_arr < 0):

                # Get indices of heat consumers who feature supply temperature lower than target reflux temperature (indicates too low mass flow)
                idxs_undercut = np.where(T_from_current - self.target_T_arr < 0)[0]
                print(f'\n### Supply temperature is lower than target temperature at heat_consumer idxs {np.array(self.heatConsumer_active_idxs)[idxs_undercut]} ###')
                
                # Get indices of heat consumers whose supply temperature is less than self.min_dT above target reflux temperature
                idxs_temp = np.where(T_from_current - self.target_T_arr < self.min_dT)[0]

                self.target_T_arr[idxs_temp] = T_from_current[idxs_temp] - self.min_dT

            mdot_error = (qext_current/(self.cp * (T_from_current - self.target_T_arr)) - qext_current/(self.cp*(T_from_current - T_to_current)))

        new_mdot = mdot_current + mdot_error * self.proportional_gain
        new_mdot = new_mdot.clip(self.min_mdot, np.inf)

        # Set new controlled mas flow at active heat consumers
        net.heat_consumer.loc[self.heatConsumer_active_idxs, 'controlled_mdot_kg_per_s'] = new_mdot

        return super(ppi_HeatConsumerSetTempCtrl, self).control_step(net)
    
    def finalize_control(self, net):
        if self.converged:
            net.heat_consumer.loc[self.heatConsumer_active_idxs, 'target_T'] = self.target_T_arr
    

class ppi_CircPumpMassPthermalCtrl(BasicCtrl):
    """
    Class of controller in pandapipes network for...
    - controlling mass flow through component circ_pump_const_mass_flow and corresponding flow_controller in such way that the defined target value of thermal power is reached.\n

    Accessing results from heat transfer calculations in pandapipes network.\n
    The controller may be used in closed loop networks (feed and reflux line) with a modelling approach using a combination of circ_pump_const_mass_flow and downstream flow_controller.\n

    Inherits from class **BasicCtrl** in pandapower.control.basic_controller    
    """

    def __init__(self, net, target_Pth, circPumpMass_idxs:list = None, flow_controller_idxs:list = None, abs_tol:float = 0.1, proportional_gain:float = 0.5, order:int = 0, level:int = 0, index:int = None, **kwargs):
        """
        Init function for the controller **ppi_CircPumpMassPthermalController**

        :param net: pandapipes network.\n
        :param target_Pth: float denoting the target value of thermal power (W) specified for the chosen circ_pump_const_mass_flow.\n
        :param circPumpMass_idxs: int denoting the indexes of the circ_pump_const_mass_flow components for which to create the controller.\n
        :param flow_controller_idxs: int denoting the index of the flow_controller components corresponding to the circ_pump_const_mass_flow components.\n
        :param abs_tol: float denoting the absolute tolerance between comutet/controlled and target value of thermal power, defaults to 0.1.\n
        :param proportional_gain: float, denoting a proportional gain for the controller, defaults to 0.5\n
        :param order: defaults to 0\n
        :param level: defaults to 0\n
        :param index: defaults to None\n
        """

        super(ppi_CircPumpMassPthermalCtrl, self).__init__(net, **kwargs)

        self.circPumpMass_idxs = circPumpMass_idxs if circPumpMass_idxs is not None else net.circ_pump_mass.index
        self.flow_controller_idxs = flow_controller_idxs if flow_controller_idxs is not None else net.flow_control.index
        self.target_Pth = target_Pth if isinstance(target_Pth, np.ndarray) else np.array(target_Pth) if isinstance(target_Pth, list) else target_Pth
        self.abs_tol = abs_tol
        self.proportional_gain = proportional_gain
        self.iterations = 0

        if index is None and "controller" in net.keys():
            index = get_free_id(net.controller)

        # Plausibility check for user-defined pipeflow options
        checkUserPFOptionsThermal(net = net, ctrlName = "ppi_CircPumpMassPthermalController")

        # Controller is created and added to network **net**
        self.index = self.add_controller_to_net(net = net, in_service = True, initial_run = True, index = index, order = order, level = level, recycle = False, overwrite = True, drop_same_existing_ctrl = True, **kwargs)

        # Save initial conditions of network
        ## During control step, circ_pump_mass mass flows and flow_control mass flows are varied
        net.circ_pump_mass['mdot_flow_kg_per_s_init'] = net.circ_pump_mass['mdot_flow_kg_per_s'].copy()
        net.flow_control['controlled_mdot_kg_per_s_init'] = net.flow_control['controlled_mdot_kg_per_s'].copy()

    
    
    def add_controller_to_net(self, net, in_service, initial_run, order, level, index, recycle,
                              drop_same_existing_ctrl, overwrite, **kwargs):
        """
        adds the controller to net['controller'] dataframe.

        INPUT:
            **in_service** (bool) - in service status

            **order** (int) - order

            **index** (int) - index

            **recycle** (bool) - if controller needs a new bbm (ppc, Ybus...) or if it can be used \
                                 with prestored values. This is mostly needed for time series \
                                 calculations

        """
        if drop_same_existing_ctrl:
            drop_same_type_existing_controllers(net, type(self), index=index, **kwargs)
        else:
            log_same_type_existing_controllers(net, type(self), index=index, **kwargs)

        # use base class method to raise an error if the object is in DF and overwrite = False
        # if the index is None, the base class is in charge of obtaining the next free index in the data frame
        fill_dict = {"in_service": in_service, "initial_run": initial_run, "recycle": recycle,
                     "order": order, "level": level}
        added_index = super().add_to_net(net=net, element='controller', index=index, overwrite=overwrite,
                           fill_dict=fill_dict, preserve_dtypes=True)
        return added_index
    

    def initialize_control(self, net):

        # Define residual mass flow in circ_pump_mass components
        self.mdotmin_circPumpMass = 1e-05


    def control_step(self, net):

        self.iterations += 1

        current_Pth, current_mdot, current_tempdiff, cp = self.get_Pth_mdot_tempdiff(net, idxs = self.circPumpMass_idxs)
        Pth_error = self.target_Pth - current_Pth

        new_mdot = np.fmax(self.mdotmin_circPumpMass, current_mdot + Pth_error / (current_tempdiff * cp) * self.proportional_gain)
        net.circ_pump_mass.loc[self.circPumpMass_idxs, 'mdot_flow_kg_per_s'] = new_mdot
        net.flow_control.loc[self.flow_controller_idxs, 'controlled_mdot_kg_per_s'] = new_mdot

        return super(ppi_CircPumpMassPthermalCtrl, self).control_step(net)

    def is_converged(self, net):

        # Extract current results
        current_Pth, _, _, _ = self.get_Pth_mdot_tempdiff(net, idxs = self.circPumpMass_idxs)
        Pth_error = self.target_Pth - current_Pth

        if all(abs(Pth_error) <= self.abs_tol):
            convergence = True
        else:
            convergence = False

        return convergence

    def get_Pth_mdot_tempdiff(self, net, idxs):
        """
        Function to extract necessary results from result DataFrame of pandapipes network component.\n

        :return: tuple of four arrays denoting (current thermal power, current mass flow (from), current temperature difference (to - from), value of specific heat capacity) at specified index positions of circ_pump_mass_components
        """

        mdot_from = net.res_circ_pump_mass.loc[idxs, 'mdot_from_kg_per_s'].values
        temp_from_junction = net.res_circ_pump_mass.loc[idxs, 't_from_k'].values
        temp_to_junction = net.res_circ_pump_mass.loc[idxs, 't_to_k'].values
        tempdiff = temp_to_junction - temp_from_junction
        cp = net.fluid.get_heat_capacity((temp_to_junction+temp_from_junction)/2)

        Pth = mdot_from * cp * tempdiff

        return Pth, mdot_from, tempdiff, cp

class ppi_CircPumpMassPthermalCtrl_limited(BasicCtrl):
    """
    Class of controller in pandapipes network for...
    - controlling mass flow through component circ_pump_const_mass_flow and corresponding flow_controller in such way that the defined target value of thermal power is reached.\n
    - Ensures that the peak load producer, modelled as circ_pump_const_pressure component, features a residual mass flow if its power is not needed.\n
    - Makes use of a cascaded priorisation order if multiple base load producers, modelled as circ_pump_const_mass_flow components, are included in the network. Components with lower priority values are first used to cover the heat demand in the network.\n
    
    **Controllers for base load producer components must be set with ascending level corresponding to their priority values!**\n
    
    Controlled variables:
    - mass flow of base load producer (circ_pump_mass and corresponding flow_control component)\n

    Target variables/Restrictions:
    - target values for thermal power of base load producer
    - restricts mass flow of peak load producer (circ_pump_pressure component) to minimum value > 0\n

    Accessing results from heat transfer calculations in pandapipes network.\n
    The controller may be used in closed loop networks (feed and reflux line) with a modelling approach using a combination of circ_pump_const_mass_flow and downstream flow_controller.\n

    Inherits from class **BasicCtrl** in pandapower.control.basic_controller    
    """

    def __init__(self, net, target_Pth, circPumpMass_idx:int, flow_controller_idx:int, circPumpPressure_idx:int = None, priority_list = None, priority:int = None, abs_tol:float = 0.1, proportional_gain:float = 0.5, order:int = 0, level:int = 0, index:int = None, **kwargs):

        """
        Init function for the controller **ppi_CircPumpMassPthermalCtrl_limited**

        :param net: pandapipes network.\n
        :param target_Pth: float denoting the target value of thermal power (W) specified for the chosen circ_pump_const_mass_flow.\n
        :param circPumpMass_idx: int denoting the index of the circ_pump_const_mass_flow component for which to create the controller.\n
        :param flow_controller_idx: int denoting the index of the flow_controller component corresponding to the circ_pump_const_mass_flow component.\n
        :param circPumpPressure_idx: int denoting the index of the (single) circ_pump_const_pressure component which serves as a peak load producer in the network.\n
        :param priority_list: list of int, denoting the priorisation order/cascade for circ_pump_const_mass_flow components if multiple components are in the network. Lower values mean higher priority and priorised operation as base load producer.\n
        :param priority: int, priority value for the current component.\n
        :param abs_tol: float denoting the absolute tolerance between computed/controlled and target value of thermal power, defaults to 0.1.\n
        :param proportional_gain: float, denoting a proportional gain for the controller, defaults to 0.5\n
        :param order: defaults to 0\n
        :param level: defaults to 0\n
        :param index: defaults to None\n
        """

        super(ppi_CircPumpMassPthermalCtrl_limited, self).__init__(net, **kwargs)

        self.circPumpMass_idx = circPumpMass_idx
        self.flow_controller_idx = flow_controller_idx
        self.target_Pth = target_Pth
        self.circPumpPressure_idx = circPumpPressure_idx
        self.priority_list = priority_list
        self.priority = priority
        self.abs_tol = abs_tol
        self.proportional_gain = proportional_gain
        self.iterations = 0
        # self.controllerType = 'ppi_CircPumpMassPthermalCtrl_limited'

        if index is None and "controller" in net.keys():
            index = get_free_id(net.controller)

        # Plausibility check for user-defined pipeflow options
        checkUserPFOptionsThermal(net = net, ctrlName = "ppi_CircPumpMassPthermalController")

        # Controller is created and added to network **net**
        self.index = self.add_controller_to_net(net = net, in_service = True, initial_run = True, index = index, order = order, level = level, recycle = False, overwrite = True, drop_same_existing_ctrl = True, **kwargs)

        # Save initial conditions of network
        ## During control step, circ_pump_mass mass flows and flow_control mass flows are varied
        net.circ_pump_mass['mdot_flow_kg_per_s_init'] = net.circ_pump_mass['mdot_flow_kg_per_s'].copy()
        net.flow_control['controlled_mdot_kg_per_s_init'] = net.flow_control['controlled_mdot_kg_per_s'].copy()

    
    
    def add_controller_to_net(self, net, in_service, initial_run, order, level, index, recycle,
                              drop_same_existing_ctrl, overwrite, **kwargs):
        """
        adds the controller to net['controller'] dataframe.

        INPUT:
            **in_service** (bool) - in service status

            **order** (int) - order

            **index** (int) - index

            **recycle** (bool) - if controller needs a new bbm (ppc, Ybus...) or if it can be used \
                                 with prestored values. This is mostly needed for time series \
                                 calculations

        """
        if drop_same_existing_ctrl:
            drop_same_type_existing_controllers(net, type(self), index=index, **kwargs)
        else:
            log_same_type_existing_controllers(net, type(self), index=index, **kwargs)

        # use base class method to raise an error if the object is in DF and overwrite = False
        # if the index is None, the base class is in charge of obtaining the next free index in the data frame
        fill_dict = {"in_service": in_service, "initial_run": initial_run, "recycle": recycle,
                     "order": order, "level": level}
        added_index = super().add_to_net(net=net, element='controller', index=index, overwrite=overwrite,
                           fill_dict=fill_dict, preserve_dtypes=True)
        return added_index
    

    def initialize_control(self, net):

        self.mdotmin_circPumpPressure = 1e-05 # min value for mass flow rate of circ_pump_pressure component (peak load producer)
        self.mdotmin_circPumpMass = (1e-05)/(len(net.circ_pump_mass) + 1) # min value for residual mass flow of circ_pump_mass component (base load producer)

        self.target_Pth_limited = self.target_Pth

        # Initialize convergence flag
        self.convergence = False

        ## Applying cascaded operation of base load producers
        # If the priority value of the current circ_pump_mass is the lowest found in priority list, it is the first to be used for demand coverage.
        # In this case, mass flows of all circ_pump_mass components are set to residual value to support convergence and avoid reverse flow in circ_pump_pressure component.

        # Get indices of all circ pump mass components and flow control components which feature the identical controller type handling
        self.ctrlType = type(net.controller.loc[self.index, 'object'])

        (self.cpm_idxs, self.fc_idxs) = np.array([(k.circPumpMass_idx, k.flow_controller_idx) for k in net.controller['object'] if isinstance(k, self.ctrlType)]).T

        # Initialize priority list if necessary
        if self.priority_list is None:
            self.priority_list = [0] * len(self.cpm_idxs)
        if self.priority is None:    
            self.priority = 0

        # Get indices of circ pump mass components with lower priority_list values
        self.higherPriority_idxs = list(np.where(np.array(self.priority_list) < self.priority)[0])

        if self.priority == self.priority_list.index(min(self.priority_list)): # -> Highest priority

            self.highest_priority = True # True if highest priorioty or same priority in all controlled components

            # If current circ pump mass features highest priority, initialize thermal power of all circ pump mass and flow control components accessed by the same type of controller in the network
            net.circ_pump_mass.loc[self.cpm_idxs, 'mdot_flow_kg_per_s'] = self.mdotmin_circPumpMass
            net.flow_control.loc[self.fc_idxs, 'controlled_mdot_kg_per_s'] = self.mdotmin_circPumpMass

            if hasattr(net, 'res_circ_pump_mass'):
                net.res_circ_pump_mass.loc[self.cpm_idxs, 'mdot_from_kg_per_s'] = self.mdotmin_circPumpMass

            if hasattr(net, 'res_flow_control'):
                net.res_flow_control.loc[self.fc_idxs, 'controlled_mdot_kg_per_s'] = self.mdotmin_circPumpMass

        if 'priorityCtrl' not in net.circ_pump_mass.columns:
            net.circ_pump_mass.loc[self.cpm_idxs, 'priorityCtrl'] = self.priority_list

        # Initialize temporary column indicating convergence for prioritized circ pump mas components in cascaded power supply control
        net.circ_pump_mass.loc[self.cpm_idxs, 'converged'] = False



    def control_step(self, net):
        """
        Function for control step of controller **ppi_CircPumpMassPthermalCtrl_limited**\n
        """

        if not self.convergence:

            # Check if all circ pump mass components with lower value in priority_list (-> higher priority) have reached convergence
            if (len(self.higherPriority_idxs) > 0):
                if all(net.circ_pump_mass.loc[self.cpm_idxs[self.higherPriority_idxs], 'converged']):
                    runCtrl = True
                else:
                    runCtrl = False

            elif self.highest_priority:
                    runCtrl = True

            if runCtrl:
                self.iterations += 1
                
                # Extract state of circ_pump_pressure component
                _, current_mdot_circPumpPressure = self.get_Pth_mdot_circPumpPressure(net = net)

                # Extract state of circ_pump_mass component
                current_Pth, current_mdot, current_tempdiff, cp = self.get_Pth_mdot_tempdiff(net, idxs = self.circPumpMass_idx)

                current_mdot_potential = current_mdot_circPumpPressure - self.mdotmin_circPumpPressure # Calculate difference between current mass flow in peak load producer and the specified residual mass flow

                Pth_error = self.target_Pth - current_Pth # Calculate difference between target thermal power and current thermal power
                mdot_target = self.target_Pth / (current_tempdiff * cp)
                mdot_error = mdot_target - current_Pth / (current_tempdiff * cp)

                if current_mdot_circPumpPressure > 0: # Ensure positive mass flow rate through circ_pump_pressure component

                    new_mdot = max(self.mdotmin_circPumpMass, current_mdot + min(current_mdot_potential, mdot_error) * self.proportional_gain)
                    net.circ_pump_mass['mdot_flow_kg_per_s'].at[self.circPumpMass_idx] = new_mdot
                    net.flow_control['controlled_mdot_kg_per_s'].at[self.flow_controller_idx] = new_mdot

                else:
                    new_mdot = max(self.mdotmin_circPumpMass, current_mdot + current_mdot_potential * 2) # Make sure circ_pump_pressure component leaves region with negative mass flow
                    net.circ_pump_mass['mdot_flow_kg_per_s'].at[self.circPumpMass_idx] = new_mdot
                    net.flow_control['controlled_mdot_kg_per_s'].at[self.flow_controller_idx] = new_mdot

        return super(ppi_CircPumpMassPthermalCtrl_limited, self).control_step(net)
    


    def is_converged(self, net):

        # Extract current results
        current_Pth, _, _, _ = self.get_Pth_mdot_tempdiff(net, idxs = self.circPumpMass_idx)

        # Extract current state of circ_pump_pressure component
        _, current_mdot_circPumpPressure = self.get_Pth_mdot_circPumpPressure(net = net)

        # Calculate potential step for altered mass flow
        current_mdot_potential = current_mdot_circPumpPressure - self.mdotmin_circPumpPressure

        Pth_error = self.target_Pth - current_Pth

        if current_mdot_circPumpPressure > 0: # Ensure positive mass flow rate through circ_pump_pressure component
            if (Pth_error <= self.abs_tol) or (abs(current_mdot_circPumpPressure/self.mdotmin_circPumpPressure - 1) <= 0.5) or (abs(current_mdot_potential/(2e-03)-1) <= 0.5):

                self.target_Pth_limited = current_Pth
                
                net.circ_pump_mass.at[self.circPumpMass_idx, 'converged'] = True
                self.convergence = True
                convergence = True
            else:
                convergence = False
        else:
            convergence = False

        return convergence
       

    def get_Pth_mdot_tempdiff(self, net, idxs):
        """
        Function to extract necessary results from result DataFrame of pandapipes network component.\n

        :return: tuple of four floats denoting (current thermal power, current mass flow (from), current temperature difference (to - from), value of specific heat capacity) at specified index positions of circ_pump_mass_components
        """

        mdot_from = net.res_circ_pump_mass['mdot_from_kg_per_s'].at[idxs]
        temp_from_junction = net.res_circ_pump_mass['t_from_k'].at[idxs]
        temp_to_junction = net.res_circ_pump_mass['t_to_k'].at[idxs]
        tempdiff = temp_to_junction - temp_from_junction
        cp = net.fluid.get_heat_capacity((temp_to_junction+temp_from_junction)/2)

        Pth = mdot_from * cp * tempdiff

        return Pth, mdot_from, tempdiff, cp
    
    def get_Pth_mdot_circPumpPressure(self, net):

        mdot_from = net.res_circ_pump_pressure['mdot_from_kg_per_s'].at[self.circPumpPressure_idx]
        temp_from_junction = net.res_circ_pump_pressure['t_from_k'].at[self.circPumpPressure_idx]
        temp_to_junction = net.res_circ_pump_pressure['t_to_k'].at[self.circPumpPressure_idx]

        cp = net.fluid.get_heat_capacity((temp_to_junction+temp_from_junction)/2)
        Pth = mdot_from * cp * (temp_to_junction - temp_from_junction)

        return Pth, mdot_from
    
class ppi_JunctionsMinAbsolutePressureCtrl(BasicCtrl):
    """
    Class of controller in pandapipes network for...
    - controlling minimum allowed absolute pressure at junctions in the network.\n
    - varying the pressure level at the flow junction of the component circ_pump_const_pressure at index **circPumpPressure_idx**.\n
    The junctions which shall be included in the overall control of minimum pressure level can be specified by their respective indices **junction_idxs**. If None, all junction in the network will be checked.\n
    Inherits from class **BasicCtrl** in pandapower.control.basic_controller.\n

    Prerequisits for the use of this controller:
    - At least one component of circ_pump_const_pressure has to be included in the network as a production site defining the network pressure (type "pt").\n
    """

    def __init__(self, net, circPumpPressure_idx:int, target_pmin_bar:float = 1, proportional_gain:float = 0.5, abs_tol:float = 0.05, junction_idxs:list = None, index = None, order = 1, level = 1, **kwargs):
        """
        Init function for controller of type ppi_JunctionsMinAbsolutePressureCtrl\n

        :param net: pandapipes network object\n
        :param circPumpPressure_idx: int, denoting the index of the circ_pump_pressure component in the pandapipes network at which the pressure level shall be varied to reach the target minimum pressure.\n
        :param target_pmin_bar: float, denoting the target value for minimum absolute pressure at junctions **junction_idxs** in the network, defaults to 1\n
        :param proportional_gain: float, denoting a proportional gain for the control step, defaults to 0.5\n
        :param abs_tol: float, denoting the absolute tolerance for convergence of the controller, defaults to 0.05\n
        :param junction_idxs: list-like, denoting the indices of junctions at which the prerssure shall be observed. If None, all junctions in the network are included, defaults to None\n
        :param index: int, defaults to None\n
        :param order: int, defaults to 1\n
        :param level: int, defaults to 1\n
        """

        super(ppi_JunctionsMinAbsolutePressureCtrl, self).__init__(net, **kwargs)

        self.target_pmin_bar = target_pmin_bar
        self.circPumpPressure_idx = circPumpPressure_idx
        self.junction_idxs = junction_idxs if junction_idxs is not None else net.junction.index
        self.iterations = 0
        self.proportional_gain = proportional_gain
        self.abs_tol = abs_tol

        if index is None and "controller" in net.keys():
            index = get_free_id(net.controller)

        self.index = self.add_controller_to_net(net = net, in_service = True, initial_run = True, index = index, order = order, level = level, recycle = False, overwrite = True, drop_same_existing_ctrl = True, **kwargs)

        # Save initial conditions of network
        ## During control step, p_flow_bar at circ_pump_pressure components is changed
        net.circ_pump_pressure['p_flow_bar_init'] = net.circ_pump_pressure['p_flow_bar'].copy()
        


    def add_controller_to_net(self, net, in_service, initial_run, order, level, index, recycle,
                              drop_same_existing_ctrl, overwrite, **kwargs):
        """
        adds the controller to net['controller'] dataframe.

        INPUT:
            **in_service** (bool) - in service status

            **order** (int) - order

            **index** (int) - index

            **recycle** (bool) - if controller needs a new bbm (ppc, Ybus...) or if it can be used \
                                 with prestored values. This is mostly needed for time series \
                                 calculations

        """
        if drop_same_existing_ctrl:
            drop_same_type_existing_controllers(net, type(self), index=index, **kwargs)
        else:
            log_same_type_existing_controllers(net, type(self), index=index, **kwargs)

        # use base class method to raise an error if the object is in DF and overwrite = False
        # if the index is None, the base class is in charge of obtaining the next free index in the data frame
        fill_dict = {"in_service": in_service, "initial_run": initial_run, "recycle": recycle,
                     "order": order, "level": level}
        added_index = super().add_to_net(net=net, element='controller', index=index, overwrite=overwrite,
                           fill_dict=fill_dict, preserve_dtypes=True)
        return added_index
    

    def control_step(self, net):

        self.iterations += 1        

        # Get current minimum absolute pressure in the network
        pmin_current, _ = self.get_pmin(net = net, idxs = self.junction_idxs)        
        pmin_error = self.target_pmin_bar - pmin_current

        current_pflow = net.circ_pump_pressure['p_flow_bar'].at[self.circPumpPressure_idx]
        new_pflow = current_pflow + pmin_error * self.proportional_gain
        net.circ_pump_pressure['p_flow_bar'].at[self.circPumpPressure_idx] = new_pflow

        return super(ppi_JunctionsMinAbsolutePressureCtrl, self).control_step(net)
    
    
    def is_converged(self, net):

        # Extract current results
        pmin_current, idx_pmin_current = self.get_pmin(net = net, idxs = self.junction_idxs)        
        pmin_error = self.target_pmin_bar - pmin_current

        if abs(pmin_error) <= self.abs_tol:
            self.junction_idx_pmin = idx_pmin_current
            convergence = True
        else:
            convergence = False

        return convergence


    def get_pmin(self, net, idxs):

        subDF = net.res_junction.loc[idxs, 'p_bar']

        zipped = list(zip(subDF.values, subDF.index))
        pmin, idx_min = min(zipped, key=lambda x: x[0])

        return pmin, idx_min    
    
class ppi_HeatConsumersMinDiffPressureCtrl(BasicCtrl):
    """
    Class of controller in pandapipes network for...
    - controlling minimum allowed pressure differnce at heat consumers in the network.\n
    - varying the pressure lift at the component circ_pump_const_pressure at index **circPumpPressure_idx**.\n
    The heat consumers which shall be included in the overall control of minimum pressure level can be specified by their respective indices **heatConsumer_idxs**. If None, all heat consumers in the network will be checked.\n
    Inherits from class **BasicCtrl** in pandapower.control.basic_controller.\n

    Prerequisits for the use of this controller:
    - At least one component of circ_pump_const_pressure has to be included in the network as a production site defining the network pressure (type "pt").\n
    """

    def __init__(self, net, circPumpPressure_idx:int, target_dpmin_bar:float = 1, proportional_gain:float = 0.5, abs_tol:float = 0.05, heatConsumer_idxs:list = None, index = None, order = 1, level = 1, **kwargs):

        """
        Init function for controller of type ppi_HeatConsumersMinDiffPressureCtrl\n

        :param net: pandapipes network object\n
        :param circPumpPressure_idx: int, denoting the index of the circ_pump_pressure component in the pandapipes network at which the pressure level shall be varied to reach the target minimum pressure.\n
        :param target_dpmin_bar: float, denoting the target value for minimum pressure difference at heat consumers **heatConsumer_idxs** in the network (positive, p_from_bar - p_to_bar), defaults to 1\n
        :param proportional_gain: float, denoting a proportional gain for the control step, defaults to 0.5\n
        :param abs_tol: float, denoting the absolute tolerance for convergence of the controller, defaults to 0.05\n
        :param heatConsumer_idxs: list-like, denoting the indices of heat consumers at which the prerssure differnce shall be observed. If None, all heat consumers in the network are included, defaults to None\n
        :param index: int, defaults to None\n
        :param order: int, defaults to 1\n
        :param level: int, defaults to 1\n
        """

        super(ppi_HeatConsumersMinDiffPressureCtrl, self).__init__(net, **kwargs)

        self.target_dpmin_bar = target_dpmin_bar
        self.circPumpPressure_idx = circPumpPressure_idx
        self.heatConsumer_idxs = heatConsumer_idxs if heatConsumer_idxs is not None else net.heat_consumer.index
        self.iterations = 0
        self.proportional_gain = proportional_gain
        self.abs_tol = abs_tol

        if index is None and "controller" in net.keys():
            index = get_free_id(net.controller)

        self.index = self.add_controller_to_net(net = net, in_service = True, initial_run = True, index = index, order = order, level = level, recycle = False, overwrite = True, drop_same_existing_ctrl = True, **kwargs)

        # Save initial conditions of network
        ## During control step, plift_bar at circ_pump_pressure components is changed
        net.circ_pump_pressure['plift_bar_init'] = net.circ_pump_pressure['plift_bar'].copy()


    def add_controller_to_net(self, net, in_service, initial_run, order, level, index, recycle,
                              drop_same_existing_ctrl, overwrite, **kwargs):
        """
        adds the controller to net['controller'] dataframe.

        INPUT:
            **in_service** (bool) - in service status

            **order** (int) - order

            **index** (int) - index

            **recycle** (bool) - if controller needs a new bbm (ppc, Ybus...) or if it can be used \
                                 with prestored values. This is mostly needed for time series \
                                 calculations

        """
        if drop_same_existing_ctrl:
            drop_same_type_existing_controllers(net, type(self), index=index, **kwargs)
        else:
            log_same_type_existing_controllers(net, type(self), index=index, **kwargs)

        # use base class method to raise an error if the object is in DF and overwrite = False
        # if the index is None, the base class is in charge of obtaining the next free index in the data frame
        fill_dict = {"in_service": in_service, "initial_run": initial_run, "recycle": recycle,
                     "order": order, "level": level}
        added_index = super().add_to_net(net=net, element='controller', index=index, overwrite=overwrite,
                           fill_dict=fill_dict, preserve_dtypes=True)
        return added_index


    def control_step(self, net):

        self.iterations += 1

        dpmin_current, idx_min_current = self.get_dpmin(net = net, idxs = self.heatConsumer_idxs)
        dpmin_error = self.target_dpmin_bar - dpmin_current
        
        plift_current = net.circ_pump_pressure['plift_bar'].at[self.circPumpPressure_idx]
        plift_new = plift_current + dpmin_error * self.proportional_gain
        net.circ_pump_pressure['plift_bar'].at[self.circPumpPressure_idx] = plift_new

        return super(ppi_HeatConsumersMinDiffPressureCtrl, self).control_step(net)
    

    def is_converged(self, net):

        # Extract current results
        dpmin_current, idx_min_current = self.get_dpmin(net = net, idxs = self.heatConsumer_idxs)
        dpmin_error = self.target_dpmin_bar - dpmin_current

        if abs(dpmin_error) <= self.abs_tol:
            self.heatConsumer_idx_min = idx_min_current
            convergence = True
        else:
            convergence = False

        return convergence
    

    def get_dpmin(self, net, idxs):

        subDF = net.res_heat_consumer.loc[idxs]

        dps = []        
        for idx, p_from, p_to, qext_w in zip(subDF.index, subDF['p_from_bar'], subDF['p_to_bar'], subDF['qext_w']):
            if qext_w != 0:
                dp = p_from - p_to
                dps.append((dp, idx))

        dpmin, idx_min = min(dps, key=lambda x: x[0])

        return dpmin, idx_min