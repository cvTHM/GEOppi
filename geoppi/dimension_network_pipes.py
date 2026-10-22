import pandapipes as ppi
import numpy as np
import pandas as pd
import logging
from geoppi.auxFunctions import _print_info
from typing import Callable, Optional

from geoppi.auxFunctions import (cycleEdges, func_swameejain, func_nikuradse, simultaneity_rh, simultaneity_DHW)


def hydraulicDimensioningNetwork_singleLoadPoint(
        net: ppi.pandapipesNet, # Pandapipes network
        networkType:str = 'KMR',
        dictNominalInnerDiameter:dict = {
            'distributionPipes':{
                    40: 41.9,
                    50: 53.9,
                    65: 69.7,
                    80: 82.5,
                    100: 107.1,
                    125: 132.5,
                    150: 160.3,
                    175: 184.7,
                    200: 210.1,
                    225: 234.5,
                    250: 263,
                    300: 321.7,
                    350: 344.4,
                    400: 393.8,
                    450: 444.4,
                    500: 495.4,
                    550: 546.2,
                    600: 595.8
                },
            'connectionPipes':{
                    10: 17.267,
                    20: 21.7,
                    25: 27.3,
                    32: 36,
                    40: 41.9,
                    50: 53.9,
                    65: 69.7,
                    80: 82.5,
                    100: 107.1,
                    125: 132.5,
                }
        },
        dictWallRoughness:dict = {
            10: 0.0469,
            20: 0.0469,
            25: 0.0469,
            32: 0.0469,
            40: 0.0469,
            50: 0.0469,
            65: 0.0469,
            80: 0.0469,
            100: 0.0469,
            125: 0.0469,
            150: 0.0469,
            175: 0.0469,
            200: 0.0469,
            225: 0.0469,
            250: 0.0469,
            300: 0.0469,
            350: 0.0469,
            400: 0.0469,
            450: 0.0469,
            500: 0.0469,
            550: 0.0469,
            600: 0.0469
        },
        subsetDistribution:list = [100, 150],
        subsetConnection:list = [25, 32, 40],
        connectionTypeAttr:str = 'connectionType',
        elongationFactorPipes:float = 1,
        nominalWidthAttr:str = 'nominalWidth',
        dpSpec:float = 80,
        thresh_percentage:list = [0.4, 0.2],
        respectHeight:bool = False,
        closedValvesIdx:list = [],
        frictionModel:str = 'swamee-jain',
        radialSystem:bool = True,
        define_loops_with_opened_valves:bool = True,
        fluidProperties:dict = {'rho':1000, 'cp':4190, 'nu':0.413e-06},
        greaterWidthFinalStep:dict = {'distributionPipes':0, 'connectionPipes':0},
        deltaT: Optional[float] = None,
        simultaneitySpaceHeating: Optional[Callable[[float], float]] = None,
        simultaneityHotWater: Optional[Callable[[float], float]] = None,
        **kwargs
    ):

    """
    Function for automatic dimensioning of pandapipes network with target value for specific pressure loss in pipes.

    :param net: denotes a pandapies network with only ONE LAYER (usually feed line)
    :param networkType: denotes the type of the network. Currently, "KMR" and "5GDHC" (for uninsulated networks with plastic pipes) are implemented.
    :param dictNominalInnerDiameter: is a nested dictionary with mapping of pipe inner diameters (values) to pipe nominal diameters (keys). Values are given in mm.
        On the first level, the keys "distributionPipes" and "connectionPipes" allow for the differentiation between distribution pipes and house connection pipes.
    :param dictWallRoughness: denotes a dictionary mapping wall roughness (mm) to nominal widths.
    :param subsetDistribution: denotes a list of nominal widhts which shall be considered when iterating through available nominal widths in *dictNominalInnnerDiameter['distributionPipes'].keys()*
    :param subsetConnection: ... the same principle applied to connection pipes.
    :param connectionTypeAttr: denotes the attribute name in the DataFrame where information on the discrimination between distribution pipes ('distribution') and connection pipes ('houseConnection') are stored.
    :param elongationFactor: is a value which is used for temporary elongatino of the pipes' length attribute. This is used to consider the influence of armatures etc. on the pressure loss which is modelled as an altered pipe length.
    :param nominalWidthAttr: denotes the attribute name in the DataFrame where values for nominal widths are stored.
    :param dpSpec: denotes the target value for length-specific pressure loss in single pipes (Pa/m) for the iterative dimensioning process.
    :param thresh_percentage: is a lower and an upper limit for the allowed relative undercutting / exceedance of the target pressure loss value in pipes.
        If a pipe's specific pressure loss falls below *1-thresh_percentage[0]* or exceeds by *1+thresh_percentage[1]* the target value for pressure losses, its nominal width is adapted in the next iteration step.
    :param respectHeight: denotes a switch whether geodetic heights shall be considered in the dimensioning process.
        The default value False should be used for closed networks as length-specific pressure losses are the target value for dimensioning.
    :param closedValvesIdx: denotes a list of indices for valves which shall remain open (controlling loops in the network).
        If it is set to [-1], all existing valves are closed!
    :param frictionModel: denotes the desired friction model which shall be used in hydraulic calculations (either 'colebrook' or 'nikuradse' are implemented).
    :param radialSystem: denotes a switch whether the algorithm shall check for existing loops in the network. 
        If a pure radial system is present, the dimensioning is reduced to a direct calculation of the matching / next higher nominal width.\n
        In mixed networks, only bridges with producers on exactly one side are sized directly.
        Producer connections and loop pipes retain iterative sizing. False disables direct sizing.
    :param define_loops_with_opened_valves: controls whether the definition of loops shall include valves that are opened.
    :param fluidProperties: denotes a dictionary with information on used fluid properties in the network (density "rho", specific heat capacity "cp" and kinematic viscosity "nu").
    :param greaterWidthFinalStep: denotes a dictionary which allows for the selection of altered nominal widths after the algorithm found best matches for the tagret pressure loss.
        The steps for distribution pipes and for connection pipes refer to steps within the provided range of keys for nominal widths.

        Directly sized radial pipes retain their smallest complying nominal width.

    Separate power columns are configured during create_ppi_network_from_gdf and read
    exclusively from net.sink via net['design_power_attrs'] ('rh' and 'DHW').
    With no configured columns, radial sizing uses simulated pipe mass flows.

    :param deltaT: float positive feed-minus-return temperature difference in K, required for separate powers.
    :param simultaneitySpaceHeating: callable taking downstream building count; default is constant 1.
    :param simultaneityHotWater: callable taking downstream building count; default is 1.264 * n**(-0.82), without clipping.
        Diversity applies only to radial pipe design flows, not to the network's sink/source simulation loads.
        Pipe columns design_power_kW, design_mdot_kg_per_s and design_Pa_per_m describe this local design case.

    :returns: pandapipes network object
    """

    ### Imports
    import copy

    net = net.deepcopy()
    net.pipe.reset_index(drop = True, inplace = True)


    ### Plausbility checks
    if networkType not in ('KMR', '5GDHC'):
        logging.warning('\n### Choose proper "networkType". Either "KMR" or "5GDHC" are implemented types. Aborting... ###')
        return net
    
    if frictionModel not in ('colebrook', 'swamee-jain', 'nikuradse'):
        logging.warning('\n### Choose proper frictionModel for hydraulic caculations. Implemented values are "colebrook", "swamee-jain" and "nikuradse". Aborting... ###')
        return net
    
    if any([nom not in dictNominalInnerDiameter['distributionPipes'].keys() for nom in subsetDistribution]):
        logging.warning('\n### Selected subset of nominal widths for distribution pipes could not be found in provided dictionary of nominal widths. Aborting... ###')
        return net

    if any([nom not in dictNominalInnerDiameter['connectionPipes'].keys() for nom in subsetConnection]):
        logging.warning('\n### Selected subset of nominal widths for house connection pipes could not be found in provided dictionary of nominal widths. Aborting... ###')
        return net
    
    ### Preparations
    # Define simulation options from input
    user_pf_options = {'mode':'hydraulics', 'friction_model':frictionModel, 'reset':True,'quit_on_inconsistency_connectivity':True}
    ppi.set_user_pf_options(net = net, **user_pf_options)

    # Get indices of selected nominal widths for distribution pipes and house connection pipes in provided dictionaries
    dictNominalInnerDiameter_red = copy.deepcopy(dictNominalInnerDiameter)

    # Create reduced nominal width dictionary for the used/selected subset of nominal widths only
    dictNominalInnerDiameter_red['distributionPipes'] = {k:v for k,v in dictNominalInnerDiameter['distributionPipes'].items() if k in subsetDistribution}
    dictNominalInnerDiameter_red['connectionPipes'] = {k:v for k,v in dictNominalInnerDiameter['connectionPipes'].items() if k in subsetConnection}

    # Initialisation of nominal widths for pipes
    initNWDistribution = subsetDistribution[len(subsetDistribution)//2]
    initNWConnection = subsetConnection[len(subsetConnection)//2]

    net.pipe[nominalWidthAttr]      = net.pipe.apply(lambda pipe: initNWDistribution if pipe[connectionTypeAttr] == 'distribution' else initNWConnection, axis = 1)
    net.pipe['diameter_m']          = net.pipe.apply(lambda x: 1/1e3 * (dictNominalInnerDiameter['distributionPipes'][x[nominalWidthAttr]] if x[connectionTypeAttr] == 'distribution' else dictNominalInnerDiameter['connectionPipes'][x[nominalWidthAttr]]), axis = 1)
    net.pipe['k_mm']                = net.pipe.apply(lambda x: dictWallRoughness[x[nominalWidthAttr]], axis = 1)

    # Initialisation of pipe indices for different connection types - get indices of distribtuion pipes and of connection pipes in net.pipe DataFrame
    distributionPipes_Idx   = np.array(net.pipe[net.pipe[connectionTypeAttr] == 'distribution'].index)
    connectionPipes_Idx     = np.array(net.pipe[net.pipe[connectionTypeAttr] == 'houseConnection'].index)


    # Control of closed valves
    if hasattr(net, 'valve'):
        if list(closedValvesIdx) == [-1]:
            net.valve['opened'] = False
        elif len(closedValvesIdx) > 0:
            net.valve.loc[closedValvesIdx, 'opened'] = False
        elif len(closedValvesIdx) == 0:
            net.valve['opened'] = True
        else:
            pass



    # Classify topology and sum consumers along geometrically weighted source paths.

    # The network owns both the consumer powers and their column-name mapping.
    # Reject the former GIS arguments instead of silently ignoring them via kwargs.
    obsolete = {'buildings', 'buildingIdAttr', 'spaceHeatingPowerAttr', 'hotWaterPowerAttr'} & kwargs.keys()
    if obsolete:
        raise TypeError('Configure building power attributes in create_ppi_network_from_gdf, not pipe sizing.')
    power_attrs = net.get('design_power_attrs', {})
    separate_powers = bool(power_attrs)
    sink_attributes = {'n_buildings_downstream': pd.Series(1.0, index=net.sink.index)}
    if separate_powers:
        if set(power_attrs) != {'rh', 'DHW'}:
            raise ValueError("net['design_power_attrs'] must contain 'rh' and 'DHW'.")
        if deltaT is None or not np.isfinite(deltaT) or deltaT <= 0:
            raise ValueError('Separate design powers require a positive finite deltaT in K.')
        spaceHeatingPowerAttr, hotWaterPowerAttr = power_attrs['rh'], power_attrs['DHW']
        power_data = net.sink  # No building layer or repeated spatial/ID matching is needed.
        for column, output in ((spaceHeatingPowerAttr, 'space_heating_kW_downstream'),
                               (hotWaterPowerAttr, 'hot_water_kW_downstream')):
            values = pd.to_numeric(power_data[column], errors='raise')
            active_values = values.loc[net.sink.in_service]
            if not np.isfinite(active_values).all() or (active_values < 0).any():
                raise ValueError('Active consumers require finite nonnegative design powers.')
            sink_attributes[output] = values

    net, bool_cycle = check_for_cycle(
        ppi_network=net, assign_cycle_affiliation=True,
        define_loops_with_opened_valves=define_loops_with_opened_valves,
        length_attr='length_km', radial_sink_attributes=sink_attributes)
    radial_mask = net.pipe['is_radial'].to_numpy() & bool(radialSystem)
    niter = 1 if radial_mask[net.pipe.in_service.to_numpy()].all() else 50

    # Calculate power at design load point for pipes outside of cycles and outside direct connections between multiple producers
    if separate_powers:
        heating_factor = simultaneitySpaceHeating or simultaneity_rh
        water_factor = simultaneityHotWater or simultaneity_DHW
        net.pipe['design_power_kW'] = np.nan
        for idx in net.pipe.index[radial_mask]:
            count = net.pipe.at[idx, 'n_buildings_downstream']
            factors = np.array([heating_factor(count), water_factor(count)]) if count > 0 else np.zeros(2)
            if not np.isfinite(factors).all() or (factors < 0).any():
                raise ValueError('Simultaneity factors must be finite and nonnegative.')
            net.pipe.at[idx, 'design_power_kW'] = (
                factors[0] * net.pipe.at[idx, 'space_heating_kW_downstream']
                + factors[1] * net.pipe.at[idx, 'hot_water_kW_downstream'])

    # Control of junction heights
    if respectHeight == False:
        heights = net.junction['height_m'].copy()
        net.junction['height_m'] = 0

    # Apply elongation factor to pipe lengths
    net.pipe['length_km'] *= elongationFactorPipes
   

    ### Calculations / Iterative dimensioning
    # Iteratively adapt nominal width and corresponding pipe features (inner diameter, roughness, heat transfer coefficient) to meet criteria based on specific pressure loss @ max. load
    for ii in range(niter):

        ## Perform simulation step
        ppi.pipeflow(net)

        net = calc_spec_dp(
            net = net, 
            colName_length = 'length_km',
            outputCol = 'Pa_per_m',
            factor_l = 1/elongationFactorPipes, 
            cp = fluidProperties['cp'], 
            rho = fluidProperties['rho']
            )
        
        new_keys = np.array(
                net.pipe.apply(lambda x: list(dictNominalInnerDiameter_red['distributionPipes'].keys()).index(x[nominalWidthAttr]) if x[connectionTypeAttr] == 'distribution'
                else list(dictNominalInnerDiameter_red['connectionPipes'].keys()).index(x[nominalWidthAttr]) if x[connectionTypeAttr] == 'houseConnection'
                else len(dictNominalInnerDiameter_red['distributionPipes'].keys()) - 2, axis=1)
            )         

        # Direct sizing runs once; iterative steps below are masked for these pipes.
        if ii == 0 and radial_mask.any():
            func_dp = func_swameejain if frictionModel in ('swamee-jain', 'colebrook') else func_nikuradse

            if dpSpec <= 0 or not np.isfinite(dpSpec):
                raise ValueError('dpSpec must be finite and positive.')
            net.pipe['design_mdot_kg_per_s'] = np.nan
            net.pipe['design_Pa_per_m'] = np.nan
            for idx in net.pipe.index[radial_mask]:
                mass_flow = (net.pipe.at[idx, 'design_power_kW'] * 1000 / (fluidProperties['cp'] * deltaT)
                             if separate_powers else abs(net.res_pipe.at[idx, 'mdot_from_kg_per_s']))
                if not np.isfinite(mass_flow):
                    raise ValueError(f'No finite design flow for radial pipe {idx}.')
                pipe_type = 'distributionPipes' if net.pipe.at[idx, connectionTypeAttr] == 'distribution' else 'connectionPipes'
                candidates = sorted(dictNominalInnerDiameter_red[pipe_type].items(), key=lambda item: item[1])
                # Evaluate each real catalogue diameter with its own roughness; no root rounding is needed.
                for nominal, diameter_mm in candidates:
                    # Calculate pressure loss as diffrence to set pressure loss of zero
                    pressure_loss = (0.0 if mass_flow <= 1e-12 else -float(func_dp(
                        diameter_mm / 1000, mass_flow, 0, fluidProperties['rho'],
                        fluidProperties['nu'], dictWallRoughness[nominal])))
                    if pressure_loss <= dpSpec:
                        break
                else:
                    raise ValueError(f'No allowed nominal width meets dpSpec for radial pipe {idx}.')
                net.pipe.at[idx, nominalWidthAttr] = nominal
                net.pipe.at[idx, 'diameter_m'] = diameter_mm / 1000
                net.pipe.at[idx, 'k_mm'] = dictWallRoughness[nominal]
                net.pipe.at[idx, 'design_mdot_kg_per_s'] = mass_flow
                net.pipe.at[idx, 'design_Pa_per_m'] = pressure_loss
                new_keys[idx] = list(dictNominalInnerDiameter_red[pipe_type]).index(nominal)

        if not radial_mask.all():
            ### Identify pipes whose resulting specific pressure drop is out of boundaries
            dp = abs(np.array(net.res_pipe['Pa_per_m']))
            np.nan_to_num(dp, copy = False, nan = 0)

            # Create index steps to greater / smaller nominal width in provided dictionaries
            steps = np.array([
                1 if ( (dp[npipe] - dpSpec) / dpSpec > thresh_percentage[1])
                else -1 if ( (dp[npipe] - dpSpec) / dpSpec < -1 * thresh_percentage[0])
                else 0
                for npipe in range(len(net.pipe))
            ])

            steps[radial_mask] = 0
            new_keys += steps

            new_keys[distributionPipes_Idx] = np.clip(new_keys[distributionPipes_Idx], 0, len(dictNominalInnerDiameter_red['distributionPipes'].keys())-1)
            new_keys[connectionPipes_Idx] = np.clip(new_keys[connectionPipes_Idx], 0, len(dictNominalInnerDiameter_red['connectionPipes'].keys())-1)

            # Assign newly determined nominal widths to pipes
            net.pipe.loc[distributionPipes_Idx, nominalWidthAttr] = [list(dictNominalInnerDiameter_red['distributionPipes'].keys())[kk] for kk in new_keys[distributionPipes_Idx]]
            net.pipe.loc[connectionPipes_Idx, nominalWidthAttr] = [list(dictNominalInnerDiameter_red['connectionPipes'].keys())[kk] for kk in new_keys[connectionPipes_Idx]]

            # Assign corresponding inner diameters to pipes
            net.pipe['diameter_m'] = net.pipe.apply(lambda x: 
                                                    1/1e3 * dictNominalInnerDiameter['distributionPipes'][x[nominalWidthAttr]] if x[connectionTypeAttr] == 'distribution'
                                                    else 1/1e3 * dictNominalInnerDiameter['connectionPipes'][x[nominalWidthAttr]] if x[connectionTypeAttr] == 'houseConnection'
                                                    else 0.05, axis= 1)      

            net.pipe['k_mm'] = net.pipe.apply(lambda x:
                                              dictWallRoughness[x[nominalWidthAttr]], axis = 1)
        

    radial_keys = new_keys.copy()
    # Apply final step to greater / lower nominal width dependent on input
    if greaterWidthFinalStep['distributionPipes'] != 0:
        new_keys[distributionPipes_Idx] += greaterWidthFinalStep['distributionPipes']
        new_keys[distributionPipes_Idx] = np.clip(new_keys[distributionPipes_Idx], 0, len(dictNominalInnerDiameter_red['distributionPipes'].keys())-1)

    if greaterWidthFinalStep['connectionPipes'] != 0:
        new_keys[connectionPipes_Idx] += greaterWidthFinalStep['connectionPipes']
        new_keys[connectionPipes_Idx] = np.clip(new_keys[connectionPipes_Idx], 0, len(dictNominalInnerDiameter_red['connectionPipes'].keys())-1)

    new_keys[radial_mask] = radial_keys[radial_mask]
    # Assign final nominal widths and diameters
    net.pipe.loc[distributionPipes_Idx, nominalWidthAttr] = [list(dictNominalInnerDiameter_red['distributionPipes'].keys())[kk] for kk in new_keys[distributionPipes_Idx]]
    net.pipe.loc[connectionPipes_Idx, nominalWidthAttr] = [list(dictNominalInnerDiameter_red['connectionPipes'].keys())[kk] for kk in new_keys[connectionPipes_Idx]]

    
    net.pipe['diameter_m'] = net.pipe.apply(lambda x: 
                                            1/1e3 * dictNominalInnerDiameter['distributionPipes'][x[nominalWidthAttr]] if x[connectionTypeAttr] == 'distribution'
                                            else 1/1e3 * dictNominalInnerDiameter['connectionPipes'][x[nominalWidthAttr]] if x[connectionTypeAttr] == 'houseConnection'
                                            else 0.05, axis= 1)
    
    # Assign final wall roughness
    net.pipe['k_mm'] = net.pipe.apply(lambda x:
                                              dictWallRoughness[x[nominalWidthAttr]], axis = 1)
    

    ## Perform final simulation step
    ppi.pipeflow(net)

    net = calc_spec_dp(
        net = net, 
        colName_length = 'length_km',
        outputCol = 'Pa_per_m',
        factor_l = 1/elongationFactorPipes, 
        cp = fluidProperties['cp'], 
        rho = fluidProperties['rho'])
    
    ### Corrections before returning values
    # Set pipes out of service that feature no mass flow
    mask0 = (abs(net.res_pipe['v_mean_m_per_s']) < 1e-12) | (net.res_pipe['v_mean_m_per_s'].isna())
    net.pipe.loc[mask0, 'in_service'] = False

    # Correct formerly altered pipe length
    net.pipe['length_km'] *= 1/elongationFactorPipes

    # Correct formerly altered junction heights
    if respectHeight == False:
        net.junction['height_m'] = heights
        
    return net

def assign_insulation_type(
        df_pipes:pd.DataFrame,
        colNameNominalWidth:str = 'Nennweite',
        network_type:str = 'KMR', 
        insulation_type:str = 'std',
        colNameHTC:str = 'alpha_w_per_m2k',
        calcFromNominalWidth:bool = True,
        colNameInsulationThickness:str = None,
        colNameWallThickness:str = None,
        colNameJacketThickness:str = None,
        colNameDiameter:str = None,
        keepInfoColumns:bool = False,
        lambdaPE:float = 0.35,
        lambdaSteel:float = 50,
        lambdaInsul:float = 0.027,
        )->pd.DataFrame:

    """
    Function that assigns heat transfer coefficient for pipe insulations.

    :param df_pipes: is a DataFrame of pipe objects (correpsonding to pandapipes network pipe DataFrame)
    :param colNameNominalWidth: denotes the colum name where the float value for nominal width is found
    :param network_type: either 'KMR' or '5GDHC' for differentiation between conventional insulated heating networks (KMR) or heating and cooling networks (5GDHC without insulation)
    :param insulation_type: either 'std', '1x' or 'None' defines the quality and thickness of insulation
    :param colNameHTC: denotes the target column name in which to store values of heattransfer coefficient
    :param calcFromNominalWidth: is a boolean selection whether heat trnasfer coefficient shall be calculated from documented insulation thickness dependent on nominal width. If False, the columns described in the following are used for calculations.
    :param colNameInsulationThickness:  denotes the optional column name where information on the insulation thickness is stored.
    :param colNameWallThickness: denotes the optional column name where information on the medium pipe wall thickness is stored.
    :param colNameJacketThickness: denotes the optional column name where information on the outer plastics jackete wall thickness is stored.
    :param colNameDiameter: denotes the column name where information on the inner diameter of the pipe is stored.
    :param keepInfoColumns: Boolena selection if additional for layer diameters should be stored in the output DataFrame

    :returns: DataFrame of pipes
    """

    # Sources: https://www.isoplus.group/fileadmin/products/IP_Planungshandbuch_DE.pdf, Handbuch Firma Isoplus (2011)

    # Imports
    from geoppi import pipe_characteristics
    import numpy as np

    # Temporary copy
    df = df_pipes.copy()

    # Loading dictionaries with pipe characteristics
    dict_nominal_to_inner_diameters_KMR         = pipe_characteristics.nominal_to_inner_diameters_KMR
    dict_wall_thickness_KMR                     = pipe_characteristics.dict_Wanddicke_KMR
    dict_jacket_thickness_std_KMR               = pipe_characteristics.dict_mantel_thickness_std_KMR
    dict_jacket_thickness_1x_KMR                = pipe_characteristics.dict_mantel_thickness_1x_KMR
    dict_insulation_thickness_std_KMR           = pipe_characteristics.dict_insulation_thickness_std_KMR
    dict_insulation_thickness_1x_KMR            = pipe_characteristics.dict_insulation_thickness_1x_KMR

    dict_nominal_to_inner_diameters_5GDHC       = pipe_characteristics.nominal_to_inner_diameters_5GDHC
    dict_wall_thicknes_5GDHC                    = pipe_characteristics.dict_Wanddicke_5GDHC


    if calcFromNominalWidth:
        try:
            NW=list(df[colNameNominalWidth])
        except:
            raise Exception('\n### No nominal width assigned to pipes! ###')
        

        if network_type == 'KMR':

            if insulation_type == 'None':
                insulation_mm=0
                jacket_mm = [dict_jacket_thickness_std_KMR[nw] for nw in NW ]

            elif insulation_type=='std':
                insulation_mm=[dict_insulation_thickness_std_KMR[nw] for nw in NW ]
                jacket_mm = [dict_jacket_thickness_std_KMR[nw] for nw in NW ]

            elif insulation_type=='1x':
                insulation_mm=[dict_insulation_thickness_1x_KMR[nw] for nw in NW]
                jacket_mm = [dict_jacket_thickness_1x_KMR[nw] for nw in NW ]

            # Wall thickness
            S = [dict_wall_thickness_KMR[nw] for nw in NW]

            # Repair pipes with missing data (diameter is null)
            D=[dict_nominal_to_inner_diameters_KMR[NW[mm]] for mm in range(len(NW))]

        elif network_type == '5GDHC':
            insulation_mm = 0
            jacket_mm = [dict_wall_thicknes_5GDHC[nw] for nw in NW]
            
            # Wall thickness = 0 since only jacket thickness is relevant in PE pipes withour insulation
            S = [0] * len(NW)
            
            # Repair pipes with missing data (diameter is null)
            D=[dict_nominal_to_inner_diameters_5GDHC[NW[mm]] for mm in range(len(NW))]

        else:
            _print_info('\n ### Choose proper network typ. Either "KMR" or "5GDHC". ###')
            return df
        
    
    else: # calcFromNominalWidth == False
        try:
            DD=list(df[colNameInsulationThickness])
        except:
            DD=[np.nan for l in range(len(df))]
            # raise Exception('No insulation thickness column in pipes DataFrame!')
        try:
            D=list(df[colNameDiameter])
        except:
            raise Exception('\n### No diameter column assigned to pipes DataFrame! ###')
            return
        try:
            S=list(df[colNameWallThickness])
        except:
            _print_info('\n### No wall thickness column assigned to pipes DataFrame! ###')


    

    # recalculate diameters from pipe thickness for better readability

    D                           = np.array(D)
    arrLamInsulation            = np.ones(len(D)) * lambdaInsul
    arrLamPipe                  = np.ones(len(D)) * lambdaSteel
    arrLamJAcket                = np.ones(len(D)) * lambdaPE

    try:
        DMO             = D + np.array(S)/1000*2
        DMIO            = DMO+np.array(insulation_mm)/1000*2
        DMOUTER         = DMIO + np.array(jacket_mm)/1000*2

    except:
        raise Exception('There has been an error with assigning insulation diameters to the pipes.')
        
    if keepInfoColumns:

        # Assign results to pipes #    
        df['insulation_mm']                 = insulation_mm
        df['lambda_insulation_W_per_mK']    = lambdaInsul
        df['lambda_mantel_W_per_mK']        = lambdaPE
        df['diameter_m_o']                  = DMO
        df['diameter_m_insul_o']            = DMIO
        df['diameter_m_outer']              = DMOUTER



    ### 3) Calculate heat transfer coefficient
    # assuming the inner wall temperature of the medium pipe is equal to the fluid temperature (alpha=infinite)
    # assuming outer wall temperature outside insulation is equal to soil temperature
    df[colNameHTC] = calculate_heat_transfer_coeff(
        (D,DMO,
        DMO,DMIO,
        DMIO, DMOUTER),
        arrLamPipe,
        arrLamInsulation,
        arrLamJAcket
        )

    # Preserve the selected output name while maintaining both pandapipes aliases.
    df['alpha_w_per_m2k'] = df[colNameHTC]
    df['u_w_per_m2k'] = df[colNameHTC]
    return(df)

def calculate_heat_transfer_coeff(
        diameters:list,
        *lambda_vals):
    
    """
    Function that calculates heat transfer coefficient through multi-layer cylindrical pipe
    ### Calculate heat transfer coefficient referred to inner diameter for each pipe ###
    *diameters* given as tuple of lists from inner diameter to outer diameter (two diameters for each layer!!!)
    *lambda_vals* given as lists (one entry for each pipe), number of lists correspondent to number of layers ###

    """

    import numpy as np
    import math

    diams=[]
    for ii in range(0,int(len(diameters)/2)):
        diams+=[[np.array(diameters[2*ii]),np.array(diameters[2*ii+1])]]

    lambdas=[np.array(m) for m in lambda_vals]

    alpha_i=math.inf*np.ones(np.shape(lambdas[0]))
    alpha_o = 40 # Estimation for outer heat transfer coefficient (W/(m²*K)) from pipe to undisturbed ground temperature
    alphas=[alpha_i,alpha_o]

    k_lambdas=np.zeros(np.shape(lambdas[0]))
    for ii in range(0,len(lambdas)):
        k_lambdas+=(np.log(diams[ii][1]/diams[ii][0])/lambdas[ii])

    k_alphas=np.zeros(np.shape(lambdas[0]))
    for ll in range(0,len(alphas)):
        if ll==0:
            k_alphas+=(1/(alphas[ll]*diams[0][0]/2))
        if ll==len(alphas)-1:
            k_alphas+=(1/(alphas[ll]*diams[-1][1]/2))

    # Reference diameter for heat trasnfer coefficient should be the innermost diameter with contact to the fluid
    k_ref_inner=((k_lambdas+k_alphas)*diams[0][0]/2)**(-1)
    k_ref_outer=((k_lambdas+k_alphas)*diams[-1][1]/2)**(-1)
        

    return(k_ref_inner)

def assign_nominal_widths_manually(
        net: ppi.pandapipesNet, # Pandapipes network
        dictNW:dict, 
        dictNominalInnerDiameter:dict = {
            10: 17.267,
            20: 21.7,
            25: 27.3,
            32: 36,
            40: 41.9,
            50: 53.9,
            65: 69.7,
            80: 82.5,
            100: 107.1,
            125: 132.5,
            150: 160.3,
            175: 184.7,
            200: 210.1,
            225: 234.5,
            250: 263,
            300: 321.7,
            350: 344.4,
            400: 393.8,
            450: 444.4,
            500: 495.4,
            550: 546.2,
            600: 595.8
        },
        dictWallRoughness:dict = {
            10: 0.0469,
            20: 0.0469,
            25: 0.0469,
            32: 0.0469,
            40: 0.0469,
            50: 0.0469,
            65: 0.0469,
            80: 0.0469,
            100: 0.0469,
            125: 0.0469,
            150: 0.0469,
            175: 0.0469,
            200: 0.0469,
            225: 0.0469,
            250: 0.0469,
            300: 0.0469,
            350: 0.0469,
            400: 0.0469,
            450: 0.0469,
            500: 0.0469,
            550: 0.0469,
            600: 0.0469
        },
        nominalWidthAttr:str = 'nominalWidth'
        ):
    
    """
    Function to manually adapt nominal width of pipes in pandapipes network by selection of pipe indices.
    Hydraulic properties 
    - wall roughness
    - inner diameter
    are updated after nominal width is assigned.

    :param net: pandapipes network object
    :param dictNW: dictionary matching keys of desired nominal widths to values of pipe indices
    :param dictNominalInnerDiameter: dictionary matching keys of nominal diameters to values of inner diameters (mm)
    :param dictWallRoughness: dictionary matching keys of nominal diamaters to values of wall roughness (mm)
    :param nominalWidthAttr: str denoting the attribute name of nominal widths found in net.pipe

    :returns: pandapipes network object
    """

    ### Imports
    import numpy as np

    ### Plausibility checks
    net = net.deepcopy()
    
    ### Updating nominal width and associated hydraulic pipe properties
    pipeIdxArr = np.array(net.pipe.index)

    for nw, idxArr in dictNW.items():
        idxArr_ed = np.intersect1d(idxArr, pipeIdxArr)
        mask = net.pipe.index.isin(idxArr_ed)

        net.pipe.loc[mask, nominalWidthAttr] = nw
        net.pipe.loc[mask, 'diameter_m'] = dictNominalInnerDiameter[nw] / (1e3)
        net.pipe.loc[mask, 'k_mm'] = dictWallRoughness[nw]
    

    return net

def calc_spec_dp(
        net: ppi.pandapipesNet, # pandapipes network object
        colName_length:str = 'length_km',
        outputCol:str = 'Pa_per_m',
        factor_l:float = 1.0, 
        rho:float = 1000, 
        g:float = 9.81, 
        **kwargs):
    
    """
    Function that calculates length-specific pressure loss (Pa/m) due to wall friction inside pipe in pandapipes network class.
    Geodetic height differences are omitted in this calculation - only pressure loss due to wall friction and additional loss coefficients in the wall are taken into account.

    :param net: Pandpapipes network class
    :param colName_length: column name where length information is stored in net.pipe DataFrame (unit: km)
    :param factor_l: Factor which is applied to the length of each pipe BEFORE referring pressure loss to pie length.
    :param rho: density of medium, necessary for subtracting pressure differneces due to height differences
    :param g: gravitational constant, necessary for subtracting pressure differneces due to height differences
    :return: Pandapipes network with pipe DataFrame extended by **outputCol** with specific pressure loss per m.
    """

    import numpy as np

    if 'constants' in kwargs:
        constants = kwargs['constant']
        g = constants.g
        rho = constants.rho

    net.res_pipe['diff_h'] =      np.array(net.junction.loc[net.pipe['from_junction'], 'height_m']) - np.array(net.junction.loc[net.pipe['to_junction'], 'height_m'])

    net.res_pipe[outputCol] =     (( (net.res_pipe['p_from_bar'] - net.res_pipe['p_to_bar']) * 10**5 - rho * g * net.res_pipe['diff_h'])) / (
        net.pipe[colName_length] * factor_l * 1e3)

    return(net)

def check_for_cycle(
        ppi_network: ppi.pandapipesNet,
        assign_cycle_affiliation:bool = True, 
        define_loops_with_opened_valves:bool = True,
        length_attr:str = 'length_km',
        radial_sink_attributes: Optional[dict[str, pd.Series]] = None
        ):
    
    """
    Function that checks given input of pandapipes network for cycles in pipe DataFrame.

    :param ppi_network: Pandapipes network
    :param assign_cycle_affiliation: Boolean if the affiliation to the cycle with the biggest summed length shall be transferred to the original pipe DataFrame as additional attribute.
    :param define_loops_with_opened_valves: Boolean if only opened valves shall be used for definition of potential cycles. If False, all valves regardless of their state are considered.
    :param length_attr: Attribute name for pipe length found in ppi_network.pipe DataFrame

    :param radial_sink_attributes: optional dict mapping output columns to Series indexed like net.sink.
        Enables bridge classification and downstream sums; stored pipe direction is never changed.
    :return: ppi_network Pandapipes network with new column inside pipe DataFrame and with ID connecting pipes to possibly existent loops/cycles
    :return: cycle either True if cycle is found or False if no cycle is found.
    """

    # Imports
    import networkx as nx
    import pandas as pd

    # Initialisations
    components = [ppi_network.pipe]

    if hasattr(ppi_network, 'valve'):
        if define_loops_with_opened_valves == True:
            components.append(ppi_network.valve[ppi_network.valve['opened'] == True])
        else:
            components.append(ppi_network.valve)

    if hasattr(ppi_network, 'heat_exchanger'):
        components.append(ppi_network.heat_exchanger)

    componentsDF = pd.concat(tuple(components), axis = 0)
    componentsDF['weight'] = componentsDF[length_attr]
    componentsDF['weight'] = componentsDF['weight'].fillna(0)

    # Check for cycles
    G = nx.from_pandas_edgelist(
        df = componentsDF,
        source = 'from_junction',
        target = 'to_junction',
        edge_attr = 'weight'
        )


    cyclelist = nx.cycle_basis(G)
    cycle = True if len(cyclelist) != 0 else False

    # Get list of summed weights within cycles and sort the list of cycles in cyclelist in corresponding ascending order
    _, cycleWeights = cycleEdges(G = G, cycleNodes = cyclelist, weights = 'weight')
    cyclelist_sorted = [x for _,x in sorted(zip(cycleWeights, cyclelist), reverse = False)]

    # Assign cycle number of the cycle with the highest summed weight to all contained edges
    if assign_cycle_affiliation == True:
        ppi_network.pipe['loop'] = None
        for n, c in enumerate(cyclelist_sorted):
            ppi_network.pipe.loc[(ppi_network.pipe['from_junction'].isin(c) & ppi_network.pipe['to_junction'].isin(c)), 'loop'] = n

    if radial_sink_attributes is not None:
        from geoppi.suitable_network_routing import sum_attrs_along_shortest_paths

        # Preserve parallel pipes and self-loops, which must never be classified as bridges.
        graph = nx.MultiGraph()
        active_nodes = set(ppi_network.junction.index[ppi_network.junction.in_service])
        graph.add_nodes_from(active_nodes)
        for component_name in ('pipe', 'valve', 'heat_exchanger', 'flow_control'):
            component = ppi_network.get(component_name, pd.DataFrame())
            for idx, row in component.iterrows():
                if not row.get('in_service', True):
                    continue
                if component_name == 'valve' and define_loops_with_opened_valves and not row.opened:
                    continue
                u, v = row.from_junction, row.to_junction
                if u not in active_nodes or v not in active_nodes:
                    continue
                # Geometric pipe length has priority over hydraulic equivalent length.
                geometry = row.get('geometry')
                length = (geometry.length if geometry is not None and hasattr(geometry, 'length')
                          else row.get(length_attr, 0) * (1000 if length_attr == 'length_km' else 1))
                if not np.isfinite(length) or length < 0:
                    raise ValueError('Graph lengths must be finite and nonnegative.')
                graph.add_edge(u, v, key=(component_name, idx), weight=length)

        producers = set()
        for component_name in ('ext_grid', 'source'):
            for _, row in ppi_network.get(component_name, pd.DataFrame()).iterrows():
                if row.get('in_service', True) and row.junction in graph:
                    producers.add(row.junction)
        if not producers:
            raise ValueError('Radial classification requires active feed-line producers.')

        ppi_network.pipe['is_radial'] = False
        ## A bridge is a producer connection precisely when both separated sides contain producers.
        # This ensures only those pipes are markes as radial if they are not inside a direct connection between multiple producers (and can be supplied by more than one producer node)
        for u, v in list(nx.bridges(graph)):
            key = next(iter(graph[u][v]))
            attributes = graph[u][v][key].copy()
            graph.remove_edge(u, v, key)
            side_u = nx.node_connected_component(graph, u)
            side_v = nx.node_connected_component(graph, v)
            graph.add_edge(u, v, key=key, **attributes)
            if key[0] == 'pipe' and bool(producers & side_u) != bool(producers & side_v):
                ppi_network.pipe.at[key[1], 'is_radial'] = True

        sinks = ppi_network.sink.loc[ppi_network.sink.in_service]
        reachable = nx.multi_source_dijkstra_path_length(graph, producers, weight='weight')
        if not sinks.junction.isin(reachable).all():
            raise ValueError('Active consumers must be connected to an active producer.')
        for output, values in radial_sink_attributes.items():
            totals = values.reindex(sinks.index).groupby(sinks.junction).sum().to_dict()
            graph, edge_values = sum_attrs_along_shortest_paths(
                graph, startNodes=sorted(producers), endNodes=list(totals), weight='weight',
                output_attr=output, dictAttrsEndNodes=totals)
            ppi_network.pipe[output] = [edge_values.get(('pipe', idx), 0) for idx in ppi_network.pipe.index]

    return(ppi_network, cycle)

def update_ppi_results(
        net: ppi.pandapipesNet, # pandapipes network
        user_pf_options:dict = None,
        respectHeight:bool = False,
        fluidProperties:dict = {'rho':1000, 'cp':4190, 'nu':0.413e-06},
        elongationFactorPipes:float = 1,
        resultAttributes:dict = {'pipe':['Pa_per_m', 'mdot_from_kg_per_s'], 'junction':['t_k', 'p_bar']},
    ):

    """
    Function that updates simulation results for pandapipes network.\n    
    A simulation with the currently entered charatceristics of the network is performed. Results are transferred to the network components' DataFrames.\n 
    :param user_pf_options: dictionary denoting the user-defined pipeflow options which shall be transfered to the pipeflow-command (key frictionModel: denotes the chosen friction model for hydraulic calculations. Either "colebrook", "swamee-jain", "nikuradse" are implemented. Key mode: Calculation mode chosen for simulations. Can either be "hydraulics" or "all".) Defaults to net.user_pf_options if None.\n
    :param respectHeight: denotes a switch whether geodetic height differences between junctions shall be considered in the simulation. Otherwise, all heights are set to zero.\n
    :param fluidProperties: denotes a dictionary with information on used fluid properties in the network (density "rho", specific heat capacity "cp" and kinematic viscosity "nu").\n
    :param elongationFactorPipes: is a value which is used for temporary elongation of the pipes' length attribute. This is used to consider the influence of armatures etc. on the pressure loss which is modelled as an altered pipe length.\n
    :param resultAttributes: is a dictionary with key, value pairs of network components and the respective result entities which shall be transferred to the DataFrame
    
    """

    ### Plausibility checks for user pipeflow options
    if user_pf_options is None:
        user_pf_options = net.user_pf_options

    else:        
        if user_pf_options['friction_model'] not in ('colebrook', 'swamee-jain', 'nikuradse'):
            logging.warning('\n### Choose proper frictionModel for hydraulic calculations. Implemented values are "colebrook", "swamee-jain" and "nikuradse". Aborting... ###')
            return net
        
        if user_pf_options['mode'] not in ('hydraulics', 'all', 'sequential', 'bidirectional'):
            logging.warning('\n### Choose proper calculation mode for pandapipes simulation. Either "hydraulics" for only hydraulic calculation or "all" for subsequent hydraulic and thermalcalculation is implemented. ###')
            return net
    
    # Store user pf options and later restore them
    orig_pf_options = net.user_pf_options.copy()    
    
    ### Preparations
    net = net.deepcopy()
    ppi.set_user_pf_options(net = net, reset = True, **user_pf_options)

    # Control of junction heights
    if respectHeight == False:
        heights = net.junction['height_m'].copy()
        net.junction['height_m'] = 0

    # Apply elongation factor to pipe lengths
    net.pipe['length_km'] *= elongationFactorPipes

    ## Perform simulation step
    ppi.pipeflow(net)

    net = calc_spec_dp(
        net = net, 
        colName_length = 'length_km',
        outputCol = 'Pa_per_m',
        factor_l = 1/elongationFactorPipes, 
        cp = fluidProperties['cp'], 
        rho = fluidProperties['rho'])
    
    # Correct formerly altered pipe length
    net.pipe['length_km'] *= 1/elongationFactorPipes

    # Correct formerly altered junction heights
    if respectHeight == False:
        net.junction['height_m'] = heights

    ### Transfer generated results to component DataFrames
    for k, v in resultAttributes.items():
        if hasattr(net, k):
            comp = net.__getattr__(k)

        if hasattr(net, 'res_'+k):
            compRes = net.__getattr__('res_'+k)
            for val in v:
                if val in compRes.columns:
                    # Correct for different namings in result DataFrame and component DataFrame
                    if val == 'Pa_per_m':
                        comp[val] = abs(compRes[val])

                    else:
                        comp[val] = compRes[val]


    # Restore user pf options
    ppi.set_user_pf_options(net = net, reset = True, **orig_pf_options)
    
    return net

