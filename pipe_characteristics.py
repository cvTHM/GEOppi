
''' Material properties '''

''' Thermal conductivities of pipe and inulsation material (W/(m*K))'''
lambdaPE = 0.35
lambdaSteel = 50
lamdbaPU = 0.027

''' Function definitions '''
def load_hydraulic_pipe_characteristics(
        networkType:str = 'KMR'
    ):
    
    """
    Function to load hydraulic pipe characteristics dependent on input "networkType"

    :param networkType: strr denoting the type of the network pipes, either 'KMR' for plastic jacket pipes with insulation and steel medium pipe or '5GDHC' for pure plastic pipes.
    :returns: dictionaries for matching of nominal widhts to inner diameters and nominal widths to wall roughness values
    """
    
    if networkType == 'KMR':
        dictNominalInnerDiameter = nominal_to_inner_diameters_KMR
        dictWallRoughness = dict_k_mm_KMR

    if networkType == '5GDHC':
        dictNominalInnerDiameter = nominal_to_inner_diameters_5GDHC
        dictWallRoughness = dict_k_mm_5GDHC


    return dictNominalInnerDiameter, dictWallRoughness


''' KMR networks '''

# Dictionary nominal diameters to inner diameters (mm to mm):
nominal_to_inner_diameters_KMR = {
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
}

# Wall thickness (inside insulation, steel etc.) (mm to mm)
dict_Wanddicke_KMR = {
    10: 2.6,
    20: 2.6,
    25: 3.2,
    32: 3.2,
    40: 3.2,
    50: 3.2,
    65: 3.2,
    80: 3.2,
    100: 3.6,
    125: 4,
    150: 4.5,
    175: 4.5,
    200: 6.3,
    225: 6.3,
    250: 6.3,
    300: 7.1,
    350: 8,
    400: 8.8,
    450: 10,
    500: 11,
    550: 11,
    600: 12.5
}


dict_mantel_thickness_std_KMR = {
        10:2.62,
        20:3,
        25:3,
        32:3,
        40:3,
        50:3,
        65:3,
        80:3,
        100:3.2,
        125:3.4,
        150:3.6,
        175:3.9,
        200:4.1,
        225:4.5,
        250:4.8,
        300:5.2,
        350:5.6,
        400:6,
        450:6.6,
        500:6.9,
        550:7.2,
        600:7.9  
}

dict_mantel_thickness_1x_KMR = {
        10:2.62,
        20:3,
        25:3,
        32:3,
        40:3,
        50:3,
        65:3,
        80:3,
        100:3.4,
        125:3.6,
        150:3.9,
        175:4.2,
        200:4.5,
        225:4.8,
        250:5.2,
        300:5.6,
        350:6,
        400:6.6,
        450:6.9,
        500:7.2,
        550:7.9,
        600:8.7  
}

# Wall roughness (mm to mm)
dict_k_mm_KMR = {
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
}

dict_k_mm_PEX = {
    10: 0.007,
    20: 0.007,
    25: 0.007,
    32: 0.007,
    40: 0.007,
    50: 0.007,
    65: 0.007,
    80: 0.007,
    100: 0.007,
    125: 0.007,
    150: 0.007,
    175: 0.007,
    200: 0.007,
    225: 0.007,
    250: 0.007,
    300: 0.007,
    350: 0.007,
    400: 0.007,
    450: 0.007,
    500: 0.007,
    550: 0.007,
    600: 0.007
}

# Insulation thickness
dict_insulation_thickness_std_KMR={
        10:25.935,
        20:28.55,
        25:25.15,
        32:30.8,
        40:27.85,
        50:29.35,
        65:28.95,
        80:32.55,
        100:39.65,
        125:39.25,
        150:37.25,
        175:39.25,
        200:43.85,
        225:50.75,
        250:57.75,
        300:57.85,
        350:66.6,
        400:70.8,
        450:79.9,
        500:94.1,
        550:68.4,
        600:87.1
        }

dict_insulation_thickness_1x_KMR={
    10:31.769,
    20:38.55,
    25:35.15,
    32:38.3,
    40:35.35,
    50:36.85,
    65:38.95,
    80:42.55,
    100:51.95,
    125:51.55,
    150:51.95,
    175:56.55,
    200:63.45,
    225:72.95,
    250:83.3,
    300:82.45,
    350:96.2,
    400:105.2,
    450:119.6,
    500:138.8,
    550:113.4,
    600:137.1
    }



# Doppelrohre: Lichter Rohrabstand
dict_pipe_distance_doubleKMR = {
    20:19,
    25:19,
    32:19,
    40:19,
    50:20,
    65:20,
    80:25,
    100:25,
    125:30,
    150:40,
    175:45,
    200:45
}

# Doppelrohre: Außendurchmesser
dict_pipe_diam_out_std_doubleKMR = {
    20:125,
    25:140,
    32:160,
    40:160,
    50:200,
    65:225,
    80:250,
    100:315,
    125:400,
    150:450,
    175:500,
    200:560   

}

dict_mantel_thickness_doppelrohr_std_doubleKMR = {
    20:3,
    25:3,
    32:3,
    40:3,
    50:3.2,
    65:3.4,
    80:3.6,
    100:4.1,
    125:4.8,
    150:5.2,
    175:5.5,
    200:6
}


# %%
''' PE SDR11 networks (5GDHC networks) '''

nw_keys = [25,32,40,50,63,75,90,110,125,140,160,180,200,225,250,280,315,355,400,450,500,560,630,710,800]

# Dictionary nominal diameters to inner diameters:
nominal_to_inner_diameters_5GDHC = {
    25: 20.4,
    32: 26.2,
    40: 32.8,
    50: 41,
    63: 51.6,
    75: 61.4,
    90: 73.6,
    110: 90,
    125: 102.2,
    140: 114.6,
    160: 131,
    180: 147.2,
    200: 163.6,
    225: 184,
    250: 204.6,
    280: 229,
    315: 257.8,
    355: 290.4,
    400: 327.2,
    450: 368.2,
    500: 409,
    560: 458.2,
    630: 515.4,
    710: 581,
    800: 654.6
}



# Wall thickness (inside insulation, steel etc.)
dict_Wanddicke_5GDHC = {
    x:x/11 for x in nw_keys
}

# Wall roughness
dict_k_mm_5GDHC = {
    x:0.007 for x in nw_keys
}
