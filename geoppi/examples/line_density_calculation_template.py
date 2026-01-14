# %%

import os
import geoppi
import geopandas as gp
from pathlib import Path

cs = 'EPSG:25832'

## Load example data
flp = os.getcwd() / Path(r'data/exampleNetwork2/')
lines = gp.read_file(flp / Path(r'streets.gpkg')).to_crs(cs)

buildings = gp.read_file(flp / Path(r'buildings.gpkg')).to_crs(cs)

hex = None#gp.read_file(flp / Path(r'hex_div.gpkg')).to_crs(cs)

# Define whether a randomized sampling in each hexagon with aim connection ratio shall take place
rand_sampling = True

# Number of randomized samples for selection of buildings to attain aimAG (sample with median heat demand is chosen)
nSamples = 1

# Define minimum number of buildings which shall be found within single hexagon below which the aimAG is ignored and ALL buildings are considered for calculation
nBuildings_min = 1

# Define aim of connection ratio within each subdivision of the regarded area in hex_div
target_CR = 1

# Define hexagon-specific connection ratio as attribute name of its value in hex
attr_hex_CR = 'connectionRatio'


# Define attribute from buildings for calculation of line density
target_attr = 'demand_use_th'

# Define additional attributes from the buildings layer which shall be summed on closest line objects
summed_target_attrs = ['demand_use_th_2045_san']

# Define maximum distance between buildings and line objects to include them into calculation of line density
# -> Example uses function-like distance calculation
def calc_distance_from_heat_demand(
    heat_demand, 
    min_line_density = 1000, # At least line density of 1 MWh/m at house connection line to assign closest line to polygon
    minLimit = 10,
    maxLimit = 100):

    return(min(maxLimit, max(minLimit, heat_demand/min_line_density)) )



# %%

### Start calculation
lines_out, polys_out = geoppi.sum_attributes_on_lines(
    polygons = buildings.head(100),
    lines = lines,
    spatial_distribution = hex,
    rand_sampling = rand_sampling,
    nSamples = nSamples,
    nPolygons_min = nBuildings_min,
    target_connection_ratio = target_CR,
    spatial_connection_ratio = attr_hex_CR,
    target_attr = 'demand_use_th',
    agg_func = 'median',
    additional_attr = summed_target_attrs,
    func_max_distance = None #calc_distance_from_heat_demand
)

lines_out[f'ld_{target_attr}_MWh_per_m'] = lines_out[f'summed_{target_attr}'] / 1e03 / lines_out.geometry.length


# lines_out.to_file(flp / Path(r'results/streets_line_density.gpkg'), driver = 'GPKG')
# polys_out.to_file(flp / Path(r'results/polygons_line_density.gpkg'), driver = 'GPKG')