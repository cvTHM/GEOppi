# %%
import geopandas as gp
import numpy as np
from pathlib import Path
import geoppi
from geoppi.create_network_topology import createUniqueJunctions


DATA_DIR = Path(__file__).resolve().parent / "geoppi" / "examples" / "data" / "exampleNetwork"


# %% Load raw example data from geoppi/examples/data/exampleNetwork
lines = gp.read_file(DATA_DIR / "streetsRaw.gpkg").to_crs("EPSG:25832")
buildings = gp.read_file(DATA_DIR / "buildings.gpkg").to_crs("EPSG:25832")
producers = gp.read_file(DATA_DIR / "startpoints_producers2.gpkg").to_crs("EPSG:25832")

# Split input lines at internal intersections so the graph topology reflects crossings.
junctions, lines = createUniqueJunctions(lines.explode(index_parts=False).reset_index(drop=True), splitLines=True)

# Use only the first producer for this simplified example.
# producers = producers.head(1).reset_index(drop=True)

lines["unique_ID_lines"] = np.arange(len(lines), dtype=int)
lines["length"] = lines.geometry.length
buildings["unique_ID_polys"] = np.arange(len(buildings), dtype=int)


# %% Example 1: count downstream buildings for a single producer using a building attribute value of 1

# Use only the first producer for this variant.
producer_single = producers.head(1).reset_index(drop=True)

# Create a building attribute that equals 1 for every building so counts can be summed.
buildings['n_buildings'] = 1

lines_downstream_single = geoppi.sum_attrs_to_closest_supplier(
    lines=lines,
    buildings=buildings,
    producers=producer_single,
    building_attrs=['n_buildings', "demand_use_th"],
    output_attr=["sum_n_buildings_downstream", "sum_demand_use_th_downstream"],
    weight='length',
    line_id='unique_ID_lines',
    building_id='unique_ID_polys',
    max_distance=100.0,
)

print('Downstream building counts for single producer:')

print(lines_downstream_single[['unique_ID_lines', 'sum_n_buildings_downstream', 'sum_demand_use_th_downstream']].head(10))

lines_downstream_single.to_file(DATA_DIR / 'lines_downstream_buildings_single.gpkg')

# %% Example 2: count downstream buildings for both producers using the same building attribute
producer_attributes_both = {idx: 1 for idx in producers.index}
buildings['n_buildings'] = 1

lines_downstream_both = geoppi.sum_attrs_to_closest_supplier(
    lines=lines,
    buildings=buildings,
    producers=producers,
    building_attrs=['n_buildings', "demand_use_th"],
    output_attr=["sum_n_buildings_downstream_both", "sum_demand_use_th_downstream_both"],
    weight='length',
    line_id='unique_ID_lines',
    building_id='unique_ID_polys',
    max_distance=100.0,
)

print('Downstream building counts for both producers:')
print(lines_downstream_both[['unique_ID_lines', 'sum_n_buildings_downstream_both', 'sum_demand_use_th_downstream_both']].head(10))

lines_downstream_both.to_file(DATA_DIR / 'lines_downstream_buildings_both.gpkg')


# %% Example 3: Gas network Kleinlinden

DATA_DIR = Path(r"Z:\13_Skizzen_Anträge_Projekte\17_VNEP_Gasnetzrückbauplanung_Entwicklung\example_gas_network")

# Load raw example data from geoppi/examples/data/exampleNetwork
lines = gp.read_file(DATA_DIR / "gas_network.gpkg").to_crs("EPSG:25832")
lines.drop(columns = [n for n in lines.columns if n not in ("DRUCKSTUFE", "NENNWEITE1", "LEITUNGSFU", "geometry")], inplace = True)

lines = lines[lines["DRUCKSTUFE"] == "Niederdruck"].reset_index(drop = True)

buildings = gp.read_file(DATA_DIR / "buildings_final.gpkg").to_crs("EPSG:25832")
buildings = buildings[buildings["demand_use_th"] > 1000].reset_index(drop = True)

producers = gp.read_file(DATA_DIR / "substations.gpkg").to_crs("EPSG:25832")

# Split input lines at internal intersections so the graph topology reflects crossings.
_, lines = createUniqueJunctions(lines.explode(index_parts=False).reset_index(drop=True), splitLines=True)

lines = lines[lines["LEITUNGSFU"] == "Versorgungsleitung"].reset_index(drop = True)
lines["unique_ID_lines"] = np.arange(len(lines), dtype=int)
lines["length"] = lines.geometry.length

buildings["unique_ID_polys"] = np.arange(len(buildings), dtype=int)

buildings['n_buildings'] = 1

lines_downstream_both = geoppi.sum_attrs_to_closest_supplier(
    lines=lines,
    buildings=buildings,
    producers=producers,
    building_attrs=['n_buildings', "demand_use_th"],
    output_attr=["sum_n_buildings_downstream_both", "sum_demand_use_th_downstream_both"],
    weight='length',
    line_id='unique_ID_lines',
    building_id='unique_ID_polys',
    max_distance=100.0,
)

print('Downstream building counts for both producers:')
print(lines_downstream_both[['unique_ID_lines', 'sum_n_buildings_downstream_both', 'sum_demand_use_th_downstream_both']].head(10))

lines_downstream_both.to_file(DATA_DIR / 'lines_downstream_buildings_both.gpkg')

# %% Example 4: Use nominal width as additional weight in graph

producers = producers.head(1).reset_index(drop=True)

buildings['n_buildings'] = 1
lines["nominal_width"] = (
    lines["NENNWEITE1"]
    .astype(str)
    .str.extract(r"(\d+)", expand=False)
    .astype("Int64")
)

lines["weight_total"] = lines["length"] * (1/lines["nominal_width"]/1000)**5

lines_downstream_both = geoppi.sum_attrs_to_closest_supplier(
    lines=lines,
    buildings=buildings,
    producers=producers,
    building_attrs=['n_buildings', "demand_use_th"],
    output_attr=["sum_n_buildings_downstream_single", "sum_demand_use_th_downstream_single"],
    weight='weight_total',
    line_id='unique_ID_lines',
    building_id='unique_ID_polys',
    max_distance=100.0,
)

print('Downstream building counts for both producers:')
print(lines_downstream_both[['unique_ID_lines', 'sum_n_buildings_downstream_single', 'sum_demand_use_th_downstream_single']].head(10))

lines_downstream_both.to_file(DATA_DIR / 'lines_downstream_buildings_single_NW.gpkg')