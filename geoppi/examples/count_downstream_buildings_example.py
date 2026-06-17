# %%
import geopandas as gp
import numpy as np
from pathlib import Path
import geoppi
from geoppi.create_network_topology import createUniqueJunctions


DATA_DIR = Path(__file__).resolve().parent / "data" / "exampleNetwork"


# %% Load raw example data from geoppi/examples/data/exampleNetwork
lines = gp.read_file(DATA_DIR / "streetsRaw.gpkg").to_crs("EPSG:25832")
buildings = gp.read_file(DATA_DIR / "buildings.gpkg").to_crs("EPSG:25832")
producers = gp.read_file(DATA_DIR / "startpoints_producers2.gpkg").to_crs("EPSG:25832")

# Split input lines at internal intersections so the graph topology reflects crossings.
junctions, lines = createUniqueJunctions(lines.explode(index_parts=False).reset_index(drop=True), splitLines=True)

# Use only the first producer for this simplified example.
producers = producers.head(1).reset_index(drop=True)

lines["unique_ID_lines"] = np.arange(len(lines), dtype=int)
buildings["unique_ID_polys"] = np.arange(len(buildings), dtype=int)

# %% Example 1: count downstream buildings for a single producer using an attribute sum per producer

# Use only the first producer for this variant.
producer_single = producers.head(1).reset_index(drop=True)
producer_attributes_single = {producer_single.index[0]: 1}

lines_downstream_single = geoppi.count_downstream_buildings(
    lines=lines,
    buildings=buildings,
    producers=producer_single,
    output_attr="n_buildings_downstream_single",
    weight="length",
    line_id="unique_ID_lines",
    building_id="unique_ID_polys",
    dict_producer_attrs=producer_attributes_single,
    max_distance=100.0,
)

print("Downstream building counts for single producer:")
print(lines_downstream_single[["unique_ID_lines", "n_buildings_downstream_single"]].head(10))

lines_downstream_single.to_file(DATA_DIR / "lines_downstream_buildings_single.gpkg")

# %% Example 2: count downstream buildings for both producers with a custom producer attribute mapping
producer_attributes_both = {idx: 1 for idx in producers.index}

lines_downstream_both = geoppi.count_downstream_buildings(
    lines=lines,
    buildings=buildings,
    producers=producers,
    output_attr="n_buildings_downstream_both",
    weight="length",
    line_id="unique_ID_lines",
    building_id="unique_ID_polys",
    dict_producer_attrs=producer_attributes_both,
    max_distance=100.0,
)

print("Downstream building counts for both producers:")
print(lines_downstream_both[["unique_ID_lines", "n_buildings_downstream_both"]].head(10))

lines_downstream_both.to_file(DATA_DIR / "lines_downstream_buildings_both.gpkg")


