# %%
import geopandas as gp
import numpy as np
from pathlib import Path
import geoppi
from geoppi.create_network_topology import createUniqueJunctions
from geoppi.line_density_calculation import closest_lines_to_polygons


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

# %% Example 1: count downstream buildings by matching building boundaries to nearest line segments

lines_downstream = geoppi.count_downstream_buildings(
    lines=lines,
    buildings=buildings,
    producers=producers,
    output_attr="n_buildings_downstream",
    weight="length",
    line_id="unique_ID_lines",
    building_id="unique_ID_polys",
    max_distance=100.0,
)

print("Downstream building counts (nearest boundary matching):")
print(lines_downstream[["unique_ID_lines", "n_buildings_downstream"]].head(10))

lines_downstream.to_file(DATA_DIR / "lines_downstream_buildings.gpkg")

# %% Example 2: precompute a matching dictionary from building IDs to line IDs
buildings_matched = closest_lines_to_polygons(
    polygons=buildings[["unique_ID_polys", "geometry"]].copy(),
    lines=lines[["unique_ID_lines", "geometry"]].copy(),
    maxDistances=100.0,
)

buildings_matched = buildings_matched.rename(columns={"nearestline": "nearest_line_index"})
dict_building2line = {
    int(row["unique_ID_polys"]): int(lines.loc[row["nearest_line_index"], "unique_ID_lines"])
    for _, row in buildings_matched.dropna(subset=["nearest_line_index"]).iterrows()
}

lines_downstream_dict = geoppi.count_downstream_buildings(
    lines=lines,
    buildings=buildings,
    producers=producers,
    output_attr="n_buildings_downstream_dict",
    weight="length",
    line_id="unique_ID_lines",
    building_id="unique_ID_polys",
    dict_polyID_lineID=dict_building2line,
    max_distance=100.0,
)

print("Downstream building counts (precomputed ID mapping):")
print(lines_downstream_dict[["unique_ID_lines", "n_buildings_downstream_dict"]].head(10))


