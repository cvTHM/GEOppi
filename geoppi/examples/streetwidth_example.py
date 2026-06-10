# %%

import time
from pathlib import Path
import geopandas as gp
import numpy as np

from geoppi import measure_street_width, detect_lines_in_narrow_passages

_input = Path(r'data/exampleStreetwidth/input')
_output = Path(r'data/exampleStreetwidth/output')

axis_lines = gp.read_file(_input / 'streets.gpkg')
street_parcels = gp.read_file(_input / 'parcels_traffic.gpkg')
parcels = gp.read_file(_input / 'parcels.gpkg')


print("=" * 60)
print("Test 1: measure_street_width")
print("=" * 60)

t0 = time.time()

axis_lines_gdf, profile_gdf, lines_gdf = measure_street_width(
    axis_lines=axis_lines,
    street_parcels=street_parcels,
    step_size=5.0,
    max_range=30.0,
    max_valid_width=25.0, # Total width
    fan_angles=np.linspace(-15, 15, 7).tolist(), # Ex.: 7 measurement line with angle of -15° to +15° rel. to street axis
    overlap_threshold=2,
    min_boundary_angle=70.0,
)


t1 = time.time()

print(f"Duration: {t1 - t0:.2f}s")
print(f"Profile points: {len(profile_gdf)}")
print(f"Measurement lines:     {len(lines_gdf)}")
axis_lines_gdf.to_file(_output / 'width_street_lines.gpkg')
profile_gdf.to_file(_output / 'width_profile.gpkg', driver = "GPKG")
lines_gdf.to_file(_output / 'width_messlinien.gpkg', driver = "GPKG")
print()


print("=" * 60)
print("Test 2: detect_lines_in_narrow_passages")
print("=" * 60)

t0 = time.time()

result_lines, shortest_lines = detect_lines_in_narrow_passages(
    lines=axis_lines,
    polygons=parcels,
    merge_touching_polygons=True,
    threshDistance=10,
    distPointsCircumference=5,
    nNeighbours=10,
    col="narrowPassage_m",
)

t1 = time.time()

print(f"Duration: {t1 - t0:.2f}s")
print(f"Lines with narrow passage: {result_lines['narrowPassage_m'].notna().sum()} / {len(result_lines)}")
print(f"Shortest connection: {len(shortest_lines)}")
result_lines = result_lines.drop(columns=['fid'], errors='ignore')
result_lines.to_file(_output / 'narrow_engstellen.gpkg')
shortest_lines.to_file(_output / 'narrow_verbindungen.gpkg')