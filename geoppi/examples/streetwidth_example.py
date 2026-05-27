import os
import sys
import time
import geopandas as gp
import numpy as np  

sys.path.append(os.path.join(os.path.dirname(__file__), "..", ".."))
from auxFunctions import measure_street_width, detect_lines_in_narrow_passages

_input = os.path.join(os.path.dirname(__file__), "..", "data", "exampleStreetwidth", "input")
_output = os.path.join(os.path.dirname(__file__), "..", "data", "exampleStreetwidth", "output")
axis_lines = gp.read_file(os.path.join(_input, "strassen_buseck.shp"))
street_parcels = gp.read_file(os.path.join(_input, "flurstuecke_strassenverkehr_buseck.gpkg"))
buildings = gp.read_file(os.path.join(_input, "flursteucke_buseck_ohne_straßenverkehr.gpkg"))
buildings = buildings[buildings.geometry.length >= 1]


print("=" * 60)
print("Test 1: measure_street_width")
print("=" * 60)

t0 = time.time()

profile_gdf, lines_gdf = measure_street_width(
    axis_lines=axis_lines,
    street_parcels=street_parcels,
    step_size=5.0,
    max_range=30.0,
    max_valid_width=25.0, # gesamte breite
    fan_angles=np.linspace(-15, 15, 7).tolist(),
    overlap_threshold=2,
    min_boundary_angle=70.0,
)

t1 = time.time()

print(f"Dauer: {t1 - t0:.2f}s")
print(f"Profile-Punkte: {len(profile_gdf)}")
print(f"Messlinien:     {len(lines_gdf)}")
profile_gdf.to_file(os.path.join(_output, "width_profile.gpkg"))
lines_gdf.to_file(os.path.join(_output, "width_messlinien.gpkg"))
print()


print("=" * 60)
print("Test 2: detect_lines_in_narrow_passages")
print("=" * 60)

t0 = time.time()

result_lines, shortest_lines = detect_lines_in_narrow_passages(
    lines=axis_lines,
    polygons=buildings,
    merge_touching_polygons=True,
    threshDistance=10,
    distPointsCircumference=5,
    nNeighbours=10,
    col="narrowPassage_m",
)

t1 = time.time()

print(f"Dauer: {t1 - t0:.2f}s")
print(f"Linien mit Engstelle: {result_lines['narrowPassage_m'].notna().sum()} / {len(result_lines)}")
print(f"Kürzeste Verbindungen: {len(shortest_lines)}")
result_lines = result_lines.drop(columns=['fid'], errors='ignore')
result_lines.to_file(os.path.join(_output, "narrow_engstellen.gpkg"))
shortest_lines.to_file(os.path.join(_output, "narrow_verbindungen.gpkg"))