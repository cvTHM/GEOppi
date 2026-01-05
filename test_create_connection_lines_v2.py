# %%
import geopandas as gp
from shapely.geometry import LineString, Point, MultiPoint
from shapely.ops import nearest_points, split
import numpy as np
from tqdm import tqdm
import shapely
import shapely.plotting

from shapely.ops import unary_union
import matplotlib.pyplot as plt

from geoppi.auxFunctions import assign_attr_by_max_intersection_area

def create_shortest_connections_optimized(polygons_gdf, lines_gdf, buffer_dist=0.1):
    # Schritt 1: Baue einen räumlichen Index für die Linien
    lines_gdf = lines_gdf.copy()
    lines_gdf['geometry'] = lines_gdf['geometry'].apply(lambda geom: geom if geom.is_valid else geom.buffer(0))
    sindex = lines_gdf.sindex

    connection_lines = []

    for poly in tqdm(polygons_gdf.geometry, total = len(polygons_gdf)):
        # Optional: Punkte entlang des Polygon-Randes (z.B. alle 1 m)
        edge = poly.exterior
        num_points = min(8, int(edge.length))

        possible_connection_lines_poly = []

        for i in range(num_points):
            pt = edge.interpolate(i / num_points, normalized=True)
            # 0.1m ins Polygon hinein verschieben
            direction = poly.centroid.coords[0]
            vec = np.array(direction) - np.array(pt.coords[0])
            vec = vec / np.linalg.norm(vec) * buffer_dist
            inward_pt = Point(np.array(pt.coords[0]) + vec)
            # Suche nur Linien, die im Umkreis von 10m liegen (Bounding Box)
            possible_matches_index = list(sindex.intersection(inward_pt.buffer(100).bounds))
            if not possible_matches_index:
                continue
            possible_lines = lines_gdf.iloc[possible_matches_index]
            # Finde die Linie mit der minimalen Distanz
            dists = possible_lines.distance(inward_pt)
            min_idx = dists.idxmin()
            nearest_line = possible_lines.loc[min_idx].geometry
            # Finde den nächsten Punkt auf der Linie
            nearest_on_line = nearest_points(inward_pt, nearest_line)[1]
            connection = LineString([inward_pt, nearest_on_line])

            possible_connection_lines_poly.append(connection)

        if len(possible_connection_lines_poly) > 0:
            allLengths = [l.length for l in possible_connection_lines_poly]
            connection_lines.append(possible_connection_lines_poly[np.argmin(allLengths)])        


    connection_gdf = gp.GeoDataFrame(geometry=connection_lines, crs=polygons_gdf.crs)

    def round_point(pt, precision=0.01):
        """Rundet die Koordinaten eines Punktes auf die angegebene Präzision."""
        x = round(pt.x / precision) * precision
        y = round(pt.y / precision) * precision
        return Point(x, y)

    def extend_line(line, length):
        coords = list(line.coords)
        if len(coords) < 2:
            return line
        # Richtung am Anfang
        dx1 = coords[1][0] - coords[0][0]
        dy1 = coords[1][1] - coords[0][1]
        norm1 = np.hypot(dx1, dy1)
        x0 = coords[0][0] - dx1 / norm1 * length
        y0 = coords[0][1] - dy1 / norm1 * length
        # Richtung am Ende
        dx2 = coords[-1][0] - coords[-2][0]
        dy2 = coords[-1][1] - coords[-2][1]
        norm2 = np.hypot(dx2, dy2)
        x1 = coords[-1][0] + dx2 / norm2 * length
        y1 = coords[-1][1] + dy2 / norm2 * length
        new_coords = [(x0, y0)] + coords[1:-1] + [(x1, y1)]
        return LineString(new_coords)

    def insert_points_in_linestring(line, points, precision=0.01, tolerance=1e-8):
        """Fügt Punkte als zusätzliche Stützpunkte in einen LineString ein, mit Rundung."""
        if not points:
            return line
        # Alle Stützpunkte der Linie als Liste
        coords = list(line.coords)
        # Für jeden Punkt: berechne den Abstand entlang der Linie
        distances = [line.project(pt) for pt in points]
        # Kombiniere bestehende Vertices (mit deren Position entlang der Linie)
        vertex_distances = [line.project(Point(c)) for c in coords]
        # Füge alle Punkte (Vertices + neue Punkte) zusammen
        all_distances = vertex_distances + distances
        all_points = [round_point(line.interpolate(d), precision) for d in all_distances]
        # Entferne Duplikate (innerhalb der Toleranz)
        unique = []
        used = []
        for d, pt in sorted(zip(all_distances, all_points)):
            if not any(pt.equals_exact(u, tolerance) for u in unique):
                unique.append(pt)
                used.append(d)
        return LineString(unique)

    def fast_split_lines_by_points(lines_gdf, connection_lines_gdf, extend_connection=1.0, precision=0.01, tolerance=1e-8):
        # 1. Verbindungslinien um 1m verlängern
        extended_conns = connection_lines_gdf.geometry.apply(
            lambda l: extend_line(l, extend_connection)
        )
        # 2. Räumlicher Index für die Verbindungslinien
        splitter_sindex = gp.GeoSeries(extended_conns, crs=connection_lines_gdf.crs).sindex

        split_lines = []
        for line in lines_gdf.geometry:
            # Nur relevante Verbindungslinien
            possible_idx = list(splitter_sindex.intersection(line.bounds))
            if not possible_idx:
                split_lines.append(line)
                continue
            splitters = [extended_conns.iloc[i] for i in possible_idx]

            # Alle Schnittpunkte sammeln
            split_points = []
            for splitter in splitters:
                inter = line.intersection(splitter)
                if inter.is_empty:
                    continue
                if inter.geom_type == 'Point':
                    split_points.append(inter)
                elif inter.geom_type == 'MultiPoint':
                    split_points.extend(list(inter.geoms))
                elif inter.geom_type == 'LineString':
                    split_points.extend([Point(c) for c in inter.coords])
                elif inter.geom_type == 'GeometryCollection':
                    for geom in inter.geoms:
                        if geom.geom_type == 'Point':
                            split_points.append(geom)
                        elif geom.geom_type == 'LineString':
                            split_points.extend([Point(c) for c in geom.coords])

            # Schnittpunkte exakt auf die Linie projizieren (snappen & runden)
            snapped_points = []
            used_proj = set()
            for pt in split_points:
                d = line.project(pt)
                if precision < d < line.length - precision and d not in used_proj:
                    snapped = round_point(line.interpolate(d), precision)
                    snapped_points.append(snapped)
                    used_proj.add(d)

            # 3. Kreuzungspunkte als Stützpunkte in die Linie einfügen
            if snapped_points:
                new_line = insert_points_in_linestring(line, snapped_points, precision, tolerance)
                # 4. Splitte an genau diesen Punkten (MultiPoint)
                splitter_geom = MultiPoint(snapped_points)
                result = split(new_line, splitter_geom)
                split_lines.extend([seg for seg in result.geoms if seg.length > precision])
            else:
                split_lines.append(line)
        return gp.GeoDataFrame(geometry=split_lines, crs=lines_gdf.crs)

    split_lines_gdf = fast_split_lines_by_points(lines_gdf = lines_gdf, connection_lines_gdf = connection_gdf, extend_connection = 0.1, precision = 0.01, tolerance = 1e-08)

    return connection_gdf, split_lines_gdf

def create_shortest_connections_optimized_v2(polygons_gdf, lines_gdf, buffer_dist=0.1):
    # Schritt 1: Baue einen räumlichen Index für die Linien
    lines_gdf = lines_gdf.copy()
    lines_gdf['geometry'] = lines_gdf['geometry'].apply(lambda geom: geom if geom.is_valid else geom.buffer(0))
    sindex = lines_gdf.sindex

    connection_lines = []

    for poly in tqdm(polygons_gdf.geometry, total = len(polygons_gdf)):
        # Optional: Punkte entlang des Polygon-Randes (z.B. alle 1 m)
        edge = poly.exterior
        num_points = min(8, int(edge.length))

        possible_connection_lines_poly = []

        for i in range(num_points):
            pt = edge.interpolate(i / num_points, normalized=True)
            # 0.1m ins Polygon hinein verschieben
            direction = poly.centroid.coords[0]
            vec = np.array(direction) - np.array(pt.coords[0])
            vec = vec / np.linalg.norm(vec) * buffer_dist
            inward_pt = Point(np.array(pt.coords[0]) + vec)
            # Suche nur Linien, die im Umkreis von 100m liegen (Bounding Box)
            possible_matches_index = list(sindex.intersection(inward_pt.buffer(100).bounds))
            if not possible_matches_index:
                continue
            possible_lines = lines_gdf.iloc[possible_matches_index]
            # Finde die Linie mit der minimalen Distanz
            dists = possible_lines.distance(inward_pt)
            min_idx = dists.idxmin()
            nearest_line = possible_lines.loc[min_idx].geometry
            # Finde den nächsten Punkt auf der Linie
            nearest_on_line = nearest_points(inward_pt, nearest_line)[1]
            connection = LineString([inward_pt, nearest_on_line])

            possible_connection_lines_poly.append(connection)

        if len(possible_connection_lines_poly) > 0:
            allLengths = [l.length for l in possible_connection_lines_poly]
            connection_lines.append(possible_connection_lines_poly[np.argmin(allLengths)])        


    connection_gdf = gp.GeoDataFrame(geometry=connection_lines, crs=polygons_gdf.crs)

    def round_point(pt, precision=0.01):
        x = round(pt.x / precision) * precision
        y = round(pt.y / precision) * precision
        return Point(x, y)

    def extend_line(line, length):
        coords = list(line.coords)
        if len(coords) < 2:
            return line
        # Richtung am Anfang
        dx1 = coords[1][0] - coords[0][0]
        dy1 = coords[1][1] - coords[0][1]
        norm1 = np.hypot(dx1, dy1)
        x0 = coords[0][0] - dx1 / norm1 * length
        y0 = coords[0][1] - dy1 / norm1 * length
        # Richtung am Ende
        dx2 = coords[-1][0] - coords[-2][0]
        dy2 = coords[-1][1] - coords[-2][1]
        norm2 = np.hypot(dx2, dy2)
        x1 = coords[-1][0] + dx2 / norm2 * length
        y1 = coords[-1][1] + dy2 / norm2 * length
        new_coords = [(x0, y0)] + coords[1:-1] + [(x1, y1)]
        return LineString(new_coords)

    def insert_points_in_linestring(line, points, precision=0.01, tolerance=1e-8):
        if not points:
            return line
        coords = list(line.coords)
        distances = [line.project(pt) for pt in points]
        vertex_distances = [line.project(Point(c)) for c in coords]
        all_distances = vertex_distances + distances
        all_points = [round_point(line.interpolate(d), precision) for d in all_distances]
        unique = []
        used = []
        for d, pt in sorted(zip(all_distances, all_points)):
            if not any(pt.equals_exact(u, tolerance) for u in unique):
                unique.append(pt)
                used.append(d)
        return LineString(unique)

    def generate_precise_connection_lines(polygons_gdf, lines_gdf, buffer_dist=0.1, precision=0.01, search_radius=100):
        """
        Schnelle Erzeugung von Verbindungslinien von Polygonen zu Liniennetz.
        Die Verbindungslinien enden exakt auf gerundeten Schnittpunkten (precision).
        """
        # 1. Räumlicher Index für Liniennetz
        lines_sindex = lines_gdf.sindex

        connection_lines = []
        crossing_points = []

        for poly in tqdm(polygons_gdf.geometry, total = len(polygons_gdf)):
            min_dist = np.inf
            best_start = None
            best_target = None

            edge = poly.exterior
            num_points = min(int(edge.length), 10)
            for i in range(num_points):
                pt = edge.interpolate(i / num_points, normalized=True)
                centroid = poly.centroid
                direction = np.array(centroid.coords[0]) - np.array(pt.coords[0])
                direction = direction / np.linalg.norm(direction) * buffer_dist
                inner_pt = Point(np.array(pt.coords[0]) + direction)

                # Nur Linien im Nahbereich betrachten (schnell!)
                bbox = inner_pt.buffer(search_radius).bounds
                possible_idx = list(lines_sindex.intersection(bbox))
                if not possible_idx:
                    continue
                possible_lines = lines_gdf.iloc[possible_idx]

                # Finde die nächste Linie und den exakten Zielpunkt
                dists = possible_lines.distance(inner_pt)
                min_idx = dists.idxmin()
                nearest_line = possible_lines.loc[min_idx].geometry
                # Exakter Zielpunkt auf der Linie
                proj = nearest_line.project(inner_pt)
                target = nearest_line.interpolate(proj)
                dist = inner_pt.distance(target)
                if dist < min_dist:
                    min_dist = dist
                    best_start = inner_pt
                    best_target = target

            if best_start is not None and best_target is not None:
                rounded_target = round_point(best_target, precision)
                connection_lines.append(LineString([best_start, rounded_target]))
                crossing_points.append(rounded_target)

        connection_gdf = gp.GeoDataFrame(geometry=connection_lines, crs=polygons_gdf.crs)
        crossing_points_gdf = gp.GeoDataFrame(geometry=crossing_points, crs=polygons_gdf.crs)
        return connection_gdf, crossing_points_gdf

    def fast_split_lines_by_points(lines_gdf, connection_lines_gdf, precision=0.01, tolerance=1e-8, extend_connection=1.0):
        """
        Schneidet lines_gdf an den Kreuzungspunkten mit connection_lines_gdf.
        Die Kreuzungspunkte werden auf 0.01m gerundet und als Stützpunkte eingefügt.
        """
        # 1. Verbindungslinien um 1m verlängern (für robustes Kreuzungsfinden)
        extended_conns = connection_lines_gdf.geometry.apply(
            lambda l: extend_line(l, extend_connection)
        )
        splitter_sindex = gp.GeoSeries(extended_conns, crs=connection_lines_gdf.crs).sindex

        split_lines = []
        for line in tqdm(lines_gdf.geometry, total = len(lines_gdf)):
            possible_idx = list(splitter_sindex.intersection(line.bounds))
            if not possible_idx:
                split_lines.append(line)
                continue
            splitters = [extended_conns.iloc[i] for i in possible_idx]
            split_points = []
            for splitter in splitters:
                inter = line.intersection(splitter)
                if inter.is_empty:
                    continue
                if inter.geom_type == 'Point':
                    split_points.append(inter)
                elif inter.geom_type == 'MultiPoint':
                    split_points.extend(list(inter.geoms))
                elif inter.geom_type == 'LineString':
                    split_points.extend([Point(c) for c in inter.coords])
                elif inter.geom_type == 'GeometryCollection':
                    for geom in inter.geoms:
                        if geom.geom_type == 'Point':
                            split_points.append(geom)
                        elif geom.geom_type == 'LineString':
                            split_points.extend([Point(c) for c in geom.coords])
            # Kreuzungspunkte exakt auf die Linie projizieren (snappen & runden)
            snapped_points = []
            used_proj = set()
            for pt in split_points:
                d = line.project(pt)
                if precision < d < line.length - precision and d not in used_proj:
                    snapped = round_point(line.interpolate(d), precision)
                    snapped_points.append(snapped)
                    used_proj.add(d)
            # Kreuzungspunkte als Stützpunkte in die Linie einfügen
            if snapped_points:
                new_line = insert_points_in_linestring(line, snapped_points, precision, tolerance)
                splitter_geom = MultiPoint(snapped_points)
                result = split(new_line, splitter_geom)
                split_lines.extend([seg for seg in result.geoms if seg.length > precision])
            else:
                split_lines.append(line)
        return gp.GeoDataFrame(geometry=split_lines, crs=lines_gdf.crs)
    
    # Aufrufe
    connection_lines_gdf, crossing_points_gdf = generate_precise_connection_lines(polygons_gdf, lines_gdf)

    split_lines_gdf = fast_split_lines_by_points(lines_gdf = lines_gdf, connection_lines_gdf = connection_lines_gdf, extend_connection = 0.1, precision = 0.01, tolerance = 1e-08)

    return connection_gdf, split_lines_gdf


def extend_line(line, offset):
    coords = list(line.coords)
    # Get the direction vector at the start of the LineString
    start_vec = np.array(coords[1]) - np.array(coords[0])
    start_dir = start_vec / np.linalg.norm(start_vec)
    # Calculate the new start point by moving backwards along the direction
    new_start = np.array(coords[0]) - offset * start_dir

    # Get the direction vector at the end of the LineString
    end_vec = np.array(coords[-1]) - np.array(coords[-2])
    end_dir = end_vec / np.linalg.norm(end_vec)
    # Calculate the new end point by moving forwards along the direction
    new_end = np.array(coords[-1]) + offset * end_dir

    # Build a new coordinate list with the extended endpoints
    new_coords = [tuple(new_start)] + coords[1:-1] + [tuple(new_end)]
    return LineString(new_coords)

from geoppi.auxFunctions import (nnearest, create_circumferential_points, create_point_splitter_npoints, split_lines_at_length, closest_objects_to_points)

cs = 'EPSG:25832'



lines = gp.read_file(r'geoppi/examples/data/lineSeparation/lines_separation_raw.gpkg').to_crs(cs)
polys = gp.read_file(r'geoppi/examples/data/lineSeparation/buildings.gpkg').to_crs(cs)



# %%

lines['ID'] = np.arange(len(lines))
polys['unnID'] = np.arange(len(polys))

polys_buff = polys.copy()
polys_buff['geometry'] = polys_buff['geometry'].buffer(-0.1)

circumferentialPoints = create_circumferential_points(polygons = polys_buff, npoints = 10, lmin = 2)
circumferentialPoints['unnID'] = polys['unnID'].copy()
circumferentialPoints = circumferentialPoints.explode(index_parts = True).reset_index(drop = True)

# lines_split, splittingPoints = split_lines_at_length(lines = lines, distance = 2.5, min_distance_last_segment= 2, return_splittingPoints = True, keep_original_line_idx = True)

allJunctionPoints = lines['geometry'].apply(lambda x: create_point_splitter_npoints(x, distance = 2.5, lmin = 2))
allJunctionPoints = allJunctionPoints[~(allJunctionPoints.is_empty)].reset_index(drop = True)
allJunctionPoints = allJunctionPoints.explode(index_parts = False).reset_index(drop = True)

import matplotlib.pyplot as plt

fig, ax = plt.subplots()
allJunctionPoints.plot(ax = ax)
plt.show()



# %% Find closest points for each polygon
import pandas as pd
from shapely import STRtree
from collections import Counter

# Find closest object to each point (indicating index of closest line segment)
idxs, distances = closest_objects_to_points(points = circumferentialPoints, geomObjects = allJunctionPoints, maxDist = np.ones(len(circumferentialPoints))*100)
    
# Assign nearest street segment to each point at buildings' boundaries
circumferentialPoints['nearestPoint']                            = idxs
circumferentialPoints['distanceToPoint']                         = distances

# Find index of points which feature the shortes distance to adjacent line segment within group of points for each polygon
idx_ps_min_distances = circumferentialPoints[~circumferentialPoints['nearestPoint'].isnull()].groupby(by='unnID')['distanceToPoint'].idxmin().values
circumferentialPoints_ed = circumferentialPoints.loc[idx_ps_min_distances]
circumferentialPoints_ed['nearestPoint_geometry'] = allJunctionPoints.iloc[circumferentialPoints_ed['nearestPoint'].values].geometry.values

junctionPoints = circumferentialPoints_ed['nearestPoint_geometry'].copy()

# Create lines between buildings circumferential points and splitting points
# connectionLines = gp.GeoDataFrame(geometry = [extend_line(LineString([p1, p2]), offset = 0.1) for p1, p2 in zip(circumferentialPoints_ed['geometry'], circumferentialPoints_ed['nearestPoint_geometry'])]).set_crs(cs)

connectionLines = gp.GeoDataFrame(geometry = [extend_line(LineString([p1, p2]), offset = 0.001) for p1, p2 in zip(circumferentialPoints_ed['geometry'], circumferentialPoints_ed['nearestPoint_geometry'])]).set_crs(cs)

# Separate lines at intersections
allLines_temp = unary_union(pd.concat((lines.geometry, connectionLines.geometry), axis = 0))
allLines = gp.GeoDataFrame(geometry = list(allLines_temp.geoms)).set_crs(cs)

# Filter only these connection lines which feature a single common start- or endpoint and have a short length
shortConnectionLines = allLines[allLines.length < 0.0011]

# 2. Endpunkte der kurzen Linien extrahieren
endpoints_short = []
for line in shortConnectionLines.geometry:
    coords = np.array(line.coords)
    endpoints_short.append((Point(coords[0]), Point(coords[-1])))

start_pts_short = pd.Series([pt[0] for pt in endpoints_short], index=shortConnectionLines.index)
end_pts_short = pd.Series([pt[1] for pt in endpoints_short], index=shortConnectionLines.index)


# 3. ALLE Punkte aus gs sammeln (Start, Ende UND Stützpunkte)
all_points = []
for line in allLines.geometry:
    if not line.is_empty:
        coords = list(line.coords)  # ALLE Koordinaten (inkl. Stützpunkte)
        all_points.extend([Point(xy) for xy in coords])

# Punkt-Häufigkeiten zählen (für exakte Matches)
point_counter = Counter(all_points)
print(f"Insgesamt {len(all_points)} Punkte, {len(point_counter)} unique")

# 4. Funktion: Häufigkeit eines Punktes ermitteln
def get_point_frequency(pt):
    return point_counter.get(pt, 0)

# 5. Häufigkeiten für Start- und Endpunkte der kurzen Linien
start_freqs = np.array([get_point_frequency(pt) for pt in start_pts_short])
end_freqs = np.array([get_point_frequency(pt) for pt in end_pts_short])

all_points_arr = np.array(all_points)
start_freqs_2 = [sum(all_points_arr == pt) for pt in start_pts_short]
end_freqs_2 = [sum(all_points_arr == pt) for pt in end_pts_short]

all_freqs_2 = start_freqs_2 + end_freqs_2

# 6. Filter: ENTWEIDER Start- ODER Endpunkt taucht ≤1 auf
condition = (start_freqs <= 1) | (end_freqs <= 1)

# Ergebnis
mask = shortConnectionLines[condition]

result = allLines[~allLines.index.isin(mask.index)]



# %% Wichtig:

allLines['id'] = np.arange(len(allLines))




# %%
lines_buff = lines.copy()
lines_buff.geometry = lines_buff.geometry.buffer(0.05)

allLines_buff = allLines.copy()
allLines_buff.geometry = allLines_buff.geometry.buffer(0.05)

allLines_buff = assign_attr_by_max_intersection_area(gp1 = allLines, gp_source = lines_buff, gp1_id = 'id', attr = 'ID')

allLines_attrs = allLines.copy()
allLines_attrs['ID'] = allLines_buff['ID']


# %%
allLines_attrs.to_file(r'C:\Git\geoppi\geoppi\examples\data\lineSeparation\lines_separated_out.shp')
lines_buff.to_file(r'C:\Git\geoppi\geoppi\examples\data\lineSeparation\lines_buff.shp')

