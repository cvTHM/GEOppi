# %%
import geopandas as gp
from shapely.geometry import LineString, Point, MultiPoint
from shapely.ops import nearest_points, split
import numpy as np
from tqdm import tqdm

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


    connection_gdf = gpd.GeoDataFrame(geometry=connection_lines, crs=polygons_gdf.crs)

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
        splitter_sindex = gpd.GeoSeries(extended_conns, crs=connection_lines_gdf.crs).sindex

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
        return gpd.GeoDataFrame(geometry=split_lines, crs=lines_gdf.crs)

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


    connection_gdf = gpd.GeoDataFrame(geometry=connection_lines, crs=polygons_gdf.crs)

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

        connection_gdf = gpd.GeoDataFrame(geometry=connection_lines, crs=polygons_gdf.crs)
        crossing_points_gdf = gpd.GeoDataFrame(geometry=crossing_points, crs=polygons_gdf.crs)
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
        splitter_sindex = gpd.GeoSeries(extended_conns, crs=connection_lines_gdf.crs).sindex

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
        return gpd.GeoDataFrame(geometry=split_lines, crs=lines_gdf.crs)
    
    # Aufrufe
    connection_lines_gdf, crossing_points_gdf = generate_precise_connection_lines(polygons_gdf, lines_gdf)

    split_lines_gdf = fast_split_lines_by_points(lines_gdf = lines_gdf, connection_lines_gdf = connection_lines_gdf, extend_connection = 0.1, precision = 0.01, tolerance = 1e-08)

    return connection_gdf, split_lines_gdf


from geoppi.auxFunctions import (nnearest, create_circumferential_points, create_point_splitter_npoints, split_lines_at_length, closest_objects_to_points)


lines = gp.read_file(r'geoppi/examples/data/lineSeparation/lines_separation_raw.gpkg')
polys = gp.read_file(r'geoppi/examples/data/lineSeparation/buildings.gpkg')


# %%

polys['unnID'] = np.arange(len(polys))
circumferential_points = create_circumferential_points(polygons = polys, npoints = 10, lmin = 5)
circumferential_points['unnID'] = polys['unnID'].copy()
circumferential_points = circumferential_points.explode(index_parts = True).reset_index(drop = True)

lines_split, splittingPoints = split_lines_at_length(lines = lines, distance = 10, min_distance_last_segment= 2, return_splittingPoints = True)



import matplotlib.pyplot as plt

fig, ax = plt.subplots()
lines_split.plot(ax = ax)
splittingPoints.plot(ax = ax)
plt.show()


# %% Find closest points for each polygon
import pandas as pd

# Find closest object to each point (indicating index of closest line segment)
idxs, distances = closest_objects_to_points(points = circumferential_points, geomObjects = splittingPoints, maxDist = np.ones(len(circumferential_points))*100)

seriesIdx = pd.Series(idxs)
distances = pd.Series(distances)      
    
# Assign nearest street segment to each point at buildings' boundaries
circumferential_points['nearestPoint']                            = seriesIdx.values
circumferential_points['distanceToPoint']                               = distances.values

# Find index of points which feature the shortes distance to adjacent line segment within group of points for each polygon
idx_ps_min_distances = circumferential_points[~circumferential_points['nearestPoint'].isnull()].groupby(by='unnID')['distanceToPoint'].idxmin().values
circumferential_points_ed = circumferential_points.loc[idx_ps_min_distances]




# %% Imports

# import momepy
# import geopandas as gp
# from pathlib import Path

# from geoppi.suitable_network_routing import (network_span_dfs_level_search)



# # %% Load data

# flp = Path('C:/GitLab/paper_generic_network_creation/data/lineDensity')
# flp_out = flp

# # Coordinate system
# cs = 'EPSG:25832'

# lines               = gp.read_file(flp / Path(r'Streets_LineDensity_bestand_renov_v1.gpkg')).explode(index_parts = False).to_crs(cs)
# lines['length']     = lines.geometry.length
# lines['inverse_ld_demand_use_th'] = 1/lines['ld_demand_use_th']

# startpoints         = gp.read_file(flp / Path(r'startpointsProducers.gpkg')).to_crs(cs)
 
# consumerConnections = gp.read_file(flp / Path(r'Points_consumerIntersection.gpkg')).to_crs(cs)

# buildings = gp.read_file(flp / Path(r'Buildings_LineDensity_renov_v1.gpkg')).explode(index_parts = False)

# buildings = buildings.to_crs(lines.crs)


# # %%
# # Aufruf
# connection_lines, lines_sep = create_shortest_connections_optimized_v2(polygons_gdf = buildings, lines_gdf = lines, buffer_dist = 0.1)