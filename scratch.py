# %%
import geopandas as gp
from pathlib import Path
from auxFunctions import nnearest


polygons = gp.read_file(Path(r'C:\Users\Voelzel\Desktop\ITEE_WS2425\Upload_heating_network_dimensioning_Python\exampleNetwork\buildings_final.gpkg'))
lines = gp.read_file(Path(r'C:\Users\Voelzel\Desktop\ITEE_WS2425\Upload_heating_network_dimensioning_Python\exampleNetwork\network_raw.gpkg'))

polygons.plot()



# %%

def detect_lines_in_narrow_passages(
    lines:gp.GeoDataFrame,
    polygons:gp.GeoDataFrame,
    threshDistance:float = 10,
    distPointsCircumference:float = 1,
    nNeighbours:int = 10,
    col:str = 'narrowPassage' 
    )->gp.GeoDataFrame:

    # Imports
    import libpysal
    import numpy as np
    from scipy.spatial import cKDTree
    from shapely import LineString
    from auxFunctions import nnearest

    # Plausibility checks
    nNeighbours = max(1, nNeighbours)

    cs = lines.crs

    # Definitions
    def qweights_wrapper(df, col : str = None):

        if col is None:
            col = 'qweights'

        # testing: This does broadcast
        # queen_weights = libpysal.weights.Queen.from_dataframe(df)
        # return (queen_weights.component_labels + df.index[0]).astype(int)

        # ... maybe this does?
        queen_weights = libpysal.weights.Queen.from_dataframe(df)
        df[col] = (queen_weights.component_labels + df.index[0]).astype(int)
        return df


    def create_points_along_boundary(polygon, dist:float = 1):

        """
        Function creating circumferential points along boundary of a polygon object.\n

        :param polygon: shapely.Polygon object.\n
        :param dist: float denoting the desired distance between circumferential points along boundary to create.\n

        :return: list of points along boundary of the polygon
        """

        boundary = polygon.boundary
        if boundary.length == 0:
            return []
        
        numPoints = int(boundary.length // dist)
        return [(np.round(boundary.interpolate(ii * dist).x, 2), np.round(boundary.interpolate(ii * dist).y, 2)) for ii in range(numPoints+1)]
           
    ## Data preparation
    # Merge all touching/overlapping polygons to reduce number of processed polygons
    print(f'\n... Merging touching/overlapping polygons\n')
    polygons = qweights_wrapper(polygons, col = 'qWeight')
    polygons = polygons.dissolve(by = 'qWeight', aggfunc = 'first')

    # Create circumferenital points along boundary of polygons
    all_points = [create_points_along_boundary(poly, dist = distPointsCircumference) for poly in polygons.geometry]

    # Initialize container for shortest line objects between neighbouring points
    shortestLines = []

    ## Start looping through all points
    for n, pts in enumerate(all_points):

        # Initialize container for lines for each polygon
        singleLines = []

        # Exclude points belonging to the same polygon from the search for neighbours
        otherpts = [pt for j in range(len(all_points)) if j != n for pt in all_points[j]]
        lOthers = len(otherpts)

        if lOthers > 0:
            idx = nnearest(A = np.array(pts), B = np.array(otherpts), distance = threshDistance, n = nNeighbours)

            for nn, p in enumerate(idx):
                singleLines.append( [(pts[nn], otherpts[ix]) for ix in p if ix < lOthers] )
                singleLines = [b for b in singleLines if b]

            shortestLines.append(singleLines)

    # Flatten list of lines
    shortestLines = [x for xs in shortestLines for x in xs]
    shortestLineObjects = gp.GeoDataFrame(geometry = [LineString(coords) for ls in shortestLines for coords in ls])
    shortestLineObjects.set_crs(cs, inplace = True)
    
    # Order all geometries in canonical form and drop duplicates (Connecting lines feature only one start- and one endpoint)
    shortestLineObjects['geometry'] = shortestLineObjects.normalize()
    shortestLineObjects.drop_duplicates()

    # Transfer results to output DF of lines
    lines = lines.copy()
    lines = lines.sjoin(shortestLineObjects, how = 'left', predicate = 'intersects')
    idxs = lines[~lines['index_right'].isna()].drop_duplicates(subset = 'geometry').index

    lines[col] = False
    lines.loc[idxs, col] = True


    return lines

# %%


lines_new = detect_lines_in_narrow_passages(
    lines = lines,
    polygons = polygons,
    nNeighbours = 5,
    distPointsCircumference = 2,
    threshDistance = 10,
    col = 'narrowPassage'
)

lines_new.to_file(Path(r'C:\Users\Voelzel\Desktop\ITEE_WS2425\Upload_heating_network_dimensioning_Python\exampleNetwork\outStreets_new.gpkg'), driver = 'GPKG')

# %%

import numpy as np
from scipy.spatial import cKDTree
from shapely.ops import unary_union
from shapely import LineString
import libpysal

import time


def qweights_wrapper(df, col : str = None):

    if col is None:
        col = 'qweights'

    # testing: This does broadcast
    # queen_weights = libpysal.weights.Queen.from_dataframe(df)
    # return (queen_weights.component_labels + df.index[0]).astype(int)

    # ... maybe this does?
    queen_weights = libpysal.weights.Queen.from_dataframe(df)
    df[col] = (queen_weights.component_labels + df.index[0]).astype(int)
    return df


def create_points_along_boundary(polygon, dist:float = 1):

    boundary = polygon.boundary
    if boundary.length == 0:
        return []
    
    numPoints = int(boundary.length // dist)
    return [(np.round(boundary.interpolate(ii * dist).x, 2), np.round(boundary.interpolate(ii * dist).y, 2)) for ii in range(numPoints+1)]

def nnearest_dist(arrA:np.array, arrB:np.array, dist:float = 1, n:int = 2):

    btree = cKDTree(arrB)
    dist, idx = btree.query(arrA, k=n, distance_upper_bound=dist)

    return dist, idx


polygons = qweights_wrapper(polygons, col = 'qWeight')
polygons = polygons.dissolve(by = 'qWeight', aggfunc = 'first')

###
starttime1p5 = time.time()

all_points = [create_points_along_boundary(poly, dist = 2) for poly in polygons.geometry]

endtime1p5 = time.time()
duration1p5 = endtime1p5 - starttime1p5
print(f'### Duration for creating points along boundaries is {duration1p5:.2}s')

shortestLines = []


# %%
###
starttime2 = time.time()

for n, pts in enumerate(all_points):

    singleLines = []

    otherpts = [pt for j in range(len(all_points)) if j != n for pt in all_points[j]]
    lOthers = len(otherpts)

    if lOthers > 0:
        dist, idx = nnearest_dist(arrA = np.array(pts), arrB = np.array(otherpts), dist = 10, n = 5)

        for nn, p in enumerate(idx):
            singleLines.append( [(pts[nn], otherpts[ix]) for ix in p if ix < lOthers] )
            singleLines = [b for b in singleLines if b]

        shortestLines.append(singleLines)

# Flatten list of lines
shortestLines = [x for xs in shortestLines for x in xs]


# %%

shortestLineObjects = gp.GeoDataFrame(geometry = [LineString(coords) for ls in shortestLines for coords in ls])
shortestLineObjects['geometry'] = shortestLineObjects.normalize()
shortestLineObjects.drop_duplicates()

shortestLineObjects.to_file(Path(r'C:\Users\Voelzel\Desktop\ITEE_WS2425\Upload_heating_network_dimensioning_Python\exampleNetwork\testLines.gpkg'), driver = 'GPKG')


# %%

lines_new = lines.sjoin(shortestLineObjects, how = 'left', predicate = 'intersects')
lines_new[~lines_new['index_right'].isna()].drop_duplicates(subset = 'geometry').to_file(Path(r'C:\Users\Voelzel\Desktop\ITEE_WS2425\Upload_heating_network_dimensioning_Python\exampleNetwork\outStreets.gpkg'), driver = 'GPKG')

# %% Suche nach Schnitten mit Straßenzügen

sindex = lines.sindex

for idx, line in shortestLineObjects.iterrows():
    possible_matches_index = list(sindex.query(line.geometry))
    possible_matches = lines.iloc[possible_matches_index]

    actual_matches = possible_matches[possible_matches.intersects(line.geometry)]

lines['narrow'] = False
# lines.loc[lines.index[actual_matches], 'narrow'] = True

# %%

lines.to_file(Path(r'C:\Users\Voelzel\Desktop\ITEE_WS2425\Upload_heating_network_dimensioning_Python\exampleNetwork\outStreets.gpkg'), driver = 'GPKG')

# %%
import shapely

shortestLineObjects.geometry.dissolve()

shapely.intersects(lines.geometry, shortestLineObjects.geometry)


# %%

# np.sqrt((all_points[0][4][0] - all_points[1][0][0])**2 + (all_points[0][4][1] - all_points[1][0][1])**2)

# %%

################################################################################
def create_point_splitter_npoints(line, npoints:int = 8, lmin:float = 5):

    """
    Create MultiPoint layer with points along idx_ps_min_distances given LineString geometry at idx_ps_min_distances predefined number of points (respecting minimum values for single line segments between points).

    Author: C. Völzel, 2023-07
    """

    import numpy as np
    from shapely.ops import unary_union

    distances = np.round(np.linspace(0, line.length, npoints+1), 1) if line.length < 200 else np.round(np.arange(0, line.length, 20), 1)
    
    # If last element is shorter than lmin, skip last splitting point
    if line.length - distances[-1] < lmin:
        distances = distances[:-1]


    points = [line.interpolate(distance) for distance in distances]
    result = unary_union(points)

    return(result)
################################################################################

# %%


