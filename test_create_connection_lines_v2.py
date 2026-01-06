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

from geoppi.auxFunctions import (nnearest, create_circumferential_points, create_point_splitter_npoints, split_lines_at_length, closest_objects_to_points)

cs = 'EPSG:25832'


lines = gp.read_file(r'geoppi/examples/data/exampleNetwork2/streets.gpkg').to_crs(cs)
polys = gp.read_file(r'geoppi/examples/data/exampleNetwork2/buildings.gpkg').to_crs(cs)

testpoly = gp.read_file(r'geoppi/examples/data/exampleNetwork2/testpoly.gpkg').to_crs(cs)

# %%

def assign_attr_by_max_intersection_area_v2(gp1:gp.GeoDataFrame, gp_source:gp.GeoDataFrame, attr:str, gp1_id:str='id'):

    '''
    Function to assign attr from GeoDataFrame gp_source to input-GeoDataFrame gp1 by maximum intersection area between each object of gp1 and gp_source.

    :param gp1: GeoDataFrame containing objects to which attributes from gp_source shall be matched.\n
    :param gp_source: GeoDataFrame containing attributes which shall be transferred to objects in gp1 with max. intersection.\n
    :param attr: str or list of strings denoting the attributes that shall be transferred.\n
    :param gp1_id: str denoting a unique identifier column name for objects in gp1.\n        
    '''

    import geopandas as gp

    dictemp = dict(zip(gp1[gp1_id], gp1.index))

    # Cut selection to those attributes which are contained in gp_source
    if not isinstance(attr, list):
        if isinstance(attr, str):
            attr = [attr]
    
    attr = [at for at in attr if at in gp_source.columns]
    
    # Drop attributes already contained in gp1
    ex_attr = [at for at in attr if at in gp1.columns]
    if len(ex_attr) > 0:
        print(f'\nAttributes {ex_attr} already contained in gp1. Columns are dropped and overwritten.')

    gp1.drop(columns = ex_attr, inplace = True)

    gtemp = gp.overlay(gp1, gp_source[attr + ['geometry']], how='intersection', keep_geom_type = False)

    gtemp['temporary'] = gtemp.length

    if (all(gp1.geom_type.isin(['LineString', 'MultiLineString'])) & all(gp_source.geom_type.isin(['Polygon', 'MultiPolygon']))) | (all(gp_source.geom_type.isin(['LineString', 'MultiLineString'])) & all(gp1.geom_type.isin(['Polygon', 'MultiPolygon']))):        
        gtemp['temporary'] = gtemp.length

    elif (all(gp1.geom_type.isin(['Polygon', 'MultiPolygon'])) & all(gp_source.geom_type.isin(['Polygon', 'MultiPolygon']))):
        gtemp['temporary'] = gtemp.area

    idxmax = gtemp.groupby(by=gp1_id)['temporary'].idxmax()
    ls = gtemp.loc[idxmax.values, [gp1_id] + attr].set_index(gp1_id, drop=True)
    ls.index = [dictemp[p] for p in ls.index]

    gp1[attr] = ls

    return (gp1)




polys = assign_attr_by_max_intersection_area_v2(gp1 = polys, gp_source = testpoly, attr = ['name', 'fixedname'], gp1_id = 'unique_ID_polys')


polys.to_file(r'geoppi/examples/data/exampleNetwork2/polys_out.gpkg', driver  = 'GPKG')