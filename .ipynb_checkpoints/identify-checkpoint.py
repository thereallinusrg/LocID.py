'''
QTM LocID identify
Version 0.2.0 - 2023-08-28
Author: linus.rueegg@geo.uzh.ch

This script is part of the LocID project.
Its code is used to identify, as in creating an identifier for, objects.
'''

from array import array
import re
import ast

from shapely.geometry import Point, Polygon, MultiPolygon, LineString, MultiLineString, mapping

import pandas as pd
import geopandas as gpd

import QTM as q

# Polygon Mapping
def poly_g2z(pgxy):
    '''
    Converts Polygon to ZOT coordinates

    Parameters: pgxy (shapely.geometry Polygon): Polygon to convert. Polygon with WGS84 coordinates
    Returns:    pzxy (shapely.geometry Polygon): Polygon with ZOT coordinates
    '''
    poly_mapped = mapping(pgxy)
    poly_coordinates = poly_mapped['coordinates'][0]

    l= list(map(lambda x:q.GeoPoint(x[1], x[0]), poly_coordinates))

    return Polygon(list(map(lambda y:q.GeoToZot(y).raw(), l)))

# Recursive Grid Polygon Creation


def tri_split_helper(tri):
    # Helper for Recursive grid Polygon Creation
    snx = [item[0] for item in tri]
    sny = [item[1] for item in tri]
    return(snx, sny)


def tri_split(snx, sny, lvl, lst, gdf, poly):
    # Recursively determine ZOT coordinates of smallest level triangles
    # Create shapely Polygons of these lowest level triangles
    # Return GeoDataFrame of these Polygons
     
    ssnx = set(snx)
    ssny = set(sny)

    xmid = sum(ssnx) / 2
    ymid = sum(ssny) / 2
    xmax = max(ssnx)
    xmin = min(ssnx)
    ymax = max(ssny)
    ymin = min(ssny)

    if lst.count(0) % 2 != 0: # Check if we need to handle a turned triangle (each 0 facet turns everything inside)
        tri2 = [[xmax, ymid], [xmid, ymid], [xmax, ymin]]
        tri0 = [[xmax, ymid], [xmid, ymax], [xmid, ymid]] 
        tri3 = [[xmid, ymax], [xmid, ymid], [xmin, ymax]]
        tri1 = [[xmax, ymax], [xmid, ymax], [xmax, ymid]]
    else:
        tri2 = [[xmid, ymin], [xmid, ymid], [xmax, ymin]]
        tri0 = [[xmin, ymid], [xmid, ymin], [xmid, ymid]]
        tri3 = [[xmin, ymid], [xmid, ymid], [xmin, ymax]]
        tri1 = [[xmin, ymin], [xmid, ymin], [xmin, ymid]]

    if (lvl > 1):   # recursive call
        lvl = lvl - 1
        snx2, sny2 = tri_split_helper(tri2); gdf2 = tri_split(snx2, sny2, lvl, lst + [2], gdf, poly)
        snx0, sny0 = tri_split_helper(tri0); gdf0 = tri_split(snx0, sny0, lvl, lst + [0], gdf, poly)
        snx3, sny3 = tri_split_helper(tri3); gdf3 = tri_split(snx3, sny3, lvl, lst + [3], gdf, poly)
        snx1, sny1 = tri_split_helper(tri1); gdf1 = tri_split(snx1, sny1, lvl, lst + [1], gdf, poly)  
        gdf = pd.concat([gdf2, gdf0, gdf3, gdf1], ignore_index=True)

    elif (lvl > 0):
        g2 = gpd.GeoDataFrame(pd.DataFrame({'qloc': [lst+[2]]}, index = [0]), geometry =[Polygon(tri2)])
        g0 = gpd.GeoDataFrame(pd.DataFrame({'qloc': [lst+[0]]}, index = [0]), geometry =[Polygon(tri0)])
        g3 = gpd.GeoDataFrame(pd.DataFrame({'qloc': [lst+[3]]}, index = [0]), geometry =[Polygon(tri3)])
        g1 = gpd.GeoDataFrame(pd.DataFrame({'qloc': [lst+[1]]}, index = [0]), geometry =[Polygon(tri1)])

        gdf = pd.concat([gdf, g1, g2, g3, g0], ignore_index=True)  

    else:
        pass

    return gdf


def stack_max(stack):
    # Returns the smallest populated level of a stack
    for n in range(1,30):
        x = stack[n].QID
        if x == 9:
            return (n-1)
            break


def facet2shape(facet):
    # Converts a Facet Object to a shapely Polygon
    x = facet.node_x
    y = facet.node_y
    return Polygon([[x[0], y[0]], [x[1], y[1]], [x[2], y[2]]])


def stack2gdf(stacks, SeT_lvl, max_lvl):
    '''
    This function takes a stack and returns a GeoDataFrame with the shapely polygons of the facets between SeT_lvl and the max_lvl
    Parameters: 
        stacks: list of stacks
        SeT_lvl: int, level of the SeT
        max_lvl: int, max level of the stack
    Returns:
        gdf: GeoDataFrame with the shapely polygons of the facets between SeT_lvl and the max_lvl
    '''

    glst = [] # geometry list
    nlst = [] # n Facet list
    lvlst = [] # lvl list
    xlst = []
    ylst = []
    qlst = []
    xnlst = []
    ynlst = []

    for n in range(0, len(stacks)):
        for lvl in range(SeT_lvl, max_lvl+1):
            glst.append(facet2shape(stacks[n][lvl]))
            xlst.append(stacks[n][lvl].node_x)
            ylst.append(stacks[n][lvl].node_y)
            qlst.append(stacks[n][lvl].QID)
            xnlst.append(stacks[n][lvl].xnode)
            ynlst.append(stacks[n][lvl].ynode)
            nlst.append(n)
            lvlst.append(lvl)


    print(n+1, "stacks with max lvl", lvl)

    gdf = gpd.GeoDataFrame(data={'n':nlst, 'lvl':lvlst, 'node_x':xlst, 'node_y':ylst, 'QID':qlst, 'xnode':xnlst, 'ynode':ynlst},geometry=glst)
    print("Polygons in gdf:", len(gdf), "--> That is the number of facets between SeT_lvl", SeT_lvl, "and the max lvl", lvl)

    return gdf


def qloc_append(row, pqxy):
    # Helper function to append qloc to a row
    return pqxy[row['n']][:row['lvl']]


def prepend_qloc(row, qloc_base):
    # Helper function to prepend the base qloc to a qloc in a row
    return qloc_base + row['qloc']


def create_missing_facets(gdf, polygon):
    '''
    Creates a GeoDataFrame of the missing facets of a polygon
        Parameters: gdf (GeoDataFrame): GeoDataFrame of Polygons with qloc column
                    polygon (Polygon): Polygon to create grid for. Shapely Polygon with WGS84 coordinates
        Returns:    grid (GeoDataFrame): GeoDataFrame of Facets with qloc column
    '''
    polygon_gdf = gpd.GeoDataFrame(geometry=[id.poly_g2z(polygon)])
    msng_pc = polygon_gdf

    grid = gpd.GeoDataFrame()
    level = gdf.lvl.max()
    # filter df to only include polygons of one QTM level
    gdf_f = gdf[gdf['lvl'] == level]
    # filter the GeoDataFrame to only include polygons that cover the polygon
    # this is a dataframe operation 
    covering_gdf = gdf_f[gdf_f['geometry'].intersects(msng_pc['geometry'].iloc[0])]

    # append covering polygons to init_grid
    grid = pd.concat([grid, covering_gdf], ignore_index=True)

    # remove pieces which are covered by covering polygons of current level from polygon, create "missing piece"
    # this is a geometric operation, not a dataframe operation
    msng_pc = msng_pc.overlay(covering_gdf, how='difference')


    # from max_level-1 on a while loop, with cover instead of intersect
    while len(msng_pc) > 0:
        level = level - 1
        gdf_f = gdf[gdf['lvl'] == level]
        covering_gdf = gdf_f[gdf_f['geometry'].covers(msng_pc['geometry'].iloc[0])]
        grid = pd.concat([grid, covering_gdf], ignore_index=True)

        msng_pc = msng_pc.overlay(covering_gdf, how='difference')

    lvl = grid.lvl.max()
    grid_f = grid[grid['lvl'] == lvl]
    grid_ex = grid[grid['lvl'] != lvl]

    msng_pc = polygon_gdf
    msng_pc = msng_pc.overlay(grid_f, how='difference')

    for index, row in grid_ex.iterrows():
        qloc_base = row.qloc
        snx = row.node_x
        sny = row.node_y
        lst = []
        gdfin = gpd.GeoDataFrame()
        lvl_dwn = lvl - row.lvl

        out_gdf = id.tri_split(snx, sny, lvl_dwn, lst, gdfin, msng_pc)
        msng_pc = msng_pc.overlay(out_gdf, how='difference')
        out_gdf['qloc'] = out_gdf.apply(lambda row: prepend_qloc(row, qloc_base), axis=1)
        grid_f = pd.concat([grid_f, out_gdf], ignore_index=True)
    return grid_f


def map_poly(polygon, tol, lvl, bytes_out=False):
    '''
    Maps Polygon nodes to QTM coordinates and returns level of smallest enclosing triangle.

    Parameters: polygon (shapely.geometry Polygon (EPSG: 4326))
                tol (int) - tolerance level
                lvl (int) - maximum lvl difference allowed
    Returns:    pqxy (list) - list of QTM coordinates of nodes
                SeT_lvl (int) - level of smallest enclosing triangle
                stacks (list) - list of stacks of Facets
                lms (list) - list of levels of last identical facet, compared to the previous QTM coordinate
    '''
    poly_mapped = mapping(polygon)
    pzxy = poly_mapped['coordinates'][0]
    pzxy = list(dict.fromkeys(pzxy))  # remove duplicates, maintain order

    SeT_lvl = 99
    pqxy = []
    stacks = []
    max_lvl = 0

    for cord in pzxy: #This is intentionally a loop, as g2q works with a memory stack, not calculating everything again       
        qid, QL, lastMatch, stack = q.g2q(cord[1], cord[0], tol, qid_only=False, bytes_out=bytes_out)

        if lastMatch == max_lvl: continue #Skipping duplicate encodings (in the same facet)
        elif lastMatch > max_lvl: max_lvl = deepcopy(lastMatch)

        stacks.append(deepcopy(stack)) ## deepcopy is important, due to the global stack --> to not modify to global stack
        pqxy.append(qid)

        if SeT_lvl > lastMatch: SeT_lvl = lastMatch

    if max_lvl > lvl:
        max_lvl = lvl

    # turn stacks into gdf
    gdf = stack2gdf(stacks, SeT_lvl, max_lvl)

    # add qloc column with qid 
    gdf["qloc"] = gdf.apply(lambda row: qloc_append(row, pqxy), axis=1)

    # filter for duplicate facets
    gdf['hash_x'] = gdf.node_x.apply(lambda x: hash(str(x)))
    gdf['hash_y'] = gdf.node_y.apply(lambda x: hash(str(x)))
    gdf = gdf[~gdf.duplicated(['hash_x', 'hash_y'])]

    # create missing facets to intersect polygon with
    grid = create_missing_facets(gdf, polygon)

    # intersect polygon with grid
    full_pc = gpd.GeoDataFrame(geometry=[id.poly_g2z(polygon)])
    qid_lst = grid.overlay(full_pc, how='intersection').qloc.to_list()

    return qid_lst


def create_tree(lst):
    '''
    Creates a tree from a list of lists

    Parameters: lst (list): list of qloc lists
    Returns:    tree (list): LocID style facet coordinate tree
    '''
    tree = {}
    for path in lst:
        current = tree
        for node in path:
            if node not in current:
                current[node] = {}
            current = current[node]

    xs = str(tree).replace("{", "").replace(":", ",").replace(" ", "").split(',')
    xs = re.sub(r"}+", "4", str(xs))
    xl = ast.literal_eval(xs)
    return list(map(int, xl[:-1]))


def encode_polygon(polygon, tol, lvl, bytes_out=False):
    '''
    Encodes a polygon into a LocID

    Parameters: polygon (Polygon): Polygon to encode
                tol (float): tolerance for polygon encoding
                lvl (int): maximum encoding level (Maximum = 30)
    Returns:    locid (bytes): LocID of polygon
    '''
    qids = map_poly(polygon, tol, lvl, bytes_out=False)  # Polygon Mapping is only supported with bytes_out = False atm
   
    # Compress QTM coordinates to polygon LocID
    locid = create_tree(qids)

    # Return LocID - as bytes or list
    if bytes_out is True:
        locid = array('B', locid)
        return locid.tobytes()
    else:
        return locid


def identify(geo_object, tol, lvl=17):
    # lvl is not a fixed parameter, but a parameter to be determined

    # Encoding of a point
    if type(geo_object) is Point:
        lat, lon = geo_object.coords[0]
        return q.g2q(lat, lon, tol)

    # Encoding of a polygon
    elif type(geo_object) is Polygon:
        return encode_polygon(geo_object, tol, lvl)
    # TODO: Implement MultiPolygon encoding
    elif type(geo_object) is MultiPolygon:
        raise NotImplementedError("MultiPolygon type is not implemented yet")

    # Encoding of a line
    # TODO: Implement line encoding
    elif type(geo_object) is LineString:
        raise NotImplementedError("LineString type is not implemented yet")
    elif type(geo_object) is MultiLineString:
        raise NotImplementedError("MultiLineString type is not implemented yet")
    else:
        raise NotImplementedError("Unsupported object type")


