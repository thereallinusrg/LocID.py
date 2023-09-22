'''
QTM LocID Operations (match & compare)
Version Version 0.2.1 - 2023-09-22
Author: linus.rueegg@geo.uzh.ch

This script is part of the LocID project.
Its code is used to perform quasi-spatial operations on LocIDs.
'''

import numpy as np


def split_locid(locid):
    '''
    Split locid into SeT and geometry parts

    Parameters:
        locid (bytes): locid
    Returns:
        SeT (bytearray): SeT part of locid
        geoM (bytearray): geometry part of locid
        geomType (int): geometry type    
    '''
    if type(locid) is not bytes:
        locid = bytes(locid)
        
    no_oct = locid[1:]
    oct_id = locid[:1]

    if b'\x08' in no_oct:
        geomType = 8; #print("Polygon")
        SeT, geoM = no_oct.split(b'\x08')
        SeT = oct_id + SeT

    elif b'\x07' in no_oct:
        geomType = 7; #print("Line")
        SeT, geoM = no_oct.split(b'\x07')
        SeT = oct_id + SeT
    else:
        geomType = 0; #print("Point")
        SeT = locid
        geoM = b''
    
    return SeT, geoM, geomType


def inflate_geom_tree(geoM, bytes_out=True):
    '''
    Inflates a geometry tree to a list of full QTM coordinates

    Parameters:
        geoM (bytearray): geometry part of a locid with geometry
        bytes_out (bool): if true, return a list of bytearrays, else return a list of lists
    Returns:
        tree (list): inflated geometry tree
    '''
    ba_list = geoM.split(b'\x04')
    branch = ba_list[0]

    tree = []

    for ba in ba_list:
        branch = branch[:-len(ba)] + ba
       
        if bytes_out == False: tree.append(list(branch))
        else: tree.append(branch)

    return tree


def compare_geom(SeT1, geoM1, SeT2, geoM2, SeT, diff_index):
    '''
    Compare two locids with geometry

    Parameters:
        SeT1 (bytearray): SeT part of locid 1
        geoM1 (bytearray): geometry part of locid 1
        SeT2 (bytearray): SeT part of locid 2
        geoM2 (bytearray): geometry part of locid 2
        SeT (bytearray): Smallest enclosing Triangle, enclosing both locids
        diff_index (int): index of first differing byte in SeT
   
    Returns:
        IoU (float): Intersection over Union, 1.0 = identical, 0.0 = no overlap
        SeT (bytearray): Smallest enclosing Triangle, enclosing both locids 
    '''
    #Geometry likenesses
    limb1 = SeT1[diff_index:]
    limb2 = SeT2[diff_index:]

    #inflate tree 
    twigs1 = inflate_geom_tree(geoM1)
    twigs2 = inflate_geom_tree(geoM2)

    #add differing limb to branching out inflated trees
    branches1 = list(map(lambda x: (bytes(limb1 + x)), twigs1))
    branches2 = list(map(lambda x: (bytes(limb2 + x)), twigs2))

    return branches1, branches2


def match_geom(SeT1, geoM1, SeT2, geoM2, SeT, diff_index):
    branches1, branches2 = compare_geom(SeT1, geoM1, SeT2, geoM2, SeT, diff_index)
    return len(np.intersect1d(branches1, branches2)) / len(np.union1d(branches1, branches2)), SeT


def contains_geom(SeT1, geoM1, SeT2, geoM2, SeT, diff_index):
    branches1, branches2 = compare_geom(SeT1, geoM1, SeT2, geoM2, SeT, diff_index)
    return set(branches1) & set(branches2) == set(branches2)


def change_geom(SeT1, geoM1, SeT2, geoM2, SeT, diff_index):
    branches1, branches2 = compare_geom(SeT1, geoM1, SeT2, geoM2, SeT, diff_index)
    return((len(set(branches1) - set(branches2)) + len(set(branches2) - set(branches1)))/len(branches1))

def match(locid1, locid2):
    ''' 
    Match two locids

    Parameters:
        locid1 (bytes): locid 1
        locid2 (bytes): locid 2

    Returns:
        IoU (float): Intersection over Union, 1.0 = identical, 0.0 = no overlap
        SeT (bytes): Smallest enclosing Triangle, enclosing both locids         
    '''
    SeT1, geoM1, geomType1 = split_locid(locid1)
    SeT2, geoM2, geomType2 = split_locid(locid2)

    if geomType1 != geomType2:
        raise ValueError("The two locids are not of the same geometry type.")
   
    #SeT likenesses
    shorter = min(len(SeT1), len(SeT2))
    diff_index = np.argmax(np.not_equal(SeT1[:shorter], SeT2[:shorter]))
   
    if diff_index == 0:
        diff_index = shorter

    SeT = SeT1[:diff_index] # the trunk


    # Matching of a point
    # TODO: Implement point matching
    if geomType1 == 0:
        raise NotImplementedError("Point matching is not supported yet")

    # Matching of a polygon
    elif geomType1 == 8:
        return match_geom(SeT1, geoM1, SeT2, geoM2, SeT, diff_index)
    
    # Matching of a line
    # TODO: Implement line matching
    elif geomType1 == 7:
        raise NotImplementedError("Line matching is not supported yet")

    else:
        raise ValueError("Unsupported object type to match or bad encoded locids")


def contains(locid1, locid2):
    '''
    Check if locid1 contains locid2

    Parameters:
        locid1 (bytes): locid 1 --> The larger geometry
        locid2 (bytes): locid 2 --> The smaller, contained geometry

    Returns:
        IoU (float): Intersection over Union, 1.0 = identical, 0.0 = no overlap
        SeT (bytes): Smallest enclosing Triangle, enclosing both locids        
    '''
    SeT1, geoM1, geomType1 = split_locid(locid1)
    SeT2, geoM2, geomType2 = split_locid(locid2)

    if geomType1 != geomType2:
        raise ValueError("The two locids are not of the same geometry type.")
    
    # SeT likenesses
    shorter = min(len(SeT1), len(SeT2))
    diff_index = np.argmax(np.not_equal(SeT1[:shorter], SeT2[:shorter]))
    
    if diff_index == 0:
        diff_index = shorter

    SeT = SeT1[:diff_index] # the trunk

    # Matching of a point
    # TODO: Implement point matching
    if geomType1 == 0:
        raise NotImplementedError("Point matching is not supported yet")

    # Matching of a polygon
    elif geomType1 == 8:
        return contains_geom(SeT1, geoM1, SeT2, geoM2, SeT, diff_index)
  
    # Matching of a line
    # TODO: Implement line matching
    elif geomType1 == 7:
        raise NotImplementedError("Line matching is not supported yet")

    else:
        raise ValueError("Unsupported object type to match or bad encoded locids")
  
def change(locid1, locid2):
    ''' 
    Measure change between two locids

    Parameters:
        locid1 (bytes): locid 1
        locid2 (bytes): locid 2

    Returns:
        IoU (float): Intersection over Union, 1.0 = identical, 0.0 = no overlap
        SeT (bytes): Smallest enclosing Triangle, enclosing both locids         
    '''
    SeT1, geoM1, geomType1 = split_locid(locid1)
    SeT2, geoM2, geomType2 = split_locid(locid2)

    if geomType1 != geomType2:
        raise ValueError("The two locids are not of the same geometry type.")
   
    #SeT likenesses
    shorter = min(len(SeT1), len(SeT2))
    diff_index = np.argmax(np.not_equal(SeT1[:shorter], SeT2[:shorter]))
   
    if diff_index == 0:
        diff_index = shorter

    SeT = SeT1[:diff_index] # the trunk


    # Matching of a point
    # TODO: Implement point matching
    if geomType1 == 0:
        raise NotImplementedError("Point change detection is not supported")

    # Matching of a polygon
    elif geomType1 == 8:
        return change_geom(SeT1, geoM1, SeT2, geoM2, SeT, diff_index)
    
    # Matching of a line
    # TODO: Implement line matching
    elif geomType1 == 7:
        raise NotImplementedError("Line change detection is not supported yet")

    else:
        raise ValueError("Unsupported object type to detect changes or bad encoded locids")