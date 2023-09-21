'''
QTM LocID Operations (match & compare)
Version Version 0.2.0 - 2023-08-28
Author: linus.rueegg@geo.uzh.ch

This script is part of the LocID project.
Its code is used to perform quasi-spatial operations on LocIDs.
'''

import numpy as np


def split_qlocid(qlocid):
    '''
    Split qlocid into SeT and geometry parts

    Parameters:
        qlocid (bytes): qlocid
    Returns:
        SeT (bytearray): SeT part of qlocid
        geoM (bytearray): geometry part of qlocid
        geomType (int): geometry type    
    '''
    if type(qlocid) is not bytes:
        qlocid = bytes(qlocid)
        
    no_oct = qlocid[1:]
    oct_id = qlocid[:1]

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
        SeT = qlocid
        geoM = b''
    
    return SeT, geoM, geomType


def inflate_geom_tree(geoM, bytes_out=True):
    '''
    Inflates a geometry tree to a list of full QTM coordinates

    Parameters:
        geoM (bytearray): geometry part of a qlocid with geometry
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
    Compare two qlocids with geometry

    Parameters:
        SeT1 (bytearray): SeT part of qlocid 1
        geoM1 (bytearray): geometry part of qlocid 1
        SeT2 (bytearray): SeT part of qlocid 2
        geoM2 (bytearray): geometry part of qlocid 2
        SeT (bytearray): Smallest enclosing Triangle, enclosing both qlocids
        diff_index (int): index of first differing byte in SeT
   
    Returns:
        IoU (float): Intersection over Union, 1.0 = identical, 0.0 = no overlap
        SeT (bytearray): Smallest enclosing Triangle, enclosing both qlocids 
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

def match(qlocid1, qlocid2):
    ''' 
    Match two qlocids

    Parameters:
        qlocid1 (bytes): qlocid 1
        qlocid2 (bytes): qlocid 2

    Returns:
        IoU (float): Intersection over Union, 1.0 = identical, 0.0 = no overlap
        SeT (bytes): Smallest enclosing Triangle, enclosing both qlocids         
    '''
    SeT1, geoM1, geomType1 = split_qlocid(qlocid1)
    SeT2, geoM2, geomType2 = split_qlocid(qlocid2)

    if geomType1 != geomType2:
        raise ValueError("The two qlocids are not of the same geometry type.")
   
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
        raise ValueError("Unsupported object type to match or bad encoded qlocids")


def contains(qlocid1, qlocid2):
    '''
    Check if qlocid1 contains qlocid2

    Parameters:
        qlocid1 (bytes): qlocid 1 --> The larger geometry
        qlocid2 (bytes): qlocid 2 --> The smaller, contained geometry

    Returns:
        IoU (float): Intersection over Union, 1.0 = identical, 0.0 = no overlap
        SeT (bytes): Smallest enclosing Triangle, enclosing both qlocids        
    '''
    SeT1, geoM1, geomType1 = split_qlocid(qlocid1)
    SeT2, geoM2, geomType2 = split_qlocid(qlocid2)

    if geomType1 != geomType2:
        raise ValueError("The two qlocids are not of the same geometry type.")
    
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
        raise ValueError("Unsupported object type to match or bad encoded qlocids")
  
def change(qlocid1, qlocid2):
    ''' 
    Measure change between two qlocids

    Parameters:
        qlocid1 (bytes): qlocid 1
        qlocid2 (bytes): qlocid 2

    Returns:
        IoU (float): Intersection over Union, 1.0 = identical, 0.0 = no overlap
        SeT (bytes): Smallest enclosing Triangle, enclosing both qlocids         
    '''
    SeT1, geoM1, geomType1 = split_qlocid(qlocid1)
    SeT2, geoM2, geomType2 = split_qlocid(qlocid2)

    if geomType1 != geomType2:
        raise ValueError("The two qlocids are not of the same geometry type.")
   
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
        raise ValueError("Unsupported object type to detect changes or bad encoded qlocids")