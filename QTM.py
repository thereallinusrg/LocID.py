'''
QTM LocID Encoder & Decoder       (Decoder not fully implemented yet)
Version Version 0.2.1 - 2023-11-22
Author: linus.rueegg@geo.uzh.ch

This is the implementation of Geoffrey Dutton's QTM algorithm in Python.
This is the main QTM file of the project. It contains the QTM class and all the functions needed to work with QTM IDs.
The code in this file is NOT necessarily held close to the original code in the paper. For that, please refer to QTM_original.py.
In this file edits are made to make the code more pythonic and to work seamlessly with the other functions of the QTM LocID project.
'''

import doctest
from array import array
from math import floor


class ZotPoint:
    '''
    Defines point object in ZOT coordinate space.
    Checks upon initialization that the coordinates are between -1 and 1.

    Attributes:
        x (float): x coordinate
        y (float): y coordinate
    '''

    def __init__(self, x, y):
        if not -1 <= x <= 1:
            raise ValueError("ZOT coordinates must be between -1 and 1 (x)")
        self.x = x
        if not -1 <= y <= 1:
            raise ValueError("ZOT coordinates  must be between -1 and 1 (y)")
        self.y = y

    def __str__(self):
        return f"x: {self.x}, y: {self.y}"

    def __repr__(self):
        return f"ZotPoint({self.x}, {self.y})"

    def raw(self):
        return (self.x, self.y)


class GeoPoint:
    '''
    Defines point object in geographic coordinate space.
    Checks upon initialization that the coordinates are between -90 and 90 for latitude and -180 and 180 for longitude.

    Attributes:
        lat (float): latitude
        lon (float): longitude
    '''

    def __init__(self, lat, lon):
        if not -90 <= lat <= 90:
            raise ValueError("Latitude must be between -90 and 90 (lat)")
        self.lat = lat
        if not -180 <= lon <= 180:
            raise ValueError("Longitude must be between -180 and 180 (lon)")
        self.lon = lon

    def __str__(self):
        return f"lat: {self.lat}, lon: {self.lon}"

    def __repr__(self):
        return f"GeoPoint({self.lat}, {self.lon} (lat, lon))"


class Facet:
    '''
    Defines empty facet object.

    Attributes:
        node_y (list): y-coordinates of the p-, y- and xnodes
        node_x (list): x-coordinates of the p-, y- and xnodes
        QID (int): QTM base facet number (1-8) -->  QID=9 indicates value not set
        xnode (int): basis number of the x-node
        ynode (int): basis number of the y-node
    '''

    def __init__(self, node_x=[None, None, None], node_y=[None, None, None], QID=None, xnode=None, ynode=None):
        self.node_y = node_y
        self.node_x = node_x
        self.QID = QID
        self.ynode = ynode
        self.xnode = xnode

    def __str__(self):
        return f"node_x: {self.node_x}, node_y: {self.node_y}, QID: {self.QID}, xnode: {self.xnode}, ynode: {self.ynode}"


def emptyStack(leng):
    '''
    Creates a "stack" of empty Facet objects with the given length.

    Parameters: leng (int): length of the stack
    Returns: list of Facet objects
    '''
    return [Facet() for _ in range(leng)]


def getOctant(zxy):
    '''
    Validates  ZOT coordinates;
    returns octant (1-8) occupied by zxy or neg error code 
    Error codes are from -1 to -15, indicating what was out of range. 
    Inequalities in range tests  on x  and y  determine  what  octant 
    points lying on the equator and principal  meridians  will  occupy. 
    Dutton, 1999 (Appendix A, p 151)

    Parameters: zxy (ZotPoint): ZOT coordinates
    Returns: oct (int): octant number
    '''
    error = 0
    oct = 0
    if zxy.x > 0:
        if zxy.x > 1:
            error = error - 1
        if zxy.y > 0:
            if zxy.y > 1:
                error = error - 8
            oct = 1
        else:
            if zxy.y < -1:
                error = error - 4
            oct = 2

    else:
        if zxy.x < -1:
            error = error - 2
        if zxy.y < 0:
            if zxy.y < -1:
                error = error - 4
            oct = 3
        else:
            if zxy.y > 1:
                error = error - 8
            oct = 4

    if error != 0:
        raise Exception(
            "Error in getOctant- possibly the ZOT coordinates are not in the range [-1,1]")

    # In  south hemisphere,  add 4  to octant number:
    if (abs(zxy.x) + abs(zxy.y)) > 1:
        oct = oct + 4
    return oct


def initStack(oct, leng=30, stack=emptyStack(30)):
    '''
    Inits stack used by QTMeval for a specified ZOT octant 
    Locations for octant vertex positions in ZOTspace.
    Dutton, 1999 (Appendix A, p 155)

    Parameters: oct (int): octant number
                leng (int): length of the stack
                stack (list): stack of Facet objects

    Returns: stack (list): stack of Facet objects
    '''
    # Locations for octant vertex positions in ZOTspace: ##[x,y] --> according to the table on page 156
    zPolenode = [[0., 0.], [0., 0.], [0., 0.], [0., 0.], [1., 1.], [
        1., -1.], [-1., -1.], [-1., 1.]]  # those have been altered into lists of lists
    zXnode = [[1., 0.], [1., 0.], [-1., 0.], [-1., 0.],
              [0., 1.], [0., -1.], [0., -1.], [0., 1.]]
    zYnode = [[0., 1.], [0., -1.], [0., -1.], [0., 1.],
              [1., 0.], [1., 0.], [-1., 0.], [-1., 0.]]

    octidx = oct-1
    stack[0].QID = oct  # First QTMdigit is octant num

    # Octant polenodes IDs always = 1,  but do not need to be stored
    # We do need to store the pole coords,  however:
    stack[0].node_x[0] = zPolenode[octidx][0]
    stack[0].node_y[0] = zPolenode[octidx][1]

    if oct < 5:                     # North hemisphere octants:
        # if oct > 5:                     # North hemisphere octants:
        stack[0].xnode = 3        # x-nodebasis num
        stack[0].ynode = 2        # y-nodebasis num
        ix = 2
        iy = 1              # indices for these

    else:                           # South hemisphere octants:
        stack[0].xnode = 2        # x-nodebasis num
        stack[0].ynode = 3        # y-nodebasis num
        ix = 1
        iy = 2             # indices for these

    # Now that nodes are numbered, install coords for them
    stack[0].node_x[ix] = zXnode[octidx][0]
    stack[0].node_y[ix] = zXnode[octidx][1]
    stack[0].node_x[iy] = zYnode[octidx][0]
    stack[0].node_y[iy] = zYnode[octidx][1]

    # Clear rest of stack;  QID=9 indicates value not set:

    for i in range(1, len(stack)):  # This should probably be replaced by an apply function
        stack[i].QID = 9
        stack[i].xnode = stack[i].ynode = 0

        stack[i].node_x = [0, 0, 0]  # Added this instead of the loop below
        stack[i].node_y = [0, 0, 0]
# for j in range(0,3): ## Loop seems to overwrite stack[0] with zeros
# stack[i].node_x[j] = stack[i].node_y[j] = 0
    return stack


# Global Variables initialized for QTMencode
lastMatch = 0
stack = initStack(0)


def QTMencode(zxy, tol, bytes_out=True):
    '''
    Encodes a ZOT point into a QTM ID.
    Assume zxy is a  valid ZOT coordinate in the range  (-1,-1)  to  (1,1) <br> 
    Assume tol is a valid ZOT distance in the range  (2^-29 , 1.)<br> 
    The stack may be initialized, pristine or contain data 
    Dutton, 1999 (Appendix A, p 152)

    Parameters: zxy (ZOT): ZOT coordinate
                tol (float): ZOT distance
                lastMatch (int): last match
                stack (list): stack of Facet objects
    Returns:    qid (list): QTM ID
                QL (int): QTM level
                lastMatch (int): last match
                stack (list): stack of Facet objects
    '''
    ## My comments or commented out original code lines marked with ##

    # use global variables
    global stack
    global lastMatch

    if not 2**-29 <= tol <= 1:  # Added this test
        raise ValueError("ZOT distances  must be between 2⁻²⁹ and 1 (tol)")

    oct = getOctant(zxy)  # Compute QTM Octant,  oct,  occupied by pxy:

    if oct != stack[0].QID:
        stack = initStack(oct)  # initStack  (oct,  MAXLEVS,  stack)

    qid = array('B', [oct])  # my addition
    QL = 0  # Init QTM level
    SameLoc = True  # Init  flag  for stack's  validity
    dz = 1.0  # Init ZOT x,y edge length  (of octant)

    # Loop from QTM lvl 0 to the highest that tol will allow:
    while (dz > tol):
        # Get IDs of X- and Y-nodes marked with   (xn, yn)  from level QL of stack  (1,2,3):
        xn = stack[QL].xnode
        yn = stack[QL].ynode
        pn = 6 - (xn + yn)  # Compute Polenode ID

        # Retrieve ZOT x,y of polenode, pn_x, pn_y from level QL of stack:
        pn_x = stack[QL].node_x[pn-1]
        pn_y = stack[QL].node_y[pn-1]
        # pn_x  =  stack[0].node_x[pn-1]
        # pn_y  =  stack[0].node_y[pn-1]
        dz = dz / 2  # Reduce closeness criterion by half

    # Compute displacement of zxy from polenode:
        dx = abs(zxy.x - pn_x)
        dy = abs(zxy.y - pn_y)
    # Identify closest node to zxy (node 0 represents central facet):
        node = 0
        if (dx + dy) <= dz:
            node = pn
        elif dx >= dz:
            node = xn
        elif dy >= dz:
            node = yn

        QL = QL + 1  # Increment QTM level

    # Is stack state still consistent at this level?
        if SameLoc == True and stack[QL].QID == node:
            lastMatch = QL
        else:  # Need to update rest of stack:
            SameLoc = False
            xn_x = stack[QL-1].node_x[xn-1]
            xn_y = stack[QL-1].node_y[xn-1]
            yn_x = stack[QL-1].node_x[yn-1]
            yn_y = stack[QL-1].node_y[yn-1]
            # xn_x    =  stack[0].node_x[xn-1]
            # xn_y    =  stack[0].node_y[xn-1]
            # yn_x    =  stack[0].node_x[yn-1]
            # yn_y    =  stack[0].node_y[yn-1]

    # Copy pn_x, pn_y,  xn_x,  xn_y, yn_x,  yn_v to next level overriding these values for the two nonselected nodes:
            if node == pn:  # !=
                stack[QL].node_x[pn-1] = pn_x
                stack[QL].node_y[pn-1] = pn_y
            else:
                stack[QL].node_x[pn-1] = (xn_x + yn_x) / 2
                stack[QL].node_y[pn-1] = (xn_y + yn_y) / 2

            if node == yn:  # !=
                stack[QL].node_x[yn-1] = yn_x
                stack[QL].node_y[yn-1] = yn_y
            else:
                stack[QL].node_x[yn-1] = (xn_x + pn_x) / 2
                stack[QL].node_y[yn-1] = pn_y

            if node == xn:  # !=
                stack[QL].node_x[xn-1] = xn_x
                stack[QL].node_y[xn-1] = xn_y
            else:
                stack[QL].node_x[xn-1] = pn_x
                stack[QL].node_y[xn-1] = (yn_y + pn_y) / 2

        # Renumber nodes:

            if node == xn:
                pn = yn
                yn = 6 - (xn + yn)
            elif node == yn:
                pn = xn
                xn = 6 - (pn + yn)
            elif node == pn:
                yn = xn
                xn = 6 - (yn + pn)

        # Check the block with the function "renumber" later in this document
        # if  node  ==  xn:       yn  =  6  -  (xn  +  yn)
        # elif  node  ==  yn:     xn  =  6  -  (xn  +  yn)
        # elif  node  ==  pn:    yn = xn ; xn  =  6  -  (yn +  pn)

        # Store xnode and ynode ids on stack  (pole ID will be derived):
            stack[QL].xnode = xn
            stack[QL].ynode = yn
            # Added check for 2/3 QID change for even oct ID's.
            if oct & 1 == 1:
                stack[QL].QID = node  # Store QID for this level
            else:
                if node == 2:
                    stack[QL].QID = 3
                elif node == 3:
                    stack[QL].QID = 2
                else:
                    stack[QL].QID = node
            # stack[QL].QID       =  node

        qid.append(stack[QL].QID)  # type: ignore # my addition
    if bytes_out == True:
        qid = qid.tobytes()  # my addition
    else:
        qid = qid.tolist()
    return qid, QL, lastMatch, stack


def GeoToZot(gxy):
    '''
    Converts geographic coordinates to Zot coordinates.
    Dutton, 1999 (Appendix A, p 158)

    Parameters: gxy (GeoPoint)
    Returns:    zxy (ZotPoint)
    '''
    # DOUBLE Ion,  londeg,  dx,  dy,  dxy,  temp
    # Assume /geoPoint.lat/  .LE. 90; Assume /geoPoint.lon/  .LE.  180

    dxy = 1. - (abs(gxy.lat) / 90.)      # Get fractional colatitude

    lon = gxy.lon
    if lon < 0:
        lon = lon + 360.                # Make longitude positive

    # get degree part of longitude
    londeg = floor(lon)

    # Normalize lon to  (0.-> 90. )
    dx = (londeg % 90) + lon - londeg
    dx = dx * (dxy / 90.)                      # Derive Manhattan x-distance
    dy = dxy - dx                                   # Derive Manhattan y-distance
    if gxy.lat < 0.:                               # Reverse dir.  in s. hemi
        temp = 1. - dx
        dx = 1. - dy
        dy = temp

    if (floor(lon) / 90) % 2 == 1:            # Reverse x,y if lon in (90, 180, 270)
        temp = dx
        dx = dy
        dy = temp

    # Negatize Y value in top half of ZOT space:
    if lon > 90. and lon <= 270.:
        dy = -dy

    # Negatize X  value in left half of ZOT space:
    if lon > 180.:
        dx = -dx

    zotPoint = ZotPoint(x=dx, y=dy)

    return zotPoint


def g2q(lat, lon, tol, qid_only=True, bytes_out=True):
    '''
    Directly converts geographic coordinates to QTM ID.

    lat, lon: Order of coordinates according to ISO 6709.
    lon, lat: Order of coordinates according to GeoPandas:

    Parameters: lat (float), lon (float), tol (float)
                lastMatch (int), stack (list)
    Returns:    qid (list), QL (int), lastMatch (int), stack (list)

    Tests:
    Zurich
        >>> g2q(47.38, 8.53, 0.0001, bytes_out = False)
        [1, 1, 3, 3, 0, 1, 3, 1, 3, 1, 3, 1, 1, 1, 0]

    Athens
        >>> g2q(37.97, 23.72, 0.0001, bytes_out = False)
        [1, 0, 3, 0, 0, 3, 0, 3, 3, 2, 2, 2, 2, 2, 2]

    Santiago MEX (ORD41)
        >>> g2q(25.42, -99.85, 0.0001, bytes_out = False)
        [3, 3, 2, 1, 2, 1, 3, 1, 1, 1, 3, 0, 3, 1, 3]

    Athens Maine
        >>> g2q(44.95, -69.67, 0.0001, bytes_out = False)
        [4, 0, 2, 2, 3, 3, 3, 1, 1, 3, 3, 2, 0, 2, 0]

    Berlin ZA
        >>> g2q(-32.90, 27.58, 0.0001, bytes_out = False)
        [5, 0, 3, 2, 0, 2, 2, 0, 3, 0, 2, 1, 2, 0, 0]

    Sydney
        >>> g2q(-33.87, 151.22, 0.0001, bytes_out = False)
        [6, 0, 3, 0, 3, 2, 2, 3, 3, 2, 3, 0, 2, 1, 2]

    Santiago Bolivia
        >>> g2q(-18.32, -59.57, 0.0001, bytes_out = False)
        [8, 3, 1, 3, 3, 0, 0, 2, 0, 1, 1, 1, 2, 1, 1]
    '''
    gxy = GeoPoint(lat, lon)
    zxy = GeoToZot(gxy)
    qid, QL, lastMatch, local_stack = QTMencode(zxy, tol, bytes_out=bytes_out)

    if qid_only == True:
        return qid
    else:
        return qid, QL, lastMatch, local_stack


doctest.testmod(verbose=False)
