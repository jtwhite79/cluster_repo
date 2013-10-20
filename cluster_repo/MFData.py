import re
import sys
import os
import math
import numpy as np
from MFInterpolators import bilinear_interpolation

#import MFArrayUtil as au

def free_u1drel(length,file):
    line = file.readline().strip().split()
    data = []
    if line[0].upper() == 'CONSTANT':
        for l in range(length): data.append(float(line[1]))
        
    elif line[0].upper() == 'INTERNAL': 
        cnst = float(line[1])       
        for l in range(length):
            dline = f.readline().strip()
            data.append(float(dline) * cnst)
    
    elif line[0].upper() == 'EXTERNAL':
        raise NameError('External not supported')
    
    elif line[0].upper() == 'OPEN/CLOSE':
        cnst = float(line[2])
        f2 = open(line[1])
        for l in range(length):
            dline = f2.readline()
            #print dline,length
            data.append(float(dline.strip()) * cnst)
        f2.close
    else:
        raise TypeError('unrecognized keyword: '+line[0].upper())
    
    return data

def load_dis_file(file): 
    f = open(file,'r')
    off = re.compile('offset',re.IGNORECASE)
    
    #--read comment lines
    #--try to get the offset
    while True:
        line = f.readline()
        if line[0] != '#': break
        if off.search(line) != None:
            try:
                raw = line.strip().split('=')[-1].split(',')
                xoff = float(raw[0])
                yoff = float(raw[1])
                rotation = float(raw[2])
                offset = [xoff,yoff,rotation]
                print 'x-offset = %g y-offset = %g rotation = %g\n' % (xoff,yoff,rotation)
            except:
                print 'offset not found in dis file header...continuing'
                offset = [-999,-999]
    
    #--parse the first line
    raw = line.split()
    nlay = int(raw[0])
    nrow = int(raw[1])
    ncol = int(raw[2])
    nper = int(raw[3])
    itmuni = int(raw[4])
    lenunit = int(raw[5])  
    
    print 'nlay = {0:3d} nrow = {1:3d}, ncol = {2:3d}'.format(nlay,nrow,ncol) 
    
    #--parse the laycbd line
    line = f.readline()
    raw = line.strip().split()
    if len(raw) != nlay:
        raise IndexError('need '+str(nlay)+' entries for dataset 2')
    laycbd = []
    for r in raw: 
        laycbd.append(float(r))
    #--read delr and delc
    delr = free_u1drel(ncol,f) 
    delc = free_u1drel(nrow,f)
    #--close DIS file
    f.close()
    #--return data                        
    return offset,nlay,nrow,ncol,np.array(delr),np.array(delc)

def cell_coordinates(nrow,ncol,delr,delc):
    xcell = np.zeros((ncol),np.float)
    ycell = np.zeros((nrow),np.float)
    #--cast single value to numpy array if necessary
    try:
        i = delr.shape[0]
    except:
        t = delr
        delr = np.empty( (ncol), np.float )
        delr.fill( t )
    try:
        j = delc.shape[0]
    except:
        t = delc
        delc = np.empty( (nrow), np.float )
        delc.fill( t )
    #--node coordinates
    xcell[0] = delr[0] / 2.0
    for i in range(1,ncol):
        xcell[i] = xcell[i-1] + ( delr[i-1] + delr[i] ) / 2.0
    ycell[0] = delc[0] / 2.0
    for i in range(1,nrow):
        ycell[i] = ycell[i-1] + ( delc[i-1] + delc[i] ) / 2.0
    ycell = ycell[::-1]
    return xcell, ycell
                                               
def edge_coordinates(nrow,ncol,delr,delc):
    xedge = np.zeros((ncol+1),np.float)
    yedge = np.zeros((nrow+1),np.float)
    #--cast single value to numpy array if necessary
    try:
        i = delr.shape[0]
    except:
        t = delr
        delr = np.empty( (ncol), np.float )
        delr.fill( t )
    try:
        j = delc.shape[0]
    except:
        t = delc
        delc = np.empty( (nrow), np.float )
        delc.fill( t )
    #--edge of cells
    xedge[0] = 0.0
    for i in range(1,ncol+1):
        xedge[i] = xedge[i-1] + delr[i-1]
    yedge[0] = 0.0
    j = nrow - 1
    for i in range(1,nrow+1):
        yedge[i] = yedge[i-1] + delc[j]
        j -= 1
    yedge = yedge[::-1]
    return xedge, yedge

def twoArrayToPointArray(x,y):
    point = []
    for [xt,yt] in zip(x,y):
        point.append( [xt,yt] )
    return np.array(point)

def pointArrayToTwoArray(point):
    x = []
    y = []
    for [xt,yt] in point:
        x.append( xt )
        y.append( yt )
    return np.array(x), np.array(y)
    

def rotate(box,angle):
    new_box = []
    sin_phi = math.sin(angle*math.pi/180.0)
    cos_phi = math.cos(angle*math.pi/180.0)
    #print sin_phi,cos_phi
    for point in box:
        new_x = (point[0] * cos_phi) - (point[1] * sin_phi)
        new_y = (point[0] * sin_phi) + (point[1] * cos_phi)
        new_box.append([new_x,new_y])
    return new_box                                   


def add_offset(box,offset):
    for point in box:
        point[0] += offset[0]
        point[1] += offset[1]
    return box

def subtract_offset(box,offset):
    for point in box:
        point[0] -= offset[0]
        point[1] -= offset[1]
    return box


def get_row_col(nrow,ncol,xedge,yedge,xi,yi):
    icol = 0
    for i in range(1,ncol+1):
        xo = xedge[i]
        icol = i
        if xo > xi:
            break
    irow = 0
    for i in range(1,nrow+1):
        yo = yedge[i]
        irow = i
        if yo <= yi:
            break
    return (irow - 1), (icol - 1)

#--processing
def makeCellEdgePointsAlongLine(x,y,xedge,yedge,returnVertices=False):
    small_value = 1.0e-1
    #--build list of points along current line
    pts = []
    npts = len(x)
    dlen = 0.
    for idx in xrange(1,npts):
        x0 = x[idx-1]
        x1 = x[idx]
        y0 = y[idx-1]
        y1 = y[idx]
        a  = x1 - x0
        b  = y1 - y0
        c  = math.sqrt( math.pow(a,2.) + math.pow(b,2.) )
        #--find cells with (x0,y0) and (x1,y1)
        irow0,jcol0 = findRowCol(x0,y0,xedge,yedge)
        irow1,jcol1 = findRowCol(x1,y1,xedge,yedge)
        #--determine direction to go in the x- and y-directions
        goRight = False
        jx   = 0
        incx =  abs( small_value * a / c )
        goDown = False
        iy   = 0
        incy = -abs(small_value * b / c )
        if a == 0.: 
            incx = 0.
        elif a > 0.: 
            goRight = True
            jx   = 1
            incx *= -1.
        if b == 0.: 
            incy = 0.
        elif b < 0.: 
            goDown = True
            iy   = 1
            incy *= -1.
        #--process data
        if irow0 >= 0 and jcol0 >= 0:
            iAdd = True
            if idx > 1 and returnVertices==True: iAdd = False
            if (iAdd==True): pts.append( [x0,y0,dlen] )
        icnt = 0
        while True:
            icnt += 1
            dx  = xedge[jcol0+jx] - x0
            dlx = 0.
            if a != 0.:
                dlx = c * dx / a
            dy  = yedge[irow0+iy] - y0
            dly = 0.
            if b != 0.:
                dly = c * dy / b
            if dlx != 0. and dly != 0.:
                if abs(dlx) < abs(dly):
                    dy = dx * b / a
                else:
                    dx = dy * a / b
            xt = x0 + dx + incx
            yt = y0 + dy + incy
            dl = math.sqrt( math.pow((xt-x0),2.) + math.pow((yt-y0),2.) )
            dlen += dl
            if (returnVertices==False): pts.append( [xt,yt,dlen] )
            x0,y0 = xt,yt
            xt = x0 - 2. * incx
            yt = y0 - 2. * incy
            dl = math.sqrt( math.pow((xt-x0),2.) + math.pow((yt-y0),2.) )
            dlen += dl
            x0,y0 = xt,yt
            irow0,jcol0 = findRowCol(x0,y0,xedge,yedge)
            if irow0 >= 0 and jcol0 >= 0:
                if (returnVertices==False): pts.append( [xt,yt,dlen] )
            elif irow1 < 0 or jcol1 < 0:
                dl = math.sqrt( math.pow((x1-x0),2.) + math.pow((y1-y0),2.) )
                dlen += dl
                break
            if irow0 == irow1 and jcol0 == jcol1:
                dl = math.sqrt( math.pow((x1-x0),2.) + math.pow((y1-y0),2.) )
                dlen += dl
                pts.append( [x1,y1,dlen] )
                break
    return np.array( pts )

def makeEqualSpacePointsAlongLine(pt_dist,x,y):
    #--build list of points along current line
    pts = []
    npts = len(x)
    dlen = 0.
    remainder = 0.
    for idx in xrange(1,npts):
        x0 = x[idx-1]
        x1 = x[idx]
        y0 = y[idx-1]
        y1 = y[idx]
        a  = x1 - x0
        b  = y1 - y0
        c  = math.sqrt( math.pow(a,2.) + math.pow(b,2.) )
        if remainder != 0.:
            px = remainder * a / c
            py = remainder * b / c
            pts.append( [x0,y0,dlen] )
            x0 += px
            y0 += py
        dlen += remainder
        px = pt_dist * a / c
        py = pt_dist * b / c
        xint = np.arange(x0,x1,px)
        yint = np.arange(y0,y1,py)
        for [xt,yt] in zip(xint,yint):
            pts.append( [xt,yt,dlen] )
            dlen += pt_dist
        xdist = math.sqrt( math.pow((xt-x1),2.) + math.pow((yt-y1),2.) )
        remainder = pt_dist - xdist
        dlen += xdist
        if idx == (npts - 1):
            pts.append( [x1,y1,dlen] )
    return np.array( pts )

def findRowCol(x,y,xedge,yedge):
    #--find the modflow cell nw of the cross-section point
    jcol = -100
    for jdx,xmf in enumerate(xedge):
        if xmf > x:
            jcol = jdx - 1
            break
    irow = -100
    for jdx,ymf in enumerate(yedge):
        if ymf < y:
            irow = jdx - 1
            break
    return irow,jcol

def bilinearInterpToPoints(x,y,vindex,xcell,ycell,vdata):
    vinterp = []
    for idx,[xt,yt,vit] in enumerate( zip(x,y,vindex) ):
        #--find the modflow cell nw of the cross-section point
        jcol = -100
        for jdx,xmf in enumerate(xcell):
            if xmf > xt:
                jcol = jdx - 1
                break
        irow = -100
        for jdx,ymf in enumerate(ycell):
            if ymf < yt:
                irow = jdx - 1
                break
        if irow >= 0 and jcol >= 0:
            d = []
            for iy in xrange(irow,irow+2):
                for jx in xrange(jcol,jcol+2):
                    t = [ xcell[jx], ycell[iy], vdata[iy,jx] ]
                    d.append( t )
            vinterp.append( [vit, bilinear_interpolation( xt, yt, d )] )
    return np.array( vinterp )

def cellValueAtPoints(x,y,vindex,xedge,yedge,vdata):
    vcell = []
    for idx,[xt,yt,vit] in enumerate( zip(x,y,vindex) ):
        #--find the modflow cell containing point
        irow,jcol = findRowCol(xt,yt,xedge,yedge)
        if irow >= 0 and jcol >= 0:
            vcell.append( [vit, vdata[irow,jcol]] )
    return np.array( vcell )
           