#--2d real array utilities 
#--for reading, writing and plotting
#--ASCII 2-D MODFLOW real arrays
#--libraries
import numpy as np
import os
import math
import sys
import gc
import time

#import blnUtil
#reload(blnUtil)

no_data = -999

def TestDirExist(ctest):
    #--Evaluate if the directories in ctest exist.
    #  Create the directories in ctest if they do not already exist.
    for f in ctest:
        sys.stdout.write( 'evaluating...\n  "{0}"\n'.format( os.path.dirname( f ) ) )
        fa = os.path.abspath(f)
        d = os.path.dirname(fa)
        if not os.path.exists(d):
            sys.stdout.write( 'creating directory path...\n  "{0}"\n'.format( os.path.dirname( f ) ) )
            os.makedirs(d)


def mapBndCellsArray(nrow,ncol,bndcells,**kwargs):  
    try:
        parm = kwargs['parm']
    except:
        parm = ''
    array = np.zeros((nrow,ncol),dtype='double')-1.0e+30
    for cells in range(0,len(bndcells)): 
        if cmp(parm,'stage') == 0:
            array[(bndcells[cells].row)-1,(bndcells[cells].column)-1] = bndcells[cells].stage
        elif cmp(parm,'rate') == 0:
            array[(bndcells[cells].row)-1,(bndcells[cells].column)-1] = bndcells[cells].rate
        elif cmp(parm,'concen') == 0:
            array[(bndcells[cells].row)-1,(bndcells[cells].column)-1] = bndcells[cells].concen
        elif cmp(parm,'cond') == 0:
            array[(bndcells[cells].row)-1,(bndcells[cells].column)-1] = bndcells[cells].cond
        elif cmp(parm,'int') == 0:
            array[(bndcells[cells].row)-1,(bndcells[cells].column)-1] = 1
        else:        
            try:
                array[(bndcells[cells].row)-1,(bndcells[cells].column)-1] = bndcells[cells].stage
            except:
                try:
                    array[(bndcells[cells].row)-1,(bndcells[cells].column)-1] = bndcells[cells].rate
                except:
                    try:
                        array[(bndcells[cells].row)-1,(bndcells[cells].col)-1] = bndcells[cells].concen
                    except:
                        raise TypeError 
    return array




def loadArrayFromFile(nrow,ncol,file,ctype=None):
    '''
    read 2darray from file
    file(str) = path and filename
    '''
    try:
        file_in = open(file,'r')
        openFlag = True
    except:
#       assert os.path.exists(file)
        file_in = file
        openFlag = False
    
    data = np.zeros((nrow*ncol),dtype=np.float) #-1.0E+10
    data.fill( -1.0E+10 )
    d = 0
    while True:
        line = file_in.readline()
        if line is None or d == nrow*ncol:break
        raw = line.strip('\n').split()
        for a in raw:
            try:
                data[d] = float(a)
            except:
                print 'error casting to float on line: ',line
                sys.exit()
            if d == (nrow*ncol)-1:
                assert len(data) == (nrow*ncol)
                #data.resize(nrow,ncol)
                #return(data) 
            d += 1  
    file_in.close()
    if ctype != None:
        if ctype.lower() == 'i':
            data = data.astype(np.int32)
    data.resize(nrow,ncol)
    return(data)

def loadRealArrayFromBinaryFile(nrow,ncol,file,realbyte=4):
    '''
    read 2darray from file
    file(str) = path and filename
    '''
    if realbyte == 4:
        dtype_np = np.float
        tsize    = np.float32
    elif realbyte == 8:
        dtype_np = np.double    
        tsize    = np.float64
    try:
        f = open(file,'rb')
    except:
#       assert os.path.exists(file)
        f = file
        return False, np.zeros( (nrow,ncol), dtype_np )
    
    data=np.fromfile(file=f, dtype=tsize, count=nrow*ncol)
    data.shape=(nrow,ncol)
    f.close()
    return True, data

def loadSurferGrdFromFile(file):
    '''
    read 2darray from file
    file(str) = path and filename
    '''
    try:
        file_in = open(file,'r')
        openFlag = True
    except:
#       assert os.path.exists(file)
        file_in = file
        openFlag = False
    
    line = file_in.readline()  #dssa
    line = file_in.readline()  #ncol nrow
    raw = line.split()
    try:
        ncol = int( raw[0] )
        nrow = int( raw[1] )
    except:
        print 'error parsing ncol and nrow on line:\n  ',line
        sys.exit()
    line = file_in.readline()  #xmin xmax
    raw = line.split()
    try:
        xmin = float( raw[0] )
        xmax = float( raw[1] )
    except:
        print 'error parsing xmin and xmax on line:\n  ',line
        sys.exit()
    line = file_in.readline()  #ymin ymax
    raw = line.split()
    try:
        ymin = float( raw[0] )
        ymax = float( raw[1] )
    except:
        print 'error parsing ymin and ymax on line:\n  ',line
        sys.exit()
    line = file_in.readline()  #zmin zmax
    raw = line.split()
    try:
        zmin = float( raw[0] )
        zmax = float( raw[1] )
    except:
        print 'error parsing zmin and zmax on line:\n  ',line
        sys.exit()
    
    #calculate cell size and offset
    dx = ( xmax - xmin ) / float( ncol - 1 )
    dy = ( ymax - ymin ) / float( nrow - 1 )
    x0 = xmin - 0.5 * dx
    y0 = ymin - 0.5 * dy
    offset   = [ x0, y0 ]
    cellsize = [ dx, dy ]
    
    data = np.zeros((nrow*ncol),dtype='double')-1.0E+10
    d = 0
    while True:
        line = file_in.readline()
        if line is None or d == nrow*ncol:break
        raw = line.strip('\n').split()
        for a in raw:
            try:
                data[d] = float(a)
            except:
                print 'error casting to float on line: ',line
                sys.exit()
            if d == (nrow*ncol)-1:
                assert len(data) == (nrow*ncol)
                data.resize(nrow,ncol)
                #return nrow,ncol,data
                exit
            d += 1  

    file_in.close()
    #print nrow, ncol

    data.resize(nrow,ncol)
    #invert data to modflow orientation
    data = np.flipud(data)
    #mask where less than zmin
    data = np.ma.masked_less( data, zmin )
    #mask where greater than zmax
    data = np.ma.masked_greater( data, zmax )
    return nrow,ncol,offset,cellsize,data


def ref2grd(file,ref,nrow,ncol,offset,delt,nodata=-999):
    f = open(file,'w')
    f.write('ncols '+str(ncol)+'\n')
    f.write('nrows '+str(nrow)+'\n')
    f.write('xllcorner '+str(offset[0])+'\n')
    f.write('yllcorner '+str(offset[1])+'\n')
    f.write('cellsize '+str(delt)+'\n')
    f.write('nodata_value {0}\n'.format( nodata ))
    writeArrayToFile(ref,f,nWriteCol=ncol) 
    f.close()
    return
        

    
def writeArrayToFile(array,file,**kwargs):
    '''
    write 2-d array to file
    nWriteCol(int) = number of columns in output file array
    oFormat(str) = output format ex: {0:12.4E}
    file(str) = path and filename
    '''
    
    #--get keyword arguments
    try:
        nWriteCol = kwargs['nWriteCol']
    except:
        nWriteCol = 10
    
    try:
        oFormat = kwargs['oFormat']
        if len(oFormat) == 1:
            if cmp(oFormat,'i') == 0 or cmp(oFormat,'I') == 0:
                oFormat = '{0:3.0f}'
    except:
        oFormat = '{0:14.6E}'

    try:
        oFormat = kwargs['cFormat']
    except:
        cFormat = False
    
    
    assert len(np.shape(array)) == 2
    nrow,ncol = np.shape(array)
    #--check for line return flag at end each column
    if ncol%nWriteCol == 0:
        lineReturnFlag = False
    else:
        lineReturnFlag = True   
        
    #--try to open file, if fail, file is an open fileobj       
    try:
        file_out = open(file,'w')
        openFlag = True
    except:
        file_out = file
        openFlag = False
            
    #--write the array                      
    for a in range(0,nrow):
        for b in range(0,ncol):
            try:
                file_out.write(oFormat.format(float(array[a][b])))
            except:
                print 'NAN at row,col: ',a,b,array[a][b]
                sys.exit()
            if (b+1)%nWriteCol == 0.0 and b != 0:
                file_out.write('\n')
        if lineReturnFlag == True:
            file_out.write('\n')
    if openFlag == True:        
        file_out.close()
    return True




#--generic plot of 2d array with grid
def plotArray(array,rowDim,colDim,**kwargs):
    import pylab as pl
    import matplotlib as mpl
    from matplotlib import colors, cm
    #--get keyword arguments
    try:
        title = kwargs['title']
    except:
        title = ''

    try:
        cBarLoc,cBarLabel = kwargs['cBarLoc'],kwargs['cBarLabel']
    except:
        cBarLoc,cBarLabel = [0.25,0.025,0.5,0.025],''
    try:
        cbticks = kwargs['cbticks']
        cbticks_len = len( cbticks )
    except:
        cbticks = None
        cbticks_len = 0
    
    
    try:
        xOff,yOff = kwargs['offset'][0],kwargs['offset'][1]
    except:
        xOff,yOff = 0.0,0.0     
    try:
        vmax = kwargs['vmax']
    except:
        vmax = np.max(array)
    
    try:
        vmin = kwargs['vmin']
    except:
        vmin = np.min(array)
    print 'min,max',vmin,vmax
    try:
        figuresize = kwargs['figuresize']
    except:
        if np.shape(array)[0] > np.shape(array)[1]:
            figuresize = (8.5,11)
        else:
            figuresize = (11,8.5)
    
    try:
        gridFlag = kwargs['gridFlag']
    except:
        gridFlag = False
        
    try:
        outputFlag = kwargs['outputFlag']
    except:
        try:
            outputFlag = kwargs['output']
        except:
            outputFlag = 'show'
    
    try:
        cmap = kwargs['cmap']
    except:
        cmap='jet'
    
    #--get contour array
    try:
        con_array = kwargs['con_array']
        con_len = len(con_array)
        con_text = []
        for c in con_array:
            ct = '{0:15.7g}'.format( c ).strip()
            con_text.append( ct )
    except:
        con_array = []
        con_len = 0
        con_text = []
        
    #--get bln lines and index array
    try:
        blnlines = kwargs['bln']
#       blnPoints,blnIdx = blnUtil.loadBlnFile(blnlines)
        blnPoints,blnIdx = loadBlnFile(blnlines)
    except:
        blnPoints = []
        blnIdx = []
    
    #-- get generic points to plot
    try:
        genPoints = kwargs['gpts']
    except:
        genPoints = []
    
    try:
        fig = kwargs['fig']
    except:
        fig = pl.figure(figsize=(figuresize))
    
    try:
        Aspect = kwargs['Aspect']
    except:
        Aspect = 'auto'

    try:
        ax = kwargs['ax']
    except:
        ax = fig.add_subplot(111,aspect=Aspect)

    try:
        isample = kwargs['isample']
    except:
        isample = 1
    
    try:
        LogTransform = kwargs['LogTransform']
    except:
        LogTransform = False
    
    try:
        plotColorBar = kwargs['plotColorBar']
    except:
        plotColorBar = True
    
    try:
        polyline_list = kwargs['polyline']
        IsPolyline = True
    except:
        polyline_list = []
        IsPolyline = False
    try:
        polyline_color = kwargs['polyline_color']
    except:
        polyline_color = 'k'
    
                        
    #--array params
    nrow,ncol = np.shape(array)
    
    #--log transform
    if LogTransform == True:
        mask = np.ma.getmask( array )
        mc = np.ma.count_masked( array )
        if mc < 1:
            mask = np.ma.make_mask_none((nrow,ncol), dtype=None)
        for i in range(0,nrow):
            for j in range(0,ncol):
                if mask[i,j] == False:
#                    array[i,j] = math.log10( array[i,j] )
                    try:
                        array[i,j] = math.log10( array[i,j] )
                    except:
                        print i,j,array[i,j]
                        return
        if con_len > 0:
            for i in range(0,con_len):
                con_array[i] = math.log10( float( con_array[i] ) )
        try:
            vmax = math.log10( kwargs['vmax'] )
        except:
            vmax = math.ceil( np.max(array) )
        try:
            vmin = math.log10( kwargs['vmin'] )
        except:
            vmin = math.floor( np.min(array) )
        if cbticks_len > 0:
            cbticks = array( cb_ticks )
            t = copy( cbticks )
            for i in range(0,cbticks_len):
                t[i] = math.log10( t[i] )
            cbticks = copy( t )
        else:
            cbticks = np.arange(vmin,vmax+0.1,1)
            if cbticks[len(cbticks)-1] != vmax:
                cbticks = np.append( cbticks, vmax )
        print 'revised min,max',vmin,vmax

    #--set x and y dimensions       
    try:
        x = np.arange(xOff,xOff+(ncol*colDim)+colDim,colDim)
        xmin,xmax = xOff,xOff+(ncol*colDim)
    except:
        try:            
            f_in = open(colDim,'r')
            x = np.zeros((1),dtype='float')
            for line in f_in:
                raw = line.split()
                for a in range(0,len(raw)):
                    x = np.append(x,float(raw[a]))
            x = np.delete(x,[0])
            xmin,xmax = np.min(x),np.max(x)                         
        except:
            assert np.shape(colDim)[0] == ncol+1
            x = colDim
            xmin,xmax = np.min(x),np.max(x)
            
    try:
        y = np.arange(yOff,yOff+(nrow*rowDim)+rowDim,rowDim)
        ymin,ymax = yOff,yOff+(nrow*rowDim)
    except:
        try:            
            f_in = open(rowDim,'r')
            y = np.zeros((1),dtype='float')
            for line in f_in:
                raw = line.split()
                for a in range(0,len(raw)):
                    y = np.append(y,float(raw[a]))
            y = np.delete(y,[0])
            ymin,ymax = np.min(y),np.max(y)                         
        except:
            assert np.shape(rowDim)[0] == nrow+1
            y = rowDim
            ymin,ymax = np.min(y),np.max(y)

    #--set x and y dimensions at cell centers, if necessary
    if con_len > 0:
        xcell = np.zeros( (ncol), np.float )
        for i in xrange(0,ncol):
            xcell[i] = 0.5 * ( x[i] + x[i+1] )
        ycell = np.zeros( (nrow), np.float )
        for i in xrange(0,nrow):
            ycell[i] = 0.5 * ( y[i] + y[i+1] )
                
    #array = np.flipud(array)
    #fig = figure(figsize=figuresize)
    #ax = subplot(1,1,1,aspect='equal')
    ax.set_title(title)     
            
    #--define meshgrid
    #X,Y = np.meshgrid(x,y)
    
    #--set up color map
    numColors = 128
    palette = cm.get_cmap(cmap,numColors)
    palette.set_over('w')
    palette.set_under('w')
    palette.set_bad('w')
    
    
    #--mask
#   array = ma.masked_where(array<vmin,array)
#   array = ma.masked_where(array>vmax,array)
    
    #-- plot array
#    c = ax.pcolor(X[::isample,::isample],Y[::isample,::isample],array[::isample,::isample],vmin=vmin,vmax=vmax,cmap=palette,alpha=0.5,edgecolors='None')
    c = ax.pcolor(x,y,np.flipud(array),vmin=vmin,vmax=vmax,cmap=palette,alpha=0.5,edgecolors='None')
    
    #-- plot grid if gridFlag
    if gridFlag:
        row = 0
        while row < nrow:
            plot([xmin,xmax],[y[row],y[row]],'k-',linewidth=0.1)
            row += 1    
        col = 0
        while col < ncol:
            plot([x[col],x[col]],[ymin,ymax],'k-',linewidth=0.1)
            col += 1
    #print xmin,xmax,ymin,ymax
    
    
    #--plot BLN lines
    if len(blnPoints) > 0:
        for a in range(1,np.shape(blnIdx)[0]):
            #print blnIdx[a-1]
            ax.plot(blnPoints[blnIdx[a-1]:blnIdx[a],0],blnPoints[blnIdx[a-1]:blnIdx[a],1],'k-',lw=2.0)

    #--plot polyline shapedata            
    if IsPolyline == True:
        polyline_plot( ax, polyline_list, polyline_color )

    #--plot generic points
    if len(genPoints) > 0:
    
        #try:
        ax.plot(genPoints[:,0],genPoints[:,1],'k+')
        #except:
        #   print 'error plotting generic points...'
    
    #--plot contours
    if con_len > 0:
        cs = ax.contour(xcell,ycell,np.flipud(array),con_array,colors='k')
        fmt = {}
        for l,s in zip( cs.levels, con_text ):
            fmt[l] = s
        ax.clabel(cs,inline=1,fmt=fmt)
        
    #--plot color bar
    if plotColorBar:
        cax=mpl.pyplot.axes(cBarLoc)
        fig.colorbar(c,cax=cax,orientation='horizontal',ticks=cbticks)                                       
        cax.set_title(cBarLabel) 
        if LogTransform == True:
            lt1 = cax.get_xticklabels()
            lt2 = []
            for i,t in enumerate( lt1 ):
#                rv = float ( t.get_text() )
                rv = cbticks[i]
                v = math.pow( 10., rv )
                if rv >= 6.0:
                    t = '{0:10.3e}'.format( v ).strip()
                elif rv >= 0.0 :
                    t = '{0:10.0f}'.format( v ).strip()
                elif rv >= -1.0:
                    t = '{0:10.1f}'.format( v ).strip()
                elif rv >= -2.0:
                    t = '{0:10.2f}'.format( v ).strip()
                elif rv >= -3.0:
                    t = '{0:10.3f}'.format( v ).strip()
                else:
                    t = '{0:10.3e}'.format( v ).strip()
                lt2.append( t )
            cax.set_xticklabels( lt2 )

    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin,ymax)

    if outputFlag == 'save':
        if title == '':
            title = str(time.time())
        plotname=title+'.png'
        fig.savefig(plotname,orientation='portrait',format='png',dpi=150)
    elif outputFlag == 'show':
        show()
    elif outputFlag == None:
        #print 'no output produced...'
        return ax
    else:
        try:
            fmt = outputFlag.split('.')[-1]
        
        except:
            print 'unable to parse outputFlag for output format type: ',outputFlag
        fig.savefig(outputFlag,orientation='portrait',format=fmt,dpi=150)
    c = []
    mpl.pyplot.close('all')
    gc.collect()
    return


def loadBlnFile(file):
    
    data = np.array([0,0,0])
    dataIdx = np.array([0])
    f_in = open(file,'r')
    while True:
        try:
            points = int(f_in.readline().strip('\n').split(',')[0])
        except:
            break               
        for a in range(0,points):
            point = f_in.readline().strip('\n').split(',')
            data = np.vstack((data,[float(point[0]),float(point[1]),float(point[2])]))
        dataIdx = np.append(dataIdx,dataIdx[-1]+points)
    data = np.delete(data,0,axis=0)
    return(data,dataIdx) 

def UpscaleArray(nrow0,ncol0,array,iscale):
    nrow1 = nrow0 / iscale
    ncol1 = ncol0 / iscale 
    if iscale == 1:
        return array 
    upscale = np.zeros( (nrow1,ncol1), np.float )
    for irow in range(0,nrow1):
        for jcol in range(0,ncol1):
            i0 = irow * iscale
            i1 = min( nrow0 - 1, (irow+1)*iscale )
            j0 = jcol * iscale
            j1 = min( ncol0 - 1, (jcol+1)*iscale )
            slice = array[i0:i1,j0:j1]
            upscale[irow,jcol] = np.mean( slice )
    return upscale

def DownscaleArray(nrow0,ncol0,array,iscale):
    nrow1 = nrow0 * iscale
    ncol1 = ncol0 * iscale
    if iscale == 1:
        return array 
    downscale = np.zeros( (nrow1,ncol1), np.float )
    for irow in range(0,nrow1,iscale):
        for jcol in range(0,ncol1,iscale):
            i0 = irow / iscale
            j0 = jcol / iscale
            downscale[irow:irow+iscale,jcol:jcol+iscale] = array[i0,j0]
    return downscale


#--generic image plotting function
def plotimage( d, **kwargs ):
    import pylab as pl
    import matplotlib as mpl
    from matplotlib import colors, cm
    try:
        figuresize = kwargs['figuresize']
    except:
        figuresize = (8.5,11)
    try:
        vmin = kwargs['vmin']
    except:
        vmin = np.min( d )
    try:
        vmax = kwargs['vmax']
    except:
        vmax = np.max( d )
    try:
        text = kwargs['text']
    except:
        text = 'None'
    try:
        fout = kwargs['fout']
    except:
        fout = 'temp.png'
    try:
        interpolation = kwargs['interpolation']
    except:
        interpolation = 'none'
    try:
        cmap = kwargs['cmap']
    except:
        cmap = 'jet'
    try:
        addcolorbar = kwargs['addcolorbar']
        if addcolorbar.lower() != 'vertical' and addcolorbar.lower() != 'horizontal':
            addcolorbar = 'None'
    except:
        addcolorbar = 'None'
    try:
        polyline_list = kwargs['polyline']
        IsPolyline = True
    except:
        polyline_list = []
        IsPolyline = False
    try:
        polyline_color = kwargs['polyline_color']
    except:
        polyline_color = 'k'
    try:
        im_extent = kwargs['extent']
    except:
        im_extent = None
        
    #--plot data
    fig = pl.figure(figsize=(figuresize))
    ax = fig.add_subplot(111,aspect=1)
    dv = ax.imshow(d,vmin=vmin,vmax=vmax,interpolation=interpolation,extent=im_extent,cmap=cmap)
    if IsPolyline == True:
        polyline_plot( ax, polyline_list, polyline_color )
    if addcolorbar != 'None':
        cb = fig.colorbar(dv,orientation=addcolorbar,ticks=[vmin,vmax])
        if addcolorbar.lower() == 'vertical':
            cb.ax.set_yticklabels( (vmin,vmax), size=6 )
        else:
            cb.ax.set_xticklabels( (vmin,vmax), size=6 )
    if im_extent != None:
        ax.set_xlim(im_extent[0],im_extent[1])
        ax.set_ylim(im_extent[2],im_extent[3])
    mpl.pyplot.xticks( size=6 )
    mpl.pyplot.yticks( size=6 )
    if text != 'None':
        ax.text(0.0, 1.01,text,size=8,horizontalalignment='left', verticalalignment='bottom', transform = ax.transAxes)
    fig.savefig(fout,dpi=150)
    mpl.pyplot.close('all')
    gc.collect()
    return

#--generic pcolor plotting function
def plotgrid( x, y, d, **kwargs ):
    import pylab as pl
    import matplotlib as mpl
    from matplotlib import colors, cm
    try:
        figuresize = kwargs['figuresize']
    except:
        figuresize = (8.5,11)
    try:
        vmin = kwargs['vmin']
    except:
        vmin = np.min( d )
    try:
        vmax = kwargs['vmax']
    except:
        vmax = np.max( d )
    try:
        text = kwargs['text']
    except:
        text = 'None'
    try:
        fout = kwargs['fout']
    except:
        fout = 'temp.png'
    try:
        cmap = kwargs['cmap']
    except:
        cmap = 'jet'
    try:
        addcolorbar = kwargs['addcolorbar']
        if addcolorbar.lower() != 'vertical' and addcolorbar.lower() != 'horizontal':
            addcolorbar = 'None'
    except:
        addcolorbar = 'None'
    try:
        polyline_list = kwargs['polyline']
        IsPolyline = True
    except:
        polyline_list = []
        IsPolyline = False
    try:
        polyline_color = kwargs['polyline_color']
    except:
        polyline_color = 'k'
    #--plot data
    fig = pl.figure(figsize=(figuresize))
    ax = fig.add_subplot(111,aspect=1)
    dv = ax.pcolor(x, y[::-1], np.flipud(d),vmin=vmin,vmax=vmax,cmap=cmap)
    if IsPolyline == True:
        polyline_plot( ax, polyline_list, polyline_color )
    if addcolorbar != 'None':
        cb = fig.colorbar(dv,orientation=addcolorbar,ticks=[vmin,vmax])
        if addcolorbar.lower() == 'vertical':
            cb.ax.set_yticklabels( (vmin,vmax), size=6 )
        else:
            cb.ax.set_xticklabels( (vmin,vmax), size=6 )
    ax.set_xlim(x.min(),x.max())
    ax.set_ylim(y.min(),y.max())
    mpl.pyplot.xticks( size=6 )
    mpl.pyplot.yticks( size=6 )
    if text != 'None':
        ax.text(0.0, 1.01,text,size=8,horizontalalignment='left', verticalalignment='bottom', transform = ax.transAxes)
    fig.savefig(fout,dpi=150)
    del dv
    mpl.pyplot.close('all')
    gc.collect()
    return

def polyline_plot( ax, polyline_list, polyline_color='black', polyline_width=0.5, line_style='solid' ):
    for line in polyline_list:
        ax.plot(line[0,:],line[1,:],linestyle=line_style,color=polyline_color,linewidth=polyline_width)
    return

def point_plot( ax, point_list, marker='o', markersize=3, markeredgecolor='blue', markerfacecolor='None' ):
    for [x,y] in point_list:
        ax.plot(x, y, marker=marker, markersize=markersize, markeredgecolor=markeredgecolor,markerfacecolor=markerfacecolor)
    return
            