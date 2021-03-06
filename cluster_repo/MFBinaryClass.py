import sys
import numpy as np
import struct as strct
from pylab import ma, flipud
import string
import math

#generic functions
def kij_from_icrl(icrl,nlay,nrow,ncol):
    'Convert the modflow node number to row, column, and layer.'
    nrc = nrow * ncol
    #k=int( icrl / nrow / ncol )+1
    #i=int( (icrl-(k-1)*nrow*ncol) / ncol )+1
    #j=icrl - (k-1)*nrow*ncol - (i-1)*ncol
    k  = int( icrl / nrc )
    if ( k * nrc ) < icrl:
        k += 1
    ij = int( icrl - ( k - 1 ) * nrc )
    i  = int( ij / ncol )
    if ( i * ncol ) < ij:
        i += 1
    j = ij - ( i - 1 ) * ncol
    #print k, i, j
    return k,i,j

def icrl_from_kij(k,i,j,nlay,nrow,ncol):
    'Convert layer, row, and column to the modflow node number.'
    nrc = nrow * ncol
    icrl=int( ( ( k - 1 ) * nrc ) + ( ( i - 1 ) * ncol ) + j  )
    return icrl

def MFarray_to_plotarray(mfarray,maskvalue,orientation,rcl):
    '''Create a 2d plotting array from a 3d modflow array.
    mfarray: a 3d modflow array
    maskvalue: the value to mask (e.g. hdry)
    orientation: 'layer' 'row' or 'column'
    rcl: the layer row or column
    '''
    rcl=rcl-1
    nlay,nrow,ncol=shape(mfarray)
    if(orientation=='layer'):
        Z=flipud(mfarray[rcl,:,:]).copy()
    elif(orientation=='row'):
        Z=flipud(mfarray[:,rcl,:]).copy()
    elif(orientation=='column'):
        Z=flipud(mfarray[:,:,rcl]).copy()
    Z=ma.masked_where(Z == maskvalue,Z)
    return Z

class SWRReadBinaryStatements:
    integer = np.int32
    real = np.float64
    character = np.uint8
    integerbyte = 4
    realbyte = 8
    textbyte = 4
    def read_integer(self):
        intvalue=strct.unpack('i',self.file.read(1*SWRReadBinaryStatements.integerbyte))[0]
        return intvalue
    def read_real(self):
        realvalue=strct.unpack('d',self.file.read(1*SWRReadBinaryStatements.realbyte))[0]
        return realvalue
    def read_text(self):
        #textvalue=strct.unpack('cccccccccccccccc',self.file.read(16*self.textbyte))
        textvalue=np.fromfile(file = self.file, dtype=SWRReadBinaryStatements.character, count=16).tostring()
        if not isinstance(textvalue, str):
            textvalue = textvalue.decode().strip()
        else:
            textvalue = textvalue.strip()
        return textvalue
    def read_obs_text(self,nchar=20):
        #textvalue=strct.unpack('cccccccccccccccc',self.file.read(16*self.textbyte))
        textvalue=np.fromfile(file = self.file, dtype=MFReadBinaryStatements.character, count=nchar).tostring()
        if not isinstance(textvalue, str):
            textvalue = textvalue.decode().strip()
        else:
            textvalue = textvalue.strip()
        return textvalue
    def read_record(self):
#        x = np.fromfile(file=self.file,dtype=SWRReadBinaryStatements.real,count=self.nrecord*self.items)
#        x.resize(self.nrecord,self.items)
        if self.skip == True:
            lpos = self.file.tell() + ( self.nrecord*self.items*SWRReadBinaryStatements.realbyte )
            self.file.seek(lpos)
            x = np.zeros((self.nrecord*self.items),SWRReadBinaryStatements.real)
        else:
            x = np.fromfile(file=self.file,dtype=SWRReadBinaryStatements.real,count=self.nrecord*self.items)
        x.resize(self.nrecord,self.items)
        return x
    def read_items(self):
#        x = np.fromfile(file=self.file,dtype=SWRReadBinaryStatements.real,count=self.items)
#        return x
        if self.skip == True:
            lpos = self.file.tell() + ( self.items * SWRReadBinaryStatements.realbyte )
            self.file.seek(lpos)
            x = np.zeros((self.items),SWRReadBinaryStatements.real)
        else:
            x = np.fromfile(file=self.file,dtype=SWRReadBinaryStatements.real,count=self.items)
        return x
    def read_1dintegerarray(self):
        i = np.fromfile(file=self.file,dtype=SWRReadBinaryStatements.integer,count=self.nrecord)
        return i


class MFReadBinaryStatements:
    'Class of methods for reading MODFLOW binary files'
    #--byte definition
    integer=np.int32
    real=np.float32
    double=np.float64
    character=np.uint8
    integerbyte=4
    realbyte=4
    doublebyte=8
    textbyte=1
    def read_integer(self):
        intvalue=strct.unpack('i',self.file.read(1*MFReadBinaryStatements.integerbyte))[0]
        return intvalue
    def read_real(self):
        realvalue=strct.unpack('f',self.file.read(1*MFReadBinaryStatements.realbyte))[0]
        return realvalue
    def read_double(self):
        doublevalue=strct.unpack('f',self.file.read(1*MFReadBinaryStatements.doublebyte))[0]
        return doublevalue
    def read_text(self):
        #textvalue=strct.unpack('cccccccccccccccc',self.file.read(16*self.textbyte))
        textvalue=np.fromfile(file = self.file, dtype=MFReadBinaryStatements.character, count=16).tostring()
        if not isinstance(textvalue, str):
            textvalue = textvalue.decode().strip()
        else:
            textvalue = textvalue.strip()
        return textvalue
    def read_hyd_text(self,nchar=20):
        #textvalue=strct.unpack('cccccccccccccccc',self.file.read(16*self.textbyte))
        textvalue=np.fromfile(file = self.file, dtype=MFReadBinaryStatements.character, count=nchar).tostring()
        if not isinstance(textvalue, str):
            textvalue = textvalue.decode().strip()
        else:
            textvalue = textvalue.strip()
        return textvalue
    def read_3drealarray(self):
        x=np.fromfile(file = self.file, dtype=MFReadBinaryStatements.real, count=self.nlay*self.nrow*self.ncol)
        x.shape=(self.nlay,self.nrow,self.ncol)
        return x
    def skip_3drealarray(self):
        lpos = self.file.tell() + self.nlay*self.nrow*self.ncol*MFReadBinaryStatements.realbyte
        self.file.seek(lpos)
        return True
    def read_2drealarray(self):
        x=np.fromfile(file = self.file, dtype=MFReadBinaryStatements.real, count=self.nrow*self.ncol)
        x.shape=(self.nrow,self.ncol)
        return x
    def skip_2drealarray(self):
        lpos = self.file.tell() + self.nrow*self.ncol*MFReadBinaryStatements.realbyte
        self.file.seek(lpos)
        return True
    def read_2dintegerarray(self):
        i=np.fromfile(file = self.file, dtype=MFReadBinaryStatements.integer, count=self.nrow*self.ncol)
        i.shape=(self.nrow,self.ncol)
        return i
    def skip_2dintegerarray(self):
        lpos = self.file.tell() + self.nrow*self.ncol*MFReadBinaryStatements.integerbyte
        self.file.seek(lpos)
        return True
    def read_1drealarray(self,i):
        x=np.fromfile(file = self.file, dtype=MFReadBinaryStatements.real, count=i)
        return x

class MF_Discretization:
    def assign_rowcollay(self,nlay,nrow,ncol):
        #initialize grid information
        self.nrow=nrow
        self.ncol=ncol
        self.nlay=nlay
    #def read_PESTGridSpecificationFile(filename):
    #def read_MFDiscretizationFile(filename)

class SWR_BinaryObs(SWRReadBinaryStatements):
    'Reads binary head output from MODFLOW head file'
    def __init__(self,filename,verbose=False):
        #initialize class information
        self.skip = False
        self.verbose = verbose
        #--open binary head file
        self.file=open(filename,'rb')
        #--NOBS
        self.nobs=self.read_integer()
        self.v = np.empty((self.nobs),dtype='float')
        self.v.fill(1.0E+32)
        #--read obsnames
        obsnames = []
        for idx in range(0,self.nobs):
            cid = self.read_obs_text()
            obsnames.append( cid )
        self.obsnames = np.array( obsnames )
        #print self.obsnames
        #--set position
        self.datastart = self.file.tell()
        #get times
        self.times = self.time_list()

    def get_time_list(self):
        return self.times
        
    def get_num_items(self):
        return self.nobs

    def get_obs_labels(self):
        return self.obsnames
    
    def rewind_file(self):    
        self.file.seek(self.datastart)
        return True
        
    def time_list(self):    
        self.skip = True
        self.file.seek(self.datastart)
        times = []
        while True:
            current_position = self.file.tell()
            totim,v,success = self.next()
            if success == True:
                times.append([totim,current_position])
            else: 
                self.file.seek(self.datastart)
                times = np.array( times )
                self.skip = False
                return times

        
    def __iter__(self):
        return self

    def read_header(self):
        try:
            totim=self.read_real()
            return totim,True
        except:
            return -999.,False 

    def next(self):
        totim,success=self.read_header()
        if(success):
            for idx in range(0,self.nobs):
                self.v[idx] = self.read_real()
        else:
            if self.verbose == True:
                print('_BinaryObs object.next() reached end of file.')
            self.v.fill(1.0E+32)
        return totim,self.v,success

    def get_values(self,idx):
        iposition = int( self.times[idx,1] )
        self.file.seek(iposition)
        totim,v,success = self.next()
        if success == True:
            return totim,v,True
        else:
            self.v.fill( 1.0E+32 )
            return 0.0,self.v,False 

    def get_time_gage(self,record):
        idx = -1
        try:
            idx = int( record ) - 1
            if self.verbose == True:
                print('retrieving SWR observation record [{0}]'.format( idx+1 ))
        except:
            for icnt,cid in enumerate(self.obsnames):
                if record.strip().lower() == cid.strip().lower():
                    idx = icnt
                    if self.verbose == True:
                        print('retrieving SWR observation record [{0}] {1}'.format( idx+1, record.strip().lower() ))
                    break
        gage_record = np.zeros((2))#tottime plus observation
        if idx != -1 and idx < self.nobs:
            #--find offset to position
            ilen = self.get_point_offset(idx)
            #--get data
            for time_data in self.times:
                self.file.seek(int(time_data[1])+ilen)
                v=self.read_real()
                this_entry = np.array([float(time_data[0])])
                this_entry = np.hstack((this_entry,v))
                gage_record = np.vstack((gage_record,this_entry))
            #delete the first 'zeros' element
            gage_record = np.delete(gage_record,0,axis=0)
        return gage_record

    def get_point_offset(self,ipos):
        self.file.seek(self.datastart)
        lpos0 = self.file.tell()
        point_offset = int(0)
        totim,success=self.read_header()
        idx = (ipos)
        lpos1 = self.file.tell() + idx*SWRReadBinaryStatements.realbyte
        self.file.seek(lpos1)
        point_offset = self.file.tell() - lpos0
        return point_offset


class SWR_Record(SWRReadBinaryStatements):
    def __init__(self,swrtype,filename,verbose=False):
        #--swrtype =  0 = stage record
        #--swrtype = -1 = reach group record
        #--swrtype = -2 = reach group connection velocity record
        #--swrtype >  0 = aq-reach exchange record type = nlay
        self.file = open(filename,'rb')
        self.type = None
        try:
            ctype = swrtype.lower()
        except:
            ctype = None
            pass
        if ctype is not None:
            self.type = ctype
        else:
            try:
                itype = int(swrtype)
            except:
                print('SWR data type not defined')
                raise
            if itype == 0:
                self.type = 'stage'
            elif itype == -1:
                self.type = 'reachgroup'
            elif itype == -2:
                self.type = 'qm'
            elif itype > 0:
                self.type = 'qaq'
        if self.type is None:
            print('undefined SWR data type')
            raise
            
        self.verbose = verbose
        self.nrgout = 0
        if self.type == 'qm':
            self.nrgout = self.read_integer()
        self.nrecord = self.read_integer()
        self.items = self.get_num_items()
        self.null_record = np.zeros((self.nrecord,self.items)) + 1.0E+32
        #
        self.missingData = -9999.9
        self.dataAvailable = True
        self.skip = False
        #read connectivity for velocity data if necessary
        if self.type == 'qm':
            self.connectivity = self.read_connectivity()
            if self.verbose == True:
                print(self.connectivity)
        #--initialize reachlayers and nqaqentries for qaq data
        if self.type == 'qaq':
            self.reachlayers = np.zeros( (self.nrecord), np.int )
            self.nqaqentries = 0
        self.qaq_dtype = np.dtype([('layer','i4'),\
                                      ('bottom','f8'),('stage','f8'),\
                                      ('depth','f8'),('head','f8'),\
                                      ('wetper','f8'),('cond','f8'),\
                                      ('headdiff','f8'),('qaq','f8')])
        
        self.datastart = self.file.tell()
        #get times
        self.times = self.time_list()

    def get_nrecords(self):
        return self.nrgout, self.nrecord
        
    def get_time_list(self):
        return self.times
    
    def read_connectivity(self):
        conn = np.zeros( (self.nrecord,3), np.int )
        icount = 0
        for nrg in range(0,self.nrgout):
            nconn = self.read_integer()
            for ic in range(0,nconn):
                conn[icount,0] = nrg
                conn[icount,1] = self.read_integer()
                conn[icount,2] = self.read_integer()
                icount += 1
        return conn
    
    def get_connectivity(self):
        if self.type == 'qm':
            return self.connectivity
        else:
            return None
        
    def get_num_items(self):
        if self.type == 'stage': 
            return  1 
        elif self.type == 'reachgroup':
            return 14 
        elif self.type == 'qm':
            return  2 
        elif self.type == 'qaq':
            return  10 
        else:
            return -1
    
    def get_header_items(self):
        return ['totim','dt','kper','kstp','swrstp','success_flag'] 

    def get_item_list(self):
        if self.type == 'stage':
            list = ['stage']
        if self.type == 'reachgroup':
            list = ['stage','qsflow','qlatflow','quzflow','rain','evap',\
                            'qbflow','qeflow','qexflow','qbcflow','qcrflow','dv','inf-out','volume']
        if self.type == 'qm':
            list = ['flow','velocity']
        if self.type == 'qaq':
            list = ['reach','layer','bottom','stage','depth','head',\
                    'wetper','cond','headdiff','qaq']
        return list
         
    def get_temporal_list(self):
        list = ['totim','dt','kper','kstp','swrstp','success']
        return list

    def get_item_number(self,value,isTimeSeriesOutput=True):
        l = self.get_item_list()
        ioff = 6
        if isTimeSeriesOutput == False:
            ioff = 0
        try:
            i = l.index(value.lower())
            i += ioff
        except ValueError:
            l = self.get_temporal_list()
            try:
                i = l.index(value.lower())
            except ValueError:
                i = -1  #-no match
                print('no match to: ', value.lower())
        return i

    def return_gage_item_from_list(self,r,citem,scale=1.0):
        ipos = self.get_item_number(citem)
        n = r.shape[0]
        if n < 1:
            return self.null_record
        v = np.zeros( (n), np.float )
        for i  in range(0,n):
            v[i]  = r[i,ipos] * scale
        return v
        
    def read_header(self):
        if self.type == 'qaq':
            try:
                self.nqaqentries = 0
                for i in range(0,self.nrecord):
                    self.reachlayers[i] = self.read_integer() 
                    self.nqaqentries += self.reachlayers[i]
                    #print i+1, self.reachlayers[i]
                #print self.nqaqentries 
            except:
                if self.verbose == True:
                    sys.stdout.write('\nCould not read reachlayers')
                return 0.0,0.0,0,0,0,False
        try: 
            totim = self.read_real()
            dt = self.read_real()
            kper = self.read_integer()
            kstp = self.read_integer()
            swrstp = self.read_integer()
            return totim,dt,kper,kstp,swrstp,True
        except:
            return 0.0,0.0,0,0,0,False
    
    def get_record(self,*args):
        #--pass a tuple of timestep,stress period
        try:
            kkspt = args[0]
            kkper = args[1]
            while True:
                totim,dt,kper,kstp,swrstp,success,r = self.next()
                if success == True:
                    if kkspt == kstp and kkper == kper:
                        if self.verbose == True:
                            print(totim,dt,kper,kstp,swrstp,True)
                        return totim,dt,kper,kstp,swrstp,True,r
                else:
                    return 0.0,0.0,0,0,0,False,self.null_record
        except:
            #--pass a scalar of target totim - 
            #--returns either a match or the first
            #--record that exceeds target totim
            try:
                ttotim = float(args[0])
                while True:
                    totim,dt,kper,kstp,swrstp,r,success = self.next()
                    if success == True:
                        if ttotim <= totim:
                            return totim,dt,kper,kstp,swrstp,True,r
                    else:
                        return 0.0,0.0,0,0,0,False,self.null_record    
            except:
                #--get the last successful record
                previous = self.next()
                while True:
                    this_record = self.next()
                    if this_record[-2] == False:
                        return previous
                    else: previous = this_record
    
    def get_gage(self,rec_num=0,iconn=0,rec_lay=1):    
        if self.type == 'qaq':
            gage_record = np.zeros((self.items+6))#items plus 6 header values, reach number, and layer value
        else:
            gage_record = np.zeros((self.items+6))#items plus 6 header values
        while True:
            totim,dt,kper,kstp,swrstp,success,r = self.next()
            if success == True:
                this_entry = np.array([totim,dt,kper,kstp,swrstp,success])
                irec = rec_num - 1
                #find correct entry for record and layer
                if self.type == 'qaq':
                    ifound = 0
                    ilay = rec_lay
                    ilen = np.shape(r)[0]
                    #print np.shape(r)
                    for i in range(0,ilen):
                        ir = int(r[i,0])
                        il = int(r[i,1])
                        if ir == rec_num and il == ilay:
                            ifound = 1
                            irec = i
                            break
                    if ifound < 1:
                        r[irec,:] = 0.0
                elif self.type == 'qm':
                    ifound = 0
                    for i in range(0,self.nrecord):
                        inode = self.connectivity[i,1]
                        ic    = self.connectivity[i,2]
                        if rec_num == inode and ic == iconn:
                            ifound = 1
                            irec = i
                            break
                    if ifound < 1:
                        r[irec,:] = 0.0
                
                #if self.type == 'qaq':
                #    print 'shape gage_record {}'.format(gage_record.shape)
                #    print gage_record
                #    print 'shape this_entry {}'.format(this_entry.shape)
                #    print this_entry
                   
                this_entry = np.hstack((this_entry,r[irec]))
                gage_record = np.vstack((gage_record,this_entry))

            else: 
                gage_record = np.delete(gage_record,0,axis=0) #delete the first 'zeros' element
                return gage_record
                     
    def next(self):
        totim,dt,kper,kstp,swrstp,success = self.read_header()        
        if success == False: 
            if self.verbose == True:
                print('SWR_Record.next() object reached end of file')
            return 0.0,0.0,0,0,0,False,self.null_record
        else:
            if self.type == 'qaq':
                r = self.read_qaq()
                return totim,dt,kper,kstp,swrstp,True,r
            else:
                r = self.read_record()
        return totim,dt,kper,kstp,swrstp,True,r

    def read_qaq(self):
        x = np.zeros((self.nqaqentries,self.items), SWRReadBinaryStatements.real)                
        if self.skip == True:
            bytes = self.nqaqentries * (SWRReadBinaryStatements.integerbyte + 8*SWRReadBinaryStatements.realbyte)
            lpos = self.file.tell() + ( bytes )
            self.file.seek(lpos)
        else:
            qaq_list = self.get_item_list()
            bd = np.fromfile(self.file,dtype=self.qaq_dtype,count=self.nqaqentries)
            ientry = 0
            for irch in range(self.nrecord):
                klay = self.reachlayers[irch]
                for k in range(klay):
                    x[ientry,0] = irch+1
                    ientry += 1
            for idx, k in enumerate(qaq_list[1:]):
                x[:,idx+1] = bd[k]
        #print 'shape x: {}'.format(x.shape)
        return x
        
    
    def rewind_file(self):    
        self.file.seek(self.datastart)
        return True
        
    def time_list(self):    
        self.skip = True
        self.file.seek(self.datastart)
        idx = 0
        sys.stdout.write('Generating SWR binary data time list\n')
        times = []
        while True:
            #--output something to screen so it is possible to determine
            #  that the time list is being created
            idx += 1
            v = divmod( float(idx), 100. )
            if v[1] == 0.0:
                sys.stdout.write('.')
            #--get current position
            current_position = self.file.tell()
            totim,dt,kper,kstp,swrstp,success,r = self.next()
            if success == True:
                times.append([totim,dt,kper,kstp,swrstp,current_position])
            else: 
                self.file.seek(self.datastart)
                times = np.array( times )
                self.skip = False
                sys.stdout.write('\n')
                return times

    def get_time_record(self,time_index=0):
        self.file.seek(int(self.times[time_index][5]))
        totim,dt,kper,kstp,swrstp,success,r = self.next()
        if success == True:
            if self.verbose == True:
                print(totim,dt,kper,kstp,swrstp,True)
            return totim,dt,kper,kstp,swrstp,True,r
        else:
            return 0.0,0.0,0,0,0,False,self.null_record
    
    def get_point_offset(self,rec_num,iconn):
        self.file.seek(self.datastart)
        lpos0 = self.file.tell()
        point_offset = int(0)
        totim,dt,kper,kstp,swrstp,success = self.read_header()
        #--qaq terms
        if self.type == 'qaq':
            sys.stdout.write('MFBinaryClass::get_point_offset can not be used to extract QAQ data')
            sys.exit( 1 )
        #--stage and reach group terms
        elif self.type == 'stage' or self.type == 'reachgroup':
            idx = (rec_num-1)*self.items
            lpos1 = self.file.tell() + idx*SWRReadBinaryStatements.realbyte
            self.file.seek(lpos1)
            point_offset = self.file.tell() - lpos0
        #--connection flux and velocity terms
        elif self.type == 'qm':
            frec = -999
            for i in range(0,self.nrecord):
                inode = self.connectivity[i,1]
                ic    = self.connectivity[i,2]
                if rec_num == inode and ic == iconn:
                    frec = i
                    break
            if frec == -999:
                self.dataAvailable = False
            else:
                self.dataAvailable = True
                idx = (frec)*self.items
                lpos1 = self.file.tell() + idx*SWRReadBinaryStatements.realbyte
                self.file.seek(lpos1)
                point_offset = self.file.tell() - lpos0
        return point_offset

    def get_time_gage(self,rec_num=0,iconn=0):    
        if self.type == 'qaq':
            sys.stdout.write('MFBinaryClass::get_time_gage can not be used to extract QAQ data\n')
            sys.exit( 1 )
        num_records = self.items+6 #items plus 6 header values
        gage_record = np.zeros((num_records),np.float)
        #--find offset to position
        ilen = int(0)
        if rec_num > 0:
            ilen = self.get_point_offset(rec_num,iconn)
        else:
            self.dataAvailable = False
        if self.dataAvailable == False:
            sys.stdout.write('  Error: data is not available for reach {0} '.format(rec_num))
            if self.type == 'qm':
                sys.stdout.write('connected to reach {0}'.format(iconn))
            sys.stdout.write('\n')
        #--get data
        if len( self.times ) > 0:
            for time_data in self.times:
                totim = time_data[0]
                dt = time_data[1]
                kper = time_data[2]
                kstp = time_data[3]
                swrstp = time_data[4]
                success = True
                #--get data
                if self.dataAvailable == True:
                    self.file.seek(int(time_data[5])+ilen)
                    r = self.read_items()
                else:
                    r = np.empty((self.items),np.float)
                    r.fill(self.missingData)
                #--push the data to the data structure
                this_entry = np.array([totim,dt,kper,kstp,swrstp,success])
                #--update this_entry and current gage_record    
                this_entry = np.hstack((this_entry,r))
                gage_record = np.vstack((gage_record,this_entry))
        #--delete first empty entry and return gage_record
        gage_record = np.delete(gage_record,0,axis=0) #delete the first 'zeros' element
        return gage_record
       
        
    
class MODFLOW_Head(MFReadBinaryStatements,MF_Discretization):
    'Reads binary head output from MODFLOW head file'
    def __init__(self,nlay,nrow,ncol,filename,verbose=False):
        #initialize grid information
        self.assign_rowcollay(nlay,nrow,ncol)
        self.h = np.zeros((self.nlay, self.nrow, self.ncol)) + 1.0E+32
        self.items = self.get_num_items()
        self.x0 = 0.0
        self.y0 = 0.0
        self.dx = 0.0
        self.dy = 0.0
        self.verbose = verbose
        self.skip = False
        self.datastart = 0
        #open binary head file
        self.file=open(filename,'rb')
        #get times
        self.times = self.time_list()

    def get_time_list(self):
        return self.times
        
    def get_num_items(self):
        return 1 #heads
    
    def set_coordinates(self,x0,y0,dx,dy):
        self.x0 = x0
        self.y0 = y0
        self.dx = dx
        self.dy = dy
        #--set x and y coordinate for each row and column
        self.x = np.empty( (self.nrow+1,self.ncol+1) )
        self.y = np.empty( (self.nrow+1,self.ncol+1) )
        xt = np.zeros( self.ncol+1 )
        yt = np.zeros( self.nrow+1 )
        xt[0] = self.x0
        for j in range(1,self.ncol+1):
                xt[j] = xt[j-1] + self.dx
        yt[self.nrow] = self.y0
        for i in range(self.nrow-1,-1,-1):
                yt[i] = yt[i+1] + self.dy
                #print i, yt[i]
        for i in range(0,self.nrow+1):
                y = yt[i]
                for j in range(0,self.ncol+1):
                        x = xt[j]
                        self.x[i,j] = x
                        self.y[i,j] = y
        return 1
        
    def get_ijfromcoordinates(self,x,y):
        i = 0
        j = 0
        iexit = 0
        jexit = 0
        for ii in range(0,self.nrow):
                y1 = self.y[ii,0]
                y2 = self.y[ii+1,0]
                if y < y1 and y >= y2:
                        i = ii + 1
                        iexit = 1
                        for jj in range(0,self.ncol):
                                x1 = self.x[ii,jj]
                                x2 = self.x[ii,jj+1]
                                if x >= x1 and x < x2:
                                        j = jj + 1
                                        exit
                if iexit > 0:
                        exit

        return i,j 

    def get_nodefromrcl(self,row,col,lay):
        inode = icrl_from_kij(lay,row,col,self.nlay,self.nrow,self.ncol)
        return inode

    def read_header(self):
        try:
            kstp=self.read_integer()
            kper=self.read_integer()
            pertim=self.read_real()
            totim=self.read_real()
            text=self.read_text()
            ncol=self.read_integer()
            nrow=self.read_integer()
            ilay=self.read_integer()
            if self.verbose == True:
                print(kstp,kper,ilay,nrow,ncol,pertim,totim,True)
            return kstp,kper,pertim,totim,ncol,nrow,ilay,True
        except:
            return 0,0,0.,0.,0,0,0,False 
        
    def read_layerheads(self):
        if self.skip == True:
            success = self.skip_2drealarray()
            hl = np.zeros((self.nrow*self.ncol),np.float)
        else:
            hl = self.read_2drealarray()
        hl.shape=(self.nrow,self.ncol)
        return hl

    def __iter__(self):
        return self

    def next(self):
        for k in range(self.nlay):
            kstp,kper,pertim,totim,ncol,nrow,ilay,success=self.read_header()
            if(success):
                assert ncol==self.ncol, 'NCOL not consistent with binary heads file.'
                assert nrow==self.nrow, 'NROW not consistent with binary heads file.'
                assert ilay==k+1, 'Layers in head file are not sequential'
                self.h[ilay - 1, :, :] = self.read_layerheads()
            else:
                if self.verbose == True:
                    print('MODFLOW_Head object.next() reached end of file.')
                return 0.,0,0, np.zeros((self.nlay, self.nrow, self.ncol),\
                dtype='float')+1.0E+32,False
        self.KSTP=kstp
        self.KPER=kper
        self.PERTIM=pertim
        self.TOTIM=totim
        if self.verbose == True:
            print('Heads read for time step ',kstp,' and stress period ',kper)
        return totim,kstp,kper,self.h,True
     
    def get_record(self,*args):
        try:
            kkspt = args[0]
            kkper = args[1]
            while True:
                totim,kstp,kper,h,success = self.next()
                if success == True:
                    if kstp == kkspt and kkper == kper:
                        if self.verbose == True:
                            print(totim,kstp,kper,True)
                        return totim,kstp,kper,h,True
                else:
                    return 0.0,0,0,np.zeros((self.nlay,self.nrow,self,ncol),dtype='float')+1.0E+32,False 
        except:
            try:    
                target_totim = float(args[0])
                while True:
                    totim,kstp,kper,h,success = self.next()
                    if success:
                        if target_totim <= totim:
                            return totim,kstp,kper,h,True
                    else:
                        return 0.0,0,0,np.zeros((self.nlay,self.nrow,self.ncol),dtype='float')+1.0E+32,False
            
            except:
                #--get the last successful record
                previous = self.next()
                while True:
                    this_record = self.next()
                    if this_record[-1] == False:
                        return previous
                    else: previous = this_record
                                
    #rec_num is modflow node number
    def get_gage(self,rec_num):    
        k, i, j = kij_from_icrl(rec_num,self.nlay,self.nrow,self.ncol)
        if self.verbose == True:
            print('node=', rec_num, 'row=', i, ' col=', j, 'lay=', k)
        gage_record = np.zeros((self.items+1))#items plus tottime
        while True:
            totim,kstp,kper,h,success = self.next()
            if success == True:
                #print totim,np.shape(h[rec_num-1])
                this_entry = np.array([totim])
                this_entry = np.hstack((this_entry,h[k-1,i-1,j-1]))
                gage_record = np.vstack((gage_record,this_entry))
            else: 
                gage_record = np.delete(gage_record,0,axis=0) #delete the first 'zeros' element
                return gage_record

    def rewind_file(self):    
        self.file.seek(self.datastart)
        return True
        
    def time_list(self):    
        self.skip = True
        self.file.seek(self.datastart)
        idx = 0
        sys.stdout.write('Generating binary head data time list\n')
        times = []
        while True:
            #--output something to screen so it is possible to determine
            #  that the time list is being created
            idx += 1
            v = divmod( float(idx), 100. )
            if v[1] == 0.0:
                sys.stdout.write('.')
            #--get current file position
            current_position = self.file.tell()
            totim,kstp,kper,h,success = self.next()
            if success == True:
                times.append([totim,kstp,kper,current_position])
            else: 
                sys.stdout.write('\n')
                self.file.seek(self.datastart)
                times = np.array( times )
                self.skip = False
                return times

    def get_array(self,iposition):
        self.file.seek(iposition)
        totim,kstp,kper,h,success = self.next()
        if success == True:
            if self.verbose == True:
                print(totim,kstp,kper,True)
            return totim,kstp,kper,h,True
        else:
            return 0.0,0,0,np.zeros((self.nlay,self.nrow,self,ncol),dtype='float')+1.0E+32,False 

    def get_time_gage(self,rec_num):
        k, i, j = kij_from_icrl(rec_num,self.nlay,self.nrow,self.ncol)
        if self.verbose == True:
            print('node=', rec_num, 'row=', i, ' col=', j, 'lay=', k)
        gage_record = np.zeros((self.items+1))#items plus tottime
        #--find offset to position
        ilen = self.get_point_offset(k,i,j)
        #--get data
        for time_data in self.times:
            self.file.seek(int(time_data[3])+ilen)
            v=self.read_real()
            this_entry = np.array([float(time_data[0])])
            this_entry = np.hstack((this_entry,v))
            gage_record = np.vstack((gage_record,this_entry))
        #delete the first 'zeros' element
        gage_record = np.delete(gage_record,0,axis=0)
        return gage_record

    def get_point_offset(self,kpos,ipos,jpos):
        self.file.seek(self.datastart)
        lpos0 = self.file.tell()
        point_offset = int(0)
        for k in range(kpos):
            kstp,kper,pertim,totim,ncol,nrow,ilay,success=self.read_header()
            if(success):
                assert ncol==self.ncol, 'NCOL not consistent with binary heads file.'
                assert nrow==self.nrow, 'NROW not consistent with binary heads file.'
                assert ilay==k+1, 'Layers in head file are not sequential'
            if k < (kpos-1):
                success = self.skip_2drealarray()
            else:
                idx = (ipos-1)*self.ncol + (jpos-1)
                lpos1 = self.file.tell() + idx*MFReadBinaryStatements.realbyte
                self.file.seek(lpos1)
                point_offset = self.file.tell() - lpos0
        return point_offset



class MODFLOW_CBB(MFReadBinaryStatements,MF_Discretization):
    'Reads binary cell by cell output from MODFLOW cbb file'
    def __init__(self,nlay,nrow,ncol,filename,verbose=False):
        #initialize grid information
        self.assign_rowcollay(nlay,nrow,ncol)
        self.flux = np.empty((self.nlay, self.nrow, self.ncol))
        self.verbose = verbose
        #open binary head file
        self.file=open(filename,'rb')
        #set skip
        self.skip = False
        self.unique_items, self.all_times = self.get_all_times()
        sys.stdout.write( 'CBC data types\n' )
        for item in self.unique_items:
            sys.stdout.write( '{0} \n'.format( item ) )
        sys.stdout.write( '\n' )
        
    def get_time_list(self,fluxtype):
        self.times = self.time_list(fluxtype)
        return self.times

    def get_flux_item_list(self):
        return (self.unique_items)

    def read_header(self):
        try:
            kstp=self.read_integer()
            kper=self.read_integer()
            text=self.read_text()
            ncol=self.read_integer()
            nrow=self.read_integer()
            nlay=self.read_integer()
            ubdsvtype=0;delt=0.;pertim=0.;totim=0.
            if (nlay < 0):
                nlay=-nlay
                ubdsvtype = self.read_integer()
                delt = self.read_real()
                pertim = self.read_real()
                totim = self.read_real()
            if self.verbose == True:
                print(kstp,kper,text,nlay,nrow,ncol,ubdsvtype,delt,pertim,totim,True)
            return kstp,kper,text,nlay,nrow,ncol,ubdsvtype,delt,pertim,totim,True
        except:
#            return kstp,kper,text,nlay,nrow,ncol,ubdsvtype,delt,pertim,totim,False
            return 0,0,'',0,0,0,0,0.0,0.0,0.0,False

    def read_cbbdata(self,nlay,nrow,ncol,ubdsvtype,text):
        temp=np.zeros((nlay,nrow,ncol))
        if(ubdsvtype < 2):
            if self.skip == True:
                success = self.skip_3drealarray()
            else:
                temp[:,:,:]=self.read_3drealarray()
        elif(ubdsvtype == 2):
            nlist = self.read_integer()
            for i in range(nlist):
                icrl=self.read_integer()
                Q=self.read_real()
                k,i,j=kij_from_icrl(icrl,nlay,nrow,ncol)
                temp[k-1,i-1,j-1] += Q
        elif (ubdsvtype == 3):
            if self.skip == True:
                success = self.skip_2dintegerarray()
                success = self.skip_2drealarray()
            else:
                il = self.read_2dintegerarray()
                hl = self.read_2drealarray()
                for i in range(0,nrow):
                    for j in range(0,ncol):
                        k = il[i,j] - 1
                        temp[k,i,j] = hl[i,j]
        elif (ubdsvtype == 4):
            if self.skip == True:
                success = self.skip_2drealarray()
            else:
                temp[0,:,:] = self.read_2drealarray()
        elif (ubdsvtype == 5):
            naux = 1 - self.read_integer()
            if (naux > 0):
                for i in range(naux):
                    dummy=self.read_text()
            nlist = self.read_integer()
            for i in range(nlist):
                icrl=self.read_integer()
                Q=self.read_real()
                if self.skip == False:
                    k,i,j=kij_from_icrl(icrl,nlay,nrow,ncol)
                    temp[k-1,i-1,j-1] = temp[k-1,i-1,j-1] + Q
                if (naux > 0):
                    if self.skip == True:
                        for j in range(naux):
                            v = self.read_real()
                    else:
                        val = []
                        for j in range(naux):
                            val.append( self.read_real() )
        self.flux=temp
        return

    def next(self):
        kstp,kper,text,nlay,nrow,ncol,ubdsvtype,delt,pertim,totim,success=self.read_header()
        if(success):
            self.read_cbbdata(nlay,nrow,ncol,ubdsvtype,text)
            if self.verbose == True:
                print(text,totim,kstp,kper)
            return text,totim,kstp,kper,True
        else:
            if self.verbose == True:
                print('MODFLOW_CBB object.read_next_cbb() reached end of file.')
            return '',0.0,0,0,False

    def read_next_fluxtype(self,fluxtype):
        while(True):
#            text,totim,kstp,kper,success=self.read_next_cbb()
            text,totim,kstp,kper,success=self.next()
            #print text,totim,kstp,kper
            if (success):
                if (string.strip(string.ljust(text,16)) == string.strip(string.ljust(fluxtype,16))):
#               if (cmp(string.strip(string.ljust(text,16)),\
#               string.strip(string.ljust(fluxtype,16)))) == 0:
                    return self.flux,totim,True
            else:
                return np.empty((self.nlay,self.nrow,self.ncol)),0.,False
                
    def get_record(self,fluxtype,*args):
        while(True):
            text,totim,kstp,kper,success=self.next()
            if (success):
#               if (cmp(string.strip(string.ljust(text,16)),\
#               string.strip(string.ljust(fluxtype,16)))) == 0:
                if ( string.strip(string.ljust(text,16)) == string.strip(string.ljust(fluxtype,16)) ):
                    try:
                        kkstp = args[0]
                        kkper = args[1]
                        if (kstp == kkstp and kper == kkper):
                            return self.flux,totim,True
                    except:
                        try:
                            target_totim = float(args[0])
                            #dt = abs( totim - target_totim )
                            if target_totim <= totim:
                                return self.flux,totim,True
                        except:
                            return self.flux,totim,True
            else:
                return np.zeros((self.nlay,self.nrow,self.ncol),dtype='float')+1.0E+32,0.,False

    def print_cbb_info(self):
        success=True
        while(success):
            text,totim,success=self.read_next_cbb()
            print(text,'totim:',totim)
        return

    def get_all_times(self):
        self.skip = True
        self.file.seek(0)
        idx = 0
        sys.stdout.write('Generating MODFLOW CBC binary data time list\n')
        all_times = []
        while True:
            current_position = self.file.tell()
            text,totim,kstp,kper,success=self.next()
            if success == True:
                idx += 1
                v = divmod( float(idx), 100. )
                if v[1] == 0.0:
                    sys.stdout.write('.')
                all_times.append([text,totim,kstp,kper,current_position])
            else: 
                self.file.seek(0)
                sys.stdout.write('\n')
                unique_text = []
                for t in all_times:
                    if t[0] not in unique_text:
                        unique_text.append( t[0] )
                self.skip = False
                return (unique_text,all_times)

    def time_list(self,fluxtype=None):
        if fluxtype == None:
            fluxtype = self.all_times[0][0]
        if fluxtype not in self.unique_items:
            return np.array([0,0,0,0])
        times = []
        for [text,totim,kstp,kper,current_position] in self.all_times:
            if text == fluxtype:
                    times.append([totim,kstp,kper,current_position])
        times = np.array( times )
        return times
#           
#        self.skip = True
#        if fluxtype==None:
#            self.file.seek(0)
#            text,totim,kstp,kper,success=self.next()
#            fluxtype = text
#        self.file.seek(0)
#        times = []
#        while True:
#            current_position = self.file.tell()
#            text,totim,kstp,kper,success=self.next()
#            if success == True:
#                if text==fluxtype:
#                    times.append([totim,kstp,kper,current_position])
#            else: 
#                self.file.seek(0)
#                times = np.array( times )
#                self.skip = False
#                return times

    def get_array(self,iposition):
        self.file.seek(iposition)
        text,totim,kstp,kper,success=self.next()
        if success == True:
            return self.flux,totim,True
        else:
            return np.empty((self.nlay,self.nrow,self.ncol)),0.,False

    
class MT3D_Concentration(MFReadBinaryStatements,MF_Discretization):
    'Reads binary concentration output from MT3D concentration file'
    def __init__(self,nlay,nrow,ncol,filename,verbose=False):
        #initialize grid information
        self.assign_rowcollay(nlay,nrow,ncol)
        self.h = np.zeros((self.nlay, self.nrow, self.ncol),dtype='float')+1.0E+32
        self.verbose = verbose
        #open binary head file
        self.file=open(filename,'rb')
        self.times = self.time_list()
        
        self.totim = -999
        self.kper = -999
        self.kstp = -999
        self.ntrans = -999



    def get_time_list(self):
        return self.times


    def get_array(self,iposition):
        self.file.seek(iposition)
        totim,h,kstp,kper,success = self.next()
        if success == True:        
            return totim,kstp,kper,h,True
        else:
            return 0.0,0,0,np.zeros((self.nlay,self.nrow,self,ncol),dtype='float')+1.0E+32,False 


    def read_header(self):
        try:

            NTRANS=self.read_integer()
            KSTP=self.read_integer()
            KPER=self.read_integer()
            TOTIM=self.read_real()
            TEXT=self.read_text()
            NCOL=self.read_integer()
            NROW=self.read_integer()
            ILAY=self.read_integer()
            return NTRANS,KSTP,KPER,TOTIM,TEXT,NCOL,NROW,ILAY,True
        except:
            return -999,-999,-999,-999,-999,-999,-999,-999,False

    def read_layerconcens(self):
        cl = self.read_2drealarray()
        cl.shape=(self.nrow,self.ncol)
        return cl

    def __iter__(self):
        return self

    
    def time_list(self):    
        self.file.seek(0)
#        current_position = self.file.tell()
        idx = 0
        sys.stdout.write('Generating MT3D UCN binary data time list\n')
        times = []
        while True:
            #--output something to screen so it is possible to determine
            #  that the time list is being created
            idx += 1
            v = divmod( float(idx), 100. )
            if v[1] == 0.0:
                sys.stdout.write('.')
            #--get current file position
            current_position = self.file.tell()
            totim,h,kstp,kper,success = self.next()
            if success == True:
                #this_time = [totim,kstp,kper,current_position]
                times.append([totim,kstp,kper,current_position])
#                current_position = self.file.tell()
            else: 
                sys.stdout.write('\n')
                self.file.seek(0)
                times = np.array( times )
                return times

    def next(self):
     for k in range(self.nlay):
         NTRANS,KSTP,KPER,TOTIM,TEXT,NCOL,NROW,ILAY,success=self.read_header()
         if success:
             assert NCOL==self.ncol, 'NCOL not consistent with binary heads file.'
             assert NROW==self.nrow, 'NROW not consistent with binary heads file.'
             self.h[ILAY-1, :, :] = self.read_layerconcens()
         else:
             if self.verbose == True:
                print('MT3DMS_Concentration object.read_next_heads() reached end of file.')
             return 0., np.zeros((self.nlay, self.nrow, self.ncol),dtype='float')+1.0E+32, 0,0,False
     #print 'MT3DMS concentration read (ntrans,kstp,kper,time): ',NTRANS,KSTP,KPER,TOTIM 
     self.kper = KPER
     self.ntrans = NTRANS
     self.kstp = KSTP
     self.totim = TOTIM
     return self.totim,self.h,self.kstp,self.kper,True
     
    def get_record(self,*args):
        try:
            kkspt = args[0]
            kkper = args[1]
            while True:
                totim,concen,kstp,kper,success = self.next()
                if success:
                    if kstp == kkspt and kkper == kper:
                        return totim,kstp,kper,concen,True
                else:
                    return 0.0,0,0,np.zeros((self.nlay,self.nrow,self,ncol),dtype='float')+1.0E+32,False 
        except:
            target_totim = args[0]
            while True:
                totim,concen,kstp,kper,success = self.next()
                if success:
                    if target_totim == totim:
                        return totim,kstp,kper,concen,True
                else:
                    return 0.0,0,0,np.zeros((self.nlay,self.nrow,self.ncol),dtype='float')+1.0E+32,False
                

class MODFLOW_HYDMOD(MFReadBinaryStatements):
    'Reads binary head output from MODFLOW head file'
    def __init__(self,filename,double=False,slurp=False,verbose=False):
        '''slurp is a short cut to read all output using numpy fromfile()
        if you use it, you don't need to read times
        '''
        #initialize class information
        self.skip = True
        self.double = bool(double)
        self.verbose = verbose
        #--open binary head file
        self.file=open(filename,'rb')
        #--NHYDTOT,ITMUNI
        self.nhydtot=self.read_integer()
        self.itmuni=self.read_integer()
        if self.nhydtot < 0:
            self.double = True
            self.nhydtot = abs( self.nhydtot )
        self.v = np.empty((self.nhydtot),dtype='float')
        self.v.fill(1.0E+32)
        ctime = self.read_hyd_text(nchar=4)
        #--read HYDLBL
        hydlbl = []
        #hydid = []
        for idx in range(0,self.nhydtot):
            cid = self.read_hyd_text()
            hydlbl.append( cid )
        self.hydlbl = np.array( hydlbl )
        if self.verbose == True:
            print(self.hydlbl)
        if not slurp:
            #--set position
            self.datastart = self.file.tell()
            #get times
            self.times = self.time_list()


    def get_time_list(self):
        return self.times
        
    def get_num_items(self):
        return self.nhydtot

    def get_hyd_labels(self):
        return self.hydlbl
    
    def rewind_file(self):    
        self.file.seek(self.datastart)
        return True
        
    def time_list(self):    
        self.skip = True
        self.file.seek(self.datastart)
        times = []
        while True:
            current_position = self.file.tell()
            totim,v,success = self.next()
            if success == True:
                times.append([totim,current_position])
            else: 
                self.file.seek(self.datastart)
                times = np.array( times )
                self.skip = False
                return times

        
    def __iter__(self):
        return self

    def slurp(self):
        if self.double:
            float_type = np.float64
        else:
            float_type = np.float32
        dtype_list = [('totim',float_type)]
        for site in self.hydlbl:
            dtype_list.append((site[6:].strip(),float_type))
        dtype = np.dtype(dtype_list)
        data = np.fromfile(self.file,dtype,count=-1)
        return data        

    def read_header(self):
        try:
            totim=self.read_real()
            return totim,True
        except:
            return -999.,False 

    def next(self):
        totim,success=self.read_header()
        if(success):
            for idx in range(0,self.nhydtot):
                if self.double==True:
                    self.v[idx] = float(self.read_double())
                else:
                    self.v[idx] = self.read_real()
        else:
            if self.verbose == True:
                print('MODFLOW_HYDMOD object.next() reached end of file.')
            self.v.fill(1.0E+32)
        return totim,self.v,success

    def get_values(self,idx):
        iposition = int( self.times[idx,1] )
        self.file.seek(iposition)
        totim,v,success = self.next()
        if success == True:
            return totim,v,True
        else:
            self.v.fill( 1.0E+32 )
            return 0.0,self.v,False 

    def get_time_gage(self,record,lblstrip=6):
        idx = -1
        try:
            idx = int( record ) - 1
            if idx >= 0 and idx < self.nhydtot:
                if self.verbose == True:
                    print('retrieving HYDMOD observation record [{0}]'.format( idx+1 ))
            else:
                print('Error: HYDMOD observation record {0} not found'.format( record.strip().lower() ))
        except:
            for icnt,cid in enumerate(self.hydlbl):
                if lblstrip > 0:
                    tcid = cid[lblstrip:len(cid)]
                else:
                    tcid = cid
                if record.strip().lower() == tcid.strip().lower():
                    idx = icnt
                    if self.verbose == True:
                        print('retrieving HYDMOD observation record [{0}] {1}'.format( idx+1, record.strip().lower() ))
                    break
            if idx == -1:
                print('Error: HYDMOD observation record {0} not found'.format( record.strip().lower() ))
        gage_record = np.zeros((2))#tottime plus observation
        if idx != -1 and idx < self.nhydtot:
            #--find offset to position
            ilen = self.get_point_offset(idx)
            #--get data
            for time_data in self.times:
                self.file.seek(int(time_data[1])+ilen)
                if self.double == True:
                    v=float(self.read_double())
                else:
                    v=self.read_real()
                this_entry = np.array([float(time_data[0])])
                this_entry = np.hstack((this_entry,v))
                gage_record = np.vstack((gage_record,this_entry))
            #delete the first 'zeros' element
            gage_record = np.delete(gage_record,0,axis=0)
        return gage_record

    def get_point_offset(self,ipos):
        self.file.seek(self.datastart)
        lpos0 = self.file.tell()
        point_offset = int(0)
        totim,success=self.read_header()
        idx = (ipos)
        if self.double == True:
            lpos1 = self.file.tell() + idx*MFReadBinaryStatements.doublebyte
        else:
            lpos1 = self.file.tell() + idx*MFReadBinaryStatements.realbyte
        self.file.seek(lpos1)
        point_offset = self.file.tell() - lpos0
        return point_offset

                    