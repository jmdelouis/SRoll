#!PYPATH

# ==============================================================================
# 
#                      iSROLL : INVERT DATA WITH IA
# 
# ==============================================================================

"""iSROLL : INVERT DATA WITH IA

"""


import time
import ctypes
import numpy as np
import Csroll as Cs

import argparse
import os
import sys
import time

from mpi4py import MPI
import psutil
from inspect import currentframe, getframeinfo  
from scipy import interpolate

def getloadavg():
    return(open("/proc/loadavg").readline().split(" ")[:3])

#==============================================================================
#from mpi4py import MPI
#==============================================================================

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

param_file=sys.argv[1]
exec(open(param_file).read())
if rank==0:
    print('TFLEARN',TFLEARN)
#==============================================================================
#
#               PARAMETER BY DEFAULT
#
#==============================================================================

WSCALE=4

if (len(sys.argv)>2):
    for k in range(len(sys.argv)-2):
        if sys.argv[2+k]=='--dotemplate':
            dotemplate=True
        if sys.argv[2+k]=='--doplot':
            DOPLOT=True
        if sys.argv[2+k]=='--doreference':
            DOREFERENCE=True
        if sys.argv[2+k]=='--onlyparam':
            PARAMETERONLY=True
        if sys.argv[2+k]=='--fullyconnected':
            DOFULLYCONNECTED=True
try:
    DOREFERENCE
except NameError:
    DOREFERENCE=False

try:
    TFLEARN
except NameError:
    TFLEARN=Out_Offset[0]
try:
    PARAMETERONLY
except NameError:
    PARAMETERONLY=False

if PARAMETERONLY==True:
    if rank==0:
        print('====================================================================')
        print('==                       PARAMETER ONLY                           ==')
        print('====================================================================')
try:
    dotemplate
except NameError:
    dotemplate=False

try:
    showdata
except NameError:
    showdata=False
try:
    DOPLOT
except NameError:
    DOPLOT=False
    
try:
    NUM_EPOCHS
except NameError:
    NUM_EPOCHS = 10000
try:
    EVAL_FREQUENCY
except NameError:
    EVAL_FREQUENCY = 500
try:
    NITT_CNN
except NameError:
    NITT_CNN = 200
try:
    PRT_FREQUENCY
except NameError:
    PRT_FREQUENCY = 10
try:
    DECAY_RATE
except NameError:
    DECAY_RATE = 1.0
try:
    LEARNING_RATE
except NameError:
    LEARNING_RATE = 0.01
try:
    NITT
except NameError:
    NITT = 10
try:
    SEED
except NameError:
    SEED = 1234
try:
    BATCH_SIZE
except NameError:
    BATCH_SIZE=-1
try:
    Nside_cnn
except NameError:
    Nside_cnn=Nside
try:
    verbose
except NameError:
    verbose = 0
try:
    WSCALE
except NameError:
    WSCALE = 4

#==============================================================================
#
#               DO THE JOB !!!
#
#==============================================================================

if DOPLOT==True:
    import matplotlib.pyplot as plt
    
start_time_total = time.time()
#==============================================================================
#
#               Tool fo neural network synthesis
#
#==============================================================================

    
def gettype(name):
    if name=='FLOAT32':
        return('float32')
    if name=='FLOAT64':
        return('float64')
    if name=='INT64':
        return('int64')
    if name=='INT32':
        return('int32')
    if name=='INT16':
        return('int16')
    return('NoType')

def getmap(name,nside,dtype='float32',rgsz=27664,step=1):

    lname=name.split(':')
    if len(lname)<3:
        print('problem to read ',name)
    ldtype=gettype(lname[1])
    if lname[0]=='BINARY':
        bytsz=sys.getsizeof(np.zeros([2],dtype=ldtype))-sys.getsizeof(np.zeros([1],dtype=ldtype))

        if showdata==True:
            print('SHOWDATA %s'%(lname[2]))
            return(np.ones([12*nside*nside],dtype=dtype))

        f = open(lname[2], "rb")
        try:
            out = (np.frombuffer(f.read(12*nside*nside*bytsz),dtype=ldtype).reshape(12*nside*nside)).astype(dtype)
        finally:
            f.close()
        return(out)
    else:
        print('Unknown type ',lname[0])
def getdata(name,beginrg,endrg,dtype='float32',rgsz=27664,step=1):
    
    lname=name.split(':')
    if len(lname)<3:
        print('problem to read ',name)
    ldtype=gettype(lname[1])
    if lname[0]=='BINARY':
        bytsz=sys.getsizeof(np.zeros([2],dtype=ldtype))-sys.getsizeof(np.zeros([1],dtype=ldtype))

        if showdata==True:
            print('SHOWDATA %s'%(lname[2]))
            return(np.ones([(endrg-beginrg+1)//step,rgsz],dtype=dtype))

        f = open(lname[2], "rb")
        if step==1:
            try:
                f.seek(beginrg*rgsz*bytsz)
                out = (np.frombuffer(f.read((endrg-beginrg+1)*rgsz*bytsz),dtype=ldtype).reshape(endrg-beginrg+1,rgsz)).astype(dtype)
            finally:
                f.close()
        else:
            nring=(endrg-beginrg+1)//step
            out=np.zeros([nring,rgsz],dtype=dtype)
            try:
                for ir in range(nring):
                    f.seek((beginrg+ir*step)*rgsz*bytsz)
                    out[ir,:] = np.frombuffer(f.read((rgsz*bytsz)),dtype=ldtype).astype(dtype)
            finally:
                f.close()
        return(out)

def main():
    
    #=============================================================================================================================
    #                                    INIT MPI
    #=============================================================================================================================
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    NUMBEROFTHREAD=NUMBER_TF
    
    if NUMBEROFTHREAD>size:
        NUMBEROFTHREAD=size
        
    if NUMBEROFTHREAD!=1:
        color = int(rank/NUMBEROFTHREAD)
        key  = int(rank%NUMBEROFTHREAD)

        ncomm = comm.Split(color,key)
        ncomm_transpose = comm.Split(key,color)
    else:
        ncomm = comm
        ncomm_transpose = comm

    if rank%NUMBEROFTHREAD==0:
        IMPORTTF1
        IMPORTTF2

    Cfunc = Cs.Csroll("S4PATH")

    def ChangeReso(idx,nsidein,nsideout):
        if nsidein==nsideout:
            return(idx)
        lndata=len(idx.flatten())
        th=np.zeros([lndata],dtype='float')
        ph=np.zeros([lndata],dtype='float')
        nhidx=np.zeros([lndata],dtype='int64')
        Cfunc.pix2ang(nsidein,(idx.flatten()).astype('int64'),th,ph,lndata)
        Cfunc.ang2pix(nsideout,th,ph,nhidx,lndata)
        del th
        del ph
        return(nhidx.astype('int32'))
        
    #=============================================================================================================================
    #                                    COMPUTE BEGIN/ENDRG MPI
    #=============================================================================================================================
    NringPerProc = ((EndRing - BeginRing+1)//size)
    Nring2share = (EndRing - BeginRing+1)-NringPerProc*size
    RingTab = np.zeros([size],dtype='int32')+NringPerProc
    RingTab[0:Nring2share]+=1

    BeginRingTab = np.zeros([size],dtype='int32')
    EndRingTab = np.zeros([size],dtype='int32')
    BeginRingTab[0]=BeginRing
    for i in range(size-1):
        BeginRingTab[i+1]=BeginRingTab[i]+RingTab[i]
    EndRingTab[0:size-1]=BeginRingTab[1:size]-1
    EndRingTab[size-1]=EndRing
    if (rank==0):
        print('BEGINRING = ',BeginRingTab)
        print('ENDRING = ',EndRingTab)
    
    #=============================================================================================================================
    #                                    READ AND EXCHANGE DATA
    #=============================================================================================================================
    nbolo=len(Calibration)
    data={}
    tmp_data={}
    tmp_func={}
    func_input={}
    add_data={}
    calib={}
    hdata={}
    hdata2={}
    hidx={}
    ridx={}
    val_ring={}
    nring={}
    coefnorm={}
    cnn_idx={}
    cnn_imref={}
    rotidx={}
    protidx={}
    maskrot={}
    map_fring={}
    
    icalib=-np.ones([nbolo],dtype='int32')
    
    if (rank==0):
        print('Nbolo ',nbolo)

    normrelative=False
    docalib=False
    ncnn=0
    DOMAIN={}
    XIMAGE_SIZE    = {}
    YIMAGE_SIZE    = {}
    KERNELSZ       = {}
    NHIDDEN        = {}
    NUM_DCONV_CHAN = {}
    NDCONV         = {}
    YCNN           = {}
    CIRCULAR       = {}
    is2D           = {}
    isSCAT         = {}
    isFSL          = {}
    ty             = {}
    tico           = {}
    cnn_ib         = {}
    cnn_ibolo      = {}
    REGULAMP       = {}
    VAR2FIT        = {}
    AMPDECOD       = {}
    INITDATA       = {}
    CNNTEMPLATE    = {}
    ROTIDX         = {}
    IMRIDX         = {}
    NSIDEFSL       = {}
    
    COMPUTE_SCATTER_OBJECTIVE   = {}
    FILTER_SCATTER              = {}
    SCATTER_CLEAN               = {}
    WSCATT                      = {}
    NSCATTER                    = {}
    NHARM_SCATTER               = {}
    NHARM_LARGE_SCALE           = {}
    SCATTER_OBJECTIVE_FILE      = {}
    SCATTER_0_MEAN_FILE         = {}
    SCATTER_1_MEAN_FILE         = {}
    SCATTER_2_MEAN_FILE         = {}
    SCATTER_0_STD_FILE          = {}
    SCATTER_1_STD_FILE          = {}
    SCATTER_2_STD_FILE          = {}
    SCATTER_NORM                = {}
    SCATTER_WALL                = {}
    NHARM                       = {}
    
    #=============================================================================================================================
    #                                    MANAGE MASK
    #=============================================================================================================================
    if DOREFERENCE==True:
        mask=np.ones([12*Nside*Nside],dtype='float32')
    else:
        mask=getmap(Mask,Nside,dtype='float32')
    
    
    cnn_ibolo=np.zeros([nbolo+1],dtype='int32')
    
    for ib in range(nbolo):
        
        cnn_ibolo[ib]=ncnn
        #=============================================================================================================================
        #                                    READ DATA
        #=============================================================================================================================
        data[ib]     = (getdata(Signal[ib],BeginRingTab[rank],EndRingTab[rank],step=RSTEP)-Monop[ib])/Calibration[ib]
        if rank==0:
            print('DATA ',Signal[ib])
        hdata[ib]    = getdata(Hit[ib],BeginRingTab[rank],EndRingTab[rank],step=RSTEP)*(Calibration[ib]/NEP[ib])**2
        if rank==0:
            print('HDATA ',Hit[ib])
        ph           = getdata(Ptg_PH[ib],BeginRingTab[rank],EndRingTab[rank],dtype='float64',step=RSTEP)
        if rank==0:
            print('PH PTG ',Ptg_PH[ib])
        th           = getdata(Ptg_TH[ib],BeginRingTab[rank],EndRingTab[rank],dtype='float64',step=RSTEP)
        if rank==0:
            print('PH PTG ',Ptg_TH[ib])
        psi          = getdata(Ptg_PSI[ib],BeginRingTab[rank],EndRingTab[rank],dtype='float64',step=RSTEP)
        if rank==0:
            print('PH PSI ',Ptg_PSI[ib])
        pidx=np.zeros(len(th.flatten()),dtype='int64')
        Cfunc.ang2pix(Nside,th.flatten(),ph.flatten(),pidx,len(th.flatten()))
        hidx[ib]     = pidx.reshape(th.shape)
        hdata2[ib]   = hdata[ib]*(mask[hidx[ib]]+1E-2)

        f_bad=getdata(Badring[ib],0,EndRing,rgsz=1,dtype='int32') #np.fromfile(Badring[ib],dtype='int32')
        nx,ny=data[ib].shape
        ridx[ib]=(BeginRingTab[rank]+np.arange(nx)*RSTEP)
        idx=np.where(f_bad[ridx[ib]]==0)
            
                                     
        #=============================================================================================================================
        #                                    MANAGE LINEAR TEMPLATE VALUE
        #=============================================================================================================================
        ntmp=0
        nfunc=0
        ntotfunc=0
        ntotMfunc=0
        ninput_template=len(Template[ib])
        for itmp in range(ninput_template):
            if Template[ib][itmp][3]!='RELATIV' or nbolo>1:
                if Template[ib][itmp][1]=='CALIB':
                    icalib[ib]=ntmp
                    docalib=True
                if Template[ib][itmp][1]=='LINEAR' or Template[ib][itmp][0]=='CALIB':
                    ntmp=ntmp+1
                if Template[ib][itmp][1]=='SPLINE1D':
                    nfunc=nfunc+1

                if Template[ib][itmp][1]=='CNN1D':
                    DOMAIN[ncnn]         = Template[ib][itmp][4]
                    YCNN[ncnn]           = Template[ib][itmp][5]
                    NHIDDEN[ncnn]        = int(Template[ib][itmp][6])
                    NUM_DCONV_CHAN[ncnn] = Template[ib][itmp][7]
                    KERNELSZ[ncnn]       = int(Template[ib][itmp][8])
                    CIRCULAR[ncnn]       = Template[ib][itmp][9]
                    XIMAGE_SIZE[ncnn]    = Template[ib][itmp][10]
                    NHARM[ncnn]          = Template[ib][itmp][11]
                    REGULAMP[ncnn]       = Template[ib][itmp][12]
                    VAR2FIT[ncnn]        = 1.0
                    AMPDECOD[ncnn]       = 1.0
                    NDCONV[ncnn]         = len(NUM_DCONV_CHAN[ncnn])-1
                    INITDATA[ncnn]       = 'NA'
                    CNNTEMPLATE[ncnn]    = -1
                    if rank==0:
                        print(NUM_DCONV_CHAN[ncnn])
                        print(NDCONV[ncnn])
                    isSCAT[ncnn]=False
                    is2D[ncnn]=False
                    isFSL[ncnn]=False
                    cnn_ib[ncnn]=ib
                    ncnn+=1
                
                if Template[ib][itmp][1]=='CNN1D_SCAT':
                    DOMAIN[ncnn]         = Template[ib][itmp][4]
                    YCNN[ncnn]           = Template[ib][itmp][5]
                    NHIDDEN[ncnn]        = int(Template[ib][itmp][6])
                    NUM_DCONV_CHAN[ncnn] = Template[ib][itmp][7]
                    KERNELSZ[ncnn]       = int(Template[ib][itmp][8])
                    CIRCULAR[ncnn]       = Template[ib][itmp][9]
                    XIMAGE_SIZE[ncnn]    = Template[ib][itmp][10]
                    NDCONV[ncnn]         = len(NUM_DCONV_CHAN[ncnn])-1
                    COMPUTE_SCATTER_OBJECTIVE[ncnn]     = Template[ib][itmp][11]
                    FILTER_SCATTER[ncnn]                = Template[ib][itmp][12]
                    SCATTER_CLEAN[ncnn]                 = Template[ib][itmp][13]
                    NSCATTER[ncnn]                      = Template[ib][itmp][14]
                    NHARM_SCATTER[ncnn]                 = 2*Template[ib][itmp][15]
                    NHARM_LARGE_SCALE[ncnn]             = Template[ib][itmp][16]
                    WSCATT[ncnn]                        = Template[ib][itmp][17]
                    NHARM[ncnn]          = 0
                    if WSCATT[ncnn]==0:
                        SCATTER_NORM[ncnn]=False
                    else:
                        SCATTER_NORM[ncnn]=True
                    SCATTER_WALL[ncnn]                  = Template[ib][itmp][18]
                    SCATTER_OBJECTIVE_FILE[ncnn]        = Template[ib][itmp][19]
                    SCATTER_0_MEAN_FILE[ncnn]           = Template[ib][itmp][20]
                    SCATTER_1_MEAN_FILE[ncnn]           = Template[ib][itmp][21]
                    SCATTER_2_MEAN_FILE[ncnn]           = Template[ib][itmp][22]
                    SCATTER_0_STD_FILE[ncnn]            = Template[ib][itmp][23]
                    SCATTER_1_STD_FILE[ncnn]            = Template[ib][itmp][24]
                    SCATTER_2_STD_FILE[ncnn]            = Template[ib][itmp][25]
                    REGULAMP[ncnn]       = Template[ib][itmp][26]
                    VAR2FIT[ncnn]        = 1.0
                    AMPDECOD[ncnn]       = 1.0
                    if rank==0:
                        print("================= CNN1D_SCAT ====================")
                        print("NUM_DCONV_CHAN "        ,NUM_DCONV_CHAN[ncnn])
                        print("NDCONV "                ,NDCONV[ncnn])
                        print("NHARM_SCATTER "         ,NHARM_SCATTER[ncnn])
                        print("SCATTER_OBJECTIVE_FILE ",SCATTER_OBJECTIVE_FILE[ncnn])
                        print("WSCATT "                ,WSCATT[ncnn])
                        print("SCATTER_0_MEAN "        ,SCATTER_0_MEAN_FILE[ncnn])
                        print("SCATTER_1_MEAN "        ,SCATTER_1_MEAN_FILE[ncnn])
                        print("SCATTER_2_MEAN "        ,SCATTER_2_MEAN_FILE[ncnn])
                        print("SCATTER_0_STD_FILE "    ,SCATTER_0_STD_FILE[ncnn])
                        print("SCATTER_1_STD_FILE "    ,SCATTER_1_STD_FILE[ncnn])
                        print("SCATTER_2_STD_FILE "    ,SCATTER_2_STD_FILE[ncnn])
                        print("=================================================")
                    is2D[ncnn]=False
                    if dotemplate==False:
                        isSCAT[ncnn]=True
                    else:
                        isSCAT[ncnn]=False
                    isFSL[ncnn]=False
                    cnn_ib[ncnn]=ib
                    INITDATA[ncnn]       = 'NA'
                    CNNTEMPLATE[ncnn]    = -1
                    ncnn+=1
                
                if Template[ib][itmp][1]=='CNN2D':
                    DOMAIN[ncnn]         = Template[ib][itmp][4]
                    YCNN[ncnn]           = Template[ib][itmp][5]
                    NHIDDEN[ncnn]        = int(Template[ib][itmp][6])
                    NUM_DCONV_CHAN[ncnn] = Template[ib][itmp][7]
                    KERNELSZ[ncnn]       = int(Template[ib][itmp][8])
                    CIRCULAR[ncnn]       = Template[ib][itmp][9]
                    XIMAGE_SIZE[ncnn]    = Template[ib][itmp][10]
                    YIMAGE_SIZE[ncnn]    = Template[ib][itmp][11]
                    NDCONV[ncnn]         = len(NUM_DCONV_CHAN[ncnn])-1
                    NHARM[ncnn]          = Template[ib][itmp][12]
                    REGULAMP[ncnn]       = Template[ib][itmp][13]
                    VAR2FIT[ncnn]        = 1.0
                    AMPDECOD[ncnn]       = 1E-3
                    
                    if rank==0:
                        print(NUM_DCONV_CHAN[ncnn])
                        print(NDCONV[ncnn])
                    is2D[ncnn]=True
                    isSCAT[ncnn]=False
                    isFSL[ncnn]=False
                    cnn_ib[ncnn]=ib
                    INITDATA[ntmp]      = 'NA'
                    CNNTEMPLATE[ncnn]   = ntmp
                    ntmp=ntmp+1
                    ncnn+=1
                if Template[ib][itmp][1]=='CNN2D_FSL':
                    DOMAIN[ncnn]         = Template[ib][itmp][4]
                    YCNN[ncnn]           = Template[ib][itmp][5]
                    NHIDDEN[ncnn]        = int(Template[ib][itmp][6])
                    NUM_DCONV_CHAN[ncnn] = Template[ib][itmp][7]
                    KERNELSZ[ncnn]       = int(Template[ib][itmp][8])
                    CIRCULAR[ncnn]       = Template[ib][itmp][9]
                    XIMAGE_SIZE[ncnn]    = Template[ib][itmp][10]
                    YIMAGE_SIZE[ncnn]    = Template[ib][itmp][11]
                    NDCONV[ncnn]         = len(NUM_DCONV_CHAN[ncnn])-1
                    NHARM[ncnn]          = Template[ib][itmp][12]
                    REGULAMP[ncnn]       = Template[ib][itmp][13]
                    VAR2FIT[ncnn]        = Template[ib][itmp][18]
                    ROTIDX[ncnn]         = {}
                    ROTIDX[ncnn]['data'] = (np.load(Template[ib][itmp][15]+'.npy')).transpose()
                    ROTIDX[ncnn]['pidx'] = (np.load(Template[ib][itmp][15]+'_pidx.npy')).transpose()
                    ROTIDX[ncnn]['mask'] = (np.load(Template[ib][itmp][15]+'_mask.npy'))
                    IMRIDX[ncnn]         = Template[ib][itmp][16]
                    NSIDEFSL[ncnn]       = int(Template[ib][itmp][17])
                    if rank==0:
                        print(NUM_DCONV_CHAN[ncnn])
                        print(NDCONV[ncnn])
                    is2D[ncnn]=True
                    isSCAT[ncnn]=False
                    isFSL[ncnn]=True
                    cnn_ib[ncnn]=ib
                    INITDATA[ntmp]      = Template[ib][itmp][14]
                    CNNTEMPLATE[ncnn]   = -1
                    AMPDECOD[ncnn]       = 1.0
                    ncnn+=1
                if Template[ib][itmp][1]=='CNN1D_FSL':
                    DOMAIN[ncnn]         = Template[ib][itmp][4]
                    YCNN[ncnn]           = Template[ib][itmp][5]
                    NHIDDEN[ncnn]        = int(Template[ib][itmp][6])
                    NUM_DCONV_CHAN[ncnn] = Template[ib][itmp][7]
                    KERNELSZ[ncnn]       = int(Template[ib][itmp][8])
                    CIRCULAR[ncnn]       = Template[ib][itmp][9]
                    XIMAGE_SIZE[ncnn]    = int(Template[ib][itmp][10])
                    NHARM[ncnn]          = Template[ib][itmp][11]
                    REGULAMP[ncnn]       = Template[ib][itmp][12]
                    VAR2FIT[ncnn]        = Template[ib][itmp][16]
                    AMPDECOD[ncnn]       = 1.0
                    ROTIDX[ncnn]         = {}
                    ROTIDX[ncnn]['data'] = (np.load(Template[ib][itmp][13]+'.npy')).transpose()
                    #ROTIDX[ncnn]['pidx'] = (np.load(Template[ib][itmp][13]+'_pidx.npy')).transpose()
                    #ROTIDX[ncnn]['mask'] = (np.load(Template[ib][itmp][13]+'_mask.npy'))
                    IMRIDX[ncnn]         = Template[ib][itmp][14]
                    NSIDEFSL[ncnn]       = int(Template[ib][itmp][15])
                    NDCONV[ncnn]         = len(NUM_DCONV_CHAN[ncnn])-1
                    INITDATA[ncnn]       = 'NA'
                    CNNTEMPLATE[ncnn]    = ntmp
                    if rank==0:
                        print(NUM_DCONV_CHAN[ncnn])
                        print(NDCONV[ncnn])
                    isSCAT[ncnn]=False
                    is2D[ncnn]=False
                    isFSL[ncnn]=True
                    cnn_ib[ncnn]=ib
                    ncnn+=1
                    ntmp=ntmp+1

                if Template[ib][itmp][3]=='RELATIV':
                    normrelative=True

        ninput_template_map=len(TemplateMap[ib])
        for itmp in range(ninput_template_map):
            if TemplateMap[ib][itmp][4]!='RELATIV' or nbolo>1:
                if TemplateMap[ib][itmp][1]=='CALIB':
                    icalib[ib]=ntmp
                    docalib=True
                if TemplateMap[ib][itmp][1]=='LINEAR' or TemplateMap[ib][itmp][1]=='CALIB' or TemplateMap[ib][itmp][1]=='POLEFF' or TemplateMap[ib][itmp][1]=='POLANG':
                    ntmp=ntmp+1
                if TemplateMap[ib][itmp][4]=='RELATIV':
                    normrelative=True

        tmp_func[ib]={}
        func_input[ib]={}
        name_func={}
        tmp_func[ib]=np.zeros([nfunc,nx,ny],dtype='float32')
        nfunc=0
        if ncnn>0:
            cnn_idx[ib]=np.zeros([ncnn-cnn_ibolo[ib],nx,ny],dtype='int32')
        ncnn=cnn_ibolo[ib]
        
        for itmp in range(ninput_template):
            if Template[ib][itmp][1]=='SPLINE1D':
                if (rank==0):
                    print(Template[ib][itmp][0],Template[ib][itmp][1],Template[ib][itmp][2],Template[ib][itmp][3])
                tmp_func[ib][nfunc,:,:]=getdata(Template[ib][itmp][2],BeginRingTab[rank],EndRingTab[rank],step=RSTEP)
                name_func[nfunc]=Template[ib][itmp][0]
                func_input[nfunc]=itmp
                nfunc=nfunc+1
                ntotfunc+=4
                ntotMfunc+=Template[ib][itmp][3]
            if Template[ib][itmp][1]=='CNN1D' or Template[ib][itmp][1]=='CNN1D_SCAT' or Template[ib][itmp][1]=='CNN1D_FSL' or Template[ib][itmp][1]=='CNN2D' or Template[ib][itmp][1]=='CNN2D_FSL':
                tmp_cnn=getdata(Template[ib][itmp][2],BeginRingTab[rank],EndRingTab[rank],step=RSTEP)
                if Template[ib][itmp][3]=='LINEAR':
                    cnn_idx[ib][ncnn-cnn_ibolo[ib],:,:]=np.fmod(((tmp_cnn-DOMAIN[ncnn][0])/(DOMAIN[ncnn][1]-DOMAIN[ncnn][0])*XIMAGE_SIZE[ncnn]).astype('int32'),XIMAGE_SIZE[ncnn])

                if Template[ib][itmp][3]=='HIST':
                    nx,ny=tmp_cnn.shape
                    tmp_cnn=((tmp_cnn-DOMAIN[ncnn][0])/(DOMAIN[ncnn][1]-DOMAIN[ncnn][0])*XIMAGE_SIZE[ncnn]*10).astype('int')
                    for kk in range(nx):
                        lhidx=np.where(hdata2[ib][kk,:]>0)[0]
                        if len(lhidx)>0:
                            hh,hx=np.histogram(tmp_cnn[kk,lhidx],range=[0.0,XIMAGE_SIZE[ncnn]*10+1.0],bins=10*XIMAGE_SIZE[ncnn]+1)
                            #plt.subplot(1,3,1)
                            #plt.plot(tmp_cnn[kk,lhidx])
                            #plt.subplot(1,3,2)
                            #plt.plot(hh)
                            hh=np.cumsum(hh)
                            #plt.subplot(1,3,3)
                            #plt.plot(hh)
                            #plt.show()
                            hh=((hh*(XIMAGE_SIZE[ncnn]-1))/hh.max()).astype('int')
                            
                            #plt.plot(tmp_cnn[kk,:],hh[tmp_cnn[kk,:]],'.')
                            #plt.show()
                            tmp_cnn[kk,:]=hh[tmp_cnn[kk,:]]
                        del lhidx
                    cnn_idx[ib][ncnn-cnn_ibolo[ib],:,:]=tmp_cnn
                if (rank==0):
                    print('INFO ',Template[ib][itmp][0],Template[ib][itmp][1],Template[ib][itmp][2],Template[ib][itmp][3],cnn_idx[ib][ncnn-cnn_ibolo[ib],:,:].min(),cnn_idx[ib][ncnn-cnn_ibolo[ib],:,:].max())
                del tmp_cnn
                
                    
                #DO SECOND AXES IF GIVEN
                #if YCNN[ncnn][0]=='BINARY':
                #    cnny_idx[ib]=np.zeros([1,nx,ny],dtype='int32')
                #    tmp_cnn=getdata(YCNN[ncnn][1],BeginRingTab[rank],EndRingTab[rank],step=RSTEP)
                #    cnny_idx[ib][ncnn-cnn_ibolo[ib],:,:]=np.fmod(((tmp_cnn-YCNN[ncnn][2])/(YCNN[ncnn][3]-YCNN[ncnn][2])*YIMAGE_SIZE[ncnn]).astype('int32'),YIMAGE_SIZE[ncnn])
                #    if (rank==0):
                #        print('INFO Y ',YCNN[ncnn][1],cnny_idx[ib][ncnn-cnn_ibolo[ib],:,:].min(),cnny_idx[ib][ncnn-cnn_ibolo[ib],:,:].max())
                #    del tmp_cnn
                
                ncnn+=1
        frelat=np.zeros(ntmp,dtype='int')
        tmp_data[ib]=np.zeros([ntmp,nx,ny],dtype='float32')
        name_tmp={}
        ntmp=0
        for itmp in range(ninput_template):
            if Template[ib][itmp][3]!='RELATIV' or nbolo>1:
                if Template[ib][itmp][1]=='LINEAR' or Template[ib][itmp][1]=='CALIB':
                    if (rank==0):
                        print(Template[ib][itmp][0],Template[ib][itmp][1],' TEMPLATE ',Template[ib][itmp][2],Template[ib][itmp][3])
                    tmp_data[ib][ntmp,:,:]=getdata(Template[ib][itmp][2],BeginRingTab[rank],EndRingTab[rank],step=RSTEP)
                    if Template[ib][itmp][1]=='CALIB':
                        calib[ib]=1.0*(tmp_data[ib][ntmp,:,:]).reshape(nx,ny)
                    if Template[ib][itmp][3]=='RELATIV':
                        frelat[ntmp]=1
                    name_tmp[ntmp]=Template[ib][itmp][0]
                    ntmp=ntmp+1

            if (Template[ib][itmp][1]=='CNN1D_FSL' or Template[ib][itmp][1]=='CNN2D') and PARAMETERONLY==True:
                if (rank==0):
                    print(Template[ib][itmp][0],Template[ib][itmp][1],' FIT DECOD ')
                nx,ny=tmp_data[ib][ntmp,:,:].shape
                tmp_data[ib][ntmp,:,:]=np.random.rand(nx,ny)
                name_tmp[ntmp]=Template[ib][itmp][0]
                ntmp=ntmp+1

        for itmp in range(ninput_template_map):
            if TemplateMap[ib][itmp][4]!='RELATIV' or nbolo>1:
                if TemplateMap[ib][itmp][1]=='LINEAR' or TemplateMap[ib][itmp][1]=='CALIB' or TemplateMap[ib][itmp][1]=='POLEFF' or TemplateMap[ib][itmp][1]=='POLANG':
                    if (rank==0):
                        print(TemplateMap[ib][itmp][0],TemplateMap[ib][itmp][1],' TEMPLATE MAP[',TemplateMap[ib][itmp][2],'] ',TemplateMap[ib][itmp][3],TemplateMap[ib][itmp][4])

                    lidx = np.zeros(len(th.flatten()),dtype='int64')
                    Cfunc.ang2pix(TemplateMap[ib][itmp][2],th.flatten(),ph.flatten(),lidx,len(th.flatten()))
                    lidx = lidx.reshape(th.shape)

                    if TemplateMap[ib][itmp][1]=='POLEFF' or TemplateMap[ib][itmp][1]=='POLANG':
                        namemap=TemplateMap[ib][itmp][3].replace("@POLAR@", "Q")
                        imq=getmap(namemap,TemplateMap[ib][itmp][2])
                        namemap=TemplateMap[ib][itmp][3].replace("@POLAR@", "U")
                        imu=getmap(namemap,TemplateMap[ib][itmp][2])
                        if TemplateMap[ib][itmp][1]=='POLEFF':
                            tmp_data[ib][ntmp,:,:]= imq[lidx]*np.cos(2*psi)+imu[lidx]*np.sin(2*psi)
                        if TemplateMap[ib][itmp][1]=='POLANG':
                            tmp_data[ib][ntmp,:,:]= -imq[lidx]*np.sin(2*psi)+imu[lidx]*np.cos(2*psi)
                    else:
                        tmp_data[ib][ntmp,:,:] = getmap(TemplateMap[ib][itmp][3],TemplateMap[ib][itmp][2])[lidx]

                    if TemplateMap[ib][itmp][1]=='CALIB':
                        calib[ib]=1.0*(tmp_data[ib][ntmp,:,:]).reshape(nx,ny)
                    if TemplateMap[ib][itmp][4]=='RELATIV':
                        frelat[ntmp]=1
                    name_tmp[ntmp]=TemplateMap[ib][itmp][0]
                    ntmp=ntmp+1
        del th
        del ph
        del psi
        comm.Barrier()
        process = psutil.Process(os.getpid())
        if (rank%16==0):
            print('MEM %s %d %.3f MB'%(getframeinfo(currentframe()).lineno,rank,(process.memory_info().rss)/(1024.*1024.)),getloadavg())
            
        #=============================================================================================================================
        #                                    MANAGE ADDED VALUE
        #=============================================================================================================================
        nadd=0
        ninput_add=len(Added[ib])
        for iadd in range(ninput_add):
            nadd=nadd+1
        add_data[ib]=np.zeros([nx,ny],dtype='float32')
        doadd=0
        for iadd in range(ninput_add):
            if (rank==0):
                print(Added[ib][iadd][0],' TEMPLATE ',Added[ib][iadd][1])
            add_data[ib][:,:]=add_data[ib][:,:]+float(Added[ib][iadd][0])*getdata(Added[ib][iadd][1],BeginRingTab[rank],EndRingTab[rank],step=RSTEP)
            doadd=1

        data[ib]     = data[ib][idx[0],:].astype('float32')
        hdata[ib]    = hdata[ib][idx[0],:].astype('float32')
        hdata2[ib]    = hdata2[ib][idx[0],:].astype('float32')
        hidx[ib]     = hidx[ib][idx[0],:]
        if ncnn>0:
            cnn_idx[ib] = cnn_idx[ib][:,idx[0],:].astype('int32')
        if ntmp>0:
            tmp_data[ib] = tmp_data[ib][:,idx[0],:].astype('float32')
        if nfunc>0:
            tmp_func[ib] = tmp_func[ib][:,idx[0],:].astype('float32')
            
        for icnn in range(ncnn):
            if isFSL[icnn]==True:
                if rank==0:
                    print(IMRIDX[icnn],12*NSIDEFSL[icnn]*NSIDEFSL[icnn])
                imrtmp=getdata(IMRIDX[icnn],BeginRingTab[rank],EndRingTab[rank],step=RSTEP,rgsz=12*NSIDEFSL[icnn]*NSIDEFSL[icnn])
                imrtmp=imrtmp[idx[0],:]
                nx,ny=imrtmp.shape
                lnring=np.zeros([size],dtype='int32')
                comm.Allgather(sendbuf=(np.array([nx],dtype='int32'),MPI.INT),recvbuf=(lnring,MPI.INT))

                # complicated exchange to avoid copy all the IMRIDX timeline in memory
                for irank in range(size):
                    if irank%NUMBEROFTHREAD==0:
                        if rank==irank:
                            IMRIDX[icnn]=np.zeros([lnring.sum(),12*NSIDEFSL[icnn]*NSIDEFSL[icnn]],dtype='float32')
                            n0=0
                            for jrank in range(size):
                                if irank!=jrank:
                                    tmp=comm.recv(source=jrank, tag=11)
                                    IMRIDX[icnn][n0:n0+lnring[jrank],:]=tmp.reshape(lnring[jrank],ny)
                                    del tmp
                                else:
                                    IMRIDX[icnn][n0:n0+lnring[jrank],:]=imrtmp
                                n0+=lnring[jrank]
                        else:
                            comm.send((imrtmp.astype('float32')).flatten(), dest=irank, tag=11)
                del imrtmp
                
                if rank%NUMBEROFTHREAD==0:
                    print(IMRIDX[icnn][0,:].sum())
                    IMRIDX[icnn][:,:]/=2.16710895E4
                    IMRIDX[icnn][:,:]*=0.00013406
                    nnx,nny=ROTIDX[icnn]['data'].shape
                    cnn_imref[icnn] = tf.constant(IMRIDX[icnn][:,:].astype('float32'))
                    print('FSL DECOD [%d] '%(rank),cnn_imref[icnn].get_shape().as_list())
                    rotidx[icnn]=tf.constant(ROTIDX[icnn]['data'].flatten())
                    #protidx[icnn]=tf.constant(ROTIDX[icnn]['pidx'].flatten())
                    #if DOREFERENCE==False:
                    #    maskrot[icnn]=tf.constant(ROTIDX[icnn]['mask'].flatten())
                    #else:
                    maskrot[icnn]=tf.constant(np.zeros([nny],dtype='float32'))
                    print('ROIDX ',ROTIDX[icnn]['data'].max())
                    
                    norient=tf.constant(np.array([128],dtype='int32'))
                    npidx=nny
                    print('npidx',npidx)
                    del ROTIDX[icnn]
                    del IMRIDX[icnn]
                    
                    
        add_data[ib] = add_data[ib][idx[0],:].astype('float32')
        if docalib:
            calib[ib]    = calib[ib][idx[0],:].astype('float32')

        nx,ny=data[ib].shape
        ridx[ib]=(np.repeat(ridx[ib][idx[0]],ny).reshape(nx,ny))

        data[ib]     = data[ib][hdata[ib]>0]
        hidx[ib]     = hidx[ib][hdata[ib]>0]
        if ncnn>0:
            cnn_idx[ib] = cnn_idx[ib][:,hdata[ib]>0]
        if ntmp>0:
            tmp_data[ib] = tmp_data[ib][:,hdata[ib]>0]
        if nfunc>0:
            tmp_func[ib] = tmp_func[ib][:,hdata[ib]>0]
        if doadd>0:
            add_data[ib] = add_data[ib][hdata[ib]>0]
        ridx[ib]     = ridx[ib][hdata[ib]>0]
        if docalib:
            calib[ib]    = calib[ib][hdata[ib]>0]
        hdata2[ib]    = hdata2[ib][hdata[ib]>0]
        hdata[ib]    = hdata[ib][hdata[ib]>0]
    
       # Find valid ring and replace it in ridx table
        fring= np.zeros([EndRing+1])
        fring[ridx[ib]]=1
        tfring = np.zeros([EndRing+1]) 
        comm.Allreduce(fring,tfring)

        val_ring[ib]=np.where(tfring>0)[0]
        nring[ib]=len(val_ring[ib])
        
        if rank==0:
            print("NRING[%d] (valid) %d"%(ib,nring[ib]),icalib,nring[ib])
        tfring[val_ring[ib]]=np.arange(nring[ib])
        ridx[ib]=tfring[ridx[ib]].astype('int')

    if showdata==True:
        exit(0)
    cnn_ibolo[nbolo]=ncnn
        
    if ncnn==0:
        l_Nside_cnn=Nside
    else:
        l_Nside_cnn=Nside_cnn
    #=============================================================================================================================
    #                                    MPI EXCHANGE DATA
    #=============================================================================================================================
    if size>0:
        if rank==0:
            rpidx=(np.random.rand(l_Nside_cnn*l_Nside_cnn*12)*size).astype('int32')
        else:
            rpidx=np.zeros([l_Nside_cnn*l_Nside_cnn*12],dtype='int32')


        comm.Bcast(rpidx, root=0)
     
        for ib in range(nbolo):

            lhidx=ChangeReso(hidx[ib],Nside,l_Nside_cnn)

            if rank==0:
                print('Transfer Data Bolo ',rank,ib)

            lpidx=rpidx[lhidx]
            
            scounts=np.zeros([size],dtype='int')
            for i in range(size):
                scounts[i]=(lpidx==i).sum()

            tidx=np.argsort(lpidx)

            rcounts = np.array(comm.alltoall(scounts))

            if rank==0:
                print(rcounts,scounts)

            nlocal=int(rcounts.sum())

            buffer1=np.zeros([nlocal],dtype='float32')
            comm.Alltoallv(sendbuf=(data[ib][tidx].astype('float32'),scounts,MPI.FLOAT),recvbuf=(buffer1,rcounts,MPI.FLOAT)) 
            del data[ib]
            data[ib]=buffer1

            buffer1=np.zeros([nlocal],dtype='float32')
            comm.Alltoallv(sendbuf=(hdata2[ib][tidx].astype('float32'),scounts,MPI.FLOAT),recvbuf=(buffer1,rcounts,MPI.FLOAT)) 
            del hdata2[ib]
            hdata2[ib]=buffer1

            if docalib:
                buffer1=np.zeros([nlocal],dtype='float32')
                comm.Alltoallv(sendbuf=(calib[ib][tidx].astype('float32'),scounts,MPI.FLOAT),recvbuf=(buffer1,rcounts,MPI.FLOAT))
                del calib[ib]
                calib[ib]=buffer1

            buffer2=np.zeros([nlocal],dtype='float32')
            comm.Alltoallv(sendbuf=(hdata[ib][tidx].astype('float32'),scounts,MPI.FLOAT),recvbuf=(buffer2,rcounts,MPI.FLOAT))
            del hdata[ib]
            hdata[ib]=buffer2
        
            buffer3=np.zeros([nlocal],dtype='int32')
            comm.Alltoallv(sendbuf=(hidx[ib][tidx].astype('int32'),scounts,MPI.INT),recvbuf=(buffer3,rcounts,MPI.INT)) 
            del hidx[ib]
            hidx[ib]=buffer3
        
            buffer4=np.zeros([nlocal],dtype='int32')
            comm.Alltoallv(sendbuf=(ridx[ib][tidx].astype('int32'),scounts,MPI.INT),recvbuf=(buffer4,rcounts,MPI.INT))
            del ridx[ib]
            ridx[ib]=buffer4

            if ntmp>0:
                buffer5=np.zeros([nlocal],dtype='float32')
                ltmp_data=np.zeros([ntmp,nlocal],dtype='float32')
                for i in range(ntmp):
                    comm.Alltoallv(sendbuf=((tmp_data[ib][i,tidx]).reshape(len(tidx)).astype('float32'),scounts,MPI.FLOAT),recvbuf=(buffer5,rcounts,MPI.FLOAT)) 
                    ltmp_data[i,:]=buffer5
                del tmp_data[ib]
                tmp_data[ib]=ltmp_data
                del buffer5

            if nfunc>0:
                buffer5=np.zeros([nlocal],dtype='float32')
                ltmp_func=np.zeros([nfunc,nlocal],dtype='float32')
                for i in range(nfunc):
                    comm.Alltoallv(sendbuf=((tmp_func[ib][i,tidx]).reshape(len(tidx)).astype('float32'),scounts,MPI.FLOAT),recvbuf=(buffer5,rcounts,MPI.FLOAT)) 
                    ltmp_func[i,:]=buffer5
                del tmp_func[ib]
                tmp_func[ib]=ltmp_func
                del buffer5

            if ncnn>0:
                buffer5=np.zeros([nlocal],dtype='int32')
                ltmp_func=np.zeros([ncnn-cnn_ibolo[ib],nlocal],dtype='int32')
                for i in range(cnn_ibolo[ib],cnn_ibolo[ib+1]):
                    comm.Alltoallv(sendbuf=((cnn_idx[ib][i-cnn_ibolo[ib],tidx]).reshape(len(tidx)).astype('int32'),scounts,MPI.INT),recvbuf=(buffer5,rcounts,MPI.INT)) 
                    ltmp_func[i-cnn_ibolo[ib],:]=buffer5
                del cnn_idx[ib]
                cnn_idx[ib]=ltmp_func
                del buffer5


            if doadd==1:
                buffer2=np.zeros([nlocal],dtype='float32')
                comm.Alltoallv(sendbuf=(add_data[ib][tidx].astype('float32'),scounts,MPI.FLOAT),recvbuf=(buffer2,rcounts,MPI.FLOAT)) 
                del add_data[ib]
                add_data[ib]=buffer2
    
   
        
    fpix= np.zeros([12*Nside*Nside],dtype='int')
    for ib in range(nbolo):
        fpix[hidx[ib]]=1
    
    val_pidx=np.where(fpix>0)[0]
    npixel=len(val_pidx)
    fpix[val_pidx]=np.arange(npixel)
    for ib in range(nbolo):
        hidx[ib]=fpix[hidx[ib]]
    
    ndata=np.zeros([nbolo],dtype='int')
    for ib in range(nbolo):
        ndata[ib]=data[ib].shape[0]
    print("Rank %d [%d]: NPIXEL = %d, valid samples = "%(rank,ntmp,npixel),ndata)

    #data[ib]=np.exp(-0.01*(cnn_idx[ib][0,:]-64)**2)
    avv=0
    navv=0
    for ib in range(nbolo):
        #data[ib]=np.cos(ridx[ib]/1000.)
        avv+=(hdata[ib]*data[ib]).sum()
        navv+=(hdata[ib]).sum()

    for ib in range(nbolo):
        data[ib]-=avv/navv
        
    if domodel:
        for ib in range(nbolo):
            #tmp_data[ib]=np.random.randn(ntmp,ndata[ib])

            model=np.zeros([nring[ib]+ntmp],dtype='float')
            model[0:nring[ib]]=np.cos((2*np.pi*(np.arange(nring[ib])+100*ib))/nring[ib])
            model[0:nring[ib]]-=model[0:nring[ib]].mean()
            data[ib]=model[ridx[ib]]
            calib[ib]=1.0*(tmp_data[ib][ntmp-1,:])

            if ntmp>0:
                vv=0.01*(1+np.arange(ntmp).reshape(1,ntmp))+0.003*ib
                if docalib:
                    vv[0,icalib[ib]]=1.02+0.003*ib
                data[ib]+=np.dot(vv,tmp_data[ib]).reshape(ndata[ib])
                model[nring[ib]:]=vv

            np.save(Out_Offset[ib]+'_mod',model)

    if rank==0:
        for ib in range(nbolo):
            print('DATA INFO',rank,data[ib].min(),data[ib].max(),hdata[ib].min(),hdata[ib].max(),ridx[ib].min(),ridx[ib].max(),hidx[ib].min(),hidx[ib].max())

    comm.Barrier()
    #=============================================================================================================================
    #                                    MANAGE AMPLITUDE TEMPLATE TO FIT
    #=============================================================================================================================
    coef2        = np.zeros([ntmp],dtype='float')
    coef         = np.zeros([ntmp],dtype='float')
    ncoef        = np.zeros([ntmp],dtype='float')

    for ib in range(nbolo):
        for i in range(ntmp):
            coef2[i]+=(tmp_data[ib][i,:]**2).sum()
            coef[i]+=(tmp_data[ib][i,:]).sum()
            ncoef[i]+=ndata[ib]
            
    tcoef2        = np.zeros([ntmp],dtype='float')
    tcoef         = np.zeros([ntmp],dtype='float')
    tncoef        = np.zeros([ntmp],dtype='float')

    comm.Allreduce((coef2,MPI.DOUBLE),(tcoef2,MPI.DOUBLE))
    comm.Allreduce((coef,MPI.DOUBLE),(tcoef,MPI.DOUBLE))
    comm.Allreduce((ncoef,MPI.DOUBLE),(tncoef,MPI.DOUBLE))

    coef2=np.sqrt(tcoef2/tncoef-(tcoef/tncoef)**2)
    for i in range(ntmp):
        if rank==0:
            print('REGUL ',name_tmp[i],coef2[i])
        for ib in range(nbolo):
            tmp_data[ib][i,:]/=coef2[i]

    #=============================================================================================================================
    #                                    DESIGN TENSORFLOW MINIMIZATION
    #=============================================================================================================================
    
    
    buffer = Cfunc.finit(Nside,npixel,nbolo,ntmp,nfunc,ntotfunc,ntotMfunc,int(docalib==1),int(doadd==1),int(ncnn>0),cnn_ibolo,rank,size)

    if nfunc>0:
        for ib in range(nbolo):
            for ifunc in range(nfunc):
                if Template[ib][func_input[ifunc]][1]=='SPLINE1D':
                    extrem=np.array([tmp_func[ib][ifunc,:].min()],dtype='float')
                    lextrem=np.zeros([size],dtype='float')
                    comm.Allgather((extrem,MPI.DOUBLE),(lextrem,MPI.DOUBLE))
                    themin=lextrem.min()
                    extrem=np.array([tmp_func[ib][ifunc,:].max()],dtype='float')
                    comm.Allgather((extrem,MPI.DOUBLE),(lextrem,MPI.DOUBLE))
                    themax=lextrem.max()

                    hh,xhh=np.histogram(tmp_func[ib][ifunc,:],range=[themin,themax*1.001],bins=20000)
                    hh=hh.astype('float')
                    lhh=0.0*hh
                    comm.Allreduce((hh,MPI.DOUBLE),(lhh,MPI.DOUBLE))
                    lhh=np.cumsum(lhh)
                    lhh/=lhh.max()

                    tmp_func[ib][ifunc,:]=lhh[np.floor((tmp_func[ib][ifunc,:]-xhh[0])/(xhh[1]-xhh[0])).astype('int')]
                    
                    Cfunc.finitSpline1D(buffer,ifunc,Template[ib][func_input[ifunc]][3],0.0,1.0)

    for ib in range(nbolo):
        Cfunc.fcount_data_per_pixel(buffer,hidx[ib].astype('int32'),ndata[ib])

    comm.Barrier()
    process = psutil.Process(os.getpid())
    if (rank%16==0):
        print('MEM %s %d %.3f MB'%(getframeinfo(currentframe()).lineno,rank,(process.memory_info().rss)/(1024.*1024.)),getloadavg())

    Cfunc.fallocbolo(buffer)

    if ncnn>0:
        if rank%NUMBEROFTHREAD==0:
            TFNULL=tf.constant(np.zeros([1],dtype='float32'))
            
        if DOFULLYCONNECTED:
            fc2_weights={}
            fc2_biases={}
        dconv_weights={}
        dconv_biases={}
        co_filt={}
        ico_filt={}
        co2_filt={}
        ico2_filt={}
        param={}
        theoffset={}
        cnn_corr={}
        circular={}

        scatter_objective_cos = {}
        scatter_0_objective = {}
        scatter_1_objective = {}
        scatter_2_objective = {}
        
        for k in range(ncnn):
            if isSCAT[k]==True:
                # Scattering Transform filters
                phi_tf={}
                psi_tf={}
                psi1_tf={}
                psi2_tf={}
                order_1_tf={}
                order_2_tf={}

                scatter_objective = {}

                scatter_model = {}

                scatter_0_std = {}
                scatter_1_std = {}
                scatter_2_std = {}

                scatter_mask = {}
                
                #Debug
                check_scatter_model_cos = {}
                check_scatter_0_model = {}
                check_scatter_1_model = {}
                check_scatter_2_model = {}
                check_cost = {}

        if (rank==0):
            print('MEM %s %d %.3f MB'%(getframeinfo(currentframe()).lineno,rank,(process.memory_info().rss)/(1024.*1024.)),getloadavg())
        var2fit={}
        for k in range(ncnn):
            np.random.seed(1234+k)
            if PARAMETERONLY:
                if rank%NUMBEROFTHREAD==0:
                    vv=np.array([VAR2FIT[k]],dtype='float32')
                    var2fit[k] = vv
                    train_var2fit = tf.compat.v1.placeholder(tf.float32,shape=1)

                    if DOFULLYCONNECTED:
                        fc2_weights[k] = tf.constant(np.load('%s_%d_fc2w.npy'%(TFLEARN,k)))
                        fc2_biases[k]  = tf.constant(np.load('%s_%d_fc2b.npy'%(TFLEARN,k)))
                        if rank==0:
                            print('READ FL WEIGHTS')
                    
                    dconv_weights[k]={}
                    dconv_biases[k]={}
                    for i in range(NDCONV[k]):
                        dconv_weights[k][i] = tf.constant(np.load('%s_%d_w_%d.npy'%(TFLEARN,k,i)))
                        dconv_biases[k][i]  = tf.constant(np.load('%s_%d_b_%d.npy'%(TFLEARN,k,i)))

                    #dconv_weights[k][NDCONV[k]-1] = tf.Variable(np.load('%s_%d_w_%d.npy'%(TFLEARN,k,NDCONV[k]-1)))
                    #dconv_biases[k][NDCONV[k]-1]  = tf.Variable(np.load('%s_%d_b_%d.npy'%(TFLEARN,k,NDCONV[k]-1)))
            else:
                if rank%NUMBEROFTHREAD==0:
                    vv=np.array([VAR2FIT[k]],dtype='float32')
                    var2fit[k] = vv
                    train_var2fit = tf.compat.v1.placeholder(tf.float32,shape=1)
                if DOFULLYCONNECTED:
                    if is2D[k]==False:
                        tmp=np.random.randn(NHIDDEN[k], XIMAGE_SIZE[k] // (WSCALE**NDCONV[k]) *  NUM_DCONV_CHAN[k][0]).astype('float32')
                    else:
                        tmp=np.random.randn(NHIDDEN[k], XIMAGE_SIZE[k] // (WSCALE**NDCONV[k]) * YIMAGE_SIZE[k] // (WSCALE**NDCONV[k]) *  NUM_DCONV_CHAN[k][0]).astype('float32')
                    comm.Bcast(tmp,root=0)
                    if rank%NUMBEROFTHREAD==0:
                        fc2_weights[k] = tf.Variable(tmp)
                    if is2D[k]==False:
                        tmp=np.zeros(XIMAGE_SIZE[k] // (WSCALE**NDCONV[k])  * NUM_DCONV_CHAN[k][0]).astype('float32')
                    else:
                        tmp=np.zeros(XIMAGE_SIZE[k] // (WSCALE**NDCONV[k])  * YIMAGE_SIZE[k] // (WSCALE**NDCONV[k])  * NUM_DCONV_CHAN[k][0]).astype('float32')
                    if rank%NUMBEROFTHREAD==0:
                        fc2_biases[k] = tf.Variable(tmp)
                dconv_weights[k]={}
                dconv_biases[k]={}
                for i in range(NDCONV[k]):
                    if is2D[k]==False:
                        tmp=np.random.randn(KERNELSZ[k], 1,NUM_DCONV_CHAN[k][i+1], NUM_DCONV_CHAN[k][i]).astype('float32')
                    else:
                        tmp=np.random.randn(KERNELSZ[k], KERNELSZ[k],NUM_DCONV_CHAN[k][i+1], NUM_DCONV_CHAN[k][i]).astype('float32')
                    comm.Bcast(tmp,root=0)
                    if rank%NUMBEROFTHREAD==0:
                        dconv_weights[k][i] = tf.Variable(tmp)
                    tmp=np.zeros(NUM_DCONV_CHAN[k][i+1]).astype('float32')
                    if rank%NUMBEROFTHREAD==0:
                        dconv_biases[k][i] = tf.Variable(tmp)

        
        nndata=0
        for ib in range(nbolo):
            idx=np.where(hdata2[ib]>0)[0]
            nndata+=len(idx)

        if NUMBEROFTHREAD!=1:
            vv=np.array([nndata],dtype='int32')
            lvv=np.array([0],dtype='int32')
        
            ncomm.Allreduce([vv,1,MPI.INT],lvv,op=MPI.MAX)
            mndata=lvv[0]
        else:
            mndata=nndata

        if rank%NUMBEROFTHREAD==0:
            train_ridx={}
            train_wridx={}
            if BATCH_SIZE==-1:
                train_data = tf.compat.v1.placeholder(tf.float32,shape=mndata)
                train_hh   = tf.compat.v1.placeholder(tf.float32,shape=mndata)
                train_pidx = tf.compat.v1.placeholder(tf.int32,shape=mndata)
            
                for icnn in range(ncnn):
                    train_ridx[icnn] = tf.compat.v1.placeholder(tf.int32,shape=mndata)
                    train_wridx[icnn] = tf.compat.v1.placeholder(tf.float32,shape=mndata)
            else:
                train_data = tf.compat.v1.placeholder(tf.float32,shape=BATCH_SIZE)
                train_hh   = tf.compat.v1.placeholder(tf.float32,shape=BATCH_SIZE)
                train_pidx = tf.compat.v1.placeholder(tf.int32,shape=BATCH_SIZE)
                for icnn in range(ncnn):
                    train_ridx[icnn] = tf.compat.v1.placeholder(tf.int32,shape=BATCH_SIZE)
                    train_wridx[icnn] = tf.compat.v1.placeholder(tf.float32,shape=BATCH_SIZE)
            
        

        if rank%NUMBEROFTHREAD==0:
            def conv_transpose_4(itest,testdec):
                dimin=itest.get_shape().as_list()
                dimin2=testdec.get_shape().as_list()

                widx2=np.zeros([4,dimin2[2],dimin[3]],dtype='int32')
                for k in range(dimin2[2]):
                  for l in range(dimin[3]):
                    widx2[:,k,l]=np.array([1,0,5,6],dtype='int32')*dimin2[2]*dimin[3]+k*dimin[3]+l
                widx2=tf.constant(widx2,dtype='int32')

                widx=tf.constant(np.repeat([1,0,5,6],dimin2[2]*dimin[3]).reshape(4,dimin2[2],dimin[3]),dtype='int32')
                wcomp=tf.constant(np.array([1.0,0,0,0,0,1.0,1.0]),dtype='float32')
                testdec3=tf.reshape(tf.gather(wcomp,widx),[1,4,dimin2[2],dimin[3]])
                testdec2=tf.reshape(tf.gather(tf.reshape(testdec,[7*dimin2[2]*dimin[3]]),widx2),[1,4,dimin2[2],dimin[3]])

                sidx=np.zeros([dimin[0],4*dimin[1],1,dimin2[2]],dtype='int32')
                for k in range(dimin[0]):
                  for l in range(dimin2[2]):
                    sidx[k,0,0,l]=(dimin2[2]*(4*k+2)+l)
                    sidx[k,1,0,l]=(dimin2[2]*(4*k+3)+l)
                    sidx[k,4*dimin[1]-1,0,l]=(dimin2[2]*(4*k+1)+l)
                sidx=tf.constant(sidx)
                sidx2=tf.constant(np.array([1,0,dimin[1]-1,dimin[1]-1]))

                ltest=tf.nn.conv2d_transpose(itest,testdec,strides=[1, WSCALE, 1, 1],padding='SAME',output_shape=[dimin[0],4*dimin[1],1,dimin2[2]])

                sval3=tf.gather(itest,sidx2,axis=1)

                sval2=tf.reduce_sum(sval3*testdec3*testdec2,axis=3)

                sval=tf.gather(tf.reshape(sval2,[dimin[0]*4*dimin2[2]]),sidx)
                return(sval+ltest)

            #=====================================================================================
            # SCATTERING TRANSFORM FUNCTIONS
            #=====================================================================================

            def fit_cosine(signal,nrings,ximage_size,filter_flag,tico,ty):

                scatter = tf.reshape(signal,[nrings,ximage_size]); 
                #scatter = scatter[:nrings,:]

                if filter_flag:
                    scatter = tf.matmul(tico,tf.matmul(tf.transpose(tico),scatter))

                scatter_cos = tf.reshape(tf.matmul(scatter,ty),[nrings,-1]); 

                return scatter_cos

            def scattering(signal,phi,psi,psi1,psi2,order_1,average=1):

                ## PAD VARIABLE
                [Nring,Nh] = signal.get_shape().as_list()

                signal = tf.transpose(signal)

                m = np.int(np.ceil(np.log2(Nring)))

                r = 2**m - Nring
                signal_padded = tf.concat([signal,np.zeros([Nh,r])],axis=1)

                # GO TO FREQ DOMAIN (APPLY FFT)
                signal_fft = tf.signal.fft((tf.cast(signal_padded,dtype=tf.complex64)))

                # SCATT 0
                scatt_0 = tf.multiply(signal_fft,tf.cast(phi,dtype = tf.complex64))
                scatt_0_time = tf.real(tf.signal.ifft(scatt_0))

                ## ORDER 1
                # FILTER WITH SCALE FILTERS
                signal_filtered = tf.multiply(tf.expand_dims(signal_fft,axis=1),tf.cast(psi,dtype = tf.complex64))

                # GO BACK TO TIME DOMAIN (APPLY FFT)
                signal_filtered_time = tf.signal.ifft(signal_filtered)

                # TAKE MODULUS
                signal_filter_time_abs = tf.abs(signal_filtered_time)

                # GO TO FREQ DOMAIN (APPLY FFT)
                signal_abs_fft = tf.signal.fft((tf.cast(signal_filter_time_abs,dtype=tf.complex64)))

                # FILTER WITH LARGE SCALE FILTER
                signal_abs_fft_large_filtered = tf.multiply(signal_abs_fft,tf.cast(phi,dtype=tf.complex64))

                # GO BACK TO TIME DOMAIN (APPLY FFT)
                scatter_complex =  tf.signal.ifft(signal_abs_fft_large_filtered)

                #TAKE REAL PART AND UNPAD
                scatter_1 =  tf.real(scatter_complex)[:,:,:Nring]

                ## ORDER 2

                # FILTER WITH SCALE 1 FILTERS
                signal_filtered_order_2 = tf.multiply(tf.expand_dims(signal_fft,axis=1),tf.cast(psi1,dtype = tf.complex64))

                # GO BACK TO TIME DOMAIN (APPLY FFT)
                signal_filtered_time_order_2 = tf.signal.ifft(signal_filtered_order_2)

                # TAKE MODULUS
                signal_filter_time_abs_order_2 = tf.abs(signal_filtered_time_order_2)

                # GO TO FREQ DOMAIN (APPLY FFT)
                signal_abs_fft_order_2 = tf.signal.fft((tf.cast(signal_filter_time_abs_order_2,dtype=tf.complex64)))

                # FILTER WITH SCALE 2 FILTERS
                signal_abs_fft_second_filtered_order_2 = tf.multiply(signal_abs_fft_order_2,tf.cast(psi2,dtype = tf.complex64))

                # GO BACK TO TIME DOMAIN (APPLY FFT)
                scatter_complex_order_2 =  tf.signal.ifft(signal_abs_fft_second_filtered_order_2)

                # TAKE MODULUS
                scatter_abs_order_2 = tf.abs(scatter_complex_order_2)

                # GO TO FREQ DOMAIN (APPLY FFT)
                scatter_abs_fft_order_2 = tf.signal.fft((tf.cast(scatter_abs_order_2,dtype=tf.complex64)))

                # FILTER WITH LARGE SCALE FILTER
                scatter_abs_fft_large_filtered_order_2 = tf.multiply(scatter_abs_fft_order_2,tf.cast(phi,dtype=tf.complex64))

                # GO BACK TO TIME DOMAIN (APPLY FFT)
                scatter_complex_order_2 =  tf.signal.ifft(scatter_abs_fft_large_filtered_order_2)

                #TAKE REAL PART AND UNPAD
                scatter_2 =  tf.real(scatter_complex_order_2)[:,:,:Nring]

                scatter_0 = tf.transpose(scatt_0_time[:,:Nring])

                if average:
                    scatter_0 = tf.reduce_mean(scatter_0,axis=0)

                    scatter_1 = tf.reduce_mean(scatter_1,axis=2)

                    scatter_2 = tf.reduce_mean(scatter_2,axis=2)

                return scatter_0,scatter_1,scatter_2

            def normalize_and_select(scatter_0,scatter_1,scatter_2,order_1,scatter_mask,wscatt):
                ## SCATTERING MOMENTS NORMALIZATION
                scatter_1= tf.transpose(tf.transpose(scatter_1)/scatter_0)
                scatter_2 = scatter_2/tf.gather(scatter_1,tf.reshape(order_1,[-1]),axis=1)

                ## SELECT ONLY VALID SCATTERING COEFFICIENTS
                scatter_2 = tf.boolean_mask(scatter_2,scatter_mask,axis=1)
                
                scatter_0 = wscatt[0]*scatter_0
                scatter_1 = wscatt[1]*scatter_1
                scatter_2 = wscatt[2]*scatter_2
                
                return scatter_1,scatter_2

            def scatter_cost(scatter_0_model,scatter_1_model,scatter_2_model,
                             scatter_0_objective,scatter_1_objective,scatter_2_objective,
                             scatter_0_std,scatter_1_std,scatter_2_std,
                             COMPUTE_SCATTER_OBJECTIVE,SCATTER_CLEAN,
                             orders=[1,1,1]):

                scatter_relu = 0.0

                print('VERIFY DIMS')
                print(scatter_2_objective.shape)
                print(scatter_2_model.shape)
                print(scatter_2_std.shape)

                if COMPUTE_SCATTER_OBJECTIVE:

                    if not SCATTER_CLEAN:

                        if orders[0]:
                          scatter_relu = scatter_relu + tf.reduce_sum(tf.square(tf.reshape(scatter_0_model-scatter_0_objective,[-1])))

                        if orders[1]:
                          scatter_relu = scatter_relu + tf.reduce_sum(tf.square(tf.reshape(scatter_1_model-scatter_1_objective,[-1])))

                        if orders[2]:
                          scatter_relu = scatter_relu + tf.reduce_sum(tf.square(tf.reshape(scatter_2_model-scatter_2_objective,[-1])))

                    else:

                        if orders[0]:
                          scatter_relu = scatter_relu + (-tf.reduce_sum(tf.log(1-tf.exp(-(tf.square(tf.reshape(scatter_0_model-scatter_0_objective,[-1])))))))     

                        if orders[1]:
                          scatter_relu = scatter_relu + (-tf.reduce_sum(tf.log(1-tf.exp(-(tf.square(tf.reshape(scatter_1_model-scatter_1_objective,[-1])))))))

                        if orders[2]:
                          scatter_relu = scatter_relu + (-tf.reduce_sum(tf.log(1-tf.exp(-(tf.square(tf.reshape(scatter_2_model-scatter_2_objective,[-1])))))))

                else: 

                    if not SCATTER_CLEAN:

                        if  orders[0]:
                          scatter_relu = scatter_relu + tf.reduce_sum(tf.square(tf.reshape((scatter_0_model-scatter_0_objective)/scatter_0_std,[-1])))

                        if orders[1]:
                          scatter_relu = scatter_relu + tf.reduce_sum(tf.square(tf.reshape((scatter_1_model-scatter_1_objective)/scatter_1_std,[-1])))

                        if orders[2]:
                          scatter_relu = scatter_relu + tf.reduce_sum(tf.square(tf.reshape((scatter_2_model-scatter_2_objective)/scatter_2_std,[-1])))

                    else:

                        if  orders[0]:
                          scatter_relu = scatter_relu + (-tf.reduce_sum(tf.log(1-tf.exp(-(tf.square(tf.reshape((scatter_0_model-scatter_0_objective)/scatter_0_std,[-1])))))))

                        if orders[1]:
                          scatter_relu = scatter_relu + (-tf.reduce_sum(tf.log(1-tf.exp(-(tf.square(tf.reshape((scatter_1_model-scatter_1_objective)/scatter_1_std,[-1]))))))) 

                        if orders[2]:
                          scatter_relu = scatter_relu + (-tf.reduce_sum(tf.log(1-tf.exp(-(tf.square(tf.reshape((scatter_2_model-scatter_2_objective)/scatter_2_std,[-1])))))))


                return scatter_relu
            
            def lrelu(x, alpha=0.3):
                return tf.maximum(x, tf.multiply(x, alpha))

            def cnn_decod_noFL(parameter,dconv_weights,dconv_biases,k,circular=True):

                hidden = parameter
                if DORELU==True:
                  hidden = lrelu(hidden)
                relu_shape = hidden.get_shape().as_list()
                if rank==0:
                    print(relu_shape,parameter.get_shape().as_list())
                if isFSL[k]==False:
                    relu = tf.reshape(hidden,[relu_shape[0], XIMAGE_SIZE[k] // (WSCALE**NDCONV[k]) ,
                                              1 , NUM_DCONV_CHAN[k][0]])
                else:
                    relu = tf.reshape(hidden,[relu_shape[0], 12*32*32 // (WSCALE**NDCONV[k]) ,
                                              1 , NUM_DCONV_CHAN[k][0]])
                # DeConvolution kernel
                for i in range(NDCONV[k]):
                    relu_shape = relu.get_shape().as_list()
                    if rank==0:
                        print('CIRCULAR :',circular)
                    if circular==True:
                        conv = conv_transpose_4(relu,dconv_weights[i])
                    else:
                        if isFSL[k]==False:
                            conv = tf.nn.conv2d_transpose(relu,dconv_weights[i],
                                                          strides=[1, WSCALE, 1, 1],padding='SAME',
                                                          output_shape=[relu_shape[0],XIMAGE_SIZE[k] // (WSCALE**(NDCONV[k]-i-1)),
                                                                        1,NUM_DCONV_CHAN[k][i+1]])
                        else:
                            conv = tf.nn.conv2d_transpose(relu,dconv_weights[i],
                                                          strides=[1, WSCALE, 1, 1],padding='SAME',
                                                          output_shape=[relu_shape[0], 12*32*32 // (WSCALE**(NDCONV[k]-i-1)),
                                                                        1,NUM_DCONV_CHAN[k][i+1]])

                    relu=tf.nn.bias_add(conv, dconv_biases[i])
                    relu_shape = relu.get_shape().as_list()

                    if DORELU==True: #and 
                        if rank==0:
                            print("Non LINEAR")
                        if i!=NDCONV[k]-1:
                            relu = lrelu(relu)
                        else:
                            relu = lrelu(relu,alpha=0)
                            
                        
                    if rank==0:
                        print(relu_shape)

                relu=tf.reshape(relu,[relu_shape[0]*relu_shape[1]])/REGULAMP[k]
                if NHARM[k]!=0:
                    wfilt=tf.matmul(tf.reshape(relu,[relu_shape[0],relu_shape[1]]),co_filt[k])
                    print(wfilt.get_shape().as_list())
                    wfilt=tf.matmul(co2_filt[k],wfilt)
                    print(wfilt.get_shape().as_list())
                    wfilt=tf.matmul(ico2_filt[k],wfilt)
                    print(wfilt.get_shape().as_list())
                    wfilt=tf.matmul(wfilt,ico_filt[k])
                    print(wfilt.get_shape().as_list())
                    filt=tf.reshape(wfilt,[relu_shape[0]*relu_shape[1]])
                    if NHARM[k]>0:
                        relu = relu-filt
                    else:
                        relu = filt
                    
                return relu
            
            def cnn_decod_noFL_2D(parameter,dconv_weights,dconv_biases,k,circular=True):

                hidden = parameter
                if DORELU==True:
                  hidden = lrelu(hidden)
                relu_shape = hidden.get_shape().as_list()
                if rank==0:
                    print(relu_shape,parameter.get_shape().as_list())
                relu = tf.reshape(hidden,[1, YIMAGE_SIZE[k] // (WSCALE**NDCONV[k]) ,
                                          XIMAGE_SIZE[k] // (WSCALE**NDCONV[k]) , NUM_DCONV_CHAN[k][0]])
                # DeConvolution kernel
                for i in range(NDCONV[k]):
                    relu_shape = relu.get_shape().as_list()
                    if rank==0:
                        print('CIRCULAR :',circular)
                    conv = tf.nn.conv2d_transpose(relu,dconv_weights[i],
                                                  strides=[1, WSCALE, WSCALE, 1],padding='SAME',
                                                  output_shape=[1,YIMAGE_SIZE[k] // (WSCALE**(NDCONV[k]-i-1)),
                                                                XIMAGE_SIZE[k] // (WSCALE**(NDCONV[k]-i-1)),
                                                                NUM_DCONV_CHAN[k][i+1]])

                    relu=tf.nn.bias_add(conv, dconv_biases[i])
                    relu_shape = relu.get_shape().as_list()

                    if DORELU==True: # and i!=NDCONV[k]-1:
                        if rank==0:
                            print("Non LINEAR")
                        if i!=NDCONV[k]-1:
                            relu = lrelu(relu)
                        else:
                            relu = lrelu(relu,alpha=0)
                        
                    if rank==0:
                        print(relu_shape)

                relu=tf.reshape(relu,[XIMAGE_SIZE[k]*YIMAGE_SIZE[k]])/REGULAMP[k]

                if NHARM[k]!=0:
                    wfilt=tf.matmul(tf.reshape(relu,[XIMAGE_SIZE[k],YIMAGE_SIZE[k]]),co_filt[k])
                    print(wfilt.get_shape().as_list())
                    wfilt=tf.matmul(co2_filt[k],wfilt)
                    print(wfilt.get_shape().as_list())
                    wfilt=tf.matmul(ico2_filt[k],wfilt)
                    print(wfilt.get_shape().as_list())
                    wfilt=tf.matmul(wfilt,ico_filt[k])
                    print(wfilt.get_shape().as_list())
                    filt=tf.reshape(wfilt,[XIMAGE_SIZE[k]*YIMAGE_SIZE[k]])
                    if NHARM[k]>0:
                        relu = relu-filt
                    else:
                        relu = filt
                return relu

            def cnn_decod(parameter,fc2_weights,
                          fc2_biases,
                          dconv_weights,
                          dconv_biases,k,circular=True):

                param = parameter
                
                # convert into image/bias
                if DOFULLYCONNECTED == False:
                    hidden = parameter
                else:
                    hidden = tf.matmul(parameter, fc2_weights) + fc2_biases
                    
                if DORELU==True:
                  hidden = lrelu(hidden)

                relu_shape = hidden.get_shape().as_list()
                if rank==0:
                    print(relu_shape,param.get_shape().as_list())
                if isFSL[k]==True:
                    relu = tf.reshape(hidden,[relu_shape[0], 12*32*32 // (WSCALE**NDCONV[k]) ,
                                              1 , NUM_DCONV_CHAN[k][0]])
                else:
                    relu = tf.reshape(hidden,[relu_shape[0], XIMAGE_SIZE[k] // (WSCALE**NDCONV[k]) ,
                                              1 , NUM_DCONV_CHAN[k][0]])
                # DeConvolution kernel
                for i in range(NDCONV[k]):
                    relu_shape = relu.get_shape().as_list()
                    if rank==0:
                        print('CIRCULAR :',circular)
                    if circular==True:
                        conv = conv_transpose_4(relu,dconv_weights[i])
                    else:
                        if isFSL[k]==False:
                            conv = tf.nn.conv2d_transpose(relu,dconv_weights[i],
                                                          strides=[1, WSCALE, 1, 1],padding='SAME',
                                                          output_shape=[relu_shape[0],XIMAGE_SIZE[k] // (WSCALE**(NDCONV[k]-i-1)),
                                                                        1,NUM_DCONV_CHAN[k][i+1]])
                        else:
                            conv = tf.nn.conv2d_transpose(relu,dconv_weights[i],
                                                          strides=[1, WSCALE, 1, 1],padding='SAME',
                                                          output_shape=[relu_shape[0],12*32*32// (WSCALE**(NDCONV[k]-i-1)),
                                                                        1,NUM_DCONV_CHAN[k][i+1]])

                    relu=tf.nn.bias_add(conv, dconv_biases[i])
                    relu_shape = relu.get_shape().as_list()

                    if DORELU==True: # and i!=NDCONV[k]-1:
                        if rank==0:
                            print("Non LINEAR")
                        if i!=NDCONV[k]-1:
                            relu = lrelu(relu)
                        else:
                            relu = lrelu(relu,alpha=0)
                        
                    if rank==0:
                        print(relu_shape)

                relu=tf.reshape(relu,[relu_shape[0]*relu_shape[1]])/REGULAMP[k]
                if NHARM[k]!=0:
                    wfilt=tf.matmul(tf.reshape(relu,[relu_shape[0],relu_shape[1]]),co_filt[k])
                    print(wfilt.get_shape().as_list())
                    wfilt=tf.matmul(co2_filt[k],wfilt)
                    print(wfilt.get_shape().as_list())
                    wfilt=tf.matmul(ico2_filt[k],wfilt)
                    print(wfilt.get_shape().as_list())
                    wfilt=tf.matmul(wfilt,ico_filt[k])
                    print(wfilt.get_shape().as_list())
                    filt=tf.reshape(wfilt,[relu_shape[0]*relu_shape[1]])
                    if NHARM[k]>0:
                        relu = relu-filt
                    else:
                        relu = filt
                return relu

            def cnn_decod_2D(parameter,fc2_weights,
                          fc2_biases,
                          dconv_weights,
                             dconv_biases,k,circular=True):

                param = parameter

                # convert into image/bias
                if DOFULLYCONNECTED == False:
                    hidden = parameter
                else:
                    hidden = tf.matmul(parameter, fc2_weights) + fc2_biases

                if DORELU==True:
                  hidden = lrelu(hidden)

                relu_shape = hidden.get_shape().as_list()
                if rank==0:
                    print(relu_shape,param.get_shape().as_list())
                relu = tf.reshape(hidden,[relu_shape[0], YIMAGE_SIZE[k] // (WSCALE**NDCONV[k]) ,
                                         XIMAGE_SIZE[k] // (WSCALE**NDCONV[k]) , NUM_DCONV_CHAN[k][0]])
                # DeConvolution kernel
                for i in range(NDCONV[k]):
                    relu_shape = relu.get_shape().as_list()
                    if rank==0:
                        print('CIRCULAR :',circular)
                    if circular==True:
                        conv = conv_transpose_4(relu,dconv_weights[i])
                    else:
                        conv = tf.nn.conv2d_transpose(relu,dconv_weights[i],
                                                      strides=[1, WSCALE, WSCALE, 1],padding='SAME',
                                                      output_shape=[relu_shape[0],YIMAGE_SIZE[k] // (WSCALE**(NDCONV[k]-i-1)),
                                                                    XIMAGE_SIZE[k] // (WSCALE**(NDCONV[k]-i-1)),
                                                                    NUM_DCONV_CHAN[k][i+1]])

                    relu=tf.nn.bias_add(conv, dconv_biases[i])
                    relu_shape = relu.get_shape().as_list()

                    if DORELU==True: # and i!=NDCONV[k]-1:
                        if rank==0:
                            print("Non LINEAR")
                        if i!=NDCONV[k]-1:
                            relu = lrelu(relu)
                        else:
                            relu = lrelu(relu,alpha=0)
                        
                    if rank==0:
                        print(relu_shape)

                relu=tf.reshape(relu,[relu_shape[0]*relu_shape[1]*relu_shape[2]])/REGULAMP[k]
                if NHARM[k]!=0:
                    wfilt=tf.matmul(tf.reshape(relu,[relu_shape[0],relu_shape[1]]),co_filt[k])
                    print(wfilt.get_shape().as_list())
                    wfilt=tf.matmul(wfilt,co2_filt[k])
                    print(wfilt.get_shape().as_list())
                    wfilt=tf.matmul(wfilt,ico2_filt[k])
                    print(wfilt.get_shape().as_list())
                    wfilt=tf.matmul(wfilt,ico_filt[k])
                    print(wfilt.get_shape().as_list())
                    filt=tf.reshape(wfilt,[relu_shape[0]*relu_shape[1]])
                    if NHARM[k]>0:
                        relu = relu-filt
                    else:
                        relu = filt
                return relu

            def cnn_model(data, parameter,ihh,ipidx,iridx,wridx,NPIX,
                          fc2_weights,
                          fc2_biases,
                          dconv_weights,
                          dconv_biases,circular,train=False):
                """The Model definition."""
                avv=TFNULL
                fsl = TFNULL
                scatter_model_cos = TFNULL
                scatter_0_model = TFNULL
                scatter_1_model = TFNULL
                scatter_2_model = TFNULL
                vdestrip  = 0
                norm_relu=0
                testscat = 0
                relu=0.0
                
                for k in range(ncnn):
                    if is2D[k]==False:
                        vdestrip1 = cnn_decod(parameter[k],
                                              fc2_weights[k],
                                              fc2_biases[k],
                                              dconv_weights[k],
                                              dconv_biases[k],k,
                                              circular=circular[k])
                    else:
                        vdestrip1 = cnn_decod_2D(parameter[k],
                                                 fc2_weights[k],
                                                 fc2_biases[k],
                                                 dconv_weights[k],
                                                 dconv_biases[k],k,
                                                 circular=circular[k])
                        
                    tsize=vdestrip1.get_shape().as_list()
                    
                    if isSCAT[k]==True:
                        testscat=1
                        
                        scatter_model_cos = fit_cosine(vdestrip1,cnn_ntotring[k],XIMAGE_SIZE[k],FILTER_SCATTER[k],tico[k],ty[k])
                        
                        scatter_0_model,scatter_1_model,scatter_2_model = scattering(scatter_model_cos,phi_tf[k],psi_tf[k],psi1_tf[k],psi2_tf[k],order_1_tf[k])
                         
                        if SCATTER_NORM[k]==True:
                            scatter_1_model,scatter_2_model = normalize_and_select(scatter_0_model,scatter_1_model,scatter_2_model,order_1_tf[k],scatter_mask[k],WSCATT[k])

                        norm_relu = norm_relu + SCATTER_WALL[k]*scatter_cost(scatter_0_model,scatter_1_model,scatter_2_model,
                                                                             scatter_0_objective[k],scatter_1_objective[k],scatter_2_objective[k],
                                                                             scatter_0_std[k],scatter_1_std[k],scatter_2_std[k],
                                                                             COMPUTE_SCATTER_OBJECTIVE[k],SCATTER_CLEAN[k])
                    
                    if isFSL[k]==True:
                        # precompute all rotation of the FSL
                        fsl = vdestrip1
                        #vdestrip1 = tf.reshape(tf.transpose(tf.reshape(tf.tile(fsl,norient),[128,npidx])),[128*npidx])
                        #vdestrip1 = tf.reshape(tf.gather(vdestrip1,rotidx[k]),[12*32*32,128])
                        vdestrip1 = tf.transpose(tf.reshape(tf.gather(fsl,rotidx[k]),[128,12*32*32]))
                        if DOREFERENCE==True:
                            vdestrip1 = vdestrip1 + 1E-3 * theoffset[k]
                        vdestrip1 = tf.matmul(cnn_imref[k],vdestrip1)
                        l_size=vdestrip1.get_shape().as_list()
                        vdestrip1=tf.reshape(vdestrip1,[l_size[0]*l_size[1]])
                        tsize=vdestrip1.get_shape().as_list()
                        #relu=relu+1E4*tf.reduce_sum(tf.square(maskrot[k]*fsl))
                        print('DECOD SIZE',vdestrip1.get_shape().as_list())
                    else:
                        avv=avv+tf.reduce_sum(tf.square(tf.reduce_sum(tf.reshape(vdestrip1,[tsize[0]//XIMAGE_SIZE[k],XIMAGE_SIZE[k]]),1)))
                    
                    vdestrip  =  vdestrip +  train_var2fit*wridx[k]*tf.gather(vdestrip1,iridx[k])

                omap=tf.math.unsorted_segment_sum(ihh*vdestrip, ipidx, NPIX)
                hmap=tf.math.unsorted_segment_sum(ihh,          ipidx, NPIX)
                imap=tf.math.unsorted_segment_sum(ihh*data   ,  ipidx, NPIX)
                nmap=tf.gather(hmap,ipidx)

                ww = ihh/((1+ihh)*nmap+1E-7)
                if DOREFERENCE==True:
                    vdestrip = ww*(vdestrip*nmap)
                    creshape = ww*(data*nmap)
                    norm_relu = TFNULL
                else:
                    odest=tf.gather(omap,ipidx)
                    vdestrip = ww*(vdestrip*nmap - odest)
                    creshape = ww*(data*nmap - tf.gather(imap,ipidx))

                    #1E6 make the Wloss efficiency around 1.0 
                    shmap=tf.reduce_sum(hmap)/(1E6)

                    norm_relu=shmap * tf.reduce_sum(tf.square(omap/(1E-6+hmap)))

                tmp=(vdestrip-creshape)
                
                relu=relu+tf.reduce_sum(tf.square(tmp))

                return relu,norm_relu,avv,scatter_model_cos,scatter_0_model,scatter_1_model,scatter_2_model,fsl #,omap,imap,hmap,vdestrip1

            def cnn_model_noFL(data, parameter,ihh,ipidx,iridx,wridx,NPIX,dconv_weights,
                               dconv_biases,circular,train=False):
                """The Model definition."""

                fsl = TFNULL
                norm_relu = 0
                avv = TFNULL
                vdestrip = 0
                scatter_model_cos = TFNULL
                scatter_0_model = TFNULL
                scatter_1_model = TFNULL
                scatter_2_model = TFNULL
                testscat = 0
                relu=0.0
                
                for k in range(ncnn):
                    if is2D[k]==False:
                        vdestrip1 = cnn_decod_noFL(parameter[k],
                                                   dconv_weights[k],
                                                   dconv_biases[k],k, 
                                                   circular=circular[k])
                    else:
                        vdestrip1 = cnn_decod_noFL_2D(parameter[k],
                                                      dconv_weights[k],
                                                      dconv_biases[k],k, 
                                                      circular=circular[k])
                    tsize=vdestrip1.get_shape().as_list()

                    if isSCAT[k]==True:
                        testscat=1
                        
                        scatter_model_cos = fit_cosine(vdestrip1,cnn_ntotring[k],XIMAGE_SIZE[k],FILTER_SCATTER[k],tico[k],ty[k])
                    
                        scatter_0_model,scatter_1_model,scatter_2_model = scattering(scatter_model_cos,phi_tf[k],psi_tf[k],psi1_tf[k],psi2_tf[k],order_1_tf[k])
                        
                        if SCATTER_NORM[k]==True:
                            scatter_1_model,scatter_2_model = normalize_and_select(scatter_0_model,scatter_1_model,scatter_2_model,order_1_tf[k],scatter_mask[k],WSCATT[k])

                        #norm_relu = norm_relu + SCATTER_WALL[k]* tf.reduce_sum(tf.square(scatter_model_cos-scatter_objective_cos[k]))
                        
                        norm_relu = norm_relu + SCATTER_WALL[k]* scatter_cost(scatter_0_model,scatter_1_model,scatter_2_model,
                                                                              scatter_0_objective[k],scatter_1_objective[k],scatter_2_objective[k],
                                                                              scatter_0_std[k],scatter_1_std[k],scatter_2_std[k],
                                                                              COMPUTE_SCATTER_OBJECTIVE[k],SCATTER_CLEAN[k])
            
                    if isFSL[k]==True:
                        # precompute all rotation of the FSL
                        fsl = vdestrip1
                        print(fsl.get_shape().as_list())
                        #vdestrip1 = tf.reshape(tf.transpose(tf.reshape(tf.tile(fsl,norient),[norient,npidx])),[norient*npidx])
                        #vdestrip1 = tf.reshape(tf.gather(vdestrip1,rotidx[k]),[12*32*32,128])
                        vdestrip1 = tf.transpose(tf.reshape(tf.gather(fsl,rotidx[k]),[128,12*32*32]))

                        if DOREFERENCE==True:
                            vdestrip1 = vdestrip1 + 1E-3 * theoffset[k]
                        
                        print(vdestrip1.get_shape().as_list())
                        vdestrip1 = tf.matmul(cnn_imref[k],vdestrip1)
                        l_size=vdestrip1.get_shape().as_list()
                        vdestrip1=tf.reshape(vdestrip1,[l_size[0]*l_size[1]])
                        tsize=vdestrip1.get_shape().as_list()
                        #relu=relu+1E4*tf.reduce_sum(tf.square(maskrot[k]*fsl))
                    else:
                        avv=avv+tf.reduce_sum(tf.square(tf.reduce_sum(tf.reshape(vdestrip1,[tsize[0]//XIMAGE_SIZE[k],XIMAGE_SIZE[k]]),1)))
                    
                    vdestrip  =  vdestrip +  train_var2fit*wridx[k]*tf.gather(vdestrip1,iridx[k])
                
                    
                omap=tf.math.unsorted_segment_sum(ihh*vdestrip, ipidx, NPIX)
                hmap=tf.math.unsorted_segment_sum(ihh,          ipidx, NPIX)
                imap=tf.math.unsorted_segment_sum(ihh*data   ,  ipidx, NPIX)
                nmap=tf.gather(hmap,ipidx)

                ww = ihh/((ihh)*nmap+1E-7)
                
                if DOREFERENCE:
                    vdestrip = ihh*vdestrip
                    creshape = ihh*data
                    norm_relu = TFNULL
                else:
                    odest=tf.gather(omap,ipidx)
                    vdestrip = ww*(vdestrip*nmap - odest)
                    idest=tf.gather(imap,ipidx)
                    creshape = ww*(data*nmap - idest)
                
                    #1E6 make the Wloss efficiency around 1.0 
                    shmap=tf.reduce_sum(hmap)/(1E6)

                    if testscat==0:
                        norm_relu = shmap * tf.reduce_sum(tf.square(omap/(1E-6+hmap)))
                    
                tmp=(vdestrip-creshape)

                relu=relu+tf.reduce_sum(tf.square(tmp))

                return relu,norm_relu,avv,scatter_model_cos,scatter_0_model,scatter_1_model,scatter_2_model,fsl #,omap,imap,hmap,vdestrip1

    if ncnn>0:
        cnn_ntotring={}
        for k in range(ncnn):
            cnn_ntotring[k]=0

    for ib in range(nbolo):
        
        if doadd:
            theadd=add_data[ib].astype('float')
        else:
            theadd=np.zeros([1],dtype='float')

        if docalib:
            thecalib=calib[ib].astype('float')
        else:
            thecalib=np.zeros([1],dtype='float')

        if ntmp>0:
            the_data=tmp_data[ib].astype('float')
        else:
            the_data=np.zeros([1],dtype='float')

        if nfunc>0:
            the_func=tmp_func[ib].astype('float')
        else:
            the_func=np.zeros([1],dtype='float')

        if ncnn>0:
            l_train_ridx={}
            l_train_wridx={}
            icorrection=np.zeros([cnn_ibolo[ib+1]-cnn_ibolo[ib],ndata[ib]],dtype='int32')
            for k in range(ncnn):
                if cnn_ib[k]==ib:
                    nnx=0
                    if YCNN[k][0]=='BINARY':
                        if is2D[k]==True:
                            icorrection[k-cnn_ibolo[ib],:]=((cnny_idx[ib][k-cnn_ibolo[ib],:]).astype('int32')*XIMAGE_SIZE[k]+cnn_idx[ib][k-cnn_ibolo[ib],:]).astype('int32')
                            cnn_ntotring[k]+=YIMAGE_SIZE[k]
                            nnx=YIMAGE_SIZE[k]
                        else:
                            icorrection[k-cnn_ibolo[ib],:]=((cnny_idx[ib][k-cnn_ibolo[ib],:]).astype('int32')*XIMAGE_SIZE[k]+cnn_idx[ib][k-cnn_ibolo[ib],:]).astype('int32')
                            cnn_ntotring[k]+=(cnny_idx[ib][k-cnn_ibolo[ib],:]).max()+1
                            nnx=ntotr
                        
                    if YCNN[k][0]=='RING':
                        if YCNN[k][1]==0:
                            ntotr=nring[ib]
                        else:
                            if YCNN[k][1]==-1:
                                ntotr=nring[ib]
                            else:
                                ntotr=YCNN[k][1]
                            
                        if is2D[k]==True and isFSL[k]==False:
                            if ntotr==nring[ib]:
                                if YCNN[k][1]==-1:
                                    iy2=np.load(YCNN[k][2])
                                    icorrection[k-cnn_ibolo[ib],:]=((iy2[val_ring[ib][ridx[ib]]]*XIMAGE_SIZE[k]).astype('int32')+cnn_idx[ib][k-cnn_ibolo[ib],:]).astype('int32')
                                else:
                                    icorrection[k-cnn_ibolo[ib],:]=((((YIMAGE_SIZE[k]*(ridx[ib]%ntotr))//ntotr)*XIMAGE_SIZE[k]).astype('int32')+cnn_idx[ib][k-cnn_ibolo[ib],:]).astype('int32')
                            else:
                                icorrection[k-cnn_ibolo[ib],:]=(((YIMAGE_SIZE[k]*((val_ring[ib][ridx[ib]]-BeginRing)%ntotr))//ntotr)*XIMAGE_SIZE[k]+cnn_idx[ib][k-cnn_ibolo[ib],:]).astype('int32')
                            cnn_ntotring[k]+=YIMAGE_SIZE[k]
                            nnx=YIMAGE_SIZE[k]
                        else:
                            icorrection[k-cnn_ibolo[ib],:]=(((ridx[ib]%ntotr).astype('int32'))*XIMAGE_SIZE[k]).astype('int32')+cnn_idx[ib][k-cnn_ibolo[ib],:].astype('int32')
                            cnn_ntotring[k]+=ntotr
                            nnx=ntotr
                            
                    if YCNN[k][0]=='ONE':
                        icorrection[k-cnn_ibolo[ib],:]=cnn_idx[ib][k-cnn_ibolo[ib],:].astype('int32')
                        cnn_ntotring[k]+=1
                        nnx=1
                    if YCNN[k][0]=='RIDX':
                        theidx=(YIMAGE_SIZE[k]*(val_ring[ib][ridx[ib]]-BeginRing)/(EndRing - BeginRing)).astype('int32')
                        icorrection[k-cnn_ibolo[ib],:]=(theidx)*XIMAGE_SIZE[k]+cnn_idx[ib][k-cnn_ibolo[ib],:].astype('int32')
                        cnn_ntotring[k]+=YIMAGE_SIZE[k]
                        del theidx
                        nnx=YIMAGE_SIZE[k]
                    if YCNN[k][0]=='HPIX':
                        l_nside=YCNN[k][1]
                        lndata=len(val_pidx[hidx[ib]])
                        th=np.zeros([lndata],dtype='float')
                        ph=np.zeros([lndata],dtype='float')
                        nhidx=np.zeros([lndata],dtype='int64')
                        Cfunc.pix2ang(Nside,val_pidx[hidx[ib]],th,ph,lndata)
                        Cfunc.ang2pix(l_nside,th,ph,nhidx,lndata,nest=True)

                        icorrection[k-cnn_ibolo[ib],:]=(nhidx)*XIMAGE_SIZE[k]+cnn_idx[ib][k-cnn_ibolo[ib],:].astype('int32')
                        cnn_ntotring[k]+=l_nside*l_nside*12
                        nnx=l_nside*l_nside*12
                        del nhidx
                        del th
                        del ph
                    
                    if rank%NUMBEROFTHREAD==0:
                        if NHARM[k]!=0:
                            nharm=int(abs(NHARM[k]))
                            co=np.zeros([nharm*2+1,XIMAGE_SIZE[k]],dtype='float32')
                            ico=np.zeros([nharm*2+1,XIMAGE_SIZE[k]],dtype='float32')
                            co2=np.zeros([nharm*2+1,nnx],dtype='float32')
                            ico2=np.zeros([nharm*2+1,nnx],dtype='float32')
                            x=np.arange(XIMAGE_SIZE[k])/float(XIMAGE_SIZE[k])*2*np.pi
                            # do the survey repeating
                            y=np.arange(nnx)/(float(11194-240)/RSTEP)*2*np.pi
                            co[0,:]=1.0/np.sqrt(float(XIMAGE_SIZE[k]))
                            ico[0,:]=1.0/np.sqrt(float(XIMAGE_SIZE[k]))
                            co2[0,:]=1.0/np.sqrt(float(nnx))
                            ico2[0,:]=1.0/np.sqrt(float(nnx))

                            for i in range(nharm):
                                co[2*i+1,:]=np.cos(x*(i+1))/np.sqrt(0.5*float(XIMAGE_SIZE[k]))
                                co[2*i+2,:]=np.sin(x*(i+1))/np.sqrt(0.5*float(XIMAGE_SIZE[k]))
                                ico[2*i+1,:]=np.cos(x*(i+1))/np.sqrt(0.5*float(XIMAGE_SIZE[k]))
                                ico[2*i+2,:]=np.sin(x*(i+1))/np.sqrt(0.5*float(XIMAGE_SIZE[k]))
                                co2[2*i+1,:]=np.cos(y*(i+1))/np.sqrt(0.5*float(nnx))
                                co2[2*i+2,:]=np.sin(y*(i+1))/np.sqrt(0.5*float(nnx))
                                ico2[2*i+1,:]=np.cos(y*(i+1))/np.sqrt(0.5*float(nnx))
                                ico2[2*i+2,:]=np.sin(y*(i+1))/np.sqrt(0.5*float(nnx))
                                if i>nharm//2:
                                    tmp=(1+np.cos(np.pi*(i-nharm//2)/float(nharm//2)))/2.0
                                    print(tmp)
                                    ico[2*i+1,:]*=tmp
                                    ico[2*i+2,:]*=tmp
                                    ico2[2*i+1,:]*=tmp
                                    ico2[2*i+2,:]*=tmp
                            co_filt[k]=tf.constant(co.transpose())
                            ico_filt[k]=tf.constant(ico)
                            co2_filt[k]=tf.constant(co2)
                            ico2_filt[k]=tf.constant(ico2.transpose())
            

                    if is2D[k]==True and isFSL[k]==False:
                        Cfunc.falloccorrection(buffer,XIMAGE_SIZE[k]*YIMAGE_SIZE[k],k)
                    else:
                        Cfunc.falloccorrection(buffer,XIMAGE_SIZE[k]*cnn_ntotring[k],k)

                    if DOFULLYCONNECTED:
                        tmp=np.random.randn(cnn_ntotring[k],NHIDDEN[k])/100.0
                        if rank==0:
                            print('Fully CONNECTED ',[cnn_ntotring[k],NHIDDEN[k]])
                    else:
                        if is2D[k]==False:
                            if isFSL[k]==False:
                                tmp=np.random.randn(cnn_ntotring[k],NUM_DCONV_CHAN[k][0]*XIMAGE_SIZE[k]//(WSCALE**(NDCONV[k])))/100.0
                                if rank==0:
                                    print('NOT FULLY CONNECTED 1D ',[cnn_ntotring[k],NUM_DCONV_CHAN[k][0]*XIMAGE_SIZE[k]//(WSCALE**(NDCONV[k]))])
                            else:
                                tmp=np.random.randn(1,NUM_DCONV_CHAN[k][0]*12*32*32//(WSCALE**(NDCONV[k])))/100.0
                                if rank==0:
                                    print('NOT FULLY CONNECTED 1D ',[1,NUM_DCONV_CHAN[k][0]*12*32*32//(WSCALE**(NDCONV[k]))])
                        else:
                            tmp=np.random.randn(YIMAGE_SIZE[k]//(WSCALE**(NDCONV[k])),NUM_DCONV_CHAN[k][0]*XIMAGE_SIZE[k]//(WSCALE**(NDCONV[k])))/100.0
                            if rank==0:
                                print('NOT FULLY CONNECTED 2D ',[YIMAGE_SIZE[k]//(WSCALE**(NDCONV[k])),NUM_DCONV_CHAN[k][0]*XIMAGE_SIZE[k]//(WSCALE**(NDCONV[k]))])

                    comm.Bcast(tmp, root=0)

                    if rank%NUMBEROFTHREAD==0:

                        if isSCAT[k]==True:
                            #=====================================================================================
                            # BUILD SCATTERING TRANSFORM FILTER BANK AND AUXILIARY BASES
                            #=====================================================================================

                            # Signal length

                            # Timeline size
                            T = 2**np.int(np.ceil(np.log2(cnn_ntotring[k])))

                            # Filters
                            x=np.arange(T)/T
                            phi=np.fft.fft(np.exp(-8*((x-T/2)**2)/(T**2)))


                            rr2=(x**2)*16
                            psi=np.zeros([NSCATTER[k],T])

                            for j in range(NSCATTER[k]):
                                psi[j,:]=(np.exp(-rr2*2**(2*j))*rr2).astype('float')
                                psi[j,:]/=(psi[j,:]).max()


                            psi1 = np.repeat(psi,NSCATTER[k],axis=0)

                            psi2 = np.tile(psi,[NSCATTER[k],1])

                            phi_tf[k]=tf.constant(phi.reshape([1,T]))
                            psi_tf[k]=tf.constant(psi.reshape([NSCATTER[k],T]))
                            psi1_tf[k]=tf.constant(psi1.reshape([NSCATTER[k]**2,T]))
                            psi2_tf[k]=tf.constant(psi2.reshape([NSCATTER[k]**2,T]))

                            order_1 = np.repeat(np.arange(NSCATTER[k]),NSCATTER[k]).reshape([-1,1]).astype('int32')
                            order_2 = np.tile(np.arange(NSCATTER[k]).reshape(-1,1),[NSCATTER[k],1]).astype('int32')

                            order_1_tf[k] = tf.constant(order_1)
                            order_2_tf[k] = tf.constant(order_2)

                            scatter_mask[k] = tf.reshape((order_2_tf[k]-order_1_tf[k])>=0,[-1])

                            #=====================================================================================
                            # COSINE BASE FOR SCATTERING TRANSFORM
                            #=====================================================================================

                            # Cosine base
                            y = np.zeros([NHARM_SCATTER[k],XIMAGE_SIZE[k]],dtype='float32')
                            for i in range(NHARM_SCATTER[k]//2):
                                y[i*2,:]=np.cos(((1+i)*np.arange(XIMAGE_SIZE[k])/(XIMAGE_SIZE[k]*1.0))*2*np.pi)/np.sqrt(XIMAGE_SIZE[k]//2)
                                y[i*2+1,:]=np.sin(((1+i)*np.arange(XIMAGE_SIZE[k])/(XIMAGE_SIZE[k]*1.0))*2*np.pi)/np.sqrt(XIMAGE_SIZE[k]//2)
                            ty[k] = tf.constant(y.transpose())


                            ## TICO
                            ico=np.zeros([cnn_ntotring[k],NHARM_LARGE_SCALE[k]*2+1],dtype='float32')
                            for i in range(NHARM_LARGE_SCALE[k]):
                                ico[:,2*i]=np.cos((i+1)*np.arange(cnn_ntotring[k])*2*np.pi/cnn_ntotring[k])*np.sqrt(2.0/cnn_ntotring[k])
                                ico[:,2*i+1]=np.sin((i+1)*np.arange(cnn_ntotring[k])*2*np.pi/cnn_ntotring[k])*np.sqrt(2.0/cnn_ntotring[k])
                            ico[:,NHARM_LARGE_SCALE[k]*2]=np.sqrt(1/cnn_ntotring[k])
                            tico[k]=tf.constant(ico)


                            #=====================================================================================
                            # COMPUTE OR LOAD OBJECTIVE SCATTERING TRANSFORM
                            #=====================================================================================
                            if COMPUTE_SCATTER_OBJECTIVE[k]:

                                scatter_objective[k] = np.load(SCATTER_OBJECTIVE_FILE[k]).astype('float32')

                                nx,ny=scatter_objective[k].shape

                                if nx!=cnn_ntotring[k] or ny!=XIMAGE_SIZE[k]:
                                    x=np.arange(nx)
                                    y=np.arange(ny)
                                    f_scipy = interpolate.interp2d(x, y, scatter_objective[k].flatten(), kind='cubic')
                                    nx2=cnn_ntotring[k]
                                    ny2=XIMAGE_SIZE[k]
                                    x2=np.arange(nx2)/nx2*nx
                                    y2=np.arange(ny2)/ny2*ny
                                    scatter_objective[k] = (f_scipy(x2, y2).transpose()).astype('float32')
                                    print('CHANGE TEMPLATE SIZE ',nx,ny,scatter_objective[k].shape)

                                scatter_objective_cos[k] = fit_cosine(scatter_objective[k],cnn_ntotring[k],XIMAGE_SIZE[k],FILTER_SCATTER[k],tico[k],ty[k])

                                scatter_0_objective[k],scatter_1_objective[k],scatter_2_objective[k] = scattering(scatter_objective_cos[k],phi_tf[k],psi_tf[k],psi1_tf[k],psi2_tf[k],order_1_tf[k])

                                # Dummy scatter std
                                scatter_0_std[k] = tf.ones(scatter_0_objective[k].shape)
                                scatter_1_std[k] = tf.ones(scatter_1_objective[k].shape)
                                scatter_2_std[k] = tf.ones(scatter_2_objective[k].shape)

                                #DEBUG
                                '''
                                scatter_0_std[k] = tf.ones(scatter_0_objective[k].shape)
                                scatter_1_std[k] = tf.ones(scatter_1_objective[k].shape)
                                scatter_2_std[k] = tf.ones(scatter_2_objective[k].shape)

                                scatter_0_objective[k] = 1E-1*tf.random.normal(scatter_0_objective[k].shape)
                                scatter_1_objective[k] = 1E-1*tf.random.normal(scatter_1_objective[k].shape)
                                scatter_2_objective[k] = 1E-1*tf.random.normal(scatter_2_objective[k].shape)
                                '''

                                ## NORMALIZED AND SELECT ONLY VALID SCATTERING COEFFICIENTS
                                if SCATTER_NORM[k]:
                                    scatter_1_objective[k],scatter_2_objective[k]= normalize_and_select(scatter_0_objective[k],scatter_1_objective[k],scatter_2_objective[k],order_1_tf[k],scatter_mask[k],WSCATT[k])

                            else:

                                scatter_0_objective[k] = tf.constant(np.fromfile(SCATTER_0_MEAN_FILE[k],dtype='float32').reshape([NHARM_SCATTER[k]]))
                                scatter_1_objective[k] = tf.constant(np.fromfile(SCATTER_1_MEAN_FILE[k],dtype='float32').reshape([NHARM_SCATTER[k],-1]))
                                scatter_2_objective[k] = tf.constant(np.fromfile(SCATTER_2_MEAN_FILE[k],dtype='float32').reshape([NHARM_SCATTER[k],-1]))

                                scatter_0_std[k] = tf.constant(np.fromfile(SCATTER_0_STD_FILE[k],dtype='float32').reshape([NHARM_SCATTER[k]]))
                                scatter_1_std[k] = tf.constant(np.fromfile(SCATTER_1_STD_FILE[k],dtype='float32').reshape([NHARM_SCATTER[k],-1]))
                                scatter_2_std[k] = tf.constant(np.fromfile(SCATTER_2_STD_FILE[k],dtype='float32').reshape([NHARM_SCATTER[k],-1]))

                                #DEBUG
                                '''
                                scatter_0_std[k] = tf.ones(scatter_0_objective[k].shape)
                                scatter_1_std[k] = tf.ones(scatter_1_objective[k].shape)
                                scatter_2_std[k] = tf.ones(scatter_2_objective[k].shape)

                                scatter_0_objective[k] = 1E-1*np.random.normal(scatter_0_objective[k].shape)
                                scatter_1_objective[k] = 1E-1*tf.random.normal(scatter_1_objective[k].shape)
                                scatter_2_objective[k] = 1E-1*tf.random.normal(scatter_2_objective[k].shape)
                                '''

                                ## SELECT ONLY VALID SCATTERING COEFFICIENTS (ALREADY NORMALIZED)
                                scatter_2_objective[k] = tf.boolean_mask(scatter_2_objective[k],scatter_mask[k],axis=1)

                            ## SELECT ONLY VALID SCATTERING COEFFICIENTS
                            scatter_2_std[k] = tf.boolean_mask(scatter_2_std[k],scatter_mask[k],axis=1)

                            '''
                            print('DEBUG LOAD')

                            print(scatter_0_objective[k])
                            print(scatter_1_objective[k])
                            print(scatter_2_objective[k])

                            print(scatter_0_std[k])
                            print(scatter_1_std[k])
                            print(scatter_2_std[k])                
                            '''

                            #scatter_1_objective[k],scatter_2_objective[k]= normalize_and_select(scatter_0_objective[k],scatter_1_objective[k],scatter_2_objective[k],order_1_tf[k],scatter_mask[k])

                            '''                               
                            print('DEBUG NORMALIZATION')

                            print(scatter_mask[k])               

                            print(scatter_0_objective[k])
                            print(scatter_1_objective[k])
                            print(scatter_2_objective[k])

                            print(scatter_0_std[k])
                            print(scatter_1_std[k])
                            print(scatter_2_std[k])
                            '''
                        else:
                            scatter_objective_cos[k]=TFNULL
                            scatter_0_objective[k]=TFNULL
                            scatter_1_objective[k]=TFNULL
                            scatter_2_objective[k]=TFNULL

                        if PARAMETERONLY==True:
                            param[k]=tf.Variable(0*np.load('%s_%d_param.npy'%(TFLEARN,k)))
                        else:
                            param[k]=tf.Variable(tmp.astype('float32'))
                            if DOREFERENCE==True:
                                ltmp=np.zeros([1],dtype='float32')
                                theoffset[k]=tf.Variable(ltmp)

                        circular[k]=(CIRCULAR[k]=='CIRCULAR')
                        if DOFULLYCONNECTED:
                            if is2D[k]==False:
                                cnn_corr_tmp=cnn_decod(param[k],
                                                      fc2_weights[k],
                                                      fc2_biases[k],
                                                      dconv_weights[k],
                                                      dconv_biases[k],k,circular=circular[k])
                            else:
                                cnn_corr_tmp=cnn_decod_2D(param[k],
                                                          fc2_weights[k],
                                                          fc2_biases[k],
                                                          dconv_weights[k],
                                                          dconv_biases[k],k,circular=circular[k])
                        else:
                            if is2D[k]==False:
                                cnn_corr_tmp=cnn_decod_noFL(param[k],
                                                           dconv_weights[k],
                                                           dconv_biases[k],k,circular=circular[k])
                            else:
                                cnn_corr_tmp=cnn_decod_noFL_2D(param[k],
                                                               dconv_weights[k],
                                                               dconv_biases[k],k,circular=circular[k])
                        if isFSL[k]==True:
                            print('FSL DECOD',cnn_corr_tmp.get_shape().as_list())
                            #cnn_corr_tmp=tf.reshape(tf.transpose(tf.reshape(tf.tile(cnn_corr_tmp,norient),[128,npidx])),[128*npidx])
                            print('FSL DECOD',cnn_corr_tmp.get_shape().as_list())
                            #vdestrip1 = tf.reshape(tf.gather(cnn_corr_tmp,rotidx[k]),[12*32*32,128])
                            vdestrip1 = tf.transpose(tf.reshape(tf.gather(cnn_corr_tmp,rotidx[k]),[128,12*32*32]))
                            print('FSL DECOD',vdestrip1.get_shape().as_list(),cnn_imref[k].get_shape().as_list())
                            vdestrip1 = tf.matmul(cnn_imref[k],vdestrip1)
                            if DOREFERENCE==True:
                                vdestrip1= vdestrip1+ 1E-3 *theoffset[k]
                            print('FSL DECOD',vdestrip1.get_shape().as_list())
                            l_size=vdestrip1.get_shape().as_list()
                            print('FSL DECOD',l_size,cnn_corr_tmp.get_shape().as_list(),rotidx[k].get_shape().as_list())
                            cnn_corr[k]=tf.reshape(vdestrip1,[l_size[0]*l_size[1]])
                        else:
                            cnn_corr[k]=cnn_corr_tmp
                            
                        if rank==0:
                           print(circular[k])	
        else:
            icorrection=np.zeros([1],dtype='float')

        Cfunc.fsetbolo(buffer,theadd,data[ib].astype('float'),hdata[ib].astype('float'),hdata2[ib].astype('float'),
                       the_data,the_func,thecalib,hidx[ib].astype('int32'), ridx[ib].astype('int32'),icorrection.astype('int32'),
                       ib,ndata[ib]) 

        del add_data[ib]
        del data[ib]
        del hdata[ib]
        del hdata2[ib]
        del tmp_data[ib]
        del ridx[ib]
        del icorrection

    if ncnn>0:
        tmp_hh=np.zeros([mndata],dtype='float32')
        tmp_pidx=np.zeros([mndata],dtype='int32')
        Cfunc.fgetdataidx(buffer,tmp_hh,tmp_pidx)

        if rank%NUMBEROFTHREAD==0:
            l_train_hh={}
            l_train_hh[0]=1*tmp_hh
            for irank in range(1,NUMBEROFTHREAD):
                ncomm.Recv([tmp_hh, mndata, MPI.FLOAT],   irank, 100)
                l_train_hh[irank]   = 1*tmp_hh
        else:
            ncomm.Send([tmp_hh, mndata, MPI.FLOAT],   0, 100)

        hidx_cnn=val_pidx[tmp_pidx]
        # Compress healpix index to the one used only
        fpix_cnn = np.zeros([12*l_Nside_cnn*l_Nside_cnn],dtype='int')
        hidx_cnn = ChangeReso(hidx_cnn,Nside,l_Nside_cnn)
        fpix_cnn[hidx_cnn]=1

        val_pidx_cnn=np.where(fpix_cnn>0)[0]
        npixel_cnn=len(val_pidx_cnn)
        fpix_cnn[val_pidx_cnn]=np.arange(npixel_cnn)
        hidx_cnn=fpix_cnn[hidx_cnn]
        cnn_npixel=(hidx_cnn.max()+1)
        print("Rank %d [%d]: NPIXEL_CNN = %d "%(rank,ncnn,npixel_cnn))

        if NUMBEROFTHREAD!=1:
            vv=np.array([cnn_npixel],dtype='int32')
            lvv=np.array([0],dtype='int32')

            ncomm.Allreduce([vv,1,MPI.INT],lvv,op=MPI.MAX)
            l_cnn_npixel=lvv[0]
        else:
            l_cnn_npixel=cnn_npixel

        ntotring=0
        for ib in range(nbolo):
            ntotring+=nring[ib]

        if rank%NUMBEROFTHREAD==0:
            l_train_pidx={}
            l_train_pidx[0]=1*hidx_cnn
            for irank in range(1,NUMBEROFTHREAD):
                ncomm.Recv([hidx_cnn, mndata, MPI.INT],   irank, 100)
                l_train_pidx[irank]   = 1*hidx_cnn
        else:
            ncomm.Send([hidx_cnn, mndata, MPI.INT],   0, 100)

        del tmp_hh
        del tmp_pidx
        del hidx_cnn

        for k in range(ncnn):
            tmp=np.zeros([mndata],dtype='int32')
            tmp2=np.zeros([mndata],dtype='float32')
            
            Cfunc.fgetidxdata(buffer,tmp,tmp2,k)
            
            if rank%NUMBEROFTHREAD==0:
                l_train_ridx[k]={}
                l_train_wridx[k]={}
                print('RIDX ',tmp)
                print('WRIDX ',tmp2)
                l_train_ridx[k][0]=1*tmp
                l_train_wridx[k][0]=AMPDECOD[k]*tmp2
                for irank in range(1,NUMBEROFTHREAD):
                    ltmp=np.zeros([mndata],dtype='int32')
                    ltmp2=np.zeros([mndata],dtype='float32')
                    ncomm.Recv([ltmp, mndata, MPI.INT],   irank, 100)
                    ncomm.Recv([ltmp2, mndata, MPI.FLOAT],   irank, 101)
                    l_train_ridx[k][irank]   = 1*ltmp
                    l_train_wridx[k][irank]  = AMPDECOD[k]*ltmp2
                    del ltmp
                    del ltmp2
            else:
                ncomm.Send([tmp, mndata, MPI.INT],   0, 100)
                ncomm.Send([tmp2, mndata, MPI.FLOAT],  0, 101)
            del tmp
            del tmp2

            


        if rank%NUMBEROFTHREAD==0:
            if DOFULLYCONNECTED:
                logits,nlogits,lavv,lscat,lscat0,lscat1,lscat2,lfsl = cnn_model(train_data,param,train_hh,train_pidx,
                                                                           train_ridx,train_wridx,l_cnn_npixel,
                                                                           fc2_weights,
                                                                           fc2_biases,
                                                                           dconv_weights,
                                                                           dconv_biases,circular)
            else:
                logits,nlogits,lavv,lscat,lscat0,lscat1,lscat2,lfsl  = cnn_model_noFL(train_data,param,train_hh,train_pidx,
                                                                                 train_ridx,train_wridx,l_cnn_npixel,
                                                                                 dconv_weights,
                                                                                 dconv_biases,circular)

            loss = logits + WLoss*(nlogits + lavv)

            batch = tf.Variable(0, dtype=tf.float32)

            # Decay once per epoch, using an exponential schedule starting at 0.01.
            learning_rate = tf.compat.v1.train.exponential_decay(
                LEARNING_RATE,       # Base learning rate.
                batch,  # Current index into the dataset.
                10,          # Decay step.
                DECAY_RATE,          # Decay rate.
                staircase=True)

            # Use simple momentum for the optimization.
            opti=tf.compat.v1.train.AdamOptimizer(learning_rate,0.9)

            optimizer = opti.minimize(loss,global_step=batch)

            gradient={}
            igrad={}

            apply_transform_op={}

            for icnn in range(ncnn):
                igrad[icnn]={}
                if PARAMETERONLY==True:
                    gradient[icnn] = opti.compute_gradients(loss,var_list=[param[icnn]]) #,dconv_weights[icnn][2],dconv_biases[icnn][2]])
                    
                    if DOFULLYCONNECTED:
                        if isFSL[icnn]==False:
                            if is2D[icnn]==True:
                                igrad[icnn][0]=tf.compat.v1.placeholder(tf.float32,shape=(1,NHIDDEN[icnn])) #par
                            else:
                                igrad[icnn][0]=tf.compat.v1.placeholder(tf.float32,shape=(cnn_ntotring[icnn],NHIDDEN[icnn])) #par
                        else:
                            igrad[icnn][0]=tf.compat.v1.placeholder(tf.float32,shape=(1,NHIDDEN[icnn])) #par
                    else:
                        if is2D[icnn]==False:
                            if isFSL[icnn]==False:
                                igrad[icnn][0]=tf.compat.v1.placeholder(tf.float32,shape=(cnn_ntotring[icnn],
                                                                                          NUM_DCONV_CHAN[k][0]*XIMAGE_SIZE[k]//(WSCALE**(NDCONV[k])))) #param
                            else:
                                igrad[icnn][0]=tf.compat.v1.placeholder(tf.float32,shape=(1,
                                                                                          NUM_DCONV_CHAN[k][0]*12*32*32//(WSCALE**(NDCONV[k])))) #param
                        else:
                            igrad[icnn][0]=tf.compat.v1.placeholder(tf.float32,shape=(YIMAGE_SIZE[k]//(WSCALE**(NDCONV[k])),
                                                                                      NUM_DCONV_CHAN[k][0]*XIMAGE_SIZE[k]//(WSCALE**(NDCONV[k])))) #param
                    #if is2D[icnn]==False:
                    #    igrad[icnn][1]=tf.compat.v1.placeholder(tf.float32,shape=(KERNELSZ[icnn], 1,NUM_DCONV_CHAN[icnn][3], NUM_DCONV_CHAN[icnn][2])) #dconv_weights[2]
                    #else:
                    #    igrad[icnn][1]=tf.compat.v1.placeholder(tf.float32,shape=(KERNELSZ[icnn], KERNELSZ[icnn],NUM_DCONV_CHAN[icnn][3], NUM_DCONV_CHAN[icnn][2])) #dconv_weights[2]
                    #igrad[icnn][2]=tf.compat.v1.placeholder(tf.float32,shape=(NUM_DCONV_CHAN[icnn][3])) #dconv_biases[2]
                    
                else:
                    if DOFULLYCONNECTED:
                        if DOREFERENCE==True:
                            gradient[icnn] = opti.compute_gradients(loss,var_list=[param[icnn],fc2_weights[icnn],fc2_biases[icnn],
                                                                                   dconv_weights[icnn][0],dconv_weights[icnn][1],dconv_weights[icnn][2],
                                                                                   dconv_biases[icnn][0],dconv_biases[icnn][1],dconv_biases[icnn][2],theoffset[icnn]])
                            igrad[icnn][9]=tf.compat.v1.placeholder(tf.float32,shape=(1)) #theoffset
                        else:
                            gradient[icnn] = opti.compute_gradients(loss,var_list=[param[icnn],fc2_weights[icnn],fc2_biases[icnn],
                                                                                   dconv_weights[icnn][0],dconv_weights[icnn][1],dconv_weights[icnn][2],
                                                                                   dconv_biases[icnn][0],dconv_biases[icnn][1],dconv_biases[icnn][2]])
                        if is2D[icnn]==False:
                            if isFSL[icnn]==False:
                                igrad[icnn][7]=tf.compat.v1.placeholder(tf.float32,shape=(NHIDDEN[icnn],
                                                                                          XIMAGE_SIZE[icnn] // (WSCALE**NDCONV[icnn]) *  NUM_DCONV_CHAN[icnn][0])) #fc2_weights
                                igrad[icnn][8]=tf.compat.v1.placeholder(tf.float32,shape=(XIMAGE_SIZE[icnn] // (WSCALE**NDCONV[icnn])  * NUM_DCONV_CHAN[icnn][0])) #fc2_biases
                            else:
                                igrad[icnn][7]=tf.compat.v1.placeholder(tf.float32,shape=(NHIDDEN[icnn],
                                                                                          12*32*32 // (WSCALE**NDCONV[icnn]) *  NUM_DCONV_CHAN[icnn][0])) #fc2_weights
                                igrad[icnn][8]=tf.compat.v1.placeholder(tf.float32,shape=(12*32*32 // (WSCALE**NDCONV[icnn])  * NUM_DCONV_CHAN[icnn][0])) #fc2_biases
                        else:
                            igrad[icnn][7]=tf.compat.v1.placeholder(tf.float32,shape=(NHIDDEN[icnn],XIMAGE_SIZE[icnn] // (WSCALE**NDCONV[icnn]) *
                                                                                      YIMAGE_SIZE[icnn] // (WSCALE**NDCONV[icnn]) *  NUM_DCONV_CHAN[icnn][0])) #fc2_weights
                            igrad[icnn][8]=tf.compat.v1.placeholder(tf.float32,shape=(XIMAGE_SIZE[icnn] // (WSCALE**NDCONV[icnn])  *
                                                                                      YIMAGE_SIZE[icnn] // (WSCALE**NDCONV[icnn])  * NUM_DCONV_CHAN[icnn][0])) #fc2_biases

                    else:
                        if DOREFERENCE==True:
                            gradient[icnn] = opti.compute_gradients(loss,var_list=[param[icnn],
                                                                                   dconv_weights[icnn][0],dconv_weights[icnn][1],dconv_weights[icnn][2],
                                                                                   dconv_biases[icnn][0],dconv_biases[icnn][1],dconv_biases[icnn][2],theoffset[icnn]])
                            igrad[icnn][7]=tf.compat.v1.placeholder(tf.float32,shape=(1)) #theoffset
                        else:
                            gradient[icnn] = opti.compute_gradients(loss,var_list=[param[icnn],
                                                                                   dconv_weights[icnn][0],dconv_weights[icnn][1],dconv_weights[icnn][2],
                                                                                   dconv_biases[icnn][0],dconv_biases[icnn][1],dconv_biases[icnn][2]])

                    if DOFULLYCONNECTED:
                        igrad[icnn][0]=tf.compat.v1.placeholder(tf.float32,shape=(cnn_ntotring[icnn],NHIDDEN[icnn])) #param
                    else:
                        if is2D[icnn]==False:
                            if isFSL[icnn]==False:
                                igrad[icnn][0]=tf.compat.v1.placeholder(tf.float32,shape=(cnn_ntotring[icnn],
                                                                                          NUM_DCONV_CHAN[k][0]*XIMAGE_SIZE[k]//(WSCALE**(NDCONV[k])))) #param
                            else:
                                igrad[icnn][0]=tf.compat.v1.placeholder(tf.float32,shape=(1,
                                                                                          NUM_DCONV_CHAN[k][0]*12*32*32//(WSCALE**(NDCONV[k])))) #param
                        else:
                            igrad[icnn][0]=tf.compat.v1.placeholder(tf.float32,shape=(YIMAGE_SIZE[k]//(WSCALE**(NDCONV[k])),
                                                                                      NUM_DCONV_CHAN[k][0]*XIMAGE_SIZE[k]//(WSCALE**(NDCONV[k])))) #param


                    if is2D[icnn]==False:
                        igrad[icnn][1]=tf.compat.v1.placeholder(tf.float32,shape=(KERNELSZ[icnn], 1,NUM_DCONV_CHAN[icnn][1], NUM_DCONV_CHAN[icnn][0])) #dconv_weights[0]
                        igrad[icnn][2]=tf.compat.v1.placeholder(tf.float32,shape=(KERNELSZ[icnn], 1,NUM_DCONV_CHAN[icnn][2], NUM_DCONV_CHAN[icnn][1])) #dconv_weights[1]
                        igrad[icnn][3]=tf.compat.v1.placeholder(tf.float32,shape=(KERNELSZ[icnn], 1,NUM_DCONV_CHAN[icnn][3], NUM_DCONV_CHAN[icnn][2])) #dconv_weights[2]
                    else:
                        igrad[icnn][1]=tf.compat.v1.placeholder(tf.float32,shape=(KERNELSZ[icnn], KERNELSZ[icnn],NUM_DCONV_CHAN[icnn][1], NUM_DCONV_CHAN[icnn][0])) #dconv_weights[0]
                        igrad[icnn][2]=tf.compat.v1.placeholder(tf.float32,shape=(KERNELSZ[icnn], KERNELSZ[icnn],NUM_DCONV_CHAN[icnn][2], NUM_DCONV_CHAN[icnn][1])) #dconv_weights[1]
                        igrad[icnn][3]=tf.compat.v1.placeholder(tf.float32,shape=(KERNELSZ[icnn], KERNELSZ[icnn],NUM_DCONV_CHAN[icnn][3], NUM_DCONV_CHAN[icnn][2])) #dconv_weights[2]
                    igrad[icnn][4]=tf.compat.v1.placeholder(tf.float32,shape=(NUM_DCONV_CHAN[icnn][1])) #dconv_biases[0]
                    igrad[icnn][5]=tf.compat.v1.placeholder(tf.float32,shape=(NUM_DCONV_CHAN[icnn][2])) #dconv_biases[1]
                    igrad[icnn][6]=tf.compat.v1.placeholder(tf.float32,shape=(NUM_DCONV_CHAN[icnn][3])) #dconv_biases[2]

            tgradient={}
            for icnn in range(ncnn):
                tgradient[icnn] = [(igrad[icnn][i],GeV[1]) for (i,GeV) in enumerate(gradient[icnn])] #list [(grad,var)]
                apply_transform_op[icnn] = opti.apply_gradients(tgradient[icnn],global_step=batch)
            if rank==0:
                print('Initialized!')

    comm.Barrier()

    process = psutil.Process(os.getpid())
    if (rank%16==0):
        print('MEM %s %d %.3f MB'%(getframeinfo(currentframe()).lineno,rank,(process.memory_info().rss)/(1024.*1024.)),getloadavg())

    nfit_bolo=[(ntmp+nring[i]+ntotMfunc) for i in nring]
    nfit=np.array(nfit_bolo).sum()

    #=============================================================================================================================
    #                                   COMPUTE INDEX NORM
    #=============================================================================================================================
    if normrelative:
        Cfunc.frelativ(buffer,frelat.astype('int32'))

    ntotdata=0
    ntotring=0
    for ib in range(nbolo):
        ntotdata+=ndata[ib]
        ntotring+=nring[ib]

    Cfunc.finitgrad(buffer)
    A0=np.zeros([1],dtype='float')
    LA0=np.zeros([1],dtype='float')
    Cfunc.calc_A0(buffer,LA0)
    comm.Allreduce((LA0,MPI.DOUBLE),(A0,MPI.DOUBLE))

    gain=np.ones([nbolo])


    # Create a local session to run the training.
    start_time = time.time()

    xk=np.zeros([nfit],dtype='float')
    nnn=0
    for ib in range(nbolo):
        xk[nnn:nring[ib]+nnn]=np.cos(np.arange(nring[ib])/1000.)
        xk[nnn:nring[ib]+nnn]-=xk[nnn:nring[ib]+nnn].mean()
        nnn+=nring[ib]
    
    if rank%NUMBEROFTHREAD==0:
        # chek if tensorflow is a GPU based one
        if 'gpu_options' in dir(tf):
            gpu_options = tf.GPUOptions(per_process_gpu_memory_fraction=(1.0*size)/NUMBEROFTHREAD)
            sess=tf.compat.v1.Session(config=tf.ConfigProto(inter_op_parallelism_threads=np.max([1,NUMBEROFTHREAD//2]),
                                                            intra_op_parallelism_threads=np.max([1,NUMBEROFTHREAD//2]),
                                                            gpu_options=gpu_options)) 
        else:
            sess=tf.compat.v1.Session(config=tf.ConfigProto(inter_op_parallelism_threads=np.max([1,NUMBEROFTHREAD//2]),
                                                            intra_op_parallelism_threads=np.max([1,NUMBEROFTHREAD//2])))

        # Run all the initializers to prepare the trainable parameters.
        #tf.global_variables_initializer().run()
        tf.compat.v1.global_variables_initializer().run(session=sess)
        print('Rank %d Tensorflow Initialized!'%(rank))

    if ncnn>0:
        if is2D[0]==True and isFSL[0]==False:
            init_decod=np.zeros([YIMAGE_SIZE[0],XIMAGE_SIZE[0]],dtype='float')
        else:
            init_decod=np.zeros([cnn_ntotring[0],XIMAGE_SIZE[0]],dtype='float')

    for itt in range(NITT):
        b=np.zeros([nfit],dtype='float')
        y=np.zeros([nfit],dtype='float')
        rk=np.zeros([nfit],dtype='float')
        Apk=np.zeros([nfit],dtype='float')

        if rank==0:
            print('=========================================')
            print('==     ITT = %4d                      =='%(itt))
            print('=========================================')
            print('GAIN',gain)
            print('A0',A0[0])
            print('NFIT',nfit)

        # DO NOT COMPUTE THE CALIBRATION THAT COULD BE DEGENERATED WITH TEMPLATE FITTING
        if itt>NITTCALIB and docalib==1:
            Cfunc.nocalib(buffer,icalib)
        else:
            Cfunc.docalib(buffer)
            
        if ncnn>0 and (itt>0 or PARAMETERONLY==True):
            for k in range(ncnn):
                if rank%NUMBEROFTHREAD==0:
                    #decod = sess.run([cnn_corr[k]])[0]
                    decod,par,lf1 = sess.run([cnn_corr[k],param[k],lfsl])
                    #decod = init_decod.astype('float32')
                else:
                    if is2D[k]==True and isFSL[k]==False:
                        decod=np.zeros([YIMAGE_SIZE[k]*XIMAGE_SIZE[k]],dtype='float32')
                    else:
                        decod=np.zeros([cnn_ntotring[k]*XIMAGE_SIZE[k]],dtype='float32')
    
                ncomm.Bcast(decod,root=0)
                if CNNTEMPLATE[k]!=-1:
                    
                    coef2[CNNTEMPLATE[k]]=1.0
                    #if itt!=0:
                    #    decod*=xk[ntotring+CNNTEMPLATE[k]]
                    
                Cfunc.fsetcorrection(buffer, AMPDECOD[k]*decod.astype('float'),int(k),int(CNNTEMPLATE[k]))
                
                #f DOPLOT==True and rank%NUMBEROFTHREAD==0 and itt==0:
                #   print(par,lpidx.shape,lf1.shape)
                #   plt.subplot(1,2,1)
                #   if isFSL[0]==True:
                #       lpidx=np.load('healpix2ap_32_nest.npy')
                #       plt.imshow(lf1[lpidx].reshape(128,128),cmap='jet')
                #   plt.subplot(1,2,2)
                #   plt.imshow(AMP2DECOD[k]*decod.reshape(cnn_ntotring[k],XIMAGE_SIZE[k]),cmap='jet')
                #   plt.show()
                #
                #   print('DECOD ',decod.std())

        Cfunc.fsetbolo_gain(buffer,gain.astype('float'),1,1)

        Cfunc.projdata(buffer,A0[0],b)

        Cfunc.proj(buffer,xk,y)

        lrk=np.zeros([nfit],dtype='float')
        lrk[:]=b-y
        
        comm.Allreduce((lrk,MPI.DOUBLE),(rk,MPI.DOUBLE))

        pk=1.0*rk
        deltak=(rk**2).sum()
        if itt==0:
            delta0=1.0*deltak
        delta=1.0*deltak
        step=0
        while(step<NUM_EPOCHS and delta>1E-20*delta0):

            delta=1.0*deltak
            
            lApk=np.zeros([nfit],dtype='float')
            
            Cfunc.proj(buffer,pk,lApk)

            comm.Allreduce((lApk,MPI.DOUBLE),(Apk,MPI.DOUBLE))

            alpha=delta/(pk*Apk).sum()

            xk += alpha*pk

            rk -= alpha*Apk
            
            deltak=(rk**2).sum()

            beta = deltak/delta

            pk = 1.0*(rk + beta*pk)
            
            if step%EVAL_FREQUENCY==0 or delta<=1E-20*delta0:
                if ncnn>0 and step!=0 and itt!=NITT-1:
                    sig=np.zeros([mndata],dtype='float')

                    Cfunc.fsetbolo_gain(buffer,gain.astype('float'),1,0)

                    #xk[:]=0.0
                    if DOREFERENCE==True or dotemplate==True:
                        test=1
                    else:
                        test=0
                        
                    Cfunc.fgetdata(buffer,xk,sig,init_decod,0,int(test),CNNTEMPLATE[0])
                     #if rank==0:
                    #   test=1
                    #   Cfunc.fgetdata(buffer,xk,sig,init_decod,0,int(test),CNNTEMPLATE[0])
                    #   np.save(Out_Offset[0]+'_INITDECOD',init_decod)
                    #   print('INITDECOD ',init_decod.std())
                    #   print(np.polyfit(decod.flatten(),init_decod.flatten(),1))
                    #   plt.subplot(1,3,1)
                    #   plt.imshow(init_decod,cmap='jet')
                    #   nx,ny=init_decod.shape
                    #   plt.subplot(1,3,2)
                    #   plt.imshow(AMPDECOD[0]*decod.reshape(nx,ny),cmap='jet')
                    #   plt.subplot(1,3,3)
                    #   plt.imshow(init_decod-AMPDECOD[0]*decod.reshape(nx,ny),cmap='jet')
                    #   plt.show()
                    #comm.Barrier()
                    #exit(0)
                    
                    if dotemplate==True:
                        if rank==0:
                            print('==========================================')
                            print('save TEMPLATE : %s_TEMPLATE.npy'%(Out_Offset[0]))
                            np.save('%s_TEMPLATE.npy'%(Out_Offset[0]),init_decod)
                            print('=                 DONE                   =')
                            print('==========================================')
                        exit(0)
                    if rank%NUMBEROFTHREAD==0:
                        l_sig={}
                        l_sig[0]=1*sig
                        for irank in range(1,NUMBEROFTHREAD):
                            ncomm.Recv([sig, mndata, MPI.DOUBLE],   irank, 100)
                            l_sig[irank]   = 1*sig
                    else:
                        ncomm.Send([sig, mndata, MPI.DOUBLE],   0, 100)

                    dt1=0
                    dt1b=0
                    dt2=0

                    if rank%NUMBEROFTHREAD==0:
                        for itt_cnn in range(NITT_CNN):
                                    
                            feed_dict = {}
                            if CNNTEMPLATE[0]!=-1:
                                vv=np.array([xk[ntotring+CNNTEMPLATE[0]]],dtype='float32')
                                feed_dict[train_var2fit] =  vv
                            if BATCH_SIZE==-1:
                                feed_dict[train_data] = l_sig[itt_cnn%NUMBEROFTHREAD]
                                feed_dict[train_hh]   = l_train_hh[itt_cnn%NUMBEROFTHREAD]
                                feed_dict[train_pidx] = l_train_pidx[itt_cnn%NUMBEROFTHREAD]
                                for icnn in range(ncnn):
                                    feed_dict[train_ridx[icnn]] = l_train_ridx[icnn][itt_cnn%NUMBEROFTHREAD]
                                    feed_dict[train_wridx[icnn]] = l_train_wridx[icnn][itt_cnn%NUMBEROFTHREAD]
                            else:
                                idx=(np.random.rand(BATCH_SIZE)*mndata).astype('int')
                                feed_dict[train_data] = l_sig[itt_cnn%NUMBEROFTHREAD][idx]
                                feed_dict[train_hh]   = l_train_hh[itt_cnn%NUMBEROFTHREAD][idx]
                                feed_dict[train_pidx] = l_train_pidx[itt_cnn%NUMBEROFTHREAD][idx]
                                for icnn in range(ncnn):
                                    feed_dict[train_ridx[icnn]] = l_train_ridx[icnn][itt_cnn%NUMBEROFTHREAD][idx]
                                    feed_dict[train_wridx[icnn]] = l_train_wridx[icnn][itt_cnn%NUMBEROFTHREAD][idx]

                            l={}
                            vgrad={}
                            for icnn in range(ncnn):
                                start_time1 = time.time()
                                l[icnn],vgrad[icnn] = sess.run([loss,gradient[icnn]], feed_dict=feed_dict)
                                dt1-= start_time1 - time.time()
                                
                            feed_dict2 = {}
                            for icnn in range(ncnn):
                                start_time1 = time.time()
                                rval={}
                                #====================================================
                                # MERGE ALL GRADIENT FROM ALL RANK
                                #====================================================

                                for i in range(len(vgrad[icnn])):
                                    g,v=vgrad[icnn][i]
                                    lg=np.zeros([len(g.flatten())],dtype='float32')
                                    ncomm_transpose.Allreduce((g.flatten(),MPI.FLOAT),(lg,MPI.FLOAT))
                                    rval[i] = lg
                                    
                                del l[icnn],vgrad[icnn]
                                
                                if PARAMETERONLY==True:
                                    if DOFULLYCONNECTED:
                                        feed_dict2[igrad[icnn][0]] =rval[0].reshape(1,NHIDDEN[icnn]) #param
                                        
                                        if is2D[k]==False:
                                            if isFSL[k]==False:
                                                feed_dict2[igrad[icnn][0]] =rval[0].reshape(cnn_ntotring[icnn],NHIDDEN[icnn]) #param
                                            else:
                                                feed_dict2[igrad[icnn][0]] =rval[0].reshape(1,NHIDDEN[icnn]) #param

                                        else:
                                            feed_dict2[igrad[icnn][0]] =rval[0].reshape(1,NHIDDEN[icnn]) #param
                                        #if is2D[k]==False:
                                        #    feed_dict2[igrad[icnn][1]] =rval[1].reshape(KERNELSZ[icnn], 1,NUM_DCONV_CHAN[icnn][3], NUM_DCONV_CHAN[icnn][2]) #dconv_weights[2]
                                        #else:
                                        #    feed_dict2[igrad[icnn][1]] =rval[1].reshape(KERNELSZ[icnn], KERNELSZ[icnn],NUM_DCONV_CHAN[icnn][3], NUM_DCONV_CHAN[icnn][2]) #dconv_weights[2]
                                        #feed_dict2[igrad[icnn][2]] =rval[2].reshape(NUM_DCONV_CHAN[icnn][3]) #dconv_biases[2]
                                    else:
                                        if is2D[k]==False:
                                            if isFSL[k]==False:
                                                feed_dict2[igrad[icnn][0]] =rval[0].reshape(cnn_ntotring[icnn],NUM_DCONV_CHAN[icnn][0]*XIMAGE_SIZE[icnn]//(WSCALE**(NDCONV[icnn]))) #param
                                            else:
                                                feed_dict2[igrad[icnn][0]] =rval[0].reshape(1,NUM_DCONV_CHAN[icnn][0]*12*32*32//(WSCALE**(NDCONV[icnn]))) #param
                                        else:
                                            feed_dict2[igrad[icnn][0]] =rval[0].reshape(YIMAGE_SIZE[icnn] // (WSCALE**NDCONV[icnn]),
                                                                                        NUM_DCONV_CHAN[icnn][0]*XIMAGE_SIZE[icnn]//(WSCALE**(NDCONV[icnn]))) #param
                                        #feed_dict2[igrad[icnn][1]] =rval[1].reshape(1) #fc2_biases

                                else:
                                    if DOFULLYCONNECTED:
                                        if is2D[k]==False:
                                            if isFSL[k]==False:
                                                feed_dict2[igrad[icnn][0]] =rval[0].reshape(cnn_ntotring[icnn],NHIDDEN[icnn]) #param
                                                feed_dict2[igrad[icnn][7]] =rval[7].reshape(NHIDDEN[icnn], XIMAGE_SIZE[icnn] // (WSCALE**NDCONV[icnn]) *  NUM_DCONV_CHAN[icnn][0]) #fc2_weights
                                                feed_dict2[igrad[icnn][8]] =rval[8].reshape(XIMAGE_SIZE[icnn]//(WSCALE**NDCONV[icnn])*NUM_DCONV_CHAN[icnn][0]) #fc2_biases
                                            else:
                                                feed_dict2[igrad[icnn][0]] =rval[0].reshape(1,NHIDDEN[icnn]) #param
                                                feed_dict2[igrad[icnn][7]] =rval[7].reshape(NHIDDEN[icnn], 12*32*32 // (WSCALE**NDCONV[icnn]) *  NUM_DCONV_CHAN[icnn][0]) #fc2_weights
                                                feed_dict2[igrad[icnn][8]] =rval[8].reshape(12*32*32//(WSCALE**NDCONV[icnn])*NUM_DCONV_CHAN[icnn][0]) #fc2_biases
                                        else:
                                            feed_dict2[igrad[icnn][0]] =rval[0].reshape(YIMAGE_SIZE[icnn] // (WSCALE**NDCONV[icnn]),NHIDDEN[0]) #param
                                            feed_dict2[igrad[icnn][7]] =rval[7].reshape(NHIDDEN[icnn], YIMAGE_SIZE[icnn] // (WSCALE**NDCONV[icnn]) * XIMAGE_SIZE[icnn] // (WSCALE**NDCONV[icnn]) *  NUM_DCONV_CHAN[icnn][0]) #fc2_weights
                                            feed_dict2[igrad[icnn][8]] =rval[8].reshape(YIMAGE_SIZE[icnn] // (WSCALE**NDCONV[icnn])* XIMAGE_SIZE[icnn]//(WSCALE**NDCONV[icnn])*NUM_DCONV_CHAN[icnn][0]) #fc2_biases
                                        if DOREFERENCE==True:
                                            feed_dict2[igrad[icnn][9]] =rval[9].reshape(1) # theoffset
                                    else:
                                        if is2D[k]==False:
                                            if isFSL[k]==False:
                                                feed_dict2[igrad[icnn][0]] =rval[0].reshape(cnn_ntotring[icnn],NUM_DCONV_CHAN[icnn][0]*XIMAGE_SIZE[icnn]//(WSCALE**(NDCONV[icnn]))) #param
                                            else:
                                                feed_dict2[igrad[icnn][0]] =rval[0].reshape(1,NUM_DCONV_CHAN[icnn][0]*12*32*32//(WSCALE**(NDCONV[icnn]))) #param
                                        else:
                                            feed_dict2[igrad[icnn][0]] =rval[0].reshape(YIMAGE_SIZE[icnn] // (WSCALE**NDCONV[icnn]),
                                                                                        NUM_DCONV_CHAN[icnn][0]*XIMAGE_SIZE[icnn]//(WSCALE**(NDCONV[icnn]))) #param
                                        if DOREFERENCE==True:
                                            feed_dict2[igrad[icnn][7]] =rval[7].reshape(1) # theoffset

                                    if is2D[k]==False:
                                        feed_dict2[igrad[icnn][1]] =rval[1].reshape(KERNELSZ[icnn], 1,NUM_DCONV_CHAN[icnn][1], NUM_DCONV_CHAN[icnn][0]) #dconv_weights[0]
                                        feed_dict2[igrad[icnn][2]] =rval[2].reshape(KERNELSZ[icnn], 1,NUM_DCONV_CHAN[icnn][2], NUM_DCONV_CHAN[icnn][1]) #dconv_weights[1]
                                        feed_dict2[igrad[icnn][3]] =rval[3].reshape(KERNELSZ[icnn], 1,NUM_DCONV_CHAN[icnn][3], NUM_DCONV_CHAN[icnn][2]) #dconv_weights[2]
                                    else:
                                        feed_dict2[igrad[icnn][1]] =rval[1].reshape(KERNELSZ[icnn], KERNELSZ[icnn],NUM_DCONV_CHAN[icnn][1], NUM_DCONV_CHAN[icnn][0]) #dconv_weights[0]
                                        feed_dict2[igrad[icnn][2]] =rval[2].reshape(KERNELSZ[icnn], KERNELSZ[icnn],NUM_DCONV_CHAN[icnn][2], NUM_DCONV_CHAN[icnn][1]) #dconv_weights[1]
                                        feed_dict2[igrad[icnn][3]] =rval[3].reshape(KERNELSZ[icnn], KERNELSZ[icnn],NUM_DCONV_CHAN[icnn][3], NUM_DCONV_CHAN[icnn][2]) #dconv_weights[2]
                                    feed_dict2[igrad[icnn][4]] =rval[4].reshape(NUM_DCONV_CHAN[icnn][1]) #dconv_biases[0]
                                    feed_dict2[igrad[icnn][5]] =rval[5].reshape(NUM_DCONV_CHAN[icnn][2]) #dconv_biases[1]
                                    feed_dict2[igrad[icnn][6]] =rval[6].reshape(NUM_DCONV_CHAN[icnn][3]) #dconv_biases[2]
                                dt1b-= start_time1 - time.time() 
                                
                            start_time2 = time.time()       
                            l      = sess.run(loss, feed_dict=feed_dict)
                            for icnn in range(ncnn):
                                result = sess.run(apply_transform_op[icnn], feed_dict=feed_dict2)

                            # TODO ONCE SUCCEED TO GO BACK TO PREVIOUS STEP
                            #lnew   = sess.run(loss, feed_dict=feed_dict)
                                
                            #vv=np.array([l,lnew],dtype='float32')
                            #nvv=np.array([0.0,0.0],dtype='float32')
                            #ncomm_transpose.Allreduce((vv,MPI.FLOAT),(nvv,MPI.FLOAT))
                            #dvv=0.1
                            #while nvv[1]>1.2*nvv[0] and dvv>1E-3:
                            #    feed_dict3={}
                            #    for imod in feed_dict2:
                            #        feed_dict3[imod]=-(1-dvv)*feed_dict2[imod]
                            #        
                            #    for icnn in range(ncnn):
                            #        result = sess.run(apply_transform_op[icnn], feed_dict=feed_dict3)
                            #        
                            #    for imod in feed_dict2:
                            #        feed_dict2[imod]=dvv*feed_dict2[imod]
                            #    if rank==0:
                            #        print('CONVERGENCE PROBLEM : GRADIENT ADJUSTED ',itt_cnn,nvv[0],nvv[1],nvv[0]/nvv[1])
                            #            
                            #    lnew   = sess.run(loss, feed_dict=feed_dict)
                            #    vv=np.array([l,lnew],dtype='float32')
                            #    nvv=np.array([0.0,0.0],dtype='float32')
                            #    ncomm_transpose.Allreduce((vv,MPI.FLOAT),(nvv,MPI.FLOAT))
                            #    if rank==0:
                            #        print('CONVERGENCE CORRECTED : GRADIENT ADJUSTED ',itt_cnn,nvv[0],nvv[1],nvv[0]/nvv[1])
                            #    dvv/=10
                            
                            dt2-= start_time2 - time.time()   

                            if itt_cnn%PRT_FREQUENCY==0:

                                elapsed_time2 = time.time() - start_time
                                start_time = time.time()
                                if rank==0:
                                    start_time3 = time.time()    
                                    if PARAMETERONLY:   
                                        decod,par,meanamp,l1 = sess.run([cnn_corr[0],param[0],lavv,nlogits],feed_dict=feed_dict)
                                        v1=var2fit[0]
                                        dt3 = start_time3 - time.time()

                                        print('Itt %d-%d %8.5f s Loss = %.6g Par[0] = %.6g'%(itt,itt_cnn,elapsed_time2/PRT_FREQUENCY,l,par[0]),
                                              par.std(),v1,decod.std(),meanamp,(l1),
                                              getloadavg(),'%.2f,%.2f,%.2f,%.2f'%(dt1,dt1b,dt2,dt3))
                                    else:
                                        decod,par,meanamp,l1 = sess.run([cnn_corr[0],param[0],lavv,nlogits],feed_dict=feed_dict)
                                        dt3 = start_time3 - time.time()

                                        print('Itt %d-%d %8.5f s Loss = %.6g'%(itt,itt_cnn,elapsed_time2/PRT_FREQUENCY,l),
                                              par.std(),decod.std(),meanamp,(l1),
                                              getloadavg(),'%.2f,%.2f,%.2f,%.2f'%(dt1,dt1b,dt2,dt3))
                                        
                                    dt1=0
                                    dt1b=0
                                    dt2=0


                    for k in range(ncnn):
                        if rank%NUMBEROFTHREAD==0:
                            decod,par,ill,il0,il1,il2,ll,l0,l1,l2,lf1 = sess.run([cnn_corr[k],param[k],scatter_objective_cos[k],scatter_0_objective[0],scatter_1_objective[0],scatter_2_objective[0],lscat,lscat0,lscat1,lscat2,lfsl],feed_dict=feed_dict)
                            
                            #decod=init_decod.astype('float32')
                            print('DECOD ',decod.shape,cnn_ntotring[k]*XIMAGE_SIZE[k],cnn_ntotring[k],XIMAGE_SIZE[k],init_decod.shape)
                        else:
                            if is2D[k]==True and isFSL[k]==False:
                                decod=np.zeros([YIMAGE_SIZE[k]*XIMAGE_SIZE[k]],dtype='float32')
                            else:
                                decod=np.zeros([cnn_ntotring[k]*XIMAGE_SIZE[k]],dtype='float32')
                        comm.Bcast(decod,root=0)
                        #if CNNTEMPLATE[k]!=-1:
                        #    decod*=xk[ntotring+CNNTEMPLATE[k]]
                        Cfunc.fsetcorrection(buffer,AMPDECOD[k]*decod.astype('float'),int(k),int(CNNTEMPLATE[k]))
                        if rank==0:
                            np.save(Out_Offset[0]+'_DECOD_%d'%(k),decod)
                            if isFSL[k]==True:
                                np.save(Out_Offset[0]+'_FSL_%d'%(k),lf1)
                                
                        if DOPLOT==True:
                            if rank==0:
                                print(decod.std())
                                l_amp=3*AMPDECOD[k]*decod.std()
                                if NHARM[k]!=0:
                                    wfilt=np.dot(init_decod,co.transpose())
                                    print(wfilt.shape)
                                    wfilt=np.dot(co2,wfilt)
                                    print(wfilt.shape)
                                    wfilt=np.dot(wfilt,ico2)
                                    print(wfilt.shape)
                                    wfilt=(np.dot(ico.transpose(),wfilt)).transpose()
                                    print(wfilt.shape)
                                plt.figure(figsize=(6,12))
                                plt.subplot(1,3,1)
                                plt.imshow(init_decod-np.median(init_decod),cmap='jet',vmin=-l_amp,vmax=l_amp)
                                plt.subplot(1,3,2)
                                plt.imshow(AMPDECOD[k]*(decod.reshape(cnn_ntotring[k],XIMAGE_SIZE[k])-np.median(decod)),cmap='jet',vmin=-l_amp,vmax=l_amp)
                                #if NHARM[k]!=0:
                                #    plt.subplot(2,2,3)
                                #    plt.imshow(wfilt,cmap='jet',vmin=-l_amp,vmax=l_amp)
                                #if isFSL[k]==True:
                                #    plt.subplot(2,2,3)
                                #    if len(lf1)==128*128:
                                #        plt.imshow(lf1.reshape(128,128),cmap='jet')
                                #    else:
                                #        lpidx=np.load('healpix2ap_32_nest.npy')
                                #        plt.imshow(lf1[lpidx].reshape(128,128),cmap='jet')
                                        
                                plt.subplot(1,3,3)
                                tmp=init_decod-AMPDECOD[k]*decod.reshape(cnn_ntotring[k],XIMAGE_SIZE[k])
                                plt.imshow(tmp-np.median(tmp),cmap='jet',vmin=-l_amp,vmax=l_amp)
                                plt.show()

                    del sig

                    Cfunc.fsetbolo_gain(buffer,gain.astype('float'),1,1)

                    Cfunc.projdata(buffer,A0[0],b)

                    Cfunc.proj(buffer,xk,y)

                    lrk[:]=b-y

                    comm.Allreduce((lrk,MPI.DOUBLE),(rk,MPI.DOUBLE))

                    pk=1.0*rk


            if step%PRT_FREQUENCY==0:
                elapsed_time = time.time() - start_time
                start_time = time.time()
                if rank==0:
                    print('Step %d, %8.5f s Delta0: %.6g, delta: %.6g' % (step,elapsed_time/PRT_FREQUENCY,delta0,delta),getloadavg())

                    #print(["%10.4lg"%(xk[ntotring+ib+np.arange(nbolo,dtype='int')*ntmp].sum()) for ib in np.arange(ntmp,dtype='int')[frelat==1]])

            deltak=(rk**2).sum()
            step+=1


        vtmp=np.zeros([ntmp,nbolo])
        for ib in range(nbolo):
            if icalib[ib]!=-1:
                gain[ib]-=xk[ntotring+icalib[ib]+ib*ntmp]/coef2[icalib[ib]]

            vtmp[:,ib]=(xk[ntotring+ib*ntmp:ntotring+(ib+1)*ntmp])/coef2

        if rank==0:
            for ib in range(ntmp):
                print('TMP %s = '%(name_tmp[ib]),['%10.4lg'%(vtmp[ib,i]) for i in range(nbolo)],' MEAN = %10.4lg'%(vtmp[ib,:].mean()),' STD = %10.4lg'%(vtmp[ib,:].std()))
            for ib in range(nbolo):
                np.save(Out_Offset[ib]+'_VALUE',vtmp[:,ib])

        #=============================================================================================================================
        #                                  SAVE DATA
        #=============================================================================================================================
        if itt==NITT-1:

            if ncnn>0:
                if rank==0 and PARAMETERONLY==False:
                    for icnn in range(ncnn):
                        if DOFULLYCONNECTED:
                            fc2w,fc2b = sess.run([fc2_weights[icnn],fc2_biases[icnn]])
                            np.save('%s_%d_fc2w.npy'%(Out_Offset[0],icnn),fc2w)
                            np.save('%s_%d_fc2b.npy'%(Out_Offset[0],icnn),fc2b)
                            
                        for i in range(NDCONV[k]):
                            w,b=sess.run([dconv_weights[k][i],dconv_biases[k][i]])
                            np.save('%s_%d_w_%d.npy'%(Out_Offset[0],icnn,i),w)
                            np.save('%s_%d_b_%d.npy'%(Out_Offset[0],icnn,i),b)
                        par=sess.run([param[k]])
                        np.save('%s_%d_param.npy'%(Out_Offset[0],icnn),par[0])
                        
                        
                            
                for k in range(ncnn):
                    if rank%NUMBEROFTHREAD==0:
                        decod = sess.run([cnn_corr[k]])[0]
                        #decod = init_decod.astype('float32')
                    else:
                        decod=np.zeros([cnn_ntotring[k],XIMAGE_SIZE[k]],dtype='float32')
                        
                    comm.Bcast(decod,root=0)
                    
                    #Cfunc.fsetcorrection(buffer,decod.astype('float'),int(k),(CNNTEMPLATE[k]))

            if OUTCLEANCALIB==True:
                Cfunc.fsetbolo_gain(buffer,gain.astype('float'),1,1)
            else:
                Cfunc.fsetbolo_gain(buffer,gain.astype('float'),0,1)


            if rank==0:
                nn=0
                for ib in range(nbolo):
                    theoff=np.zeros([EndRing+1],dtype='float32')
                    theoff[:]=Cfunc.UNSEEN
                    theoff[val_ring[ib]]=xk[nn:nn+nring[ib]]
                    np.save(Out_Offset[ib],theoff)

                    nn+=nring[ib]

            for imap in range(len(MAP)):
                comm.Barrier()
                process = psutil.Process(os.getpid())
                if rank%16==0:
                    print('MEM %s %d %.3f MB'%(getframeinfo(currentframe()).lineno,rank,(process.memory_info().rss)/(1024.*1024.)),getloadavg())

                if rank==0:
                    map_fin=np.zeros([12*Nside*Nside],dtype='float32')
                    hmap_fin=np.zeros([12*Nside*Nside],dtype='float32')
                else:
                    map_fin=np.zeros([1],dtype='float32')
                    hmap_fin=np.zeros([1],dtype='float32')

                omap=np.zeros([npixel],dtype='float32')
                ohmap=np.zeros([npixel],dtype='float32')

                ffring=np.zeros([ntotring],dtype='int32')
                nn=0
                for ib in range(nbolo):
                    fidx=np.where((val_ring[ib]>=MAP[imap][3])*(val_ring[ib]<=MAP[imap][4]))[0]
                    ffring[nn+fidx]=1
                    nn+=nring[ib]

                Cfunc.domap(buffer,xk,omap,ohmap,np.array(MAP[imap][2],dtype='int32'),ffring)

                vv=np.array([npixel],dtype='int32')
                lvv=np.zeros([size],dtype='int32')

                comm.Allgather((vv,MPI.INT),(lvv,MPI.INT))

                if rank==0:
                    print('NPIXELS ',lvv.sum())
                    map_fin[val_pidx]=omap
                    hmap_fin[val_pidx]=ohmap
                    del omap
                    del ohmap

                    for i in range(1,size):
                        outidx=np.zeros([lvv[i]],dtype='int32')
                        outdata=np.zeros([lvv[i]],dtype='float32')

                        comm.Recv([outidx,  lvv[i], MPI.INT],   i, 100)
                        comm.Recv([outdata, lvv[i], MPI.FLOAT], i, 101)
                        map_fin[outidx]=outdata
                        comm.Recv([outdata, lvv[i], MPI.FLOAT], i, 102)
                        hmap_fin[outidx]=outdata

                        del outidx
                        del outdata
                else:
                        comm.Send([val_pidx.astype('int32'), npixel, MPI.INT], dest=0, tag=100)
                        comm.Send([omap.astype('float32'),   npixel, MPI.FLOAT], dest=0, tag=101)
                        comm.Send([ohmap.astype('float32'),  npixel, MPI.FLOAT], dest=0, tag=102)

                        del omap
                        del ohmap

                process = psutil.Process(os.getpid())
                if (rank%16==0):
                    print('MEM %s %d %.3f MB'%(getframeinfo(currentframe()).lineno,rank,(process.memory_info().rss)/(1024.*1024.)),getloadavg())

                if rank==0:
                    map_fin[hmap_fin>0]/=hmap_fin[hmap_fin>0]
                    map_fin[hmap_fin==0]=Cfunc.UNSEEN

                    name='!'+MAP[imap][0]+'_I.fits'
                    Cfunc.writemap(map_fin,Nside,name.encode('utf-8'),0,b'G')
                    name='!'+MAP[imap][0]+'_H.fits'
                    Cfunc.writemap(hmap_fin,Nside,name.encode('utf-8'),0,b'G')
                    print(MAP[imap][0],' Done')

                del map_fin
                del hmap_fin

            comm.Barrier()
        else:
            for ib in range(nbolo):
                if icalib[ib]!=-1:
                    xk[ntotring+icalib[ib]+ib*ntmp]=0.0


    if rank==0:
        print('TOTAL DURATION ',time.time() - start_time_total)

#=============================================================================================================================
#                                  MAIN : flag to be improved and used
#=============================================================================================================================
if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  
  parser.add_argument(
      '--dotemplate',
      default=False,
      help='Compute Scattering Transform Template',
      action='store_true')

  FLAGS, unparsed = parser.parse_known_args()
  main()

