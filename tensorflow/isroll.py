#!//opt/software/occigen/tools/python/2.7.12/intel/17.0/bin/python

# ==============================================================================
# 
#                      iSROLL : INVERT DATA WITH IA
# 
# ==============================================================================

"""iSROLL : INVERT DATA WITH IA

"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import matplotlib.pyplot as plt
import argparse
import os
import sys
import time

import numpy
from six.moves import xrange
import tensorflow as tf
import scipy.signal as sig
import sys
import healpy as hp
from mpi4py import MPI

#==============================================================================
#from mpi4py import MPI
#==============================================================================

param_file=sys.argv[1]
exec(open(param_file).read())

domodel=False
#==============================================================================
#
#               TEST ON Bflux (Does not work)
#
#==============================================================================

start_time_total = time.time()

def getmap(name,nside,dtype='float32',rgsz=27664,step=1):

    bytsz=sys.getsizeof(numpy.zeros([2],dtype=dtype))-sys.getsizeof(numpy.zeros([1],dtype=dtype))
    f = open(name, "rb")
    try:
        out = numpy.frombuffer(f.read(12*nside*nside*bytsz),dtype=dtype).reshape(12*nside*nside)
    finally:
        f.close()
    return(out)

def getdata(name,beginrg,endrg,dtype='float32',rgsz=27664,step=1):

    bytsz=sys.getsizeof(numpy.zeros([2],dtype=dtype))-sys.getsizeof(numpy.zeros([1],dtype=dtype))
    f = open(name, "rb")
    if step==1:
        try:
            f.seek(beginrg*rgsz*bytsz)
            out = numpy.frombuffer(f.read((endrg-beginrg+1)*rgsz*bytsz),dtype=dtype).reshape(endrg-beginrg+1,rgsz)
        finally:
            f.close()
    else:
        nring=(endrg-beginrg+1)//step
        out=numpy.zeros([nring,rgsz])
        try:
            for ir in range(nring):
                f.seek((beginrg+ir*step)*rgsz*bytsz)
                out[ir,:] = numpy.frombuffer(f.read((rgsz*bytsz)),dtype=dtype)
        finally:
            f.close()
    return(out)

def main(_):

    #=============================================================================================================================
    #                                    INIT MPI
    #=============================================================================================================================
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    #=============================================================================================================================
    #                                    COMPUTE BEGIN/ENDRG MPI
    #=============================================================================================================================
    NringPerProc = ((EndRing - BeginRing+1)//size)
    Nring2share = (EndRing - BeginRing+1)-NringPerProc*size
    RingTab = numpy.zeros([size],dtype='int32')+NringPerProc
    RingTab[0:Nring2share]+=1

    BeginRingTab = numpy.zeros([size],dtype='int32')
    EndRingTab = numpy.zeros([size],dtype='int32')
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
    add_data={}
    calib={}
    hdata={}
    hdata2={}
    hidx={}
    ridx={}
    val_ring={}
    nring={}
    icalib=-numpy.ones([nbolo],dtype='int32')
    
    if (rank==0):
        print('Nbolo ',nbolo)

    normrelative=False
    docalib=False

    #=============================================================================================================================
    #                                    MANAGE MASK
    #=============================================================================================================================
    mask=getmap(Mask,Nside,dtype='int32')

    for ib in range(nbolo):

        #=============================================================================================================================
        #                                    READ DATA
        #=============================================================================================================================
        data[ib]  = (getdata(Signal[ib],BeginRingTab[rank],EndRingTab[rank],step=RSTEP)-Monop[ib])/Calibration[ib]
        hdata[ib] = getdata(Hit[ib],BeginRingTab[rank],EndRingTab[rank],step=RSTEP)*(Calibration[ib]/NEP[ib])**2
        ph        = getdata(Ptg_PH[ib],BeginRingTab[rank],EndRingTab[rank],dtype='float64',step=RSTEP)
        th        = getdata(Ptg_TH[ib],BeginRingTab[rank],EndRingTab[rank],dtype='float64',step=RSTEP)
        psi       = getdata(Ptg_PSI[ib],BeginRingTab[rank],EndRingTab[rank],dtype='float64',step=RSTEP)
        hidx[ib]  =hp.ang2pix(Nside,th,ph)
        hdata2[ib]=hdata[ib]*mask[hidx[ib]]

        f_bad=numpy.fromfile(Badring[ib],dtype='int32')
        nx,ny=data[ib].shape
        ridx[ib]=(BeginRingTab[rank]+numpy.arange(nx)*RSTEP)
        idx=numpy.where(f_bad[ridx[ib]]==0)

        #=============================================================================================================================
        #                                    MANAGE LINEAR TEMPLATE VALUE
        #=============================================================================================================================
        ntmp=0
        ninput_template=len(Template[ib])
        for itmp in range(ninput_template):
            if Template[ib][itmp][3]!='RELATIV' or nbolo>1:
                if Template[ib][itmp][1]=='CALIB':
                    icalib[ib]=ntmp
                    docalib=True
                if Template[ib][itmp][1]=='LINEAR' or Template[ib][itmp][0]=='CALIB':
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
        
        frelat=numpy.zeros(ntmp,dtype='int')
        tmp_data[ib]=numpy.zeros([ntmp,nx,ny],dtype='float32')
        name_tmp={}
        ntmp=0
        for itmp in range(ninput_template):
            if Template[ib][itmp][3]!='RELATIV' or nbolo>1:
                if Template[ib][itmp][1]=='LINEAR' or Template[ib][itmp][1]=='CALIB':
                    if (rank==0):
                        print(Template[ib][itmp][0],Template[ib][itmp][1],' TEMPLATE ',Template[ib][itmp][2])
                    tmp_data[ib][ntmp,:,:]=getdata(Template[ib][itmp][2],BeginRingTab[rank],EndRingTab[rank],step=RSTEP)
                    if Template[ib][itmp][1]=='CALIB':
                        calib[ib]=1.0*(tmp_data[ib][ntmp,:,:]).reshape(nx,ny)
                    if Template[ib][itmp][3]=='RELATIV':
                        frelat[ntmp]=1
                    name_tmp[ntmp]=Template[ib][itmp][0]
                    ntmp=ntmp+1

        for itmp in range(ninput_template_map):
            if TemplateMap[ib][itmp][4]!='RELATIV' or nbolo>1:
                if TemplateMap[ib][itmp][1]=='LINEAR' or TemplateMap[ib][itmp][1]=='CALIB' or TemplateMap[ib][itmp][1]=='POLEFF' or TemplateMap[ib][itmp][1]=='POLANG':
                    if (rank==0):
                        print(TemplateMap[ib][itmp][0],TemplateMap[ib][itmp][1],' TEMPLATE MAP[',TemplateMap[ib][itmp][2],'] ',TemplateMap[ib][itmp][3])

                    lidx=hp.ang2pix(TemplateMap[ib][itmp][2],th,ph)

                    if TemplateMap[ib][itmp][1]=='POLEFF' or TemplateMap[ib][itmp][1]=='POLANG':
                        namemap=TemplateMap[ib][itmp][3].replace("@POLAR@", "Q")
                        imq=getmap(namemap,TemplateMap[ib][itmp][2])
                        namemap=TemplateMap[ib][itmp][3].replace("@POLAR@", "U")
                        imu=getmap(namemap,TemplateMap[ib][itmp][2])
                        if TemplateMap[ib][itmp][1]=='POLEFF':
                            tmp_data[ib][ntmp,:,:]= imq[lidx]*numpy.cos(2*psi)+imu[lidx]*numpy.sin(2*psi)
                        if TemplateMap[ib][itmp][1]=='POLANG':
                            tmp_data[ib][ntmp,:,:]= -imq[lidx]*numpy.sin(2*psi)+imu[lidx]*numpy.cos(2*psi)
                    else:
                        tmp_data[ib][ntmp,:,:] = getmap(TemplateMap[ib][itmp][3],TemplateMap[ib][itmp][2])[lidx]

                    if TemplateMap[ib][itmp][1]=='CALIB':
                        calib[ib]=1.0*(tmp_data[ib][ntmp,:,:]).reshape(nx,ny)
                    if TemplateMap[ib][itmp][4]=='RELATIV':
                        frelat[ntmp]=1
                    name_tmp[ntmp]=TemplateMap[ib][itmp][0]
                    ntmp=ntmp+1
        th=0
        ph=0
        psi=0

        #=============================================================================================================================
        #                                    MANAGE ADDED VALUE
        #=============================================================================================================================
        nadd=0
        ninput_add=len(Added[ib])
        for iadd in range(ninput_add):
            nadd=nadd+1
        add_data[ib]=numpy.zeros([nadd,nx,ny],dtype='float32')
        nadd=0
        for iadd in range(ninput_add):
            if (rank==0):
                print(Added[ib][iadd][0],' TEMPLATE ',Added[ib][iadd][1])
            add_data[ib][nadd,:,:]=getdata(Added[ib][iadd][1],BeginRingTab[rank],EndRingTab[rank],step=RSTEP)
            nadd=nadd+1

        

        data[ib]     = data[ib][idx[0],:].astype('float32')
        hdata[ib]    = hdata[ib][idx[0],:].astype('float32')
        hdata2[ib]    = hdata2[ib][idx[0],:].astype('float32')
        hidx[ib]     = hidx[ib][idx[0],:]
        tmp_data[ib] = tmp_data[ib][:,idx[0],:].astype('float32')
        add_data[ib] = add_data[ib][:,idx[0],:].astype('float32')
        if docalib:
            calib[ib]    = calib[ib][idx[0],:].astype('float32')

        nx,ny=data[ib].shape
        ridx[ib]=(numpy.repeat(ridx[ib][idx[0]],ny).reshape(nx,ny))

        data[ib]     = data[ib][hdata[ib]>0]
        hidx[ib]     = hidx[ib][hdata[ib]>0]
        if ntmp>0:
            tmp_data[ib] = tmp_data[ib][:,hdata[ib]>0]
        if nadd>0:
            add_data[ib] = add_data[ib][:,hdata[ib]>0]
        ridx[ib]     = ridx[ib][hdata[ib]>0]
        if docalib:
            calib[ib]    = calib[ib][hdata[ib]>0]
        hdata2[ib]    = hdata2[ib][hdata[ib]>0]
        hdata[ib]    = hdata[ib][hdata[ib]>0]
    
    
        # Find valid ring and replace it in ridx table
        fring= numpy.zeros([EndRing+1])
        fring[ridx[ib]]=1
        tfring = numpy.zeros([EndRing+1]) 
        if size>0:
            comm.Allreduce(fring,tfring)
        else:
            tfring=fring

        val_ring[ib]=numpy.where(tfring>0)[0]
        nring[ib]=len(val_ring[ib])
        
        if rank==0:
            print("NRING[%d] (valid) %d"%(ib,nring[ib]),icalib,nring[ib])
        tfring[val_ring[ib]]=numpy.arange(nring[ib])
        ridx[ib]=tfring[ridx[ib]].astype('int')

        
    #=============================================================================================================================
    #                                    MPI EXCHANGE DATA
    #=============================================================================================================================
    if size>0:
        if rank==0:
            rpidx=(numpy.random.rand(Nside*Nside*12)*size).astype('int32')
        else:
            rpidx=numpy.zeros([Nside*Nside*12],dtype='int32')

        comm.Bcast(rpidx, root=0) 
     
        for ib in range(nbolo):
            if rank==0:
                print('Transfer Data Bolo ',rank,ib)

            lpidx=rpidx[hidx[ib]]
            
            scounts=numpy.zeros([size],dtype='int')
            for i in range(size):
                scounts[i]=(lpidx==i).sum()

            tidx=numpy.argsort(lpidx)

            rcounts = numpy.array(comm.alltoall(scounts))

            if rank==0:
                print(rcounts,scounts)

            nlocal=int(rcounts.sum())

            buffer1=numpy.zeros([nlocal],dtype='float32')
            comm.Alltoallv(sendbuf=(data[ib][tidx].astype('float32'),scounts,MPI.FLOAT),recvbuf=(buffer1,rcounts,MPI.FLOAT)) 
            data[ib]=buffer1

            buffer1=numpy.zeros([nlocal],dtype='float32')
            comm.Alltoallv(sendbuf=(hdata2[ib][tidx].astype('float32'),scounts,MPI.FLOAT),recvbuf=(buffer1,rcounts,MPI.FLOAT)) 
            hdata2[ib]=buffer1

            if docalib:
                buffer1=numpy.zeros([nlocal],dtype='float32')
                comm.Alltoallv(sendbuf=(calib[ib][tidx].astype('float32'),scounts,MPI.FLOAT),recvbuf=(buffer1,rcounts,MPI.FLOAT)) 
                calib[ib]=buffer1

            buffer2=numpy.zeros([nlocal],dtype='float32')
            comm.Alltoallv(sendbuf=(hdata[ib][tidx].astype('float32'),scounts,MPI.FLOAT),recvbuf=(buffer2,rcounts,MPI.FLOAT)) 
            hdata[ib]=buffer2
        
            buffer3=numpy.zeros([nlocal],dtype='int32')
            comm.Alltoallv(sendbuf=(hidx[ib][tidx].astype('int32'),scounts,MPI.INT),recvbuf=(buffer3,rcounts,MPI.INT)) 
            hidx[ib]=buffer3
        
            buffer4=numpy.zeros([nlocal],dtype='int32')
            comm.Alltoallv(sendbuf=(ridx[ib][tidx].astype('int32'),scounts,MPI.INT),recvbuf=(buffer4,rcounts,MPI.INT)) 
            ridx[ib]=buffer4

            buffer5=numpy.zeros([nlocal],dtype='float32')
            ltmp_data=numpy.zeros([ntmp,nlocal],dtype='float32')
            for i in range(ntmp):
                comm.Alltoallv(sendbuf=((tmp_data[ib][i,tidx]).reshape(len(tidx)).astype('float32'),scounts,MPI.FLOAT),recvbuf=(buffer5,rcounts,MPI.FLOAT)) 
                ltmp_data[i,:]=buffer5
            tmp_data[ib]=ltmp_data

            ladd_data=numpy.zeros([nadd,nlocal],dtype='float32')
            for i in range(nadd):
                comm.Alltoallv(sendbuf=((add_data[ib][i,tidx]).reshape(len(tidx)).astype('float32'),scounts,MPI.FLOAT),recvbuf=(buffer5,rcounts,MPI.FLOAT)) 
                ladd_data[i,:]=buffer5
            add_data[ib]=ladd_data
            buffer5=0

    # Compress healpix index to the one used only
    fpix= numpy.zeros([12*Nside*Nside],dtype='int')
    
    for ib in range(nbolo):
        fpix[hidx[ib]]=1

    val_pidx=numpy.where(fpix>0)[0]
    npixel=len(val_pidx)
    fpix[val_pidx]=numpy.arange(npixel)
    for ib in range(nbolo):
        hidx[ib]=fpix[hidx[ib]]
    
    ndata=numpy.zeros([nbolo],dtype='int')
    for ib in range(nbolo):
        ndata[ib]=data[ib].shape[0]
    print("Rank %d : NPIXEL = %d, valid samples = "%(rank,npixel),ndata)

    if domodel:
        for ib in range(nbolo):
            data[ib]=numpy.cos((2*numpy.pi*(ridx[ib]+100*ib))/nring[ib])
            if ntmp>0:
                vv=numpy.array([0.003,0.002,0.0015,0.001,0.98]).reshape(1,ntmp)
                data[ib]+=numpy.dot(vv,tmp_data[ib]).reshape(ndata[ib])
    
            model=numpy.zeros([nring[ib]+ntmp])
            model[0:nring[ib]]=numpy.cos((2*numpy.pi*(numpy.arange(nring[ib])+100*ib))/nring[ib])
            model[nring[ib]:]=vv
            numpy.save(Out_Offset[ib]+'_mod',model)
        

    if rank==0:
        for ib in range(nbolo):
            print('DATA INFO',rank,data[ib].min(),data[ib].max(),hdata[ib].min(),hdata[ib].max(),ridx[ib].min(),ridx[ib].max(),hidx[ib].min(),hidx[ib].max())

    comm.Barrier()

    #=============================================================================================================================
    #                                    DESIGN TENSORFLOW MINIMIZATION
    #=============================================================================================================================
    
    train_data     = {}
    train_tmp_data = {}
    train_hdata    = {}
    train_hdata2    = {}
    train_pidx     = {}
    train_ridx     = {}

    for ib in range(nbolo):
        train_tmp_data[ib] = tf.constant(tmp_data[ib].reshape(ntmp,ndata[ib]).astype('float64'))
        train_hdata[ib]    = tf.constant(hdata[ib].reshape(ndata[ib]).astype('float64'))
        train_hdata2[ib]   = tf.constant(hdata2[ib].reshape(ndata[ib]).astype('float64'))
        train_pidx[ib]     = tf.constant(hidx[ib].reshape(ndata[ib]).astype('int32'))
        train_ridx[ib]     = tf.constant(ridx[ib].reshape(ndata[ib]).astype('int32'))

    nfit_bolo=[(ntmp+nring[i]) for i in nring]
    nfit=numpy.array(nfit_bolo).sum()

    #=============================================================================================================================
    #                                   COMPUTE INDEX NORM
    #=============================================================================================================================
    if normrelative:
        nnorm=frelat.sum()
        relat_idx=numpy.zeros([nfit],dtype='int32')+ntmp
        nidx=0
        for ib in range(nbolo):
            nidx+=nring[ib]
            for itmp in range(ntmp):
                if frelat[i]==1:
                    relat_idx[nidx+itmp]=itmp
                else:
                    relat_idx[nidx+itmp]=ntmp
            nidx+=ntmp

    train_xk = tf.compat.v1.placeholder(tf.float64,shape=(nfit))
    if normrelative:
        lnorm_idx = tf.constant(relat_idx.astype('int32'))

    ntotdata=0
    for ib in range(nbolo):
        ntotdata+=ndata[ib]
    train_data = tf.compat.v1.placeholder(tf.float64,shape=(ntotdata))
        
    def domap(xval,the_data,tbolo=[0]):
        ldata=tf.split(the_data,ndata)
        xval_bolo=tf.split(xval,nfit_bolo)
        for ib in tbolo:
            off,vtmp=tf.split(xval_bolo[ib],[nring[ib],ntmp])
            vdestrip=tf.gather(off,train_ridx[ib])
            if ntmp>0:
                vdestrip = vdestrip + tf.reshape(tf.matmul(tf.reshape(vtmp,[1,ntmp]),train_tmp_data[ib]),[ndata[ib]])
            if ib==tbolo[0]:
                omap=tf.math.unsorted_segment_sum(train_hdata[ib]*(ldata[ib]-vdestrip),train_pidx[ib],npixel)
                homap=tf.math.unsorted_segment_sum(train_hdata[ib],train_pidx[ib],npixel)
            else:
                omap=omap+tf.math.unsorted_segment_sum(train_hdata[ib]*(ldata[ib]-vdestrip),train_pidx[ib],npixel)
                homap=homap+tf.math.unsorted_segment_sum(train_hdata[ib],train_pidx[ib],npixel)
        return(omap,homap)
    
    def proj(xval):
        xval_bolo=tf.split(xval,nfit_bolo)
        vdestrip={}
        for ib in range(nbolo):
            off,vtmp=tf.split(xval_bolo[ib],[nring[ib],ntmp])
            vdestrip[ib]=tf.gather(off,ridx[ib])
            if ntmp>0:
                vdestrip[ib] = vdestrip[ib] + tf.reshape(tf.matmul(tf.reshape(vtmp,[1,ntmp]),train_tmp_data[ib]),[ndata[ib]])
            if ib==0:
                omap=tf.math.unsorted_segment_sum(train_hdata2[ib]*vdestrip[ib],train_pidx[ib],npixel)
                ohmap=tf.math.unsorted_segment_sum(train_hdata2[ib],train_pidx[ib],npixel)
                lnorm=tf.reduce_sum(off)
            else:
                omap=omap+tf.math.unsorted_segment_sum(train_hdata2[ib]*vdestrip[ib],train_pidx[ib],npixel)
                ohmap=ohmap+tf.math.unsorted_segment_sum(train_hdata2[ib],train_pidx[ib],npixel)
                lnorm=lnorm+tf.reduce_sum(off)
                
        #UN PEU VIOLENT SI PAS DE FIT DE MONOPOLE
        #lnorm       = tf.reduce_sum(omap)
        if normrelative:
            lnorm_relat,garbage = tf.split(tf.math.unsorted_segment_sum(xval,lnorm_idx,ntmp+1),[ntmp,1])

        lohmap={}
        for ib in range(nbolo):
            lomap=tf.gather(omap,train_pidx[ib])
            lohmap[ib]=tf.gather(ohmap,train_pidx[ib])
            vdestrip[ib]= (lohmap[ib]*vdestrip[ib]-lomap)
            if ib==0:
                vmap=tf.math.unsorted_segment_sum(vdestrip[ib],train_pidx[ib],npixel)
            else:
                vmap=vmap+tf.math.unsorted_segment_sum(vdestrip[ib],train_pidx[ib],npixel)

        for ib in range(nbolo):
            lvmap=tf.gather(vmap,train_pidx[ib])
            dvdestrip=(lohmap[ib]*vdestrip[ib]-train_hdata2[ib]*lvmap)
            part1=tf.reshape(tf.nn.bias_add(tf.reshape(tf.math.unsorted_segment_sum(dvdestrip,ridx[ib],nring[ib]),[nring[ib],1]),
                                            tf.reshape(lnorm,[1])),[nring[ib]])
            if ntmp>0:
                part2=tf.reshape(tf.matmul(train_tmp_data[ib],tf.reshape(dvdestrip,[ndata[ib],1])),[ntmp])
                if normrelative:
                    part2=part2+lnorm_relat

            if ib==0:
                if ntmp>0:
                    res=tf.concat([part1,part2],0)
                else:
                    res=part1
            else:
                if ntmp>0:
                    res=tf.concat([res,part1,part2],0)
                else:
                    res=tf.concat([res,part1],0)
        return(res)
  
    def projdata(the_data):
        ldata=tf.split(the_data,ndata)
        for ib in range(nbolo):
            if ib==0:
                omap=tf.math.unsorted_segment_sum(train_hdata2[ib]*ldata[ib],train_pidx[ib],npixel)
                ohmap=tf.math.unsorted_segment_sum(train_hdata2[ib],train_pidx[ib],npixel)
            else:
                omap=omap+tf.math.unsorted_segment_sum(train_hdata2[ib]*ldata[ib],train_pidx[ib],npixel)
                ohmap=ohmap+tf.math.unsorted_segment_sum(train_hdata2[ib],train_pidx[ib],npixel)

        destrip={}
        lohmap={}
        for ib in range(nbolo):
            lomap=tf.gather(omap,train_pidx[ib])
            lohmap[ib]=tf.gather(ohmap,train_pidx[ib])
            destrip[ib]=(lohmap[ib]*ldata[ib]-lomap)
            if ib==0:
                vmap=tf.math.unsorted_segment_sum(destrip[ib],train_pidx[ib],npixel)
            else:
                vmap=vmap+tf.math.unsorted_segment_sum(destrip[ib],train_pidx[ib],npixel)

        for ib in range(nbolo):
            lvmap=tf.gather(vmap,train_pidx[ib])
            ddestrip=(lohmap[ib]*destrip[ib]-train_hdata2[ib]*lvmap)
            
            part1=tf.math.unsorted_segment_sum(ddestrip,train_ridx[ib],nring[ib])
            if ntmp>0:
                part2=tf.reshape(tf.matmul(train_tmp_data[ib],tf.reshape(ddestrip,[ndata[ib],1])),[ntmp])

            if ib==0:
                if ntmp>0:
                    res=tf.concat([part1,part2],0)
                else:
                    res=part1
            else:
                if ntmp>0:
                    res=tf.concat([res,part1,part2],0)
                else:
                    res=tf.concat([res,part1],0)
        return(res)

    themap={}
    thehmap={}
    for imap in range(len(MAP)):
        themap[imap],thehmap[imap]=domap(train_xk,train_data,tbolo=MAP[imap][2])
  
    b=projdata(train_data)
    y=proj(train_xk)

    gain=numpy.ones([nbolo])

    # Create a local session to run the training.
    start_time = time.time()
    with tf.compat.v1.Session() as sess:
        # Run all the initializers to prepare the trainable parameters.
        tf.compat.v1.global_variables_initializer().run()
        if rank==0:
            print('Initialized!')
        # Loop through training steps.

        xk=numpy.zeros([nfit],dtype='float')
        if size>0:
            rk=numpy.zeros([nfit],dtype='float')
            Apk=numpy.zeros([nfit],dtype='float')

        for itt in range(NITT):

            if rank==0:
                print('=========================================')
                print('==     ITT = %4d                      =='%(itt))
                print('=========================================')
                print('GAIN',gain)

            totdata=numpy.array([])
            for ib in range(nbolo):
                if docalib:
                    tmp=gain[ib]*numpy.array(data[ib])-numpy.array(calib[ib])
                else:
                    tmp=numpy.array(data[ib])

                for iadd in range(nadd):
                    tmp+=Added[ib][iadd][0]*add_data[ib][iadd,:]
                    
                totdata=numpy.append(totdata,[tmp])

            lb,ly = sess.run([b,y],feed_dict={train_xk:xk,
                                              train_data:totdata.flatten()})
            lrk=lb-ly
        
            if size>0:
                comm.Allreduce(lrk,rk)
            else:
                rk=lrk

            pk=1.0*rk
            deltak=(rk**2).sum()
            if itt==0:
                delta0=1.0*deltak
            delta=1.0*deltak

            step=0
            while(step<NUM_EPOCHS and delta>1E-30*delta0):
            
                delta=1.0*deltak

                lApk=sess.run([y],feed_dict={train_xk:pk})[0]
                if size>0:
                    comm.Allreduce(lApk,Apk)
                else:
                    Apk=lApk
            
                alpha=delta/(pk*Apk).sum()
            
                xk = xk + alpha*pk
                rk = rk - alpha*Apk

                deltak=(rk**2).sum()
                
                beta = deltak/delta

                pk=rk + beta*pk

                if step%EVAL_FREQUENCY==0:
                    elapsed_time = time.time() - start_time
                    start_time = time.time()
                    if rank==0:
                        print('Step %d, %.3f s Delta0: %.6g, delta: %.6g' % (step,elapsed_time/EVAL_FREQUENCY,delta0,delta))

                step+=1
        

            nn=0
            vtmp=numpy.zeros([ntmp,nbolo])
            for ib in range(nbolo):
                if icalib[ib]!=-1:
                    gain[ib]-=xk[nn+nring[ib]+icalib[ib]]

                vtmp[:,ib]=xk[nn+nring[ib]:nn+nring[ib]+ntmp]
                nn+=nring[ib]+ntmp

            if rank==0:
                for ib in range(ntmp):
                    print('TMP %s = '%(name_tmp[ib]),[vtmp[ib,i] for i in range(nbolo)],' MEAN = ',vtmp[ib,:].mean(),' STD = ',vtmp[ib,:].std())

    #=============================================================================================================================
    #                                  SAVE DATA
    #=============================================================================================================================
            if itt==NITT-1:
                if rank==0:
                    nn=0
                    for ib in range(nbolo):
                        theoff=numpy.zeros([EndRing+1],dtype='float32')
                        theoff[:]=hp.UNSEEN
                        theoff[val_ring[ib]]=xk[nn:nn+nring[ib]]
                        numpy.save(Out_Offset[ib],theoff)

                        nn+=nring[ib]+ntmp

                for imap in range(len(MAP)):
                    map_fin,hmap_fin=sess.run([themap[imap],thehmap[imap]],
                                              feed_dict={train_xk:xk,
                                              train_data:totdata.flatten()})
                    lthemap=numpy.zeros([12*Nside*Nside],dtype='float32')
                    lthehmap=numpy.zeros([12*Nside*Nside],dtype='float32')
                    if size==0:
                        lthemap[:]=hp.UNSEEN

                    lthemap[val_pidx]=map_fin
                    lthehmap[val_pidx]=hmap_fin

                    if size>0:
                        allmap=numpy.zeros([12*Nside*Nside],dtype='float32')
                        allhmap=numpy.zeros([12*Nside*Nside],dtype='float32')
                        comm.Allreduce((lthemap.astype('float32'),MPI.FLOAT),(allmap,MPI.FLOAT))
                        comm.Allreduce((lthehmap.astype('float32'),MPI.FLOAT),(allhmap,MPI.FLOAT))
                    else:
                        allmap=lthemap
                        allhmap=lthehmap

                    if rank==0:
                        allmap[allhmap>0]/=allhmap[allhmap>0]
                        allmap[allhmap==0]=hp.UNSEEN

                        hp.write_map(MAP[imap][0]+'_I',allmap,overwrite=True)
                        hp.write_map(MAP[imap][0]+'_H',allhmap,overwrite=True)
                        print(MAP[imap][0],' Done')

                comm.Barrier()

        if rank==0:
            print('TOTAL DURATION ',time.time() - start_time_total)

#=============================================================================================================================
#                                  MAIN : flag to be improved and used
#=============================================================================================================================
if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument(
      '--use_fp16',
      default=False,
      help='Use half floats instead of full floats if True.',
      action='store_true')
  parser.add_argument(
      '--self_test',
      default=False,
      action='store_true',
      help='True if running a self test.')

  FLAGS, unparsed = parser.parse_known_args()
  tf.compat.v1.app.run(main=main, argv=[sys.argv[0]] + unparsed)

