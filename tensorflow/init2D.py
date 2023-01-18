#!/venv/py3-phyocean/bin/python

from scipy.ndimage import gaussian_filter
import numpy as np
import sys
import time
import logging
import os
import matplotlib.pyplot as plt
import argparse
import healpy as hp

DATAPATH='/scratch/cnt0028/ias1717/SHARED/RD12_data/dmc/MISS03/DATA/PBR_JMD'
OUTNAME='545_FSL'
OUTNAME='857_FSL'

begr=240
edr=26050
RSTEP=10
nshift=5

bolo=['545-1','545-4']

bolo=['857-1','857-2','857-3','857-4','545-1','545-2','545-4']
nx=512
ny=1024
vdata=np.zeros([len(bolo),ny,nx],dtype='float32')
hdata=np.zeros([len(bolo),ny,nx],dtype='float32')

dosquare=False
dolinear=False

for ib in range(len(bolo)):
    print(bolo[ib])
    
# ==================================================================================================================
#
#   HERE IS THE PLACE WHERE YOU CAN CHANGE THE HPR TO FIT FOR THE MODEL
#
# ==================================================================================================================
 
    fsl=np.fromfile(DATAPATH+'/%s_REP7_2'%(bolo[ib]),dtype='float32',count=(edr+1)*27664)
    hh=np.fromfile(DATAPATH+'/%s_REP6_hit'%(bolo[ib]),dtype='float32',count=(edr+1)*27664)
    
# ==================================================================================================================
#
#   THE NEXT LINES ARE THERE TO COMPUT THE INDEX OF THE Y AXIS IN SUCH WAY THAT THE VARIATION OF THE SIGNAL
#
# ==================================================================================================================
    y2=np.bincount(10*nshift+np.arange(len(fsl))//27664,weights=hh*fsl*fsl,minlength=edr+20*nshift+1)
    ny2=np.bincount(10*nshift+np.arange(len(fsl))//27664,weights=hh,minlength=edr+20*nshift+1)
    y2[ny2>0]/=ny2[ny2>0]
    y2[edr+10*nshift+1:]=y2[edr+10*nshift]
    if dosquare==False:
        y2=np.sqrt(y2)

    y3=1.0*y2[0:edr+1]
    y4=1.0*y2[20*nshift:edr+20*nshift+1]
    y2=1.0*y2[10*nshift:edr+10*nshift+1]
    
    for i in range(len(y2)-1):
        y2[i+1]+=y2[i]
        y3[i+1]+=y3[i]
        y4[i+1]+=y4[i]
        
    y2*=ny/y2[edr]
    iy2=y2.astype('int')
    y3*=ny/y3[edr]
    iy3=y3.astype('int')
    y4*=ny/y4[edr]
    iy4=y4.astype('int')
    
    if dolinear==True:
        iy2=(ny*(np.arange(edr+1)-begr))//(edr+1-begr)
    iy2[iy2<0]=0
    iy2[iy2==ny]=ny-1
    iy3[iy3<0]=0
    iy3[iy3==ny]=ny-1
    iy4[iy4<0]=0
    iy4[iy4==ny]=ny-1

    if dosquare==False:
        if dolinear==True:
            np.save(DATAPATH+'/CODEIDX1_%s_%s.npy'%(OUTNAME,bolo[ib]),iy2)
        else:
            np.save(DATAPATH+'/CODEIDX_%s_%s.npy'%(OUTNAME,bolo[ib]),iy2)
    else:
        np.save(DATAPATH+'/CODEIDX2_%s_%s.npy'%(OUTNAME,bolo[ib]),iy2)
# ==================================================================================================================
#
#   END OF THE INDEX COMPUTATION AND SAVE IT FOR SROLL4
#
# ==================================================================================================================
# ==================================================================================================================
#
#   AND NOW COMPUTE THE PROJECTION WEIGTHED BY THE HIT COUNT
#
# ==================================================================================================================
    
    iy2=iy2[begr+RSTEP*np.arange((edr-begr+1)//RSTEP)]
    iy3=iy3[begr+RSTEP*np.arange((edr-begr+1)//RSTEP)]
    iy4=iy4[begr+RSTEP*np.arange((edr-begr+1)//RSTEP)]
    
    pa=np.fromfile(DATAPATH+'/%s_REP6_phregul'%(bolo[ib]),dtype='float32',count=(edr+1)*27664)
    ph=np.fromfile(DATAPATH+'/%s_REP6_ptg'%(bolo[ib]),dtype='float',count=(edr+1)*27664)
    th=np.fromfile(DATAPATH+'/%s_REP6_ptg_TUPLE_1'%(bolo[ib]),dtype='float',count=(edr+1)*27664)
    hh=hh[begr*27664:].reshape(edr-begr+1,27664)[RSTEP*np.arange((edr-begr+1)//RSTEP),:].flatten()
    th=th[begr*27664:].reshape(edr-begr+1,27664)[RSTEP*np.arange((edr-begr+1)//RSTEP),:].flatten()
    ph=ph[begr*27664:].reshape(edr-begr+1,27664)[RSTEP*np.arange((edr-begr+1)//RSTEP),:].flatten()
    fsl=fsl[begr*27664:].reshape(edr-begr+1,27664)[RSTEP*np.arange((edr-begr+1)//RSTEP),:].flatten()
    pa=pa[begr*27664:].reshape(edr-begr+1,27664)[RSTEP*np.arange((edr-begr+1)//RSTEP),:].flatten()
    pidx=hp.ang2pix(64,th,ph)
    
    nring=len(pa)//27664
    idx=np.where(hh>0)[0]
    
    x1=(nx*(pa[idx]/np.pi/2)).astype('int')+nx*(iy2[idx//27664]).astype('int')
    
    nfres=np.bincount(x1,weights=hh[idx],minlength=ny*nx)
    res=np.bincount(x1,weights=hh[idx]*(fsl[idx]),minlength=ny*nx)

    print(fsl[idx].std())
    
    res[nfres>0]/=nfres[nfres>0]
    nfres=nfres.reshape(ny,nx)
    vdata[ib,:,:]=res.reshape(ny,nx)
    hdata[ib,nfres>0]=nfres[nfres>0]

    #vdata[ib*5+1,:,nshift:nx]=vdata[ib*5,:,0:nx-nshift]
    #vdata[ib*5+1,:,0:nshift]=vdata[ib*5,:,nx-nshift:]
    #hdata[ib*5+1,:,nshift:nx]=hdata[ib*5,:,0:nx-nshift]
    #hdata[ib*5+1,:,0:nshift]=hdata[ib*5,:,nx-nshift:nx]
    #vdata[ib*5+2,:,0:nx-nshift]=vdata[ib*5,:,nshift:nx]
    #vdata[ib*5+2,:,nx-nshift:nx]=vdata[ib*5,:,0:nshift]
    #hdata[ib*5+2,:,0:nx-nshift]=hdata[ib*5,:,nshift:nx]
    #hdata[ib*5+2,:,nx-nshift:nx]=hdata[ib*5,:,0:nshift]
    
    fsl2=(vdata[ib,:,:].reshape(ny*nx))[x1]
    im1=np.bincount(pidx[idx],weights=hh[idx]*fsl[idx],minlength=12*64*64)
    im2=np.bincount(pidx[idx],weights=hh[idx]*(fsl[idx]-fsl2),minlength=12*64*64)
    imh=np.bincount(pidx[idx],weights=hh[idx],minlength=12*64*64)
    im1/=imh
    im2/=imh
    hp.mollview(im1-np.median(im1),cmap='jet',hold=False,sub=(2,len(bolo),1+ib),min=-1E-3,max=1e-3)
    hp.mollview(im2,cmap='jet',hold=False,sub=(2,len(bolo),1+len(bolo)+ib),min=-1E-4,max=1e-4)
    
    #x1=(nx*(pa[idx]/np.pi/2)).astype('int')+nx*(iy3[idx//27664]).astype('int')
    #res=np.bincount(x1,weights=hh[idx]*(fsl[idx]),minlength=ny*nx)
    #nfres=np.bincount(x1,weights=hh[idx],minlength=ny*nx)
    #res[nfres>0]/=nfres[nfres>0]
    #nfres=nfres.reshape(ny,nx)
    #vdata[ib*5+3,:,:]=res.reshape(ny,nx)
    #hdata[ib*5+3,nfres>0]=nfres[nfres>0]
    #
    #x1=(nx*(pa[idx]/np.pi/2)).astype('int')+nx*(iy4[idx//27664]).astype('int')
    #res=np.bincount(x1,weights=hh[idx]*(fsl[idx]),minlength=ny*nx)
    #nfres=np.bincount(x1,weights=hh[idx],minlength=ny*nx)
    #res[nfres>0]/=nfres[nfres>0]
    #nfres=nfres.reshape(ny,nx)
    #vdata[ib*5+4,:,:]=res.reshape(ny,nx)
    #hdata[ib*5+4,nfres>0]=nfres[nfres>0]

    #plt.figure(figsize=(8,8))
    #for i in range(4):
    #    plt.subplot(2,2,1+i)
    #    plt.imshow((vdata[ib*5+1+i,:,:]-vdata[ib*5,:,:])*(hdata[ib*5+1+i,:,:]>0)*(hdata[ib*5,:,:]>0),
    #               cmap='jet',vmin=-0.0003,vmax=0.0003)
    #
plt.show()

if dosquare==False:
    if dolinear==True:
        np.save(DATAPATH+'/INITCODE1_%s_ALL.npy'%(OUTNAME),vdata)
        np.save(DATAPATH+'/INITCODE1_%s_HALL.npy'%(OUTNAME),hdata)
    else:
        np.save(DATAPATH+'/INITCODE_%s_ALL.npy'%(OUTNAME),vdata)
        np.save(DATAPATH+'/INITCODE_%s_HALL.npy'%(OUTNAME),hdata)
else:
    np.save(DATAPATH+'/INITCODE2_%s_ALL.npy'%(OUTNAME),vdata)
    np.save(DATAPATH+'/INITCODE2_%s_HALL.npy'%(OUTNAME),hdata)
