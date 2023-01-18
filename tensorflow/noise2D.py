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

DATAPATH='/scratch/cnt0028/ias1717/SHARED/RD12_data/dmc/MISS03/DATA/PBR_JMD/'
OUTNAME='545_FSL'

begr=240
edr=26050
RSTEP=10
nshift=5

bolo=['545-1','545-4']

refbolo=[sys.argv[1]]
doplot=int(sys.argv[2])==1
OUTNAME='857_FSL'
if sys.argv[1]=='857-1':
    bolo=['857-2','857-3','857-4']
    lamp=27.5
if sys.argv[1]=='857-2':
    bolo=['857-1','857-3','857-4']
    lamp=2.0
if sys.argv[1]=='857-3':
    bolo=['857-1','857-2','857-4']
    lamp=18.0
if sys.argv[1]=='857-4':
    bolo=['857-1','857-2','857-3'] 
    lamp=8.0


nx=512
ny=1024
vdata=np.zeros([len(bolo)+1,ny,nx],dtype='float32')
hdata=np.zeros([len(bolo)+1,ny,nx],dtype='float32')

dosquare=False
dolinear=False

omasx=hp.smoothing(hp.ud_grade(np.fromfile('/scratch/cnt0028/ias1717/SHARED/bware/sroll22_trick/SrollEx/tensorflow//mask_RD_545_10DEG',dtype='float32'),64),3/180.*np.pi)
imref=hp.ud_grade(hp.read_map('/scratch/cnt0028/ias1717/SHARED/bware/sroll2r94_jmd/545/MAP/%s_RD12_REP6_DATA_NORM_FSL_4_-2_Full_I.fits'%(refbolo[0])),64)
    
# ==================================================================================================================
#
#   HERE IS THE PLACE WHERE YOU CAN CHANGE THE HPR TO FIT FOR THE MODEL
#
# ==================================================================================================================

fsl=np.fromfile(DATAPATH+'/%s_REP7_2'%(refbolo[0]),dtype='float32',count=(edr+1)*27664)
hh=np.fromfile(DATAPATH+'/%s_REP6_hit'%(refbolo[0]),dtype='float32',count=(edr+1)*27664)

iy2=np.load(DATAPATH+'/CODEIDX_%s_%s.npy'%(OUTNAME,refbolo[0]))
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
    
pa=np.fromfile(DATAPATH+'/%s_REP6_phregul'%(refbolo[0]),dtype='float32',count=(edr+1)*27664)
ph=np.fromfile(DATAPATH+'/%s_REP6_ptg'%(refbolo[0]),dtype='float',count=(edr+1)*27664)
th=np.fromfile(DATAPATH+'/%s_REP6_ptg_TUPLE_1'%(refbolo[0]),dtype='float',count=(edr+1)*27664)
hh=hh[begr*27664:].reshape(edr-begr+1,27664)[RSTEP*np.arange((edr-begr+1)//RSTEP),:].flatten()
th=th[begr*27664:].reshape(edr-begr+1,27664)[RSTEP*np.arange((edr-begr+1)//RSTEP),:].flatten()
ph=ph[begr*27664:].reshape(edr-begr+1,27664)[RSTEP*np.arange((edr-begr+1)//RSTEP),:].flatten()
fsl=fsl[begr*27664:].reshape(edr-begr+1,27664)[RSTEP*np.arange((edr-begr+1)//RSTEP),:].flatten()
pa=pa[begr*27664:].reshape(edr-begr+1,27664)[RSTEP*np.arange((edr-begr+1)//RSTEP),:].flatten()
pidx=hp.ang2pix(64,th,ph)
    
nring=len(pa)//27664
idx=np.where(hh>0)[0]

x1=(nx*(pa[idx]/np.pi/2)).astype('int')+nx*(iy2[idx//27664]).astype('int')
    
nfres = np.bincount(x1,weights=hh[idx],minlength=ny*nx)
res   = np.bincount(x1,weights=hh[idx]*(fsl[idx]),minlength=ny*nx)
res[nfres>0]/=nfres[nfres>0]
vdata[0,:,:]=res.reshape(ny,nx)
hdata[0,:,:]=nfres.reshape(ny,nx)
    
co=np.fromfile('/scratch/cnt0028/ias1717/SHARED/RD12_data/dmc/MISS03/DATA/sroll_templates/MAP_JMD_128/NEW_cleaned_12CO_K.float32.bin',dtype='float32')
co=hp.ud_grade(co,64)

for ib in range(len(bolo)):
    im=hp.ud_grade(hp.read_map('/scratch/cnt0028/ias1717/SHARED/bware/sroll2r94_jmd/545/MAP/%s_RD12_REP6_DATA_NORM_FSL_4_-2_Full_I.fits'%(bolo[ib])),64)
    mat=np.array([[(omasx*im*im).sum(),(omasx*co*im).sum(),(omasx*im).sum()],
                  [(omasx*co*im).sum(),(omasx*co*co).sum(),(omasx*co).sum()],
                  [(omasx*im).sum(),   (omasx*co).sum(),   (omasx).sum()]])
    vec=np.array([(omasx*imref*im).sum(),(omasx*co*imref).sum(),(omasx*imref).sum()])
    rr=np.linalg.solve(mat,vec)
    print(bolo[ib])
    print(rr)
    refim=rr[0]*im+rr[1]*co+rr[2]-imref
    vx,vy,vz=hp.pix2vec(64,np.arange(12*64*64))
    v0,v1,v2=hp.ang2vec((90+20.42)/180.*np.pi,(-150.62)/180.*np.pi)
    w0,w1,w2=hp.ang2vec((90+17.2)/180.*np.pi,(-0.6)/180.*np.pi)
    o0,o1,o2=hp.ang2vec((90+52.53)/180.*np.pi,(107.52)/180.*np.pi)
    p0,p1,p2=hp.ang2vec((90+47.21)/180.*np.pi,(36.4)/180.*np.pi)
    q0,q1,q2=hp.ang2vec((90-2.94)/180.*np.pi,(-157.1)/180.*np.pi)
    hp.mollview(omasx*refim,cmap='jet',hold=False,sub=(2,1,1))
    amp=64
    amp=128
    refim*=1-np.exp(-amp*((vx-v0)**2+(vy-v1)**2+(vz-v2)**2))-np.exp(-amp*((vx-w0)**2+(vy-w1)**2+(vz-w2)**2))-np.exp(-amp*((vx-o0)**2+(vy-o1)**2+(vz-o2)**2))-np.exp(-amp*((vx-p0)**2+(vy-p1)**2+(vz-p2)**2))-np.exp(-amp*((vx-q0)**2+(vy-q1)**2+(vz-q2)**2))
    refim=omasx*hp.smoothing(refim,1/180.*np.pi)
    if doplot==True:
        hp.mollview(refim,cmap='jet',hold=False,sub=(2,1,2))
        plt.show()
    
    res2  = np.bincount(x1,weights=hh[idx]*(refim[pidx[idx]]),minlength=ny*nx)

    print(fsl[idx].std(),res.max(),res2.max())
    
    res2[nfres>0]/=nfres[nfres>0]
    
    vdata[ib+1,:,:]=0*res.reshape(ny,nx)+res2.reshape(ny,nx) #/lamp
    hdata[ib+1,:,:]=nfres.reshape(ny,nx)
    if doplot==True:
        plt.subplot(1,3,1)
        plt.imshow(res.reshape(ny,nx),cmap='jet')
        plt.subplot(1,3,2)
        plt.imshow(res2.reshape(ny,nx),cmap='jet')
        plt.subplot(1,3,3)
        plt.imshow(vdata[ib+1,:,:],cmap='jet')
        plt.show()
    
##vdata[ib*5+1,:,nshift:nx]=vdata[ib*5,:,0:nx-nshift]
##vdata[ib*5+1,:,0:nshift]=vdata[ib*5,:,nx-nshift:]
##hdata[ib*5+1,:,nshift:nx]=hdata[ib*5,:,0:nx-nshift]
##hdata[ib*5+1,:,0:nshift]=hdata[ib*5,:,nx-nshift:nx]
##vdata[ib*5+2,:,0:nx-nshift]=vdata[ib*5,:,nshift:nx]
##vdata[ib*5+2,:,nx-nshift:nx]=vdata[ib*5,:,0:nshift]
##hdata[ib*5+2,:,0:nx-nshift]=hdata[ib*5,:,nshift:nx]
##hdata[ib*5+2,:,nx-nshift:nx]=hdata[ib*5,:,0:nshift]
#
#fsl2=(vdata[ib,:,:].reshape(ny*nx))[x1]
#im1=np.bincount(pidx[idx],weights=hh[idx]*fsl[idx],minlength=12*64*64)
#im2=np.bincount(pidx[idx],weights=hh[idx]*(fsl[idx]-fsl2),minlength=12*64*64)
#imh=np.bincount(pidx[idx],weights=hh[idx],minlength=12*64*64)
#im1/=imh
#im2/=imh
#hp.mollview(im1-np.median(im1),cmap='jet',hold=False,sub=(2,len(bolo),1+ib),min=-1E-3,max=1e-3)
#hp.mollview(im2,cmap='jet',hold=False,sub=(2,len(bolo),1+len(bolo)+ib),min=-1E-4,max=1e-4)
#
##x1=(nx*(pa[idx]/np.pi/2)).astype('int')+nx*(iy3[idx//27664]).astype('int')
##res=np.bincount(x1,weights=hh[idx]*(fsl[idx]),minlength=ny*nx)
##nfres=np.bincount(x1,weights=hh[idx],minlength=ny*nx)
##res[nfres>0]/=nfres[nfres>0]
##nfres=nfres.reshape(ny,nx)
##vdata[ib*5+3,:,:]=res.reshape(ny,nx)
##hdata[ib*5+3,nfres>0]=nfres[nfres>0]
##
##x1=(nx*(pa[idx]/np.pi/2)).astype('int')+nx*(iy4[idx//27664]).astype('int')
##res=np.bincount(x1,weights=hh[idx]*(fsl[idx]),minlength=ny*nx)
##nfres=np.bincount(x1,weights=hh[idx],minlength=ny*nx)
##res[nfres>0]/=nfres[nfres>0]
##nfres=nfres.reshape(ny,nx)
##vdata[ib*5+4,:,:]=res.reshape(ny,nx)
##hdata[ib*5+4,nfres>0]=nfres[nfres>0]
#
##plt.figure(figsize=(8,8))
##for i in range(4):
##    plt.subplot(2,2,1+i)
##    plt.imshow((vdata[ib*5+1+i,:,:]-vdata[ib*5,:,:])*(hdata[ib*5+1+i,:,:]>0)*(hdata[ib*5,:,:]>0),
##               cmap='jet',vmin=-0.0003,vmax=0.0003)
##
#show()

np.save(DATAPATH+'/ERRCODE_%s_%s_ALL.npy'%(OUTNAME,refbolo[0]),vdata)
np.save(DATAPATH+'/ERRCODE_%s_%s_HALL.npy'%(OUTNAME,refbolo[0]),hdata)
