import numpy as np
import matplotlib.pyplot as plt
import healpy as hp

nout2=2048
nout=256
bolo=['857-1','857-2','857-3','857-4']

ref=hp.ud_grade(np.fromfile('/home1/scratch/jmdeloui/DATA4SROLL4/map_857_2018.float32.bin',dtype='float32'),nout)
Calibration = [3.30076826046e-16,3.55811287601e-16,3.18681631353e-16,2.219187708e-16]
"""
im={}
for ib in range(len(bolo)):
    print(bolo[ib])
    map=np.zeros([12*nout2*nout2])
    hmap=np.zeros([12*nout2*nout2])
    for i in range(240,26005,100):
        ss=np.fromfile('/export/home1/jmdeloui/DATA4SROLL4/%s_REP6_SIGTEST'%(bolo[ib]),dtype='float32',offset=i*27664*4,count=27664)/Calibration[ib]
        hh=np.fromfile('/export/home1/jmdeloui/DATA4SROLL4/%s_REP6_hit'%(bolo[ib]),dtype='float32',offset=i*27664*4,count=27664)
        ph=np.fromfile('/export/home1/jmdeloui/DATA4SROLL4/%s_REP6_ptg'%(bolo[ib]),dtype='float',offset=i*27664*8,count=27664)
        th=np.fromfile('/export/home1/jmdeloui/DATA4SROLL4/%s_REP6_ptg_TUPLE_1'%(bolo[ib]),dtype='float',offset=i*27664*8,count=27664)
    
        pidx=hp.ang2pix(nout2,th,ph)
        map[pidx]+=hh*ss
        hmap[pidx]+=hh

    map[hmap>0]/=hmap[hmap>0]
    map[hmap==0]=hp.UNSEEN
    im[ib]=hp.ud_grade(map,nout)


bolo=['857-1','857-2','857-3','857-4']
for ib in range(len(bolo)):
    off=np.fromfile('/home1/scratch/jmdeloui/VEC/temp_Nside2048_%s_offset_OFF'%(bolo[ib]),dtype='float')
    plt.plot(off[off>-1000])
"""    
im2={}
for ib in range(len(bolo)):
    im2[ib]=hp.ud_grade(hp.read_map('/home1/scratch/jmdeloui/MAP/temp_Nside2048_%s_full_0.fits'%(bolo[ib])),nout)
"""
plt.figure(figsize=(16,8))
for i in range(4):
    hp.mollview(im2[i],cmap='jet',norm='hist',hold=False,sub=(2,2,1+i),title='%d'%(i+1))
    
hp.mollview(ref,cmap='jet',min=-2,max=2,hold=False)
plt.figure(figsize=(16,8))
for i in range(4):
    tmp=im[i]-ref
    tmp[im[i]==hp.UNSEEN]=hp.UNSEEN
    hp.mollview(hp.ud_grade(tmp,32),cmap='jet',min=-2,max=2,hold=False,sub=(2,2,1+i),title='%d - REF'%(i+1))
plt.figure(figsize=(16,8))
for i in range(4):
    tmp=im2[i]-ref
    tmp[im2[i]==hp.UNSEEN]=hp.UNSEEN
    hp.mollview(hp.ud_grade(tmp-np.median(tmp[tmp!=hp.UNSEEN]),32),min=-2,max=2,cmap='jet',hold=False,sub=(2,2,1+i),title='%d - REF'%(i+1))
"""
plt.figure(figsize=(16,8))
for i in range(3):
    print(im2[i].min(),im2[i].max())
    for j in range(3-i):
        tmp=im2[i]-im2[i+j+1]
        tmp[im2[i]==hp.UNSEEN]=hp.UNSEEN
        tmp[im2[i+j+1]==hp.UNSEEN]=hp.UNSEEN
        
        hp.mollview(hp.ud_grade(tmp,64),cmap='jet',min=-0.03,max=0.03,hold=False,sub=(3,3,1+3*i+i+j),title='%d - %d'%(i+1,i+j+2))
plt.show()
