import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import sys
nout =32


def down(im,nside):
    nin=int(np.sqrt(im.shape[0]/12))
    if nin==nside:
        return(im)
    nn=(nin//nside)**2
    return(np.mean(im.reshape(12*nin**2//nn,nn),1))

plt.figure()
corr2=np.load('/export/home/jmdeloui/heal_cnn/FSLDATA.npy')[0:4*256*256].reshape(256,256,4)
corr=np.load('corrfsl.npy').reshape(256,256,4)
for i in range(4):
    plt.subplot(3,4,1+i)
    plt.imshow(corr[:,:,i],cmap='jet')
    plt.subplot(3,4,5+i)
    plt.imshow(corr2[:,:,i],cmap='jet',vmin=-0.01,vmax=0.3)
    plt.subplot(3,4,9+i)
    plt.imshow(corr[:,:,i]-corr2[:,:,i],cmap='jet',vmin=-0.1,vmax=0.1)

plt.figure()

poldustmol=np.zeros([12*nout*nout,2])
poldustmol[:,0] = down(np.load('/export/home/jmdeloui/heal_cnn/map_353_1_256_nest.npy'),nout)*2000*1.85
poldustmol[:,1] = down(np.load('/export/home/jmdeloui/heal_cnn/map_353_2_256_nest.npy'),nout)*2000*1.85

map=np.load('myMap_out.npy').reshape(12*nout*nout,2)
print(map.shape)
hp.mollview(map[:,0],cmap='jet',hold=False,sub=(3,2,1),nest=True,title='MAP_%d.npy'%(nout))
hp.mollview(map[:,1],cmap='jet',hold=False,sub=(3,2,2),nest=True,title='MAP_%d.npy'%(nout))
hp.mollview(poldustmol[:,0],min=-10,max=10,cmap='jet',hold=False,sub=(3,2,3),nest=True,title='MAP_%d.npy'%(nout))
hp.mollview(poldustmol[:,1],min=-10,max=10,cmap='jet',hold=False,sub=(3,2,4),nest=True,title='MAP_%d.npy'%(nout))
hp.mollview(map[:,0]-poldustmol[:,0],min=-0.3,max=0.3,cmap='jet',hold=False,sub=(3,2,5),nest=True,title='MAP_%d.npy'%(nout))
hp.mollview(map[:,1]-poldustmol[:,1],min=-0.3,max=0.3,cmap='jet',hold=False,sub=(3,2,6),nest=True,title='MAP_%d.npy'%(nout))








plt.show()
