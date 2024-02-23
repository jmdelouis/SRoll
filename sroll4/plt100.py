import healpy as hp
import numpy as np
import matplotlib.pyplot as plt


imc=hp.read_map('/travail/jdelouis/sroll4_map/100ds1_Full_COND.fits')
hp.mollview(imc,cmap='jet')

vv=np.fromfile('/travail/jdelouis/sroll4_map/100ds1_VEC_100-1a_OFF',dtype='float')
plt.figure()
plt.plot(vv[vv>-10000])

imi=1E6*hp.read_map('/travail/jdelouis/sroll4_map/100ds1_Full_0.fits')
imq=1E6*hp.read_map('/travail/jdelouis/sroll4_map/100ds1_Full_1.fits')
imu=1E6*hp.read_map('/travail/jdelouis/sroll4_map/100ds1_Full_2.fits')

imir=1E6*hp.read_map('/travail/jdelouis/sroll4_map/100ds1_Full_0_RAW.fits')
imqr=1E6*hp.read_map('/travail/jdelouis/sroll4_map/100ds1_Full_1_RAW.fits')
imur=1E6*hp.read_map('/travail/jdelouis/sroll4_map/100ds1_Full_2_RAW.fits')

imir=imir-np.median(imi)
imi=imi-np.median(imi)

plt.figure(figsize=(8,12))
hp.mollview(imi,cmap='jet',hold=False,sub=(3,2,1),min=-300,max=300)
hp.mollview(imq,cmap='jet',hold=False,sub=(3,2,3),min=-30,max=30)
hp.mollview(imu,cmap='jet',hold=False,sub=(3,2,5),min=-30,max=30)

hp.mollview(imir,cmap='jet',hold=False,sub=(3,2,2),min=-300,max=300)
hp.mollview(imqr,cmap='jet',hold=False,sub=(3,2,4),min=-30,max=30)
hp.mollview(imur,cmap='jet',hold=False,sub=(3,2,6),min=-30,max=30)
plt.show()
