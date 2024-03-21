import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
import foscat.Spline1D as spl 


imc=hp.read_map('/travail/jdelouis/sroll4_map/100-1A_Full_COND.fits')
hp.mollview(imc,cmap='jet')

vv=np.fromfile('/travail/jdelouis/sroll4_map/100-1A_VEC_100-1a_X2',dtype='float')

# le premir terme est le fit de gain
nspline=vv.shape[0]-1
splinetime=spl.Spline1D(nspline)
ntot=16*nspline
res=np.zeros([ntot])
for i in range(ntot):
    res[i]=(splinetime.calculate(i/ntot)*vv[1:]).sum()

plt.figure()
plt.plot(res,color='b',label='Spline Val')
plt.legend()

imi=1E6*hp.read_map('/travail/jdelouis/sroll4_map/100-1A_Full_0.fits')

imir=1E6*hp.read_map('/travail/jdelouis/sroll4_map/100-1A_Full_0_RAW.fits')

imir=imir-np.median(imi)
imi=imi-np.median(imi)

plt.figure(figsize=(8,8))
hp.mollview(imi,cmap='jet',hold=False,sub=(2,1,1),min=-1000,max=1000)
hp.mollview(imir,cmap='jet',hold=False,sub=(2,1,2),min=-1000,max=1000)
plt.show()
