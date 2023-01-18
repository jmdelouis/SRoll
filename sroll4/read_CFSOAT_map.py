

import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
#Read offset correction 
off = np.fromfile('/export/home/tfoulquier/VEC/CFOSAT_LS08_offset_OFF',dtype = 'float')
plt.plot(off[off>-10000])



#Read output sroll Q
imq = hp.read_map("/export/home/tfoulquier/MAP/CFOSAT_CFOSAT_full_0.fits")
#im[im<-1e10]=hp.UNSEEN
imq = hp.ud_grade(imq,64) # nside 64

print(hp.UNSEEN)
hp.mollview(imq,cmap="jet",min = -0.02,max = 0.02,rot=[180,0],title ='sroll result Q') # projection 


#Read output sroll U
imu = hp.read_map("/export/home/tfoulquier/MAP/CFOSAT_CFOSAT_full_1.fits")
#im[im<-1e10]=hp.UNSEEN
imu = hp.ud_grade(imu,64) # nside 64

hp.mollview(imu,cmap="jet",min = -0.02,max = 0.02,rot=[180,0],title ='sroll result U') # projection 


#Read output sroll Q templates
map_temp_q = hp.read_map("/export/home/tfoulquier/MAP/temp_Nside128_CFOSAT_full_0.fits")
#im[im<-1e10]=hp.UNSEEN
map_temp_q = hp.ud_grade(map_temp_q,64) # nside 64

print(hp.UNSEEN)
hp.mollview(map_temp_q,cmap="jet",min = -0.01,max = 0.01,rot=[180,0],title ='sroll result Q templates') # projection 


#Read output sroll U templates
map_temp_u = hp.read_map("/export/home/tfoulquier/MAP/temp_Nside128_CFOSAT_full_1.fits")
#im[im<-1e10]=hp.UNSEEN
map_temp_u = hp.ud_grade(map_temp_u,64) # nside 64

hp.mollview(map_temp_u,cmap="jet",min = -0.001,max = 0.001,rot=[180,0],title ='sroll result U templates') # projection 

#Calcul direction
nout=int(np.sqrt(imq.shape[0]//12))
amp=np.ones([imq.shape[0]])*hp.UNSEEN

idx=np.where(imq>-1E20)[0]
amp[idx]=np.sqrt(imq[idx]**2+imu[idx]**2)
#hp.cartview(amp,cmap="plasma",min = 0.0,max = 0.02,title ='Wave amplitude') # projection
hp.cartview(amp,cmap='plasma',min = 0.0,max = 0.02,rot=[180,0],title ='Wave amplitude') # projection
lout=16
th,ph=hp.pix2ang(lout,np.arange(12*lout*lout))
pidx=hp.ang2pix(nout,th,ph)

lidx=np.where(imq[pidx]>-1E20)[0]
phi=2*np.arctan2(imq,imu)
for ii in lidx:
  dx=amp[pidx[ii]]
  x=[th[ii]+np.sin(phi[pidx[ii]])*dx,th[ii]-np.sin(phi[pidx[ii]])*dx]
  y=[ph[ii]+np.cos(phi[pidx[ii]])*dx,ph[ii]-np.cos(phi[pidx[ii]])*dx]

  hp.projplot(x,y,color='yellow')


"""
TF0=[-0.000292936,6.89498e-05,2.50308e-05,0.00018745]
TF1=[-9.21165e-05,-9.16179e-07,-8.73221e-05,0.000184181]
TF2=[0.00782424,-0.00010907,-0.00313055,-0.00446979]
TF3=[-0.000163835,-1.43946e-05,3.80121e-05,0.000137207]
"""
TF0=[-0.0150525,0.00425564,0.00444107,0.00592793,]
TF1=[-0.0028925,0.0030586,-0.00208858,0.00187862,]
TF2=[0.385975,-0.00643398,-0.149772,-0.226291,]
TF3=[-0.00634197,-0.000265213,0.00213657,0.00442346,]
TF4=[-0.0207998,0.00427556,0.00574678,0.0103579,]
TF5=[-0.0280448,0.0142401,0.00402221,0.00952115,]
TF6=[0.351053,-0.00378559,-0.13075,-0.204799,]
TF7=[-0.017052,0.000979414,0.00559888,0.0098503,]



plt.figure()
x = np.arange(128)/128*2*np.pi
for i in range(4):
  plt.plot(x,TF0[i]*np.cos(x)+TF1[i]*np.sin(x)+TF2[i]*np.cos(2*x)+TF3[i]*np.sin(2*x)+TF4[i]*np.cos(3*x)+TF5[i]*np.sin(3*x)+TF6[i]*np.cos(4*x)+TF7[i]*np.sin(4*x),label='LS%02d'%(4+2*i)) 
plt.legend()


#Print data without denoising
map1 = hp.read_map("/export/home1/jmdeloui/CFOSAT/CFOSAT_Full_Q.fits")
hp.mollview(map1,cmap="jet",min = -0.02,max = 0.02,rot=[180,0],title ='Map no denoised Q') # projection 

map2 = hp.read_map("/export/home1/jmdeloui/CFOSAT/CFOSAT_Full_U.fits")
hp.mollview(map2,cmap="jet",min = -0.02,max = 0.02,rot=[180,0],title ='Map no denoised U') # projection 


"""
map_fake = hp.read_map("/export/home/jmdeloui/reduced_FAKE_U")
map_fake = hp.ud_grade(map_fake,64) # nside 64
hp.mollview(map_fake,cmap="jet",min = -0.0003,max = 0.0003,title ='reduced_FAKE_U') # projection 


map_fake = hp.read_map("/export/home/jmdeloui/reduced_FAKE_Q")
map_fake = hp.ud_grade(map_fake,64) # nside 64
hp.mollview(map_fake,cmap="jet",min = -0.0003,max = 0.0003,title ='reduced_FAKE_Q') # projection 
"""

plt.show()


