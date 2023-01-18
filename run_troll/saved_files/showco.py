import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
import pickle
from numpy import linalg as LA


nside=128
coref=np.fromfile('/scratch/cnt0028/ias1717/SHARED/RD12_data/dmc/MISS03/DATA/sroll_templates/MAP_JMD_128/COSIMU.float32.bin',dtype='float32')
#cl,alm=hp.anafast(coref,alm=True)
#m,l=hp.Alm.getlm(len(cl)-1)
#alm=alm*np.exp(8*(l/float(len(cl)))**2)*np.exp(-4*(l/float(len(cl)))**4)
##coref2=hp.alm2map(alm,128)*(coref>1E-6)
#coef1=hp.alm2map(alm,128,fwhm=0.1)
#coef2=hp.alm2map(alm,128,fwhm=0.03)
#coref2=abs(hp.synfast(cl,128)*coef2)
#coref2=coef1+coref2*coef1.std()/coref2.std()
#
#coref13=coref2
#cldiff1=hp.anafast((coref13)*1E6)
#tcldiff1=hp.anafast((coref)*1E6)
#plt.plot(np.log(cldiff1),'--',color='blue')
#plt.plot(np.log(tcldiff1),'--',color='red')
#plt.show()
#exit()
#coref2=abs(hp.synfast(cl*(1+np.arange(384)/384)**2,128)*(coref>1E-6))
#coref2=coref2*coref.std()/coref2.std()
#(coref2.astype('float32')).tofile('/scratch/cnt0028/ias1717/SHARED/RD12_data/dmc/MISS03/DATA/sroll_templates/MAP_JMD_128/NEW_cleaned_13CO_SIMU.float32.bin')

coref13=np.fromfile('/scratch/cnt0028/ias1717/SHARED/RD12_data/dmc/MISS03/DATA/sroll_templates/MAP_JMD_128/NEW_cleaned_13CO_SIMU.float32.bin',dtype='float32')
mat=np.fromfile('matrix.dat',dtype='float').reshape(12*nside*nside,5,5)
vec=np.fromfile('vec.dat',dtype='float').reshape(12*nside*nside,5)
print('SAVE DONE')

cond=np.zeros([12*nside*nside])
for i in range(12*nside*nside):
    cond[i]=LA.cond(mat[i,:,:])

vec2=1*vec
for i in range(12*nside*nside):
    if cond[i]<1E10:
        vec2[i,:]=LA.solve(mat[i,:,:],vec[i,:])
    else:
        vec2[i,:]=hp.UNSEEN
    
vec2[:,0]-=np.median(vec2[:,0])

cldiff=hp.anafast((vec2[:,4]-coref13)*1E6)
cldiff1=hp.anafast((coref13)*1E6)
tcldiff=hp.anafast((vec2[:,3]-coref)*1E6)
tcldiff1=hp.anafast((coref)*1E6)
plt.plot(np.log(cldiff1),'--',color='blue')
plt.plot(np.log(cldiff),color='blue')
plt.plot(np.log(tcldiff1),'--',color='red')
plt.plot(np.log(tcldiff),color='red')
plt.show()

vec2.tofile('/scratch/cnt0028/ias1717/SHARED/RD12_data/dmc/MISS03/DATA/sroll_templates/MAP_JMD_128/100x143_RESULT.float32.bin')
tab1=[-300,-30,-30,-5,-5]
tab2=[300,30,30,30,30]
for i in range(5):
    hp.mollview(vec2[:,i]*1E6,min=tab1[i],max=tab2[i],cmap='jet')
plt.show()
exit()
