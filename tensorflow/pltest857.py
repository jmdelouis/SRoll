import healpy as hp
import numpy as np
import matplotlib.pyplot as plt

nout=64
amp=0.03
beam=1/180.*np.pi

co=np.fromfile('/scratch/cnt0028/ias1717/SHARED/RD12_data/dmc/MISS03/DATA/sroll_templates/MAP_JMD_128/NEW_cleaned_12CO_K.float32.bin',dtype='float32')
co=hp.ud_grade(co,nout)

#new=['CNN2D_TL_FSL_1_-32','NORM_FSL_4_-2','NORM_FSL_4_-2','NORM_FSL_4_-2']
new=['CNN2D_TL_FSL_1_-32','CNN2D_TL_FSL_1_-32','CNN2D_TL_FSL_1_-32','CNN2D_TL_FSL_1_-32']
#new=['CNN2D_TL_FSL_1_-32','NORM_FSL_4_-2','CNN2D_TL_FSL_1_-32','NORM_FSL_4_-2']

plt.figure(figsize=(8,6))
for i in range(4):
    im1=hp.ud_grade(hp.read_map('/scratch/cnt0028/ias1717/SHARED/bware/sroll2r94_jmd/545/MAP/857-%d_RD12_REP6_DATA_%s_Full_I.fits'%(i+1,new[i])),nout)
    im2=hp.ud_grade(hp.read_map('/scratch/cnt0028/ias1717/SHARED/bware/sroll2r94_jmd/545/MAP/857-%d_RD12_REP6_DATA_NORM_FSL_4_-2_Full_I.fits'%(i+1)),nout)
    hp.mollview(im1-im2-np.median(im1-im2),cmap='jet',min=-amp,max=amp,hold=False,sub=(2,2,1+i),title='CNN-NORM 857-%d'%(i+1))
plt.savefig('FSLCORRECTION.png')
im={}
for i in range(4):
    im[i]=hp.ud_grade(hp.read_map('/scratch/cnt0028/ias1717/SHARED/bware/sroll2r94_jmd/545/MAP/857-%d_RD12_REP6_DATA_%s_Full_I.fits'%(i+1,new[i])),nout)

cl1={}
plt.figure(figsize=(8,6))
for i in range(3):
    for j in range(3-i):
        mat=np.zeros([3,3])
        vec=np.zeros([3])

        mat[0,0]=np.sum(im[i+j+1]*im[i+j+1])
        mat[1,0]=np.sum(co*im[i+j+1])
        mat[2,0]=np.sum(   im[i+j+1])
        mat[0,1]=np.sum(im[i+j+1]*co)
        mat[1,1]=np.sum(co*co)
        mat[2,1]=np.sum(   co)
        mat[0,2]=np.sum(im[i+j+1])
        mat[1,2]=np.sum(co)
        mat[2,2]=np.sum(1+0*co)

        vec[0]=np.sum(im[i+j+1]*im[i])
        vec[1]=np.sum(co*im[i])
        vec[2]=np.sum(im[i])

        aa=np.linalg.solve(mat,vec)
        print(i,i+j+1,aa)
        tmp=hp.smoothing(im[i]-aa[0]*im[i+j+1]-aa[1]*co-aa[2],beam)
        hp.mollview(tmp-np.median(tmp),cmap='jet',min=-amp,max=amp,hold=False,sub=(3,3,1+3*i+j+i),title='CNN %d-%d'%(i+1,i+j+2))
        cl1[1+3*i+j+i]=hp.anafast(tmp-np.median(tmp),gal_cut=20)
plt.savefig('FSLMAPCNN.png')

im={}
im[0]=hp.ud_grade(hp.read_map('/scratch/cnt0028/ias1717/SHARED/bware/sroll2r94_jmd/545/MAP/857-1_RD12_REP6_DATA_NORM_FSL_4_-2_Full_I.fits'),nout)
im[1]=hp.ud_grade(hp.read_map('/scratch/cnt0028/ias1717/SHARED/bware/sroll2r94_jmd/545/MAP/857-2_RD12_REP6_DATA_NORM_FSL_4_-2_Full_I.fits'),nout)
im[2]=hp.ud_grade(hp.read_map('/scratch/cnt0028/ias1717/SHARED/bware/sroll2r94_jmd/545/MAP/857-3_RD12_REP6_DATA_NORM_FSL_4_-2_Full_I.fits'),nout)
im[3]=hp.ud_grade(hp.read_map('/scratch/cnt0028/ias1717/SHARED/bware/sroll2r94_jmd/545/MAP/857-4_RD12_REP6_DATA_NORM_FSL_4_-2_Full_I.fits'),nout)
cl2={}
plt.figure(figsize=(8,6))
for i in range(3):
    for j in range(3-i):
        mat=np.zeros([3,3])
        vec=np.zeros([3])

        mat[0,0]=np.sum(im[i+j+1]*im[i+j+1])
        mat[1,0]=np.sum(co*im[i+j+1])
        mat[2,0]=np.sum(   im[i+j+1])
        mat[0,1]=np.sum(im[i+j+1]*co)
        mat[1,1]=np.sum(co*co)
        mat[2,1]=np.sum(   co)
        mat[0,2]=np.sum(im[i+j+1])
        mat[1,2]=np.sum(co)
        mat[2,2]=np.sum(1+0*co)

        vec[0]=np.sum(im[i+j+1]*im[i])
        vec[1]=np.sum(co*im[i])
        vec[2]=np.sum(im[i])

        aa=np.linalg.solve(mat,vec)
        print(i,i+j+1,aa)
        tmp=hp.smoothing(im[i]-aa[0]*im[i+j+1]-aa[1]*co-aa[2],beam)
        hp.mollview(tmp-np.median(tmp),cmap='jet',min=-amp,max=amp,hold=False,sub=(3,3,1+3*i+j+i),title='TEMP %d-%d'%(i+1,i+j+2))
        cl2[1+3*i+j+i]=hp.anafast(tmp-np.median(tmp),gal_cut=20)
plt.savefig('FSLMAPREF.png')

plt.figure(figsize=(8,6))
for i in range(3):
    for j in range(3-i):
        plt.subplot(3,3,1+3*i+j+i)
        plt.plot(cl1[1+3*i+j+i],color='red',label='CNN')
        plt.plot(cl2[1+3*i+j+i],color='blue',label='TEMP')
        plt.yscale('log')
        plt.xscale('log')
        plt.legend()
plt.savefig('FSLSPEC.png')
plt.show()
        
