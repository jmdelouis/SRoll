
import healpy as hp
import matplotlib.pyplot as plt
import numpy as np

#===========================================================================
# SCRIPT TO GENERATE TEMPLATE FOR SROLL
#===========================================================================


def up_grade(im,onside,FWHM=20./60./180.*np.pi):
        th,ph=hp.pix2ang(onside,np.arange(12*onside*onside))
        val=hp.get_interp_val(im,th,ph)
        val=hp.smoothing(val,FWHM)
        return(val)


iname=['/scratch/cnt0028/ias1717/SHARED/RD12_data/dmc/MISS03/DATA/MAP_JMD_128/freefree_NORM_I',
       '/scratch/cnt0028/ias1717/SHARED/RD12_data/dmc/MISS03/DATA/sroll_templates/MAP_JMD_128/353_Q.float32.bin',
       '/scratch/cnt0028/ias1717/SHARED/RD12_data/dmc/MISS03/DATA/sroll_templates/MAP_JMD_128/353_U.float32.bin',
       '/scratch/cnt0028/ias1717/SHARED/RD12_data/dmc/MISS03/DATA/sroll_templates/MAP_JMD_128/NEW_cleaned_12CO.float32.bin',
       '/scratch/cnt0028/ias1717/SHARED/RD12_data/dmc/MISS03/DATA/sroll_templates/MAP_JMD_128/NEW_cleaned_13CO.float32.bin',
       '/scratch/cnt0028/ias1717/SHARED/RD12_data/dmc/MISS03/DATA/MAP_JMD_128/Dust_New_I',
       '/scratch/cnt0028/ias1717/SHARED/RD12_data/dmc/MISS03/DATA/MAP_JMD_128/Dust_New_Q',
       '/scratch/cnt0028/ias1717/SHARED/RD12_data/dmc/MISS03/DATA/MAP_JMD_128/Dust_New_U']

oname=['/scratch/cnt0028/ias1717/SHARED/bware/JMD_TEMPLATE/MAP_JMD_1024/freefree_NORM_I',
       '/scratch/cnt0028/ias1717/SHARED/bware/JMD_TEMPLATE/MAP_JMD_1024/353pol_Q.float32.bin',
       '/scratch/cnt0028/ias1717/SHARED/bware/JMD_TEMPLATE/MAP_JMD_1024/353pol_U.float32.bin',
       '/scratch/cnt0028/ias1717/SHARED/bware/JMD_TEMPLATE/MAP_JMD_1024/NEW_cleaned_12CO.float32.bin',
       '/scratch/cnt0028/ias1717/SHARED/bware/JMD_TEMPLATE/MAP_JMD_1024/NEW_cleaned_13CO.float32.bin',
       '/scratch/cnt0028/ias1717/SHARED/bware/JMD_TEMPLATE/MAP_JMD_1024/Dust_New_I',
       '/scratch/cnt0028/ias1717/SHARED/bware/JMD_TEMPLATE/MAP_JMD_1024/Dust_New_Q',
       '/scratch/cnt0028/ias1717/SHARED/bware/JMD_TEMPLATE/MAP_JMD_1024/Dust_New_U']


for i in range(8):
        print('Start compute ',iname[i])
        im=np.fromfile(iname[i],dtype='float32')
        rim=up_grade(im,1024)
        rim=rim.astype('float32')
        rim.tofile(oname[i])
        print('End compute ',oname[i])
        
