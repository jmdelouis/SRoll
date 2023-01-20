import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
import sys

#nfile_I = '/export/home/tfoulquier/MAP/857_maps/temp_OUT2048_test_rstep10_input_REP6_857GHz_I_crab_full_0.fits'
#nfile_Q = '/export/home/tfoulquier/MAP/857_maps/temp_OUT2048_test_rstep10_input_REP6_857GHz_I_crab_full_1.fits'
#nfile_U = '/export/home/tfoulquier/MAP/857_maps/temp_OUT2048_test_rstep10_input_REP6_857GHz_I_crab_full_2.fits'
nfile_I = '/export/home/tfoulquier/MAP/857_maps/temp_OUT2048_rstep10_857GHz_IQU_fullsky_full_0.fits'
nfile_Q = '/export/home/tfoulquier/MAP/857_maps/temp_OUT2048_rstep10_857GHz_IQU_fullsky_full_1.fits'
nfile_U = '/export/home/tfoulquier/MAP/857_maps/temp_OUT2048_rstep10_857GHz_IQU_fullsky_full_2.fits'

nfile2_I = '/export/home/tfoulquier/MAP/857_maps/temp_OUT2048_test2_rstep1_input_REP6_857GHz_IQU_ODUST_crab_full_0.fits'
nfile2_Q = '/export/home/tfoulquier/MAP/857_maps/temp_OUT2048_test2_rstep1_input_REP6_857GHz_IQU_ODUST_crab_full_1.fits'
nfile2_U = '/export/home/tfoulquier/MAP/857_maps/temp_OUT2048_test2_rstep1_input_REP6_857GHz_IQU_ODUST_crab_full_2.fits'

nfile_COND = '/export/home/tfoulquier/MAP/857_maps/temp_OUT2048_test2_rstep1_input_REP6_857GHz_I_crab_full_COND.fits'

im353_I = hp.read_map('/export/home1/tfoulquier/SRoll20_SkyMap_353psb_full.fits')
im353_Q = hp.read_map('/export/home1/tfoulquier/SRoll20_SkyMap_353psb_full.fits',1)
im353_U = hp.read_map('/export/home1/tfoulquier/SRoll20_SkyMap_353psb_full.fits',2)
#im353ds2 = hp.read_map('/export/home1/tfoulquier/SRoll20_SkyMap_353ds2_full.fits')




print(nfile_I)
"""
valmin = sys.argv[1]
valmax = sys.argv[2]
"""
name_I = nfile_I.split("/")
name_I = name_I[len(name_I)-1]

name_Q = nfile_Q.split("/")
name_Q = name_Q[len(name_Q)-1]

name_U = nfile_U.split("/")
name_U = name_U[len(name_U)-1]

#ph_gal=287.60158
#th_gal=-0.64452
ph_gal=184.557
th_gal=-5

im_I = hp.read_map(nfile_I)
im_Q = hp.read_map(nfile_Q)
im_U = hp.read_map(nfile_U)

im2_I = hp.read_map(nfile2_I)
im2_Q = hp.read_map(nfile2_Q)
im2_U = hp.read_map(nfile2_U)


im_COND = hp.read_map(nfile_COND)

hp.mollview(im_I,cmap='jet',min = 0,max=100,title =name_I )
hp.mollview(im_Q,cmap='jet',min = -3,max=3,title =name_Q )
hp.mollview(im_U,cmap='jet',min = -3,max=3,title =name_U )

"""
hp.gnomview(im_I,cmap='jet',rot=[ph_gal,th_gal],min = 0,max = 100,title =name_I )
hp.gnomview(im_Q,cmap='jet',rot=[ph_gal,th_gal],min = -3,max = 3,title =name_Q )
hp.gnomview(im_U,cmap='jet',rot=[ph_gal,th_gal],min = -3,max = 3,title =name_U )

hp.gnomview(im353_I,cmap='jet',rot=[ph_gal,th_gal],min = -0.001,max = 0.001,title ='im353_I' )
hp.gnomview(im353_Q,cmap='jet',rot=[ph_gal,th_gal],min = -0.001,max = 0.001,title ='im353_Q' )
hp.gnomview(im353_U,cmap='jet',rot=[ph_gal,th_gal],min = -0.001,max = 0.001,title ='im353_U' )



hp.gnomview(im_COND,cmap='jet',rot=[ph_gal,th_gal],min = -10,max = 10,title ='COND' )
"""

plt.show()
