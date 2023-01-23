import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
import sys

nfileQ = '/home1/datawork/tfoulqui/map_run_857/NsideGV2048_test_rstep2_857GHz_ODUST_IQU_full_1.fits'
nfileU = '/home1/datawork/tfoulqui/map_run_857/NsideGV2048_test_rstep2_857GHz_ODUST_IQU_full_2.fits'


nfileQ2 = '/home1/datawork/tfoulqui/map_run_857/NsideGV2048_test_rstep2_857GHz_PDUST_ODUST_GRAD_IQU_full_1.fits'
nfileU2 = '/home1/datawork/tfoulqui/map_run_857/NsideGV2048_test_rstep2_857GHz_PDUST_ODUST_GRAD_IQU_full_2.fits'

nfileQ3 = '/home1/datawork/tfoulqui/map_run_857/NsideGV2048_test_rstep2_857GHz_ODUST_GRAD_IQU_full_1.fits'
nfileU3 = '/home1/datawork/tfoulqui/map_run_857/NsideGV2048_test_rstep2_857GHz_ODUST_GRAD_IQU_full_2.fits'

nfileQ4 = '/home1/datawork/tfoulqui/map_run_857/NsideGV2048_test_rstep2_857GHz_GRAD_IQU_full_1.fits'
nfileU4 = '/home1/datawork/tfoulqui/map_run_857/NsideGV2048_test_rstep2_857GHz_GRAD_IQU_full_2.fits'

valmin = -1
valmax = 1

nside = 64

im_Q = hp.ud_grade(hp.read_map(nfileQ),nside)
im2_Q =hp.ud_grade(hp.read_map(nfileQ2),nside)

im_U = hp.ud_grade(hp.read_map(nfileU),nside)
im2_U =hp.ud_grade(hp.read_map(nfileU2),nside)

im3_Q =hp.ud_grade(hp.read_map(nfileQ3),nside)
im3_U =hp.ud_grade(hp.read_map(nfileU3),nside)


im4_Q =hp.ud_grade(hp.read_map(nfileQ4),nside)
im4_U =hp.ud_grade(hp.read_map(nfileU4),nside)

hp.mollview(im_Q,cmap='jet',norm ='hist',title = 'Q_ODUST' )
hp.mollview(im2_Q,cmap='jet',norm ='hist',title = 'Q_GRAD_PDUST_ODUST' )

hp.mollview(im_U,cmap='jet',norm ='hist',title = 'U_ODUST' )
hp.mollview(im2_U,cmap='jet',norm ='hist',title = 'U_GRAD_PDUST_ODUST' )

hp.mollview(im3_U,cmap='jet',norm ='hist',title = 'U_GRAD_ODUST' )
hp.mollview(im3_Q,cmap='jet',norm ='hist',title = 'U_GRAD_ODUST' )


hp.mollview(im4_U,cmap='jet',norm ='hist',title = 'U_GRAD' )
hp.mollview(im4_Q,cmap='jet',norm ='hist',title = 'U_GRAD' )

plt.show()


## saved map path :
# result  rstep2_857GHz_GRAD_ODUST /home1/datawork/tfoulqui/map_run_857/NsideGV2048_test_rstep1_857GHz_GRAD_PDUST_ODUST_IQU_full_0.fits



