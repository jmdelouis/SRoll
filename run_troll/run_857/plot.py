import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
import sys

nfile = sys.argv[1]
print(nfile)
valmin = sys.argv[2]
valmax = sys.argv[3]
name = nfile.split("/")
name = name[len(name)-1]

ph_gal=287.60158
th_gal=-0.64452

ph_gal_1=184.557
th_gal_1=-5

im = hp.read_map(nfile)
im2 = hp.ud_grade(im,64)
hp.mollview(im2,cmap='jet',min = valmin,max=valmax,title =name )
hp.gnomview(im,cmap='jet',rot=[ph_gal,th_gal],min = valmin,max = valmax,title ='Carina')
hp.gnomview(im,cmap='jet',rot=[ph_gal_1,th_gal_1],min = valmin,max = valmax,title ='Crab')

plt.show()


## saved map path :
# result  rstep2_857GHz_GRAD_ODUST /home1/datawork/tfoulqui/map_run_857/NsideGV2048_test_rstep1_857GHz_GRAD_PDUST_ODUST_IQU_full_0.fits
