
import healpy as hp
import numpy as np
import matplotlib.pyplot as plt

im = hp.read_map("/export/home/tfoulquier/MAP/857_maps/temp_OUT2048_857GHz_full_0.fits")
#im = hp.read_map("/export/home/tfoulquier/MAP/857_maps/temp_OUT2048_carina_857GHz_full_0.fits")
hp.mollview(im,cmap="jet",min = -10,max = 10,title ='carina result 0') 

mask = np.fromfile('/export/home1/jmdeloui/CARINA_MASK.float32.bin',dtype=np.float32)


ph_gal=287.60158
th_gal=-0.64452

hp.gnomview(im,cmap='jet',rot=[ph_gal,th_gal],min=-0.01,max=400)
hp.mollview(im,cmap='jet',title="Carina",norm='hist')


hp.gnomview(mask,cmap='jet',rot=[ph_gal,th_gal])
hp.mollview(mask,cmap='jet',title="Mask Carina")


plt.show()
