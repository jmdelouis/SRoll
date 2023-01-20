import numpy as np 
import healpy as hp 

# th ph position - get_inter_val ->  healpy 

in_template_map = np.fromfile("/export/home1/jmdeloui/DATA4SROLL4/map_857_2018.float32.bin",dtype='float32')
mask = np.fromfile("/export/home1/jmdeloui/CARINA_MASK.float32.bin",dtype='float32')
th,ph = hp.pix2ang(4096,np.arange(12*4096**2))

im = hp.get_interp_val(in_template_map,th,ph).astype('float32')
im2 = hp.get_interp_val(mask,th,ph).astype('float32')
null_map = np.arange(0,12*4096**2,dtype='float32')

im.tofile('/export/home1/tfoulquier/data_sroll/MAP/template_857_4096.float32.bin')
null_map.tofile('/export/home1/tfoulquier/data_sroll/MAP/map_null_857_4096.float32.bin')
im2.tofile('/export/home1/tfoulquier/data_sroll/MAP/mask_carina_4096.float32.bin')




