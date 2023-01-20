
import healpy as hp
import numpy as np
import matplotlib.pyplot as plt

carina = True

ph_gal=287.60158
th_gal=-0.64452

#mask = np.fromfile('/export/home1/jmdeloui/CARINA_MASK.float32.bin',dtype=np.float32)

#353Ghz planck maps
"""
im353ds1 = hp.read_map('/export/home1/tfoulquier/SRoll20_SkyMap_353ds1_full.fits')
im353ds2 = hp.read_map('/export/home1/tfoulquier/SRoll20_SkyMap_353ds2_full.fits')
"""
#I REP6
#im0 = hp.read_map('/export/home/tfoulquier/MAP/857_maps/temp_OUT2048_rstep10_input_REP6_857GHz_I_full_0.fits')


pdust_im = hp.read_map('/export/home/tfoulquier/MAP/857_maps/temp_OUT2048_test_rstep10_input_PDUST_857GHz_IQU_full_0.fits')
pdust_im1 = hp.read_map('/export/home/tfoulquier/MAP/857_maps/temp_OUT2048_test_rstep10_input_PDUST_857GHz_IQU_full_1.fits')
pdust_im2 = hp.read_map('/export/home/tfoulquier/MAP/857_maps/temp_OUT2048_test_rstep10_input_PDUST_857GHz_IQU_full_2.fits')


rep6_im = hp.read_map('/export/home/tfoulquier/MAP/857_maps/temp_OUT2048_test_rstep1_input_REP6_857GHz_IQU_carina_delta_full_0.fits')
rep6_im1 = hp.read_map('/export/home/tfoulquier/MAP/857_maps/temp_OUT2048_test_rstep1_input_REP6_857GHz_IQU_carina_delta_full_1.fits')
rep6_im2 = hp.read_map('/export/home/tfoulquier/MAP/857_maps/temp_OUT2048_test_rstep1_input_REP6_857GHz_IQU_carina_delta_full_2.fits')
rep6_im3 = hp.read_map('/export/home/tfoulquier/MAP/857_maps/temp_OUT2048_test_rstep1_input_REP6_857GHz_IQU_carina_delta_full_COND.fits')

rep6_im1 -= np.median(rep6_im1[rep6_im1>-1E10])
rep6_im2 -= np.median(rep6_im2[rep6_im2>-1E10])


hp.gnomview(pdust_im,cmap='jet',rot=[ph_gal,th_gal],min = 0,max = 10,title ='pdust gnomview result I') 
hp.gnomview(pdust_im1,cmap='jet',rot=[ph_gal,th_gal],min = -0.3,max = 0.3,title ='pdust gnomview result Q') 
hp.gnomview(pdust_im2,cmap='jet',rot=[ph_gal,th_gal],min = -0.3,max = 0.3,title ='pdust gnomview result U') 

hp.gnomview(rep6_im,cmap='jet',rot=[ph_gal,th_gal],min = -10,max = 1000,title ='gnomview result I') 
hp.gnomview(rep6_im1,cmap='jet',rot=[ph_gal,th_gal],min = -0.3,max = 0.3,title ='gnomview result Q') 
hp.gnomview(rep6_im2,cmap='jet',rot=[ph_gal,th_gal],min = -0.3,max = 0.3,title ='gnomview result U') 
hp.gnomview(rep6_im3,cmap='jet',rot=[ph_gal,th_gal],min = -0.3,max = 0.3,title ='gnomview result COND') 



"""
rep6grad2_im1 = hp.read_map('/export/home/tfoulquier/MAP/857_maps/temp_OUT2048_test_rstep1_input_REP6_857GHz_IQU_carina_GRAD2_full_1.fits')
rep6grad2_im2 = hp.read_map('/export/home/tfoulquier/MAP/857_maps/temp_OUT2048_test_rstep1_input_REP6_857GHz_IQU_carina_GRAD2_full_2.fits')

rep6grad3_im1 = hp.read_map('/export/home/tfoulquier/MAP/857_maps/temp_OUT2048_test_rstep1_input_REP6_857GHz_IQU_carina_GRAD3_full_1.fits')
rep6grad3_im2 = hp.read_map('/export/home/tfoulquier/MAP/857_maps/temp_OUT2048_test_rstep1_input_REP6_857GHz_IQU_carina_GRAD3_full_2.fits')

rep6grad2_im1 -= np.median(rep6grad2_im1[rep6grad2_im1>-1E10])
rep6grad2_im2 -= np.median(rep6grad2_im2[rep6grad2_im2>-1E10])

rep6grad3_im1 -= np.median(rep6grad3_im1[rep6grad3_im1>-1E10])
rep6grad3_im2 -= np.median(rep6grad3_im2[rep6grad3_im2>-1E10])


hp.gnomview(rep6grad2_im1,cmap='jet',rot=[ph_gal,th_gal],min = -0.3,max = 0.3,title ='gnomview result Q GRAD2') 
hp.gnomview(rep6grad2_im2,cmap='jet',rot=[ph_gal,th_gal],min = -0.3,max = 0.3,title ='gnomview result U GRAD2')

hp.gnomview(rep6grad3_im1,cmap='jet',rot=[ph_gal,th_gal],min = -0.3,max = 0.3,title ='gnomview result Q GRAD3') 
hp.gnomview(rep6grad3_im2,cmap='jet',rot=[ph_gal,th_gal],min = -0.3,max = 0.3,title ='gnomview result U GRAD3')  
"""

plt.show()
exit(0)
if(carina == True):
  #IQU REP6 mask carina + badring2
  im = hp.read_map('/export/home/tfoulquier/MAP/857_maps/temp_OUT2048_rstep10_input_REP6_857GHz_IQU_maskCarina_badring2_full_0.fits')
  im1 = hp.read_map('/export/home/tfoulquier/MAP/857_maps/temp_OUT2048_rstep10_input_REP6_857GHz_IQU_maskCarina_badring2_full_1.fits')
  im2 = hp.read_map('/export/home/tfoulquier/MAP/857_maps/temp_OUT2048_rstep10_input_REP6_857GHz_IQU_maskCarina_badring2_full_2.fits')
  im3 = hp.read_map('/export/home/tfoulquier/MAP/857_maps/temp_OUT2048_rstep10_input_REP6_857GHz_IQU_maskCarina_badring2_full_COND.fits')

  im_dust = hp.read_map('/export/home/tfoulquier/MAP/857_maps/temp_OUT2048_rstep1_input_PDUST_857GHz_IQU_maskCarina_badring2_full_0.fits')
  im_dust2 = hp.read_map('/export/home/tfoulquier/MAP/857_maps/temp_OUT2048_rstep1_input_PDUST_857GHz_IQU_maskCarina_badring2_full_1.fits')
  im_dust3 = hp.read_map('/export/home/tfoulquier/MAP/857_maps/temp_OUT2048_rstep1_input_PDUST_857GHz_IQU_maskCarina_badring2_full_2.fits')
  im_dust4 = hp.read_map('/export/home/tfoulquier/MAP/857_maps/temp_OUT2048_rstep1_input_PDUST_857GHz_IQU_maskCarina_badring2_full_COND.fits')


  
  
  hp.mollview(im0,cmap="jet",min = -10,max = 10,title ='result I - rstep10_input_REP6_857GHz_I_full_0') 
  hp.mollview(im,cmap="jet",min = -10,max = 10,title ='result I - rstep10_input_REP6_857GHz_IQU_maskCarina_badring2_full_0') 
  hp.mollview(im1,cmap="jet",min = -10,max = 10,title ='result Q - rstep10_input_REP6_857GHz_IQU_maskCarina_badring2_full_1') 
  hp.mollview(im2,cmap="jet",min = -10,max = 10,title ='result U - rstep10_input_REP6_857GHz_IQU_maskCarina_badring2_full_2') 
  hp.mollview(im3,cmap="jet",min = -10,max = 10,title ='result COND - rstep10_input_REP6_857GHz_IQU_maskCarina_badring2_full_COND')

  
  hp.mollview(im_dust,cmap="jet",min = -1E14,max = 1E14,title ='result I - rstep10_input_PDUST_857GHz_IQU_full_0') 
  hp.mollview(im_dust2,cmap="jet",min = -1E14,max = 1E14,title ='result Q - rstep10_input_PDUST_857GHz_IQU_full_1') 
  hp.mollview(im_dust3,cmap="jet",min = -1E14,max = 1E14,title ='result U - rstep10_input_PDUST_857GHz_IQU_full_2') 
  hp.mollview(im_dust4,cmap="jet",min = -1E14,max = 1E14,title ='result COND - rstep10_input_PDUST_857GHz_IQU_full_COND')


  hp.gnomview(im,cmap='jet',rot=[ph_gal,th_gal],min = 0,max = 10,title ='gnomview result I') 
  hp.gnomview(im1,cmap='jet',rot=[ph_gal,th_gal],min = -10,max = 10,title ='gnomview result Q') 
  hp.gnomview(im2,cmap='jet',rot=[ph_gal,th_gal],min = -10,max = 10,title ='gnomview result U') 


else:

  #IQU PDUST mask 857 + normal badring
  im_dust = hp.read_map('/export/home/tfoulquier/MAP/857_maps/temp_OUT2048_rstep10_input_PDUST_857GHz_IQU_full_0.fits')
  im_dust2 = hp.read_map('/export/home/tfoulquier/MAP/857_maps/temp_OUT2048_rstep10_input_PDUST_857GHz_IQU_full_1.fits')
  im_dust3 = hp.read_map('/export/home/tfoulquier/MAP/857_maps/temp_OUT2048_rstep10_input_PDUST_857GHz_IQU_full_2.fits')
  im_dust4 = hp.read_map('/export/home/tfoulquier/MAP/857_maps/temp_OUT2048_rstep10_input_PDUST_857GHz_IQU_full_COND.fits')

  
  hp.mollview(im_dust,cmap="jet",min = -1E14,max = 1E14,title ='result I - rstep10_input_PDUST_857GHz_IQU_full_0') 
  hp.mollview(im_dust2,cmap="jet",min = -1E14,max = 1E14,title ='result Q - rstep10_input_PDUST_857GHz_IQU_full_1') 
  hp.mollview(im_dust3,cmap="jet",min = -1E14,max = 1E14,title ='result U - rstep10_input_PDUST_857GHz_IQU_full_2') 
  hp.mollview(im_dust4,cmap="jet",min = -1E14,max = 1E14,title ='result COND - rstep10_input_PDUST_857GHz_IQU_full_COND')




hp.mollview(im353ds1,cmap="jet",min = -0.01,max = 0.01,title ='im353ds1') 
hp.mollview(im353ds2,cmap="jet",min = -0.01,max = 0.01,title ='im353ds2') 



plt.show()



#hp.mollview(im3,cmap="jet",min = -10,max = 10,title ='I 857-input rstep = 10') 
#hp.mollview(im4,cmap="jet",min = -1E13,max = 1E14,title ='I PDUST-input rstep = 10') 
#hp.mollview(im353ds1,cmap="jet",min = -0.01,max = 0.01,title ='im353ds1') 
#hp.mollview(im353ds2,cmap="jet",min = -0.01,max = 0.01,title ='im353ds2') 
#hp.mollview(im,cmap='jet',title="Carina",norm='hist')
#hp.gnomview(mask,cmap='jet',rot=[ph_gal,th_gal])
#hp.mollview(mask,cmap='jet',title="Mask Carina")

