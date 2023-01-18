

import healpy as hp
import numpy as np
import matplotlib.pyplot as plt

off = np.fromfile('/export/home/tfoulquier/VEC/sroll2_debug_857ghz_857-1_offset_OFF',dtype = 'float')
off1 = np.fromfile('/export/home/tfoulquier/VEC/sroll2_debug_857ghz_857-2_offset_OFF',dtype = 'float')
off2 = np.fromfile('/export/home/tfoulquier/VEC/sroll2_debug_857ghz_857-3_offset_OFF',dtype = 'float')
off3 = np.fromfile('/export/home/tfoulquier/VEC/sroll2_debug_857ghz_857-4_offset_OFF',dtype = 'float')

plt.plot(off[off>-10000])
plt.plot(off1[off1>-10000])
plt.plot(off2[off2>-10000])
plt.plot(off3[off3>-10000])

im = hp.read_map("/export/home/tfoulquier/MAP/sroll2_debug_857ghz_857ghz_full_0.fits")
#im[im<-1e10]=hp.UNSEEN
im = hp.ud_grade(im,64) # nside 64

print(hp.UNSEEN)
hp.mollview(im,cmap="jet",min = -0.0003,max = 0.0003,title ='sroll result 0') # projection 

im = hp.read_map("/export/home/tfoulquier/MAP/sroll2_debug_857ghz_857ghz_full_1.fits")
#im[im<-1e10]=hp.UNSEEN
im = hp.ud_grade(im,64) # nside 64

hp.mollview(im,cmap="jet",min = -0.0003,max = 0.0003,title ='sroll result 1') # projection 


map_fake = hp.read_map("/export/home/jmdeloui/reduced_FAKE_U")
map_fake = hp.ud_grade(map_fake,64) # nside 64
hp.mollview(map_fake,cmap="jet",min = -0.0003,max = 0.0003,title ='reduced_FAKE_U') # projection 


map_fake = hp.read_map("/export/home/jmdeloui/reduced_FAKE_Q")
map_fake = hp.ud_grade(map_fake,64) # nside 64
hp.mollview(map_fake,cmap="jet",min = -0.0003,max = 0.0003,title ='reduced_FAKE_Q') # projection 


plt.show()


