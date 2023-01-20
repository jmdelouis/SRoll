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
"""
ph_gal=184.557
th_gal=-5
"""
im = hp.read_map(nfile)
hp.mollview(im,cmap='jet',min = valmin,max=valmax,title =name )
hp.gnomview(im,cmap='jet',rot=[ph_gal,th_gal],min = valmin,max = valmax,title =name) 
plt.show()
