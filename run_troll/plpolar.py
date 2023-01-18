import numpy as np
import matplotlib.pyplot as plt
import healpy as hp

map=np.load('POLARMAP.npy').reshape(12*32*32,2)

hp.mollview(map[:,0],cmap='jet',hold=False,sub=(2,1,1),nest=True,min=-40,max=40)
hp.mollview(map[:,1],cmap='jet',hold=False,sub=(2,1,2),nest=True,min=-40,max=40)
plt.show()
