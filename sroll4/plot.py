import healpy as hp
import matplotlib.pyplot as plt
import sys

f = sys.argv[1]

min_im = sys.argv[2]
max_im = sys.argv[3]

im = hp.read_map(f)
hp.mollview(im,min = min_im ,max =max_im)
plt.show()

