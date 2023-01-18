#!/usr/bin/env python

from __future__ import division, print_function
import os
import time

import numpy as np
import healpy as hp

# les donnees MAPxxx sont sur renater file sender:
# https://filesender.renater.fr/?s=download&token=438ddabd-393d-4724-b01b-b179471cbd61

# cosed=np.load('coefbolo.npy').tolist()

cosed = {
  u'100-1a': [10.87, 16.96],
  u'100-1b': [12.61, 16.4 ],
  u'100-2a': [14.69, 14.08],
  u'100-2b': [12.01, 17.5 ],
  u'100-3a': [16.36, 14.52],
  u'100-3b': [11.78, 13.78],
  u'100-4a': [19.09, 18.64],
  u'100-4b': [16.11, 17.57],
  u'143-1a': [0.0613, 0.0022],
  u'143-1b': [0.0437, 0.0017],
  u'143-2a': [0.0523, 0.002 ],
  u'143-2b': [0.0557, 0.0022],
  u'143-3a': [0.0881, 0.003 ],
  u'143-3b': [0.0737, 0.0023],
  u'143-4a': [0.0489, 0.0018],
  u'143-4b': [0.0493, 0.0019],
  u'143-5':  [0.021 , 0.0012],
  u'143-6':  [0.0579, 0.002 ],
  u'143-7':  [0.0099, 0.0005],
  u'143-8':  [0.0404, 0.0018],
  u'217-1':  [50.22, 34.42],
  u'217-2':  [42.47, 32.73],
  u'217-3':  [51.23, 37.37],
  u'217-4':  [47.75, 30.87],
  u'217-5a': [43.97, 35.85],
  u'217-5b': [43.68, 38.54],
  u'217-6a': [38.92, 41.21],
  u'217-6b': [40.75, 33.33],
  u'217-7a': [45.5 , 41.57],
  u'217-7b': [43.58, 33.19],
  u'217-8a': [45.3 , 41.48],
  u'217-8b': [41.78, 34.16],
  u'353-1':  [170.3,  82.5],
  u'353-2':  [174. , 130.8],
  u'353-3a': [185.4, 133.3],
  u'353-3b': [200.7, 166.6],
  u'353-4a': [172.9, 121. ],
  u'353-4b': [140.9, 125.2],
  u'353-5a': [150.3, 138.1],
  u'353-5b': [159.8, 143.9],
  u'353-6a': [148.9, 143. ],
  u'353-6b': [166.4, 167.1],
  u'353-7':  [196.9, 110.9],
  u'353-8':  [185.3,  99.9],
  u'545-1':  [141143.4258,  68060.0609],
  u'545-2':  [150154.9882,  88693.9566],
  u'545-3':  [150919.4925, 110276.8587],
  u'545-4':  [157008.9537,  88696.7343],
  u'857-1':  [554873.58501, 664832.95365],
  u'857-2':  [517745.71884, 618752.26455],
  u'857-3':  [555657.15421, 634683.25204],
  u'857-4':  [529623.60368, 577309.67588],
}

im12co = hp.read_map('/bluevan1/symottet/FFP10_apr17_PSM_maps/psmmap2toi/MAP_MODEL_12CO.fits')
im13co = hp.read_map('/bluevan1/symottet/FFP10_apr17_PSM_maps/psmmap2toi/MAP_MODEL_13CO.fits')

# clean model maps for unwanted artefacts and negative values
im12co[im12co < 1e-10] = 0.0
im13co[im13co < 1e-10] = 0.0

hdr = list()
hdr += [("BEAMSIZE", 0.0, "Beam size in arcminutes")] # CO templates are in fact convolved with a 100ghz beam (9.68 arcmin), but we'll neglect it for compatibility with the old PSM maps
hdr += [("CREATOR",  os.path.realpath( __file__), "File created by")]
hdr += [("CREADATE", time.ctime(), "Creation date")]
print( hdr)

for pixname in cosed.keys():
  print( "%s\t12CO: %10g, 13CO: %10g" % (pixname, cosed[pixname][0], cosed[pixname][1]))
  bolo12co = im12co * cosed[pixname][0]
  bolo13co = im13co * cosed[pixname][1]
  psmdet = pixname.replace("-", "_")
  mapname = "/bluevan1/symottet/FFP10_apr17_PSM_maps/co12_2020_map_detector_%s.fits" % psmdet
  print( "writing", mapname)
  hp.write_map( mapname, bolo12co, column_units="K_CMB", extra_header=hdr, overwrite=True) 
  mapname = "/bluevan1/symottet/FFP10_apr17_PSM_maps/co13_2020_map_detector_%s.fits" % psmdet
  print( "writing", mapname)
  hp.write_map( mapname, bolo13co, column_units="K_CMB", extra_header=hdr, overwrite=True)
