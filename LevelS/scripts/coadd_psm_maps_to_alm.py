#!/usr/bin/env python

from __future__ import division, print_function

import os
import sys
import healpy as hp
import numpy as np
import stools
import time
from IMO_4_27 import * # BOLOID, DETSETS, CALIB, NEP, XPOL, ELECWNOISE, PHOTWNOISE


# on M4, due to the nside 4096 "*ps" maps, this script produces a memory error on login nodes
# must be run in "qsub (-I) -l nodes=1:ppn=24,walltime=24:00:00"


def map_to_alm( map_in, alm_out=None, map_out=None, nsideout=2048, overwrite=False):
  if (alm_out is not None) and (os.path.exists( alm_out)) and (not overwrite):
    print( "reading alm", alm_out)
    try:
      alm = hp.read_alm( alm_out, hdu=(1,2,3))
    except:
      alm = hp.read_alm( alm_out)
    if len( alm.shape) == 2:
      print( "polarised")
    else:
      print( "not polarised")
    return alm
  print( "reading map", map_in)
  m = hp.read_map( map_in, field=None, verbose=False)
  if len( m.shape) == 2:
    ns = hp.npix2nside( len( m[0]))
    print( "polarised, nside=%d" % ns)
  else:
    ns = hp.npix2nside( len( m))
    print( "not polarised, nside=%d" % ns)
  print( "converting to alm")
  alm = hp.map2alm( m, lmax=3*nsideout-1)
  if map_out is not None:
    print( "writing map", map_out)
    m2 = hp.alm2map( alm, nside=nsideout)
    hp.write_map( map_out, m2)
  if alm_out is not None:
    print( "writing alm", alm_out)
    hp.write_alm( alm_out, alm)
  return alm


################################################################################

NSIDEOUT = 2048

COMPLIST = [
  # from JMD /bluevan1/symottet/FFP10_apr17_PSM_maps/psmmap2toi/make_co2020_maps.py
  # nside=2048, nopol, beam=0 (9.68arcmin, actually), units=KCMB, dtype=float32
  "co12_2020", "co13_2020",                 
  # from /project/projectdirs/planck/data/ffp10/sky/PSM/PSM_OUTPUT_FFP10_GALAXY_NEWDUST
  "freefree", "synchrotron", "thermaldust", 
  # from /project/projectdirs/planck/data/ffp10/sky/PSM/FFP10_EXTRA_GALACTIC
  "firb",                                   
  # from /project/projectdirs/planck/data/ffp10/sky/PSM/FFP10_SZ
  "kineticsz", "thermalsz",                 
  # from /project/projectdirs/planck/data/ffp10/sky/PSM/FFP10_PS
  # !!! nside=4096 !!!
  "faintirps", "faintradiops",              
]

COMP_TPL = "/bluevan1/symottet/FFP10_apr17_PSM_maps/{comp}_map_detector_{psmdet}.fits"
CMB_TPL  = "/bluevan1/symottet/FFP10_apr17_PSM_maps/ffp10_lensed_scl_cmbfid_nodip_{freq}_map.fits"
ALM_TPL  = "/bluevan1/symottet/APR20_sky_alm/APR20_{comp}_alm_{pixname}.fits"
MOUT_TPL = "/bluevan1/symottet/APR20_sky_alm/APR20_{comp}_map_{pixname}.fits"

assert len( sys.argv) == 2, "usage: %s <freq>" % sys.argv[0]

freq = sys.argv[1]
assert freq in ("100", "143", "217", "353", "545", "857")

detset = freq+"ghz"

print( "****", time.ctime())

comp = "cmbfid"
cmb_alm = map_to_alm( CMB_TPL.format( **locals()),
                      overwrite = True)
#                      alm_out = ALM_TPL.format( **locals()),
#                      map_out = MOUT_TPL.format( **locals()),
#                      overwrite = False)
print( "cmb_alm.shape=", cmb_alm.shape)

for pixname in DETSETS[detset]:
  print( "****", time.ctime())
  psmdet = pixname.replace("-", "_")
  sky_alm = cmb_alm.copy()
  for comp in COMPLIST:
    comp_alm = map_to_alm( COMP_TPL.format( **locals()),
                           overwrite = True)
#                           alm_out = ALM_TPL.format( **locals()),
#                           map_out = MOUT_TPL.format( **locals()),
#                           overwrite = False)

    if len( comp_alm.shape) == 1:
      # unpolarized component only add T alm
      sky_alm[0] += comp_alm
    else:
      # add T, E, B alm
      sky_alm += comp_alm

  comp = "sky"
  almname = ALM_TPL.format( **locals())
  print( "writing", almname)
  hp.write_alm( almname, sky_alm)

  map_out = MOUT_TPL.format( **locals())
  m = hp.alm2map( sky_alm, nside=NSIDEOUT)
  hp.write_map( map_out, m)
  print( "")

print( "****", time.ctime(), "\n")
