#!/usr/bin/env python

# 21 March 2020
# authors:
# - Luca Pagano
# - Sylvain Mottet

import os
import sys
import numpy as np
import healpy as hp
import time


# command line management
try:
  assert len( sys.argv) == 4
  assert sys.argv[1] in ["32", "64"]
  assert sys.argv[2] in ["100", "143", "217", "353"]
  assert sys.argv[3] in ["0", "1", "2", "3", "4"]
  first_sim = int( sys.argv[3]) * 100
except:
  print "usage: %s <nsideout> <freq> <hundred>" % sys.argv[0]
  print "  where:"
  print "  - <nsideout> is 32 or 64,"
  print "  - <freq> is 100, 143, 217, or 353 "
  print "  - <hundred> is 0, 1, 2, 3, or 4"
  print "downgrade JAN18r60 sims for channel <freq>, from iter <hundred>*100 to <hundred>*100+99"
  print
  exit(-1)


OVERWRITE = False # if True, reprocess and overwrite FITS file if it already exists

# Effective Beam FWHM in arcmin
# LFI from L02: https://arxiv.org/pdf/1807.06206.pdf, table 1 p.4
# HFI from L03: https://arxiv.org/pdf/1807.06207, table 12 p.47
FWHM = { "33": 32.29,  "44": 26.99,  "70": 13.22, 
        "100":  9.68, "143":  7.30, "217":  5.02,
        "353":  4.94, "545":  4.83, "857":  4.64}

nsideout = int( sys.argv[1])
nsidein  = 2048
npixin   = hp.nside2npix( nsidein)

mapin_pfix =  '/scratch/cnt0028/ias1717/SHARED/tmp_sim_dgrade/JAN18r60_'
mapin_pfix =  '/scratch/cnt0028/ias1717/SHARED/bware/JAN18r60_MAP/JAN18r60_'
if nsideout == 32:
  mapout_pfix = '/scratch/cnt0028/ias1717/SHARED/bware/sroll20_sim_nside32/sroll20_sim_ns32_'
else:
  mapout_pfix = '/scratch/cnt0028/ias1717/SHARED/bware/sroll20_sim_nside64/sroll20_sim_ns64_'

channels = ['100ghz', '143psb', '217psb', '353psb']
if sys.argv[2] == "100":
  channel = "100ghz"
else:
  channel = sys.argv[2] + "psb"
assert channel in channels
channels = [channel]

# takes about 40min for 100 sim files on one occigen node
cuts = [ "ghz_full", "ghz_fulleven", "ghz_fullodd", "ghz_hm1", "ghz_hm2",
         "ds1_full", "ds1_fulleven", "ds1_fullodd", "ds1_hm1", "ds1_hm2",
         "ds2_full", "ds2_fulleven", "ds2_fullodd", "ds2_hm1", "ds2_hm2",]
if nsideout == 64:
  cuts = [ "ghz_full", "ds1_full", "ds2_full"]

wpixin  = hp.pixwin( 2048,     pol=True) #, lmax=nsideout*4)
wpixout = hp.pixwin( nsideout, pol=True) #, lmax=nsideout*4)

l = np.arange( np.size( wpixout[0]))
if nsideout <= 32:
  print "smooth map with ns=%d cosinus window function" % nsideout
  beamout = 0.5 * (1.0 + np.cos( np.pi * (l-nsideout) / 2.0 / nsideout))
  beamout[0:nsideout+1] = 1.0
  beamout[3*nsideout+1:] = 0.0
  beamout = np.array( [beamout, beamout, beamout])
else:
  # smooth the map with a gaussian beam of FWHM the size of a nsideout pixel
  fwhm_out = hp.max_pixrad( nsideout, degrees=True) * 2 * 60
  print "smooth map with %.1f arcmin FWHM gaussian beam" % fwhm_out
  beamout = hp.gauss_beam( np.deg2rad( fwhm_out / 60.0), lmax=nsideout*4, pol=True)

#print "wpixin.shape=", len(wpixin), len(wpixin[0])
#print "wpixout.shape=", len(wpixout), len(wpixout[0])
#print "beamout.shape=", beamout.shape

for channel in channels:
  beamin = hp.gauss_beam( np.deg2rad( FWHM[channel[0:3]]/60.0), lmax=nsideout*4, pol=True)
  for cut in cuts:
    cut = channel[0:3] + cut
    if channel.endswith( "psb"):
      cut = cut.replace( "ghz", "psb")
    for isim in np.arange( first_sim, first_sim+100):
      strsim = str(isim).zfill(3)
      mapin = mapin_pfix + strsim + '_' + channel + '_' + cut + '_IQU.fits'
      mapout = mapout_pfix + strsim + '_' + cut + '.fits'
      if os.path.exists( mapout) and not OVERWRITE:
        print mapout + " exists, skipping"
        continue
      if channel == "143psb":
        mapin = mapin.replace( "JAN18r60", "JAN18r79")
      print time.asctime() + ": reading " + mapin
      I,Q,U = hp.read_map( mapin, [0,1,2], verbose=False)
      assert len( I) == npixin
      med = np.median( I[I != hp.UNSEEN])
      I[I==hp.UNSEEN] = med
      med = np.median( Q[Q != hp.UNSEEN])
      Q[Q==hp.UNSEEN] = med
      med = np.median( U[U != hp.UNSEEN])
      U[U==hp.UNSEEN] = med
      alm = hp.map2alm([I,Q,U], lmax=nsideout*4)
      hp.almxfl( alm[0], beamout[:,0] / beamin[:,0] / wpixin[0][0:nsideout*4+1] * wpixout[0], inplace=True)
      hp.almxfl( alm[1], beamout[:,1] / beamin[:,1] / wpixin[1][0:nsideout*4+1] * wpixout[1], inplace=True)
      hp.almxfl( alm[2], beamout[:,2] / beamin[:,2] / wpixin[1][0:nsideout*4+1] * wpixout[1], inplace=True)
      m = hp.alm2map( alm, nsideout, lmax=nsideout*4, pixwin=False, pol=True)
      m = hp.reorder( m, r2n=True)
      print time.asctime() + ": writing " + mapout
      hp.write_map( mapout, m, nest=True, overwrite=True)
      
#      exit(0)
