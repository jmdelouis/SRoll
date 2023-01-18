#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import stools
import sys


def usage():
  print """
usage: %s <detset> <hundred>
where <detset> in [100ghz, 143psb, 217psb, 353psb]
and <hundred> in [0, 1, 2, 3, 4]
""" % sys.argv[0]
  exit(0)
  

# deltapsi simulations mask
MASK = stools.NODMCDATA + "/MASK_2048_GALACTIC_fits/HFI_Mask_GalPlane-apo2_fsky52_2048_R2.00.fits"
# JMD SRoll2 70/64% fsky mask with PS apodized
MASK = stools.NODMCDATA + "/MASK_2048_GALACTIC_fits/MASK_S2PAPER_0_64_PS.fits"
# Julia DPC paper 43% fsky mask with PS ont apodized
MASK = stools.NODMCDATA + "/MASK_2048_GALACTIC_fits/mask_fsky043_PS_100-353.fits"

if MASK.endswith( "MASK_S2PAPER_0_64_PS.fits"):
  FSKY = 64
elif MASK.endswith( "mask_fsky043_PS_100-353.fits"):
  FSKY = 43
else:
  raise Exception("unhandled mask (%s)" % MASK)

SCRIPTDIR = os.path.dirname( os.path.realpath( __file__))

SPICEMAX = 3000

DETSETS = ["100ghz", "143psb", "217psb", "353psb"]

if len( sys.argv) != 3:
  usage()

detset = sys.argv[1]
if not detset in DETSETS:
  usage()

hundred = int( sys.argv[2])
if not hundred in range(5):
  usage()

if detset == "143psb":
  FITSTPL  = "/scratch/cnt0028/ias1717/SHARED/bware/JAN18r79_MAP/JAN18r79_{niter:03d}_{detset}_{detset}_{ringcut}_IQU.fits"
else:
  FITSTPL  = "/scratch/cnt0028/ias1717/SHARED/bware/JAN18r60_MAP/JAN18r60_{niter:03d}_{detset}_{detset}_{ringcut}_IQU.fits"

for niter in range( hundred*100, (hundred+1)*100):
#for niter in range( 1):
  fits1 = FITSTPL.format( ringcut = "hm1", **locals())
  fits2 = FITSTPL.format( ringcut = "hm2", **locals())
  clname = SCRIPTDIR + "/cl_JAN18_fsky{FSKY}/cl_JAN18_fsky{FSKY}_hmdiff_{detset}_{niter:03d}.txt".format( **locals())
  if os.path.exists( fits1):
    print fits1
    if not os.path.exists( clname):
      stools.spice_diff( fits1, fits2, weight=MASK, lmax = SPICEMAX, cl_filename = clname)
