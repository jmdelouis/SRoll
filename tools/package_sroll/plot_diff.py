#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import stools
import sys

import healpy
from matplotlib import pyplot

# deltapsi simulations mask
MASK = stools.NODMCDATA + "/MASK_2048_GALACTIC_fits/HFI_Mask_GalPlane-apo2_fsky52_2048_R2.00.fits"
# JMD SRoll2 70/64% fsky mask with PS apodized
MASK = stools.NODMCDATA + "/MASK_2048_GALACTIC_fits/MASK_S2PAPER_0_64_PS.fits"
# Julia DPC paper 43% fsky mask with PS ont apodized
MASK = stools.NODMCDATA + "/MASK_2048_GALACTIC_fits/mask_fsky043_PS_100-353.fits"

map1 = "/scratch/cnt0028/ias1717/smottet/sylvain_sroll3_MAP/sroll3_s21_512MPI_100ghz_100ghz_full_I.fits"
map2 = "/scratch/cnt0028/ias1717/smottet/sylvain_sroll3_MAP/sroll3_s21_100ghz_100ghz_full_I.fits"
title = "100GHz sroll21 done with sroll3 on occigen - MPI512-MPI1024"
nside = 2048
minmax = 0.001


map2 = "/scratch/cnt0028/ias1717/smottet/sylvain_sroll3_MAP/sroll3_s21_512MPI_100ghz_100ghz_full_I.fits"
map1 = "/scratch/cnt0028/ias1717/SHARED/bware/SRoll20_REP6_M4_512N_Luca/100GHz_RD12_REP6_data_sroll2r79_dipHFI17_full_I.fits"
title = "100GHz sroll20 done with sroll3 on occigen minus Luca's sroll20"
nside = 2048
minmax = 0.001


map1 = "/scratch/cnt0028/ias1717/smottet/sylvain_sroll3_MAP/sroll3_s21_512MPI_143ghz_143ghz_full_I.fits"
map2 = "/scratch/cnt0028/ias1717/smottet/sylvain_sroll3_MAP/sroll3_s21_512MPI_143psb_143psb_full_I.fits"
title = "143GHz minus 143psb"
nside = 2048
minmax = 0.001
minmax = 10



m1 = healpy.read_map( map1)
m2 = healpy.read_map( map2)
mdiff = m1-m2
mdiff[m1 == healpy.UNSEEN] = healpy.UNSEEN
mdiff[m2 == healpy.UNSEEN] = healpy.UNSEEN

if nside < 2048:
  mdiff = healpy.ud_grade( mdiff, nside)

mdiff = healpy.remove_monopole( mdiff)

if minmax != 0:
  healpy.mollview( mdiff*1e6, cmap="jet", unit="microK", min=-minmax, max=minmax, title=title + ", nside=%d"%nside)
else:
  healpy.mollview( mdiff*1e6, cmap="jet", unit="microK", title=title + ", nside=%d"%nside)

pyplot.show()
