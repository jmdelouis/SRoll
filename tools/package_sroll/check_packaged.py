#!/usr/bin/env python

################################################################################
#
# check_packaged.py
# compare packaged sroll maps with PLA 2018 when possible
#
################################################################################


import healpy
import numpy


################################################################################

HFIFREQS = ["100", "143", "217", "353", "545", "857"]
PLACUTS  = ["full", "halfmission-1", "halfmission-2", "full-oddring", "full-evenring"]

for freq in HFIFREQS:
  for cut in PLACUTS:
    if freq <= "353":
      COMPLIST = ["I", "Q", "U", "H", "II", "IQ", "IU", "QQ", "QU", "UU"]
    else:
      COMPLIST = ["I", "H", "II"]
    pname = "/scratch/cnt0028/ias1717/SHARED/RD12_data/dmc/MISS03/DATA/PLA2018_HFI/HFI_SkyMap_{freq}_2048_R3.01_{cut}.fits".format(**locals())
    sname = "/scratch/cnt0028/ias1717/smottet/sroll21_datarelease/final_fits/SRoll21_SkyMap_{freq}ghz_{cut}.fits".format(**locals())
    print "\n" + sname + "\n" + pname
    for field in range( len( COMPLIST)):
      pmap = healpy.read_map( pname, field=field, verbose=False)
      smap = healpy.read_map( sname, field=field, verbose=False)
      smap[pmap == healpy.UNSEEN] = healpy.UNSEEN
      smap[pmap == 0] = 0
      pmap[smap == healpy.UNSEEN] = healpy.UNSEEN
      pmap[smap == 0] = 0
      ndiff = numpy.sum( numpy.not_equal( pmap, smap))
      print "  %s diff=%.5f%%" % (COMPLIST[field], ndiff*100.0/len( smap))


