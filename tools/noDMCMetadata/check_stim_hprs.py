#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, sys
import numpy
import stools
from IMO_4_27 import * # BOLOID, DETSETS, CALIB, NEP, XPOL, ELECWNOISE, PHOTWNOISE

BEGRING = 240
ENDRING = 26050
HPRSIZE = 27664

try:
  SROLLHOST  = os.environ["SROLLHOST"]
  SROLLDIR   = os.environ["SROLLDIR"]
  NODMCDATA  = os.environ["NODMCDATA"]
  SCRATCHDIR = os.environ["SCRATCHDIR"]
except:
  raise Exception( "You must 'source ./srollex_setenv.sh' before using this script")

if len( sys.argv) != 5:
  print "usage: %s <hprtemplate> <pixname/detset> <firstiter> <lastiter>" % sys.argv[0]
  print "where:"
  print "  - <hprtemplate> is the HPR name python template"
  print "                  eg. {NODMCDATA}/JAN18_stimHPR/{pixname}_JAN18_stimHPR_{iternum:03d}.float32.bin"
  print "                  or  {SCRATCHDIR}/stimr56_VEC/{pixname}_JAN18_stimHPR_{iternum:03d}.float32.bin"
  print "  - <pixname/detset> is either a single pixname or a detset"
  print "                  eg. 100-1a, 143psb, 353ghz, lofreq (100GHz-353GHz), JAN18r60 (100ghz+143psb+353psb)"
  print "  - <firstiter> is the first iteration number to check (eg. 0)"
  print "  - <lastiter> is the last iteration number to check (eg. 99)"
  exit(0)

hprtpl    = sys.argv[1]
detset    = sys.argv[2]
firstiter = int( sys.argv[3])
lastiter  = int( sys.argv[4])
verbose = True

if detset in DETSETS.keys():
  detlist = DETSETS[detset]
elif detset in BOLOID.keys():
  detlist = [detset]
elif detset == "JAN18r60":
  detlist = DETSETS["100ghz"] + DETSETS["143psb"] + DETSETS["353psb"]
else:
  raise Exception("Unkown bolometer/detset: " + detset)

print detlist

for pixname in detlist:
  boloid = BOLOID[pixname]
  badrings = numpy.fromfile( NODMCDATA + "/calROIs/%s_discarded_rings_dx11" % boloid, dtype="int32")
  iternum = firstiter
  hprname = hprtpl.format( **locals())
  print hprname
  refhpr = numpy.fromfile( hprname, dtype="float32")
  refmean = numpy.zeros( ENDRING+1, dtype="float32")
  for ring in range( BEGRING, ENDRING+1):
    refmean[ring] = numpy.mean( refhpr[ring*HPRSIZE:(ring+1)*HPRSIZE])
  refmean = numpy.abs( refmean)
  for iternum in range( firstiter+1, lastiter+1):
    hprname = hprtpl.format( **locals())
    print hprname
    if not os.path.exists( hprname):
      print "!!! file not found !!!"
      continue
    checkhpr = numpy.fromfile( hprname, dtype="float32")
    badhitrings = list()
    badmeanrings = list()
    badmeanvals = list()
    badspl = 0
    badmean = 0
    for ring in range( BEGRING, ENDRING+1):
      if badrings[ring] == 0:
        ndiffhits = numpy.sum( (refhpr[ring*HPRSIZE:(ring+1)*HPRSIZE] == 0) != (checkhpr[ring*HPRSIZE:(ring+1)*HPRSIZE] == 0))
        ringmean = numpy.abs( numpy.mean( checkhpr[ring*HPRSIZE:(ring+1)*HPRSIZE]))
        if ndiffhits != 0:
          badhitrings.append( ring)
          badspl += ndiffhits
        if (ringmean < refmean[ring]*0.8) or (ringmean > refmean[ring]*1.2):
          badmean += 1
          badmeanrings.append( ring)
          badmeanvals.append( ringmean*100.0/refmean[ring])
    if verbose:
      if badspl != 0:
        print "%d differing hits in %d rings" % (badspl, len( badhitrings)), badhitrings[:10]
      for bm in range( min( badmean, 10)):
        print "discrepant mean in ring %d (out of %d): %d%%" % (badmeanrings[bm], badmean, badmeanvals[bm])

