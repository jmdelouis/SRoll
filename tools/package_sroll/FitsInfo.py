#!/usr/bin/env python

import sys
try:
  import astropy.io.fits as pyfits
except:
  import pyfits


if len(sys.argv) != 2:
  print( "display information about a FITS file")
  print( "usage: %s <fitsfilename>" % sys.argv[0])
  exit(1)

hdulist = pyfits.open( sys.argv[1], mode="readonly", memmap=True)
print( "========== file header ==========")
hdulist.info()
for i, hdu in enumerate( hdulist):
  print( "\n========== HDU[%d] header ==========" % i)
  for keyword in hdu.header.keys():
#    if keyword == "COMMENT": continue
    print( "%15s = %15s / %s" % (keyword, hdu.header[keyword], hdu.header.comments[keyword]))
  if (i == 1) and (len( hdulist) > 2):
    goon = input("Continue? [Y/n]")
    if goon != "" and goon[0].lower() == "n":
      break
hdulist.close()

