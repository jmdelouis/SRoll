#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import healpy
import stools
import IMO_4_27


def get_map_units( mapname):
  mapname = mapname.split("/")[-1]
  mapds = None
  for detset in IMO_4_27.DETSETS:
    if detset in mapname:
      mapds = detset
      break
  if mapds == None:
    return "UNKNOWN"
  elif "545" in mapds or "857" in mapds:
    return "MJy/sr"
  else:
    return "KCMB"

  
def merge_fits( filename, suf_in, suf_out):
  # return True if filename has been processed, successfully or not
  suffix = "_%s.fits" % suf_in[0]
  if not filename.endswith( suffix):
    return False
  file_prefix = filename[:-len(suffix)]
  if not os.path.exists( file_prefix + "_%s.fits" % suf_in[1]):
    return False
  fits_out = file_prefix + "_%s.fits" % suf_out
  print "[%d/%d:%s] writing %s: " % (stools.MPI_RANK, stools.MPI_SIZE, stools.HOSTNAME, fits_out)
  maps =list()
  try:
    for sfx in suf_in:
      maps.append( healpy.read_map( file_prefix + "_%s.fits" % sfx))
    healpy.write_map( fits_out, 
                      maps, 
                      column_names=suf_in)
#                      column_units=get_map_units( file_prefix),
#                      overwrite=True)
    for sfx in suf_in:
      os.remove( file_prefix + "_%s.fits" % sfx)
    print "done\n"
  except Exception, e:
    print "failed: %s\n" % str( e)
  return True


################################################################################
################################################################################

if (__name__ == '__main__'):

  if len( sys.argv) < 2:
    print "merge I, Q and U fits maps produced by Sroll in one IQU.fits file"
    print "usage: merge_fits_maps.py <fitsfilename(s)>"

#  print "[%d/%d:%s] starting " % (stools.MPI_RANK, stools.MPI_SIZE, stools.HOSTNAME)
  
  file_list = sys.argv[1:]
  for objname in file_list[stools.MPI_RANK::stools.MPI_SIZE]:
    if os.path.exists( objname):
      if not merge_fits( objname, ("I", "Q", "U"), "IQU"):
        if not merge_fits( objname, ("I", "C", "II"), "ICII"):
          merge_fits( objname, ("II", "IQ", "IU", "QQ", "QU", "UU"), "IIQQUU")
  
