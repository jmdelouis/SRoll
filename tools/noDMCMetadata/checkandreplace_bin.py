#!/usr/bin/env python

import os
import sys
import glob
import nodmclib
import stools
import numpy

BUFLEN = int( 1e8)+100

mpi_size = int(os.getenv('OMPI_COMM_WORLD_SIZE', 1))
mpi_rank = int(os.getenv('OMPI_COMM_WORLD_RANK', 0))

def usage_and_abort():
  print "%s: check whether a bin object is identical to the nodmc one" % sys.argv[0]
  print "usage: %s [--replace] <object(s)>" % sys.argv[0]
  print "       --replace: if set, the new binary object will have the same name as the DMC object"
  print "                  and the DMC object will be deleted"
  exit(1)
  
if len(sys.argv) == 1:
  usage_and_abort()

params = sys.argv[1:]

optReplace = False
for i, objname in enumerate( params):
  if objname.strip().lower() == "--replace":
    del params[i]
    optReplace = True
  elif objname.startswith("-"):
    print "error: unkown option " + objname
    print
    usage_and_abort()

for param in params:
  for objname in glob.glob( param):
    if not os.path.exists( objname + "/no_dmc_metadata.txt"):
      print "not a DMC object: " + objname
      continue
    print "checking " + objname
    PIOType, DataType, BegIdx, EndIdx, BegRing, EndRing, Author, Date = nodmclib.getObjectInfo( objname)
    binname  = objname + ".bin"
    metaname = objname + ".meta"
    if not os.path.exists( metaname):
      print "not converted, skipped"
      continue
    identical = True
    for idx in range( BegIdx, EndIdx, BUFLEN):
      readlen = BUFLEN
      if idx + readlen - 1 > EndIdx:
        readlen = EndIdx - idx + 1
      if DataType == "PIOFLOAT":
        buf = nodmclib.read_PIOFLOAT( objname, idx, readlen)
        bintype = "float32"
      elif DataType == "PIODOUBLE":
        buf = nodmclib.read_PIODOUBLE( objname, idx, readlen)
        bintype = "float64"
      elif DataType == "PIOINT":
        buf = nodmclib.read_PIOINT( objname, idx, readlen)
        bintype = "int32"
      elif DataType == "PIOLONG":
        buf = nodmclib.read_PIOLONG( objname, idx, readlen)
        bintype = "int64"

      binbuf = stools.readBIN( binname, bintype, begin=idx, end=idx+readlen-1)
      if numpy.all( buf == binbuf):
        print "o",
      else:
        if (numpy.all(numpy.isnan( buf) == numpy.isnan( binbuf))) and (numpy.all(numpy.isfinite( buf) == numpy.isfinite( binbuf))):
          print "o",
        else:
          identical = False
          print "X",
      sys.stdout.flush()

    if identical:
      print( "OK")
      if optReplace:
        os.system("./rm_nodmc_object.py " + objname)
        os.system("mv %s %s" % (binname, objname))
        print( "replaced\n")
    else:
      print( "DIFFER\n")
