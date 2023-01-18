#!/usr/bin/env python

import os
import sys
import glob
import nodmclib
import array


################################################################################

# this function comes from stools.py, but importing stools fails on CORI where
# there is not yet any healpy installation, so I'm just copying it here to remove
# the stools dependency

def pwrite( filename, offset, data, datatype="float32"):
# python 3.3: os.pwrite( filedesc, bytestring, offset)
# file must exist, you can non-destructively create it with "open( filename, 'a').close()"
# TODO: consider locking the file with fcntl.lockf()
  if (datatype == "float32") or (datatype == "f32") or (datatype == "PIOFLOAT"):
    array_datatype = "f"
    datasize = 4
  elif (datatype == "float64") or (datatype == "f64") or (datatype == "PIODOUBLE"):
    datasize = 8
    array_datatype = "d"
  elif (datatype == "int32") or (datatype == "i32") or (datatype == "PIOINT"):
    datasize = 4
    array_datatype = "i"
  elif (datatype == "int64") or (datatype == "i64") or (datatype == "PIOLONG"):
    datasize = 8
    array_datatype = "l"
  outfile = open( filename, 'r+b')
  outfile.seek( offset * datasize)
  out_array = array.array( array_datatype, data)
  assert out_array.itemsize == datasize # itemsize may vary depending on the platform
  out_array.tofile( outfile)
  outfile.close()


################################################################################

def usage_and_abort():
  print "%s: convert a dmc object to a binary file using nodmclib.py" % sys.argv[0]
  print "    takes about 1 hour for a PIOFLOAT TOI on occigen"
  print "usage: %s <object(s)>" % sys.argv[0]
  print "       <object(s)>: full object name(s), accepts wildcards"
  print "use 'checkandreplace_bin.py' to replace the nodmc object by its binary version"
  exit(1)


################################################################################

try:
  SROLLHOST = os.environ["SROLLHOST"]
  DMCDATA   = os.environ["DMCDATA"]
except:
  raise Exception( "You must 'source ./srollex_setenv.sh' before using this script")

if SROLLHOST == "M3":
  BINOUTDIR = "/redtruck/SimuData"
elif SROLLHOST == "CC":
  BINOUTDIR = "/sps/planck/SimuData/nodmc"
elif SROLLHOST == "CORI":
  BINOUTDIR = "/global/cscratch1/sd/smottet"
else:
  # there should not be objects left in DMC format on other clusters
  raise Exception( "Unsupported SROLLHOST (%s)" % SROLLHOST)

BUFLEN = int( 1e8)
#mpi_size = int(os.getenv('OMPI_COMM_WORLD_SIZE', 1))
#mpi_rank = int(os.getenv('OMPI_COMM_WORLD_RANK', 0))

if len( sys.argv) == 1:
  usage_and_abort()

params = sys.argv[1:]

for param in params:
  for objname in glob.glob( param):
    if not os.path.exists( objname + "/no_dmc_metadata.txt"):
      print objname + ": no_dmc_metadata.txt not found, skipping\n"
      continue
    if not objname.startswith( DMCDATA):
      print objname + ": not in DMCDATA path, skipping\n"
      continue
    outpfix = objname.replace( DMCDATA, BINOUTDIR)
    metaname = outpfix + ".meta"
    if os.path.exists( metaname):
      print objname + ": already converted, skipping\n"
      continue

    print "converting " + objname
    PIOType, DataType, BegIdx, EndIdx, BegRing, EndRing, Author, Date = nodmclib.getObjectInfo( objname)
    if not (DataType in ["PIOINT", "PIOLONG", "PIOFLOAT", "PIODOUBLE"]):
      print objname + ": %s data type not implemented, skipping\n" % DataType
    if DataType == "PIOFLOAT":
      typtag = ".f32"
    elif DataType == "PIODOUBLE":
      typtag = ".f64"
    elif DataType == "PIOINT":
      typtag = ".i32"
    elif DataType == "PIOLONG":
      typtag = ".i64"

    binname = outpfix + typtag + ".bin"
    os.system( "mkdir -p " + os.path.dirname( binname))
    if Author == 68:
      Author = "smottet"
    print "PIOType:", PIOType
    print "DataType:", DataType
    print "BegIdx:", BegIdx
    print "EndIdx:", EndIdx
    print "BegRing:", BegRing
    print "EndRing:", EndRing
    print "Author:", Author
    print "Date:", Date
    print "deleting/creating " + binname
    os.system("rm -f {binname}; touch {binname}".format( **locals()))
    for idx in range( BegIdx, EndIdx, BUFLEN):
      readlen = BUFLEN
      if idx + readlen - 1 > EndIdx:
        readlen = EndIdx - idx + 1
      print "r",
      sys.stdout.flush()
      if DataType == "PIOFLOAT":
        buf = nodmclib.read_PIOFLOAT( objname, idx, readlen)
      elif DataType == "PIODOUBLE":
        buf = nodmclib.read_PIODOUBLE( objname, idx, readlen)
      elif DataType == "PIOINT":
        buf = nodmclib.read_PIOINT( objname, idx, readlen)
      elif DataType == "PIOLONG":
        buf = nodmclib.read_PIOLONG( objname, idx, readlen)
      print "w",
      sys.stdout.flush()
      pwrite( binname, idx, buf, DataType)
    os.system("cp {objname}/no_dmc_metadata.txt {metaname}".format( **locals()))
    print("\n")
