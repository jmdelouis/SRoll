##############################################################################
# example_no_dmc_python.py
#
# This python script aims to demonstrate a simple example for using the python
# wrapper of NO_DMC_LIB.
# It demonstrate the most usual function for a specified object:
# - how to access the metadata
# - how to retrieve data from object
# - how to retrieve flag (written) information from object
#
# Author  : Christian Madsen
# Date    : 2015-03-06
# version : Initial
##############################################################################

import nodmclib


#-----------------------------------------------------------------------------

def printMetadata(obj):
  # Retrieve meta using NO DMC LIB (python wrapper)
  meta = nodmclib.getObjectInfo(obj)

  print("metadata: "), meta

#-----------------------------------------------------------------------------

def readDataAsPIODOUBLE(obj, offset, nbSample):
  # Retrieve data
  data = nodmclib.read_PIODOUBLE(obj, offset, nbSample)

  print("data: "), data

#-----------------------------------------------------------------------------

def readFlagWritten(obj, offset, nbSample):
  # Retrieve Flag
  flags = nodmclib.readFlag_Written(obj, offset, nbSample)

  print("flags: "), flags



##############################################################################
# MAIN
# This code allow to execute the module as a standalone script.
##############################################################################
if __name__ == "__main__":

  obj = "/data/dmc/MISS03/DATA/RingInfo/SpinPhase"
  offset = 10
  nbSample = 12
  

  print("")
  print("****************************************************************")
  print("This is the example program for the python wrapper of NO_DMC_LIB")
  print("****************************************************************")
  print("")

  print("Data used for this example:")
  print("  >> DMC Object fullpath = %s" % obj)
  print("  >> offset = %d" % offset)
  print("  >> number of samples to be read = %d" % nbSample)

  print("")

  printMetadata(obj)
  readDataAsPIODOUBLE(obj, offset, nbSample)
  readFlagWritten(obj, offset, nbSample)

  print("")
  print("ALL SUCCESS !")
