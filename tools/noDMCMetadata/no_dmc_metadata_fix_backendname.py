#!/usr/bin/env python

##############################################################################
# no_dmc_metadata_fix_backendname.py
#
# This python script aims to fix 'backendname' entry in "no_dmc_metadata.txt"
# files according to a specified pattern, so that only files containing the
# matching path will be updated.
#
# You may use this script when you move object data (and metadata) to another
# directory so that NO_DMC_LIB can still access to them. 
#
# NOTE: this script does NOT require DMC stuff!
#
# Author  : Christian Madsen
# Date    : 2015-06-02
# version : Initial
##############################################################################

import os


### Parameters
METADATA_FILENAME = "no_dmc_metadata.txt"

NODMCMETADATA_BACKENDNAME = "Backendname"

# Debug and verbose mode
DEBUG = False # Set to True if you want debug message to be print
VERBOSE = True # Set to True if you want more log message


#-----------------------------------------------------------------------------

def logDebug(msg):
  """Print message only if DEBUG is set to True.
  """
  if DEBUG:
    print("[DBG] %s" % msg)

#-----------------------------------------------------------------------------

def log(msg):
  """Print message only if VERBOSE is set to True.
  """
  if VERBOSE:
    print(msg)

#-----------------------------------------------------------------------------

def printUsage(cmd):
  print("Usage: %s <workingDir> <old_pattern> <new_pattern>" % cmd)

#-----------------------------------------------------------------------------

def getAllMetadataFiles(workingDir):
  import os.path

  res_list = []
  for dirpath, dirnames, filenames in os.walk(workingDir):
    for filename in filenames:
      if filename == METADATA_FILENAME:
        res_list.append(os.path.join(dirpath, filename))

  return res_list

#-----------------------------------------------------------------------------

def getSelectedMetadataFiles(files, old_pattern):
  res_selectedFiles = []

  for file in files:
    if old_pattern in open(file).read():
      res_selectedFiles.append(file)

  return res_selectedFiles

#-----------------------------------------------------------------------------

def replacePatternInMetadataFiles(selectedMetadataFiles, old_pattern, new_pattern):
  import fileinput

  # Loop on all files to be proceed
  for file in selectedMetadataFiles:
    log("  . Proceed file %s" % file)

    # Open file for reading
    s = open(file).read()

    # Replace
    s = s.replace(old_pattern, new_pattern)

    # write result
    f = open(file, 'w')
    f.write(s)
    f.flush()
    f.close()

#-----------------------------------------------------------------------------


##############################################################################
# MAIN
# This code allow to execute the module as a standalone script.
##############################################################################
if __name__ == "__main__":
  import sys

  shortName = sys.argv[0].split("/")[-1]

  log("")
  log("*******************************************************************************")
  log("* %s" % shortName)
  log("* This is the NO_DMC_LIB metadata fix backendname script")
  log("*******************************************************************************")
  log("")

  if len(sys.argv) != 4:
    printUsage(shortName)
    exit(1)

  # Retrieve the working dir
  workingDir = sys.argv[1]

  # Retrieve the old_pattern (the one to be replaced)
  old_pattern = sys.argv[2]

  # Retrieve the new_pattern (the one that will replace the old one)
  new_pattern = sys.argv[3]

  log("workingDir  = %s" % workingDir)
  log("old_pattern = '%s'" % old_pattern)
  log("new_pattern = '%s'" % new_pattern)
  log("")

  ### 1) Find all metadata file
  allMetadataFiles = getAllMetadataFiles(workingDir)
  #logDebug("allMetadataFiles = %s" % allMetadataFiles)
  log("> Found a total of %d metadata files" % len(allMetadataFiles))

  ### 2) Select those that match old_pattern for 'Backendname'
  selectedMetadataFiles = getSelectedMetadataFiles(allMetadataFiles, old_pattern)
  log("> %d contains the <old_pattern>" % len(selectedMetadataFiles))

  ### 3) Apply change
  replacePatternInMetadataFiles(selectedMetadataFiles, old_pattern, new_pattern)

  # Successfull end
  print("")
  print("Done.")
