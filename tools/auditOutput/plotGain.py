#!/usr/bin/env python
##############################################################################
# plotGain.py
#
# !!!WARNING!!! This script need some fix (path update...)
#
# This python script aims to plot gain obtain after a sroll process and
# optionnaly compare it to a reference.
#
# NOTE: If diff is requested, user shall have a valid NO_DMC_LIB env.
#       see instructions in NO_DMC_LIB/doc/README.
#
# Usage:
# $> ./plotGain.py <myProcessName>
#
# Author  : Christian Madsen
# Date    : 2015-07-27
# version : Initial
##############################################################################

import numpy as np
from matplotlib.pyplot import * # for plot
import healpy as hp


### Debug and verbose mode
DEBUG = True
VERBOSE = True

### Flag controlling process
flag_display = True # True if we want to display plot
flag_compare = True # Allow to compare the new gain with the reference
flag_save = False # True to save result in png file 

configHost = "M4" # can be one of "M4", "EDISON"


### Reference directories
if configHost == "M4":
  baseDir = "/pscratch1/delouis/tmp_CM_SROLL_OUT/sroll_v1.12/" # For M4
else:
  baseDir = "/scratch2/scratchdirs/cmadsen/SROLL_OUT/" # For EDISON (NERSC) 

# File to be processed: list of bolo to be take in count
# ["143-%d%s" % (n,s) for n in range(1,5) for s in ['a','b']]
bolo_pol = ['143-1a', '143-1b', '143-2a', '143-2b', '143-3a', '143-3b', '143-4a', '143-4b']
bolo_all = bolo_pol + ['143-5', '143-6', '143-7']

#RDLabel_REF="RD12full2_"
RDLabel_NEW="RD13beta2_"

#gain_REF_dir = "/data/dmc/MISS03/DATA/VECT_JMD_12/" # on M3
gain_REF_dir = "/pscratch1/delouis/SIM_data/dmc/MISS03/DATA/VECT_JMD_12/" # on m4
gain_REF_filename_GEN = "%s_offsets_PROD_REP6_RD12full2_%sGHz_GAIN" # % (bolo, bolo[:3])

# Generique filename associated to GAIN file
gain_NEW_filename_GEN = "%s_offsets_PROD_REP6_%s%sGHz_GAIN" # % (bolo, RDLabel, bolo[:3])


### output (png config)
pngFilename_GEN = "%s.png"


### External program
#montage_bin = "montage"


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


##############################################################################
# MAIN
# This code allow to execute the module as a standalone script.
##############################################################################
if __name__ == "__main__":
  import sys
  import time
  import os
  import subprocess


  ### A - Parse command line

  assert(len(sys.argv) == 3)
  
  argpos = 1

  # 1- Retrieve the prod name to be used
  currentProdName = sys.argv[argpos]
  argpos+=1
  log("Prod name to be used: '%s'" % currentProdName)

  # Construct input gain dir = where to find new gain data files
  inputGainDir = baseDir + currentProdName + "/VECT/"
  if not os.path.exists(inputGainDir):
    log("ERROR: Directory '%s' does not exist! PLEASE check your input name." % inputGainDir)
    exit(1)

  if flag_save:
    # Create output dir for png
    pngOutputDir = baseDir + currentProdName + "/png_GAIN/"
    # Create directory (if required)
    if not os.path.exists(pngOutputDir):
      os.makedirs(pngOutputDir)

  # 2- Retrieve the requested bolo subset. Can be either "POL" or "ALL".
  bolo_subset = sys.argv[argpos]
  argpos+=1
  log("Bolo subset: '%s'" % bolo_subset)

  if bolo_subset == "POL":
    currentBoloList = bolo_pol
  elif bolo_subset == "ALL":
    currentBoloList = bolo_all
  else:
    log("ERROR: unknown bolo subset!")
    exit(1)


  ### B - Loop on all gain for bolo_all (obtain using SWB data)
  for i_bolo, bolo in enumerate(currentBoloList):
    log("")
    log("** Processing bolo %s" % bolo)

    currentGainFile = gain_NEW_filename_GEN % (bolo, RDLabel_NEW, bolo[:3])
    currentGainFile_fullpath = inputGainDir + currentGainFile
    log("   Reading gain from '%s'" % currentGainFile_fullpath)
    if not os.path.isfile(currentGainFile_fullpath):
      log("WARN: Missing file... (just skip it)")
      continue
    
    gainNEW = np.fromfile(currentGainFile_fullpath, dtype='f8') # Values are 'double'

    logDebug("gainNEW.size = %d" % gainNEW.size)
    logDebug(gainNEW)

    # Compute the mean (TOSEE XXX if we need to remove zeros from the mean...)
    gainNEW_mean = gainNEW.mean()
    logDebug("gainNEW.mean = %g" % gainNEW_mean)

    # Plot data around the mean
    delta = 5e-3
    #ylim(gainNEW_mean-delta, gainNEW_mean+delta)
    ylim(1.001, 1.011)
    #title(currentGainFile)
    title("%s GAIN from '%s'" % (bolo, currentProdName))
    #plot(gainNEW, label="%s GAIN from '%s'" % (bolo, currentProdName))
    plot(gainNEW, label="new")
    legend(loc='lower right')

    # if requested we read reference gain and make the diff
    if flag_compare and bolo in bolo_pol:
      import nodmclib
      currentGainREFFile = gain_REF_filename_GEN % (bolo, bolo[:3])
      currentGainREFFile_fullpath = gain_REF_dir + currentGainREFFile
      log("   Reading gain REF from '%s'" % currentGainREFFile_fullpath)
      gainREF = nodmclib.read_PIODOUBLE(currentGainREFFile_fullpath, 240, gainNEW.size)
      logDebug("gainREF.size = %d" % gainREF.size)
      logDebug(gainREF)
      plot(gainREF, label="ref")

      gainDIFF = gainREF - gainNEW
      plot(gainDIFF+gainNEW_mean+0.004, label="DIFF")
      legend(loc='lower right')
      

    # Saving png file
    if flag_save:
      log("   Saving png file...")
      tmpFilename = pngOutputDir + pngFilename_GEN % currentGainFile
      savefig(tmpFilename)

    if flag_display:
      #show(block=False)
      show() # blocking display

    # Close to ensure that only one plot will be draw at a time
    close()
  
  # Before exiting check if we must keep plot displayed
  if flag_display:
    show() # this one is blocking
