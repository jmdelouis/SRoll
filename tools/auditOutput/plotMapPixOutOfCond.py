##############################################################################
# plotMapPixOutOfCond.py
#
# !!!WARNING!!! This script need some fix (path update...)
#
# This python script allow to display all pixels that have been take appart
# during sroll process due to cond<10 not meet.
# This script also extract indexes of problematic pixel from the COND map.
#
#
# This script is to be used either on M3 or M4, since there is no dep.
# But be careful to path to data!
#
# Author  : Christian Madsen
# Date    : 2015-09-04
# version : Initial
##############################################################################

import numpy as np
from matplotlib.pyplot import * # for plot
import healpy as hp


### Debug and verbose mode
DEBUG = True
VERBOSE = True

### Flag controlling process
flag_display = True # True if we want to display the last map
flag_save_png = True # True to save result in png file 

### Reference directories
# Where to find the log containing pix indexes that are out of cond
inputLog_filename = "/pscratch1/cmadsen/SROLL_OUT/143_SWB_TEST/DBG_sroll_143GHz_512_proc.log_REF_FULL_NAN_PIX"
# Where to find the cond map
inputCondMap_filename = "/pscratch1/cmadsen/SROLL_OUT/143_SWB_V2/MAP_COND/143GHz_RD13beta2_REP6_survey_0_COND"
outputPNG_dir = "/pscratch1/cmadsen/SROLL_OUT/143_SWB_TEST/"

# Image filename
pngFilename_GEN = "%s.png"


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


  MAP_SIZE = 50331648

  # A - Extract pix indexes from the log file
  pixOut = np.zeros(50331648)

  # Use external command to extract list from log file
  # XXX to update !!! since new label...
  cmd="cat %s | grep Pix# | cut -d '#' -f2 | cut -d ' ' -f1 | cut -d '[' -f1" % inputLog_filename
  res = subprocess.check_output(cmd, shell=True)
  logDebug("len(res) = %d" % len(res))
  listOfIndex = res.split()
  logDebug("len(listOfIndex) = %d" % len(listOfIndex))

  
  # Loop on each pix indexes
  count = 0
  for ind_s in listOfIndex:
    # Convert string index to int
    ind = int(ind_s)
    # Debug
    if count < 10:
      logDebug("ind = %d" % ind)
      count+=1
    # Update array accordingly
    pixOut[ind] = 1


  ### Display result
  hp.mollview(pixOut, title="pixOut")
  show(block=False)


  # Final end
  log("")
  log("All done.")

  if flag_display:
    show() # Final blocking show
