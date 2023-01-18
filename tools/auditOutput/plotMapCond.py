##############################################################################
# plotMapCond.py
#
# !!!WARNING!!! This script need some fix (path update...)
#
# This python script allow to generate a MAP corresponding to the cond()
# of the corresponding matrix obtain with the pattern:
# M = [[II,IQ,IU],
#      [IQ,QQ,QU],
#      [IU,QU,UU]]
# See http://docs.scipy.org/doc/numpy/reference/generated/numpy.linalg.cond.html
# An optional png file can be generated.
#
#
# This script is to be used either on M3 or M4, since there is no dep.
# But be careful to path to data!
#
# Author  : Christian Madsen
# Date    : 2015-09-03
# version : Initial
##############################################################################

import numpy as np
from matplotlib.pyplot import * # for plot
import healpy as hp
from numpy import linalg as LA # module required for cond() function


### Debug and verbose mode
DEBUG = True
VERBOSE = True

### Flag controlling process
flag_display = True # True if we want to display the last map
#flag_diffDownscaleWithUdGrade = True # True to downscale diff map using ud_grade()
flag_save_png = True # True to save result in png file 


# Config for downscale (see ud_grade())
downScaleNSIDE = 64
downScalePOWER = 0
downScaleMAPLabel = "downscale%d_power%d_" % (downScaleNSIDE, downScalePOWER)


### Reference directories
# Where to find ref MAPs
#ref_MapDir = "/data/dmc/MISS03/DATA/MAP_JMD_2048_PROD13/"
# Where to find new MAPs
#new_MapDir = "/redtruck/SimuData/CM_temp/" # for M3
new_MapDir = "/pscratch1/cmadsen/SROLL_OUT/" # for M4


### The range of data to be processed
freq_list = ["100", "143"]
yearOrSurvey = ["survey_0", "year_0", "year_1"]
#extensions = ["I", "Q", "U"]
cond_extensions = ["II","IQ","IU","QQ","QU","UU"]
#DEBUG XXX
freq_list = ["143"]
yearOrSurvey = [yearOrSurvey[0]]
yos = yearOrSurvey[0] # force to "survey_0"

### Dir/File label (depending on freq)
dirNames_GEN = {"100":"100",
                "143":"143_SWB_V2"}

# Note that labels are to be specialized with % (freq, yos, type)
# Input MAP file are like: 143GHz_RD12full2_REP6_corr_survey_0_[II,....,QU]
fileLabels_GEN = {"100":["%sGHz_RD12full2_REP6_corr_%s_%s", "%sGHz_RD12full2_REP6_corr_%s_%s"], # [ref,new]   ex: 100GHz_RD12full2_REP6_corr_survey_0_I
                  "143":["%sGHz_RD12full2_REP6_corr_%s_%s", "%sGHz_RD13beta2_REP6_%s_%s"]} # ex: new = 143GHz_RD13beta2_REP6_survey_0_II

# image filename
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


  # Loop on all required freq
  for freq in freq_list:
    log("")
    log("Processing freq '%s'" % freq)

    currentNewDir = new_MapDir + dirNames_GEN[freq] + "/"
    logDebug("  > Directory for input MAPs: \"%s\"" % currentNewDir)

    # Where the cond map will be saved
    currentMapRefCondDIR = currentNewDir + "MAP_COND/"
    if not os.path.exists(currentMapRefCondDIR):
      logDebug("currentMapRefCondDIR = '%s'" % currentMapRefCondDIR)
      os.makedirs(currentMapRefCondDIR)
    log("  * 'MAP ref cond' will be saved in '%s'" % currentMapRefCondDIR)
    log("")

#    if flag_save_png:
#      currentPngDir = currentNewDir +  "png_DIFF_MAP" + "/" # where all images will be saved!
#      # Ensure directory exist
#      if not os.path.exists(currentPngDir):
#        logDebug("currentPngDir = '%s'" % currentPngDir)
#        os.makedirs(currentPngDir)
#      log("  * Images will be saved in '%s'" % currentPngDir)
#      log("")

    ### 1) Load all input maps (x6)
    all_maps = {}

    for type in cond_extensions: # Loop on all input MAP to load them

      # Construct the filename
      filename = fileLabels_GEN[freq][1] % (freq, yos, type)
      filename_fullpath = currentNewDir + "MAP/" + filename
      log("  -> Input MAP filename for %s: '%s'" % (type, filename_fullpath))

      # Read input map data
      logDebug("loading map...")
      dataMap = np.fromfile(filename_fullpath, dtype='f') # input map data are float!
      logDebug("     Map Data size = %d   -> nside %d" % (dataMap.size, hp.npix2nside(dataMap.size)))
      dataMap[dataMap == hp.UNSEEN] = np.nan

      # Store it temporary in the appropriate object
      all_maps[type] = dataMap

    logDebug("All maps have been loaded!")

    ### 2) Conpute the cond
    size = len(all_maps['II'])
    logDebug("len of map data: %s" % size)
    dataCond = np.empty(size)

    # Loop on all pixels
    for i in range(dataCond.size):
      # For each pixel create the matrix and compute the corresponding cond.
      mat = [[all_maps['II'][i], all_maps['IQ'][i], all_maps['IU'][i]],
             [all_maps['IQ'][i], all_maps['QQ'][i], all_maps['QU'][i]],
             [all_maps['IU'][i], all_maps['QU'][i], all_maps['UU'][i]]]

      ##### DEBUG ######
      #log("%d / %d" % (i, dataCond.size))
      #percent = (float(i)/dataCond.size)*100
      #logDebug("percent = %d" % percent)
      #log("%d%%" % percent)
      #log("mat = %s" % mat)
      if i % 503316 == 0:
        log("%d%%" % ((float(i)/dataCond.size)*100))
      ##################
      
      # Store value in result
      try:
        dataCond[i] = LA.cond(mat)
      except RuntimeWarning:
        log("!!!! EXCEPTION !!!!")
        log("i = %d" % i)
        log("mat = %s" % mat)
        log("")

    dataCond_filename = fileLabels_GEN[freq][1] % (freq, yos, 'COND')

    # Save result in file
    dataCond.tofile(currentMapRefCondDIR + dataCond_filename)

    ### 3) Optional: save the map as png

    if flag_display or flag_save_png:
      hp.mollview(dataCond, title=dataCond_filename)

    if flag_save_png:
      logDebug("   . Saving figure to image file...")
      tmpFilename = currentMapRefCondDIR + pngFilename_GEN % dataCond_filename
      savefig(tmpFilename)

    # Handle display

    if flag_display:
      show(block=False)
    else: # free memory by closing mollview
      close()


  # Final end
  log("")
  log("All done.")

  if flag_display:
    show() # Final blocking show
