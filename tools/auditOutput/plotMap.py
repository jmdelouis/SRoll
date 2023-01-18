##############################################################################
# plotMap.py
#
# !!!WARNING!!! This script need some fix (path update...)
#
# This python script allow to generate the png images corresponding to
# generated map from sroll.
# Additionaly it can produce "diff" maps between two sroll executions.
#
#
# Author  : Christian Madsen
# Date    : 2015-XX-XX
# version : Initial
##############################################################################

import numpy as np
from matplotlib.pyplot import * # for plot
import healpy as hp


### Debug and verbose mode
DEBUG = False
VERBOSE = True

### Flag controlling process
flag_display = True # True if we want to display the last map
flag_compare = False # True to compare new map with the REF map (DIFF).
flag_diffDownscaleWithUdGrade = True # True to downscale diff map using ud_grade()
flag_saveMontage_diff = True # True to save all results in a png grid file 
flag_save = False # True to save result in png file 

configHost = "M4" # can be one of "M4", "EDISON"
configType = "REAL" # can be one of "REAL", "SIM"


# Config for downscale (see ud_grade())
downScaleNSIDE = 64
downScalePOWER = 0
downScaleMAPLabel = "downscale%d_power%d_" % (downScaleNSIDE, downScalePOWER)


### Reference directories
if configHost == "M4":
  #baseDir = "/pscratch1/cmadsen/SROLL_OUT/" # For M4
  baseDir = "/pscratch1/delouis/tmp_CM_SROLL_OUT/sroll_vlast/" # For M4
else:
  baseDir = "/scratch2/scratchdirs/cmadsen/SROLL_OUT/" # For EDISON (NERSC) 
  if flag_saveMontage_diff:
    print("WARNING: 'flag_saveMontage_diff' option is not available on this host!")
    exit(1)

MAP_SUBDIR = "MAP/" # where to find map from mapBaseDir

### The list of maps to be proceed
freq_list = ["100", "143", "217", "353", "545", "857"]
detOrFull = ["-Detset1", "-Detset2", "GHz", "-Detset3"]
yearOrSurvey = ["survey_0", "year_0", "year_1"]
extensions = ["I", "Q", "U"]

# **DEBUG**: below is config for debug or particular purpose...
freq_list = ["857"]
#detOrFull = ["-1", "-2", "-3", "-4", "GHz"]
detOrFull = ["GHz"]
extensions = ["I"]
#yearOrSurvey = [yearOrSurvey[0]]
# **END DEBUG**

### The map to be processed
# The generic directory where to find map base directory
#mapBaseDir_GEN = "%s_HM2_REP8" # % freq
mapBaseDir_GEN = "%s" # % freq
# Some label in name of file map
mapRDLabel="RD13beta2_REP6"

### (optional) the reference map to be compared with
DIFF_LABEL = "_DIFF_REMDIP" # Label to be used in png name (should express idea of REF-NEW)
mapBaseDir_REF_GEN = "%s_HM2_REP8_REMDIP32" # % freq
#mapRDLabel_REF="RD12full2_"
mapRDLabel_REF="RD13beta2_"

### Mollview config
mollview_param_xsize = 2000 # Default is 800

### Image filename extention
pngFilename_GEN = "%s.png"

### External program (can be used on M4, but not available on M3)
montage_bin = "montage"



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


  # Init values
  globalMax = None

  # Loop on all required freq
  for freq in freq_list:
    log("")
    log("Processing freq '%s'" % freq)

    # Construct the mapBaseDir
    mapBaseDir = baseDir + mapBaseDir_GEN % freq + "/"
    mapDir = mapBaseDir + MAP_SUBDIR

    # Construct the mapBaseDirREF
    mapBaseDir_REF = baseDir + mapBaseDir_REF_GEN % freq + "/"
    mapDir_REF = mapBaseDir_REF + MAP_SUBDIR

    # Define dir to output png
    currentPngDir = mapBaseDir +  "png_MAP" + "/" # where all images will be saved!
    if flag_save:
      # Ensure directory exist
      if not os.path.exists(currentPngDir):
        #logDebug("currentPngDir = '%s'" % currentPngDir)
        os.makedirs(currentPngDir)
      log("  * Images will be saved in '%s'" % currentPngDir)
      log("")

    # Loop on all file matching the current freq
    for dof in detOrFull:
      for yos in yearOrSurvey:
        for type in extensions:
          # Remember that there is no Q no U for detset3!
          # And no detset3 fro 100GHz
          if dof=="-Detset3" and (type != "I" or freq == "100"):
            continue

          m_input = []

          # Construct the filename
          if configType == "SIM":
            filename = "%s%s_00000_SIM_corr_%s_%s" % (freq, dof, yos, type)
          else:
            filename_baseGEN = "%s%s_%s%s_corr_%s_%s" % (freq, dof, "%s", "%s", yos, type)
            if ("HR1" in mapBaseDir_GEN):
              halfRingLabel = "_1st"
            elif "HR2" in mapBaseDir_GEN:
              halfRingLabel = "_2nd"
            else:
              halfRingLabel = "" # By default empty, unless we are using HR data

            filename_NEW = filename_baseGEN % (mapRDLabel,halfRingLabel)

            if ("HR1" in mapBaseDir_REF_GEN):
              halfRingLabel = "_1st"
            elif "HR2" in mapBaseDir_REF_GEN:
              halfRingLabel = "_2nd"
            else:
              halfRingLabel = "" # By default empty, unless we are using HR data
            filename_REF = filename_baseGEN % (mapRDLabel_REF,halfRingLabel)
            #filename_NEW = "%s%s_%s_corr_%s_%s" % (freq, dof, mapRDLabel, yos, type)
            #filename_REF = "%s%s_%s_corr_%s_%s" % (freq, dof, mapRDLabel_REF, yos, type)
          log("  -> MAP filename_NEW: '%s'" % filename_NEW)

          # Read data
          ts = time.clock()
          data = np.fromfile(mapDir + filename_NEW, dtype='f') # new data are float!
          te = time.clock()
          logDebug("     (Loaded in %gs)" % (te-ts))
          logDebug("     Data size = %d   -> nside %d" % (data.size, hp.npix2nside(data.size)))

          # Force unseen pixel to be treat as nan...
          data[data == hp.UNSEEN] = np.nan

          #hp.mollview(data, title=filename_NEW + " (" + labelMap + ")")
          hp.mollview(data, xsize=mollview_param_xsize, title=filename_NEW)

          if flag_save:
            logDebug("   . Saving figure to image file...")
            tmpFilename = currentPngDir + pngFilename_GEN % filename_NEW
            savefig(tmpFilename)
          
            # Update m_input
            #m_input.append(baseDir + "png_" + labelMap_REF + "/" + pngFilename_GEN % filename_REF) # ref map
            m_input.append(tmpFilename) # new map

          if flag_compare and dof!="-Detset3":
            log("  -> MAP filename_REF: '%s'" % filename_REF)
            #dataREF = np.fromfile(mapDir_REF + filename_REF, dtype='f8') # ref read as double!
            dataREF = np.fromfile(mapDir_REF + filename_REF, dtype='f') # new ref read as float!
            dataREF[dataREF == hp.UNSEEN] = np.nan

            mapDiff = dataREF - data
            #print "Diff map: ", mapDiff
            currentMax = np.nanmax(np.abs(mapDiff))
            log("     >> Max diff (abs val): %g" % currentMax)
            # Update global max
            if not globalMax or currentMax > globalMax:
              globalMax = currentMax

            if flag_diffDownscaleWithUdGrade:
              mapDiff_64 = hp.ud_grade(mapDiff, nside_out=downScaleNSIDE, pess=False, power=downScalePOWER)
                
            # NOTE: the hist is only available is not all nan!
            if np.isnan(currentMax):
              log("     >> Since Diff is all Nan we do not plot anything!")
            else:
              # norm="HIST" allow to see easily the diff
              for mollviewNorm in [None, "HIST"]:
                mollviewNormLabel = "(HIST)" if mollviewNorm=="HIST" else ""
                mollviewNormLabelFile = "HIST_" if mollviewNorm=="HIST" else ""
                hp.mollview(mapDiff, xsize=mollview_param_xsize, title="mapDiff (REF-NEW) %s" % mollviewNormLabel, norm=mollviewNorm)

                if flag_save:
                  log("   . Saving figure to image file...")
                  tmpFilename = currentPngDir + "DIFF_" + DIFF_LABEL + mollviewNormLabelFile + pngFilename_GEN % filename_baseGEN % ("","")
                  savefig(tmpFilename)
                  m_input.append(tmpFilename) # new diff map

                ### Now also produce diff downscale map (if required)
                if flag_diffDownscaleWithUdGrade:
                  hp.mollview(mapDiff_64, xsize=mollview_param_xsize, title="mapDiff (REF-NEW) %s %s" % (mollviewNormLabel, downScaleMAPLabel), norm=mollviewNorm)
                  if flag_save:
                    log("   . Saving figure to image file...")
                    tmpFilename = currentPngDir + "DIFF_" + DIFF_LABEL + mollviewNormLabelFile + downScaleMAPLabel + pngFilename_GEN % filename_baseGEN % ("","")
                    savefig(tmpFilename)
                    m_input.append(tmpFilename) # new diff downscale map

          if flag_display:
            show(block=False)
          else: # free memory by closing mollview
            close()

          # Do a montage with related images (if available)
          if flag_saveMontage_diff and dof!="-Detset3":
            m_output = currentPngDir + "grid_" + pngFilename_GEN % filename_baseGEN % ("","")
            cmd = "%s -mode concatenate -tile 2x3 %s %s" % (montage_bin, ' '.join(m_input), m_output)
            print("Do: %s" % cmd)
            subprocess.call(cmd, shell=True)


  # Display globalMax
  if flag_compare and dof!="-Detset3":
    log("")
    log("!! globalMaxDiff (abs) = %g" % globalMax)


  # Final end
  log("")
  log("All done.")

  if flag_display:
    show() # Final blocking show
