#!/usr/bin/env python
##############################################################################
# doQsub.py
#
# This python script allow to configure then launch qsub job accordingly to
# parameters pass on command line (preventing the user to edit manually the
# parameter file and the qsub script...).
#
# IMPORTANT: Some required config must be load from "doQsub.cfg", therefore,
#            ensure to have it in the current directory.
#            If none are present please start from copying the file
#            "doQsub.cfg_TO_SPECIALIZE":
#            $> cp doQsub.cfg_TO_SPECIALIZE doQsub.cfg
#            Then update values as you need.
#            Note: "doQsub.cfg" must not be commited to CVS (since it may vary
#                  for each individual user environment).
#
# Usage:
# $> python doQsub.py <outputDirName> <freqGhz> <type> <subType> <procNumber> [<walltime>] [<specializedParam> ...]
#
# Exemples:
# $> python doQsub.py myProcessName 143GHz REAL "" 512
# or
# $> python doQsub.py myProcessName 143GHz REAL_HR1 "" 512
# or
# $> python doQsub.py myProcessName 143GHz REAL_HR1_detset1 "" 512
# or
# $> python doQsub.py myProcessName 143GHz SIM SWB_HALFDET1 512
# or optionnaly (setting a walltime):
# $> python doQsub.py myProcessName 143GHz SIM SWB_HALFDET1 512 04:00:00
# or optionnaly (setting a walltime and setting some specialized params):
# $> python doQsub.py myProcessName 143GHz SIM SWB_HALFDET1 512 04:00:00 toto=123 zozo="tu/ti"
# Where:
#   - 'myProcessName' is the name used in the output directory.
#   - '143GHz' is the desired bolo freq.
#   - 'SIM' specify we want to use the SIM parameters file. Other option is 'REAL'.
#     In addition a suffix of '_HR1' respectively '_HR2' force to use HalfRing!
#     In addition a suffix of '_detsetX' (where X is a number) force to use only bolos in the specified detset!
#   - 'SWB_HALFDET1' specify the parameter "sub type". Can be "" (for empty).
#   - '512' is the total number of procs you want to use for this run.
#   - '04:00:00' is the desired max walltime
#   - 'toto=123' express a specialized parameter in the form key=value so that
#                the corresponding value will be substitute in generated param file.
#                You can specify as many as you like.
#                LIMITATION: for now only the parameter for which the value has
#                            been already defined as [[%value]] can be specialized.
#
# How it works:
# The script will detect on which host you are executing this script (mainly
# either M4 or EDISON (at NERSC)) then it will load appropriate generic
# parameter file (based on your command arguments) then it will specialized it
# and finally execute the qsub command.
#
# !!! IMPORTANT !!!:
# Remember to have an updated config file "doQsub.cfg".
#
# Author  : Christian Madsen
# Date    : 2015-06-02
# version : Initial
##############################################################################

import ConfigParser


### Debug and verbose mode
DEBUG = True
VERBOSE = True

### Flag controlling process
#flag_XXX = False # True if we want to 
#flag_XXX = False # True to 
#flag_XXX = True # True to 


### Constants
DEFAULT_WALLTIME = "10:00:00"
DEFAULT_DETSET = "" # "" mean 'full' (ie. use all available bolo), otherwise user can use "detset[123]"...

# DEFAULT Specializable parameters: can be override by specific def
# 'specializableParams_SPE' or by user on commande line.
# Only the parameter for which the value is define like [[%value]] in
# the generic parameter file can be specialized!
specializableParams = { "[[%BeginRing]]": "240",
                        "[[%EndRing]]": "26050",
                        "[[%GAINSTEP]]": "128",  # Default is "128"
                        "[[%REMDIP]]": "1",      # Default is "8"
                        "[[%NADU]]": "1",        # Default is "50"
                        "[[%SEED]]": "1234",     # Default is "1234"
                        "[[%REPID]]": "REP6",    # Default is "REP6"
                        "[[%DOGAINDIP]]": "1",   # Default is "1"
                        "[[%OddEven]]": ""       # Default is "" (empty mean no parity = full)
                      }

# SPECIFIC values for each freq
specializableParams_SPE = { 
                            "353": {
                              "[[%GAINSTEP]]": "32"
                              #"[[%NADU]]": "32" # just for test instead of 50!
                            },
                            "545": {
                              "[[%GAINSTEP]]": "32",
                              #"[[%REMDIP]]": "1",
                              "[[%DOGAINDIP]]": "0",
                              #"[[%NADU]]": "32"
                            },
                            "857": {
                              "[[%GAINSTEP]]": "32",
                              #"[[%REMDIP]]": "1",
                              "[[%DOGAINDIP]]": "0",
                              #"[[%NADU]]": "32"
                            }
                          }

# Directories
# The subdirectory where are saved MAPs
MAP_DIRNAME = "MAP"
# The subdirectory where are saved other data
VECT_DIRNAME = "VECT"


# The structure to store all config
class HostConfig:
  """Class representing the configuation of a server (host for executing sroll)."""

# Here are define the corresponding qsub required info regarding the host on
# which we want to launch sroll

### M3 config
Host_M3 = HostConfig()
Host_M3.id = "ln" # Node name identifying this host
Host_M3.name = "M3" # Human readable host name. See related section in "doQsub.cfg".
Host_M3.nbProcPerNode = 8 # Number of procs per node
Host_M3.qsubScript = "do.qsub_GEN_M3" # The script to be used
Host_M3.qsubQueueName = "" # Not used for M3!

### M4 config
Host_M4 = HostConfig()
Host_M4.id = "log" # Node name identifying this host
Host_M4.name = "M4" # Human readable host name. See related section in "doQsub.cfg".
Host_M4.nbProcPerNode = 24 # Number of procs per node
Host_M4.qsubScript = "do.qsub_GEN_M4" # The script to be used
Host_M4.qsubQueueName = "" # Not used for M4!

### EDISON (NERSC) config
Host_EDISON = HostConfig()
Host_EDISON.id = "edison" # Node name identifying this host
Host_EDISON.name = "Edison (NERSC)" # Human readable host name. See related section in "doQsub.cfg".
Host_EDISON.nbProcPerNode = 24 # Number of procs per node
Host_EDISON.qsubScript = "do.qsub_GEN_EDISON" # The script to be used
Host_EDISON.qsubQueueName = "regular" # Option for the qsub process (may be used to specify the queue as 'debug' or 'regular')

### Definition of all hosts to be take in count
host_list = [Host_M3, Host_M4, Host_EDISON]


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

def query_yes_no(question, default="yes"):
  """Ask a yes/no question via raw_input() and return the answer.

  "question" is a string that is presented to the user.
  "default" is the presumed answer if the user just hits <Enter>.
      It must be "yes" (the default), "no" or None (meaning
      an answer is required of the user).

  The "answer" return value is True for "yes" or False for "no".
  """
  valid = {"yes": True, "y": True, "ye": True,
           "no": False, "n": False}
  if default is None:
    prompt = " [y/n] "
  elif default == "yes":
    prompt = " [Y/n] "
  elif default == "no":
    prompt = " [y/N] "
  else:
    raise ValueError("Invalid default answer: '%s'" % default)

  while True:
    sys.stdout.write(question + prompt)
    choice = raw_input().lower()
    if default is not None and choice == '':
      return valid[default]
    elif choice in valid:
      return valid[choice]
    else:
      sys.stdout.write("Please respond with 'yes' or 'no' "
                       "(or 'y' or 'n').\n")

#-----------------------------------------------------------------------------

# TODO to check/complete then used!
def myParseArg():
  """Parse the command line.
  """
  import argparse

  parser = argparse.ArgumentParser()

  parser.add_argument("currentProcessLabel",
                      help="Name to be used for the process. It will be used in output dir name.")
  parser.add_argument("bolo",
                      choices=["143GHz"],
                      help="The bolo id. ex: '143GHz'")
  parser.add_argument("data",
                      choices=["REAL", "SIM"],
                      help="Input data type: 'REAL' or 'SIM'")
  parser.add_argument("nbProc", type=int,
                      help="Number of proc required.")
  parser.add_argument("-w", "--walltime",
                      help="The max walltime. ex '10:00:00'")
  parser.add_argument("-s", "--swb", action="store_true",
                      help="This option allow to add the SWB bolometers (using appropriate parameter file).")
                      
  args = parser.parse_args()

  if args.walltime:
    logDebug("OPTION walltime activated with: %s" % args.walltime)
  

  #TODO DEBUG
  exit(0)

  pass

#-----------------------------------------------------------------------------

def getBoloId(boloIdFull):
  """Convert "50_143_3a" to "143-3a".
  """

  return boloIdFull[3:].replace("_","-")

#-----------------------------------------------------------------------------

def preProcessParameterFile(paramGEN_Filename, freq, currentDetset, flag_SubTypeSWB):
  """Allow to Pre-Process the parameter file to add specific entry according to
  user request (detset...).
  Mainly it load the paramFile.txt then substitute each {} entry by the
  corresponding code, and return the result as a string var for next step (specialization).

  Parameters:
  ----------
  paramGEN_Filename : Filename (fullpath) to the generic parameter file (the one to be proceceed)
  freq: The bolo freq to be considered
  currentDetset: The current detset to be used. "" mean no detset -> use of all available bolo
  flag_SubTypeSWB: True if SWB mode is activated, False otherwise.
  """

  logDebug(">> freq = %s" % freq)
  logDebug(">> currentDetset = %s" % currentDetset)
  logDebug(">> flag_SubTypeSWB = %s" % flag_SubTypeSWB)

  # This dictionary define position of each component in the array:
  boloComponents = {
    "detset": 0,
    "Calibration": 1,
    "CrossPol": 2,
    "FSLCOEF": 3,
    "Monop": 4,
    "NEP": 5
  }

  bolo_list = {
    "100": (
      ("00_100_1a", ["detset1", "1.00340513024e-13", "0.0272", "1.0", "2.2605358109e-15", "2.4767234259e-16"]),
      ("01_100_1b", ["detset1", "1.23240311873e-13", "0.0293", "1.0", "6.33346357328e-15", "2.56729781043e-16"]),
      ("80_100_4a", ["detset1", "1.46826316998e-13", "0.0219", "1.0", "6.29208169223e-16", "2.26421345033e-16"]),
      ("81_100_4b", ["detset1", "1.17469574177e-13", "0.0402", "1.0", "3.36970326506e-15", "2.45551518841e-16"]),
      ("20_100_2a", ["detset2", "1.51966997388e-13", "0.0195", "1.0", "2.33868290562e-15", "1.84486488103e-16"]),
      ("21_100_2b", ["detset2", "1.5885889296e-13", "0.0513", "1.0", "-1.13641813483e-15", "2.68522224152e-16"]),
      ("40_100_3a", ["detset2", "1.38652193689e-13", "0.0521", "1.0", "5.53657151251e-15", "1.45183188975e-16"]),
      ("41_100_3b", ["detset2", "1.1492356927e-13", "0.0339", "1.0", "4.45056923088e-15", "1.46602765997e-16"])
    ),
    "143": (
      ("02_143_1a", ["detset1", "1.88008633662e-13", "0.0915", "1.0", "7.3995195043e-15", "1.43494770343e-16"]),
      ("03_143_1b", ["detset1", "1.61678010691e-13", "0.0835", "1.0", "5.32325733232e-15", "1.91452982801e-16"]),
      ("50_143_3a", ["detset1", "1.79339046635e-13", "0.0875", "1.0", "5.75907091845e-15", "1.49002842269e-16"]),
      ("51_143_3b", ["detset1", "1.61409996829e-13", "0.053", "1.0", "6.44604950677e-15", "1.29914524043e-16"]),
      ("30_143_2a", ["detset2", "1.78260067e-13", "0.067", "1.0", "5.29193929445e-15", "1.34662861911e-16"]),
      ("31_143_2b", ["detset2", "1.82131500475e-13", "0.0572", "1.0", "2.9348297266e-15", "1.47483374152e-16"]),
      ("82_143_4a", ["detset2", "1.66035966523e-13", "0.0357", "1.0", "1.99766939061e-15", "1.46058872792e-16"]),
      ("83_143_4b", ["detset2", "1.55447052831e-13", "0.0371", "1.0", "2.69626805994e-15", "1.49667609571e-16"]),
      ("10_143_5" , ["detset3", "2.65238565653e-13", "0.882", "1.0", "5.3087867885e-15", "1.8148147328e-16"]),
      ("42_143_6" , ["detset3", "2.381597318e-13", "0.92", "1.0", "4.95338905598e-15", "1.65716991561e-16"]),
      ("60_143_7" , ["detset3", "2.57719266261e-13", "0.974", "1.0", "4.67336488715e-15", "1.57075016643e-16"])
    ),
    "217": (
      ("11_217_5a", ["detset1", "1.14058850217e-13", "0.0256", "1.0", "2.80761299651e-15", "1.67331426436e-16"]),
      ("12_217_5b", ["detset1", "1.15245126451e-13", "0.0246", "1.0", "8.9805619758e-15", "1.46628673336e-16"]),
      ("61_217_7a", ["detset1", "1.24037931807e-13", "0.0307", "1.0", "5.25596376897e-15", "1.54985747982e-16"]),
      ("62_217_7b", ["detset1", "1.19572741647e-13", "0.0327", "1.0", "4.02970261138e-15", "1.46628673336e-16"]),
      ("43_217_6a", ["detset2", "1.15733750737e-13", "0.026", "1.0", "3.04610354346e-15", "1.65906925076e-16"]),
      ("44_217_6b", ["detset2", "1.17133951737e-13", "0.0236", "1.0", "6.28524515573e-15", "1.59924019363e-16"]),
      ("71_217_8a", ["detset2", "1.20172971416e-13", "0.0298", "1.0", "5.61649675718e-15", "1.69800562127e-16"]),
      ("72_217_8b", ["detset2", "1.1391920766e-13", "0.0303", "1.0", "1.42911330698e-15", "1.65906925076e-16"]),
      ("04_217_1", ["detset3", "1.68882666765e-13", "0.926", "1.0", "4.19798716956e-15", "1.52801512563e-16"]),
      ("22_217_2", ["detset3", "1.62797073142e-13", "0.961", "1.0", "7.21865509517e-15", "1.52706545805e-16"]),
      ("52_217_3", ["detset3", "1.71515434481e-13", "0.924", "1.0", "9.17652991993e-15", "1.52136745261e-16"]),
      ("84_217_4", ["detset3", "1.66584754679e-13", "0.916", "1.0", "1.12013312388e-14", "1.3931623302e-16"])
    ),
    "353": (
      ("23_353_3a", ["detset1", "3.505451069e-14", "0.0583", "1.0", "5.56785361889e-15", "1.9971509069e-16"]),
      ("24_353_3b", ["detset1", "3.46195375785e-14", "0.0413", "1.0", "2.18112475466e-15", "1.6324785587e-16"]),
      ("53_353_5a", ["detset1", "3.31816651442e-14", "0.0829", "1.0", "1.68104224305e-15", "1.56030382312e-16"]),
      ("54_353_5b", ["detset1", "3.30720800648e-14", "0.0661", "1.0", "4.36286171304e-15", "1.57169983401e-16"]),
      ("32_353_4a", ["detset2", "2.89764285963e-14", "0.0683", "1.0", "1.23924007679e-15", "1.40075967079e-16"]),
      ("33_353_4b", ["detset2", "2.89508885946e-14", "0.0439", "1.0", "5.55355317474e-15", "1.4577397252e-16"]),
      ("63_353_6a", ["detset2", "2.38419917429e-14", "0.0664", "1.0", "5.95276400648e-15", "1.62393155054e-16"]),
      ("64_353_6b", ["detset2", "2.16742504299e-14", "0.0595", "1.0", "1.04699754947e-14", "1.40550800866e-16"]),
      ("05_353_1", ["detset3", "5.81896764062e-14", "0.938", "1.0", "6.95575630847e-15", "1.39791066807e-16"]),
      ("13_353_2", ["detset3", "6.31788126648e-14", "0.91", "1.0", "5.18389125391e-15", "1.55840448798e-16"]),
      ("45_353_7", ["detset3", "4.79439720697e-14", "0.852", "1.0", "8.20371783001e-15", "1.41975302226e-16"]),
      ("85_353_8", ["detset3", "4.50246808415e-14", "0.855", "1.0", "5.30152463514e-15", "1.46628673336e-16"])
    ),
    "545": (
      ("14_545_1", ["detset1", "3.29932298836e-16", "0.913", "1.0", "3.82437622616e-14", "2.06267796946e-16"]),
      ("34_545_2", ["detset1", "3.08792864906e-16", "0.894", "1.0", "2.57192447573e-14", "1.7863247056e-16"]),
      ("73_545_4", ["detset2", "2.64855405172e-16", "0.89", "1.0", "2.82843952366e-14", "1.63152889113e-16"])
    ),
    "857": (
      ("25_857_1", ["detset1", "3.30076826046e-16", "0.887", "1.0", "2.7228011807e-14", "2.17663807827e-16"]),
      ("35_857_2", ["detset1", "3.55811287601e-16", "0.883", "1.0", "3.20205460445e-14", "2.36372259024e-16"]),
      ("65_857_3", ["detset2", "3.18681631353e-16", "0.856", "1.0", "2.14043411258e-14", "2.08546999122e-16"]),
      ("74_857_4", ["detset2", "2.219187708e-16", "0.897", "1.0", "3.3436884405e-14", "2.01804359351e-16"])
    )
  } 

  import collections
  import re

  tmpOrderedDict = collections.OrderedDict(bolo_list[freq])
  #logDebug("tmpOrderedDict = %s" % tmpOrderedDict)

  # 0) Prepare data to be used (only a subset of bolo in case of detset!)
  if currentDetset: # Must retrict to the specified detset
    currentBoloSet = collections.OrderedDict((k, v) for k, v in tmpOrderedDict.iteritems() if v[0] == currentDetset or (flag_SubTypeSWB and v[0] == "detset3"))
  else: # Case of full bolo
    currentBoloSet = tmpOrderedDict

  #logDebug("currentBoloSet = %s" % currentBoloSet)

  if freq in ("545", "857"):
    flag_freqHF = True
    # Number of maps are equal to: nb bolo + 1 (one per bolo and one containing all bolo)
    numberofMaps = len(currentBoloSet.keys()) + 1
  else:
    flag_freqHF = False
    numberofMaps = 3
    if flag_SubTypeSWB: # Add a map set with only SWB activated
      numberofMaps += 1

  # 1) Load generic file
  with open (paramGEN_Filename, "r") as myfile:
    param_GEN = myfile.read()

  # 2) Loop on all file and replace each {{xxx}} by its corresponding content
  componentsCounter = collections.Counter()
  newparam_GEN = ""
  for line in param_GEN.splitlines():
    #logDebug("line = %s" % line)

    # a) Avoid comment (line starting with #) (avoid any false pattern detection)
    if line.strip().startswith("#") or not line.strip():
      if not line.strip().startswith("#%%"): # Allow to remove old Header
        newparam_GEN+=line+"\n"
      #logDebug("Comments or empty line detected!!!")
      continue

    # b) Replace any {{}} with appropriate content
    match = re.search(r'{{(\w+)}}', line)
    if match: # we found something between {{ and }}
      componentLabel = match.group(0)[2:-2]
      #logDebug("*** found a pattern: '%s'" % componentLabel)
      # Case of componentLabel is one of boloComponents
      if componentLabel in boloComponents.keys():
        # boloComponents are to be extracted from the corresponding entry in bolo dictionnary
        # Loop on each bolo
        for bolo_k,bolo_v in currentBoloSet.iteritems():
          #logDebug("******")
          #logDebug("bolo_k = %s" % bolo_k)
          #logDebug("componentLabel: %s" % componentLabel)
          #logDebug("index pos for this componentLabel: %d" % boloComponents[componentLabel])
          componentsCounter[componentLabel] += 1
          newparam_GEN += "%s%d = %s\n" % (componentLabel, componentsCounter[componentLabel], bolo_v[boloComponents[componentLabel]])

      #--- Case of special components ---

      elif componentLabel == "OUT_NOPOL":
        # OLD:...
        # Test to know in which situation we are. If one bolo in list ends with no letter then SWB!
        #flag_SWB = False
        #for bolo_k in currentBoloSet.keys():
        #  if not bolo_k.endswith("a") and not bolo_k.endswith("b"):
        #    flag_SWB = True
        #    break

        # Note that for 545 and 857 there is no polarized bolo so every flag must be set at 1.
        # Furthermore the number of output maps correspond to (see also bolomask) :
        #   map_1 = all bolo
        #   map_n = only Nth bolo

        # Note2: For all other freq we have only 2 situations: the first when no SWB and the other with.

        if flag_freqHF:
          for ind in range(numberofMaps):
            componentsCounter[componentLabel] += 1
            newparam_GEN += "%s%d = %s\n" % (componentLabel, componentsCounter[componentLabel], "1")
        else:
          # In all case add these lines (corresponding to PSB)
          componentsCounter[componentLabel] = 3
          newparam_GEN += "%s%d = %s\n" % (componentLabel, 1, "0")
          newparam_GEN += "%s%d = %s\n" % (componentLabel, 2, "0")
          newparam_GEN += "%s%d = %s\n" % (componentLabel, 3, "0")
          # Last line only if SWB
          if flag_SubTypeSWB:
            componentsCounter[componentLabel] += 1
            newparam_GEN += "%s%d = %s\n" % (componentLabel, componentsCounter[componentLabel], "1")

      elif componentLabel == "bolomask":
        for ind in range(numberofMaps):
          newparam_GEN += "# Map%d\n" % (ind+1)
          for bolo_k,bolo_v in currentBoloSet.iteritems():
            componentsCounter[componentLabel] += 1
            val = "0" # Default value
            boloPos = currentBoloSet.keys().index(bolo_k) # bolo position in 'currentDetset' (which may be full bolo set)
            if flag_freqHF: # In case of High Freq
              if ind == 0: # First map is all one
                val = "1" 
              elif ind == boloPos + 1: # Otherwise only the current corresponding bolo is activated
                val = "1" 
            else:
              boloDetset = int(bolo_v[0][len("detset"):])
              if ind == 0: # all bolo to 1 (but with one exception in case of detset and SWB)
                val = "1"
                # But EXCEPTION! no SWB in detset MAP!
                if currentDetset and flag_SubTypeSWB and boloPos > 3:
                  val = "0"
              else: # otherwise it depends of user request...
                if currentDetset: # Case of a specific detset requested by user
                  if ind == 3 and boloDetset == 3:
                    val = "1"
                  elif ind == 1 and (boloPos == 0 or boloPos == 1):
                    val = "1"
                  elif ind == 2 and (boloPos == 2 or boloPos == 3):
                    val = "1"
                else: # Case of full bolo
                  if ind == boloDetset:
                    val = "1"
            newparam_GEN += "%s%d = %s\n" % (componentLabel, componentsCounter[componentLabel], val)
        
      elif componentLabel == "MAP":
        line_value = line.split("=")[1].strip()
        for ind in range(numberofMaps):
          # Update counter for each line added!
          componentsCounter[componentLabel] += 1
          if flag_freqHF:
            if componentsCounter[componentLabel] == 1:
              mapLabel = "%sGHz" % freq
            else:
              mapLabel = getBoloId(currentBoloSet.keys()[ind-1]) # Retrieve boloid for current map
          else:
            if currentDetset:
              if componentsCounter[componentLabel] == 1:
                if flag_SubTypeSWB:
                  mapLabel = "%sGHz_SWBHD%s" % (freq, currentDetset[-1:])
                else:
                  mapLabel = "%sGHz_HD%s" % (freq, currentDetset[-1:])
              elif componentsCounter[componentLabel] == 2:
                mapLabel = "%sx" % getBoloId(currentBoloSet.keys()[0])[:-1] # retrieve id number of bolo
              elif componentsCounter[componentLabel] == 3:
                mapLabel = "%sx" % getBoloId(currentBoloSet.keys()[2])[:-1] # retrieve id number of bolo
              elif componentsCounter[componentLabel] == 4:
                mapLabel = "%s-SWB" % freq
            else:
              if componentsCounter[componentLabel] == 1:
                  mapLabel = "%sGHz" % freq
              elif componentsCounter[componentLabel] == 2:
                  mapLabel = "%s-Detset%d" % (freq, ind)
              elif componentsCounter[componentLabel] == 3:
                  mapLabel = "%s-Detset%d" % (freq, ind)
              elif componentsCounter[componentLabel] == 4:
                  mapLabel = "%s-Detset%d" % (freq, ind)

          # Replace line accordingly
          #logDebug("!! Found MAP component in line: %s" % line)
          tmpline_value = line_value.replace("{mapLabel}", mapLabel)
          newline = "%s%d = %s" % (componentLabel, componentsCounter[componentLabel], tmpline_value)
          #logDebug("   >> newline: '%s'" % newline)
          newparam_GEN += newline + "\n"
        
      #--- End of special components ---

      else: # Case of metaComponents
        # metaComponents are of the form:  {{metaLabel}} = blabla{boloId}
        # Loop on each bolo
        for bolo_k in currentBoloSet.keys():
          # Update counter for each line added!
          componentsCounter[componentLabel] += 1
          # Replace line accordingly
          #logDebug("!! Found metaComponents in line: %s" % line)
          line_value = line.split("=")[1].strip()
          line_value = line_value.replace("{boloID}", getBoloId(bolo_k))
          line_value = line_value.replace("{boloIDFull}", bolo_k)
          newline = "%s%d = %s" % (componentLabel, componentsCounter[componentLabel], line_value)
          #logDebug("   >> newline: '%s'" % newline)
          newparam_GEN += newline + "\n"
        
    else: # Just keep line
      newparam_GEN+=line+"\n"

  # 3) Now Replace counter for each entry {%N}
  tmpparam_GEN = newparam_GEN
  for line in newparam_GEN.splitlines():
    #logDebug(">> line: %s" % line)
    # a) Avoid comment (line starting with #) (avoid any false pattern detection)
    if line.strip().startswith("#") or not line.strip():
      continue

    # Look for the counter pattern
    match = re.search(r'number_of_(\w+).*=.*{%N}', line.strip())
    if match: # we found something between {{ and }}
      counterLabel = match.group(1)
      #logDebug("!! Found {%%N} in line: %s" % line)
      #logDebug("   >> counterLabel = %s" % counterLabel)
      # Therefore replace {%N} by the corresponding counter
      newLine = line.replace(r"{%N}", "%d" % componentsCounter[counterLabel])
      #logDebug("   >> oldLine: '%s'" % line)
      #logDebug("   >> newLine: '%s'" % newLine)
      tmpparam_GEN = tmpparam_GEN.replace(line, newLine)

  # 4) Finally return generated content (after having prepend new header)
  header="###############################################################################\n"+ \
    "# This parameter file has been auto-generated by doQsub.py script.\n" + \
    "###############################################################################\n"

  return header+tmpparam_GEN

#-----------------------------------------------------------------------------

def specializeParameterFile(param_GEN, host, specializableParams):
  """Allow to specialize the generic parameter file, save it to output dir and
  return its pathname.
  
  Parameters:
  ----------
  param_GEN : Current content of parameter file (already pre-proceceed)
  host : host config
  specializableParams : param to be update in generic parameter file

  Return:
  ------
  The specialized parameter file path (to be used by qsub job).
  """

  import subprocess
  import os

  # 2) Specialize it accordingly
  base_dict = {"[[%srollIn]]" : host.srollIn,
               # NO MORE USED... "[[%srollOut]]" : host.srollOut + "/" +  currentProcessLabel,
               "[[%srollOut_MapDir]]" : host.srollOut_MapDir,
               "[[%srollOut_VectDir]]" : host.srollOut_VectDir,
              }
  # Concatenate both dictionaries
  dic_for_replacing = {}
  dic_for_replacing.update(base_dict)
  dic_for_replacing.update(specializableParams)

  # 3) Specialize with current bunch info
  tmp_GEN = param_GEN
  for key, val in dic_for_replacing.iteritems():
    # Check key is effectively existing in the file! (otherwise error...)
    if key in tmp_GEN:
      tmp_GEN = tmp_GEN.replace(key, val)
    else:
      print("ERROR: Invalid specializable parameter name ('%s'): Not existing in input file!" % key)
      exit(1)

  # Save specialized file
  temp_filename = host.srollOut_BaseDir+"/param.txt"
  with open(temp_filename, "w") as temp_file_object:
    temp_file_object.writelines(tmp_GEN)
  logDebug("Write TMP parameters file: '%s'" % temp_filename)

  # Finaly return name of the specialized file
  return temp_filename

#-----------------------------------------------------------------------------

def specializeQsubAndLaunchIt(host, label, nbNode, nbTotProc, logfilename, paramFile, walltime):
  """Allow to specialize the qsub script and launch it.
  """
  import subprocess

  # Command to be used for launching qsub script
  qsub_bin = "qsub"

  qsub_script_filename = host.qsubScript

  # 1) Load generic qsub script
  with open (qsub_script_filename, "r") as myfile:
    qsub_script_GEN = myfile.read()

  # 2) Specialize it accordingly
  dic_for_replacing = {"[[%label]]" : label,
                       "[[%nbNode]]" : str(nbNode),
                       "[[%nbTotProc]]" : str(nbTotProc),
                       "[[%srollDir]]" : host.srollBaseDir,
                       "[[%logfilename]]" : host.srollOut_BaseDir + "/" + logfilename, # prepend dir to get fullpath
                       "[[%paramFile]]" : paramFile, # This one is full path
                       "[[%walltime]]" : walltime,
                       "[[%queueName]]" : host.qsubQueueName # only for NERSC
                       }

  # 2.a) Specialize qsub script with current bunch info
  tmp_script_GEN = qsub_script_GEN
  for key, val in dic_for_replacing.iteritems():
    tmp_script_GEN = tmp_script_GEN.replace(key, val)
  #print("================================") 
  #print("%s" % tmp_script_GEN)
  #print("================================") 

  # Save temp script file
  temp_filename = host.srollOut_BaseDir+"/do.qsub"
  with open(temp_filename, "w") as temp_file_object:
    temp_file_object.writelines(tmp_script_GEN)
  logDebug("Write qsub script file: '%s'" % temp_filename)

  # 2.b) Launch it (note that 'qsub' command return immediatly)
  log("Launching qsub '%s' using %d nodes (%d tot proc) -> log file '%s'" % (label, nbNode, nbTotProc, logfilename))
  #cmd = "%s %s" % (qsub_bin, temp_filename)
  # Retrieve the output (= JOBID) to save it in output dir (=host.srollOut_BaseDir)
  jobId = subprocess.check_output([qsub_bin, "-W umask=0022", temp_filename])
  # Rmk: the -W umask=0022 allow to have qsub output to be readable by every one...

  # Save jobId file
  temp_filename = host.srollOut_BaseDir+"/jobID"
  with open(temp_filename, "w") as temp_file_object:
    temp_file_object.writelines(jobId)
  logDebug("Write jobID file: '%s'" % temp_filename)

#-----------------------------------------------------------------------------


##############################################################################
# MAIN
# This code allow to execute the module as a standalone script.
##############################################################################
if __name__ == "__main__":
  import sys
  import os

  print("*************")
  print("* doQsub.py *")
  print("*************")

  ### First of all detect the host on which the qsub shall be launch
  currentHostName = os.environ['HOSTNAME']
  logDebug("Hostname = '%s'" % currentHostName)

  currentHost = None
  for host in host_list:
    if currentHostName.startswith(host.id):
      log("Host detected: %s" % host.name)
      currentHost = host
      logDebug("qsub filename to be used: '%s'" % host.qsubScript)
      break

  if not currentHost:
    print("ERROR: Unknown host!")
    exit(1)

  ### Read config file and update internal object accordingly
  # Ensure config file exist
  CONFIG_FILE = "doQsub.cfg"
  if not os.path.isfile(CONFIG_FILE):
    print("ERROR: Please create and fill config file '%s'" % CONFIG_FILE)
    exit(1)

  # Parse config file
  config = ConfigParser.RawConfigParser()
  config.read(CONFIG_FILE)
  currentHost.srollBaseDir = config.get(currentHost.name, "srollBaseDir")
  currentHost.srollParamDir = config.get(currentHost.name, "srollParamDir")
  currentHost.srollIn = config.get(currentHost.name, "srollIn")
  currentHost.srollOut = config.get(currentHost.name, "srollOut")

  ### Parse command line
  #myParseArg()

  
  # Make some assertion about the minimal number of arguments
  assert(len(sys.argv) >= 6)
  argpos = 1

  # 0- Process name
  currentProcessLabel = sys.argv[argpos]
  argpos+=1
  logDebug(">> Read from command line: process label = %s" % currentProcessLabel)

  # Update accordingly all related variables
  currentHost.srollOut_BaseDir = currentHost.srollOut + "/" + currentProcessLabel
  currentHost.srollOut_MapDir  = currentHost.srollOut_BaseDir + "/" + MAP_DIRNAME
  currentHost.srollOut_VectDir = currentHost.srollOut_BaseDir + "/" + VECT_DIRNAME

  # 1- Bolo freq
  currentBoloFreq = sys.argv[argpos]
  argpos+=1
  logDebug(">> Read from command line: bolo freq = %s" % currentBoloFreq)

  # 2- Process type (REAL or SIM)
  currentProcessType = sys.argv[argpos]
  argpos+=1
  logDebug(">> Read from command line: process type = %s" % currentProcessType)
  # Compute the corresponding label: in case of 'REAL' -> '' ; and in case of 'SIM' -> '_SIM'
  if currentProcessType.startswith("REAL"):
    currentProcessTypeLabel = ""
    suffix_Type = currentProcessType[len("REAL"):] # Retrieve suffix
  elif currentProcessType.startswith("SIM"):
    currentProcessTypeLabel = "_SIM"
    suffix_Type = currentProcessType[len("SIM"):] # Retrieve suffix
  else:
    print("ERROR: unknown process type '%s'" % currentProcessType)
    exit(1)

  # 2a- Check for optional HalfRing suffix ("_HR1" or "_HR2")
  HALF_RING_pattern_replacement = {"_HR1": "REP6_1st",
                                   "_HR2": "REP6_2nd"}
  suffix_HR_found = False
  for hr_pattern in HALF_RING_pattern_replacement.keys():
    if suffix_Type.startswith(hr_pattern):
      suffix_HR_found = True
      logDebug("***** DETECT halfring config!!!")
      log("  -> Input parameter will be specialized for HALF RING (using '%s')" % HALF_RING_pattern_replacement[hr_pattern])
      specializableParams["[[%REPID]]"] = HALF_RING_pattern_replacement[hr_pattern]
      suffix_Type = suffix_Type[len(hr_pattern):] # Update suffix var to keep only the remaining str (if any)
      break # since only one match allowed

  # 2b- Check for optional Detset suffix ("_detsetX" where X is a number)
  # Note that default is DEFAULT_DETSET. If "" (=full) all available bolos are
  #      used, otherwise only those in the specified detset.
  currentDetset = DEFAULT_DETSET
  for detset_pattern in ("_detset1", "_detset2", "_detset3"):
    if suffix_Type.startswith(detset_pattern):
      logDebug("***** DETECT detset config!!!")
      log("  -> Input parameter will be proceceed using '%s'" % detset_pattern)
      currentDetset = detset_pattern[1:]
      suffix_Type = suffix_Type[len(detset_pattern):] # Update suffix var to keep only the remaining str (if any)
      break # since only one match allowed

  # Ensure that the suffixType was ok (totaly fully parsed), otherwise error
  if len(suffix_Type) > 0:
    print("ERROR: Invalid <type> option: '%s'" % suffix_Type)
    exit(1)

  # 2b- Parameter "sub type" (ex: "SWB", "SWB_HALFDET2" or "")
  currentParamSubType = sys.argv[argpos]
  argpos+=1
  logDebug(">> Read from command line: parameter sub type = %s" % currentParamSubType)
  if currentParamSubType != "":
    currentParamSubType = "_" + currentParamSubType
  flag_SubTypeSWB = False
  if "SWB" in currentParamSubType:
    flag_SubTypeSWB = True

  # Construct the corresponding parameter filename
  #paramFilename = "Param_%s%s%s__GEN_v2__.txt" % (currentBoloFreq[:3], currentProcessTypeLabel, currentParamSubType)
  paramFilename = "Param_%s%s__GEN_v2__.txt" % (currentBoloFreq[:3], currentProcessTypeLabel)
  log("  -> input parameter file: '%s'" % paramFilename)

  # 3- Update default value accordingly to current freq
  # For the special case of 353 we must set NADU and GAINSTEP to 32 (instead of the default 128)
  # !!! XXX NO MORE TRUE !!!
#  if "353" in currentBoloFreq:
#    specializableParams["[[%NADU]]"] = "32"
#    specializableParams["[[%GAINSTEP]]"] = "32"

  # Check if it exist some specific param for this frequency
  if currentBoloFreq[:3] in specializableParams_SPE.keys():
    # Therefore we simply update default dict (and if required, create new entry in the final dictionary)
    specializableParams.update(specializableParams_SPE[currentBoloFreq[:3]])

  # 4- Number of proc to be used
  nbTotProc = int(sys.argv[argpos])
  argpos+=1
  logDebug(">> Read from command line: nbTotProc = %d" % nbTotProc)

  assert(nbTotProc > 0)

  # Compute nb of nodes required:
  nb_proc_per_node = currentHost.nbProcPerNode
  nbNodeRequired = nbTotProc / nb_proc_per_node # integer division
  if (nbTotProc % nb_proc_per_node) != 0:
    nbNodeRequired += 1

  log("  -> nbNodeRequired = %d" % nbNodeRequired)

  # 5- (optional) Walltime
  if argpos < len(sys.argv) and "=" not in sys.argv[argpos]: # Ensure that walltime is not confused with
                                  # any folowing specialized param
    currentWalltime = sys.argv[argpos]
    argpos+=1
  else:
    currentWalltime = DEFAULT_WALLTIME
  log("  -> walltime : %s" % currentWalltime)

  # 6- (optional) Parse any specialized parameter.
  while argpos < len(sys.argv):
    tmpArg = sys.argv[argpos]
    logDebug("tmpArg = %s" % tmpArg)

    # Ensure validity of argument
    if "=" not in tmpArg:
      print("Error: argument '%s' is not of the form key=value")
      exit(1)

    # Extract key and value
    tmpKey,tmpValue = tmpArg.split("=")

    # Store it for further use
    specializableParams["[[%%%s]]"%tmpKey] = tmpValue # Note that we store the key directly in its form of "[[%key]]"
    logDebug("  -> specialized param detected: %s=%s" % (tmpKey,tmpValue))
    argpos+=1

    # Check since not yet supported!
    if tmpKey == "REPID" and suffix_HR_found:
      print("Changing the REPid while requesting HR's is NOT YET SUPPORTED in this script, please update!")
      exit(0)

  ### Finally make some assumption to help user avoid errors...

  # In case of 353 ensure the user make use of correct values for NADU and GAINSTEP
  # !!! XXX NO MORE TRUE !!!
#  if "353" in currentBoloFreq:
#    #logDebug("NADU = %s" % specializableParams["[[%NADU]"])
#    if (specializableParams["[[%NADU]]"] != "32") or (specializableParams["[[%GAINSTEP]]"] != "32"):
#      print("***************")
#      print("*** Detection of a possible error! (see detail below)")
#      print("*** For %s, GAINSTEP and NADU should be set to 32 (current values are respectively %s and %s)" % (currentBoloFreq, specializableParams["[[%GAINSTEP]]"], specializableParams["[[%NADU]]"]))
#      print("***************")
#      if not query_yes_no("Are you sure you want to continue?", default="no"):
#        exit(0)


  ### Prepare env

  # Create directories (if required)
  # a) base dir
  outDir = currentHost.srollOut_BaseDir
  # Ensure the directory exist (create it if required)
  if not os.path.exists(outDir):
    os.makedirs(outDir)
  else:
    print("!!!!!!!!!!!!!!")
    print("!!! WARNING: directory '%s' already exist, DATA WILL BE OVERWRITTEN!" % outDir)
    print("!!!!!!!!!!!!!!")
    if not query_yes_no("Are you sure you want to continue?", default="no"):
      exit(0)
  # b) MAP dir
  outDir = currentHost.srollOut_MapDir
  if not os.path.exists(outDir):
    os.makedirs(outDir)
  # c) VECT dir
  outDir = currentHost.srollOut_VectDir
  if not os.path.exists(outDir):
    os.makedirs(outDir)


  ### Retrieve CVS sroll.c revision and save it to output dir
  # WARNING: this is the cvs revision at qsub launch time not running time!! MAY DIFFER!!!

  # Path to sroll.c file:
  fullPathToSrollSrcFile = currentHost.srollIn + "/sroll.c"

  import subprocess
  cvsFileREV = subprocess.check_output(["cvs", "status", "../sroll/sroll.c"])
  #print("***** cvsFileREV = %s" % cvsFileREV)

  # Save srollREV file
  temp_filename = currentHost.srollOut_BaseDir+"/srollREV"
  with open(temp_filename, "w") as temp_file_object:
    temp_file_object.writelines(cvsFileREV)
  logDebug("Write cvsFileREV file: '%s'" % temp_filename)


  ### Pre-Process and Specializations of the parameter file....

  # 1) Pre-Process
  currentParamFileContent = preProcessParameterFile(currentHost.srollParamDir+paramFilename, currentBoloFreq[:3], currentDetset, flag_SubTypeSWB)

#  print("")
#  print("")
#  print("")
#  print("=============================================")
#  print(currentParamFileContent)
#  print("=============================================")
#  print("")
#  print("")
#  print("")


  # 2) Specialization

  # Make a log to user about specializableParams
  log("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
  log("%%%%%%% Specializable Params %%%%%%%")
  log("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
  for spe_param,spe_param_val in specializableParams.iteritems():
    log("%s = %s" % (spe_param, spe_param_val))
  log("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")

  # Specialize parameter file (according to host configuration)
  currentParamFile = specializeParameterFile(currentParamFileContent, currentHost, specializableParams)

  # Specialize Qsub and launch it

  # Process label
  label = "sroll_%s%s" % (currentBoloFreq, currentProcessTypeLabel)

  # Log file
  logFilename = "DBG_sroll_%s%s_%d_proc.log" % (currentBoloFreq, currentProcessTypeLabel, nbTotProc)
  log("  -> logFilename: '%s'" % logFilename)

  #DEBUG XXX
  #exit(1)

  # Launch process with appropriate values
  specializeQsubAndLaunchIt(currentHost, label, nbNodeRequired, nbTotProc, logFilename, currentParamFile, currentWalltime)

  # Final end
  print("")
  print("All done.")
