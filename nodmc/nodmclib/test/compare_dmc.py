##############################################################################
# compare_dmc.py
# This python script aims to compare the former legacy DMC with result obtain
# using the new NO_DMC_LIB (result shall be identical).
# Additionnaly the script allow to display performance result.
#
# NOTE: some test suite has been commented out (preventing failed status)
#       since the behaviour of the former DMC was to restrict access on 'MAP'
#       data when request is outside valid range. This behaviour is NO MORE
#       TRUE for the new NO_DMC_LIB!
#
# Author  : Christian Madsen
# Date    : 2015-01-27
# version : Initial
##############################################################################

import time
from optparse import OptionParser
import numpy as np


### Configuring the test suite

class TestSuite:
  """Class representing a test suite"""
  pass

# PIODOUBLE MAP TestSuite
T_double_MAP = TestSuite()
T_double_MAP.label = "Test on PIODOUBLE MAP object (0 -> 50331647)"
T_double_MAP.dataType = "PIODOUBLE"
T_double_MAP.obj = "/data/dmc/MISS03/DATA/PSM190_FFP8_MAP/100-1a_fg_map_I"
T_double_MAP.dsc = "This is a MAP object"
T_double_MAP.samples = [[0, 10, "Read the first ten samples"],
                        [2048, 2, "Read at a cross chunk offset"],
                        #SEE NOTE: [50331647, 10, "Read inside then outside range"],
                        [0, 50331648, "Read the whole product"]]

# PIODOUBLE ROI TestSuite
T_double_ROI = TestSuite()
T_double_ROI.label = "Test on PIODOUBLE ROI object (0 -> 46417)"
T_double_ROI.dataType = "PIODOUBLE"
T_double_ROI.obj = "/data/dmc/MISS03/DATA/SREM_ROI/TC1_countrate"
T_double_ROI.samples = [[0, 20, "Read the first values"],
                        [46410, 30, "Read inside then outside values"],
                        [0, 46418, "Read whole values"],
                        [0, 50000, "Read whole values and after"],
                        [0, 1000000, "Read huge amount of values (inside then outside)"]]

# PIODOUBLE ROI TestSuite _2
T_double_ROI_2 = TestSuite()
T_double_ROI_2.label = "Test on PIODOUBLE ROI object (_v2) (13 -> 46416)"
T_double_ROI_2.dataType = "PIODOUBLE"
T_double_ROI_2.obj = "/data/dmc/MISS03/DATA/RingInfo/SpinPhase"
T_double_ROI_2.samples = [[7, 6, "Read outside values (before)"],
                          [5, 15, "Read outside then inside values"],
                          [13, 11, "Read the first values"],
                          [46410, 20, "Read inside then outside values"],
                          [13, 46404, "Read whole values"],
                          [15000, 20000, "Read inside values"],
                          [0, 50000, "Read before and after"],
                          [0, 1000000, "Read huge amount of values (including all valid values)"]]

# PIODOUBLE ROI TestSuite _3
T_double_ROI_3 = TestSuite()
T_double_ROI_3.label = "Test on PIODOUBLE ROI object (_v3) (4005 -> 28486)"
T_double_ROI_3.dataType = "PIODOUBLE"
T_double_ROI_3.obj = "/data/dmc/MISS03/DATA/calROIs_SUB90MIN/11_217_5a_jumpAmplitude_bigPlanets_v61"
T_double_ROI_3.samples = [[3000, 1000, "Read outside values (before)"],
                          [4000, 15, "Read outside then inside values"],
                          [4005, 20, "Read the first values"],
                          [28480, 20, "Read inside then outside values"],
                          [4005, 24482, "Read whole values"],
                          [20000, 7000, "Read inside values"],
                          [1, 50000, "Read before and after"],
                          [3, 1000000, "Read huge amount of values (including all valid values)"]]

# PIOFLOAT MAP TestSuite
T_float_MAP = TestSuite()
T_float_MAP.label = "Test on PIOFLOAT MAP object (0 -> 50331647)"
T_float_MAP.dataType = "PIOFLOAT"
T_float_MAP.obj = "/data/dmc/MISS03/DATA/DR2_MAP_GALACTIC_2048/mask_galaxy_0.80_apodised"
T_float_MAP.samples = [[0, 10000, "Read the first valid values"],
                       [25600000, 136000, "Read part of valid values (more than one chunk)"],
                       [123456, 987456, "Read huge amount of valid values"],
                       [0, 50331648, "Read the whole product"],
                       #SEE NOTE: [50331547, 120, "Read last data values and over"],
                       #SEE NOTE: [50331648, 50, "Read only outside values"]
                       ]

# PIOFLOAT MAP TestSuite _2
T_float_MAP_2 = TestSuite()
T_float_MAP_2.label = "Test on PIOFLOAT MAP object (_v2) (0 -> 12582911)"
T_float_MAP_2.dataType = "PIOFLOAT"
T_float_MAP_2.obj = "/data/dmc/MISS03/DATA/MAP1024_LFI_DX7/TEMPERATURE_44"
T_float_MAP_2.samples = [[0, 64, "Read the first valid values"],
                         [16984, 16984, "Read part of valid values"],
                         [10000000, 2582911, "Read huge amount of values"],
                         #SEE NOTE: [12581911, 2000, "Read last data values and over"],
                         #SEE NOTE: [99999999, 100000, "Read only outside values"]
                         ]

# PIOFLOAT MAP TestSuite _3
T_float_MAP_3 = TestSuite()
T_float_MAP_3.label = "Test on PIOFLOAT MAP object (_v3) (0 -> 786431)"
T_float_MAP_3.dataType = "PIOFLOAT"
T_float_MAP_3.obj = "/data/dmc/MISS03/DATA/MAP0256_LFI_DX10/LFI_SkyMap-bandpassCorrection_070_0256_DX10_full_I"
T_float_MAP_3.samples = [[0, 64, "Read the first valid values"],
                         [86431, 10000, "Read part of valid values"],
                         [0, 786432, "Read whole values"],
                         #SEE NOTE: [786430, 10, "Read last data values and over"],
                         #SEE NOTE: [786432, 124, "Read only outside values"]
                         ]

# PIOFLOAT TOI TestSuite
T_float_TOI = TestSuite()
T_float_TOI.label = "Test on PIOFLOAT TOI object"
T_float_TOI.dataType = "PIOFLOAT"
T_float_TOI.obj = "/data/dmc/MISS03/DATA/e2eJan15_TOI/143-1a_00"
T_float_TOI.samples = [[0, 10, "Read the first ten samples (that does not exist)"],
                       [1371347725, 10, "Read valid range values (first 10)"],
                       [15147406766, 5, "Read last data values and over offset"]]

# PIOFLOAT TOI TestSuite _2
T_float_TOI_2 = TestSuite()
T_float_TOI_2.label = "Test on PIOFLOAT TOI object (_v2) (2564511940 -> 12631815970)"
T_float_TOI_2.dataType = "PIOFLOAT"
T_float_TOI_2.obj = "/data/dmc/MISS03/DATA/PlanetSimu640_TOI/100-1a_background"
T_float_TOI_2.samples = [[0, 10, "Read the first ten (outside product def)"],
                         [2564511930, 20, "Read outside (before) and inside product def"],
                         [12631815960, 20, "Read inside product def and outside (after)"],
                         [12631815971, 20, "Read outside (after)"],
                         [2564511940, 20, "Read first valid values"],
                         #Too huge memory required! [2564511940, 10067304031, "Read whole set of valid values"],
                         [2564554750, 10, "Cross chunk (valid values)"],
                         [2564511940, 100001, "Read a quite huge amount of valid values (from beggining)"],
                         [3564572201, 50000000, "Read a huge amount of valid values (inside product) but some missing chunk!"],
                         [3564511940, 100000001, "Read a huge amount of valid values (inside product)"],
                         [3564511941, 100000001, "Read a huge amount of valid values (inside product). Slice from 1 comparing to previous test"],
                         [2564511940, 500000008, "Read a huge amount of valid values"]]
# PIOSTRING VECT TestSuite
T_string_VECT = TestSuite()
T_string_VECT.label = "Test on PIOSTRING object"
T_string_VECT.dataType = "PIOSTRING"
T_string_VECT.obj = "/data/dmc/MISS03/DATA/ECC_v44/NAME"
T_string_VECT.samples = [[0, 10, "Read the first ten samples"],
                         [910, 120, "Read plain values then empty values in chunk then outside chunk info"],
                         [0, 915, "Read whole valid data"]]

# PIOBYTE VECT TestSuite
T_byte_VECT = TestSuite()
T_byte_VECT.label = "Test on PIOBYTE VECT object"
T_byte_VECT.dataType = "PIOBYTE"
T_byte_VECT.obj = "/data/dmc/MISS03/DATA/IMO_DATA/MB_LFI18_TAB_Y_FM_1-0_commit-2251799814551416"
T_byte_VECT.samples = [[0, 10, "Read the first ten samples"],
                       [6614350, 120, "Read 6 valid values then outside values"],
                       [0, 6614355, "Read whole valid data"]]

# PIOFLAG VECT TestSuite
T_flag_VECT = TestSuite()
T_flag_VECT.label = "Test on PIOFLAG VECT object"
T_flag_VECT.dataType = "PIOFLAG"
T_flag_VECT.obj = "/data/dmc/MISS03/DATA/IMO_DATA/HFI__bc00_100_1a_Apod5_v202_Flag_Y_commit-34432597"
T_flag_VECT.samples = [[0, 10, "Read the first ten samples"],
                       [12309, 120, "Read 6 valid values then outside values"],
                       [0, 12314, "Read whole valid data"]]

# PIOFLAG ROI TestSuite
T_flag_ROI = TestSuite()
T_flag_ROI.label = "Test on PIOFLAG ROI object (empty object!)"
T_flag_ROI.dataType = "PIOFLAG"
T_flag_ROI.obj = "/data/dmc/MISS03/DATA/cleanROIs/KeepMe"
T_flag_ROI.samples = [[0, 10, "Read outside values"],
                      [309, 111120, "Read outside values"]]

# PIOFLAG TOI TestSuite
T_flag_TOI = TestSuite()
T_flag_TOI.label = "Test on PIOFLAG TOI object"
T_flag_TOI.dataType = "PIOFLAG"
T_flag_TOI.obj = "/data/dmc/MISS03/DATA/calTOIs/84_217_4_Flag_LFER6_JC_bigPlanets_v67"
T_flag_TOI.samples = [[2452095789, 10, "Read first 10 values"],
                      [2452095788, 123, "Read 1 before then inside values"],
                      [2452095789, 1000000001, "Read huge amount inside values (starting from the first valid value)"],
                      [15147406766, 1000000001, "Read huge amount values (last valid then outside)"]]

# PIOINT TOI TestSuite
T_int_TOI = TestSuite()
T_int_TOI.label = "Test on PIOINT TOI object (0 -> 25195079333)"
T_int_TOI.dataType = "PIOINT"
T_int_TOI.obj = "/data/dmc/MISS03/DATA/Sa_HFI_C_Bolo/HFI_00_C"
T_int_TOI.samples = [[0, 10, "Read the first ten samples"],
                     [5079333, 100, "Read valid range values"],
                     [25195079333, 5, "Read last data values and over offset"]]

# PIOINT ROI TestSuite
T_int_ROI = TestSuite()
T_int_ROI.label = "Test on PIOINT ROI object (13 -> 46416)"
T_int_ROI.dataType = "PIOINT"
T_int_ROI.obj = "/data/dmc/MISS03/DATA/RingInfo/pointingId"
T_int_ROI.samples = [[1, 10, "Read outside values"],
                     [10, 30, "Read outside then inside values"],
                     [13, 46404, "Read whole values"],
                     [5, 123456789, "Read over all values"]]

# PIOINT ROI TestSuite _2
T_int_ROI_2 = TestSuite()
T_int_ROI_2.label = "Test on PIOINT ROI object (_v2) (240->24305)"
T_int_ROI_2.dataType = "PIOINT"
T_int_ROI_2.obj = "/data/dmc/MISS03/DATA/calROIs/NumNull_62_217_7b_UnvalidDataFlag_v50"
T_int_ROI_2.samples = [[200, 39, "Read outside values"],
                       [150, 150, "Read outside then inside values"],
                       [240, 24066, "Read whole values"],
                       [239, 222222222, "Read over all values"]]

# PIOINT ROI TestSuite _3
T_int_ROI_3 = TestSuite()
T_int_ROI_3.label = "Test on PIOINT ROI object (_v3) (240->27005)"
T_int_ROI_3.dataType = "PIOINT"
T_int_ROI_3.obj = "/data/dmc/MISS03/DATA/ROI_JMD/30_143_2a_discarded_rings_H2"
T_int_ROI_3.samples = [[200, 1, "Read outside values"],
                       [0, 27000, "Read outside then inside values"],
                       [240, 26766, "Read whole values"],
                       [120, 30000, "Read over all values"]]

# PIOLONG ROI TestSuite
T_long_ROI = TestSuite()
T_long_ROI.label = "Test on PIOLONG ROI object (13 -> 46416)"
T_long_ROI.dataType = "PIOLONG"
T_long_ROI.obj = "/data/dmc/MISS03/DATA/RingInfo/RingStartTime"
T_long_ROI.samples = [[3, 11, "Read outside values"],
                     [12, 10000, "Read outside then inside values"],
                     [13, 46404, "Read whole values"],
                     [2, 50000, "Read over all values"]]

# PIOLONG ROI TestSuite _2
T_long_ROI_2 = TestSuite()
T_long_ROI_2.label = "Test on PIOLONG ROI object (_v2) (4005 -> 28486)"
T_long_ROI_2.dataType = "PIOLONG"
T_long_ROI_2.obj = "/data/dmc/MISS03/DATA/calROIs_SUB90MIN/11_217_5a_jumpRing_bigPlanets_v61"
T_long_ROI_2.samples = [[300, 378, "Read outside values"],
                     [28100, 500, "Read inside then outside values"],
                     [4005, 24482, "Read whole values"],
                     [3000, 400000, "Read over all values"]]

# PIOLONG ROI TestSuite _3
T_long_ROI_3 = TestSuite()
T_long_ROI_3.label = "Test on PIOLONG ROI object (_v3) (240 -> 27005)"
T_long_ROI_3.dataType = "PIOLONG"
T_long_ROI_3.obj = "/data/dmc/MISS03/DATA/calROIs/41_100_3b_statNb_valid_WTD_v72"
T_long_ROI_3.samples = [[189, 27020, "Read over all values"]]


# Gather all testsuites in the list
testSuiteList = [T_double_MAP, T_double_ROI, T_double_ROI_2, T_double_ROI_3,
                 T_float_MAP, T_float_MAP_2, T_float_MAP_3, T_float_TOI, T_float_TOI_2,
                 T_string_VECT,
                 T_byte_VECT,
                 T_flag_VECT, T_flag_ROI, T_flag_TOI,
                 T_int_TOI, T_int_ROI, T_int_ROI_2, T_int_ROI_3,
                 T_long_ROI, T_long_ROI_2, T_long_ROI_3]


#-----------------------------------------------------------------------------

# This function is responsible to check that metadat extract from the legacy
# DMC (see function pio.InfoObject()) and the one of new library
# (see nodmclib.getObjectInfo()) are identical.
#
# Note that the new metadata function returns 2 additional info: 'BeginRing'
# and 'EndRing'.
def checkMetadata(obj):
  print("Checking Metadata...")
  # a) Retrieve meta using legacy DMC
  meta_ref = pio.InfoObject(obj)

  # b) Retrieve meta using new lib
  meta = nodmclib.getObjectInfo(obj)

  if options.FLAG_printExtraDebug:
    print("meta_ref = ", meta_ref)
    print("meta     = ", meta)

  if meta_ref[0] == meta[0] and \
     meta_ref[1] == meta[1] and \
     meta_ref[2] == meta[2] and \
     meta_ref[3] == meta[3] and \
     meta_ref[4] == meta[6] and \
     meta_ref[5] == meta[7]:
    print("--> Metadata OK")
    return True
  else:
    print("--> Metadata FAILED!!!!!!!!!")
    return False

#-----------------------------------------------------------------------------

def parseCmdline():
  # Init default options
  ts_pos = -1 # No test suite specified by default
  ts_samp_pos = -1 # No test suite sample specified by default

  #usage = "usage: %prog [options] arg"
  #parser = OptionParser(usage)
  parser = OptionParser()
  parser.set_defaults(ts_pos=ts_pos, ts_samp_pos=ts_samp_pos)
  parser.add_option("-d", "--debug", action="store_true", dest="FLAG_printExtraDebug",
                    help="Print extra debug info (can print huge amount of text!)")
  parser.add_option("-t", "--test", dest="t_select",
                    help="Select only one specified test suite")
  parser.add_option("-m", "--metadata", action="store_true", dest="FLAG_only_metadata",
                    help="Allow user to proceed only metadata check")

  (options, args) = parser.parse_args()
#    if len(args) != 1:
#        parser.error("incorrect number of arguments")
#    if options.verbose:
#        print "reading %s..." % options.filename

  if options.FLAG_printExtraDebug:
    print("[User option]: DEBUG mode activated!")
  if options.t_select:
    if "-" in options.t_select:
      options.ts_pos, options.ts_samp_pos = [int(x) for x in options.t_select.split('-')]
      print("[User option]: AS REQUESTED, ONLY ONE TEST SUITE SAMPLE WILL BE PROCEED: #%d-%d" % (options.ts_pos, options.ts_samp_pos))
    else:
      options.ts_pos = int(options.t_select)
      print("[User option]: AS REQUESTED, ONLY ONE FULL TEST SUITE WILL BE PROCEED: #%d" % options.ts_pos)
  if options.FLAG_only_metadata:
    print("[User option]: AS REQUESTED, ONLY METADATA WILL BE PROCEED")

  print("")

  return options

#-----------------------------------------------------------------------------

def printDiffBetweenArrays(tab_ref, tab):
  # If size diff we do not go further.
  tmpSizeDiff = tab_ref.size - tab.size
  if (tmpSizeDiff != 0):
    print("Result arrays are not of the same size!!! : size_ref=%d size=%d" % (tab_ref.size, tab.size))
    return

  # Find the diff between the 2 arrays
  itemindex = np.where(~((tab_ref == tab) | (np.isnan(tab_ref) & np.isnan(tab))))

  print("Number of diff: %d" % len(itemindex[0]))

  print("#Index #Ref #New (only first 10!)")
  for index in itemindex[0][:10]:
    print("%d: %r %r" % (index, tab_ref[index], tab[index]))

#-----------------------------------------------------------------------------

def isArraysEquals(tab_ref, tab):
  """Return True if arrays are of the same size and each element has the same values (note that np.nan in this case is considered to be equal values). Return False otherwise.
  Note: there is a special test for string array since the beahve quite differently.
  """
  if tab_ref.dtype.kind == 'S':
    return (len(tab_ref) == len(tab)) and (tab_ref == tab).all()
  else:
    return (len(tab_ref) == len(tab)) and (((tab_ref == tab) | (np.isnan(tab_ref) & np.isnan(tab))).all())

#-----------------------------------------------------------------------------

def printPerf(i, j, typeDataTest, dataType, nbSample, duration_ref, duration):
  """Allow to print the performance result.
  Parameters:
  ----------
  i the test suite index number (starting from 1)
  j the sample index in test suite (starting from 1)
  typeDataTest the type of test data (flag or data)
  duration_ref the reference duration (legacy DMC)
  duration the new duration obtain with new lib code
  """

  s_testId = "#%d-%d" % (i, j)
  s_dataType_NbSample = "[%s - %d]" % (dataType, nbSample)
  print("PERFORMANCE %-6s (%s) %-25s : legacy DMC = %.2fs    no DMC lib = %.2fs    GAIN: %.2f" % (s_testId, typeDataTest, s_dataType_NbSample, duration_ref, duration, duration_ref/duration))

#-----------------------------------------------------------------------------


##############################################################################
# MAIN
# This code allow to execute the module as a standalone script.
##############################################################################
if __name__ == "__main__":
  import sys
  import piolib as pio
  import nodmclib

  print("")
  print("************************************************")
  print("* %s" % sys.argv[0])
  print("* Goal: compare legacy DMC result with the new NO_DMC_LIB")
  print("* Total test suite available: %d" % len(testSuiteList))
  print("************************************************")
  print("")


  # Parse command line options
  options = parseCmdline()

  # init count var
  nbSuccess = 0
  nbFailed = 0
  nbMetadataError = 0
  testFailedList = []

  # Loop on all test suite
  for i, ts in enumerate(testSuiteList, 1):
    # Allow to skip test suite if required by user
    if options.ts_pos != -1 and options.ts_pos != i:
      continue

    print("")
    print("================================================")
    print(">> Test suite #%d: \"%s\"" % (i, ts.label))
    print(">>   obj: %s" % ts.obj)
    print("================================================")

    # Check metadata
    if not checkMetadata(ts.obj):
      nbMetadataError += 1
      testFailedList.append("#%d (METADATA)" % i)
    if options.FLAG_only_metadata:
      continue
    
    # Loop on all samples to be tested for the current test suite
    for j, sample in enumerate(ts.samples, 1):
      # Allow to skip test suite sample if required by user
      if options.ts_samp_pos != -1 and options.ts_samp_pos != j:
        continue

      print("")
      print("-------------------------------------------------")
      print("Sample #%d-%d to be tested: offset = %d   nbSample = %d" % (i, j, sample[0], sample[1]))
      print("Description: \"%s\"" % sample[2])
      print("-------------------------------------------------")

      print("** DATA **")
      print("- Reading with legacy DMC...")
      cmd = "begin=%d;end=%d" % (sample[0], sample[0]+sample[1]-1)
      print("  using '%s'" % cmd)
      t_ref_start = time.time()
      data_ref = pio.read(ts.obj, command="%s" % cmd)
      t_ref_stop = time.time()

      print("- Reading with NO DMC...")
      func_name = "read_%s" % ts.dataType
      print("  using %s()" % func_name)
      # Call the corresponding "read" function from module "nodmclib"
      t_start = time.time()
      data = getattr(nodmclib, func_name)(ts.obj, sample[0], sample[1])
      t_stop = time.time()

      if options.FLAG_printExtraDebug:
        print("data_ref = %s" % data_ref)
        print("data     = %s" % data)

      duration_ref = t_ref_stop - t_ref_start
      duration = t_stop - t_start
      printPerf(i, j, "DATA", ts.dataType, sample[1], duration_ref, duration)

      # Print result of comparison
      print("_________________________________________________")
      # Handle comparision including case of nan
      if isArraysEquals(data_ref, data):
        print("---> SUCCESS!")
        nbSuccess += 1
      else:
        if options.FLAG_printExtraDebug:
          printDiffBetweenArrays(data_ref, data)

        print("---> !!!!! FAILED !!!!!")
        testFailedList.append("#%d-%d (DATA) : '%s'" % (i,j, sample[2]))
        nbFailed += 1

      # Try to free mem (since data is no more required)
      del data_ref
      del data

      print("")
      print("** FLAG **")
      print("- Reading flag with legacy DMC...")
      flagname = "Written"
      print("  using: flgname='%s' command='%s'" % (flagname, cmd))
      t_ref_start = time.time()
      data_flag_ref = pio.readFlg(ts.obj, flgname="%s" % flagname, command="%s" % cmd)
      t_ref_stop = time.time()

      print("- Reading flag with NO DMC...")
      t_start = time.time()
      data_flag = nodmclib.readFlag_Written(ts.obj, sample[0], sample[1])
      t_stop = time.time()

      if options.FLAG_printExtraDebug:
        print("data_flag_ref = %s" % data_flag_ref[1])
        print("data_flag     = %s" % data_flag)

      duration_ref = t_ref_stop - t_ref_start
      duration = t_stop - t_start
      printPerf(i, j, "FLAG", ts.dataType, sample[1], duration_ref, duration)

      # Print result of comparison
      print("_________________________________________________")
      if isArraysEquals(data_flag_ref[1], data_flag):
        print("---> SUCCESS!")
        nbSuccess += 1
      else:
        if options.FLAG_printExtraDebug:
          printDiffBetweenArrays(data_flag_ref[1], data_flag)

        print("---> !!!!! FAILED !!!!!")
        testFailedList.append("#%d-%d (FLAG) : '%s'" % (i,j, sample[2]))
        nbFailed += 1

      # Try to free mem (since data is no more required)
      del data_flag_ref
      del data_flag

  print("")
  print("")
  print("")
  print("************************************************")
  print("Final count:  nbSuccess=%d   nbFailed=%d" % (nbSuccess, nbFailed))
  print("              nbMetadataError = %d" % (nbMetadataError))
  if nbFailed == 0 and nbMetadataError == 0:
    print(">>> SUCCESS <<<")
  else:
    print(">>> FAILED <<<")
    print("Test failed list:")
    for tfail in testFailedList:
      print("  %s" % tfail)
  print("************************************************")
