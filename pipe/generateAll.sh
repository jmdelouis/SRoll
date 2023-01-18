#!/bin/bash

###############################################################################
# generateAll.sh
#
# This script allow to quickly launch several jobs (that rely on doQsub.py)
# using a quite simple flag configuration.
#
# Author  : Christian Madsen
# Date    : 2015-11-18
# version : Initial
###############################################################################

NB_PROC=256 # default is 512

# Prefix and suffix in the output name dir
PREFIX=""
SUFFIX=""
repid="REP6" # "REP6" (default)

walltime="15:00:00" # "10:00:00" (default)

# Flags are working by using the values: 1 for true, 0 for false

flag_SWB=1 # For the SWB mode.
# Note that SWB flag has no effet for 100GHz.
# Furthermore in case of detset (non empty) this flag is temporary force to 0!
flag_DETSET=0 # If set to true activate the production of detset1 and detset2 in addition to the full map.
flag_HM=0 # For the HalfMission
flag_FM=0 # For the FullMission
flag_PARITY=0 # For additionnal map with Odd and Even
flag_YEAR=0 # For additionnal map for Year1 and Year2 
flag_HR=0 # For additionnal map HalfRing

# Static list of freq not allowing SWB
FREQ_NO_SWB=(100 545 857)

# The NEW definition of Half Mission regarding ring index, as:
#   HM1: 240 -> HALF_MISSION_RING_INDEX 
#   HM2: (HALF_MISSION_RING_INDEX+1) -> 26050
HALF_MISSION_RING_INDEX=13144

# List of freq to be proceed
#freq_list=(100 143 217 353 545 857)
freq_list=(100 143 217 353)


# Prepare label and cmd regarding REP
if [ "${repid}" == "REP6" ]
then
  REP_LABEL=
  REP_CMD=
else
  REP_LABEL="_${repid}"
  REP_CMD="REPID=${repid}"
fi


# Loop on all required freq
for freq in "${freq_list[@]}"
do
  # Check if current freq support SWB
  if [[ " ${FREQ_NO_SWB[@]} " =~ " ${freq} " ]]
    then
      tmp_no_SWB=1
    else
      tmp_no_SWB=0
  fi

  # Prepare label and cmd regarding SWB
  if [ "${tmp_no_SWB}" == 1 ] || [ ${flag_SWB} -eq 0 ]
  then
    SWB_LABEL=""
    SWB_CMD=""
  else
    SWB_LABEL="_SWB"
    SWB_CMD="SWB"
  fi

  #============================================================================
  # Case of the HRs
  #============================================================================

  if [ ${flag_HR} -eq 1 ]
    then
      for HR in "HR1" "HR2"
      do
        #echo "Launching qsub job for ${freq} ${HR} with SWB_LABEL=${SWB_LABEL} SWB_CMD=${SWB_CMD}"
        python doQsub.py ${PREFIX}${freq}${SWB_LABEL}_${HR}${REP_LABEL}${SUFFIX} ${freq}GHz REAL_${HR} "${SWB_CMD}" ${NB_PROC} ${walltime} ${REP_CMD}
      done
  fi

#  # Case of the HRs
#  for HR in "HR1" "HR2"
#  do
#    echo "Launching qsub job for ${freq} ${HR} with SWB_LABEL=${SWB_LABEL} SWB_CMD=${SWB_CMD}"
#    python doQsub.py ${PREFIX}${freq}${SWB_LABEL}_HM1_${HR}${REP_LABEL}${SUFFIX} ${freq}GHz REAL_${HR} "${SWB_CMD}" ${NB_PROC} ${walltime} BeginRing=240 EndRing=${HALF_MISSION_RING_INDEX} ${REP_CMD}
#    python doQsub.py ${PREFIX}${freq}${SWB_LABEL}_HM2_${HR}${REP_LABEL}${SUFFIX} ${freq}GHz REAL_${HR} "${SWB_CMD}" ${NB_PROC} ${walltime} BeginRing=$((HALF_MISSION_RING_INDEX+1)) EndRing=26050 ${REP_CMD}
#  done


  #============================================================================
  # Case of Detset
  #============================================================================

  if [ ${flag_DETSET} -eq 1 ]
    then

      saved_SWB_LABEL=${SWB_LABEL}
      saved_SWB_CMD=${SWB_CMD}

      for detset in "_detset1" "_detset2"
      do
        #echo "**DETSET = '${detset}'"

        # Set SWB config to none in case of DETSET (detset1 or detset2)!!!
        if [ -n "${detset}" ]
          then
            SWB_LABEL=""
            SWB_CMD=""
        fi

        python doQsub.py ${PREFIX}${freq}${SWB_LABEL}${detset}${REP_LABEL}${SUFFIX} ${freq}GHz REAL${detset} "${SWB_CMD}" ${NB_PROC} ${walltime} ${REP_CMD}

        # DETSET end ----- Set back the SWB config
        SWB_LABEL=${saved_SWB_LABEL}
        SWB_CMD=${saved_SWB_CMD}

      done
  fi

  #============================================================================
  # Case of FM
  #============================================================================
  # Case of standard *FullMission* (optionnaly with SWB)
  if [ ${flag_FM} -eq 1 ]
    then
      python doQsub.py ${PREFIX}${freq}${SWB_LABEL}${REP_LABEL}${SUFFIX} ${freq}GHz REAL "${SWB_CMD}" ${NB_PROC} ${walltime} ${REP_CMD}
  fi

  #============================================================================
  # Case of HM
  #============================================================================

  # Half Mission - Full Ring - SWB (except for 100GHz and detset)
  if [ ${flag_HM} -eq 1 ]
    then
      python doQsub.py ${PREFIX}${freq}${SWB_LABEL}_HM1${REP_LABEL}${SUFFIX} ${freq}GHz REAL "${SWB_CMD}" ${NB_PROC} ${walltime} BeginRing=240 EndRing=${HALF_MISSION_RING_INDEX} ${REP_CMD}
      python doQsub.py ${PREFIX}${freq}${SWB_LABEL}_HM2${REP_LABEL}${SUFFIX} ${freq}GHz REAL "${SWB_CMD}" ${NB_PROC} ${walltime} BeginRing=$((HALF_MISSION_RING_INDEX+1)) EndRing=26050 ${REP_CMD}
  fi

  #============================================================================
  # Case of Odd/Even
  #============================================================================

  if [ ${flag_PARITY} -eq 1 ]
    then
      for parity in "_Odd" "_Even"
      do
        python doQsub.py ${PREFIX}${freq}${SWB_LABEL}${parity}${REP_LABEL}${SUFFIX} ${freq}GHz REAL "${SWB_CMD}" ${NB_PROC} ${walltime} ${REP_CMD} OddEven=${parity}
      done
  fi

  #============================================================================
  # Case of YR1 / YR2
  #============================================================================
  if [ ${flag_YEAR} -eq 1 ]
    then
      python doQsub.py ${PREFIX}${freq}${SWB_LABEL}_YR1${REP_LABEL}${SUFFIX} ${freq}GHz REAL "${SWB_CMD}" ${NB_PROC} ${walltime} ${REP_CMD} BeginRing=240 EndRing=11194
      python doQsub.py ${PREFIX}${freq}${SWB_LABEL}_YR2${REP_LABEL}${SUFFIX} ${freq}GHz REAL "${SWB_CMD}" ${NB_PROC} ${walltime} ${REP_CMD} BeginRing=11195 EndRing=21720
  fi

done
