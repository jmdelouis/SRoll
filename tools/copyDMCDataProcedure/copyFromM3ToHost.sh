#!/bin/bash

###############################################################################
# copyFromM3ToHost.sh
#
# NOTE: This script is to be used from a dedicated host (M4, EDISON or CC)!
#       Automatic host detection is proceed to adapt path configuration...
#       For now ONLY "M4", "CORI" (NERSC) and "CC" are supported!
#
# The script will:
# 1) Copy required data and metadata to host (in the M4 case, local mount of M3
#    filesystem will be used for improved copy perf).
# 2) Update metadata info so that it comply with the new localization of
#    objects on the host system. (cf 'backendname')
#
# Usage:
# $> ./copyFromM3ToHost.sh <to_be_copied_M3_DATA.txt> <to_be_copied_M3_PIO.txt>
# or
# $> ./copyFromM3ToHost.sh --only-pio <to_be_copied_M3_PIO.txt>
#
# Author  : Christian Madsen
# Date    : 2015-07-23
# version : Initial
###############################################################################

#------------------------------------------------------------------------------

# Display usage information about this script and exit.
usage() {
  echo "Usage: $0 [<--only-pio>] [<data.txt>] <pio.txt>"
  echo "  <--only-pio> Allow to specify only one pio text file."
  exit 1
}

#------------------------------------------------------------------------------

# Test existance of the filename specified as argument and return the result. 
is_file_exits() {
  local f="$1"
  [[ -f "$f" ]] && return 0 || return 1
}

#------------------------------------------------------------------------------

###############################################################################
#################################### MAIN #####################################
###############################################################################

### Some constants
# The directory where are stored the DMC object on M3 (the source)
FORMER_PREFIX_M3_DATA="/data/dmc/MISS03/DATA/"
FORMER_PREFIX_M3_PIO="/redtruck/dmc/dmc_objects/MISS03/"

# Some info for M4 copy optimization
PREFIX_MOUNT_M3onM4_DATA="/m3gpfs3/datadmc/dmc/MISS03/DATA/" # replace FORMER_PREFIX_M3_DATA
PREFIX_MOUNT_M3onM4_PIO="/m3gpfs3/dmc/dmc_objects/MISS03/" # replace FORMER_PREFIX_M3_PIO

# Script for updating metadata
pyScript4UpdateMetadata="../noDMCMetadata/no_dmc_metadata_fix_backendname.py"


### Detect current host (to update accordingly the paths to be used)
### 0) Detect host for autoconfig
HOST=$( hostname )
#echo "[DBG] HOST=${HOST}"
flag_M4_optim=false
case "${HOST}" in
  log[0-3]*)
    echo "M4 detected"
    targetDir_DATA="/pscratch1/RD12_data/dmc/MISS03/DATA/"
    targetDir_PIO="/pscratch1/RD12_data/dmc/dmc_objects/MISS03/"
    flag_M4_optim=true
    ;;
  edison*)
    echo "Edison detected"
    targetDir_DATA="/scratch3/scratchdirs/paganol/RD12_data/dmc/MISS03/DATA/"
    targetDir_PIO="/scratch3/scratchdirs/paganol/RD12_data/dmc/dmc_objects/MISS03/"
    if [ "${LOGNAME}" = "smottet" ]
      then
        targetDir_DATA="/scratch2/scratchdirs/smottet/RD12_data/dmc/MISS03/DATA/"
        targetDir_PIO="/scratch2/scratchdirs/smottet/RD12_data/dmc/dmc_objects/MISS03/"
      fi
    ;;
  cca*)
    echo "CC detected"
    targetDir_DATA="/sps/planck/SimuData/RD12_data/dmc/MISS03/DATA/"
    targetDir_PIO="/sps/planck/SimuData/RD12_data/dmc/dmc_objects/MISS03/"
    ;;
  cori*)
    echo "CORI detected"
    targetDir_DATA="/project/projectdirs/planck/data/hfi/RD12_data/dmc/MISS03/DATA/"
    targetDir_PIO="/project/projectdirs/planck/data/hfi/RD12_data/dmc/dmc_objects/MISS03/"
    ;;
  login[0-9]*)
    echo "OCCIGEN detected"
    targetDir_DATA="/scratch/cnt0028/ias1717/SHARED/RD12_data/dmc/MISS03/DATA/"
    targetDir_PIO="/scratch/cnt0028/ias1717/SHARED/RD12_data/dmc/dmc_objects/MISS03/"
    ;;
  *)
    echo "ERROR: unsupported host"
    exit 1;;
esac
echo "targetDir_DATA=${targetDir_DATA}"
echo "targetDir_PIO =${targetDir_PIO}"


### Check and Parse command line parameters
if [[ $# -ne 2 ]]
  then
    usage
fi

# Parse option '--only-pio' (if any)
flag_only_pio=false
if [ "$1" = "--only-pio" ]
  then
    flag_only_pio=true
fi

m3login=""
# add your username here for automatic login translation...
if [ "${LOGNAME}" = "smottet" ]
  then
    m3login="symottet"
fi

if [ ! "${flag_M4_optim}" = true ] && [ -z ${m3login} ]
  then
    echo "Please enter your 'login' name on magique3: "
    read m3login
fi

echo "m3login='${m3login}'"


# Adapt input depending on user option (cf. --only-pio)
if [ "${flag_only_pio}" = true ]
  then
    if !( is_file_exits "$2" )
      then
        echo "ERROR: You must specify an existing file!"
        exit 1
    fi
    file_PIO=$( cat $2 )
  else
    if !( is_file_exits "$1" ) || !( is_file_exits "$2" )
      then
        echo "ERROR: You must specify an existing files!"
        exit 1
    fi
    file_DATA=$( cat $1 )
    file_PIO=$( cat $2 )
fi

#echo "[DBG] file_DATA='${file_DATA}'"
#echo "[DBG] file_PIO ='${file_PIO}'"

RSYNC_OPTION="-av --no-perms --omit-dir-times -h"
#RSYNC_OPTION_DEBUG="--progress"

# Prepare filemode for the mkdir
umask 0000

#==============================================================================
# COPY PIO
#==============================================================================

### 1) Copy obj
echo
echo "* Copying objects PIO..."
echo
# Init counter
currentPos=1
totPos=$( wc -l <<< "${file_PIO}" )
for obj in ${file_PIO};
do
  echo "  Processing object [${currentPos}/${totPos}]: ${obj}"

  if [ "${flag_M4_optim}" = true ]
    then
      echo "*** Optimization for M4 copy! (flag is true)"
      sourceObj=${obj/${FORMER_PREFIX_M3_PIO}/${PREFIX_MOUNT_M3onM4_PIO}}
    else
      echo "*** NO optimization! (flag is false)"
      sourceObj="${m3login}@planck31.iap.fr:${obj}"
  fi
  
  targetObj=${obj/${FORMER_PREFIX_M3_PIO}/${targetDir_PIO}}

  echo "sourceObj=${sourceObj}"
  echo "targetObj=${targetObj}"

  # First create directory
  mkdir -p "${targetObj}"

  # Copy using rsync
  rsync ${RSYNC_OPTION} ${RSYNC_OPTION_DEBUG} "${sourceObj}/" "${targetObj}"
  if [ $? -ne 0 ]
  then
    echo "ERROR during rsync!"
    exit 1
  fi 

  if [ "${SROLLHOST}" = "CC" ]
  then
    chmod -R 700 "${targetObj}"
  fi

  # Update counter
  (( currentPos++ ))
done


#==============================================================================
# COPY DATA (optional)
#==============================================================================

if [ "${flag_only_pio}" != true ]
  then
    echo
    echo "* Copying objects DATA..."
    echo
    # Init counter
    currentPos=1
    totPos=$( wc -l <<< "${file_DATA}" )
    for obj in ${file_DATA};
    do
      echo "  Processing object [${currentPos}/${totPos}]: ${obj}"
      
      if [ "${flag_M4_optim}" = true ]
        then
          echo "*** Optimization for M4 copy! (flag is true)"
          sourceObj=${obj/${FORMER_PREFIX_M3_DATA}/${PREFIX_MOUNT_M3onM4_DATA}}
        else
          echo "*** NO optimization! (flag is false)"
          sourceObj="${m3login}@planck31.iap.fr:${obj}"
      fi
      targetObj=${obj/${FORMER_PREFIX_M3_DATA}/${targetDir_DATA}}

      echo "sourceObj=${sourceObj}"
      echo "targetObj=${targetObj}"

      # First create directory
      mkdir -p "${targetObj}"

      # Copy using rsync
      rsync ${RSYNC_OPTION} ${RSYNC_OPTION_DEBUG} "${sourceObj}/" "${targetObj}"
      if [ $? -ne 0 ]
      then
        echo "ERROR during rsync!"
        exit 1
      fi 

      ### 2) Updating backendname
      python "${pyScript4UpdateMetadata}" "${targetObj}" "${FORMER_PREFIX_M3_PIO}" "${targetDir_PIO}"
      if [ $? -ne 0 ]
      then
        echo "ERROR during metadata update!"
        exit 1
      fi 

      if [ "${SROLLHOST}" = "CC" ]
      then
        chmod -R 700 "${targetObj}"
      fi

      # Update counter
      (( currentPos++ ))
    done
fi


### Final user message
echo ""
echo "ALL DONE :)"

