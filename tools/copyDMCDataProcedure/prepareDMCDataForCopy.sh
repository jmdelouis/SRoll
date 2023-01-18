#!/bin/bash

###############################################################################
# prepareDMCDataForCopy.sh
#
# NOTE:
#   This script is to be executed on M3 (with hl2 env set), AFTER having exec
#   on target host the script "checkParamDataAvailability.sh".
#
# The script will prepare data on M3 for the copy process. Mainly it update or
# create the metadata for corresponding objects.
# Note that this script will take care to handle object TUPLE sub dirs!
#
# Usage:
# $> ./prepareDMCDataForCopy.sh [--no-metadata-update] <to_be_copied.txt>
# The "--no-metadata-update" allow to generate output files without regenerating
# the metadata.
#
# Author  : Christian Madsen
# Date    : 2015-07-23
# version : Initial
###############################################################################

### Some constants
# The directory where are stored the DMC object on M3 (the source)
PREFIX_DIR_M3="/data"


### Check and Parse command line parameters
if [[ $# -gt 2 ]]
  then
    echo "You must specify one file containing DMC object path!"
    echo "Usage:"
    echo '$> ./prepareDMCDataForCopy.sh [--no-metadata-update] <to_be_copied.txt>'
    exit 1
fi

# Handle optional parameter
option=$1
flag_noMetadataUpdate=0
if [ "${option}" == "--no-metadata-update" ]
  then
    echo ">> No metadata update!"
    flag_noMetadataUpdate=1
    shift
fi


inputObjs=$( cat $1 )
#echo "[DBG] inputObjs=${inputObjs}"


### First we need to replace prefix string with appropriate path dir for M3
inputObjs_M3=${inputObjs//__PREFIX__/${PREFIX_DIR_M3}}
#echo "[DBG] inputObjs=${inputObjs}"

### 0) Take care of possible TUPLE_x sub dirs!
# For this we loop over each entry and try to find any existing TUPLE sub dir,
# if found it is added to the list of object to be processed
suffix_Tuple="_TUPLE_"
echo
echo "* Trying to find TUPLE_x sub dir associated to objects..."
echo
declare -a additionalTUPLEobj
on_error=false
# Init counter
currentPos=1
totPos=$( wc -l <<< "${inputObjs_M3}" )
for obj in ${inputObjs_M3};
do
  tmpDirname=$( dirname "${obj}" )
  tmpBasename=$( basename "${obj}" )
  nameToCheck="${tmpBasename}${suffix_Tuple}"

  #echo "[DBG] tmpDirname=$tmpDirname"
  #echo "[DBG] tmpBasename=$tmpBasename"
  #echo "[DBG] nameToCheck=$nameToCheck"

  echo "  -> Checking for \"${tmpBasename}\" [${currentPos}/${totPos}]..."

  # Make a precheck to ensure that even the base file exist! If not this is an ERROR!
  if [ ! -d "${obj}" ]
    then
      echo "ERROR: Missing file on M3 '${obj}'"
      on_error=true
      (( currentPos++ ))
      continue
  fi

  count_tuple=1
  search_end=false
  res="" # init with empty result
  until [ "${search_end}" = true ]; do
      currentTupleName="${tmpDirname}/${nameToCheck}${count_tuple}"
      #echo "[DBG] name to be checked: '${currentTupleName}'"
      # Try to detect if there is a tuple "_TUPLE_X" (where 'X' is the current index)
      if [ -d "${currentTupleName}" ]
        then
          #echo "[DBG] found TUPLE_${count_tuple}"
          # Add it to the result list
          if [[ ! -z "${res}" ]]
            then
              res="$res\n${currentTupleName}"
            else
              res="${currentTupleName}"
          fi
          (( count_tuple++ ))
        else
          search_end=true
      fi
  done
  # Update counter to have a human readable value
  (( count_tuple-- ))

#OLD slowww function...  res=$( find "${tmpDirname}" -maxdepth 1 -type d -name "${nameToCheck}*" )

  #echo -e "[DBG] res='$res'"

  # If any tuple found add them to the tuple list
  if [[ ! -z "${res}" ]]
  then
    echo "     - ${count_tuple} tuples found!"
    additionalTUPLEobj+=("$res")
  else
    echo "     - No tuple"
  fi

  # Update counter
  (( currentPos++ ))
done

# Finally add any founded TUPLE to the list
if [[ ${#additionalTUPLEobj[@]} -gt 0 ]]
then
  tmpTuples=$( printf '%s\n' "${additionalTUPLEobj[@]}" )
  inputObjs_M3=$( printf "${inputObjs_M3}\n${tmpTuples}" )
fi

#echo "[DBG] inputObjs_M3=${inputObjs_M3}"


# Exit if error occured!
if [ "${on_error}" = true ]
  then
    echo "EXIT due to previously encountered errors!"
    exit 1
fi


### 1) Loop on each obj in order to update/create metadata on M3 (if required)
echo
echo "* Processing metadata for object on M3..."
echo
output_pio=""
# Reset counter
currentPos=1
totPos=$( wc -l <<< "${inputObjs_M3}" )
for obj in ${inputObjs_M3};
do
  # Only force update if necessary
  if [ ${flag_noMetadataUpdate} -ne 1 ]
    then
      echo "  Processing object: ${obj}"
      echo "    -> Forcing metadata update [${currentPos}/${totPos}]..."
      python ../noDMCMetadata/no_dmc_metadata.py -O ${obj}
      if [ $? -ne 0 ]
      then
        echo "Error during generation of metadata!"
        exit 1
      fi
  fi

  echo "    -> Extracting backendname..."
  # Extracting backendname from metadata
  backendname=`grep Backendname ${obj}/no_dmc_metadata.txt | cut -d" " -f3`
  # Add new line to output (if required)
  if [[ ! -z "${output_pio}" ]]
    then
      output_pio="${output_pio}\n"
  fi
  # Then add the current entry to be saved
  output_pio=$(printf "${output_pio}${backendname}")
  #echo "[DBG] output_pio=${output_pio}"

  # Update counter
  (( currentPos++ ))
done

### 2) Saving path to data and pio in 2 separate files
echo
echo "* Saving results..."
echo
# Retrieve input filename to specialize output name
tmpInputBasename=$( basename "$1" )
# Construct output filenames
filename_DATA="From_M3_DATA_For_${tmpInputBasename}"
filename_PIO="From_M3_PIO_For_${tmpInputBasename}"
# Save files
echo "${inputObjs_M3}" > "${filename_DATA}"
echo -e "${output_pio}" > "${filename_PIO}"
# Output info for user
echo "    -> DATA path saved to \"${filename_DATA}\""
echo "    -> PIO path saved to \"${filename_PIO}\""
