#!/bin/bash

###############################################################################
# checkParamDataAvailability.sh
#
# This script allow to check the availability, on the current host (ex: M4 or
# EDISON), of input data specified in a parameter file.
# Mainly, using this script before calling sroll avoid to failed due to missing
# data.
# Note that depending of the option given, the script produce full or partial
# product dependency list. 
# Note also that the script can handle generic parameter files or already
# specialized ones.
#
# WARNING: order of options matter!
#
# Usage: see function "usage()" for more info.
#
# Author  : Christian Madsen
# Date    : 2015-07-23
# version : Initial
###############################################################################


#------------------------------------------------------------------------------

### Define bolo set per freq
declare -A boloList=(
                      ["100"]="00_100_1a 01_100_1b 80_100_4a 81_100_4b 20_100_2a 21_100_2b 40_100_3a 41_100_3b"
                      ["143"]="02_143_1a 03_143_1b 50_143_3a 51_143_3b 30_143_2a 31_143_2b 82_143_4a 83_143_4b 10_143_5 42_143_6 60_143_7"
                      ["217"]="11_217_5a 12_217_5b 61_217_7a 62_217_7b 43_217_6a 44_217_6b 71_217_8a 72_217_8b 04_217_1 22_217_2 52_217_3 84_217_4"
                      ["353"]="23_353_3a 24_353_3b 53_353_5a 54_353_5b 32_353_4a 33_353_4b 63_353_6a 64_353_6b 05_353_1 13_353_2 45_353_7 85_353_8"
                      ["545"]="14_545_1 34_545_2 73_545_4"
                      ["857"]="25_857_1 35_857_2 65_857_3 74_857_4"
    )

#------------------------------------------------------------------------------

# Covert fullboloid to boloid. "14_545_1" -> "545-1"
getBoloId() {
  # Remove first 3 char
  boloId=${1:3}

  # Replace '_' by '-'
  boloId=`echo ${boloId} | tr '_' '-'`
}

#------------------------------------------------------------------------------

# Display usage information about this script and exit.
usage() {
  echo "Usage: $0 <--option> [<--repX>] [<--hrX>] [<--freqXXX>] <filename>"
  echo "  <--option> is one of \"--full\" or \"--only-missing\""
  echo "  [<--SIM>] (optional) allow to handle SIM par files"
  echo "  [<--repX>] (optional) allow to specify the exact REP id version (default \"--rep6\")"
  echo "  [<--hrX>] (optional) is one of \"--hr1\" or \"--hr2\" (default is no half ring)"
  echo "  [<--parityX>] (optional) is one of \"--parityOdd\" or \"--parityEven\" (default is no parity)"
  echo "  [<--freqXXX>] (optional) is of type \"--freq100\" where 100 is the desired freq"
  echo "  <filename> is the parameter file to be processed"
  exit 1
}

#------------------------------------------------------------------------------

# Test existance of the filename specified as argument and return the result. 
is_file_exits() {
  local f="$1"
  [[ -f "$f" ]] && return 0 || return 1
}

#------------------------------------------------------------------------------

# Test existance of the directory specified as argument and return the result. 
is_dir_exits() {
  local d="$1"
  [[ -d "$d" ]] && return 0 || return 1
}

#------------------------------------------------------------------------------

### MAIN ###

hostDataDir_M4="/pscratch1/RD12_data"
hostDataDir_EDISON="/scratch2/scratchdirs/delouis"

DEFAULT_REPID="REP6"

### 0) Detect host for autoconfig
HOST=$( hostname )
#echo "[DBG] HOST=${HOST}"
case "${HOST}" in
  log*)
    hostDataDir=${hostDataDir_M4}
    echo "[DBG] M4 detected !";;
  edison*)
    hostDataDir=${hostDataDir_EDISON}
    echo "[DBG] EDISON detected!";;
  *)
    echo "ERROR: unsupported host!"
    exit 1;;
esac
#echo "[DBG] hostDataDir=${hostDataDir}"


### A) Check command line input

# We required user to specify at least 2 parameters
if [[ $# -lt 2 || $# -gt 5 ]]
  then
    usage
fi

# Check Option
option=$1

case "${option}" in
  --full)
    echo "Option --full activated"
    ;;
  --only-missing)
    echo "Option --only-missing activated"
    ;;
  *)
    echo "ERROR: unsupported option! ($1)"
    usage
    ;;
esac

# Check for optional --SIM flag
srollIn_pattern="srollIn"

opt_SIM=$2
if [[ "${opt_SIM}" == --SIM ]] # Test if option is --SIM 
  then
    srollIn_pattern="/pscratch1/RD12_data" # Make use of the path used in SIM par files!
    # Update positional argument
    shift
fi


# Check for optional --repX flag
repID=${DEFAULT_REPID}

opt_rep=$2
if [[ "${opt_rep}" == --rep* ]] # Test if option start with --rep
  then
    # Retrieve the specified value
    repID="REP${opt_rep:5}"
    # Update positional argument
    shift
fi
REP_label="__${repID}__"

# Check for optional half ring flag
opt_hr=$2

if [[ "${opt_hr}" == --hr* ]] # Test if option start with --hr
  then
    case "${opt_hr}" in
      --hr1)
        echo "Option *HalfRing 1st* activated"
        replace_pattern_for_HR="${repID}_1st"
        HR_label="__HR1st__"
        ;;
      --hr2)
        echo "Option *HalfRing 2nd* activated"
        replace_pattern_for_HR="${repID}_2nd"
        HR_label="__HR2nd__"
        ;;
      *)
        echo "ERROR: unsupported option! (${opt_hr})"
        usage
        ;;
    esac
    shift
fi

# Check for optional --parityX flag
parity_Label=""

opt_parity=$2
if [[ "${opt_parity}" == --parity* ]] # Test if option start with --parity
  then
    # Retrieve the specified value
    #parityID="${opt_parity:8}"
    case "${opt_parity}" in
      --parityOdd)
        echo "Option *Parity Odd* activated"
        parity_Label="_Odd"
        ;;
      --parityEven)
        echo "Option *Parity Even* activated"
        parity_Label="_Even"
        ;;
      *)
        echo "ERROR: unsupported option! (${opt_parity})"
        usage
        ;;
    esac

    # Update positional argument
    shift
fi


# Check for optional freq flag
opt_freq=$2

if [[ "${opt_freq}" == --freq* ]] # Test if option start with --freq
  then
    freq="${opt_freq:6}"
    echo ">> User specify freq: ${freq}GHz"
    shift

#    echo "[BDG] bolo list for ${freq}GHz: ${boloList[${freq}]}"
#    for i in ${boloList[${freq}]}
#    do
#      getBoloId "$i"
#      echo "TMP: $i   boloid: ${boloId}"
#    done
fi


# Read next parameter
parameterFile=$2


# Check parameter filename

if !( is_file_exits "${parameterFile}" )
  then
    echo "You must specify an existing parameter file!"
    exit 1
fi

echo ">> Parameter file to be tested: '${parameterFile}'"

### B) Extract all input data from param file

# First try with generic pattern "${srollIn_pattern}"
# Extract info to tmp file
cat ${parameterFile} | grep "${srollIn_pattern}" > tmp.txt

# If previous result is empty try with already speciliazed input path (${hostDataDir})
if [ ! -s tmp.txt ]
  then
    # Special case of already specialize parameter file (named "param.txt") -> nothing to specialize
    #echo "[DBG] special case of param.txt"
    cat ${parameterFile} | grep "${hostDataDir}" > tmp.txt
fi

# Only keep the path
cat tmp.txt | cut -f2 -d "=" > tmp2.txt

echo ">> We specialize object name..."

### Specialize parameter file

# i) REPID
if [ -n "${replace_pattern_for_HR}" ]
  then # (optional) HR replacement
    echo "   - Using HR pattern: '${replace_pattern_for_HR}'"
    sed -i "s/\[\[%REPID\]\]/${replace_pattern_for_HR}/" tmp2.txt
  else
    echo "   - Using REP pattern: '${repID}'"
    sed -i "s/\[\[%REPID\]\]/${repID}/" tmp2.txt
fi

# ii) Parity
sed -i "s/\[\[%OddEven\]\]/${parity_Label}/" tmp2.txt

# iii) "{boloID}" and "{boloIDFull}"
# TODO ...
# - 1) Retrieve lines where there is "{boloID}" and "{boloIDFull}"
genericLines=`grep -E '{boloID}|{boloIDFull}' tmp2.txt`
# - 2) Secondly iterate throw each required bolo for each these generic pattern
for line in ${genericLines}
do
  #echo "Process line: $line"
  tmpNewLines=""
  for i in ${boloList[${freq}]}
  do
    getBoloId "$i"
    #echo "TMP: $i   boloid: ${boloId}"
    tmp=`echo ${line} | sed "s/{boloID}/${boloId}/" -`
    tmp=`echo ${tmp} | sed "s/{boloIDFull}/$i/" -`
    tmpNewLines=`echo -e "${tmpNewLines}\n$tmp"`
  done
  #echo "[DBG] tmpNewLines=${tmpNewLines}"
  echo -e "${tmpNewLines}" > tmpGen.txt
  # Escaping char for using expression in sed
  line_ESC=$(sed 's/\[/\\&/g' <<<"${line}")
  line_ESC=$(sed 's/\]/\\&/g' <<<"${line_ESC}")
  line_ESC=$(sed 's/\//\\&/g' <<<"${line_ESC}")
  #echo "[DBG] line_ESC=${line_ESC}"

  # Replace the generic line by the set of all bolo lines
  sed -i "/^ ${line_ESC}$/r tmpGen.txt" tmp2.txt

  # Remove the generic line
  sed -i "/^ ${line_ESC}$/d" tmp2.txt
done

# Remove temporary files
rm -f "tmp.txt" "tmp2.txt" "tmpGen.txt"


# Loop on all entry and check if it exist on the local host data dir
count_OK=0
count_MISS=0
count_BAD_METADATA=0
declare -a list_miss
declare -a list_bad
declare -a list_all
for i in `cat tmp2.txt`
  do
    # Replace [[%srollIn]] by the effectiv directory on current host.
    v="${i/\[\[%srollIn\]\]/${hostDataDir}}"
    # Print the corresponding path
    echo -n $v
    # Store it in the full list
    list_all+=("$v")
    # a) Check that the related data directory exist
    if [ -d "$v" ]
      then
        # b) Check that the metadata is ok (ie. backendname path exist)
        backendname=`grep Backendname ${v}/no_dmc_metadata.txt | cut -d" " -f3`
        #echo "[DBG] backendname='${backendname}'"
        if !( is_dir_exits "${backendname}" )
          then
            list_bad+=("$v")
            ((count_BAD_METADATA++))
            echo "  --> BAD METADATA"
            echo "[DBG] backendname='${backendname}'"
          else
            ((count_OK++))
            echo "  --> OK"
        fi
      else
        ((count_MISS++))
        list_miss+=("$v")
        echo "  --> MISSING"
    fi
  done
echo
echo "TOTAL (${#list_all[@]}): ${count_OK} OK, ${count_MISS} MISSING, ${count_BAD_METADATA} BAD METADATA"
echo

### C) Save result in file (if required)

baseParam=$( basename ${parameterFile} )
outputFile="TO_BE_COPIED_FOR_${REP_label}${HR_label}${parity_Label}_${baseParam}"

# Switch behaviour depending on input user option
declare -a list_to_be_saved
if [[ ${option} == "--full" ]]
  then
    list_to_be_saved=( ${list_all[@]} )
  else
    list_to_be_saved=( ${list_miss[@]} )
fi

if [[ ${#list_to_be_saved[@]} -gt 0 ]]
  then
    #echo "[DBG] tab has ${#list_to_be_saved[@]} elts!"
    # Convert array so that it can be saved one entry per line in output text file
    toto=$( printf '%s\n' "${list_to_be_saved[@]}" )

    # Update path with generic prefix (so that it can be used on M3 easily)
    patternToBeReplaced=${hostDataDir}
    #echo "[DBG] patternToBeReplaced=${patternToBeReplaced}"
    t=${toto//${patternToBeReplaced}/__PREFIX__}
    echo "$t" > ${outputFile}

    # Print message
    echo "Result saved to: \"${outputFile}\""
    echo -e "\n\n>>> Next step is to go on M3 and execute script 'prepareDMCDataForCopy.sh' with this file as parameter."
  else
    echo "Nothing to be copied! (no file saved)"
fi

if [[ ${#list_bad[@]} -gt 0 ]]
  then
    echo -e "\n\nIMPORTANT: Also consider to update metadata for mentioned objects!"
fi

exit 0
