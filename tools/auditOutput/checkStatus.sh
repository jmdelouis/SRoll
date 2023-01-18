#!/bin/bash

###############################################################################
# checkStatus.sh
#
# This script allow to see quickly which maps are Done and other useful info.
# One can pass a directory as argument to get info for maps in this dir, by
# default it focus on the LAST_PROD_DIR.
#
# Author  : Christian Madsen
# Date    : 2015-10-21
# version : Initial
###############################################################################


# Last production dir
LAST_PROD_DIR="/pscratch1/delouis/tmp_CM_SROLL_OUT/sroll_vlast" # The reference directory for which the script print result by default

# Where to find the sroll_*.o* qsub log files
REF_DIR_FOR_STATUS="/wrk/cmadsen/SrollEx/pipe"
REF_DIR_FOR_STATUS_BIS="${REF_DIR_FOR_STATUS}/old"


# Init counters
count_Done=0
count_Queued=0
count_Error=0
count_TOT=0


#------------------------------------------------------------------------------

# Try to retrieve the walltime associated to the generated data in dir pass in
# arg.
retrieveWalltime() {

  local pathDir="$1"
  local filename="jobID"

  # Ensure the "jobID" file exist (otherwise return 1)
  [[ ! -f "${pathDir}/${filename}" ]] && return 1

  # Extract jobID info from file
  local jobId=`cat "${pathDir}/${filename}" | cut -d . -f 1`

  # Determine the appropriate directory where to find log file (REF_DIR_FOR_STATUS or REF_DIR_FOR_STATUS_BIS)
  if [ -f ${REF_DIR_FOR_STATUS}/*.o${jobId} ]
    then
      currentRefDir="${REF_DIR_FOR_STATUS}"
    elif [ -f ${REF_DIR_FOR_STATUS_BIS}/*.o${jobId} ]
      then
        currentRefDir="${REF_DIR_FOR_STATUS_BIS}"
      else
        # Unable to find log file
        return 1
  fi

  local line=`grep RESOURCESUSED ${currentRefDir}/*.o${jobId} 2> /dev/null`

  if [ $? -ne 0 ]
    then
      walltime="???"
      return 2
  fi

  walltime=`echo "${line}" | cut -d = -f5`

  return 0
}

#------------------------------------------------------------------------------

# Try to retrieve the job exit status associated to the generated data in dir pass in
# arg.
retrieveJobExitStatus() {

  local pathDir="$1"
  local filename="jobID"

  # Ensure the "jobID" file exist (otherwise return 1)
  [[ ! -f "${pathDir}/${filename}" ]] && return 1

  # Extract jobID info from file
  local jobId=`cat "${pathDir}/${filename}" | cut -d . -f 1`

  # Determine the appropriate directory where to find log file (REF_DIR_FOR_STATUS or REF_DIR_FOR_STATUS_BIS)
  if [ -f ${REF_DIR_FOR_STATUS}/*.o${jobId} ]
    then
      currentRefDir="${REF_DIR_FOR_STATUS}"
    elif [ -f ${REF_DIR_FOR_STATUS_BIS}/*.o${jobId} ]
      then
        currentRefDir="${REF_DIR_FOR_STATUS_BIS}"
      else
        # Unable to find log file
        return 1
  fi

  local line=`grep "JOB EXIT STATUS" ${currentRefDir}/*.o${jobId} 2> /dev/null`

  if [ $? -ne 0 ]
    then
      jobExitStatus="???"
      return 2
  fi

  jobExitStatus=`echo "${line}" | cut -d " " -f 5`

  return 0
}

#------------------------------------------------------------------------------

# Handler user optional directory
if [ -n "$1" ]
  then
    targetDir="$1"
  else
    targetDir="${LAST_PROD_DIR}"
fi

# Check that dir exist! If not print error message and exit
if [ ! -d "${targetDir}" ]
  then
    echo "ERROR: Directory '${targetDir}' does not exist!"
    exit 1
fi

echo ""
echo "Printing stats for maps in directory: '${targetDir}'"
echo ""

cd "${targetDir}"

# Loop on all maps dirs
for idir in 100* 143* 217* 353* 545* 857*
do
  # 0) check that dir exist... (if not, nothing to do!)
  if [ ! -d "${idir}" ]
    then
      continue
  fi

  ((count_TOT++))
  # Reset var walltime
  unset walltime
  unset jobExitStatus

  retrieveJobExitStatus "${idir}"
  res_jobExitStatus=$?

  # Check status
  grep "Finished sucessfully" ${idir}/DBG*.log &>/dev/null
  case $? in
    0) ((count_Done++))
       retrieveWalltime "${idir}"
       # Just to be sure that qsub log is coherent with log
       if [ ${jobExitStatus} -eq 0 ]
         then 
           result="DONE"
         else
           result="WEIRD! (job exit status: ${jobExitStatus})"
       fi
       ;;
    1)
       flag_error_tmp=0
       # Try to better guess situation
       grep -i "error" ${idir}/DBG*.log &>/dev/null
       res_error=$?
       grep -i "refused to die" ${idir}/DBG*.log &>/dev/null
       res_die=$?
       if [ ${res_error} -eq 0 ] || [ ${res_die} -eq 0 ]  # One of the "error" patterns has been found!
         then
           flag_error_tmp=1
       elif [ ${res_jobExitStatus} -eq 0 ] && [ ${jobExitStatus} -ne 0 ]
         then
           flag_error_tmp=1
       fi

       if [ ${flag_error_tmp} -eq 1 ]
         then
           result="** ERROR **"
           ((count_Error++))
         else
           result="*R*  Currently processed..."
       fi
       ;;
    2) result="_._  Queued"
       ((count_Queued++))
       ;;
    *) result="???"
       ;;
  esac

  # Print result to user
  printf "%-25s ---> %4s" "${idir}" "${result}"
  # Try to retrieve the corresponding wall time (if any)
  if [ -n "${walltime}" ] 
    then
      printf " (%s)" "${walltime}"
    else
      printf "%4s" ""
  fi

  # Retrieve sroll.c CVS revision
  cvsSrollSrcRevision=$( cat ${idir}/srollREV &>/dev/null )
  case $? in
    0) revInfo=$( cat ${idir}/srollREV | grep "Working revision" | cut -f2 )
       ;;
    *) revInfo="??"
  esac

  printf "  [v ${revInfo}]\n"

done

### Print some additional info
echo ""

if [ ${count_TOT} -eq 0 ]
  then
    echo "Empty production dir... (nothing to say)"
    exit 0
fi

printf "TOTAL:\n"
printf "  * DONE: %s/%s (%d%%)\n" "${count_Done}" "${count_TOT}" "$((count_Done*100/count_TOT))"
if [ $((${count_Done}+${count_Error})) -ne ${count_TOT} ]
  then
    estimatedRemainingTime=$(((count_Queued+1)*3))
    printf "  * Estimated end in %dh --> \"%s\"\n" "${estimatedRemainingTime}" "$(date '+%a %d %b %Y %R' -d "+${estimatedRemainingTime} hours")"
fi
