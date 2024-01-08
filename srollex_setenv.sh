#!/bin/bash

###############################################################################
# srollex_setenv.sh
#
# The script aims to guaranty that SrollEx compile and execute as expected on
# the host it is used on.
# Supported host are: M3, M4.
#
# IMPORTANT:
# To use this script you have to **source** it!
# $> source srollex_setenv.sh
#
# Author  : Christian Madsen
# Date    : 2015-10-01
# version : Initial
#
# 2016-03-24 - mottet@iap.fr: added SROLLDIR, SROLLHOST, PYTHONPATH,
#                             LD_LIBRARY_PATH and NODMCDATA environment
#                             variables definition/update
###############################################################################

module () {
  eval $( $(which modulecmd) bash $*)
}

if [ ! -z ${SROLLDIR} ]
  then
    echo "WARNING!"
    echo "setenv.sh already sourced in ${SROLLDIR}"
    echo "to use a different working directory please log out and log back in to clean your environment"
    return 1
fi




### Detect host for autoconfig
HOST=$( hostname )

case "${HOST}" in    
  br146-050*)
 	echo " br146-050 detected " 
	export SROLLHOST=br146-050
	 
	export PYTHONPATH=/export/home/jmdeloui/SROLL/py_sroll/ 
	export LD_LIBRARY_PATH=/export/home/jmdeloui/SROLL/py_sroll/:$LD_LIBRARY_PATH 
 ;;   
  garoupe*)
 	echo " garoupe detected " 
	export SROLLHOST=garoupe
	 
	export PYTHONPATH=/export/home/jmdeloui/SROLL/py_sroll/ 
	export LD_LIBRARY_PATH=/export/home/jmdeloui/SROLL/py_sroll/:$LD_LIBRARY_PATH 
 ;;
  *)
    echo "WARNING (setenv.sh): unknown host '${HOST}'! (Nothing done)"
    return 1
    ;;
esac

#Remove Tensorflow warning  
export TF_CPP_MIN_LOG_LEVEL="2"
    
    
# sroll root directory is where the srollex_setenv.sh script is
SROLLDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export SROLLDIR=${SROLLDIR}

