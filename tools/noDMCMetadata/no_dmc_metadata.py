#!/usr/bin/env python

##############################################################################
# no_dmc_metadata.py
#
# This python script aims to generate the metadata files that will be required
# for the noDMCLib to work.
#
# IMPORTANT: Execution of this script required to have previously set an
#            hl2_env since the script make use of DMC queries!
#
# Author  : Christian Madsen
# Date    : 2014-12-18
# version : Initial
##############################################################################

import os

### Ensure the user have previously set an hl2_env (ex: "DMCNAME" var should
#   exist)
if "DMCNAME" not in os.environ:
  print("ERROR: You shall have set an 'hl2_user_init' before calling this script!")
  exit(1)


import sys
import piolib as pio
import subprocess
import psycopg2 # for psql queries

# Debug and verbose mode
DEBUG = False # Set to True if you want debug message to be print
VERBOSE = False # Set to True if you want more log message

### Parameters
metadata_filename = "no_dmc_metadata.txt"
VERSION = "0.0.2" # metadata version (used to detect deprecated metadata)
VERSION = "0.0.3" # SM 16March2017: added md5sum
dbpath = pio.GetDataDB() # ex: '/data/dmc/MISS03/DATA'

#dbname = dbpath.split('dmc/')[1].split('/DATA')[0] # Extract data basename ex: MISS03
if dbpath == "/data/dmc/MISS03/DATA":
  dbname = "MISS03"
elif dbpath == "/sps/planck/DMC/MISS03cc/DATA":
  dbname = "MISS03cc"
else:
  raise Exception("unknown DMC database: " + dbpath)

### Binaries (external tools)
psql_bin = "psql"
dmctool_bin = "dmctool"


DMC_TABLE_OBJECT_NAME_GEN = "dmc_piolib_%sobject" # % object type (in lowercase)

NODMCMETADATA_DATE_TAG = "noDMCmetadata_date"
NODMCMETADATA_VERSION_TAG = "noDMCmetadata_version"


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

def getMetadataIdentifierContent():
  """This function return some entries allowing to identify the metadata process.
  """
  import time

  metadata_process_info = "%s : %s" % (NODMCMETADATA_DATE_TAG, time.strftime("%Y-%m-%d %H:%M:%S")) # some comment
  metadata_process_info += "\n%s : %s" % (NODMCMETADATA_VERSION_TAG, VERSION)

  return metadata_process_info

#-----------------------------------------------------------------------------

# THIS IS THE OLD WAY OF DOING... the lib way is almost 2x faster!!!
def OLD_psqlRequestByShell(dmc_obj_name):
  #cmd_args = dbname + " --pset pager=off -F ' ' -A -X -t -c \"select backendname from dmc_backendobject where backendmetaname='%s'\";" % (dmc_obj_name)
  #logDebug("DEBUG CMD: %s" % cmd_args)
  #retcode = subprocess.call(cmd, shell=True)
  #logDebug("DEBUG: retcode=%d" % retcode)

  p = subprocess.Popen([psql_bin,
                        dbname,
                        "--pset",
                        "pager=off",
                        "-F ' '",
                        "-A",
                        "-X",
                        "-t",
                        "-c select backendname from dmc_backendobject where backendmetaname='%s';" % (dmc_obj_name)
                        ],
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE)
  out, err = p.communicate()
  #logDebug("OUT: %s" % out)
  #logDebug("ERR: %s" % err)

  # Check for error (but do not use the return code since useless for this command that always return 0)
  if (len(out.strip()) == 0) or (len(err.strip()) != 0):
    raise Exception("Error while retrieving the 'backendname' for object '%s'" % dmc_obj_name)

  return out.strip()

#-----------------------------------------------------------------------------

# Same as psqlRequestByShell, but almost 2x faster!
def psqlRequestByLib(requests_list):
  dbinfos = ""

  # Connecting to database...
  connection = psycopg2.connect("dbname=%s" % dbname)
  curs = connection.cursor()

  for request in requests_list:
    curs.execute(request[1], request[2])
    request_result = curs.fetchone()
    if (request_result == None):
      raise Exception("ERROR: Unable to retrieve '%s' parameter for object '%s'" % (request[0], request[2][0]))
    logDebug("Request for '%s'" % request[0])
    logDebug("> command = '%s'" % request[1])
    logDebug("> param   = '%s'" % request[2])
    logDebug("> result  = '%s'" % str(request_result[0]))
    dbinfos += "\n%s : %s" % (request[0], str(request_result[0]))

  curs.close()
  connection.close()

  return dbinfos

#-----------------------------------------------------------------------------

def retrieveInfoFromDMCForObject(dmc_obj_name):
  log("    - Retrieving info...")

  metadata = "" # this variable will contains all the data retrieve

  # 1) Retrieve standard info (those from "dmctool getobjectinfo")
  p = subprocess.Popen(["%s getobjectinfo %s" % (dmctool_bin, dmc_obj_name)],
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        shell=True)
  out, err = p.communicate()
  #logDebug("OUT: '%s'" % out)
  #logDebug("ERR: '%s'" % err)

  # Check for error (but do not use the return code since useless for this command that always return 0)
  if (len(out.strip()) == 0) or (len(err.strip()) != 0):
    raise Exception("Error while retrieving the 'dmctool getobjectinfo' for object '%s'" % dmc_obj_name)

  metadata += out.strip()


  # 2) Retrieve several infos directly from database using sql
  # For some request we need to know the object type...
  objectTypeStr = out.strip().split("\n",1)[0].split(":")[1].strip().lower() # ex: "map"
  dmcTableName = DMC_TABLE_OBJECT_NAME_GEN % objectTypeStr
  requests_list = [ ("Backendname", "SELECT backendname from dmc_backendobject where backendmetaname=%s", [dmc_obj_name]),
               ("Flagchunksize", "select flagchunksize from " + dmcTableName + " where backendname=%s", [dmc_obj_name]),
               ("Iooffset", "select iooffset from " + dmcTableName + " where backendname=%s", [dmc_obj_name])
  ]
  
  # If everything goes well just add the new info to the metadata
  metadata += psqlRequestByLib(requests_list)

  return metadata

#-----------------------------------------------------------------------------

def saveMetadataForObject(dmc_obj_name, metadata_str):

  metadata_fullpath = "%s/%s" % (dmc_obj_name, metadata_filename)
  log("    - Saving metadata...")

  # 1) Create the file (overwriting if previously exist)
  cmd = "%s %s --pset pager=off -F ' ' -A -X -t -c \"select PIOPLFCreate('%s',777);\"" % (psql_bin, dbname, metadata_fullpath)
  #logDebug("CMD:\n   %s" % cmd)
  p = subprocess.Popen([cmd],
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        shell=True)
  out, err = p.communicate()
  
  #logDebug("OUT: '%s'" % out)
  #logDebug("ERR: '%s'" % err)

  # Check for error (but do not use the return code since useless for this command that always return 0)
  # Note: in case of success the command return "3" (as output) !
#  if (out.strip() != "3") or (len(err.strip()) != 0):
  if not os.path.exists( metadata_fullpath):
    raise Exception("Error while creating metadata file: '%s'" % metadata_fullpath)

  # 2) Save content to the file
  with open(metadata_fullpath, "w") as text_file:
    text_file.write(metadata_str)

#-----------------------------------------------------------------------------

def getAllDMCObjectForGroup(dmc_group):
  p = subprocess.Popen(["%s getobjectlist --native --nosoftdeleted %s/%s" % (dmctool_bin,dbpath,dmc_group)],
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        shell=True)
  out, err = p.communicate()
  #logDebug("OUT: '%s'" % out)
  #logDebug("ERR: '%s'" % err)

  # Check for error (but do not use the return code since useless for this command that always return 0)
  if (len(out.strip()) == 0) or (len(err.strip()) != 0):
    raise Exception("Error while retrieving list of object for group: %s" % dmc_group)

  return out.strip().split("\n")


#-----------------------------------------------------------------------------

def getAllDMCGroup():
  p = subprocess.Popen(["%s getgrouplist %s" % (dmctool_bin,dbpath)],
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        shell=True)
  out, err = p.communicate()
  #logDebug("OUT: '%s'" % out)
  #logDebug("ERR: '%s'" % err)

  # Check for error (but do not use the return code since useless for this command that always return 0)
  if (len(out.strip()) == 0) or (len(err.strip()) != 0):
    raise Exception("Error while retrieving all DMC groups")

  return out.strip().split("\n")

#-----------------------------------------------------------------------------

def _countNumberOfObjectInDMC():
  mysum = 0
  for ind, grp in enumerate(getAllDMCGroup()):
    p = subprocess.Popen("find %s -type d | wc -l" % (dbpath+"/"+grp), stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    out, err = p.communicate()
    mysum += int(out)
    print("%d   MYSUM = %d" % (ind, mysum))
  print("FINAL SUM = %d" % mysum)

#-----------------------------------------------------------------------------

def printUsage(cmd):
  print("Usage: %s -h | [-f|-F] filename | [-o|-O] objectName | [-g|-G] group | [-a|-A]" % cmd)
  print("  -h            : this help message")
  print("  -f filename   : apply to the DMC objects listed in filename")
  print("  -F filename   : same as '-f' but force re-generation of metadata")
  print("  -o objectName : apply to the specified DMC object")
  print("  -O objectName : same as '-o' but force re-generation of metadata")
  print("  -g group      : apply to all objects in the specified DMC group")
  print("  -G group      : same as '-g' but force re-generation of metadata")
  print("  -a            : apply to all groups and object (CAUTION: this is a very long process...")
  print("  -A            : same as '-a' but force re-generation of metadata")
  print("  -YES          : answer YES to all prompt mean that there is no more user confirmation!\n\
                  CARREFUL: this option must be the first option!")

#-----------------------------------------------------------------------------

def getDMCObjectMD5SUM( metadata_str):
  # computes the md5sum of each .pio and .pio.flag file and then the md5sum of the sorted list of md5sums
  backendname = None
  mdlist = metadata_str.split("\n")
  for mdline in mdlist:
    if mdline.startswith("Backendname"):
      backendname = mdline.strip().split()[-1]
      break
  if backendname == None:
    raise Exception( "can't find Backendname in:\n" + metadata_str)
  cmd = "find %s -type f -exec md5sum {} + | awk '{print $1}' | sort | md5sum" % backendname
  process = subprocess.Popen( cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  (stdout, stderr) = process.communicate() # wait for subprocess to end and get its ouputs
  if (process.returncode != 0):
    raise Exception( "error when computing md5sum:\n" + stderr)
  md5sum = stdout.split()[0]
  assert len( md5sum) == 32, "unexpected result in md5sum: %s" % toto
  return md5sum
  

#-----------------------------------------------------------------------------

def applyMetadataForObject(fullpath_to_dmc_obj, force=False):
  """fullpath_to_dmc_obj name must be a fullpath to DMC object on which we want to apply the metadata.

  By default 'force' is set to False so that metadata are not re-generated if already existing with the good version.
  """

  print( fullpath_to_dmc_obj)

  # Check that metadata does not already exist (with the correct version)
  metadata_fullpath = "%s/%s" % (fullpath_to_dmc_obj, metadata_filename)
  if os.path.isfile(metadata_fullpath):
    # Metadata file already exist
    # Now check that the version is appropriate
    if "%s : %s" % (NODMCMETADATA_VERSION_TAG, VERSION) in open(metadata_fullpath).read():
      log("    Metadata already exist and has a valid version number...")
      if force:
        log("    |--> User ask to force regeneration!")
      else:
        log("    |--> Skip")
        return
    else:
      log("    Metadata exist BUT has an old version number...")
      log("    |--> Will be regenerated!")
  
  # Generate the metadata
  metadata_str = metadata_prefix + "\n"
  metadata_str += retrieveInfoFromDMCForObject(fullpath_to_dmc_obj) + "\n"
  metadata_str += "md5sum : " + getDMCObjectMD5SUM( metadata_str) + "\n"

  #logDebug("")
  #logDebug("#"*30)
  #logDebug("Metadata for '%s':" % fullpath_to_dmc_obj)
  #logDebug(metadata_str)
  
  saveMetadataForObject(fullpath_to_dmc_obj, metadata_str)

#-----------------------------------------------------------------------------

def applyMetadataForGroupSet(groupSet, force=False):
  """groupSet must be a group list, even if it contains only one group.

  By default 'force' is set to False so that metadata are not re-generated if already existing with the good version.
  """
  objectInError = []
  nbObjProcessed = 0

  for dmc_group in groupSet:
    log("Processing DMC group: %s" % dmc_group)

    # Loop on all objects of the current group
    for dmc_obj in getAllDMCObjectForGroup(dmc_group):
      fullpath_to_dmc_obj = dbpath + "/" + dmc_group + "/" + dmc_obj

      try:
        applyMetadataForObject(fullpath_to_dmc_obj, force)
        nbObjProcessed += 1
      except Exception:
        log("    ---> Object SKIP due to error")
        objectInError.append(fullpath_to_dmc_obj)

  # At end of process, print the list of object in error (if any)
  log("")
  if len(objectInError) == 0:
    log(">>> SUCCESS <<< (All %d objects process successfully)" % nbObjProcessed)
  else:
    log(">>> %d object(s) proceed successfully" % nbObjProcessed)
    log(">>> WARNING: %d object(s) were not proceed <<< (list below)" % len(objectInError))
    for obj_in_error in objectInError:
      log("%s" % obj_in_error)

#-----------------------------------------------------------------------------

def applyMetadataForFileList(filename, force=False):
  """filename must be a fullpath to a text file containing a list of DMC
  objects on which we want to apply the metadata.

  By default 'force' is set to False so that metadata are not re-generated if
  already existing with the good version.
  """

  objectInError = []
  nbObjProcessed = 0

  log("Processing file: \"%s\"" % filename)

  with open(filename) as listObj:
    for obj in listObj:
      # Remove any extra space
      obj = obj.strip()

      # Ignore empty lines and lines starting with '#'
      if obj == '' or obj.startswith('#'):
        continue

      try:
        applyMetadataForObject(obj, force)
        nbObjProcessed += 1
      except Exception:
        log("    ---> Object SKIP due to error")
        objectInError.append(obj)

  # At end of process, print the list of object in error (if any)
  log("")
  if len(objectInError) == 0:
    log(">>> SUCCESS <<< (All %d objects process successfully)" % nbObjProcessed)
  else:
    log(">>> %d object(s) proceed successfully" % nbObjProcessed)
    log(">>> WARNING: %d object(s) were not proceed <<< (list below)" % len(objectInError))
    for obj_in_error in objectInError:
      log("    - \"%s\"" % obj_in_error)

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

def multiprocess(force=False):
  """THIS IS A BETA VERSION!
  """
  from multiprocessing import Pool

  pool = Pool(processes=32)
  pool.map(applyMetadataForGroupSet, [getAllDMCGroup()], force)
  pool.close()
  pool.join()
  log("All jobs done!")

#-----------------------------------------------------------------------------


##############################################################################
# MAIN
# This code allow to execute the module as a standalone script.
##############################################################################
if __name__ == "__main__":
  import sys

  log("")
  log("*******************************************************************************")
  log("* %s" % sys.argv[0])
  log("* This is the metadata generator for the NO DMC LIB")
  log("*******************************************************************************")
  log("")


  ### Start parsing command line

  if len(sys.argv) == 1 or (len(sys.argv) == 2 and sys.argv[1] == "-h"):
    printUsage(sys.argv[0])
    exit(1)

  metadata_prefix = getMetadataIdentifierContent()

  # Handle the '-YES' option
  YES_answer = False
  argpos = 1
  if sys.argv[argpos] == "-YES":
    if len(sys.argv) > 2:
      YES_answer = True
      argpos = 2
    else:
      print("-YES option must be associated with another option!")
      exit(1)

  if sys.argv[argpos] == "-f":
    assert(len(sys.argv) == argpos+2)
    applyMetadataForFileList(sys.argv[argpos+1])
  elif sys.argv[argpos] == "-F":
    assert(len(sys.argv) == argpos+2)
    applyMetadataForFileList(sys.argv[argpos+1], force=True)
  elif sys.argv[argpos] == "-o":
#    assert(len(sys.argv) == argpos+2)
    for objidx in range( argpos+1, len(sys.argv)):
      applyMetadataForObject(sys.argv[objidx])
  elif sys.argv[argpos] == "-O":
#    assert(len(sys.argv) == argpos+2)
    for objidx in range( argpos+1, len(sys.argv)):
      applyMetadataForObject(sys.argv[objidx], force=True)
  elif sys.argv[argpos] == "-g":
    assert(len(sys.argv) == argpos+2)
    applyMetadataForGroupSet([sys.argv[argpos+1]])
  elif sys.argv[argpos] == "-G":
    assert(len(sys.argv) == argpos+2)
    print("!!CAUTION!!: This option will overwrite ALL objects metadata in group %s" % sys.argv[argpos+1])
    if YES_answer or query_yes_no("Are you sure you want to proceed?", default="no"):
      applyMetadataForGroupSet([sys.argv[argpos+1]], force=True)
    else:
      print("Nothing done!")
      exit(0)
  elif sys.argv[argpos] == "-a":
    assert(len(sys.argv) == argpos+1)
    print("!!CAUTION!!: This option will process ALL objects of database %s" % dbpath)
    if YES_answer or query_yes_no("Are you sure you want to proceed?", default="no"):
      multiprocess()
    else:
      print("Nothing done!")
      exit(0)
  elif sys.argv[argpos] == "-A":
    assert(len(sys.argv) == argpos+1)
    print("!!CAUTION!!: This option will process ALL objects of database %s, and overwrite any pre-existing metadata" % dbpath)
    if YES_answer or query_yes_no("Are you sure you want to proceed?", default="no"):
      multiprocess(force=True)
    else:
      print("Nothing done!")
      exit(0)
  else:
    print("'%s' option NOT supported! (try -h for detail about available options)" % sys.argv[argpos])
    exit(1)

  # Successfull end
  print("Done.")
