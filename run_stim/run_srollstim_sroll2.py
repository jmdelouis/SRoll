#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, sys
import subprocess
import time, calendar
import math
import numpy
import copy
import stools
from   IMO_4_27 import * # BOLOID, DETSETS, CALIB, NEP, XPOL, ELECWNOISE, PHOTWNOISE


BEGMISS = 240
ENDMISS = 26050
BEGHM1  = BEGMISS
ENDHM1  = 13143
BEGHM2  = ENDHM1+1
ENDHM2  = ENDMISS


# constants from sroll.c for MAPRINGS parameter
FULL        = 1
HM12        = 2
S12345      = 4
YEAR12      = 8
FULLODDEVEN = 16
HM12ODDEVEN = 32


try:
  SROLLHOST = os.environ["SROLLHOST"]
  SROLLDIR  = os.environ["SROLLDIR"]
  NODMCDATA = os.environ["NODMCDATA"]
  USER      = os.environ["USER"]
except:
  raise Exception( "You must 'source ./srollex_setenv.sh' before using this script")


################################################################################
# SrollParam: method-less class to hold all parameters needed to write
# sroll/stim parameter files and submission scripts

class SrollParam( object):

  def __init__( self, detset, srollver="s2"):

    # "s1": sroll1 executable (HFI2018 data / FFP10 sims)
    # "s2": sroll2 executable (bware/SRoll2 data / JAN18 sims)
    # "s2s1": sroll2 executable with sroll1 parameters
    # "s3": troll executable
    assert srollver in ("s1", "s2", "s2s1", "s3")

    # init environment variables and paths
    self.srollhost      = SROLLHOST
    self.srolldir       = SROLLDIR
    self.user           = USER
    self.scriptfilename = os.path.realpath( __file__)
    self.directory      = os.path.dirname( self.scriptfilename) + "/" # run_srollstim.py directory

    if self.srollhost == "M3":
      self.dbpath       = NODMCDATA # /data/dmc/MISS03/DATA
      self.brimo_path   = "/wrk/symottet/brimo_4_27"
      self.scratch      = "/redtruck/SimuData"
      self.batch_suffix = ".qsub"

    elif self.srollhost == "M4":
      self.dbpath       = NODMCDATA # /pscratch1/RD12_data/dmc/MISS03/DATA
      self.brimo_path   = "/wrk/symottet/brimo_4_27"
      self.scratch      = self.dbpath
      self.batch_suffix = ".qsub"

    elif self.srollhost == "EDISON":
      self.dbpath       = NODMCDATA # /scratch3/scratchdirs/paganol/RD12_data/dmc/MISS03/DATA
      self.brimo_path   = "/project/projectdirs/planck/data/hfi/RD12_data/brimo/brimo_4_27"
      self.scratch      = os.environ["SCRATCH"] # /scratch*/scratchdirs/$USER
      self.batch_suffix = ".slurm"

    elif self.srollhost == "CORI":
      self.dbpath       = NODMCDATA # /project/projectdirs/planck/data/hfi/RD12_data/dmc/MISS03/DATA
      self.brimo_path   = "/project/projectdirs/planck/data/hfi/RD12_data/brimo/brimo_4_27"
      self.scratch      = os.environ["SCRATCH"] # /scratch*/scratchdirs/$USER
      self.batch_suffix = ".slurm"

    elif self.srollhost == "OCCIGEN":
      self.dbpath       = NODMCDATA # /scratch/cnt0028/ias1717/SHARED/RD12_data/dmc/MISS03/DATA
      self.brimo_path   = "/scratch/cnt0028/ias1717/SHARED/RD12_data/brimo_4_27"
      self.scratch      = os.environ["SCRATCHDIR"] # /scratch/cnt0028/ias1717/$USER
      self.batch_suffix = ".slurm"      
    else:
      raise Exception( "Unknown platform: <%s>" % self.srollhost)

    # default parameters values
    self.bolometers     = None # list of single detector names
    self.freq           = detset[0:3]
    self.begring        = BEGMISS
    self.endring        = ENDMISS
    self.run_name       = "default_run_name"
    self.output_prefix  = "runsrollstim"
    self.halfring       = "" # user's halfring tag, eg. "1st" or "HR1", etc.
    self.RD12calib      = None # set to True to use RD12 calib factors, to False for DX11 calib and None to try to guess

    # batch parameters
    self.joblabel       = None
    self.parfile_name   = None
    self.nodes          = None
    self.cancel_run     = False
    self.walltime       = None  # decimal number of hours (eg: 1.5) to override the value computed by the script

    # stim parameters
    self.iterations     = 1
    self.random_seed    = 0
    self.rings_per_read = -1
    self.input_toi      = None
    self.signal_out     = 0
    self.hpr_out        = 0
    self.add_noise      = 1
    self.despike_flag   = None
    self.conv_deconv    = None
    self.sim_adu        = 1
    self.sim_inl_name   = None
    self.corr_inl_name  = None
    self.fourk_lines    = 1 # 1 for APR20, 0 for JAN18
    self.save_sim_inl   = False
    self.add_oof        = None
    self.add_4Kres      = 0 # 0 for APR20, None (auto) for JAN18
    self.stim_calib     = 0
    self.save_toi       = False
    self.do_SHPR        = 1 # 1 for APR20 and SRoll3, 0 for JAN18 and SRoll1 and SRoll2.x
    self.save_hpr       = False

    # sroll parameters
    self.srollver       = srollver # "s2": sroll2, "s1": sroll1, "s2s1": sroll2 running with sroll1 parameters
    self.dipole         = None # set to "RD12" for RD12/FFP10, or "HFI17" for 2018+
    self.rstep          = 1
    self.gainstep       = 1
    self.nadu           = 32
    self.nadustep       = 1
    self.n_sroll_iter   = 20
    self.seuilcond      = 3e-4
    self.ffp10_maps     = False
    self.map_per_bolo   = False # requires self.ffp10_maps = True
    self.all_surveys    = False
    self.map_detsets    = None # list of detsets for which to produce a map
    self.map_ringsets   = None # list of ringsets (ringcuts) for each map_detsets
    self.input_hpr      = None
    self.theonops       = True # set to False to not fit the transfer function with Theo_noPS parmeters
    self.fslcoef        = None # set to 1 to remove RD12 FSL / RD12 reproduction, 0 to force no removal
    self.kcmbin         = 0
    self.adddip         = 0
    self.fitangle       = 1
    self.fitpoleff      = 1
    self.saveCOV        = 1
    self.delta_psi      = 0.0
    self.verbose        = 0

    # frequency dependant parameters
    assert self.freq in ("100", "143", "217", "353", "545", "857")
    if self.freq == "100":
      self.seuilcond = 3e-4
      self.gainstep  = 1
      self.nadu      = 32
      self.nadustep  = 1
    elif self.freq == "143":
      self.seuilcond = 3e-4
      self.gainstep  = 1
      self.nadu      = 8
      self.nadustep  = 128
    elif self.freq == "217":
      self.seuilcond = 3e-4
      self.gainstep  = 1
      self.nadu      = 32
      self.nadustep  = 1
    elif self.freq == "353":
      self.seuilcond = 3e-3
      self.gainstep  = 1
      self.nadu      = 32
      self.nadustep  = 1
    elif self.freq in ("545", "857"):
      self.seuilcond = 1e-6
      self.gainstep  = 32
      self.nadu      = 32
      self.nadustep  = 1

    # sroll 1 reproduction values
    if self.srollver in ("s1", "s2s1"):
      self.nadu         = 1
      self.nadustep     = 1
      self.n_sroll_iter = 7
      self.fitangle     = 1
      self.fitpoleff    = 1
      if self.freq <= "217":
        self.gainstep = 128
      else:
        self.gainstep = 32

  def checkinit( self):
    assert len( self.bolometers) > 0
    if self.map_detsets != None and self.map_ringsets != None:
      assert len( self.map_detsets) == len( self.map_ringsets)

    self.freq = self.bolometers[0][0:3]
    for b in self.bolometers[1:]:
      if b[0:3] != self.freq:
        self.freq = None # multifrequency detsets
        break

    if self.output_prefix[0] != "/":
      self.output_prefix = self.scratch + "/" + self.output_prefix

    if self.joblabel == None:
      self.joblabel = self.run_name
    self.joblabel += self.halfring

    # parameter files are created in the run_srollstim.py script directory
    self.parfile_name = self.directory + self.joblabel
    while os.path.isfile( self.parfile_name + self.batch_suffix):
      # if parameter files already exist, add date and time to prefix to build a unique one
      time.sleep(1)
      lt = time.localtime()
      self.parfile_name = self.directory \
                        + "%d%s%02d_%02d%02d%02d_" % (lt.tm_year, calendar.month_name[lt.tm_mon][:3], lt.tm_mday, lt.tm_hour, lt.tm_min, lt.tm_sec) \
                        + self.joblabel

    # create output directories if needed
    os.system( "mkdir -p %s_MAP" % self.output_prefix)
    os.system( "mkdir -p %s_VEC" % self.output_prefix)


  def todict( self, toto):
    # return a dictionary containing the <toto> dictionary and self variables, for format()
    return dict( self.__dict__.items() + toto.items())


################################################################################

# check that all input objects and files in a parameter file actually exist
def check_param_inputs( parfile):
  ignore_keys = ["MAP", "Out_", "signal_out", "hpr_out", "save_sim_inl"]
  checkok = True
  parlines = open( parfile).readlines()
  for line in parlines:
    line = line.strip()
    if (len(line) == 0) or (line[0] == "#"):
      continue
    equalpos = line.find( "=")
    if equalpos != -1:
      key = line[:equalpos].strip()
      val = line[equalpos+1:].strip()
      for ignore in ignore_keys:
        if key.startswith( ignore):
          val = ""
      if (len(val) == 0) or (val[0] != "/"):
        continue
      if os.path.exists( val):
        if os.path.isdir( val):
          if not os.path.exists( val + "/no_dmc_metadata.txt"):
            print "\n**WARNING - file not found: '" + line + "/no_dmc_metadata.txt'" + "\n"
            checkok = False
      else:
        print "\n** WARNING ** file not found: '" + line + "'\n"
        checkok = False
  return checkok


################################################################################

def get_subsets( detset, noswb = False):
  freq = detset[0:3]
  if (freq == "545") or (freq == "857"):
    return [detset]
  if detset in ["143ghz", "217ghz", "353ghz"]:
    if noswb:
      return [detset, freq+"ds1", freq+"ds2", freq+"psb"]
    else:
      return [detset, freq+"ds1", freq+"ds2", freq+"psb", freq+"swb"]
  if detset.endswith("psb") or (detset == "100ghz"):
    return [detset, freq+"ds1", freq+"ds2"]
  else:
    return [detset]


################################################################################

def make_submission_script( param, exe_name, submit):

# full frequency sroll requires about 2.4TB RAM (2TB on 512 procs on M3)
# M3 nodes have 2x4 cores, 32GB RAM (4GB/core), no hyperthreading (Intel Xeon E5472 3.00GHz PassMark: 4254)
# M4 nodes have 2x12 cores, 128GB RAM (5.3GB/core), x2 hyperthreading (Intel Xeon E5-2680v3 2.5GHz PassMark: 18775)
# Edison nodes have 2x12 cores, 64GB RAM (2.6GB/core), x2 hyperthreading (Intel Ivy Bridge Xeon E5-2695v2 2.4GHz PassMark: 15708)
# Cori2 nodes have 1x68 cores, 96GB RAM (1.4GB/core), x4 hyperthreading (Intel Xeon Phi 7250 1.4GHz)
# CINES Occigen Haswell Xeon E5-2690V3 @ 2.6GHz PassMark: 19333
# CINES Occigen Broadwell Xeon E5-2690V4 @ 2.6GHz PassMark: 21801

  lpar = copy.copy( param) # local working copy of parameters structure
  batch_file  = lpar.parfile_name + lpar.batch_suffix
  logfilename = lpar.parfile_name + ".log"
  merge_fits = "{directory}merge_fits_maps.py {output_prefix}_MAP/{run_name}*I.fits >>{logfilename} 2>&1".format( **lpar.todict( locals()))

  if lpar.srollhost == "M3":
    ppn = 8    # processors per node
    rpn = 32.0 # ram per node
  elif lpar.srollhost == "M4":
    ppn = 24
    rpn = 128.0
  elif lpar.srollhost == "EDISON":
    ppn = 24
    rpn = 64.0
  elif lpar.srollhost == "CORI":
    ppn = 68
    rpn = 96.0
  elif lpar.srollhost == "OCCIGEN":
    ppn = 24
    rpn = 64.0

  if exe_name.startswith( "sroll"):
    if lpar.srollver == "s1":
      exe_dir = "sroll1"
    else:
      exe_dir = "sroll"
    exe_parfile = lpar.parfile_name + ".sroll.par"

    if lpar.iterations > 1:
      if lpar.freq >= "545":
        exe_name = "srollit545"
      else:
        exe_name = "srollit"

    if lpar.nodes == None:
      # default number of nodes for 8-12 bolometers
      if lpar.srollhost == "M3":
        lpar.nodes = 64
        if lpar.freq in ["217", "353"]:
          lpar.nodes = 74
      elif lpar.srollhost == "M4":
        lpar.nodes = 22
      elif lpar.srollhost == "EDISON":
        lpar.nodes = 43
      elif lpar.srollhost == "CORI":
        lpar.nodes = 24
      elif lpar.srollhost == "OCCIGEN":
#        lpar.nodes = 22  # 512 cores
        lpar.nodes = 43  # 1024 cores
#        lpar.nodes = 86  # 2048 cores
#        lpar.nodes = 171 # 4096 cores

      if len( lpar.bolometers) <= 4:
        lpar.nodes /= 2
      elif len( lpar.bolometers) > 12:
        lpar.nodes *= 2

    # sroll mpi size must always be a power of 2
    # get the power of two that fits in the given number of nodes
    mpisize = 2 ** int( math.log( int( lpar.nodes) * ppn, 2))

  elif exe_name == "stimexe":
    exe_dir = "stim"
    exe_parfile = lpar.parfile_name + ".stim.par"
    if lpar.nodes == None:
      lpar.nodes = 1
      mpisize = lpar.nodes * ppn
    else:
      # if nodes is given for stimexe, apply the power of 2 rule fo mpi size
      mpisize = 2 ** int( math.log( int( lpar.nodes) * ppn, 2))

  elif exe_name == "troll":
    exe_dir = "sroll"
    exe_parfile = lpar.parfile_name + ".troll.par"
    if lpar.nodes == None:
      # default number of nodes for 8-12 bolometers
      if lpar.srollhost == "M4":
        lpar.nodes = 22
      elif lpar.srollhost == "OCCIGEN":
        lpar.nodes = 43  # 1024 cores
    mpisize = 2 ** int( math.log( int( lpar.nodes) * ppn, 2))

  else:
    raise Exception( "Unknown exe_name: " + exe_name)


#-----------------------------------------------------------------------

  if lpar.srollhost == "M3":
    lpar.joblabel = lpar.joblabel[0:15]
    if exe_name == "stimexe":
      # 15 min/iteration on 128 nodes (1024 cores)
      wallhours = 1.0 + 0.25 * lpar.iterations * 128.0 / lpar.nodes
    else:
      wallhours = 6 * (lpar.iterations+1)
    if (lpar.nodes > 64) and (wallhours > 48): # paralT queue
      wallhours = 48
    if lpar.walltime != None:
      wallhours = lpar.walltime
    wallhours = int( math.ceil( wallhours))

    batch_script = """#!/bin/sh
#PBS -S /bin/sh
#PBS -N {joblabel}
#PBS -j oe
#PBS -l select={nodes}:ncpus={ppn},walltime={wallhours}:00:00

cd {srolldir}
source ./srollex_setenv.sh

export OMP_NUM_THREADS=1

cd {exe_dir}
mpiexec -np {mpisize} -x TMPDIR=/tmp --bynode ./{exe_name} {exe_parfile} &> {logfilename}

# enable the following line to automatically merge fits IQU after sroll. IT CONSUMES CPU TIME
#mpiexec -np {nodes} -x TMPDIR=/tmp --bynode {merge_fits}

exit 0
""".format( **lpar.todict( locals()))

    f = open( batch_file, "wt")
    f.write( batch_script)
    f.close()
    shell_cmd = "qsub " + batch_file

#-----------------------------------------------------------------------

  if lpar.srollhost == "M4":
    wallhours = 4 * (lpar.iterations+1)
    if lpar.walltime != None:
      wallhours = lpar.walltime
    wallhours = int( math.ceil( wallhours))

    batch_script = """#!/bin/sh
#PBS -S /bin/sh
#PBS -N {joblabel}
#PBS -j oe
#PBS -l nodes={nodes}:ppn={ppn},walltime={wallhours}:00:00

cd {srolldir}
source ./srollex_setenv.sh

export OMP_NUM_THREADS=1

cd {exe_dir}
mpiexec -np {mpisize} -x TMPDIR=/tmp --bynode ./{exe_name} {exe_parfile} &> {logfilename}

# enable the following line to automatically merge fits IQU after sroll. IT CONSUMES CPU TIME
#{merge_fits}

exit 0
""".format( **lpar.todict( locals()))

    f = open( batch_file, "wt")
    f.write( batch_script)
    f.close()
    shell_cmd = "qsub " + batch_file

#-----------------------------------------------------------------------

  if lpar.srollhost == "EDISON":
    nersc_partition = NERSC_PART
    qos_planck = ""
    if nersc_partition == "debug":
      wallhours = 0.5
    elif lpar.iterations == 5:
      if lpar.freq == "100":
        wallhours = 13.0
      elif lpar.freq == "143":
        wallhours = 17.0
      elif lpar.freq in ("217", "353") or (lpar.freq == None):
        wallhours = 20.0
      elif lpar.freq in ("545", "857"):
        wallhours = 3.0
    else:
      if lpar.freq == "100":
        wallhours = 4.5 + 2.5 * (lpar.iterations-1)
      elif lpar.freq == "143":
        wallhours = 4.0 + 3.5 * (lpar.iterations-1)
      elif lpar.freq in ("217", "353") or (lpar.freq == None):
        wallhours = 5.5 + 4.0 * (lpar.iterations-1)
      elif lpar.freq in ("545", "857"):
        wallhours = 0.5 * lpar.iterations

      if QOS_PLANCK:
        if wallhours > 12:
          wallhours = 12
        qos_planck = "#SBATCH --qos=special_planck"

    if lpar.walltime != None:
      wallhours = lpar.walltime
    walldec, wallint = math.modf( wallhours)
    walltime = "%d:%d:00" % ( wallint, walldec*60)

    batch_script = """#!/bin/bash -l
#SBATCH -p {nersc_partition}
#SBATCH -N {nodes}
#SBATCH -t {walltime}
#SBATCH -J {joblabel}
{qos_planck}

cd {srolldir}
source ./srollex_setenv.sh

# distribute mpi ranks by node
export MPICH_RANK_REORDER_METHOD=0

cd {exe_dir}
srun -n {mpisize} ./{exe_name} {exe_parfile} &> {logfilename}

# enable the following line to automatically merge fits IQU after sroll. IT CONSUMES CPU TIME
#{merge_fits}

exit 0
""".format( **lpar.todict( locals()))

    f = open( batch_file, "wt")
    f.write( batch_script)
    f.close()
    shell_cmd = "sbatch " + batch_file

#-----------------------------------------------------------------------

  if lpar.srollhost == "CORI":
    nersc_partition = NERSC_PART
    if nersc_partition == "debug":
      wallhours = 0.5
    else:
      wallhours = 2.0
    if lpar.walltime != None:
      wallhours = lpar.walltime
    walldec, wallint = math.modf( wallhours)
    walltime = "%d:%d:00" % ( wallint, walldec*60)
    ompthreads = 4

    batch_script = """#!/bin/bash -l
#SBATCH -p {nersc_partition}
#SBATCH -N {nodes}
#SBATCH -t {walltime}
#SBATCH -J {joblabel}
#SBATCH -L SCRATCH,project
#SBATCH -C knl,quad,cache

cd {srolldir}
source ./srollex_setenv.sh

export OMP_NUM_THREADS={ompthreads}

cd {exe_dir}
srun -n {mpisize} -c {ompthreads} --cpu_bind=cores ./{exe_name} {exe_parfile} &> {logfilename}

# enable the following line to automatically merge fits IQU after sroll. IT CONSUMES CPU TIME
#{merge_fits}

exit 0
""".format( **lpar.todict( locals()))

    f = open( batch_file, "wt")
    f.write( batch_script)
    f.close()
    shell_cmd = "sbatch " + batch_file

#-----------------------------------------------------------------------

  if lpar.srollhost == "OCCIGEN":
    memreq = "50GB"
    if exe_name == "stimexe":
      # 12-15+ min/iteration on 43 nodes (1024 cores)
      wallhours = 0.3 * (lpar.iterations) * 43.0 / lpar.nodes
    else:
      if lpar.nadustep >= 32:
        memreq = "100GB"
      wallhours = 2.0 * (lpar.iterations)
      if len( lpar.bolometers) > 12:
        wallhours *= 2.0
      wallhours *= 1.0 + lpar.nadustep/64.0 # 143ghz sroll2 nadustep=128, requires 2x to 3x more walltime
    if lpar.walltime != None:
      wallhours = lpar.walltime
    if wallhours < 0.5:
      wallhours = 0.5
    walldec, wallint = math.modf( wallhours)
    walltime = "%d:%d:00" % ( wallint, walldec*60)

    batch_script = """#!/bin/bash -l
#SBATCH -J {joblabel}
#SBATCH --nodes={nodes}
#SBATCH --constraint=HSW24
#SBATCH --ntasks={mpisize}
#SBATCH --ntasks-per-node=24
#SBATCH --threads-per-core=1
#SBATCH --mem={memreq}
#SBATCH --time={walltime}

cd {srolldir}/{exe_dir}
# source ./srollex_setenv.sh is not needed because environment variable are propagated to the job
# use Mellanox reliable connection http://www.mellanox.com/related-docs/prod_software/Mellanox_MXM_README_v2.0.pdf
export MXM_TLS=self,shm,rc
srun --mpi=pmi2 -K1 -n $SLURM_NTASKS --distribution=cyclic ./{exe_name} {exe_parfile} &> {logfilename}

# enable the following line to automatically merge fits IQU after sroll. IT CONSUMES CPU TIME
#srun --mpi=pmi2 -K1 -n {nodes} --distribution=cyclic {merge_fits}

exit 0
""".format( **lpar.todict( locals()))

    f = open( batch_file, "wt")
    f.write( batch_script)
    f.close()
    shell_cmd = "sbatch " + batch_file

#-----------------------------------------------------------------------

  print "submission command: " + shell_cmd

# submit job if requested
  if submit:
    if param.cancel_run:
      print "\nERROR: SOME INPUTS DON'T EXIST, SUBMISSION CANCELED\n"
      print "manually run the above submission command to submit the job"
      return
    process = subprocess.Popen(shell_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (stdout_data, stderr_data) = process.communicate()
    if process.returncode != 0:
      print "\nerror stdout:" + stdout_data
      print "error stderr:" + stderr_data
    else:
      pbs_id = stdout_data.strip().split("\n")[-1].strip()
      print "job id: " + pbs_id
  else:
    print "manually run the above submission command to submit the job"
  print


################################################################################

def make_stim_parfile( param, pixname):

  lpar = copy.copy( param) # local working copy of parameters structure

  boloid   = BOLOID[pixname]
  beltchan = boloid[0:2]
  stim_parfile = lpar.parfile_name + "_%s.stim.par" % pixname

  assert lpar.sim_adu in (0, 1)
  assert lpar.add_noise in (0, 1)

  assert lpar.input_toi != None
#  lpar.input_toi = "{dbpath}/e2e_common_TOI/{pixname}_ffp10_fgdip_cmbfid"
#  lpar.input_toi = "{dbpath}/e2e_common_TOI/{pixname}_dustapr17"
#  lpar.input_toi = "{dbpath}/e2e_common_TOI/{pixname}_JAN18sky"

  # RD12calib is used to adjust white noise levels initially computed with DX11 calib factors
  if lpar.RD12calib == None:
    if "ffp10" in lpar.input_toi.lower():
      lpar.RD12calib = False
    elif "JAN18" in lpar.input_toi:
      lpar.RD12calib = True
    else:
      raise Exception( "Unable to guess RD12calib value from stim input_toi (%s)" % lpar.input_toi)

  if lpar.RD12calib:
    rd12calib = "RD12calib = 1"
  else:
    rd12calib = "RD12calib = 0"

  if lpar.corr_inl_name == None:
    # use DX11 correction INL
    # set to "0" to use a linear INL
    lpar.corr_inl_name = "{dbpath}/RAW_VECT_FIT_INFO/HFI_{beltchan}_DNL_DX11"

  if lpar.sim_inl_name == None:
    # use 2016 simulation INL
    # set to "0" to use a linear INL
    if (lpar.freq >= "353") or not (pixname[-1] in "ab"):
      lpar.sim_inl_name = "{dbpath}/RAW_VECT_FIT_INFO/{pixname}_SIMINL1406_offsets"
    else:
      lpar.sim_inl_name = "{dbpath}/RAW_VECT_FIT_INFO/{pixname}_SIMINL2803_offsets"

  if lpar.save_sim_inl:
    lpar.save_sim_inl = "save_sim_inl = " + lpar.output_prefix + "_VEC/USEDSIMINL"
  else:
    lpar.save_sim_inl = "#save_sim_inl="

  if lpar.save_toi == True:
    if lpar.signal_out == 0:
      lpar.signal_out = lpar.output_prefix + "_VEC/{pixname}_%s%s_stimTOI" % (lpar.run_name, lpar.halfring)
  else:
    lpar.signal_out = 0

  if lpar.save_hpr == True:
    if lpar.hpr_out == 0:
      if lpar.do_SHPR:
        lpar.hpr_out = lpar.output_prefix + "_VEC/{pixname}_%s%s_stimSHPR" % (lpar.run_name, lpar.halfring)
      else:
        lpar.hpr_out = lpar.output_prefix + "_VEC/{pixname}_%s%s_stimHPR" % (lpar.run_name, lpar.halfring)
  else:
    lpar.hpr_out = 0

  if lpar.freq >= "545":
    lpar.sim_adu = 0
    lpar.despike_flag = 0
    lpar.add_oof = 0
    lpar.stim_calib = -1
    lpar.add_4Kres = 0 # {pixname}_REP6_4k_residu exist for 545/857 but were not used for FFP10

  if lpar.despike_flag == None:
    lpar.despike_flag = lpar.add_noise
  if lpar.add_oof == None:
    lpar.add_oof = lpar.add_noise
  if lpar.add_4Kres == None:
    lpar.add_4Kres = lpar.add_noise
  if lpar.conv_deconv == None:
    lpar.conv_deconv = lpar.add_noise

  if lpar.fourk_lines == 1:
    # APR20 sims
    lpar.fourk_name = "{dbpath}/RAW_4K/{beltchan}_HARMS_LM4K_DX11_PIOFLOAT"
  else:
    # JAN18 sims
    lpar.fourk_name = 0

  if lpar.add_4Kres == 1:
    lpar.add_hpr = ""
  else:
    # comment add_hpr_* items when number_of_add_hpr_* is 0
    lpar.add_hpr = "#"

  # halfring is the suffix used in output object names, can be anything containing either 1 or 2 eg. "", "_HR1" or "_HR2".
  # hrtag is the suffix used in input object names, is either "", "_1st" or "_2nd"
  if "1" in lpar.halfring:
    hrtag = "_1st"
  elif "2" in lpar.halfring:
    hrtag = "_2nd"
  else:
    hrtag = ""

  STRNOW = time.strftime("%c")

  stim_parameters = """
# stim parameter file generated by {scriptfilename}
# for {user}
# on {srollhost}
# at {STRNOW}

bolometer       = {pixname}
brimo_filename  = {brimo_path}
stay_in_memory  = -1

iterations      = {iterations}
random_seed     = {random_seed}

BeginRing       = {begring}
EndRing         = {endring}
rings_per_read  = {rings_per_read}

{rd12calib}
signal_in       = {input_toi}
signal_out      = {signal_out}
hpr_out         = {hpr_out}
TOI_pixel_index = {dbpath}/e2e_common_TOI/{pixname}_HPRIDX_ABER{hrtag}_TotalFlag_dx11
phase           = {dbpath}/Sa_HFI_C_Bolo/Phase_dx11


do_photonic_noise = {add_noise}


do_LFER_convolve = {conv_deconv}


do_electronic_noise = {add_noise}


do_despike_flag = {despike_flag}


do_LFER_deconvolve = {conv_deconv}
# 0 = use default values
deconv_LPF_FILTER = 0
deconv_R_FILTER   = 0


do_sim_adu = {sim_adu}
sim_inl_name = {sim_inl_name}
{save_sim_inl}
rawgain_name = {dbpath}/RAW_GAIN/{beltchan}_INL_GP41N
rawcst_name  = {dbpath}/RAW_CST/{beltchan}_LM4K_DX11
fourk_name   = {fourk_name}
add_baseline = {dbpath}/e2e_common_TOI/{pixname}_dsn_baseline_Feb16
electronic_noise_adu = 4.0


do_compress_decompress = {sim_adu}


do_correct_adc = {sim_adu}
corr_inl_name  = {corr_inl_name}


do_adu_to_volts = {sim_adu}
use_bolometer_nonlinearity = 0


do_bl_demod = {sim_adu}


do_gaindecorr    = {sim_adu}
thermal_baseline = {dbpath}/calTOIs/T90_sm10823_v61


nharm_oof     = {add_oof}
amp_oof_param = 0


do_calibrate = {stim_calib}


add_final_toi = 0
add_final_toi_factor = 0

do_SHPR = {do_SHPR}

number_of_add_hpr_name = {add_4Kres}
{add_hpr}add_hpr_name1 = {dbpath}/PBR_JMD/{pixname}_REP6{hrtag}_4k_residu

number_of_add_hpr_factor = {add_4Kres}
{add_hpr}add_hpr_factor1 = 1.0

""".format( **lpar.todict( locals()))

  stim_parameters = stim_parameters.format( **lpar.todict( locals())) # interpret formats in formats
  f = open( stim_parfile, "wt")
  f.write( stim_parameters)
  f.close()
  print "stim parameter file: " + stim_parfile
  if not check_param_inputs( stim_parfile):
    param.cancel_run = True


################################################################################

def make_sroll_parfile( param):

  lpar = copy.copy( param) # local working copy of parameters structure
  sroll_parfile = lpar.parfile_name + ".sroll.par"

  # try to guess version of input sroll data or simulations, for default parameter values
  DATAVER = ""
  if (lpar.input_hpr != None):
    # try to guess DATAVER from input HPR name
    hprname = lpar.input_hpr.split("/")[-1].upper()
    if "REP6" in hprname:
      DATAVER = "REP6" # activate default parameters for 2015 data TOI projection
    elif "JAN18" in hprname:
      DATAVER = "JAN18"
    elif ("DEC16" in hprname) or ("MAR17" in hprname) or ("FFP10" in hprname):
      DATAVER = "DEC16" # FFP10 parameters
    else:
      raise Exception( "Unkown input signal HPR: " + lpar.input_hpr)
  else:
    # try to guess DATAVER from input TOI name
    assert lpar.input_toi != None
    toiname = lpar.input_toi.split("/")[-1].upper()
    if ("JAN18" in toiname):
      DATAVER = "JAN18"
    elif ("DEC16" in toiname) or ("MAR17" in toiname) or ("FFP10" in hprname) or ("DUSTOCT" in hprname) or ("DUSTAPR" in hprname):
      DATAVER = "DEC16" # FFP10 parameters
    else:
      raise Exception( "Unkown input signal TOI: " + lpar.input_toi)

  if lpar.RD12calib == None:
    if DATAVER == "REP6":
      lpar.RD12calib = True
    elif DATAVER == "DEC16":
      lpar.RD12calib = False
    elif DATAVER == "JAN18":
      lpar.RD12calib = True
    else:
      raise Exception( "Unable to guess RD12calib switch value from sroll input_hpr (%s)" % lpar.input_hpr)

  if DATAVER == "DEC16":
    TPL_CO     = "{dbpath}/MAP_JMD_128/COMAP"
    TPL_CO_Q   = "{dbpath}/MAP_JMD_128/CO_ORT_Q"
    TPL_CO_U   = "{dbpath}/MAP_JMD_128/CO_ORT_U"
    TPL_13CO   = ""
    TPL_MAP_I  = ""
    TPL_MAP_Q  = ""
    TPL_MAP_U  = ""
    TPL_FIT_Q  = ""
    TPL_FIT_U  = ""
    TPL_DUST_I = "{dbpath}/MAP_JMD_128/Dust_I"
    TPL_DUST_Q = "{dbpath}/MAP_JMD_128/Dust_Q"
    TPL_DUST_U = "{dbpath}/MAP_JMD_128/Dust_U"

  elif DATAVER == "REP6":
    TPL_CO     = "{dbpath}/sroll_templates/MAP_JMD_128/NEW_cleaned_12CO.float32.bin"
    TPL_CO_Q   = "{dbpath}/MAP_JMD_128/CO_ORT_Q"
    TPL_CO_U   = "{dbpath}/MAP_JMD_128/CO_ORT_U"
    TPL_13CO   = "{dbpath}/sroll_templates/MAP_JMD_128/NEW_cleaned_13CO.float32.bin"
    TPL_MAP_I  = "{dbpath}/sroll_templates/MAP_JMD_2048/{freq}_I.float32.bin"
    TPL_MAP_Q  = "{dbpath}/sroll_templates/MAP_JMD_2048/{freq}_Q.float32.bin"
    TPL_MAP_U  = "{dbpath}/sroll_templates/MAP_JMD_2048/{freq}_U.float32.bin"
    TPL_FIT_Q  = "{dbpath}/sroll_templates/MAP_JMD_128/{freq}_Q.float32.bin"
    TPL_FIT_U  = "{dbpath}/sroll_templates/MAP_JMD_128/{freq}_U.float32.bin"
    TPL_DUST_I = "{dbpath}/MAP_JMD_128/Dust_New_I"
    TPL_DUST_Q = "{dbpath}/MAP_JMD_128/Dust_New_Q"
    TPL_DUST_U = "{dbpath}/MAP_JMD_128/Dust_New_U"

  elif DATAVER == "JAN18":
    TPL_CO     = "{dbpath}/sroll_templates/MAP_JMD_128/COSIMU.float32.bin"
    TPL_CO_Q   = "{dbpath}/MAP_JMD_128/CO_ORT_Q" # ???
    TPL_CO_U   = "{dbpath}/MAP_JMD_128/CO_ORT_U" # ???
    TPL_13CO   = ""
    TPL_MAP_I  = "{dbpath}/sroll_templates/MAP_JMD_2048/template_sims_dustapr17_aligned_{freq}_000_full_1deg_ns2048_v1_I.float32.bin"
    TPL_MAP_Q  = "{dbpath}/sroll_templates/MAP_JMD_2048/template_sims_dustapr17_aligned_{freq}_000_full_1deg_ns2048_v1_Q.float32.bin"
    TPL_MAP_U  = "{dbpath}/sroll_templates/MAP_JMD_2048/template_sims_dustapr17_aligned_{freq}_000_full_1deg_ns2048_v1_U.float32.bin"
    TPL_FIT_Q  = "{dbpath}/sroll_templates/MAP_JMD_128/template_sims_dustapr17_aligned_{freq}_000_full_1deg_ns128_v1_Q.float32.bin"
    TPL_FIT_U  = "{dbpath}/sroll_templates/MAP_JMD_128/template_sims_dustapr17_aligned_{freq}_000_full_1deg_ns128_v1_U.float32.bin"
    TPL_DUST_I = "{dbpath}/sroll_templates/MAP_JMD_128/template_sims_dustapr17_FGproj_353GHz_full_1deg_ns128_v1_Dust_I.float32.bin"
    TPL_DUST_Q = "{dbpath}/sroll_templates/MAP_JMD_128/template_sims_dustapr17_FGproj_353GHz_full_1deg_ns128_v1_Dust_Q.float32.bin"
    TPL_DUST_U = "{dbpath}/sroll_templates/MAP_JMD_128/template_sims_dustapr17_FGproj_353GHz_full_1deg_ns128_v1_Dust_U.float32.bin"

  if lpar.fslcoef == None:
    if DATAVER == "REP6":
      lpar.fslcoef = 1
    else:
      lpar.fslcoef = 0

  # disable FITPOLEFF and FITANGLE for mono bolometer
  if len( lpar.bolometers) == 1:
    lpar.fitpoleff = 0
    lpar.fitangle = 0

#  if (lpar.freq == "353") and (len( lpar.bolometers) == 1):
#    raise Exception( "error: sroll doesn't support single bolometer runs at 353GHz (or comment this line...)")

  # subsets is the list of all detsets/single bolometers for which a map will be produced
  if lpar.map_detsets == None:
    if lpar.ffp10_maps and detset.endswith("ghz"):
      lpar.map_detsets  = [detset,
                           freq + "ds1",
                           freq + "ds2"]
      lpar.map_ringsets = [FULL + HM12 + S12345 + FULLODDEVEN + HM12ODDEVEN,
                           FULL,
                           FULL]
      if lpar.freq == "353":
        lpar.map_detsets  += [freq + "psb",
                              freq + "swb"]
        lpar.map_ringsets += [FULL + HM12 + S12345 + FULLODDEVEN + HM12ODDEVEN,
                              FULL + HM12]
      elif lpar.freq > "100":
        lpar.map_detsets  += [freq + "psb",
                              freq + "swb"]
        lpar.map_ringsets += [FULL + HM12,
                              FULL + HM12]
      if lpar.map_per_bolo:
        lpar.map_detsets  += DETSETS[detset]
        lpar.map_ringsets += [(FULL + HM12 + FULLODDEVEN) for b in DETSETS[detset]]
    else:
      lpar.map_detsets = ["all"] # produce a map with all bolometers

  if lpar.map_ringsets == None:
    lpar.map_ringsets = [FULL for ds in lpar.map_detsets]

  # halfring is the suffix used in output object names, can be anything containing either 1 or 2 eg. "", "_HR1" or "_HR2".
  # hrtag is the suffix used in input object names, is either "", "_1st" or "_2nd"
  if "1" in lpar.halfring:
    hrtag = "_1st"
  elif "2" in lpar.halfring:
    hrtag = "_2nd"
  else:
    hrtag = ""

  if lpar.all_surveys:
    domaxvraie = 0
  else:
    domaxvraie = 1

  if lpar.freq >= "545":
    remhdip = "REMHDIP = 1"
  else:
    remhdip = ""

  STRNOW = time.strftime("%c")

  sroll_parameters = """
# SRoll parameter file generated by {scriptfilename}
# for {user}
# on {srollhost}
# at {STRNOW}

# PARAMETERS

BeginRing  = {begring}
EndRing    = {endring}
RSTEP      = {rstep}
AVGCO      = 0
AVGDUST    = 0
CALCODUST  = 0
CUTRG      = 1
DODIPCAL   = 1
DOGAINDIP  = 1
DOMAXVRAIE = {domaxvraie}
D_NOPOL    = 0
GAINSTEP   = {gainstep}
NADU       = {nadu}
Nside      = 2048
REMDIP     = 1
{remhdip}
ADDDIP     = {adddip}
KCMBIN     = {kcmbin}
SAVEINTMAP = 0
TESTPOL    = 1
XI2STOP    = 1.0
seuilcond  = {seuilcond}
verbose    = {verbose}
number_of_SEED = 1
SEED1      = 1234
"""

  if lpar.srollver in ("s2", "s2s1"):
    sroll_parameters += """
# sroll2 specific parameters
NADUSTEP   = {nadustep}
NITT       = {n_sroll_iter}
FITANGLE   = {fitangle}
FITPOLEFF  = {fitpoleff}
saveCOV    = {saveCOV}
delta_psi  = {delta_psi}
"""

  if lpar.srollver == "s1":
    # sroll1 parameters renamed/removed in sroll2
    sroll_parameters += """
# sroll1 specific parameters
FITTHETA   = 0
"""

  sroll_parameters += "\nnumber_of_Calibration = %d\n" % len( lpar.bolometers)
  for pidx, pixname in enumerate( lpar.bolometers):
    if lpar.RD12calib: # running RD12 reproduction
      sroll_parameters += "Calibration%d = %.12g\n" % (pidx+1, RD12CALIB[pixname])
    else: # running FFP10 simulation
      sroll_parameters += "Calibration%d = %.12g\n" % (pidx+1, DX11CALIB[pixname])

  sroll_parameters += "\nnumber_of_CrossPol = %d\n" % len( lpar.bolometers)
  for pidx, pixname in enumerate( lpar.bolometers):
    sroll_parameters += "CrossPol%d = %g\n" % (pidx+1, XPOL[pixname])

  sroll_parameters += "\nnumber_of_NEP = %d\n" % len( lpar.bolometers)
  for pidx, pixname in enumerate( lpar.bolometers):
    sroll_parameters += "NEP%d = %.12g\n" % (pidx+1, RD12NEP[pixname])

  sroll_parameters += "\nnumber_of_FSLCOEF = %d\n" % len( lpar.bolometers)
  for pidx, pixname in enumerate( lpar.bolometers):
    sroll_parameters += "FSLCOEF%d = %f\n" % (pidx+1, lpar.fslcoef)

  sroll_parameters += "\nnumber_of_Monop = %d\n" % len( lpar.bolometers)
  for pidx, pixname in enumerate( lpar.bolometers):
    sroll_parameters += "Monop%d = %.12g\n" % (pidx+1, RD12MONOP[pixname])

  sroll_parameters += "\nnumber_of_OUT_NOPOL = %d\n" % len( lpar.map_detsets)
  for pidx, subset in enumerate( lpar.map_detsets):
    if subset == "all":
      # use all bolometers to build the map
      tempset = lpar.bolometers
    elif subset in DETSETS.keys():
      # use a single detset to build the map
      tempset = DETSETS[subset]
    else:
      # single bolometre map
      assert subset in BOLOID.keys()
      tempset = [subset]
    psb = 0
    for b in tempset:
      if b[-1] in "ab":
        psb += 1
    nopol = (psb < 2) # produce polarised map if detset contains at least 2 PSBs
    sroll_parameters += "OUT_NOPOL%d = %d\n" % (pidx+1, nopol)

  nbolomask = len( lpar.bolometers) * len( lpar.map_detsets)
  sroll_parameters += "\nnumber_of_bolomask = %d\n" % nbolomask
  for sidx, subset in enumerate( lpar.map_detsets):
    sroll_parameters += "# bolomask for %s map\n" % subset
    for pidx, pixname in enumerate( lpar.bolometers):
      if subset in DETSETS.keys():
        bolomask = int(pixname in DETSETS[subset])
      elif pixname == subset: # single bolometer map
        bolomask = 1
      elif subset.lower() == "all": # all bolometers map
        bolomask = 1
      else:
        bolomask = 0
      sroll_parameters += "bolomask%d = %d\n" % (sidx * len( lpar.bolometers) + pidx + 1, bolomask)

  sroll_parameters += "\nnumber_of_MAPRINGS = %d\n" % len( lpar.map_ringsets)
  for idx, mapring in enumerate( lpar.map_ringsets):
    sroll_parameters += "MAPRINGS%d = %d\n" % (idx+1, mapring)


  if lpar.freq != None:
    sroll_parameters += """
# INPUTS

Mask = {dbpath}/MASK_2048_GALACTIC/mask_RD_{freq}

"""
  else:
    sroll_parameters += """
# INPUTS

Mask = {dbpath}/MASK_2048_GALACTIC/mask_RD_353

"""

  # SRoll2 parameters
  if lpar.srollver in ("s2", "s2s1"):
    sroll_parameters += "in_template_map_I = " + TPL_MAP_I + "\n"
    sroll_parameters += "in_template_map_Q = " + TPL_MAP_Q + "\n"
    sroll_parameters += "in_template_map_U = " + TPL_MAP_U + "\n"

    if lpar.freq != "143":
      sroll_parameters += "Theo_CO           = " + TPL_CO + "\n"
      if TPL_13CO != "":
        sroll_parameters += "Theo_13CO         = " + TPL_13CO + "\n"

    sroll_parameters += "in_polar_fit_Q    = " + TPL_FIT_Q + "\n"
    sroll_parameters += "in_polar_fit_U    = " + TPL_FIT_U + "\n"

    sroll_parameters += "Theo_Dust_I       = " + TPL_DUST_I + "\n"
    sroll_parameters += "Theo_Dust_Q       = " + TPL_DUST_Q + "\n"
    sroll_parameters += "Theo_Dust_U       = " + TPL_DUST_U + "\n"

  # SRoll1 parameters
  elif lpar.srollver == "s1":
    if len( lpar.bolometers) > 1:
      # multi-bolometer
      if lpar.freq != "143":
        sroll_parameters += "Theo_CO     = " + TPL_CO + "\n"
        sroll_parameters += "Theo_CO_Q   = " + TPL_CO_Q + "\n"
        sroll_parameters += "Theo_CO_U   = " + TPL_CO_U + "\n"

      sroll_parameters += "Theo_Dust_I = " + TPL_DUST_I + "\n"
      sroll_parameters += "Theo_Dust_Q = " + TPL_DUST_Q + "\n"
      sroll_parameters += "Theo_Dust_U = " + TPL_DUST_U + "\n"

    elif pixname[-1] in "ab":
      # single PSB
      sroll_parameters += "Theo_CO     = " + TPL_DUST_Q + "\n"
      sroll_parameters += "Theo_CO_Q   = " + TPL_DUST_Q + "\n"
      sroll_parameters += "Theo_CO_U   = " + TPL_DUST_Q + "\n"

      sroll_parameters += "Theo_Dust_I = " + TPL_DUST_U + "\n"
      sroll_parameters += "Theo_Dust_Q = " + TPL_DUST_U + "\n"
      sroll_parameters += "Theo_Dust_U = " + TPL_DUST_U + "\n"
    else:
      # no Theo_CO and Theo_Dust for single SWB
      pass


  if ((lpar.freq == "100") and (len( lpar.bolometers) > 1)) or (lpar.freq == None):
    sroll_parameters += "Theo_FREEFREE = {dbpath}/MAP_JMD_128/freefree_NORM_I\n"

  if lpar.freq >= "545":
    sroll_parameters += "TEMPLATEMAP = {dbpath}/MAP_DX11d_noZodi_2048_GALACTIC_0240_27005/{freq}GHz_W_TauDeconv_planet_I\n"

  if lpar.input_hpr != None:
    sroll_parameters += "\nnumber_of_Signal_noPS = %d\n" % len( lpar.bolometers)
    for pidx, pixname in enumerate( lpar.bolometers):
      signops = lpar.input_hpr.format( **lpar.todict( locals()))
      sroll_parameters += "Signal_noPS%d = %s\n" % (pidx+1, signops)
  else:
    sroll_parameters += "\nnumber_of_stim_paramfiles = %d\n" % len( lpar.bolometers)
    for pidx, pixname in enumerate( lpar.bolometers):
      stim_parfile = lpar.parfile_name + "_%s.stim.par" % pixname
      sroll_parameters += "stim_paramfiles%d = %s\n" % (pidx+1, stim_parfile)

  sroll_parameters += "\nnumber_of_ADU = %d\n" % len( lpar.bolometers)
  for pidx, pixname in enumerate( lpar.bolometers):
    sroll_parameters += "ADU%d = {dbpath}/PBR_JMD/%s_REP6{hrtag}_adutot\n" % (pidx+1, pixname)

  sroll_parameters += "\nnumber_of_Badring = %d\n" % len( lpar.bolometers)
  for pidx, pixname in enumerate( lpar.bolometers):
    sroll_parameters += "Badring%d = {dbpath}/calROIs/%s_discarded_rings_dx11\n" % (pidx+1, BOLOID[pixname])

  sroll_parameters += "\nnumber_of_DipOrb_noPS = %d\n" % len( lpar.bolometers)
  for pidx, pixname in enumerate( lpar.bolometers):
    if DATAVER == "REP6":
      sroll_parameters += "DipOrb_noPS%d = {dbpath}/PBR_JMD/%s_REP6_diporb_quat\n" % (pidx+1, pixname)
    elif lpar.dipole == "RD12":
      sroll_parameters += "DipOrb_noPS%d = {dbpath}/HPRREPSM_PBR/%s_v64_MAY16{hrtag}_diporb\n" % (pidx+1, pixname)
    elif lpar.dipole == "HFI17":
      sroll_parameters += "DipOrb_noPS%d = {dbpath}/PBR_JMD/%s_dipHFI17_quad_hprbin\n" % (pidx+1, pixname)
    else:
      raise Exception( "you must set the <dipole> parameter.")

  sroll_parameters += "\nnumber_of_Hit_noPS = %d\n" % len( lpar.bolometers)
  for pidx, pixname in enumerate( lpar.bolometers):
    sroll_parameters += "Hit_noPS%d = {dbpath}/PBR_JMD/%s_REP6{hrtag}_hit\n" % (pidx+1, pixname)

  sroll_parameters += "\nnumber_of_Ptg_noPS = %d\n" % len( lpar.bolometers)
  for pidx, pixname in enumerate( lpar.bolometers):
    sroll_parameters += "Ptg_noPS%d = {dbpath}/PBR_JMD/%s_REP6{hrtag}_ptg\n" % (pidx+1, pixname)

  if lpar.theonops:
    theolist = ["PBR_JMD/{pixname}_REP6{hrtag}_H%d" % i for i in range( 4)]
#    if lpar.freq == "353":
#      theolist += ["PBR_JMD/{pixname}_REP6{hrtag}_R%d" % i for i in range( 1, 4)]
#    if lpar.freq in ["545", "857"]:
#      theolist += ["PBR_JMD/{pixname}_REP6{hrtag}_R%d" % i for i in range( 4)]
#      theolist += ["PBR_JMD/{pixname}_REP7{hrtag}_2"]

    sroll_parameters += "\nnumber_of_Theo_noPS = %d\n" % (len( lpar.bolometers) * len( theolist))
    for i in range( len( theolist)):
      for pidx, pixname in enumerate( lpar.bolometers):
        sroll_parameters += ("Theo_noPS%d = {dbpath}/%s\n" % (i*len( lpar.bolometers) + pidx+1, theolist[i])).format( **lpar.todict( locals()))

  if lpar.fslcoef:
    sroll_parameters += "\nnumber_of_fsl = %d\n" % len( lpar.bolometers)
    for pidx, pixname in enumerate( lpar.bolometers):
      sroll_parameters += "fsl%d = {dbpath}/PBR_JMD/%s_REP6_corrfslfgmodQ\n" % (pidx+1, pixname)

  sroll_parameters += """

# OUTPUT
"""

  sroll_parameters += "\nnumber_of_MAP = %d\n" % len( lpar.map_detsets)
  for pidx, subset in enumerate( lpar.map_detsets):
    sroll_parameters += "MAP%d = {output_prefix}_MAP/{run_name}{halfring}_%s\n" % (pidx+1, subset)

  sroll_parameters += "\nnumber_of_Out_Offset = %d\n" % len( lpar.bolometers)
  for pidx, pixname in enumerate( lpar.bolometers):
    sroll_parameters += "Out_Offset%d = {output_prefix}_VEC/{run_name}{halfring}_%s_offset\n" % (pidx+1, pixname)

  sroll_parameters += "\nnumber_of_Out_Offset_corr = %d\n" % len( lpar.bolometers)
  for pidx, pixname in enumerate( lpar.bolometers):
    sroll_parameters += "Out_Offset_corr%d = {output_prefix}_VEC/{run_name}{halfring}_%s_offset_corr\n" % (pidx+1, pixname)

  sroll_parameters += "\nnumber_of_Out_xi2 = %d\n" % len( lpar.bolometers)
  for pidx, pixname in enumerate( lpar.bolometers):
    sroll_parameters += "Out_xi2%d = {output_prefix}_VEC/{run_name}{halfring}_%s_xi2\n" % (pidx+1, pixname)

  sroll_parameters += "\nnumber_of_Out_xi2_corr = %d\n" % len( lpar.bolometers)
  for pidx, pixname in enumerate( lpar.bolometers):
    sroll_parameters += "Out_xi2_corr%d = {output_prefix}_VEC/{run_name}{halfring}_%s_corr\n" % (pidx+1, pixname)

  sroll_parameters = sroll_parameters.format( **lpar.todict( locals()))
  f = open( sroll_parfile, "wt")
  f.write( sroll_parameters)
  f.close()

  print "sroll parameter file: " + sroll_parfile
  if not check_param_inputs( sroll_parfile):
    param.cancel_run = True


################################################################################

def make_troll_parfile( param):

  lpar = copy.copy( param) # local working copy of parameters structure
  troll_parfile = lpar.parfile_name + ".troll.par"

  # try to guess version of input troll data or simulations, for default parameter values
  DATAVER = ""
  if (lpar.input_hpr != None):
    # try to guess DATAVER from input HPR name
    hprname = lpar.input_hpr.split("/")[-1].upper()
    if "REP6" in hprname:
      DATAVER = "REP6" # activate default parameters for 2015 data TOI projection
    elif "JAN18" in hprname:
      DATAVER = "JAN18"
    elif ("DEC16" in hprname) or ("MAR17" in hprname) or ("FFP10" in hprname):
      DATAVER = "DEC16" # FFP10 parameters
    else:
      raise Exception( "Unkown input signal HPR: " + lpar.input_hpr)
  else:
    # try to guess DATAVER from input TOI name
    assert lpar.input_toi != None
    toiname = lpar.input_toi.split("/")[-1].upper()
    if ("JAN18" in toiname):
      DATAVER = "JAN18"
    elif ("DEC16" in toiname) or ("MAR17" in toiname) or ("FFP10" in hprname) or ("DUSTOCT" in hprname) or ("DUSTAPR" in hprname):
      DATAVER = "DEC16" # FFP10 parameters
    else:
      raise Exception( "Unkown input signal TOI: " + lpar.input_toi)

  if lpar.RD12calib == None:
    if DATAVER == "REP6":
      lpar.RD12calib = True
    elif DATAVER == "DEC16":
      lpar.RD12calib = False
    elif DATAVER == "JAN18":
      lpar.RD12calib = True
    else:
      raise Exception( "Unable to guess RD12calib switch value from troll input_hpr (%s)" % lpar.input_hpr)


  STRNOW = time.strftime("%c")

  troll_parameters = """
# TRoll parameter file generated by {scriptfilename}
# for {user}
# on {srollhost}
# at {STRNOW}

# PARAMETERS

BeginRing  = {begring}
EndRing    = {endring}
# sum of the bolometer coefficients 1.0*8 = 8 @100GHz
AVGFREEFREE = 8.0
# sum of the bolometer coefficients 0.02*8 = 0.8 @100GHz
AVGDUST100 = 0.16
AVGCO143 = 0.000
AVGCO = 0
AVGDUST = 0
CALCODUST = 0
CUTRG = 1
DODIPCAL = 1
DOGAINDIP = 1
DOMAXVRAIE = 0
D_NOPOL = 0
FITANGLE = 0
FITPOLEFF = 0
GAINSTEP = 1
NADU = 1
NADUSTEP = 1
Nside = 2048
NITT = 7
REMDIP = 1
RSTEP = {rstep}
SAVEINTMAP = 0
TESTPOL = 1
XI2STOP = 1.0
seuilcond = 3e-4
verbose = {verbose}

number_of_SEED = 1
SEED1 = 1234

"""

  bolofreq = lpar.bolometers[0][0:3]
  troll_parameters += "\nnumber_of_Calibration = %d\n" % len( lpar.bolometers)
  for pidx, pixname in enumerate( lpar.bolometers):
    if pixname[0:3] != bolofreq: bolofreq = pixname[0:3]; troll_parameters += "\n"
    if lpar.RD12calib: # running RD12 reproduction
      troll_parameters += "Calibration%d = %.12g\n" % (pidx+1, RD12CALIB[pixname])
    else: # running FFP10 simulation
      troll_parameters += "Calibration%d = %.12g\n" % (pidx+1, DX11CALIB[pixname])

  bolofreq = lpar.bolometers[0][0:3]
  troll_parameters += "\nnumber_of_CrossPol = %d\n" % len( lpar.bolometers)
  for pidx, pixname in enumerate( lpar.bolometers):
    if pixname[0:3] != bolofreq: bolofreq = pixname[0:3]; troll_parameters += "\n"
    troll_parameters += "CrossPol%d = %g\n" % (pidx+1, XPOL[pixname])

  bolofreq = lpar.bolometers[0][0:3]
  troll_parameters += "\nnumber_of_NEP = %d\n" % len( lpar.bolometers)
  for pidx, pixname in enumerate( lpar.bolometers):
    if pixname[0:3] != bolofreq: bolofreq = pixname[0:3]; troll_parameters += "\n"
    troll_parameters += "NEP%d = %.12g\n" % (pidx+1, RD12NEP[pixname])

  bolofreq = lpar.bolometers[0][0:3]
  troll_parameters += "\nnumber_of_FSLCOEF = %d\n" % len( lpar.bolometers)
  for pidx, pixname in enumerate( lpar.bolometers):
    if pixname[0:3] != bolofreq: bolofreq = pixname[0:3]; troll_parameters += "\n"
    troll_parameters += "FSLCOEF%d = 0\n" % (pidx+1)

  bolofreq = lpar.bolometers[0][0:3]
  troll_parameters += "\nnumber_of_Monop = %d\n" % len( lpar.bolometers)
  for pidx, pixname in enumerate( lpar.bolometers):
    if pixname[0:3] != bolofreq: bolofreq = pixname[0:3]; troll_parameters += "\n"
    troll_parameters += "Monop%d = %.12g\n" % (pidx+1, RD12MONOP[pixname])

  troll_parameters += "\nnumber_of_OUT_NOPOL = 1\n"
  troll_parameters += "OUT_NOPOL1 = 0\n"

  bolofreq = lpar.bolometers[0][0:3]
  troll_parameters += "\nnumber_of_bolomask = %d\n" % len( lpar.bolometers)
  for pidx, pixname in enumerate( lpar.bolometers):
    if pixname[0:3] != bolofreq: bolofreq = pixname[0:3]; troll_parameters += "\n"
    troll_parameters += "bolomask%d = 1\n" % (pidx+1)


  troll_parameters += """
# INPUTS

Mask =    {dbpath}/MASK_2048_GALACTIC/mask_RD_100

# component templates

Theo_Dust_I =      {dbpath}/sroll_templates/MAP_JMD_128/SIMDUST_128_I.float32.bin
Theo_Dust_Q =      {dbpath}/sroll_templates/MAP_JMD_128/SIMDUST_128_Q.float32.bin
Theo_Dust_U =      {dbpath}/sroll_templates/MAP_JMD_128/SIMDUST_128_U.float32.bin

Theo_TDust_I =     {dbpath}/sroll_templates/MAP_JMD_128/SIMTDUST_128_I.float32.bin
Theo_TDust_Q =     {dbpath}/sroll_templates/MAP_JMD_128/SIMTDUST_128_Q.float32.bin
Theo_TDust_U =     {dbpath}/sroll_templates/MAP_JMD_128/SIMTDUST_128_U.float32.bin

Theo_FREEFREE =    {dbpath}/sroll_templates/MAP_JMD_128/SIMFREEFREE_128_I.float32.bin

Theo_CO =          {dbpath}/sroll_templates/MAP_JMD_128/COSIMU.float32.bin
Theo_13CO =        {dbpath}/sroll_templates/MAP_JMD_128/NEW_cleaned_13CO.float32.bin

in_synchro_map_I = {dbpath}/sroll_templates/MAP_JMD_128/SIMSYNCHRO_128_I.float32.bin
in_synchro_map_Q = {dbpath}/sroll_templates/MAP_JMD_128/SIMSYNCHRO_128_Q.float32.bin
in_synchro_map_U = {dbpath}/sroll_templates/MAP_JMD_128/SIMSYNCHRO_128_U.float32.bin

# templates per bolometer
"""

  bolofreq = lpar.bolometers[0][0:3]
  troll_parameters += "\nnumber_of_in_polar_fit_Q = %d\n" % len( lpar.bolometers)
  for pidx, pixname in enumerate( lpar.bolometers):
    if pixname[0:3] != bolofreq: bolofreq = pixname[0:3]; troll_parameters += "\n"
    troll_parameters += "in_polar_fit_Q%d = {dbpath}/sroll_templates/MAP_JMD_128/%s_SIMU_Q.float32.bin\n" % (pidx+1, bolofreq)

  bolofreq = lpar.bolometers[0][0:3]
  troll_parameters += "\nnumber_of_in_polar_fit_U = %d\n" % len( lpar.bolometers)
  for pidx, pixname in enumerate( lpar.bolometers):
    if pixname[0:3] != bolofreq: bolofreq = pixname[0:3]; troll_parameters += "\n"
    troll_parameters += "in_polar_fit_U%d = {dbpath}/sroll_templates/MAP_JMD_128/%s_SIMU_U.float32.bin\n" % (pidx+1, bolofreq)

  bolofreq = lpar.bolometers[0][0:3]
  troll_parameters += "\nnumber_of_in_template_map_I = %d\n" % len( lpar.bolometers)
  for pidx, pixname in enumerate( lpar.bolometers):
    if pixname[0:3] != bolofreq: bolofreq = pixname[0:3]; troll_parameters += "\n"
    troll_parameters += "in_template_map_I%d = {dbpath}/sroll_templates/MAP_JMD_2048/%s_SIMU_I.float32.bin\n" % (pidx+1, bolofreq)

  bolofreq = lpar.bolometers[0][0:3]
  troll_parameters += "\nnumber_of_in_template_map_Q = %d\n" % len( lpar.bolometers)
  for pidx, pixname in enumerate( lpar.bolometers):
    if pixname[0:3] != bolofreq: bolofreq = pixname[0:3]; troll_parameters += "\n"
    troll_parameters += "in_template_map_Q%d = {dbpath}/sroll_templates/MAP_JMD_2048/%s_SIMU_Q.float32.bin\n" % (pidx+1, bolofreq)

  bolofreq = lpar.bolometers[0][0:3]
  troll_parameters += "\nnumber_of_in_template_map_U = %d\n" % len( lpar.bolometers)
  for pidx, pixname in enumerate( lpar.bolometers):
    if pixname[0:3] != bolofreq: bolofreq = pixname[0:3]; troll_parameters += "\n"
    troll_parameters += "in_template_map_U%d = {dbpath}/sroll_templates/MAP_JMD_2048/%s_SIMU_U.float32.bin\n" % (pidx+1, bolofreq)

  assert lpar.input_hpr is not None
  bolofreq = lpar.bolometers[0][0:3]
  troll_parameters += "\nnumber_of_Signal_noPS = %d\n" % len( lpar.bolometers)
  for pidx, pixname in enumerate( lpar.bolometers):
    if pixname[0:3] != bolofreq: bolofreq = pixname[0:3]; troll_parameters += "\n"
    signops = lpar.input_hpr.format( **lpar.todict( locals()))
    troll_parameters += "Signal_noPS%d = %s\n" % (pidx+1, signops)

  bolofreq = lpar.bolometers[0][0:3]
  troll_parameters += "\nnumber_of_ADU = %d\n" % len( lpar.bolometers)
  for pidx, pixname in enumerate( lpar.bolometers):
    if pixname[0:3] != bolofreq: bolofreq = pixname[0:3]; troll_parameters += "\n"
    troll_parameters += "ADU%d = {dbpath}/PBR_JMD/%s_REP6_adutot\n" % (pidx+1, pixname)

  bolofreq = lpar.bolometers[0][0:3]
  troll_parameters += "\nnumber_of_Badring = %d\n" % len( lpar.bolometers)
  for pidx, pixname in enumerate( lpar.bolometers):
    if pixname[0:3] != bolofreq: bolofreq = pixname[0:3]; troll_parameters += "\n"
    troll_parameters += "Badring%d = {dbpath}/calROIs/%s_discarded_rings_dx11\n" % (pidx+1, BOLOID[pixname])

  bolofreq = lpar.bolometers[0][0:3]
  troll_parameters += "\nnumber_of_DipOrb_noPS = %d\n" % len( lpar.bolometers)
  for pidx, pixname in enumerate( lpar.bolometers):
    if pixname[0:3] != bolofreq: bolofreq = pixname[0:3]; troll_parameters += "\n"
    if lpar.dipole == "RD12":
      troll_parameters += "DipOrb_noPS%d = {dbpath}/HPRREPSM_PBR/%s_v64_MAY16_diporb\n" % (pidx+1, pixname)
    elif lpar.dipole == "HFI17":
      troll_parameters += "DipOrb_noPS%d = {dbpath}/PBR_JMD/%s_dipHFI17_quad_hprbin\n" % (pidx+1, pixname)
    else:
      raise Exception( "you must set the <dipole> parameter.")

  bolofreq = lpar.bolometers[0][0:3]
  troll_parameters += "\nnumber_of_Hit_noPS = %d\n" % len( lpar.bolometers)
  for pidx, pixname in enumerate( lpar.bolometers):
    if pixname[0:3] != bolofreq: bolofreq = pixname[0:3]; troll_parameters += "\n"
    troll_parameters += "Hit_noPS%d = {dbpath}/PBR_JMD/%s_REP6_hit\n" % (pidx+1, pixname)

  bolofreq = lpar.bolometers[0][0:3]
  troll_parameters += "\nnumber_of_Ptg_noPS = %d\n" % len( lpar.bolometers)
  for pidx, pixname in enumerate( lpar.bolometers):
    if pixname[0:3] != bolofreq: bolofreq = pixname[0:3]; troll_parameters += "\n"
    troll_parameters += "Ptg_noPS%d = {dbpath}/PBR_JMD/%s_REP6_ptg\n" % (pidx+1, pixname)

  theolist = ["PBR_JMD/{pixname}_REP6_H%d" % i for i in range( 4)]
  bolofreq = lpar.bolometers[0][0:3]
  troll_parameters += "\nnumber_of_Theo_noPS = %d\n" % (len( lpar.bolometers) * len( theolist))
  for i in range( len( theolist)):
    for pidx, pixname in enumerate( lpar.bolometers):
      if pixname[0:3] != bolofreq: bolofreq = pixname[0:3]; troll_parameters += "\n"
      troll_parameters += ("Theo_noPS%d = {dbpath}/%s\n" % (i*len( lpar.bolometers) + pidx+1, theolist[i])).format( **lpar.todict( locals()))

  bolofreq = lpar.bolometers[0][0:3]
  troll_parameters += "\nnumber_of_fsl = %d\n" % len( lpar.bolometers)
  for pidx, pixname in enumerate( lpar.bolometers):
    if pixname[0:3] != bolofreq: bolofreq = pixname[0:3]; troll_parameters += "\n"
    troll_parameters += "fsl%d = {dbpath}/PBR_JMD/%s_REP6_corrfslfgmodQ\n" % (pidx+1, pixname)

  troll_parameters += """

# OUTPUT
"""

  troll_parameters += "\nnumber_of_MAP = 1\n"
  troll_parameters += "MAP1 = {output_prefix}_MAP/{run_name}\n"

  bolofreq = lpar.bolometers[0][0:3]
  troll_parameters += "\nnumber_of_Out_Offset = %d\n" % len( lpar.bolometers)
  for pidx, pixname in enumerate( lpar.bolometers):
    if pixname[0:3] != bolofreq: bolofreq = pixname[0:3]; troll_parameters += "\n"
    troll_parameters += "Out_Offset%d = {output_prefix}_VEC/{run_name}_%s_offset\n" % (pidx+1, pixname)

  bolofreq = lpar.bolometers[0][0:3]
  troll_parameters += "\nnumber_of_Out_Offset_corr = %d\n" % len( lpar.bolometers)
  for pidx, pixname in enumerate( lpar.bolometers):
    if pixname[0:3] != bolofreq: bolofreq = pixname[0:3]; troll_parameters += "\n"
    troll_parameters += "Out_Offset_corr%d = {output_prefix}_VEC/{run_name}_%s_offset_corr\n" % (pidx+1, pixname)

  bolofreq = lpar.bolometers[0][0:3]
  troll_parameters += "\nnumber_of_Out_xi2 = %d\n" % len( lpar.bolometers)
  for pidx, pixname in enumerate( lpar.bolometers):
    if pixname[0:3] != bolofreq: bolofreq = pixname[0:3]; troll_parameters += "\n"
    troll_parameters += "Out_xi2%d = {output_prefix}_VEC/{run_name}_%s_xi2\n" % (pidx+1, pixname)

  bolofreq = lpar.bolometers[0][0:3]
  troll_parameters += "\nnumber_of_Out_xi2_corr = %d\n" % len( lpar.bolometers)
  for pidx, pixname in enumerate( lpar.bolometers):
    if pixname[0:3] != bolofreq: bolofreq = pixname[0:3]; troll_parameters += "\n"
    troll_parameters += "Out_xi2_corr%d = {output_prefix}_VEC/{run_name}_%s_xi2_corr\n" % (pidx+1, pixname)

  troll_parameters = troll_parameters.format( **lpar.todict( locals()))
  f = open( troll_parfile, "wt")
  f.write( troll_parameters)
  f.close()

  print "troll parameter file: " + troll_parfile
  if not check_param_inputs( troll_parfile):
    param.cancel_run = True


################################################################################

def run_sroll( param, submit=True):

  lpar = copy.copy( param) # local working copy of parameters structure
  lpar.checkinit()

  # sroll parameter file
  make_sroll_parfile( lpar)
  # job script and submission
  make_submission_script( lpar, "sroll", submit)


################################################################################

def run_troll( param, submit=True):

  lpar = copy.copy( param) # local working copy of parameters structure
  lpar.checkinit()

  # sroll parameter file
  make_troll_parfile( lpar)
  # job script and submission
  make_submission_script( lpar, "troll", submit)


################################################################################

def run( param, submit=True):

  lpar = copy.copy( param) # local working copy of parameters structure
  lpar.checkinit()

  # stim parameter files (one per bolometer)
  for pixname in lpar.bolometers:
    make_stim_parfile( lpar, pixname)
  # sroll parameter file
  make_sroll_parfile( lpar)
  # job script and submission
  make_submission_script( lpar, "sroll", submit)


################################################################################

def run_stimexe( param, submit=True):

  lpar = copy.copy( param) # local working copy of parameters structure
  lpar.checkinit()
  parfile_name = lpar.parfile_name

  for pixname in lpar.bolometers:
    # create stimexe parameter file
    lpar.parfile_name = parfile_name
    make_stim_parfile( lpar, pixname)
    # job script and submission
    lpar.parfile_name = parfile_name + "_" + pixname
    make_submission_script( lpar, "stimexe", submit)


################################################################################

def job_count( prefix):
  if SROLLHOST == "OCCIGEN":
    return int( subprocess.check_output("squeue -o '%%.40j' -u %s | grep %s | wc -l" % (USER, prefix), shell=True))


################################################################################

def check_stimHPR( hprtpl, detset, firstiter=0, lastiter=100):
# eg. check_stimHPR( "/scratch/cnt0028/ias1717/smottet/stimr56_VEC/{pixname}_JAN18_stimHPR_{iternum:03d}.float32.bin", ["100-1a"])
  for pixname in DETSETS[detset]:
    iternum = firstiter
    hprname = hprtpl.format( **locals())
    print hprname
    refhpr = numpy.fromfile( hprname, dtype="float32")
    for iternum in range( firstiter+1, lastiter+1):
      hprname = hprtpl.format( **locals())
      print hprname
      if not os.path.exists( hprname):
        print "!!! file not found !!!"
        continue
      checkhpr = numpy.fromfile( hprname, dtype="float32")
      stools.compare_HPR_hits( refhpr, checkhpr, verbose=True)


################################################################################

def run_sroll2_data_2018( detset, submit=False):
  param = SrollParam( detset) # get default values for all parameters
  param.input_hpr  = "{dbpath}/PBR_JMD/{pixname}_REP6"
  param.dipole     = "RD12"
  param.bolometers = [b for b in DETSETS[detset]]
  param.run_name   = "sroll2r71_data_" + detset
  param.output_prefix = param.dbpath + "sroll2r71_data"
  param.map_detsets   = [detset, param.freq+"ds1", param.freq+"ds2"] + param.bolometers
  param.map_ringsets  = [FULL + HM12 + FULLODDEVEN + S12345] * len( param.map_detsets)
  param.saveCOV = 1 # save IIQQUU
  print( "param.bolometers=" + str(param.bolometers))
  run_sroll( param, submit=submit)


################################################################################

def make_JAN18_stimexe_HPR( detset, run_name=None, first_iter=0, itercount=1, nonoise=False, linadc=False, submit=False):
  for pixname in DETSETS[detset]:
    param = SrollParam( pixname) # get default values for all parameters
#    assert param.srollhost == "OCCIGEN"
    param.input_toi = "{dbpath}/e2e_common_TOI/{pixname}_JAN18sky_TOI"
# THIS IS A KNOWN BUG, JAN18 STIMHPR HAVE BEEN PRODUCED WITH RD12CALIB=FALSE (bug in RD12calib value autodetection)
# BUT JAN18R60 SROLL HAS RUN WITH RD12CALIB=TRUE
# THE RIGHT VALUE TO USE FOR NEXT SIMULATIONS (2019?) IS RD12CALIB=TRUE FOR BOTH
    param.RD12calib = False
    param.begring = BEGMISS
    param.endring = ENDMISS
    param.bolometers = [pixname]
    param.random_seed = first_iter
    param.iterations  = itercount
    if run_name == None:
      param.run_name = "JAN18"
    else:
      param.run_name = run_name
    if nonoise:
      param.add_noise = 0
      if run_name == None:
        param.run_name += "_nonoise"
    if linadc:
      param.sim_inl_name = 0
      param.corr_inl_name = 0
      if run_name == None:
        param.run_name += "_linadc"
#    param.nodes = 43  # occigen 1024 cores
    param.joblabel = param.run_name + "_%s_%03d" % (pixname, param.random_seed)
    param.output_prefix = "stimr56"
    param.save_hpr = True
    print( "param.bolometers=" + str( param.bolometers))
    run_stimexe( param, submit=submit)


################################################################################

def make_stimexe_projonly_HPR( detset, in_toi, out_prefix, do_SHPR=1, run_name=None, first_iter=0, submit=False):
  for pixname in DETSETS[detset]:
    param = SrollParam( pixname) # get default values for all parameters
    param.input_toi = in_toi
#    param.input_toi = "{dbpath}/e2e_common_TOI/{pixname}_JAN18sky_TOI"
#    param.input_toi = "{dbpath}/e2e_common_TOI/APR20_sky_{pixname}_mmTOI.fits"
    param.RD12calib   = True
    param.begring     = BEGMISS
    param.endring     = ENDMISS
    param.bolometers  = [pixname]
    param.random_seed = first_iter
    param.iterations  = 1
    param.do_SHPR     = do_SHPR
    if run_name == None:
      param.run_name = "APR20_projonly"
    else:
      param.run_name = run_name
    param.add_noise = 0
    param.sim_adu = 0
    if param.srollhost == "M4":
      param.nodes = 11 # M4 256 cores
    elif param.srollhost == "OCCIGEN":
      param.nodes = 11  # occigen 256 cores
    param.joblabel = param.run_name + "_%s_%03d" % (pixname, param.random_seed)
    param.output_prefix = out_prefix
    param.save_hpr = True
    print( "param.bolometers=" + str( param.bolometers))
    run_stimexe( param, submit=submit)


################################################################################

def run_sroll2_JAN18r60( detset, isim, submit=False):
  strsim = "%03d" % isim
  param = SrollParam( detset) # get default values for all parameters
  param.input_hpr  = "{dbpath}/JAN18_stimHPR/{pixname}_JAN18_stimHPR_" + strsim + ".float32.bin"
  param.saveCOV    = 0
  param.dipole     = "HFI17"
  param.bolometers = [b for b in DETSETS[detset]]
  param.run_name   = "JAN18r60_" + strsim + "_" + detset
  param.output_prefix = "JAN18r60"
  param.map_detsets   = [detset, param.freq+"ds1", param.freq+"ds2"] # + param.bolometers
  param.map_ringsets  = [FULL + HM12 + FULLODDEVEN] * len( param.map_detsets)
  print( "param.bolometers=" + str( param.bolometers))
  run_sroll( param, submit=submit)


################################################################################

def run_sroll1_JAN18( detset, isim, submit=False):
  strsim = "%03d" % isim
  param = SrollParam( detset, srollver="s1") # use SRoll1 instead of SRoll2
  param.input_hpr  = "{dbpath}/JAN18_stimHPR/{pixname}_JAN18_stimHPR_" + strsim + ".float32.bin"
  param.saveCOV    = 0
  param.dipole     = "HFI17"
  param.bolometers = [b for b in DETSETS[detset]]
  param.run_name   = "s1JAN18_" + strsim + "_" + detset
  param.output_prefix = "s1"
  param.map_detsets   = [detset, param.freq+"ds1", param.freq+"ds2"] # + param.bolometers
  param.map_ringsets  = [FULL + HM12 + FULLODDEVEN] * len( param.map_detsets)
  print( "param.bolometers=" + str( param.bolometers))
  run_sroll( param, submit=submit)


################################################################################

def run_test_troll_JAN18( detset, isim, submit=False):
  assert detset in ["psb", "ds1", "ds2"]
  nonoise = ""
  nonoise = "_nonoise_linadc"
  strsim = "%03d" % isim
  param = SrollParam( "353ghz", srollver="s3") # get default values for all parameters
  param.input_hpr  = "{dbpath}/JAN18_stimHPR/{pixname}_JAN18"+nonoise+"_stimHPR_" + strsim + ".float32.bin"
  param.saveCOV    = 1
  param.rstep      = 10
  param.verbose    = 1
  param.walltime   = 1.0
  param.dipole     = "HFI17"
  param.bolometers = DETSETS["100"+detset]+DETSETS["143"+detset]+DETSETS["217"+detset]+DETSETS["353"+detset]
  param.run_name   = "JAN18_all"+detset+nonoise+"_rstep%d" % (param.rstep)
  param.output_prefix = "/scratch/cnt0028/ias1717/SHARED/bware/sroll3test"
  print( "param.bolometers=" + str( param.bolometers))
  run_troll( param, submit=submit)


################################################################################
################################################################################

if (__name__ == '__main__'):

  print( "%s: Please edit the '__main__' section of this script to submit sroll/stim jobs with the command line.\n" % sys.argv[0])
  exit( 0)

  detset = sys.argv[1]
  assert detset in DETSETS

# use this to project TOI to HPR or SHPR using stim with no effects
#  make_stimexe_projonly_HPR( detset,
#                             "{dbpath}/e2e_common_TOI/APR20_sky_{pixname}_mmTOI.fits",
#                             "tst_SHPR",
#                             first_iter=0, submit=False)

# use this to reproduce sroll2 2018 maps on data:
#  run_sroll2_data_2018( detset)

# use this to produce JAN18 ("the 500 bware sims") stim HPR (sroll inputs)
#  make_JAN18_stimexe_HPR( detset, 0, itercount=1, nonoise=True, linadc=False, submit=True)

# use this to produce JAN18r60 maps ("the 500 bware sims")
#  run_sroll2_JAN18r60( detset, 0)

# use this to test SRoll1 on JAN18 simulation HPRs
#  run_sroll1_JAN18( detset, 0)

# use this to test Troll on JAN18 simulation HPRs without noise with linear ADC and all 32 polarised bolometers (100psb+143psb+217psb+353psb)
#  run_test_troll_JAN18( "psb", 0)
