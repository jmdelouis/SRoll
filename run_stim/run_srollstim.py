#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, sys
import subprocess
import time, calendar
import math
import numpy
import copy
from IMO_4_27 import * # BOLOID, DETSETS, CALIB, NEP, XPOL, ELECWNOISE, PHOTWNOISE


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

################################################################################
# SrollParam: method-less class to hold all parameters needed to write
# sroll/stim parameter files and submission scripts

class SrollParam( object):

  def __init__( self):
    
    # init environment variables and paths
    try:
      self.srollhost = os.environ["SROLLHOST"]
    except:
      raise Exception( "You must 'source ./srollex_setenv.sh' before using this script")
    self.srolldir  = os.environ["SROLLDIR"]
    self.directory = os.path.dirname( os.path.realpath( __file__)) + "/" # run_srollstim.py directory

    if self.srollhost == "M3":
      self.dbpath = os.environ["DMCDATA"] # /data/dmc/MISS03/DATA
      self.brimo_path = "/wrk/symottet/brimo_4_27"
      self.scratch = "/redtruck/SimuData"
      self.batch_suffix = ".qsub"
    
    elif self.srollhost == "M4":
      self.dbpath = os.environ["NODMCDATA"] # /pscratch1/RD12_data/dmc/MISS03/DATA
      self.brimo_path = "/wrk/symottet/brimo_4_27"
      self.scratch = self.dbpath
      self.batch_suffix = ".qsub"
    
    elif self.srollhost == "EDISON":
      self.dbpath = os.environ["NODMCDATA"] # /scratch3/scratchdirs/paganol/RD12_data/dmc/MISS03/DATA
      self.brimo_path = "/project/projectdirs/planck/data/hfi/RD12_data/brimo/brimo_4_27"
      self.scratch = os.environ["SCRATCH"] # /scratch*/scratchdirs/<userlogin>
      self.batch_suffix = ".slurm"
    
    elif self.srollhost == "CORI":
      self.dbpath = os.environ["NODMCDATA"] # /project/projectdirs/planck/data/hfi/RD12_data/dmc/MISS03/DATA
      self.brimo_path = "/project/projectdirs/planck/data/hfi/RD12_data/brimo/brimo_4_27"
      self.scratch = os.environ["SCRATCH"] # /scratch*/scratchdirs/<userlogin>
      self.batch_suffix = ".slurm"
    
    elif self.srollhost == "OCCIGEN":
      self.dbpath = os.environ["NODMCDATA"] # /scratch/cnt0028/ias1717/SHARED/RD12_data/dmc/MISS03/DATA
      self.brimo_path = "/scratch/cnt0028/ias1717/SHARED/RD12_data/brimo_4_27"
      self.scratch = os.environ["SCRATCHDIR"] # /scratch/cnt0028/ias1717/<userlogin>
      self.batch_suffix = ".slurm"
    
    else:
      raise Exception( "Unknown platform: " + self.srollhost)
  
    # default parameters values
    self.bolometers     = None # list of single detector names
    self.freq           = None
    self.begring        = 240
    self.endring        = 26050
    self.run_name       = "default_run_name"
    self.output_prefix  = "runsrollstim"
    self.halfring       = "" # user's halfring tag, eg. "1st" or "HR1", etc.
    self.RD12calib      = None # set to True to use RD12 calib factors, to False for DX11 calib and None to throw an exception

    # batch parameters
    self.joblabel       = None
    self.parfile_name   = None
    self.nodes          = None

    # stim parameters
    self.iterations     = 1
    self.random_seed    = 0
    self.rings_per_read = -1
    self.input_toi      = None
    self.signal_out     = 0
    self.hpr_out        = 0
    self.add_noise      = 1
    self.despike_flag   = 1
    self.conv_deconv    = 1
    self.sim_adu        = 1
    self.sim_inl_name   = None
    self.corr_inl_name  = None
    self.add_oof        = 1
    self.stim_calib     = 0
    self.save_toi       = False
    self.save_hpr       = False

    # sroll parameters
    self.srollver       = 10 # 10: sroll 1.0
    self.dipole         = None # set to "RD12" for RD12/FFP10, or "HFI17" for 2018+
    self.rstep          = 1
    self.gainstep       = None
    self.ffp10_maps     = False
    self.map_per_bolo   = False # requires self.ffp10_maps = True
    self.all_surveys    = False
    self.map_detsets    = None # list of detsets for which to produce a map
    self.map_ringsets   = None # list of ringsets (ringcuts) for each map_detsets
    self.input_hpr      = None
    self.theonops       = True # set to False to not fit the transfer function with Theo_noPS parmeters
    self.kcmbin         = 0
    self.adddip         = 0


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
    if lpar.srollver == 10:
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

  else:
    raise Exception( "Unknown exe_name: " + exe_name)

  batch_file  = lpar.parfile_name + lpar.batch_suffix
  logfilename = lpar.parfile_name + ".log"

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
exit 0
""".format( **lpar.todict( locals()))

    f = open( batch_file, "wt")
    f.write( batch_script)
    f.close()
    shell_cmd = "qsub " + batch_file

#-----------------------------------------------------------------------

  if lpar.srollhost == "M4":
    wallhours = 4 * (lpar.iterations+1)
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
exit 0
""".format( **lpar.todict( locals()))

    f = open( batch_file, "wt")
    f.write( batch_script)
    f.close()
    shell_cmd = "sbatch " + batch_file

#-----------------------------------------------------------------------

  if lpar.srollhost == "OCCIGEN":
    wallhours = 2.0 * (lpar.iterations)
    if len( lpar.bolometers) > 12:
      wallhours *= 2.0
    walldec, wallint = math.modf( wallhours) 
    walltime = "%d:%d:00" % ( wallint, walldec*60)

    batch_script = """#!/bin/bash -l
#SBATCH -J {joblabel}
#SBATCH --nodes={nodes}
#SBATCH --constraint=HSW24
#SBATCH --ntasks={mpisize}
#SBATCH --ntasks-per-node=24
#SBATCH --threads-per-core=1
#SBATCH --mem=50GB            # --mem=100GB to force use of big RAM nodes
#SBATCH --time={walltime}

cd {srolldir}
source ./srollex_setenv.sh

cd {exe_dir}
srun --mpi=pmi2 -K1 -n $SLURM_NTASKS --distribution=cyclic ./{exe_name} {exe_parfile} &> {logfilename}

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


################################################################################

def make_stim_parfile( param, pixname):

  lpar = copy.copy( param) # local working copy of parameters structure

  boloid   = BOLOID[pixname]
  beltchan = boloid[0:2]
  pixfreq  = pixname[0:3]

  stim_parfile = lpar.parfile_name + "_%s.stim.par" % pixname

  assert lpar.sim_adu in (0, 1)
  assert lpar.add_noise in (0, 1)

  assert lpar.input_toi != None
#  lpar.input_toi = "{dbpath}/e2e_common_TOI/{pixname}_ffp10_fgdip_cmbfid"
#  lpar.input_toi = "{dbpath}/e2e_common_TOI/{pixname}_dustapr17"
#  lpar.input_toi = "{dbpath}/e2e_common_TOI/{pixname}_DEC17sky"

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
    if (pixfreq >= "353") or not (pixname[-1] in "ab"):
      lpar.sim_inl_name = "{dbpath}/RAW_VECT_FIT_INFO/{pixname}_SIMINL1406_offsets"
    else:
      lpar.sim_inl_name = "{dbpath}/RAW_VECT_FIT_INFO/{pixname}_SIMINL2803_offsets"

  if lpar.save_toi == True:
    if lpar.signal_out == 0:
      lpar.signal_out = lpar.output_prefix + "_VEC/{pixname}_%s%s_stimTOI" % (lpar.run_name, lpar.halfring)
  else:
    lpar.signal_out = 0

  if lpar.save_hpr == True:
    if lpar.hpr_out == 0:
      lpar.hpr_out = lpar.output_prefix + "_VEC/{pixname}_%s%s_stimHPR" % (lpar.run_name, lpar.halfring)
  else:
    lpar.hpr_out = 0

  if pixfreq >= "545":
   lpar.sim_adu = 0
   lpar.despike_flag = 0
   lpar.add_oof = 0
   lpar.stim_calib = -1
  else:
   lpar.despike_flag = lpar.add_noise    
   lpar.add_oof = lpar.add_noise    

  if lpar.sim_adu == 1:
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

  stim_parameters = """
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
rawgain_name = {dbpath}/RAW_GAIN/{beltchan}_INL_GP41N
rawcst_name  = {dbpath}/RAW_CST/{beltchan}_LM4K_DX11
#fourk_name   = {dbpath}/RAW_4K/{beltchan}_HARMS_LM4K_DX11_PIOFLOAT
fourk_name   = 0
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


number_of_add_hpr_name = {sim_adu}
{add_hpr}add_hpr_name1 = {dbpath}/PBR_JMD/{pixname}_REP6{hrtag}_4k_residu

number_of_add_hpr_factor = {sim_adu}
{add_hpr}add_hpr_factor1 = 1.0

""".format( **lpar.todict( locals()))

  stim_parameters = stim_parameters.format( **lpar.todict( locals())) # interpret formats in formats
  f = open( stim_parfile, "wt")
  f.write( stim_parameters)
  f.close()
  print "stim parameter file: " + stim_parfile


################################################################################

def make_sroll_parfile( param):

  lpar = copy.copy( param) # local working copy of parameters structure
  sroll_parfile = lpar.parfile_name + ".sroll.par"
  if lpar.RD12calib == None:
    if (lpar.input_hpr != None):
      hprname = lpar.input_hpr.split("/")[-1]
      if "REP6" in hprname:
        RD12 = True
      elif ("DEC16" in hprname) or ("MAR17" in hprname) or ("ffp10" in hprname):
        RD12 = False
      elif "DEC17" in hprname:
        RD12 = True
      else:
        raise Exception( "Unable to guess RD12calib switch value from sroll input_hpr")
    else:
      raise Exception( "RD12calib *must* now be set when using stim and sroll")
  else:
    RD12 = lpar.RD12calib
    
  if lpar.freq == None:
    seuilcond = 3e-4
    par_gainstep = 32
  elif lpar.freq <= "217":
    seuilcond = 3e-4
    par_gainstep = 128
  elif lpar.freq == "353":
    seuilcond = 3e-3
    par_gainstep = 32
  elif lpar.freq >= "545":
    seuilcond = 1e-6
    par_gainstep = 32
  else:
    assert False

  if lpar.gainstep != None:
    par_gainstep = lpar.gainstep

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
    remhdip = "REMHDIP=1"
  else:
    remhdip = ""


  sroll_parameters = """
# SRoll parameter file generated by Srollex/run_stim/run_srollstim.py

# PARAMETERS

BeginRing = {begring}
EndRing   = {endring}
RSTEP = {rstep}
AVGCO = 0
AVGDUST = 0
CALCODUST = 0
CUTRG = 1
DODIPCAL = 1
DOGAINDIP = 1
DOMAXVRAIE = {domaxvraie}
D_NOPOL = 0
FITTHETA = 0
GAINSTEP = {par_gainstep}
NADU = 1
Nside = 2048
REMDIP = 1
{remhdip}
ADDDIP = {adddip}
KCMBIN = {kcmbin}
SAVEINTMAP = 0
TESTPOL = 1
XI2STOP = 1.0
seuilcond = {seuilcond}
verbose = 0
number_of_SEED = 1
SEED1 = 1234

"""

  sroll_parameters += "\nnumber_of_Calibration = %d\n" % len( lpar.bolometers)
  for pidx, pixname in enumerate( lpar.bolometers):
    if RD12: # running RD12 reproduction
      sroll_parameters += "Calibration%d = %g\n" % (pidx+1, RD12CALIB[pixname])
    else: # running FFP10 simulation
      sroll_parameters += "Calibration%d = %g\n" % (pidx+1, DX11CALIB[pixname])

  sroll_parameters += "\nnumber_of_CrossPol = %d\n" % len( lpar.bolometers)
  for pidx, pixname in enumerate( lpar.bolometers):
    sroll_parameters += "CrossPol%d = %g\n" % (pidx+1, XPOL[pixname])

  sroll_parameters += "\nnumber_of_NEP = %d\n" % len( lpar.bolometers)
  for pidx, pixname in enumerate( lpar.bolometers):
    sroll_parameters += "NEP%d = %g\n" % (pidx+1, RD12NEP[pixname])

  sroll_parameters += "\nnumber_of_FSLCOEF = %d\n" % len( lpar.bolometers)
  for pidx, pixname in enumerate( lpar.bolometers):
    sroll_parameters += "FSLCOEF%d = 0.0\n" % (pidx+1)

  sroll_parameters += "\nnumber_of_Monop = %d\n" % len( lpar.bolometers)
  for pidx, pixname in enumerate( lpar.bolometers):
    sroll_parameters += "Monop%d = %g\n" % (pidx+1, RD12MONOP[pixname])

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


  if len( lpar.bolometers) > 1:
    # polarised detset
    if lpar.freq != "143":
      sroll_parameters += """
Theo_CO     = {dbpath}/MAP_JMD_128/COMAP
Theo_CO_Q   = {dbpath}/MAP_JMD_128/CO_ORT_Q
Theo_CO_U   = {dbpath}/MAP_JMD_128/CO_ORT_U
"""
    sroll_parameters += """
Theo_Dust_I = {dbpath}/MAP_JMD_128/Dust_I
Theo_Dust_Q = {dbpath}/MAP_JMD_128/Dust_Q
Theo_Dust_U = {dbpath}/MAP_JMD_128/Dust_U
"""

  elif pixname[-1] in "ab":
    # single PSB
    sroll_parameters += """
Theo_CO     = {dbpath}/MAP_JMD_128/Dust_Q
Theo_CO_Q   = {dbpath}/MAP_JMD_128/Dust_Q
Theo_CO_U   = {dbpath}/MAP_JMD_128/Dust_Q
Theo_Dust_I = {dbpath}/MAP_JMD_128/Dust_U
Theo_Dust_Q = {dbpath}/MAP_JMD_128/Dust_U
Theo_Dust_U = {dbpath}/MAP_JMD_128/Dust_U
"""
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
    if lpar.dipole == "RD12":
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
    if lpar.freq == "353":
      theolist += ["PBR_JMD/{pixname}_REP6{hrtag}_R%d" % i for i in range( 1, 4)]
    if lpar.freq in ["545", "857"]:
      theolist += ["PBR_JMD/{pixname}_REP6{hrtag}_R%d" % i for i in range( 4)]
      theolist += ["PBR_JMD/{pixname}_REP7{hrtag}_2"]
      
    sroll_parameters += "\nnumber_of_Theo_noPS = %d\n" % (len( lpar.bolometers) * len( theolist))
    for i in range( len( theolist)):
      for pidx, pixname in enumerate( lpar.bolometers):
        sroll_parameters += ("Theo_noPS%d = {dbpath}/%s\n" % (i*len( lpar.bolometers) + pidx+1, theolist[i])).format( **lpar.todict( locals()))

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


################################################################################

def run_sroll( param, submit=True):

  lpar = copy.copy( param) # local working copy of parameters structure
  lpar.checkinit()

  # sroll parameter file
  make_sroll_parfile( lpar)
  # job script and submission
  make_submission_script( lpar, "sroll", submit)


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
################################################################################

if (__name__ == '__main__'):

  print( "Please edit the '__main__' section of this script to submit sroll jobs with the command line.\n")
  exit( 0)


  # full quick stim+sroll on bolometer test
  # M3 32 nodes: 25min.
  if 0:
    param = SrollParam() # get default values for all parameters
    pixname = "353-3b" # 353-3b has the biggest difference between DX11 and RD12 (2.5%)
#    param.input_toi = "{dbpath}/e2e_common_TOI/{pixname}_ffp10_fgdip_cmbfid"
    param.input_toi = "{dbpath}/e2e_common_TOI/{pixname}_DEC17sky"
    if "ffp10" in param.input_toi:
      param.run_name = "tstFFP10_" + pixname
      param.RD12calib = False # reproduce FFP10 sims
      param.dipole = "RD12" # reproduce FFP10 sims
    elif "DEC17" in param.input_toi:
      param.run_name = "tstDEC17_" + pixname
      param.RD12calib = True # test 2018 sims
      param.dipole = "HFI17" # test 2018 sims
    else:
      raise Exception( "please set <input_toi>")
    param.save_toi = True
    param.save_hpr = True
    param.srollver = 10
    param.rstep    = 1
    param.gainstep = 1
    param.nodes    = 32
    param.begring  = BEGMISS
    param.endring  = ENDMISS
    param.bolometers = [pixname]
    param.output_prefix = "sr55rc"
    param.map_detsets = [pixname,]
    param.map_ringsets = [FULL + HM12 + FULLODDEVEN] * len( param.map_detsets)
    print( "param.bolometers=" + str(param.bolometers))
    run( param, submit=False)



  # sroll input HPR production with stimexe
  if 0:
    param = SrollParam() # get default values for all parameters
    param.RD12claib = False # reproduce FFP10 sims
    param.begring = BEGMISS
    param.endring = ENDMISS
    param.bolometers = ["143-1a", "143-1b"]
    param.random_seed = 75
    param.iterations  = 2
    param.nodes = 20 # the whole timeline must fit in memory, needs more than 480GB, less than 640GB
    param.run_name = "MAR17v1_JAN18"
    param.output_prefix = "stimNOV17"
    param.save_hpr = True
    print( "param.bolometers=" + str(param.bolometers))
    run_stimexe( param, submit=False)



  if 0:
    # M4: 100-1a, 256 cores, rstep 10: 2min
    # M4: 100ds1, 512 cores, rstep 10: 20min
    # M4: 100ds1, 512 cores, rstep 1: 1h15, 1TB RAM
    # M4: 100ghz, 512 cores, rstep 1: 2h15, 1.4TB RAM
    # occigen: 100ds1+143ds1+217ds1+353ds1, 2048 cores, rstep 1: 1h45
    for freq in [100, 143, 217, 353, 545, 857]:
      detset = "%dghz" % freq
      param = SrollParam() # get default values for all parameters
      param.rstep = 1
      param.gainstep = 1
#      param.input_hpr = "{dbpath}/PBR_JMD/{pixname}_REP6"
      param.input_hpr = "/redtruck/SimuData/FFP10_HFI_input_maps/Projected/HPR/{pixname}_ffp10_fgdip_hprbin"
#      param.input_hpr = "/redtruck/SimuData/FFP10_HFI_input_maps/Projected/HPR/{pixname}_dustapr17_fgdip_hprbin"
      param.dipole = "RD12"
      param.kcmbin = 1
      param.bolometers = [b for b in DETSETS[detset]]
      if "dustapr17" in param.input_hpr:
        dusttag = "dustapr17_"
      else:
        dusttag = ""
      param.run_name = "psmprojsroll_" + dusttag + detset
      param.output_prefix = "stimNOV17"
      param.map_detsets = [detset,]
      if freq == 353:
        param.map_detsets += ["353psb",]
      param.map_ringsets = [FULL + HM12 + FULLODDEVEN] * len( param.map_detsets)
      print( "param.bolometers=" + str(param.bolometers))
      run_sroll( param, submit=False)
  
