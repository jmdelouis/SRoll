#!/usr/bin/env python

import os
import numpy
import healpy
import tempfile
import subprocess
import array
import time

try:
  import pyfits
except:
  import astropy.io.fits as pyfits

from IMO_4_27 import * # BOLOID, DETSETS, CALIB, NEP, XPOL, ELECWNOISE, PHOTWNOISE
from beginringindex import BRI

# guess task count and task number with OpenMPI environment variables
# http://www.open-mpi.de/faq/?category=running#mpi-environmental-variables
MPI_SIZE = int( os.getenv( 'OMPI_COMM_WORLD_SIZE', 1))
MPI_RANK = int( os.getenv( 'OMPI_COMM_WORLD_RANK', 0))
if MPI_SIZE == 1:
  # srun -mpi=pmi2 environment variables
  # https://slurm.schedmd.com/srun.html
  MPI_SIZE = int( os.getenv( 'SLURM_NTASKS', 1))
  MPI_RANK = int( os.getenv( 'SLURM_PROCID', 0))

HPRSIZE = 27664
PBRSIZE = 10822

BEGMISS = 240
ENDMISS = 26050
BEGHM1  = BEGMISS
ENDHM1  = 13143
BEGHM2  = ENDHM1+1
ENDHM2  = ENDMISS

LOGNAME   = os.environ["LOGNAME"] # system user login, eg. symottet
DMCNAME   = os.getenv('DMCNAME', None) # only defined on magique3
DMCDATA   = os.getenv('DMCDATA', None) # defined on magique3 and magique4
NODMCDATA = os.getenv('NODMCDATA', None)

SROLLHOST = os.environ["SROLLHOST"]
SROLLDIR  = os.environ["SROLLDIR"]

# subprocess.check_outptut() doesn't exist on M3...
def subprocess_check_output( cmd):
  process = subprocess.Popen( cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  (stdout, stderr) = process.communicate() # wait for subprocess to end and get its ouputs
  if (process.returncode != 0):
    raise Exception( "error in %s:\n%s\n" % (cmd, stderr))
  return stdout

HOSTNAME = subprocess_check_output( "hostname").strip()

if SROLLHOST == "M4":
  DBPATH   = os.environ["NODMCDATA"] # /pscratch1/RD12_data/dmc/MISS03/DATA
  DX12DIR  = "/m3gpfs3/delouis/PROD_RD12ll/DATA/PROD_RD12_RC4"
  MONTAGE  = "/softs/ImageMagick/7/bin/montage"
  SPICEEXE = "export OMP_NUM_THREADS=24; /home/hivon/softs/spice/extra/PolSpice_v03-02-00_m4/src/spice"

elif SROLLHOST == "M3":
  DBPATH   = os.environ["DMCDATA"]   # /data/dmc/MISS03/DATA
  DX12DIR  = "/redtruck/delouis/PROD_RD12ll/DATA/PROD_RD12_RC4"
  MONTAGE  = "/softs/ImageMagick/6.8.8-1/bin/montage"
  SPICEEXE = "export OMP_NUM_THREADS=8; /wrk/hivon/tmp/Task_pkg/HL2_Spice/v03-01-000000_CVSHEAD/src/spice"
  DMCNAME  = None # deactivate DMC
  NODMCDATA = DMCDATA

elif SROLLHOST == "CC":
  DBPATH  = os.environ["DMCDATA"]   # /sps/planck/DMC/MISS03cc/DATA
  DX12DIR = "/sps/planck/SimuData/RD12_data/dmc/MISS03/DATA/DX12_MAP"

elif SROLLHOST == "EDISON":
  DBPATH  = os.environ["NODMCDATA"] # /scratch3/scratchdirs/paganol/RD12_data/dmc/MISS03/DATA
  DX12DIR = "/project/projectdirs/planck/data/compsep/exchanges/dx12/maps/hfi"

elif SROLLHOST == "CORI":
  DBPATH  = os.environ["NODMCDATA"] # /project/projectdirs/planck/data/hfi/RD12_data/dmc/MISS03/DATA
  DX12DIR = "/project/projectdirs/planck/data/compsep/exchanges/dx12/maps/hfi"

elif SROLLHOST == "OCCIGEN":
  DBPATH   = os.environ["NODMCDATA"] # /scratch/cnt0028/ias1717/SHARED/RD12_data/dmc/MISS03/DATA
  DX12DIR  = "/scratch/cnt0028/ias1717/SHARED/RD12_data/dmc/MISS03/DATA/DX12_MAP"
  MONTAGE  = "/usr/bin/montage"
  if HOSTNAME.startswith("n"):
    SPICEEXE = "export OMP_NUM_THREADS=24" # use all cores on compute nodes
  else:
    SPICEEXE = "export OMP_NUM_THREADS=8" # use 1/3 cores on login nodes
  SPICEEXE += "; /scratch/cnt0028/ias1717/SHARED/PolSpice_v03-05-01/bin/spice"

else:
  raise Exception( "Unsupported platform (yet): " + SROLLHOST)

BINTAB = numpy.array(( 0, 1, 2,    3,    4,    5,    6,    7,    8,    9,   10,   12,   14,
                16,   20,   25,   30,   35,   40,   45,   50,   60,   70,   80,
                90,  100,  110,  120,  130,  140,  150,  160,  180,  200,  220,
               240,  260,  280,  300,  320,  340,  360,  380,  400,  420,  440,
               460,  480,  500,  520,  540,  560,  580,  600,  620,  640,  660,
               680,  700,  720,  740,  760,  780,  800,  820,  840,  860,  880,
               900,  920,  940,  960,  980, 1000, 1050, 1100, 1150, 1200, 1250,
              1300, 1350, 1400, 1450, 1500, 1600, 1700, 1800, 1900, 2000, 2250,
              2500, 2750, 3000, 3250, 3500, 3750, 4000, 4250, 4500, 4750, 5000,
              5250, 5500, 5750, 6000, 6250))

PAGANOSCRATCH3 = "/scratch3/scratchdirs/paganol/RD12_data"
PLANCKPROJDIR  = "/project/projectdirs/planck/data/hfi/RD12_data"
MOTTETSCRATCH2 = "/scratch2/scratchdirs/smottet/RD12_data"
M3MISS03       = "/data"
M4MISS03       = "/pscratch1/RD12_data"

BEGINRINGINDEX = None
isDMC = (DMCNAME != None)

if isDMC:
  import piolib
else:
  import nodmclib


def ERI( ring):
  return BRI[ ring+1]-1


def badrings_dx11( pixname):
  boloid = BOLOID[pixname]
  freq   = int( pixname[0:3])
  if freq in (100, 143, 217):
    badrings_suffix = "_v64_extended"
  elif freq in (353, 545, 857):
    badrings_suffix = "_v61ter"
  return nodmclib.read_PIOINT( NODMCDATA + "/calROIs/%s_discarded_rings%s" % (boloid, badrings_suffix), 0, 27010)


def bad_bolos_per_ring( detset = ""):
  # return an array containing for each ring the number of bolometers for which the ring is bad
  badbolos = numpy.zeros(27006)
  if detset == "":
    pixnames = BOLOID.keys()
  else:
    pixnames = DETSETS[detset]
  for pixname in pixnames:
    badr = badrings_dx11( pixname)[0:27006]
    for i in range(27006):
      if badr[i] > 0:
        badbolos[i] += 1
  return badbolos


def totalflag_dx11_name( pixname):
  boloid = BOLOID[pixname]
  freq   = int( pixname[0:3])
  if freq in (100, 143, 217):
    flag_suffix = "_v64"
  elif freq == 353:
    flag_suffix = "_v68"
  elif freq in (545, 857):
    flag_suffix = "_v61"
  return NODMCDATA + "/calTOIs/%s_TotalFlag%s" % (boloid, flag_suffix)


def readBIN( objname, datatype, begin=-1, end=-1, begin_ring=-1, end_ring=-1, begin_HPR=-1, end_HPR=-1):
  if begin_ring != -1:
    first_sample = BRI[begin_ring]
  elif begin_HPR != -1:
    first_sample = begin_HPR * HPRSIZE
  else:
    first_sample = begin

  if end_ring != -1:
    last_sample = BRI[end_ring+1]-1
  elif end_HPR != -1:
    last_sample = (end_HPR+1) * HPRSIZE - 1
  else:
    last_sample = end

  sample_count = last_sample - first_sample + 1

  if datatype.endswith("32"):
    datasize = 4
  elif datatype.endswith("64"):
    datasize = 8
  elif datatype == "byte":
    datasize = 1
  else:
    raise Exception("Unknown datatype %s" % datatype)

  if not os.path.exists(objname):
    objname = "%s.%s.bin" % (objname, datatype)
    assert os.path.exists(objname)

  fobj = open( objname, "r")
  fobj.seek( first_sample * datasize)
  objdata = numpy.fromfile( fobj, dtype=datatype, count=sample_count)
  fobj.close()
  return objdata


def readHPRfloat( objname, begin_ring=-1, end_ring=-1) :
  fobj = open( objname, "r")
  fobj.seek( begin_ring * HPRSIZE * 4)
  objdata = numpy.fromfile( fobj, dtype = "float32", count = (end_ring-begin_ring+1) * HPRSIZE)
  fobj.close()
  return objdata


def readHPRdouble( objname, begin_ring=-1, end_ring=-1) :
  fobj = open( objname, "r")
  fobj.seek(begin_ring * HPRSIZE * 8)
  objdata = numpy.fromfile( fobj, dtype = "float64", count = (end_ring-begin_ring+1) * HPRSIZE)
  fobj.close()
  return objdata


def pbroll( phase, signal, flag):
  binidx = phase * PBRSIZE / 2 / numpy.pi # transform phase in radians to bin index
  sigpbr = numpy.zeros(PBRSIZE, dtype="float32")
  hitpbr = numpy.zeros(PBRSIZE, dtype="int32")
  for i in range(len(signal)):
    if flag[i] == 0:
      sigpbr[binidx[i]] += signal[i]
      hitpbr[binidx[i]] += 1
  for i in range(PBRSIZE):
    if hitpbr[i] != 0:
      sigpbr[i] /= hitpbr[i]
    else:
      return numpy.zeros(PBRSIZE, dtype="float32")
  return sigpbr


def pbunroll( phase, sigpbr):
  binidx = phase * PBRSIZE / 2 / numpy.pi # transform phase to bin index
  signal = numpy.zeros(len(phase), dtype="float32")
  for i in range(len(signal)):
    signal[i] = sigpbr[binidx[i]]
  return signal


def hpr2map( sigobj, pixname = None, begin_ring = 240, end_ring = 26050, remove_dipole = False, verbose = True):
  nside     = 2048
  ringbunch = 200
  ptgobj = NODMCDATA + "/PBR_JMD/%s_REP6_ptg" % pixname
  dipobj = NODMCDATA + "/PBR_JMD/%s_REP6_diporb" % pixname
  if pixname == None:
    pixname = get_pixname( sigobj)
  sigmap = numpy.zeros( 12 * nside * nside, dtype="float64")
  hitmap = numpy.zeros( 12 * nside * nside, dtype="int64")
  for ring in range( begin_ring, end_ring, ringbunch):
    if ring + ringbunch > end_ring:
      rings_to_read = end_ring - ring + 1
    else:
      rings_to_read = ringbunch
    if verbose:
      print( "%d-%d" % (ring, ring+rings_to_read-1))
    phi   = nodmclib.read_PIODOUBLE( ptgobj, ring * HPRSIZE, rings_to_read * HPRSIZE)
    theta = nodmclib.read_PIODOUBLE( ptgobj + "_TUPLE_1", ring * HPRSIZE, rings_to_read * HPRSIZE)
    sig = nodmclib.read_PIOFLOAT( sigobj, ring * HPRSIZE, rings_to_read * HPRSIZE)
    if remove_dipole:
      dip = nodmclib.read_PIOFLOAT( dipobj, ring * HPRSIZE, rings_to_read * HPRSIZE)
      sig[dip != 0] -= dip[dip != 0]
#    assert numpy.all((theta!=0) == (phi!=0))
    pixnums = healpy.ang2pix( nside, theta[theta != 0], phi[phi != 0])
    sigmap[pixnums] += sig[theta != 0]
    hitmap[pixnums] += 1
  sigmap[hitmap != 0] /= hitmap[hitmap != 0]
  sigmap[hitmap == 0] = numpy.nan
  return sigmap


# computes the noise power spectra per ring of a TOI, use meanVEC() to compute the mean detector noise power spectrum
# invoke (in a batch script) with something like:
# mpiexec -np 8 python -c "import stools; stools.detnoiseMPI('{toiname}')"
# takes about 9 hours for 1 full TOI on 1 M3 node with 8 MPI ranks
def detnoiseMPI( objname, pixname=None, begin_ring=240, end_ring=26050, outfile=None, flag=None, bad_rings=True, FFTSIZE=2**16):
  if pixname == None:
    pixname = get_pixname( objname)
  if outfile == None:
    outfile = objname + "_dtnoise"
  powspec_size = FFTSIZE/2+1
  if bad_rings:
    badr = badrings_dx11( pixname)
  good_rings = 0
  bad_rings = 0
  if MPI_RANK == 0:
    os.system("rm -f {outfile}; touch {outfile}".format( **locals()))
  else:
    time.sleep(1)
  for ring in range( begin_ring + MPI_RANK, end_ring + 1, MPI_SIZE):
    print "%d/%d: ring %d" % (MPI_RANK, MPI_SIZE, ring)
    if bad_rings and (badr[ring] != 0):
      bad_rings += 1
      powspec = numpy.zeros( powspec_size, dtype="float64")
      pwrite( outfile, powspec_size*ring, powspec, datatype="float64")
      continue
    phase = nodmclib.read_PIODOUBLE( DMCDATA + "/Sa_HFI_C_Bolo/Phase_dx11", BRI[ring], BRI[ring+1] - BRI[ring])
    signal = nodmclib.read_PIOFLOAT( objname, BRI[ring], BRI[ring+1] - BRI[ring])
    if flag != None:
      sigflag = nodmclib.read_PIOFLAG( flag, BRI[ring], BRI[ring+1] - BRI[ring])
    else:
      sigflag = numpy.zeros( len( signal), dtype="byte")
    sigpbr = pbroll( phase, signal, sigflag)
    if numpy.sum( sigpbr) != 0:
      noise  = signal - pbunroll( phase, sigpbr)
      powspec = numpy.abs( numpy.fft.rfft( noise[:FFTSIZE]))
      good_rings += 1
    else:
      powspec = numpy.zeros( powspec_size, dtype="float64")
      bad_rings += 1
    pwrite( outfile, powspec_size*ring, powspec, datatype="float64")
  print( "good rings = %d ; bad rings = %d" % (good_rings, bad_rings))
  return good_rings


def meanVEC( vecname, ringsize, begin_ring=240, end_ring=26050, no0=True, datatype="float64"):
  mimine = numpy.zeros( ringsize, dtype="float64")
  nrings = 0
  for ring in range( begin_ring, end_ring+1):
    if datatype == "float64":
      toto = nodmclib.read_PIODOUBLE( vecname, ring*ringsize, ringsize)
    elif datatype == "float32":
      toto = nodmclib.read_PIOFLOAT( vecname, ring*ringsize, ringsize)
    if no0 and numpy.sum( toto) == 0.0:
      continue
    nrings += 1
    mimine += toto
  print "meanVEC() used %d/%d rings" % (nrings, end_ring-begin_ring+1)
  return mimine / nrings, nrings


def lincombfit( vec, ref, offset=1):
  nvec   = len(vec)
  veclen = len(ref)

  a = numpy.zeros( (nvec + offset, nvec + offset), dtype="float64")
  b = numpy.zeros( nvec + offset, dtype="float64")

  for i in range( 0, nvec):
    assert len( vec[i]) == veclen
    for j in range( 0, nvec):
      a[i,j] = numpy.sum( vec[i] * vec[j])
    b[i] = numpy.sum( vec[i] * ref)
    if offset:
      a[i, nvec] = numpy.sum( vec[i])
      a[nvec, i] = a[i, nvec]
  if offset:
    a[nvec, nvec] = veclen
    b[nvec] = numpy.sum( ref)
  return numpy.linalg.solve( a, b)


def spice_diff( fits1, fits2, weight=None, mask=None, lmax=None, cl_filename=None, thetamax=False):
  """save fits1-fits2 in a temporary fits file and spice it,
if <weight> is given, diff map is multiplied by <weight>,
if <mask> is given, it is passed to spice"""
  f1, fitsname = tempfile.mkstemp()
  os.close( f1)

  hdulist = pyfits.open( fits1)
  field_count = len(hdulist[1].data[0])
  hdulist.close()
  assert (field_count == 1) or (field_count == 3)

  if weight != None:
    wmap = healpy.read_map( weight)

  maps = list()
  for field in range( field_count):
    map1 = healpy.read_map( fits1, field=field)
    map2 = healpy.read_map( fits2, field=field)
    mapdiff = map1-map2
    if weight != None:
      mapdiff *= wmap
    mapdiff[map1 == healpy.UNSEEN] = healpy.UNSEEN
    mapdiff[map2 == healpy.UNSEEN] = healpy.UNSEEN
    maps.append( mapdiff)

  if field_count == 3:
    healpy.write_map( fitsname, maps)
    pol = "YES"
  elif field_count == 1:
    healpy.write_map( fitsname, maps[0])
    pol = "NO"

  try:
    l_cl = spice_fits( fitsname, mask=mask, lmax=lmax, cl_filename=cl_filename, pol=pol, thetamax=thetamax)
  except:
    os.remove( fitsname) # cleaning before failing
    raise
  os.remove( fitsname)
  return l_cl


def spice_data( mapI, mapQ=None, mapU=None, mask=None, lmax=None, cl_filename=None, thetamax=False):
  """save mapI, mapQ, mapU in a temporary fits file and spice it"""
  f1, fitsname = tempfile.mkstemp()
  os.close( f1)
  if mapQ != None:
    assert mapU != None
    healpy.write_map( fitsname, (mapI, mapQ, mapU))
    pol = "YES"
  else:
    assert mapU is None
    healpy.write_map( fitsname, mapI)
    pol = "NO"

  try:
    l_cl = spice_fits( fitsname, mask=mask, lmax=lmax, cl_filename=cl_filename, pol=pol, thetamax=thetamax)
  except:
    os.remove( fitsname) # cleaning before failing
    raise
  os.remove( fitsname)
  return l_cl


def spice_fits( fits1, fits2=None, mask=None, lmax=None, cl_filename=None, pol=None, thetamax=False):
  """http://www2.iap.fr/users/hivon/software/PolSpice/"""
  if cl_filename is None:
    f2, txtclname = tempfile.mkstemp()
    os.close( f2)
  else:
    txtclname = cl_filename

  if pol == None:
    hdulist = pyfits.open( fits1)
    if len(hdulist[1].data[0]) == 1:
      pol = "NO"
    elif len(hdulist[1].data[0]) == 3:
      pol = "YES"
    else:
      raise Exception( "spice_fits(): fits file should contain either 1 or 3 columns (%d found)" % len(hdulist[1].data[0]))

  if (cl_filename is None) or (not os.path.exists( txtclname)):
    if not lmax is None:
      nlmax = " -nlmax %d " % lmax
    else:
      nlmax = ""

    SPICECMD = SPICEEXE + " -mapfile %s -clfile %s -polarization %s -subdipole YES %s -verbosity 0" % (fits1, txtclname, pol, nlmax)
    if fits2 != None:
      SPICECMD += " -mapfile2 %s " % fits2
    if mask != None:
      SPICECMD += " -weightfile %s " % mask
      if thetamax:
        # http://www2.iap.fr/users/hivon/software/PolSpice/faq.html
        # apodizesigma and thetamax should be set to 180 when fsky > 50%
        # apodizesigma and thetamax should be set to 20 when fsky = 1%
        SPICECMD += " -apodizesigma %s -thetamax %s" % (str(thetamax), str(thetamax))
      if pol == "YES":
        SPICECMD += " -decouple YES -tolerance 1e-6 "

    print SPICECMD
    try:
      os.system( SPICECMD)
    except:
      print "spice error"
      if cl_filename is None:
        os.remove( txtclname) # cleaning before failing
      raise
  else:
    # if C(l) text file already exists, don't run spice
    print "spice_fits(): reading " + txtclname

  # if 2 polarised fits files, return (l, tt, ee, bb, te, tb, eb, et, bt, be)
  # if 1 polarised fits file,  return (l, tt, ee, bb, te, tb, eb)
  # if 1 not polarised, fits file return (l, cl)
  lcl = numpy.loadtxt( txtclname, comments="#", unpack=True)

  if cl_filename is None:
    os.remove( txtclname)

  return lcl


# smooth an array with a convolved box
def conv_smooth(y, box_pts):
  box = numpy.ones(box_pts)/box_pts
  y_smooth = numpy.convolve(y, box, mode='same')
  return y_smooth


# fill map_in NANs and UNSEENs with values of map_in downgraded to nside_down
def inpaint_map( map_in, nside_down=128):
  nside_in = healpy.npix2nside( len( map_in))
  map_in[map_in == healpy.UNSEEN] = numpy.nan
  unseen = numpy.where( numpy.isnan( map_in))[0]
  while len( unseen) > 0:
    print( "inpaint_map(): filling with map downgraded to nside=%d" % nside_down)
    mapfill = healpy.ud_grade( map_in, nside_down)
    mapfill = healpy.ud_grade( mapfill, nside_in)
    map_in[unseen] = mapfill[unseen]
    unseen = numpy.where( numpy.isnan( map_in))[0]
    nside_down /= 2
  return map_in


def fix_planck_permissions( directory):
  cmd = """find {directory} -type d -exec chmod 750 {{}} \;
           find {directory} -type f -exec chmod 640 {{}} \;
           chgrp -R planck {directory}""".format( **locals())
  # {} will be replaced by find with each filename found,
  # {{}} is needed to escape curved braces in python strings,
  # \; executes one chmod per file,
  # use + instead of \; to execute one chmod for all the files
  os.system( cmd)


def md5sum_dmcobj( backendname):
  # computes the md5sum of each .pio and .pio.flag file and then the md5sum of the sorted list of md5sums
  cmd = "find %s -type f -exec md5sum {} + | awk '{print $1}' | sort | md5sum" % backendname
  toto = subprocess_check_output( cmd, shell=True)
  titi = toto.split()[0]
  assert len( titi) == 32, "unexpected result in md5sum: %s" % toto
  return titi


def ring_number( sample_number):
  for ring in range( len( BRI)):
    if BRI[ring] > sample_number:
      return ring-1
  return -1


def read_metadata_txt( mdname):
  if not os.path.exists( mdname):
    print( "stools.read_metadata_txt(): metadata text file not found: " + mdname)
    return
  with open( mdname, 'r') as f:
    metadata = f.readlines()
  piotype = ""
  datatype = ""
  begidx = -1
  endidx = -1
  author = ""
  creadate = ""
  backendname = ""
  for md in metadata:
    if md.startswith("TOItype"):
      piotype = md.split(":")[-1].strip()
    if md.startswith("Datatype"):
      datatype = md.split(":")[-1].strip()
    if md.startswith("BeginIndex"):
      begidx = int( md.split(":")[-1].split()[0].strip())
    if md.startswith("EndIndex"):
      endidx = int( md.split(":")[-1].split()[0].strip())
    if md.startswith("Author"):
      author = md.split(":")[-1].strip()
    if md.startswith("Date"):
      creadate = md.split(":")[-1].strip()
    if md.startswith("Backendname"):
      backendname = md.split(":")[-1].strip()
    if piotype == "POI":
      bri = begidx / HPRSIZE
      eri = endidx / HPRSIZE
    elif piotype == "TOI":
      bri = ring_number( begidx)
      eri = ring_number( endidx)
    else:
      bri = -1
      eri = -1
  return (  piotype, datatype, begidx, endidx, bri, eri, author, creadate, backendname)


def get_dmc_metadata( objname):
  if os.path.isdir( objname):
    if os.path.exists( objname + "/no_dmc_metadata.txt"):
      return read_metadata_txt( objname + "/no_dmc_metadata.txt")
    else:
      print "stools.get_dmc_metadata(): no_dmc_metadata.txt not found for DMCobject: " + objname
      return
  elif os.path.exists( objname + ".meta"):
    return read_metadata_txt( objname + ".meta")
  else:
    print "stools.get_dmc_metadata(): " + objname + ".meta not found"
    return


def get_backendname( objname):
  if not os.path.exists( objname):
    return ""
  mdname = objname + "/no_dmc_metadata.txt"
  if not os.path.exists( mdname):
    return ""
  with open( mdname, 'r') as f:
    metadata = f.readlines()
  for md in metadata:
    if md.startswith("Backendname"):
      return md.strip().split()[-1]
  return ""


def fix_backendname( objname, oldpath, newpath):
  if not os.path.exists( objname):
    print( "object not found: " + objname)
    return
  mdname = objname + "/no_dmc_metadata.txt"
  if not os.path.exists( mdname):
    print( "no_dmc_metadata not found: " + objname)
    return
  f = open( mdname, 'r')
  metadata = f.read()
  f.close()
  if metadata.find( oldpath) != -1:
    metadata = metadata.replace( oldpath, newpath)
    f = open( mdname, 'w')
    f.write( metadata)
    f.flush()
    f.close()
  else:
    print "  backendname prefix not found: " + objname
  return


def bin_cl(cl):
  # return an array whith same length as input cl, but with C(l) averaged per BINTAB bin
  assert len(cl)<BINTAB[-1]
  binned = numpy.zeros(len(cl))
  i=0
  while BINTAB[i+1] < len(cl):
    binned[BINTAB[i]:BINTAB[i+1]] = numpy.sum(cl[BINTAB[i]:BINTAB[i+1]]) / (BINTAB[i+1]-BINTAB[i])
    i += 1
  binned[BINTAB[i]:len(cl)] = numpy.sum(cl[BINTAB[i]:len(cl)]) / (len(cl)-BINTAB[i])
  return binned


def bin_l_cl(cl):
  # return one l and one C(l) per BINTAB bin
  assert len(cl)<BINTAB[-1]
  binned = numpy.zeros(len(BINTAB))
  i=0
  while BINTAB[i+1] < len(cl):
    binned[i] = numpy.sum(cl[BINTAB[i]:BINTAB[i+1]]) / (BINTAB[i+1]-BINTAB[i])
    i += 1
  binned[i] = numpy.sum(cl[BINTAB[i]:len(cl)]) / (len(cl)-BINTAB[i])
  return BINTAB[0:i+1]+(BINTAB[1:i+2]-BINTAB[0:i+1])/2, binned[:i+1]


def lin_bin_cl( cl, binsize):
  binnedlen = int( numpy.ceil( len( cl) / binsize))
  binned = numpy.zeros( binnedlen)
  for i in range( binnedlen):
    binned[i] = numpy.sum( cl[ i * binsize:(i+1) * binsize]) / binsize
  return numpy.arange( binnedlen) * binsize + binsize/2, binned


# find a pixname in a file name / string
def get_pixname( objname):
  for pixname in BOLOID.keys():
    if pixname in objname:
      return pixname
  for boloid in BOLOID.values():
    if boloid in objname:
      return boloid[3:].replace("_", "-")
  for pixname in BOLOID.values():
    # PSM style...
    if pixname.replace("-", "_") in objname:
      return pixname
  return ""


# minimises (x - gain * y)
# see also: slope, intercept, r_value, p_value, std_err = scipy.stats.linregress( x, y) # y = intercept + slope*x
def gain( x, y):
  return numpy.dot( x, y ) / numpy.dot( y, y )


def pwrite( filename, offset, data, datatype="float32"):
# python 3.3: os.pwrite( filedesc, bytestring, offset)
# file must exist, you can non-destructively create it with "open( filename, 'a').close()"
# TODO: consider locking the file with fcntl.lockf()
  if (datatype == "float32") or (datatype == "f32") or (datatype == "PIOFLOAT"):
    array_datatype = "f"
    datasize = 4
  elif (datatype == "float64") or (datatype == "f64") or (datatype == "PIODOUBLE"):
    datasize = 8
    array_datatype = "d"
  elif (datatype == "int32") or (datatype == "i32") or (datatype == "PIOINT"):
    datasize = 4
    array_datatype = "i"
  elif (datatype == "int64") or (datatype == "i64") or (datatype == "PIOLONG"):
    datasize = 8
    array_datatype = "l"
  outfile = open( filename, 'r+b')
  outfile.seek( offset * datasize)
  out_array = array.array( array_datatype, data)
  assert out_array.itemsize == datasize # itemsize may vary depending on the platform
  out_array.tofile( outfile)
  outfile.close()


def hpr2toi( hprname, toiname, begring=240, endring=26050, delete_before=False):
  pixname = get_pixname( hprname)
  pixidxname = DBPATH + "/e2e_common_TOI/{pixname}_HPRIDX_ABER_TotalFlag_dx11".format( **locals())
  if delete_before:
    try: piolib.DeleteObject( toiname)
    except: pass
  if not piolib.CheckObject( toiname):
    piolib.CreateTOIObject( toiname, "PIOFLOAT")
  print toiname

  piohpr = piolib.CheckObject( hprname)

  for ring in range( begring, endring):
    begsamp = BRI[ring]
    endsamp = BRI[ring+1]-1
    ringcmd = "ring=%d" % (ring)
    sampcmd = "begin=%d;end=%d" % (begsamp, endsamp)
    if (ring==begring) or (ring==endring) or (ring%1000==0):
      print ringcmd
    if piohpr:
      hpr = piolib.read( hprname, command=ringcmd)
    else:
      hpr = readHPRfloat( hprname, begin_ring=ring, end_ring=ring)
    pixidx = piolib.read( pixidxname, command=sampcmd)
    toi = hpr[pixidx]
    toi[pixidx == -1] = 0.0
    piolib.WriteTOIObject( toi, toiname, "PIOFLOAT", sampcmd)


def compare_HPR_hits( hpr1, hpr2, verbose=False):
  badrings = list()
  badspl = 0
  for ring in range( BEGMISS, ENDMISS+1):
    diff0 = numpy.sum( (hpr1[ring*HPRSIZE:(ring+1)*HPRSIZE] == 0) != (hpr2[ring*HPRSIZE:(ring+1)*HPRSIZE] == 0))
    if diff0 != 0:
      badrings.append( ring)
      badspl += diff0
  if verbose:
    if badspl == 0:
      print "same hits"
    else:
      print "%d differing hits in %d rings" % (badspl, len( badrings)), badrings[:10]
  return (badspl, badrings)


def map2hpr( pixname): # ~10 minutes per bolo on a M3 compute node
  #assert not "Comment me"
  map_name = "/redtruck/SimuData/FFP10_gaussbeam_sky/%s%s.fits" % (pixname, GBTAG)
  (mapi, mapq, mapu) = healpy.read_map( map_name, (0, 1, 2))
  epsilon = XPOL[pixname]
  objname = "/data/dmc/MISS03/DATA/cleanbeam_VEC/%s%s_HPR" % (pixname, GBTAG)
  if piolib.CheckObject( objname):
    piolib.DeleteObject( objname)
  piolib.CreateVECTObject( objname, "PIOFLOAT")
  print objname
  lastring = 26050
  ringbunch = 100
  for ring in range( 240, lastring+1, ringbunch):
    end_ring = min( lastring, ring+ringbunch-1)
    command = "begin=%d;end=%d" % (ring*HPRSIZE, (end_ring+1)*HPRSIZE-1)
    print "begin_ring=%d;end_ring=%d" % (ring, end_ring)
    phi = piolib.read( "/data/dmc/MISS03/DATA/PBR_JMD/%s_REP6_ptg" % pixname,         command=command)
    the = piolib.read( "/data/dmc/MISS03/DATA/PBR_JMD/%s_REP6_ptg_TUPLE_1" % pixname, command=command)
    psi = piolib.read( "/data/dmc/MISS03/DATA/PBR_JMD/%s_REP6_ptg_TUPLE_2" % pixname, command=command)
    badpix = numpy.where( numpy.logical_and(the == 0.0, phi == 0.0))[0]
#    print "badpixels: %d (%f%%)" % (len(badpix), float(len(badpix))/len(phi)*100)
    pixels  = healpy.ang2pix( 2048, the, phi)
    cos2psi = numpy.cos( 2.0 * psi)
    sin2psi = numpy.sin( 2.0 * psi)
    s = mapi[pixels] + (1.0 - epsilon) / (1.0 + epsilon) * (cos2psi * mapq[pixels] + sin2psi * mapu[pixels])
    s[badpix] = 0.0
    assert len(s) == len(phi)
    piolib.WriteVECTObject( s, objname, "PIOFLOAT", command)


def dipole2map( amplitude=3362.08*1e-6, deg_ra=264.021, deg_dec=48.253, nside=2048):
  # default values are HFI 2018 DPC paper in KCMB, see http://pmwiki.sylvainmottet.fr/index.php?n=Main.DipoleValues
  ra  = numpy.radians( deg_ra)
  dec = numpy.radians( deg_dec)
  dipref = numpy.array([numpy.cos(ra) * numpy.cos(dec) , numpy.sin(ra) * numpy.cos(dec), numpy.sin(dec)])
  vec = healpy.pix2vec( nside, numpy.arange( healpy.nside2npix(nside)))
  return numpy.dot( dipref, vec) * amplitude


def fitmap( data, model, mask=None):
  """ return the gain g that minimizes (data-g*model)"""
  npix = len( data)
  assert numpy.log2( numpy.sqrt( npix/12.0))%1.0 == 0.0, "input <data> has an invalid number of pixels (not a Healpix map?)"
  assert len( model) == npix
  if mask is None:
    mask = numpy.ones( npix)
  else:
    assert len( mask) == npix
  # replace UNSEEN and NaN with median of map finite values
  data[numpy.logical_not( data>healpy.UNSEEN)] = numpy.median( data[data>healpy.UNSEEN])
  mat  = [[numpy.sum( mask*model*model), numpy.sum( mask*model)], [numpy.sum( mask*model), numpy.sum( mask)]]
  imat = numpy.linalg.inv( mat)
  vec  = [numpy.sum( mask*model*data), numpy.sum( mask*data)]
  rr   = numpy.dot( vec, imat)
  return rr[0]


def inv_mat3x3( m):
  """return the inverse of matrix m"""
  if hasattr( m, "type"):
    dtype = m.type
  else:
    dtype = "float32"
  m1, m2, m3, m4, m5, m6, m7, m8, m9 = m.ravel()
  det = m1*m5*m9 + m4*m8*m3 + m7*m2*m6 - m1*m6*m8 - m3*m5*m7 - m2*m4*m9
  return numpy.array([m5*m9-m6*m8, m3*m8-m2*m9, m2*m6-m3*m5, m1*m9-m3*m7, m3*m4-m1*m6, m1*m5-m2*m4], dtype=dtype) / det


def inv_sym_mat3x3( m11, m12, m13, m22, m23, m33):
  """return the upper triangular half of the inverse of the symmetrical matrix:
m11 m12 m13
m12 m22 m23
m13 m23 m33
  example: (covii, coviq, coviu, covqq, covqu, covuu) = inv_sym_mat3x3( ii, iq, iu, qq, qu, uu)"""
  if hasattr( m11, "type"):
    dtype = m11.type
  else:
    dtype = "float32"
  print dtype
  det = m11*m22*m33 + m12*m23*m13 + m13*m12*m23 - m11*m23*m23 - m13*m22*m13 - m12*m12*m33
  # caller should catch divide by 0 errors...
  return numpy.array([m22*m33-m23*m23, m13*m23-m12*m33, m12*m23-m13*m22, m11*m33-m13*m13, m13*m12-m11*m23, m11*m22-m12*m12], dtype=dtype) / det

def downgrade_map( i, fwhm_deg, nsideout):
  """fwhm_deg: FWHM in degrees of the input maps effective beam, it will be deconvolved from the input map"""
  assert nsideout in [32, 64]
  npix = len(i)
  nsidein = healpy.npix2nside( npix)
  wpixin  = healpy.pixwin( nsidein,  pol=False) #, lmax=nsideout*4)
  wpixout = healpy.pixwin( nsideout, pol=False) #, lmax=nsideout*4)
  beamin  = healpy.gauss_beam( numpy.deg2rad( fwhm_deg), lmax=nsideout*4, pol=False)

  l = numpy.arange( numpy.size( wpixout[0]))
  if nsideout <= 32:
#    print "smooth map with ns=%d cosinus window function" % nsideout
    beamout = 0.5 * (1.0 + numpy.cos( numpy.pi * (l-nsideout) / 2.0 / nsideout))
    beamout[0:nsideout+1] = 1.0
    beamout[3*nsideout+1:] = 0.0
  else:
    # smooth the map with a gaussian beam of FWHM the size of a nsideout pixel
    fwhm_out = healpy.max_pixrad( nsideout, degrees=True) * 2 * 60
#    print "smooth map with %.1f arcmin FWHM gaussian beam" % fwhm_out
    beamout = healpy.gauss_beam( numpy.deg2rad( fwhm_out / 60.0), lmax=nsideout*4, pol=False)

  med = numpy.median( i[i != healpy.UNSEEN])
  i[i == healpy.UNSEEN] = med

  alm = healpy.map2alm( i, lmax=nsideout*4)
  healpy.almxfl( alm, beamout / beamin / wpixin[0:nsideout*4+1] * wpixout, inplace=True)
  m = healpy.alm2map( alm, nsideout, lmax=nsideout*4, pixwin=False, pol=True)
  return m


def downgrade_polmap( (i, q, u), fwhm_deg, nsideout):
  """fwhm_deg: FWHM in degrees of the input maps effective beam, it will be deconvolved from the input map"""
  assert nsideout in [32, 64]
  npix = len(i)
  assert len(q) == npix
  assert len(u) == npix
  nsidein = healpy.npix2nside( npix)
  wpixin  = healpy.pixwin( nsidein,  pol=True) #, lmax=nsideout*4)
  wpixout = healpy.pixwin( nsideout, pol=True) #, lmax=nsideout*4)
  beamin  = healpy.gauss_beam( numpy.deg2rad( fwhm_deg, lmax=nsideout*4, pol=True))

  l = numpy.arange( numpy.size( wpixout[0]))
  if nsideout <= 32:
#    print "smooth map with ns=%d cosinus window function" % nsideout
    beamout = 0.5 * (1.0 + numpy.cos( numpy.pi * (l-nsideout) / 2.0 / nsideout))
    beamout[0:nsideout+1] = 1.0
    beamout[3*nsideout+1:] = 0.0
    beamout = numpy.array( [beamout, beamout, beamout])
  else:
    # smooth the map with a gaussian beam of FWHM the size of a nsideout pixel
    fwhm_out = healpy.max_pixrad( nsideout, degrees=True) * 2 * 60
#    print "smooth map with %.1f arcmin FWHM gaussian beam" % fwhm_out
    beamout = healpy.gauss_beam( numpy.deg2rad( fwhm_out / 60.0), lmax=nsideout*4, pol=True)

  med = numpy.median( i[i != healpy.UNSEEN])
  i[i == healpy.UNSEEN] = med
  med = numpy.median( q[q != healpy.UNSEEN])
  q[q == healpy.UNSEEN] = med
  med = numpy.median( u[u != healpy.UNSEEN])
  u[u == healpy.UNSEEN] = med

  alm = healpy.map2alm( [i, q, u], lmax=nsideout*4)
  healpy.almxfl( alm[0], beamout[:,0] / beamin[:,0] / wpixin[0][0:nsideout*4+1] * wpixout[0], inplace=True)
  healpy.almxfl( alm[1], beamout[:,1] / beamin[:,1] / wpixin[1][0:nsideout*4+1] * wpixout[1], inplace=True)
  healpy.almxfl( alm[2], beamout[:,2] / beamin[:,2] / wpixin[1][0:nsideout*4+1] * wpixout[1], inplace=True)
  m = healpy.alm2map( alm, nsideout, lmax=nsideout*4, pixwin=False, pol=True)
  return m


def change_coord(m, coord):
    """ Change coordinates of a HEALPIX map
    (https://stackoverflow.com/questions/44443498/how-to-convert-and-save-healpy-map-to-different-coordinate-system)

    Parameters
    ----------
    m : map or array of maps
      map(s) to be rotated
    coord : sequence of two character
      First character is the coordinate system of m, second character
      is the coordinate system of the output map. As in HEALPIX, allowed
      coordinate systems are 'G' (galactic), 'E' (ecliptic) or 'C' (equatorial)

    Example
    -------
    The following rotate m from galactic to equatorial coordinates.
    Notice that m can contain both temperature and polarization.
    >>>> change_coord(m, ['G', 'C'])
    """
    # Basic HEALPix parameters
    npix = m.shape[-1]
    nside = healpy.npix2nside(npix)
    ang = healpy.pix2ang(nside, numpy.arange(npix))

    # Select the coordinate transformation
    rot = healpy.Rotator(coord=reversed(coord))

    # Convert the coordinates
    new_ang = rot(*ang)
    new_pix = healpy.ang2pix(nside, *new_ang)

    return m[..., new_pix]
