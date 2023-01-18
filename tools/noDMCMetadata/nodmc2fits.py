#!/usr/bin/env python

import os
import sys
import glob
import numpy
import tempfile
try:
  # pyfits is not installed at CC
  import astropy.io.fits as pyfits
except:
  import pyfits

import nodmclib
import stools

# numpy dtype for conversion to FITS data. '>' means big endian, which is the endiannes in FITS files
NPDTYPE = { "PIOFLOAT":  ">f4",
            "PIODOUBLE": ">f8",
            "PIOINT":    ">i4",
            "PIOLONG":   ">i8",
            "PIOFLAG":   "B",
}


# FITS TFORM codes: https://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/docs/cfitsio.pdf p.16
TYPETFORM = { "PIOFLOAT":  "1E",
              "PIODOUBLE": "1D",
              "PIOINT":    "1J",
              "PIOLONG":   "1K",
              "PIOFLAG":   "1B",
}

################################################################################

def create_streaming_fits( fits_name, fits_hdu):

  hdu0 = pyfits.PrimaryHDU()
  hdu0.header.tofile( fits_name, clobber=True)
  f, tmpname = tempfile.mkstemp()
  os.close( f)
  fits_hdu.header.tofile( tmpname, clobber=True)

  fits_file = open( fits_name, 'rb+')
  with open( tmpname, 'rb') as tmpfits:
    fits_file.seek( 0, 2)
    fits_file.write( tmpfits.read())
  os.remove( tmpname)
  return fits_file


################################################################################

def write_streaming_fits( fits_file, data, piotype):

  arr = numpy.array( data, dtype=NPDTYPE[ piotype])
  arr.tofile( fits_file)
  return


################################################################################

def close_streaming_fits( fits_file):

  if fits_file.tell() % 2880 != 0:
    s_temp = '\00'*(2880 - (fits_file.tell() % 2880))
    fits_file.write( s_temp)
  fits_file.close()
  return


################################################################################

def usage_and_abort():
  print "%s: convert a dmc object to a fits file using nodmclib.py" % sys.argv[0]
  print "    takes about 30min for a PIOFLOAT TOI on magique4 bluevan"
  print "usage: %s <object(s)>" % sys.argv[0]
  print "       <object(s)>: full object name(s), accepts wildcards"
  exit(1)


################################################################################

try:
  SROLLHOST = os.environ["SROLLHOST"]
except:
  raise Exception( "You must 'source ./srollex_setenv.sh' before using this script")

if SROLLHOST == "M3":
  BINOUTDIR = "/redtruck/SimuData/nodmcfits"
if SROLLHOST == "M4":
  BINOUTDIR = "/bluevan1/symottet"
elif SROLLHOST == "CC":
  BINOUTDIR = "/sps/planck/SimuData/nodmcfits"
elif SROLLHOST == "CORI":
  BINOUTDIR = "/global/cscratch1/sd/smottet"
else:
  # there should not be objects left in DMC format on other clusters
  raise Exception( "Unsupported SROLLHOST (%s)" % SROLLHOST)

DMCDATA = os.environ["DMCDATA"]

# computing the md5sum of a stream requires chunks of 4096 bytes
BUFLEN = 25000 * 4096

if len( sys.argv) == 1:
  usage_and_abort()

params = sys.argv[1:]

for param in params:
  for objname in glob.glob( param):
#    fitsname = objname.replace( DMCDATA, BINOUTDIR) + ".fits"
    fitsname = objname + ".fits"
    if os.path.exists( fitsname):
      print objname + ": already converted, skipping\n"
      continue

    # check if binary file
    if os.path.isfile( objname):
      if os.path.isfile( objname + ".meta"):
        print "reading " + objname
        PIOType, DataType, BegIdx, EndIdx, BegRing, EndRing, Author, Date, backname = stools.read_metadata_txt( objname + ".meta")
      else:
        print objname + ": no metadata found for binary file, sikipping\n"
        continue
    else:
      # check if DMC file
      if not objname.startswith( DMCDATA):
        print objname + ": not in DMCDATA path, skipping\n"
        continue
      if objname.endswith( "Written"):
        print objname + " skipping\n"
        continue
      if os.path.exists( objname + "/no_dmc_metadata.txt"):
        print "reading " + objname
        PIOType, DataType, BegIdx, EndIdx, BegRing, EndRing, Author, Date = nodmclib.getObjectInfo( objname)
      else:
        print objname + ": no_dmc_metadata.txt not found, skipping\n"
        continue

    if not (DataType in ["PIOINT", "PIOLONG", "PIOFLOAT", "PIODOUBLE", "PIOFLAG"]):
      print objname + ": %s data type not implemented, skipping\n" % DataType
      continue

    # preparing FITS header
    col = [pyfits.Column( name=PIOType+"data", format=TYPETFORM[DataType])]
    if hasattr( pyfits, "new_table"):
      hdu = pyfits.new_table( col)
    else:
      hdu = pyfits.BinTableHDU.from_columns(col)

    hdu.header['NAXIS2'] = EndIdx - BegIdx + 1
    hdu.header['HIERARCH first_sample'] = BegIdx
    if PIOType == "PIOFLAG":
      hdu.header['PIOTYPE']= "PIOBYTE"
    else:
      hdu.header['PIOTYPE']= PIOType
    hdu.header['DATATYPE'] = DataType
    hdu.header['BEGIDX']   = BegIdx
    hdu.header['ENDIDX']   = EndIdx
    hdu.header['BEGRING']  = BegRing
    hdu.header['ENDRING']  = EndRing
    hdu.header['AUTHOR']   = Author
    hdu.header['DMCOBJ']   = objname
    hdu.header['DATE']     = Date

    # creating FITS file
    print "writing " + fitsname
    os.system( "mkdir -p " + os.path.dirname( fitsname))
    fits_file = create_streaming_fits( fitsname, hdu)

    # copying data
    for idx in range( BegIdx, EndIdx, BUFLEN):
      readlen = BUFLEN
      if idx + readlen - 1 > EndIdx:
        readlen = EndIdx - idx + 1
      print "r",; sys.stdout.flush()
      if DataType == "PIOFLOAT":
        buf = nodmclib.read_PIOFLOAT( objname, idx, readlen)
      elif DataType == "PIODOUBLE":
        buf = nodmclib.read_PIODOUBLE( objname, idx, readlen)
      elif DataType == "PIOINT":
        buf = nodmclib.read_PIOINT( objname, idx, readlen)
      elif DataType == "PIOLONG":
        buf = nodmclib.read_PIOLONG( objname, idx, readlen)
      elif DataType == "PIOFLAG":
        buf = nodmclib.read_PIOFLAG( objname, idx, readlen)
      print "w",; sys.stdout.flush()
      write_streaming_fits( fits_file, buf, DataType)
    close_streaming_fits( fits_file)
    print "done\n"



"""
login0:$ $SCRATCHDIR/FitsInfo.py /scratch/cnt0028/ias1717/SHARED/bware/wrk_sroll21/LevelS/tmpdat/mmtoi_ffp10_lensed_scl_cmb_100.fits
========== file header ==========
Filename: mmtoi_ffp10_lensed_scl_cmb_100.fits
No.    Name      Ver    Type      Cards   Dimensions   Format
  0  PRIMARY       1 PrimaryHDU      45   ()
  1  xtension      1 BinTableHDU     19   2934314859R x 1C   [1E]

========== HDU[0] header ==========
         SIMPLE =            True / file does conform to FITS standard
         BITPIX =               8 / number of bits per data pixel
          NAXIS =               0 / number of data axes
         EXTEND =            True / FITS dataset may contain extensions
        COMMENT =   FITS (Flexible Image Transport System) format is defined in 'Astronomy
  and Astrophysics', volume 376, page 359; bibcode: 2001A&A...376..359H /
        COMMENT =   FITS (Flexible Image Transport System) format is defined in 'Astronomy
  and Astrophysics', volume 376, page 359; bibcode: 2001A&A...376..359H /
           DATE = 2019-08-23T05:11:08 / file creation date (YYYY-MM-DDThh:mm:ss UT)
 beam_file_type =               0 /
beam_radius_max =             2.0 /
 bypass_sampler =               T /
calibrate_signal =               T /
    detector_id =       00_100_1a /
detpt_aberration =               T /
       dip_norm =             1.0 /
   dipole_speed =           TOTAL /
dipole_thermotemp =               T /
    dipole_type =               2 /
 first_pointing =             240 /
  focalplane_db = /scratch/cnt0028/ias1717/SHARED/RD12_data/dmc/MISS03/DATA/LevelS_mission_fits/fpdb_IMO_4_27.fits /
interpol_galactic =               T /
 interpol_order =              11 /
  last_pointing =            6000 /
       map_file = /scratch/cnt0028/ias1717/SHARED/bware/wrk_sroll21/LevelS/tmpdat/mmmap_ffp10_lensed_scl_cmb_100.fits /
nominal_pointing =               F /
          nside =             128 /
    output_type =          SIGNAL /
oversampling_factor =             1.0 /
        ringset = /scratch/cnt0028/ias1717/SHARED/bware/wrk_sroll21/LevelS/tmpdat/ringset_ffp10_lensed_scl_cmb_100.fits /
sampler_timeshift =             0.0 /
    satinfo_ctr = /scratch/cnt0028/ias1717/SHARED/RD12_data/dmc/MISS03/DATA/LevelS_mission_fits/ctr.fits /
satinfo_ephemeris = /scratch/cnt0028/ias1717/SHARED/RD12_data/dmc/MISS03/DATA/LevelS_mission_fits/ephemeris.fits /
satinfo_indexobject = /scratch/cnt0028/ias1717/SHARED/RD12_data/dmc/MISS03/DATA/LevelS_mission_fits/ringindex_begin_end.fits /
satinfo_quaternions = /scratch/cnt0028/ias1717/SHARED/RD12_data/dmc/MISS03/DATA/LevelS_mission_fits/quaternions_27100.fits /
   satinfo_type =             HFI /
single_precision_detpt =               T /
   solar_dipole =         HFI2018 /
  source_dipole =               F /
   source_fsldp =               F /
   source_mixed =               T /
  source_mixed2 =               F /
     source_oof =               F /
  source_pntsrc =               F /
       tod_file = /scratch/cnt0028/ias1717/SHARED/bware/wrk_sroll21/LevelS/tmpdat/mmtoi_ffp10_lensed_scl_cmb_100.fits /
variable_pntsrc_factor =             1.0 /
        verbose =               1 /

========== HDU[1] header ==========
       XTENSION =        BINTABLE / binary table extension
         BITPIX =               8 / 8-bit bytes
          NAXIS =               2 / 2-dimensional binary table
         NAXIS1 =               4 / width of table in bytes
         NAXIS2 =      2934314859 / number of rows in table
         PCOUNT =               0 / size of special data area
         GCOUNT =               1 / one data group (required keyword)
        TFIELDS =               1 / number of fields in each row
         TTYPE1 =          signal / label for field   1
         TFORM1 =              1E / data format of field: 4-byte REAL
         TUNIT1 =      K(Antenna) / physical unit of field
        EXTNAME =        xtension / name of this binary table extension
        objType =      toi.LS_toi /
firstPointingID =             240 /
 lastPointingID =            6000 /
   first_sample =      1368455701 /
start_time_mission =    1621174818.0 /
     start_time =   1628761591.45 /
         f_samp =       -180.3737 /
"""
