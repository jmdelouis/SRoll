#!/usr/bin/env python

################################################################################
#
# package_sroll.py
# script used to convert maps produced with sroll v3.0 (SVN r115) into
# PLA-like FITS files.
#
# This happens in two steps:
#   - first the individual sroll maps (SROLL_TPL) are read separately, processed,
#     and copied into a temporary location (TMPFITS_TPL)
#   - then the temporary processed files are put together in multi-columns FITS
#     files with PLA format (FITSOUT_TPL)
#
# To control the processing, some flags are used:
#   - do_hits:    set to 1 to copy hitcounts (and check them against PLA2018 when possible)
#   - do_cov:     set to 1 to invert precision matrix to produce covaraince matrix
#   - do_monodip: set to 1 to set monopole value and adjust removed solar dipole value
#   - do_fits:    set to 1 to build final FITS files in PLA format
#
# Remarks:
#   - do_cov is the most time consuming processing (about 10 minutes per map)
#   - do_hits, do_cov and do_monodip must have been done before doing do_fits
#   - you can set all 4 flags to 1 for a single run
#
# Other flags:
#   - do_monobolo:  set to 1 to process monobolometer maps
#   - do_multibolo: set to 1 to process mutlibolometer maps
#
# Other tools:
#   - do_check: set to 1 to check full frequency sroll packaged FITS maps against
#               corresponding PLA2018 files
#
################################################################################


import os
import sys
import time
import healpy
import numpy
import stools
from   IMO_4_27 import * # BOLOID, DETSETS, CALIB, NEP, XPOL, ELECWNOISE, PHOTWNOISE


################################################################################
# command line management

# list of HFI frequencies
HFIFREQS = ["100", "143", "217", "353", "545", "857"]

try:
  assert len( sys.argv) == 2
  if sys.argv[1] == "all":
    FREQS = HFIFREQS
  else:
    assert sys.argv[1] in HFIFREQS
    FREQS = [sys.argv[1]]
except:
  print "usage: %s <freq>" % sys.argv[0]
  print "  where <freq> is 100, 143, 217, 353, 545, 857 or 'all'"
  print "package sroll21 data for channel <freq>"
  print
  exit(-1)


################################################################################
# processing controlling flags

do_hits    = 0 # set to 1 to copy hitcounts (and check them against PLA2018 when possible)
do_cov     = 0 # set to 1 to invert precision matrix to produce covaraince matrix
do_monodip = 1 # set to 1 to set monopole value and adjust removed solar dipole value
do_fits    = 1 # set to 1 to build final FITS files in PLA format

do_monobolo  = 1 # set to 1 to process monobolometer maps
do_multibolo = 1 # set to 1 to process mutlibolometer maps

do_check = 0 # set to 1 to check sroll packaged maps against corresponding PLA2018 maps


################################################################################

def inv_cov( ii, iq, iu, qq, qu, uu):
  """return ii, iq, iu, qq, qu, uu of inverted matrix built with ii, iq, iu, qq, qu, uu
def invert_3x3( m):
  m1, m2, m3, m4, m5, m6, m7, m8, m9 = m.ravel()
  det = m1*m5*m9 + m4*m8*m3 + m7*m2*m6 - m1*m6*m8 - m3*m5*m7 - m2*m4*m9
  return numpy.array([m5*m9-m6*m8, m3*m8-m2*m9, m2*m6-m3*m5, m1*m9-m3*m7, m3*m4-m1*m6, m1*m5-m2*m4]) / det
"""
  det = ii*qq*uu + iq*qu*iu + iu*iq*qu - ii*qu*qu - iu*qq*iu - iq*iq*uu
  if det != 0:
    return numpy.array([qq*uu-qu*qu, iu*qu-iq*uu, iq*qu-iu*qq, ii*uu-iu*iu, iu*iq-ii*qu, ii*qq-iq*iq], dtype="float32") / det
  else:
    return numpy.zeros(6, dtype="float32")


################################################################################
# initialisations

# map size
NSIDE = 2048
NPIX  = 12*NSIDE*NSIDE

# list of extensions in polarised multibolo FITS files
POLEXTS   = ["I", "Q", "U", "H", "II", "IQ", "IU", "QQ", "QU", "UU"]
# list of extensions in non-polarised multibolo FITS files
NOPOLEXTS = ["I", "H", "II"]
# list of extensions in un-polarised monobolo FITS files
MONOEXTS  = ["C", "D", "H"]
# list of ring cuts, from sroll2/3 to PLA
PLACUTS   = {"full": "full", "hm1": "halfmission-1", "hm2": "halfmission-2", "fullodd":"full-oddring", "fulleven": "full-evenring"}
#PLACUTS   = {"full": "full"}

# directory containing PLA 2018 HFI FITS files, for checks
PLA18_DIR   = "/scratch/cnt0028/ias1717/SHARED/RD12_data/dmc/MISS03/DATA/PLA2018_HFI"
# root dir for intermediate files
WRKDIR      = "/scratch/cnt0028/ias1717/smottet/sroll21_datarelease"
# sroll output maps filename template
SROLL_TPL   = "/scratch/cnt0028/ias1717/smottet/s3_sroll21_MAP/sroll21_{run}_{dset}_{cut}_{comp}.fits"
SROLL_TPL   = "/scratch/cnt0028/ias1717/jmdelouis/s3_sroll21b_MAP/sroll21b_{run}_{dset}_{cut}_{comp}.fits"
SROLL_TPL   = "/scratch/cnt0028/ias1717/smottet/s3_sroll21_MAP/sroll211_{run}_{dset}_{cut}_{comp}.fits"
# temporay processed files name template
TMPFITS_TPL = WRKDIR + "/processed_maps/tmp_{PROCVER}_{run}_{dset}_{cut}_{comp}.fits"
# out FITS filename
PROCVER     = "SRoll22"
FITSOUT_TPL = WRKDIR + "/final_fits/{PROCVER}_SkyMap_{dset}_{placut}.fits"

# some extra FITS keywords (aladin?)
RESTRQ = { "100": 100.89,
           "143": 142.88,
           "217": 221.16,
           "353": 357.5,
           "545": 555.2,
           "857": 866.8,}
BNWID = {  "100": 33,
           "143": 46,
           "217": 65,
           "353": 102,
           "545": 171,
           "857": 245,}


if do_monodip:
  print( "computing dipole maps and reading masks")
  # Planck2015, HFI17 and HFI18 solar dipole maps
  PLK15DIP = stools.dipole2map( amplitude=3364.5*1e-6,  deg_lon=264.00,  deg_lat=48.253, nside=NSIDE)
  HFI17DIP = stools.dipole2map( amplitude=3362.71*1e-6, deg_lon=264.021, deg_lat=48.24,  nside=NSIDE)
  HFI18DIP = stools.dipole2map( amplitude=3362.08*1e-6, deg_lon=264.021, deg_lat=48.253, nside=NSIDE)

#  mask = healpy.read_map( "/scratch/cnt0028/ias1717/SHARED/RD12_data/dmc/MISS03/DATA/MASK_2048_GALACTIC_fits/MASK_S2PAPER_0_70_PS.fits")
#  assert len( mask) == NPIX

  galmask_fsky80 = healpy.read_map( PLA18_DIR + "/HFI_Mask_GalPlane-apo0_2048_R2.00.fits", field=4, verbose=False) # fields are fsky: 20, 40, 60, 70, 80, 90, 97, and 99%
  psmask_857 = healpy.read_map( PLA18_DIR + "/HFI_Mask_PointSrc_2048_R2.00.fits", field=5, verbose=False) # fields are: 100ghz, 143ghz, 217ghz, 353ghz, 545ghz, 857ghz
  mask = galmask_fsky80 * psmask_857
  print( "galmask*psmak, fsky=%.1f%%" % ((mask==1).sum()*100.0/NPIX))


################################################################################

def read_pla18_map( run, cut, comp):
  # PLA2018 only contains full frequency maps, but for all cuts
  freq = run[0:3]
  assert freq in HFIFREQS
  cut = PLACUTS[cut]
  if freq <= "353":
    field = POLEXTS.index( comp)
  else:
    field = NOPOLEXTS.index( comp)
  fitsname = PLA18_DIR + "/HFI_SkyMap_{freq}_2048_R3.01_{cut}.fits".format(**locals())
  return healpy.read_map( fitsname, field=field, verbose=False)


################################################################################

def read_sroll_map( run, dset, cut, comp):
  if comp == "H": comp = "HH"
  fitsname = SROLL_TPL.format(**locals())
  return healpy.read_map( fitsname, verbose=False)


################################################################################
# main processing loop

for freq in FREQS:
  run = freq+"ghz"
  if freq in ["100", "143", "217", "353"]:
    SIG_UNIT = "Kcmb"
    COV_UNIT = "Kcmb^2"
    is_pol   = True
    is_submm = False
  elif freq in ["545", "857"]:
    SIG_UNIT = "MJy/sr"
    COV_UNIT = "(MJy/sr)^2"
    is_pol   = False
    is_submm = True
  else:
    raise Exception( "Unknown frequency: " + freq)

  if do_monodip:
    p18i = read_pla18_map( run, "full", "I")
    tmpmap = p18i.copy()
    tmpmap[mask==0] = healpy.UNSEEN
    p18_monop = healpy.fit_monopole( tmpmap)
    print( "\nPLA2018 {freq}ghz masked monopole = {p18_monop}".format(**locals()))

  for dset in [freq+"ghz", freq+"psb", freq+"ds1", freq+"ds2"] + list( DETSETS[freq+"ghz"]):
    # skip invalid and unwanted detector sets
    if (dset == "100psb"): continue
    if (freq in ["545", "857"]) and (dset.endswith( "psb") or dset.endswith( "ds1") or dset.endswith( "ds2")): continue
    is_monobolo = ("-" in dset)
    if (is_monobolo) and (not do_monobolo): continue
    if (not is_monobolo) and (not do_multibolo): continue

    for cut in PLACUTS.keys():
      if ((cut != "full") and (not dset.endswith( "ghz"))): continue
      print( "\nProcessing %s %s" % (dset, cut))

      if do_monodip:
        if is_monobolo:
          if is_submm:
            COMP = ["I", "D"]
          else:
            COMP = ["I", "C", "D"]
        else:
          COMP = ["I"]
        for comp in COMP:
          mapi = read_sroll_map( run, dset, cut, comp)
          tmpmap = mapi.copy()
          tmpmap[mask==0] = healpy.UNSEEN
          sroll_monop = healpy.fit_monopole( tmpmap)
          print( "  monopole adjustment: %g" % (sroll_monop-p18_monop))
          # adjust monopole
          mapi = mapi - sroll_monop + p18_monop
          # adjust removed dipole
          mapi = mapi + HFI17DIP - HFI18DIP
          # write processed temperature map to temp dir
          print( "  adjusted dipole, writing " + comp)
          healpy.write_map( TMPFITS_TPL.format(**locals()), mapi, overwrite=True)
        if is_pol and (not is_monobolo):
          print "  copying Q and U"
          for comp in ["Q", "U"]:
            map = read_sroll_map( run, dset, cut, comp)
            healpy.write_map( TMPFITS_TPL.format( **locals()), map, overwrite=True)

      if do_hits:
        comp = "H"
        maph = numpy.around( read_sroll_map( run, dset, cut, comp)).astype( "int32")
        assert len( maph) == NPIX
        if (dset.endswith( "ghz")):
          p18h = read_pla18_map( run, cut, comp)
          diffmap = p18h-maph
          diff_hits = (diffmap != 0).sum()
          print( "  %d HITcount pixels differ from PLA2018 (%.3f%%)" % (diff_hits, 100.0*diff_hits/NPIX))
          diff_hits = numpy.abs( diffmap).sum()
          print( "  %d absolute different hits (%.3f%%) (PLA18-sroll)" % (diff_hits, 100.0*diff_hits/p18h.sum()))
          diff_hits = diffmap.sum()
          print( "  %d cumulated different hits (%.3f%%) (PLA18-sroll)" % (diff_hits, 100.0*diff_hits/p18h.sum()))
        else:
          print( "  HIT: no corresponding PLA2018 map")
        # write processed temperature map to temp dir
        print( "  writing H")
        healpy.write_map( TMPFITS_TPL.format(**locals()), maph, overwrite=True)

      if do_cov:
        print( "  reading precision matrix")
        if is_monobolo or (not is_pol):
          mapii = read_sroll_map( run, dset, cut, "II")
          covii = 1.0 / mapii
          covii[mapii == 0] = 0
          covii[covii < 0] = 0
          covii[numpy.logical_not( numpy.isfinite(covii))] = 0
          print( "  writing covariance matrix")
          comp="II"; healpy.write_map( TMPFITS_TPL.format(**locals()), covii, overwrite=True)

        else:
          mapi  = read_sroll_map( run, dset, cut, "I") # I is needed to identify UNSEEN pixels
          mapii = read_sroll_map( run, dset, cut, "II")
          mapiq = read_sroll_map( run, dset, cut, "IQ")
          mapiu = read_sroll_map( run, dset, cut, "IU")
          mapqq = read_sroll_map( run, dset, cut, "QQ")
          mapqu = read_sroll_map( run, dset, cut, "QU")
          mapuu = read_sroll_map( run, dset, cut, "UU")

          # invalid COV pixels are set to 0
          covii = numpy.zeros( NPIX, dtype="float32")
          coviq = numpy.zeros( NPIX, dtype="float32")
          coviu = numpy.zeros( NPIX, dtype="float32")
          covqq = numpy.zeros( NPIX, dtype="float32")
          covqu = numpy.zeros( NPIX, dtype="float32")
          covuu = numpy.zeros( NPIX, dtype="float32")

          print "  inverting precision matrix",
          sys.stdout.flush()
          negII = 0
          for i in range( NPIX):
            if (i % (NPIX/10)) == 0:
              print (int(round(float(i)/NPIX*10))),
              sys.stdout.flush()
            # UNSEEN pixels in I are set to 0 in II, IQ, IU, QQ, QU, UU
            if (mapi[i] != healpy.UNSEEN):
              covii[i], coviq[i], coviu[i], covqq[i], covqu[i], covuu[i] = inv_cov( mapii[i], mapiq[i], mapiu[i], mapqq[i], mapqu[i], mapuu[i])
              if (covii[i] < 0):
                covii[i] = 0
                coviq[i] = 0
                coviu[i] = 0
                covqq[i] = 0
                covqu[i] = 0
                covuu[i] = 0
                negII += 1

          print
          if (negII > 0):
            print "  %d negative values in II set to 0 in COV" % negII

          # write covariance matrix maps to temp dir
          print( "  writing covariance matrix")
          comp="II"; healpy.write_map( TMPFITS_TPL.format(**locals()), covii, overwrite=True)
          comp="IQ"; healpy.write_map( TMPFITS_TPL.format(**locals()), coviq, overwrite=True)
          comp="IU"; healpy.write_map( TMPFITS_TPL.format(**locals()), coviu, overwrite=True)
          comp="QQ"; healpy.write_map( TMPFITS_TPL.format(**locals()), covqq, overwrite=True)
          comp="QU"; healpy.write_map( TMPFITS_TPL.format(**locals()), covqu, overwrite=True)
          comp="UU"; healpy.write_map( TMPFITS_TPL.format(**locals()), covuu, overwrite=True)

      if do_fits:
        print( "  packaging FITS: reading all temporary maps")
        placut = PLACUTS[cut]
        fitsout = FITSOUT_TPL.format(**locals())

        # prepare FITS columns
        if is_monobolo:
          if freq <= "353":
            # 100-353 monobolo
            COMPLIST = ["C",        "D",         "I",       "H",    "II"]
            COLNAMES = ['I_STOKES', 'I_NOBPCOR', 'I_BPCOR', 'HITS', 'II_COV']
            COLUNITS = ['Kcmb',     'Kcmb',      'Kcmb',    ' ',    'Kcmb^2']
          else:
            # 545-857 monobolo
            COMPLIST = ["I",        "D",         "H",    "II"]
            COLNAMES = ['I_STOKES', 'I_NOBPCOR', 'HITS', 'II_COV']
            COLUNITS = ['MJy/sr',   'MJy/sr',    ' ',    '(MJy/sr)^2']
        else:
          if freq <= "353":
            # 100-353 multibolo
            COMPLIST = ["I",        "Q",        "U",        "H",    "II",     "IQ",     "IU",     "QQ",     "QU",     "UU"]
            COLNAMES = ['I_STOKES', 'Q_STOKES', 'U_STOKES', 'HITS', 'II_COV', 'IQ_COV', 'IU_COV', 'QQ_COV', 'QU_COV', 'UU_COV']
            COLUNITS = ['Kcmb',     'Kcmb',     'Kcmb',     ' ',    'Kcmb^2', 'Kcmb^2', 'Kcmb^2', 'Kcmb^2', 'Kcmb^2', 'Kcmb^2']
          else:
            # 545-857 multibolo
            COMPLIST = ["I",        "H",    "II"]
            COLNAMES = ['I_STOKES', 'HITS', 'II_COV']
            COLUNITS = ['MJy/sr',   ' ',    '(MJy/sr)^2']

        # read maps to put in FITS columns
        maps = list()
        dtypes = list()
        for comp in COMPLIST:
          tmpmap = healpy.read_map( TMPFITS_TPL.format(**locals()), verbose=False)
          if comp == "H":
            maps += [healpy.reorder( tmpmap.astype( "int32"), r2n=True)]
          else:
            maps += [healpy.reorder( tmpmap.astype( "float32"), r2n=True)]
          dtypes += [maps[-1].dtype]

        # extra header
        hdr = list()
        hdr += [('FILENAME', os.path.basename( fitsout), 'FITS filename')]
        hdr += [('DATE',     time.strftime("%Y-%m-%d"), 'Creation date')]
        hdr += [('EXTNAME',  'FREQ-MAP', 'Extension name')]
        hdr += [('POLCCONV', 'COSMO', 'Polarization convention')]
        hdr += [('COORDSYS', 'GALACTIC', 'Coordinate system')]
        hdr += [('BAD_DATA', '-1.6375e+30', 'HEALPIX bad pixel value')]
        hdr += [('FREQ',     str(freq), 'reference frequency')]
        hdr += [('PROCVER',  PROCVER, 'Product version')]
        hdr += [('UNITFREQ', "GHz", 'frequency units')]
        hdr += [('BNDCTR',   str(freq), 'band center, same as FREQ')]
        hdr += [('RESTFRQ',  RESTRQ[freq], 'effective frequency')]
        hdr += [('BNDWID',   BNWID[freq], 'effective bandwidth (approximate)')]
        hdr += [('PROCDSET', run, 'processed detector set')]
        hdr += [('PROCRSET', 'full', 'processed ring set')]
        hdr += [('PLOTDSET', dset, 'plotted detector set')]
        hdr += [('PLOTRSET', PLACUTS[cut], 'plotted ring set')]
        hdr += [('URL',      'http://sroll20.ias.u-psud.fr', 'SRoll2 web page')]
        print hdr
        print "  writing final output FITS file: " + fitsout
        healpy.write_map( fitsout, maps, nest=True, column_names=COLNAMES, column_units=COLUNITS, extra_header=hdr, dtype=dtypes, fits_IDL=False, overwrite=True)

#      exit(0) # exit after first cut
#    exit(0) # exit after all cuts of first detset
#  exit(0) # exit after all detset of first freq/run


################################################################################

if do_check:
  for freq in FREQS:
    for cut in sorted( PLACUTS.values()):
      if freq <= "353":
        COMPLIST = ["I", "Q", "U", "H", "II", "IQ", "IU", "QQ", "QU", "UU"]
      else:
        COMPLIST = ["I", "H", "II"]
      pname = PLA18_DIR + "/HFI_SkyMap_{freq}_2048_R3.01_{cut}.fits".format(**locals())
      sname = WRKDIR + "/final_fits/{PROCVER}_SkyMap_{freq}ghz_{cut}.fits".format(**locals())
      print "\n" + sname + "\n" + pname
      for field in range( len( COMPLIST)):
        pmap = healpy.read_map( pname, field=field, verbose=False)
        smap = healpy.read_map( sname, field=field, verbose=False)
        smap[pmap == healpy.UNSEEN] = healpy.UNSEEN
        smap[pmap == 0] = 0
        pmap[smap == healpy.UNSEEN] = healpy.UNSEEN
        pmap[smap == 0] = 0
        ndiff = numpy.sum( numpy.not_equal( pmap, smap))
        print "  %s diff=%.5f%%" % (COMPLIST[field], ndiff*100.0/len( smap))
