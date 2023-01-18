#define _GNU_SOURCE

#include <unistd.h>

#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <fcntl.h>
#include <assert.h>

#include "brimo.h"

char DetectorList[NUMBEROFBOLO][2][10] = {
  {"100-1a",  "00_100_1a"},
  {"100-1b",  "01_100_1b"},
  {"143-1a",  "02_143_1a"},
  {"143-1b",  "03_143_1b"},
  {"217-1",   "04_217_1"},
  {"353-1",   "05_353_1"},
  {"143-5",   "10_143_5"},
  {"217-5a",  "11_217_5a"},
  {"217-5b",  "12_217_5b"},
  {"353-2",   "13_353_2"},
  {"545-1",   "14_545_1"},
  {"Dark1",   "15_Dark1"},
  {"100-2a",  "20_100_2a"},
  {"100-2b",  "21_100_2b"},
  {"217-2",   "22_217_2"},
  {"353-3a",  "23_353_3a"},
  {"353-3b",  "24_353_3b"},
  {"857-1",   "25_857_1"},
  {"143-2a",  "30_143_2a"},
  {"143-2b",  "31_143_2b"},
  {"353-4a",  "32_353_4a"},
  {"353-4b",  "33_353_4b"},
  {"545-2",   "34_545_2"},
  {"857-2",   "35_857_2"},
  {"100-3a",  "40_100_3a"},
  {"100-3b",  "41_100_3b"},
  {"143-6",   "42_143_6"},
  {"217-6a",  "43_217_6a"},
  {"217-6b",  "44_217_6b"},
  {"353-7",   "45_353_7"},
  {"143-3a",  "50_143_3a"},
  {"143-3b",  "51_143_3b"},
  {"217-3",   "52_217_3"},
  {"353-5a",  "53_353_5a"},
  {"353-5b",  "54_353_5b"},
  {"545-3",   "55_545_3"},
  {"143-7",   "60_143_7"},
  {"217-7a",  "61_217_7a"},
  {"217-7b",  "62_217_7b"},
  {"353-6a",  "63_353_6a"},
  {"353-6b",  "64_353_6b"},
  {"857-3",   "65_857_3"},
  {"143-8",   "70_143_8"},
  {"217-8a",  "71_217_8a"},
  {"217-8b",  "72_217_8b"},
  {"545-4",   "73_545_4"},
  {"857-4",   "74_857_4"},
  {"Dark2",   "75_Dark2"},
  {"100-4a",  "80_100_4a"},
  {"100-4b",  "81_100_4b"},
  {"143-4a",  "82_143_4a"},
  {"143-4b",  "83_143_4b"},
  {"217-4",   "84_217_4"},
  {"353-8",   "85_353_8"}
};


// detector sets, see getDetset()
char *ds_100ds1[] = {"100-1a", "100-1b", "100-4a", "100-4b"};
char *ds_100ds2[] = {"100-2a", "100-2b", "100-3a", "100-3b"};
char *ds_100ghz[] = {"100-1a", "100-1b", "100-2a", "100-2b", "100-3a", "100-3b", "100-4a", "100-4b"};

char *ds_143ds1[] = {"143-1a", "143-1b", "143-3a", "143-3b"};
char *ds_143ds2[] = {"143-2a", "143-2b", "143-4a", "143-4b"};
char *ds_143swb[] = {"143-5",  "143-6",  "143-7"};
char *ds_143psb[] = {"143-1a", "143-1b", "143-2a", "143-2b", "143-3a", "143-3b", "143-4a", "143-4b"};
char *ds_143ghz[] = {"143-1a", "143-1b", "143-2a", "143-2b", "143-3a", "143-3b", "143-4a", "143-4b", "143-5",  "143-6",  "143-7"};

char *ds_217ds1[] = {"217-5a", "217-5b", "217-7a", "217-7b"};
char *ds_217ds2[] = {"217-6a", "217-6b", "217-8a", "217-8b"};
char *ds_217swb[] = {"217-1",  "217-2",  "217-3",  "217-4"};
char *ds_217psb[] = {"217-5a", "217-5b", "217-6a", "217-6b", "217-7a", "217-7b", "217-8a", "217-8b"};
char *ds_217ghz[] = {"217-1",  "217-2",  "217-3",  "217-4",  "217-5a", "217-5b", "217-6a", "217-6b", "217-7a", "217-7b", "217-8a", "217-8b"};

char *ds_353ds1[] = {"353-3a", "353-3b", "353-5a", "353-5b"};
char *ds_353ds2[] = {"353-4a", "353-4b", "353-6a", "353-6b"};
char *ds_353swb[] = {"353-1",  "353-2",  "353-7",  "353-8"};
char *ds_353psb[] = {"353-3a", "353-3b", "353-4a", "353-4b", "353-5a", "353-5b", "353-6a", "353-6b"};
char *ds_353ghz[] = {"353-1",  "353-2",  "353-3a", "353-3b", "353-4a", "353-4b", "353-5a", "353-5b", "353-6a", "353-6b", "353-7",  "353-8"};

char *ds_545ghz[] = {"545-1",  "545-2",            "545-4"};
char *ds_857ghz[] = {"857-1",  "857-2",  "857-3",  "857-4"};


char *PIXNAMES( int i) {
  return( DetectorList[i][0]);
}

char *BOLOIDS( int i) {
  return( DetectorList[i][1]);
}

// return the BoloID of a given detector name (pixname or boloid), exit() on failure
// !!!DON'T MODIFY THE RETURNED STRING, COPY IT FOR YOUR OWN NEEDS!!!

char *toBoloID( char *detname) {

  for (int i = 0; i < NUMBEROFBOLO; i++) {
    if (!strcmp( detname, PIXNAMES(i)) || !strcmp( detname, BOLOIDS(i))) {
      return( BOLOIDS(i));
    }
  }
  fprintf( stderr, "%s::%s() ERROR: '%s' is not a valid pixname/boloid\n", __FILE__, __FUNCTION__, detname);
  exit( -1);
}


int PixnameToBoloid(       PIOSTRING outBoloid,
                     const PIOSTRING inPixname)
{
  int i;

  memset( outBoloid, 0, PIOSTRINGMAXLEN);

  for (i = 0; i < NUMBEROFBOLO; i++) {
    if (!strcmp( inPixname, PIXNAMES(i))) {
      strcpy(outBoloid, BOLOIDS(i));
      return( IDX2BC( i));
    }
  }

  /* tests if inPixName is already a boloid, and if so returns it */
  for (i = 0; i < NUMBEROFBOLO; i++) {
    if (!strcmp( inPixname, BOLOIDS(i))) {
      strcpy( outBoloid, BOLOIDS(i));
      return( IDX2BC( i));
    }
  }

  fprintf( stderr, "PixnameToBoloid ERROR: '%s' (strlen=%ld) is not a valid pixname/boloid\n", inPixname, strlen(inPixname));
  exit( -1);
}


/*==========================================================================
 PIOBoloIDToPixName : returns the PixName of a given BoloID                #
  -------------------------------------------------------------------------#
  Date            Author            Comment                                #
  18-jan-2008     mottet@iap.fr     initial version                        #
  24-jan-2008        "              returns inBoloID if inBoloID           #
                                    is already a PixName                   #
===========================================================================*/

int BoloidToPixname(       PIOSTRING outPixname,
                     const PIOSTRING inBoloid)
{
  int i;

  memset( outPixname, 0, PIOSTRINGMAXLEN);

  for (i = 0; i < NUMBEROFBOLO; i++) {
    if (!strcmp(inBoloid, BOLOIDS(i))) {
      strcpy(outPixname, PIXNAMES(i));
      return( IDX2BC( i));
    }
  }

  /* tests if inBoloid is already a pixname, and if so returns it */
  for (i = 0; i < NUMBEROFBOLO; i++) {
    if (!strcmp(inBoloid, PIXNAMES(i))) {
      strcpy(outPixname, PIXNAMES(i));
      return( IDX2BC( i));
    }
  }

  fprintf( stderr, "BoloidToPixname ERROR: '%s' (strlen=%ld) is not a valid pixname/boloid\n", inBoloid, strlen(inBoloid));
  exit( -1);
}


/*******************************************************************************
getCopsb(): fill the ouptut char array with the pixname's co-psb
*******************************************************************************/

void getCoPSB( char *out_copsb, char *in_pixname)
{
  int l = strlen( in_pixname);
  strncpy( out_copsb, in_pixname, l);
  if (out_copsb[l-1] == 'a') {
    out_copsb[l-1] = 'b';
  }
  else if (out_copsb[l-1] == 'b') {
    out_copsb[l-1] = 'a';
  }
}


/*******************************************************************************
getBC(): return an int equal to the BC number of bolometer
*******************************************************************************/

int getBC( const char *bolometer)
{
  PIOSTRING blablabla;
  return( PixnameToBoloid( blablabla, bolometer));
}


/*******************************************************************************
GetFrequency(): return an int equal to the bolometer's frequency
*******************************************************************************/

int GetFrequency( const char* bolometer) {

  PIOSTRING bolo_id;
  PIOSTRING bolo_freq;

  PixnameToBoloid(bolo_id, bolometer);
  sprintf( bolo_freq, "%.3s", bolo_id + 3);
  return atoi( bolo_freq);
}


/*******************************************************************************
isDark(): return 1 is bolometer is Dark, else 0
*******************************************************************************/

int isDark( int bc)
{
  // 15_Dark1, 75_Dark2
  if ((bc == 15) || (bc == 75)) {
    return( 1);
  }
  return( 0);
}


/*******************************************************************************
isPopcorned(): return 1 is bolometer is popcorned, else 0
*******************************************************************************/

int isPopcorned( int bc)
{
  // 55_545_3, 70_143_8
  if ((bc == 55) || (bc == 70)) {
    return( 1);
  }
  return( 0);
}


/*******************************************************************************
BC2IDX(): transform a BC into the bolometer's index in DetectorList[] (and in BRIMO file)
*******************************************************************************/

int BC2IDX( int bc)
{
  return( bc / 10 * 6 + bc % 10);
}


/*******************************************************************************
IDX2BC(): inverse of BC2IDX...
*******************************************************************************/

int IDX2BC( int idx)
{
  return( idx / 6 * 10 + idx % 6);
}


/*******************************************************************************
getDetset(): return the list of bolometers inside a detset
*******************************************************************************/

void getDetset( char *detset_name, char ***bolo_list, int *n_bolo) {

  int i;
  char dslow[7]; // lower case of detset name
  for (i = 0; i<7; i++){
    dslow[i] = tolower( detset_name[i]);
  }

  if (strncmp( dslow, "100ds1", 6) == 0) { *bolo_list = ds_100ds1; *n_bolo =  4; return;}
  if (strncmp( dslow, "100ds2", 6) == 0) { *bolo_list = ds_100ds2; *n_bolo =  4; return;}
  if (strncmp( dslow, "100psb", 6) == 0) { *bolo_list = ds_100ghz; *n_bolo =  8; return;}
  if (strncmp( dslow, "100ghz", 6) == 0) { *bolo_list = ds_100ghz; *n_bolo =  8; return;}
  if (strncmp( dslow, "143ds1", 6) == 0) { *bolo_list = ds_143ds1; *n_bolo =  4; return;}
  if (strncmp( dslow, "143ds2", 6) == 0) { *bolo_list = ds_143ds2; *n_bolo =  4; return;}
  if (strncmp( dslow, "143swb", 6) == 0) { *bolo_list = ds_143swb; *n_bolo =  3; return;}
  if (strncmp( dslow, "143psb", 6) == 0) { *bolo_list = ds_143psb; *n_bolo =  8; return;}
  if (strncmp( dslow, "143ghz", 6) == 0) { *bolo_list = ds_143ghz; *n_bolo = 11; return;}
  if (strncmp( dslow, "217ds1", 6) == 0) { *bolo_list = ds_217ds1; *n_bolo =  4; return;}
  if (strncmp( dslow, "217ds2", 6) == 0) { *bolo_list = ds_217ds2; *n_bolo =  4; return;}
  if (strncmp( dslow, "217swb", 6) == 0) { *bolo_list = ds_217swb; *n_bolo =  4; return;}
  if (strncmp( dslow, "217psb", 6) == 0) { *bolo_list = ds_217psb; *n_bolo =  8; return;}
  if (strncmp( dslow, "217ghz", 6) == 0) { *bolo_list = ds_217ghz; *n_bolo = 12; return;}
  if (strncmp( dslow, "353ds1", 6) == 0) { *bolo_list = ds_353ds1; *n_bolo =  4; return;}
  if (strncmp( dslow, "353ds2", 6) == 0) { *bolo_list = ds_353ds2; *n_bolo =  4; return;}
  if (strncmp( dslow, "353swb", 6) == 0) { *bolo_list = ds_353swb; *n_bolo =  4; return;}
  if (strncmp( dslow, "353psb", 6) == 0) { *bolo_list = ds_353psb; *n_bolo =  8; return;}
  if (strncmp( dslow, "353ghz", 6) == 0) { *bolo_list = ds_353ghz; *n_bolo = 12; return;}
  if (strncmp( dslow, "545ghz", 6) == 0) { *bolo_list = ds_545ghz; *n_bolo =  3; return;}
  if (strncmp( dslow, "857ghz", 6) == 0) { *bolo_list = ds_857ghz; *n_bolo =  4; return;}

  // detset not found, try single bolometer...
  for (i=0; i<NUMBEROFBOLO; i++) {
    if (strcmp( dslow, PIXNAMES(i)) == 0) { *bolo_list = (char**) DetectorList[i]; *n_bolo = 1; return;}
  }

  //really really not found
  fprintf( stderr, "%s/%s(): unknown detset name: %s\n", __FILE__, __FUNCTION__, detset_name);
  exit( 1);
}

/*******************************************************************************
read a BRIMO structure from a binary file, for bolo bc
*******************************************************************************/

int readBRIMO( char *brimo_filename, int bc, BRIMO *brimo) {

  char bfn[260];

  if (brimo_filename == NULL) {
    char *srollhost = NULL;
    srollhost = getenv("SROLLHOST");
    if (srollhost == NULL) {
      fprintf( stderr, "brimo.c::readBRIMO error: SROLLHOST environment variable not set.\n");
      exit(1);
    }
    if (strcmp( srollhost, "M3") == 0) {
      strcpy( bfn, "/wrk/symottet/brimo_4_27");
    }
    else if (strcmp( srollhost, "M4") == 0) {
      strcpy( bfn, "/wrk/symottet/brimo_4_27");
    }
    else if (strcmp( srollhost, "EDISON") == 0) {
      strcpy( bfn, "/project/projectdirs/planck/data/hfi/RD12_data/brimo/brimo_4_27");
    }
    else {
      fprintf( stderr, "brimo.c::readBRIMO error: unkonwon SROLLHOST (%s).\n", srollhost);
      exit(1);
    }
  }
  else {
    strcpy( bfn, brimo_filename);
  }

  // need only read right to read the BRIMO
  int brimo_filedesc = open( bfn, O_RDONLY);
  assert( brimo_filedesc >= 0);

  // check BRIMO file contains up-to-date BRIMO structures
  struct stat brimo_stat;
  memset( &brimo_stat, 0, sizeof( brimo_stat));
  assert( fstat( brimo_filedesc, &brimo_stat) != -1);
  assert( brimo_stat.st_size == NUMBEROFBOLO * sizeof( BRIMO));

  // read BRIMO structure associated with BC
  assert( pread( brimo_filedesc, brimo, sizeof( BRIMO), BC2IDX( bc) * sizeof( BRIMO)) == sizeof( BRIMO));
  assert( close( brimo_filedesc) == 0);

  // check we've read what we think...
  assert( brimo->bc == bc);
  return 0;
}


/*******************************************************************************
write a bolometers's BRIMO structure to a binary file
*******************************************************************************/

int writeBRIMO( char *brimo_filename, BRIMO *brimo) {

  int brimo_filedesc = open( brimo_filename, O_CREAT | O_RDWR, S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH);
  assert( brimo_filedesc >= 0);

  assert( pwrite( brimo_filedesc, brimo, sizeof( BRIMO), BC2IDX( brimo->bc) * sizeof( BRIMO)) == sizeof( BRIMO));

  assert( close( brimo_filedesc) == 0);
  return 0;
}

/*******************************************************************************
print a BRIMO structure to screen
*******************************************************************************/

void printBRIMO( BRIMO *brimo) {
  fprintf( stderr, "\n");
  fprintf( stderr, "printBRIMO( %s):\n",     brimo->pixname);
  fprintf( stderr, "  boloid         = %s\n", brimo->boloid);
  fprintf( stderr, "  bc             = %d\n", brimo->bc);
  fprintf( stderr, "  freq           = %d\n", brimo->freq);
  fprintf( stderr, "  KCMB2WATT      = %g\n", brimo->KCMB2WATT);
  fprintf( stderr, "  DSN2V          = %g\n", brimo->DSN2V);
  fprintf( stderr, "  g0             = %g\n", brimo->g0);
  fprintf( stderr, "  v0             = %g\n", brimo->v0);
  fprintf( stderr, "  polar_angle    = %g\n", brimo->polar_angle_rad);
  fprintf( stderr, "  polar_leakage  = %g\n", brimo->polar_leakage);
  fprintf( stderr, "  NET            = %g\n", brimo->NET);
  fprintf( stderr, "  DX11ADC_GPver  = %d\n", brimo->DX11ADC_GPver);
  fprintf( stderr, "  sphase         = %d\n", brimo->sphase);
  fprintf( stderr, "  compstep       = %d\n", brimo->compstep);
  fprintf( stderr, "  LFER_A1        = %g\n", brimo->LFER_A1);
  fprintf( stderr, "  LFER_tau1      = %g\n", brimo->LFER_tau1);
  fprintf( stderr, "  LFER_A2        = %g\n", brimo->LFER_A2);
  fprintf( stderr, "  LFER_tau2      = %g\n", brimo->LFER_tau2);
  fprintf( stderr, "  LFER_A3        = %g\n", brimo->LFER_A3);
  fprintf( stderr, "  LFER_tau3      = %g\n", brimo->LFER_tau3);
  fprintf( stderr, "  LFER_A4        = %g\n", brimo->LFER_A4);
  fprintf( stderr, "  LFER_tau4      = %g\n", brimo->LFER_tau4);
  fprintf( stderr, "  LFER_A5        = %g\n", brimo->LFER_A5);
  fprintf( stderr, "  LFER_tau5      = %g\n", brimo->LFER_tau5);
  fprintf( stderr, "  LFER_A6        = %g\n", brimo->LFER_A6);
  fprintf( stderr, "  LFER_tau6      = %g\n", brimo->LFER_tau6);
  fprintf( stderr, "  LFER_A7        = %g\n", brimo->LFER_A7);
  fprintf( stderr, "  LFER_tau7      = %g\n", brimo->LFER_tau7);
  fprintf( stderr, "  LFER_A8        = %g\n", brimo->LFER_A8);
  fprintf( stderr, "  LFER_tau8      = %g\n", brimo->LFER_tau8);
  fprintf( stderr, "  LFER_tau_stray = %g\n", brimo->LFER_tau_stray);
  fprintf( stderr, "  LFER_sphase    = %g\n", brimo->LFER_sphase);
  fprintf( stderr, "  global_offset  = %g\n", brimo->LFER_global_offset);
  fprintf( stderr, "  LFER_tauhp     = %g\n", brimo->LFER_tauhp);
  fprintf( stderr, "  lpf_fgauss     = %g\n", brimo->lpf_fgauss);
  fprintf( stderr, "  lpf_fc         = %g\n", brimo->lpf_fc);
  fprintf( stderr, "  lpf_ffact      = %g\n", brimo->lpf_ffact);
  fprintf( stderr, "  photon_noise   = %g\n", brimo->photonic_whitenoise);
  fprintf( stderr, "  electro_noise  = %g\n", brimo->electronic_whitenoise);
  fprintf( stderr, "  oof_slope      = %g\n", brimo->oof_slope);
  fprintf( stderr, "  oof_fknee      = %g\n", brimo->oof_fknee);
  fprintf( stderr, "\n");

}

void printBRIMOtotable( BRIMO *brimo) {

  if (brimo == NULL) {
    printf( "pixname\tboloid\tbc\tfreq\tKCMB2WATT\tDSN2V\t"
            "g0\tv0\tpolar_angle\tpolar_leakage\tNET\t"
            "DX11ADC_GPver\tsphase\tcompstep\t"
            "LFER_A1\tLFER_tau1\tLFER_A2\tLFER_tau2\t"
            "LFER_A3\tLFER_tau3\tLFER_A4\tLFER_tau4\t"
            "LFER_A5\tLFER_tau5\tLFER_A6\tLFER_tau6\t"
            "LFER_A7\tLFER_tau7\tLFER_A8\tLFER_tau8\t"
            "LFER_tau_stray\tLFER_sphase\tLFER_global_offset\tLFER_tauhp\t"
            "lpf_fgauss\tlpf_fc\tlpf_ffact\t"
            "photonic_whitenoise\telectronic_whitenoise\toof_slope\toof_fknee\n");
  } else {
    printf( "%s\t%s\t%02d\t%03d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%d\t%d\t%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
    brimo->pixname, brimo->boloid, brimo->bc, brimo->freq, brimo->KCMB2WATT, brimo->DSN2V,
    brimo->g0, brimo->v0, brimo->polar_angle_rad, brimo->polar_leakage, brimo->NET,
    brimo->DX11ADC_GPver, brimo->sphase, brimo->compstep,
    brimo->LFER_A1, brimo->LFER_tau1, brimo->LFER_A2, brimo->LFER_tau2,
    brimo->LFER_A3, brimo->LFER_tau3, brimo->LFER_A4, brimo->LFER_tau4,
    brimo->LFER_A5, brimo->LFER_tau5, brimo->LFER_A6, brimo->LFER_tau6,
    brimo->LFER_A7, brimo->LFER_tau7, brimo->LFER_A8, brimo->LFER_tau8,
    brimo->LFER_tau_stray, brimo->LFER_sphase, brimo->LFER_global_offset,  brimo->LFER_tauhp,
    brimo->lpf_fgauss, brimo->lpf_fc, brimo->lpf_ffact,
    brimo->photonic_whitenoise, brimo->electronic_whitenoise, brimo->oof_slope, brimo->oof_fknee);
  }
}


////////////////////////////////////////////////////////////////////////////////
// conversion factor from KCMB to MJyPerSr for 545GHz-857GHz

/*
import LSCtools
imo_id = LSCtools.Get_IMO_4_27_ID()
for bolo in sorted( [b for b in LSCtools.BOLOID.keys() if b.startswith( "545") or b.startswith( "857")]):
  print( '  if (!strcmp( pixname, "%s")) return( %g);' % ( bolo, LSCtools.GetKCMB2MJyPerSr( bolo, imo_id)))
*/

float get_KCMB2MJyPerSr( char* pixname) {
  if (!strcmp( pixname, "545-1")) return( 57.0831);
  if (!strcmp( pixname, "545-2")) return( 58.8824);
  if (!strcmp( pixname, "545-3")) return( 57.8793);
  if (!strcmp( pixname, "545-4")) return( 58.0595);
  if (!strcmp( pixname, "857-1")) return( 2.18904);
  if (!strcmp( pixname, "857-2")) return( 2.34559);
  if (!strcmp( pixname, "857-3")) return( 2.21322);
  if (!strcmp( pixname, "857-4")) return( 2.40212);
  return( 1.0);
}


double get_RD12calib( char *pixname) {
  if (!strcmp( pixname, "100-1a")) return( 1.00340513024e-13);
  if (!strcmp( pixname, "100-1b")) return( 1.23240311873e-13);
  if (!strcmp( pixname, "100-2a")) return( 1.51966997388e-13);
  if (!strcmp( pixname, "100-2b")) return( 1.5885889296e-13);
  if (!strcmp( pixname, "100-3a")) return( 1.38652193689e-13);
  if (!strcmp( pixname, "100-3b")) return( 1.1492356927e-13);
  if (!strcmp( pixname, "100-4a")) return( 1.46826316998e-13);
  if (!strcmp( pixname, "100-4b")) return( 1.17469574177e-13);
  if (!strcmp( pixname, "143-1a")) return( 1.88008633662e-13);
  if (!strcmp( pixname, "143-1b")) return( 1.61678010691e-13);
  if (!strcmp( pixname, "143-2a")) return( 1.78260067e-13);
  if (!strcmp( pixname, "143-2b")) return( 1.82131500475e-13);
  if (!strcmp( pixname, "143-3a")) return( 1.79339046635e-13);
  if (!strcmp( pixname, "143-3b")) return( 1.61409996829e-13);
  if (!strcmp( pixname, "143-4a")) return( 1.66035966523e-13);
  if (!strcmp( pixname, "143-4b")) return( 1.55447052831e-13);
  if (!strcmp( pixname, "143-5"))  return( 2.65238565653e-13);
  if (!strcmp( pixname, "143-6"))  return( 2.381597318e-13);
  if (!strcmp( pixname, "143-7"))  return( 2.57719266261e-13);
  if (!strcmp( pixname, "217-1"))  return( 1.68882666765e-13);
  if (!strcmp( pixname, "217-2"))  return( 1.62797073142e-13);
  if (!strcmp( pixname, "217-3"))  return( 1.71515434481e-13);
  if (!strcmp( pixname, "217-4"))  return( 1.66584754679e-13);
  if (!strcmp( pixname, "217-5a")) return( 1.14058850217e-13);
  if (!strcmp( pixname, "217-5b")) return( 1.15245126451e-13);
  if (!strcmp( pixname, "217-6a")) return( 1.15733750737e-13);
  if (!strcmp( pixname, "217-6b")) return( 1.17133951737e-13);
  if (!strcmp( pixname, "217-7a")) return( 1.24037931807e-13);
  if (!strcmp( pixname, "217-7b")) return( 1.19572741647e-13);
  if (!strcmp( pixname, "217-8a")) return( 1.20172971416e-13);
  if (!strcmp( pixname, "217-8b")) return( 1.1391920766e-13);
  if (!strcmp( pixname, "353-1"))  return( 5.81896764062e-14);
  if (!strcmp( pixname, "353-2"))  return( 6.31788126648e-14);
  if (!strcmp( pixname, "353-3a")) return( 3.505451069e-14);
  if (!strcmp( pixname, "353-3b")) return( 3.46195375785e-14);
  if (!strcmp( pixname, "353-4a")) return( 2.89764285963e-14);
  if (!strcmp( pixname, "353-4b")) return( 2.89508885946e-14);
  if (!strcmp( pixname, "353-5a")) return( 3.31816651442e-14);
  if (!strcmp( pixname, "353-5b")) return( 3.30720800648e-14);
  if (!strcmp( pixname, "353-6a")) return( 2.38419917429e-14);
  if (!strcmp( pixname, "353-6b")) return( 2.16742504299e-14);
  if (!strcmp( pixname, "353-7"))  return( 4.79439720697e-14);
  if (!strcmp( pixname, "353-8"))  return( 4.50246808415e-14);
  if (!strcmp( pixname, "545-1"))  return( 3.29932298836e-16);
  if (!strcmp( pixname, "545-2"))  return( 3.08792864906e-16);
  if (!strcmp( pixname, "545-4"))  return( 2.64855405172e-16);
  if (!strcmp( pixname, "857-1"))  return( 3.30076826046e-16);
  if (!strcmp( pixname, "857-2"))  return( 3.55811287601e-16);
  if (!strcmp( pixname, "857-3"))  return( 3.18681631353e-16);
  if (!strcmp( pixname, "857-4"))  return( 2.219187708e-16);
  fprintf( stderr, "get_RD12calib(): Unkown pixname (%s)", pixname);
  exit( 1);
}
