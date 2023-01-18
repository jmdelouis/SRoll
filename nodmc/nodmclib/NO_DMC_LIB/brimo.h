#ifndef _BRIMO_H_
#define _BRIMO_H_


#include "no_dmc_piolib_type_def.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>


// bolometer list tools

#define DETNAMELEN 10
typedef char DETNAME[DETNAMELEN];

#define NUMBEROFBOLO (54)

char *PIXNAMES( int i);

char *BOLOIDS( int i);

char *toBoloID( char *inPixname);

int PixnameToBoloid(       PIOSTRING outBoloid,
                     const PIOSTRING inPixname);

int BoloidToPixname(       PIOSTRING outPixname,
                     const PIOSTRING inBoloid);

void getCoPSB( char *out_copsb, char *in_pixname);

int getBC( const char *bolometer);

int GetFrequency( const char* bolometer);

int isDark( int bc);

int isPopcorned( int bc);

void getDetset( char *detset_name, char ***bolo_list, int *n_bolo);


// some global HFI IMO parameters

#define SAMPLING_FREQ (180.3737)


// brimo: Binary Reduced Instrument MOdel

typedef struct {
  char  pixname[10];
  char  boloid[10];
  int   bc;
  int   freq;
  float KCMB2WATT; // gain
  float DSN2V;
  float g0; // g0 and v0 are used to convert from Volt to Watt with: Watt = g0 * Volt * (1.0 + Volt / v0)
  float v0;
  float coeffT90;
  float dg_dT_PAU;
  float dg_dT_REU;
  float polar_angle_rad;
  float polar_leakage;
  float NET;
  int   DX11ADC_GPver;
  int   sphase;
  int   compstep; // compression step
  // time constant convolution / deconvolution parameters
  double LFER_A1;
  double LFER_A2;
  double LFER_A3;
  double LFER_A4;
  double LFER_tau1;
  double LFER_tau2;
  double LFER_tau3;
  double LFER_tau4;
  double LFER_tau_stray;
  double LFER_sphase;
  double LFER_global_offset;
  double LFER_A5;
  double LFER_A6;
  double LFER_tau5;
  double LFER_tau6;
  double LFER_A7;
  double LFER_A8;
  double LFER_tau7;
  double LFER_tau8;
  double LFER_tauhp;
  double fwidth_rfilter;
  // lowpass filter parameters
  double lpf_fgauss; // fgauss_lpf (1sigma)
  double lpf_fc;     // fc_lpf (cosine cutoff freq.)
  double lpf_ffact;  // ffact_lpf (cosine factor)
  // noise power spectrum parameters
  double photonic_whitenoise;
  double electronic_whitenoise;
  double oof_slope;
  double oof_fknee;
} BRIMO;

// BC2IDX(): transform a BC into the bolometer's index in DetectorList[] (and in BRIMO file)
int BC2IDX( int bc);

// IDX2BC(): inverse of BC2IDX...
int IDX2BC( int idx);

// read a BRIMO structure from a binary file, for bolo bc
int readBRIMO( char *brimo_filename, int bc, BRIMO *brimo);

// write a bolometer's BRIMO structure to a binary file,
int writeBRIMO( char *brimo_filename, BRIMO *brimo);

// print a BRIMO structure to screen
void printBRIMO( BRIMO *brimo);

void printBRIMOtotable( BRIMO *brimo);

// conversion factor from KCMB to MJyPerSr for 545GHz-857GHz
float get_KCMB2MJyPerSr( char* pixname);

// conversion factor from KCMB to Watts used in RD12 productions
double get_RD12calib( char *pixname);

#endif
