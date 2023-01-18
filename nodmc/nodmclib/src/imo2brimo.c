
#include <assert.h>

#include "HL2_PIOLIB/PIOLib.h"

#ifdef ISDMC
#include "HL2_DMC/PIODB.h"
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "PioLib/HL2_PIOLIB/PIOErr.h"

#include "brimo.h"

/*
 */


double GetIMODouble( char* bolo_id, char* imo_path, PIOIMO *imo_ptr) {

  PIOErr    pioerr;
  PIOSTRING full_path;
  IMOValue  imo_val;

  sprintf( full_path, "IMO:HFI:DET:Phot_Pixel Name='%s'%s", bolo_id, imo_path);
  pioerr = PIOGetValue( &imo_val, full_path, "PIODOUBLE", imo_ptr);
  if (pioerr != 0) {
    fprintf( stderr, "Error in PIOGetValue(%s)\n", full_path);
    exit( -1);
  }

  return( *(PIODOUBLE *)(imo_val.data));
}


float GetPolarAngleInRadians( char* bolo_id, PIOIMO *imo_ptr) {

  PIOErr    pioerr;
  PIOSTRING imo_path;
  IMOQuat   tempQuat;

  sprintf( imo_path, "IMO:HFI:DET:Phot_Pixel Name='%s':OpticalProperties:BeamI:Ellipse:Quaternion", bolo_id);
  pioerr = PIOGetQuaternion(&tempQuat, imo_path, imo_ptr);
  if (pioerr != 0) {
    fprintf( stderr, "Error in PIOGetQuaternion(%s)\n", imo_path);
    exit( -1);
  }

  return (2.0 * atan2f(tempQuat.Q[3], tempQuat.Q[0]));
}


double GetDSN2V( char* bolo_id, PIOIMO *imo_ptr) {

  // http://prof.planck.fr/IMG/pdf/HFI_TransferFunctionv6.2.1.pdf
  // http://wiki.planck.fr/index.php/DataDictionary/ElectronicChainHFIHSK

  PIOErr    pioerr;
  PIOSTRING imo_path;
  IMOValue  imo_val;
  PIODOUBLE F1_bc;
  PIODOUBLE GC_REU_of_Gamp;

  sprintf( imo_path, "IMO:HFI:REU:HFI_REU_bc bc_hexa='%.2s':F1", bolo_id);
  pioerr = PIOGetValue( &imo_val, imo_path, "PIODOUBLE", imo_ptr);
  if (pioerr != 0) {
    fprintf( stderr, "Error in PIOGetValue(%s)\n", imo_path);
    exit( -1);
  }
  assert( !strcmp( imo_val.Unit, "ADU/V"));
  F1_bc = *(PIODOUBLE*)(imo_val.data);

  if (GetFrequency( bolo_id) >= 545) {
    sprintf( imo_path, "IMO:HFI:REU:HFI_REU_bc bc_hexa='%.2s':GainC:PGain0", bolo_id);
    pioerr = PIOGetValue( &imo_val, imo_path, "PIODOUBLE", imo_ptr);
    if (pioerr != 0) {
      fprintf( stderr, "Error in PIOGetValue(%s)\n", imo_path);
      exit( -1);
    }
    GC_REU_of_Gamp = *(PIODOUBLE*)(imo_val.data);
  } else {
    GC_REU_of_Gamp = 1.0;
  }

  return( 1.0 / GC_REU_of_Gamp / F1_bc / 40.0);
}


double GetLFERParam( char* bolo_id, char* lfer, char* par_name, PIOIMO *imo_ptr) {

  PIOSTRING lfer_path;

  sprintf( lfer_path, ":NoiseAndSyst:TimeResp:%s:%s", lfer, par_name);
  return( GetIMODouble( bolo_id, lfer_path, imo_ptr));
}


double GetKCMB2WATT( char* bolo_id, PIOIMO *imo_ptr) {

  PIOErr    pioerr;
  PIOSTRING MajorVersion;
  PIOSTRING MinorVersion;
  PIOSTRING imo_leaf;
  PIOSTRING imo_path;
  IMOValue  imo_val;
  PIODOUBLE calib_factor_to_watt;

  pioerr = PIOInfoIMO( MajorVersion, MinorVersion, imo_ptr);
  if ((atoi( MajorVersion) > 4) || ((atoi( MajorVersion) == 4) && (atoi( MinorVersion) >= 16))) {
    sprintf( imo_path, "IMO:HFI:DET:Phot_Pixel Name='%s':RadiativeResp:Default", bolo_id);
    pioerr = PIOGetValue( &imo_val, imo_path, "PIOSTRING", imo_ptr);
    if (pioerr != 0) {
      fprintf( stderr, "Error in PIOGetValue(%s)\n", imo_path);
      exit( -1);
    }
    sprintf( imo_leaf, "%s", imo_val.data);
  } else {
    if (GetFrequency( bolo_id) <= 353) {
      sprintf( imo_leaf, "%s", "Dipole");
    } else {
      sprintf( imo_leaf, "%s", "Galactic");
    }
  }

  sprintf( imo_path, "IMO:HFI:DET:Phot_Pixel Name='%s':RadiativeResp:%s", bolo_id, imo_leaf);
  pioerr = PIOGetValue( &imo_val, imo_path, "PIODOUBLE", imo_ptr);
  if (pioerr != 0) {
    fprintf( stderr, "Error in PIOGetValue(%s)\n", imo_path);
    exit( -1);
  }

  if (strstr( imo_val.Unit, "W/K_CMB")) {
    calib_factor_to_watt = 1.0;
  } else if (strstr( imo_val.Unit, "W/MJy/sr")) {
//    calib_factor_to_watt = GetKCMB2MJyPerSr( bolo_id, imo_ptr);
    calib_factor_to_watt = GetIMODouble( bolo_id, ":SpectralResp:SpecTransmissions:KCMB_to_MJy_per_sr_IRAS", imo_ptr);
  } else {
    fprintf( stderr, "Unkown calibration factor unit in %s\n", imo_path);
    exit( -1);
  }

  return( *(PIODOUBLE*)(imo_val.data) * calib_factor_to_watt);
}

int GetCompStep( char* pixname) {
  // from http://cvs.planck.fr/cvs/Level2/Pipe_pkg/HL2_LSCorePipe/src/LSCtools.py
  if (!strcmp( pixname, "100-1a")) return( 38);
  if (!strcmp( pixname, "100-1b")) return( 38);
  if (!strcmp( pixname, "100-2a")) return( 37);
  if (!strcmp( pixname, "100-2b")) return( 38);
  if (!strcmp( pixname, "100-3a")) return( 39);
  if (!strcmp( pixname, "100-3b")) return( 44);
  if (!strcmp( pixname, "100-4a")) return( 36);
  if (!strcmp( pixname, "100-4b")) return( 38);
  if (!strcmp( pixname, "143-1a")) return( 40);
  if (!strcmp( pixname, "143-1b")) return( 44);
  if (!strcmp( pixname, "143-2a")) return( 42);
  if (!strcmp( pixname, "143-2b")) return( 55);
  if (!strcmp( pixname, "143-3a")) return( 47);
  if (!strcmp( pixname, "143-3b")) return( 45);
  if (!strcmp( pixname, "143-4a")) return( 42);
  if (!strcmp( pixname, "143-4b")) return( 48);
  if (!strcmp( pixname, "143-5" )) return( 45);
  if (!strcmp( pixname, "143-6" )) return( 43);
  if (!strcmp( pixname, "143-7" )) return( 45);
  if (!strcmp( pixname, "217-1" )) return( 43);
  if (!strcmp( pixname, "217-2" )) return( 42);
  if (!strcmp( pixname, "217-3" )) return( 43);
  if (!strcmp( pixname, "217-4" )) return( 46);
  if (!strcmp( pixname, "217-5a")) return( 40);
  if (!strcmp( pixname, "217-5b")) return( 50);
  if (!strcmp( pixname, "217-6a")) return( 43);
  if (!strcmp( pixname, "217-6b")) return( 45);
  if (!strcmp( pixname, "217-7a")) return( 50);
  if (!strcmp( pixname, "217-7b")) return( 49);
  if (!strcmp( pixname, "217-8a")) return( 45);
  if (!strcmp( pixname, "217-8b")) return( 47);
  if (!strcmp( pixname, "353-1" )) return( 55);
  if (!strcmp( pixname, "353-2" )) return( 56);
  if (!strcmp( pixname, "353-3a")) return( 47);
  if (!strcmp( pixname, "353-3b")) return( 44);
  if (!strcmp( pixname, "353-4a")) return( 50);
  if (!strcmp( pixname, "353-4b")) return( 43);
  if (!strcmp( pixname, "353-5a")) return( 43);
  if (!strcmp( pixname, "353-5b")) return( 43);
  if (!strcmp( pixname, "353-6a")) return( 45);
  if (!strcmp( pixname, "353-6b")) return( 46);
  if (!strcmp( pixname, "353-7" )) return( 46);
  if (!strcmp( pixname, "353-8" )) return( 48);
  if (!strcmp( pixname, "545-1" )) return( 15);
  if (!strcmp( pixname, "545-2" )) return( 14);
  if (!strcmp( pixname, "545-4" )) return( 15);
  if (!strcmp( pixname, "857-1" )) return( 26);
  if (!strcmp( pixname, "857-2" )) return( 26);
  if (!strcmp( pixname, "857-3" )) return( 31);
  if (!strcmp( pixname, "857-4" )) return( 24);
  return( 0);
}


int GetSPhase( char* pixname) {
  // from http://cvs.planck.fr/cvs/Level2/Pipe_pkg/HL2_LSCorePipe/src/LSCtools.py
  if (!strcmp( pixname, "100-1a")) return( 10);
  if (!strcmp( pixname, "100-1b")) return( 10);
  if (!strcmp( pixname, "100-2a")) return(  9);
  if (!strcmp( pixname, "100-2b")) return(  9);
  if (!strcmp( pixname, "100-3a")) return(  9);
  if (!strcmp( pixname, "100-3b")) return(  9);
  if (!strcmp( pixname, "100-4a")) return(  9);
  if (!strcmp( pixname, "100-4b")) return( 10);
  if (!strcmp( pixname, "143-1a")) return(  9);
  if (!strcmp( pixname, "143-1b")) return(  9);
  if (!strcmp( pixname, "143-2a")) return(  9);
  if (!strcmp( pixname, "143-2b")) return(  9);
  if (!strcmp( pixname, "143-3a")) return(  9);
  if (!strcmp( pixname, "143-3b")) return(  6);
  if (!strcmp( pixname, "143-4a")) return(  9);
  if (!strcmp( pixname, "143-4b")) return(  9);
  if (!strcmp( pixname, "143-5" )) return( 10);
  if (!strcmp( pixname, "143-6" )) return(  8);
  if (!strcmp( pixname, "143-7" )) return( 10);
  if (!strcmp( pixname, "217-1" )) return(  8);
  if (!strcmp( pixname, "217-2" )) return(  9);
  if (!strcmp( pixname, "217-3" )) return(  9);
  if (!strcmp( pixname, "217-4" )) return(  8);
  if (!strcmp( pixname, "217-5a")) return(  8);
  if (!strcmp( pixname, "217-5b")) return(  9);
  if (!strcmp( pixname, "217-6a")) return(  9);
  if (!strcmp( pixname, "217-6b")) return(  8);
  if (!strcmp( pixname, "217-7a")) return( 10);
  if (!strcmp( pixname, "217-7b")) return( 10);
  if (!strcmp( pixname, "217-8a")) return(  8);
  if (!strcmp( pixname, "217-8b")) return(  9);
  if (!strcmp( pixname, "353-1" )) return(  7);
  if (!strcmp( pixname, "353-2" )) return(  7);
  if (!strcmp( pixname, "353-3a")) return(  9);
  if (!strcmp( pixname, "353-3b")) return(  8);
  if (!strcmp( pixname, "353-4a")) return(  9);
  if (!strcmp( pixname, "353-4b")) return(  8);
  if (!strcmp( pixname, "353-5a")) return(  8);
  if (!strcmp( pixname, "353-5b")) return(  8);
  if (!strcmp( pixname, "353-6a")) return(  9);
  if (!strcmp( pixname, "353-6b")) return(  8);
  if (!strcmp( pixname, "353-7" )) return(  9);
  if (!strcmp( pixname, "353-8" )) return(  8);
  if (!strcmp( pixname, "545-1" )) return(  8);
  if (!strcmp( pixname, "545-2" )) return(  7);
  if (!strcmp( pixname, "545-4" )) return(  8);
  if (!strcmp( pixname, "857-1" )) return(  8);
  if (!strcmp( pixname, "857-2" )) return(  9);
  if (!strcmp( pixname, "857-3" )) return(  8);
  if (!strcmp( pixname, "857-4" )) return(  6);
  return( 0);
}

int GetDX11ADC_GPver( int bc) {
// from http://cvs.planck.fr/cvs/Level2/Task_pkg/HL2_LSMC/src/correctadc_versions.py
  if (bc == 00) return( 20);
  if (bc == 01) return( 20);
  if (bc == 02) return( 21);
  if (bc == 03) return( 20);
  if (bc == 04) return( 20);
  if (bc == 05) return( 21);
  if (bc == 10) return( 20);
  if (bc == 11) return( 21);
  if (bc == 12) return( 12);
  if (bc == 13) return( 20);
  if (bc == 14) return( 20);
  if (bc == 15) return( 00); // Dark1
  if (bc == 20) return( 21);
  if (bc == 21) return( 20);
  if (bc == 22) return( 20);
  if (bc == 23) return( 20);
  if (bc == 24) return( 20);
  if (bc == 25) return( 21);
  if (bc == 30) return( 21);
  if (bc == 31) return( 20);
  if (bc == 32) return( 20);
  if (bc == 33) return( 21);
  if (bc == 34) return( 20);
  if (bc == 35) return( 20);
  if (bc == 40) return( 20);
  if (bc == 41) return( 20);
  if (bc == 42) return( 20);
  if (bc == 43) return( 20);
  if (bc == 44) return( 20);
  if (bc == 45) return( 21);
  if (bc == 50) return( 21);
  if (bc == 51) return( 21);
  if (bc == 52) return( 20);
  if (bc == 53) return( 21);
  if (bc == 54) return( 20);
  if (bc == 55) return( 20); // 545-3, popcorned
  if (bc == 60) return( 21);
  if (bc == 61) return( 21);
  if (bc == 62) return( 20);
  if (bc == 63) return( 21);
  if (bc == 64) return( 20);
  if (bc == 65) return( 21);
  if (bc == 70) return( 00); // 143-8, popcorned
  if (bc == 71) return( 20);
  if (bc == 72) return( 21);
  if (bc == 73) return( 20);
  if (bc == 74) return( 20);
  if (bc == 75) return( 00); // Dark2
  if (bc == 80) return( 20);
  if (bc == 81) return( 20);
  if (bc == 82) return( 21);
  if (bc == 83) return( 20);
  if (bc == 84) return( 21);
  if (bc == 85) return( 21);
  return( 0);
}


int main( int argc, char *argv[]) {

  int i;
  char *endptr;
  PIOSTRING imo_id;
  PIOIMO *imo_ptr = NULL;
  PIOSTRING brimo_filename;
  PIOSTRING lfer;
  BRIMO brimo;


  // process command line arguments
  if (argc == 1) {
    brimo_filename[0] = 0; // don't write BRIMO, just print IMO contents in a table
    printf( "usage: imo2brimo <output_filename>\n");
    printf( "       if <output_filename> is omitted, brimo will be printed to screen but not saved to a file.\n\n");
    printf( "!! BRIMO NOT WRITTEN TO FILE !!\n\n");
  } else {
//    snprintf( brimo_filename, PIOSTRINGMAXLEN, "/wrk/symottet/srollex/brimo_4_27_test");
    snprintf( brimo_filename, PIOSTRINGMAXLEN, argv[1]);
    printf( "Writing BRIMO to %s\n", brimo_filename);
  }

  // open IMO
  snprintf( imo_id, PIOSTRINGMAXLEN, "/data/dmc/MISS03/METADATA%%1731634056892899949"); // IMO 4.27
  imo_ptr = PIOOpenIMOFile( imo_id, "r");
  if (imo_ptr == NULL) {
    fprintf( stderr, "Error in PIOOpenIMOFile(%s)\n", imo_id);
    exit( -1);
  }

  printf( "sizeof(BRIMO)=%ld\n", (long)sizeof(BRIMO));
  printBRIMOtotable( NULL);
  for (i = 0; i < NUMBEROFBOLO; i++) {
//    fprintf( stdout, "reading IMO for %s\n", BOLOIDS(i));
    memset( &brimo, 0, sizeof( brimo));
    brimo.bc = IDX2BC( i);
    snprintf( brimo.pixname, 10, PIXNAMES(i));
    snprintf( brimo.boloid,  10, BOLOIDS(i));
    if (isDark( brimo.bc) || isPopcorned( brimo.bc)) {
      // dark and popcorned bolos
      brimo.freq = 0;
      brimo.KCMB2WATT = 0.0;
    } else {
      brimo.freq = strtol( brimo.pixname, &endptr, 10);
      assert( endptr != brimo.pixname);
      brimo.KCMB2WATT = GetKCMB2WATT( brimo.boloid, imo_ptr);
    }
    brimo.DSN2V = GetDSN2V( brimo.boloid, imo_ptr);
    // g0 and v0 are used to convert from Volt to Watt with: Watt = g0 * Volt * (1.0 + Volt / v0)
    brimo.g0 =       GetIMODouble( brimo.boloid, ":NoiseAndSyst:NonLinearity:g0", imo_ptr);
    brimo.v0 =       GetIMODouble( brimo.boloid, ":NoiseAndSyst:NonLinearity:v0", imo_ptr);
    brimo.coeffT90 = GetIMODouble( brimo.boloid, ":NoiseAndSyst:NoiseTechData:CorrelTemp:CorrTemp Name='CorrTemp_1'", imo_ptr);

    brimo.polar_angle_rad = GetPolarAngleInRadians( brimo.boloid, imo_ptr);
    brimo.polar_leakage   = GetIMODouble( brimo.boloid, ":OpticalProperties:PolarProperty:PolarLeakage", imo_ptr);

    if (brimo.KCMB2WATT != 0.0) {
      // see http://wiki.planck.fr/index.php/Proc/DetectorWeightingInMaps#toc1 for conversion between NET and NEP
      brimo.NET = GetIMODouble( brimo.boloid, ":NoiseAndSyst:FlightNoise:TotalNoise:NEP", imo_ptr);
      brimo.NET = brimo.NET / brimo.KCMB2WATT * sqrt( SAMPLING_FREQ / 2);
    } else {
      brimo.NET = 0.0;
    }

    brimo.DX11ADC_GPver = GetDX11ADC_GPver( brimo.bc);
    brimo.sphase =        GetSPhase( brimo.pixname);
    brimo.compstep =      GetCompStep( brimo.pixname);
    if (brimo.freq == 0) {
      lfer == NULL;
    } else if (brimo.freq >= 545) {
      sprintf( lfer, "LFER4");
    } else if (brimo.freq == 353) {
      sprintf( lfer, "LFER8");
    } else {
      sprintf( lfer, "LFER6");
    }
    if (brimo.freq != 0) { // LFER4
      brimo.LFER_A1            = GetLFERParam( brimo.boloid, lfer, "par1", imo_ptr);
      brimo.LFER_A2            = GetLFERParam( brimo.boloid, lfer, "par2", imo_ptr);
      brimo.LFER_A3            = GetLFERParam( brimo.boloid, lfer, "par3", imo_ptr);
      brimo.LFER_A4            = GetLFERParam( brimo.boloid, lfer, "par9", imo_ptr);
      brimo.LFER_tau1          = GetLFERParam( brimo.boloid, lfer, "par4", imo_ptr);
      brimo.LFER_tau2          = GetLFERParam( brimo.boloid, lfer, "par5", imo_ptr);
      brimo.LFER_tau3          = GetLFERParam( brimo.boloid, lfer, "par6", imo_ptr);
      brimo.LFER_tau4          = GetLFERParam( brimo.boloid, lfer, "par10", imo_ptr);
      brimo.LFER_tau_stray     = GetLFERParam( brimo.boloid, lfer, "par7", imo_ptr);
      brimo.LFER_sphase        = GetLFERParam( brimo.boloid, lfer, "par8", imo_ptr);
      brimo.LFER_global_offset = GetLFERParam( brimo.boloid, lfer, "global_offset", imo_ptr);
      if (brimo.freq <= 353) { // LFER6
        brimo.LFER_A5   = GetLFERParam( brimo.boloid, lfer, "par11", imo_ptr);
        brimo.LFER_A6   = GetLFERParam( brimo.boloid, lfer, "par13", imo_ptr);
        brimo.LFER_tau5 = GetLFERParam( brimo.boloid, lfer, "par12", imo_ptr);
        brimo.LFER_tau6 = GetLFERParam( brimo.boloid, lfer, "par14", imo_ptr);
      } else {
        brimo.LFER_A5    = 0.0;
        brimo.LFER_A6    = 0.0;
        brimo.LFER_tau5  = 0.0;
        brimo.LFER_tau6  = 0.0;
      }
      if (brimo.freq == 353) { // LFER8
        brimo.LFER_A7    = GetLFERParam( brimo.boloid, lfer, "par15", imo_ptr);
        brimo.LFER_A8    = GetLFERParam( brimo.boloid, lfer, "par17", imo_ptr);
        brimo.LFER_tau7  = GetLFERParam( brimo.boloid, lfer, "par16", imo_ptr);
        brimo.LFER_tau8  = GetLFERParam( brimo.boloid, lfer, "par18", imo_ptr);
        brimo.LFER_tauhp = GetLFERParam( brimo.boloid, "SallenKeyHPF", "tauhp1", imo_ptr);
      } else {
        brimo.LFER_A7    = 0.0;
        brimo.LFER_A8    = 0.0;
        brimo.LFER_tau7  = 0.0;
        brimo.LFER_tau8  = 0.0;
        brimo.LFER_tauhp = 0.051;
      }
    } else {
      brimo.LFER_A1    = 0.0;
      brimo.LFER_A2    = 0.0;
      brimo.LFER_A3    = 0.0;
      brimo.LFER_A4    = 0.0;
      brimo.LFER_tau1  = 0.0;
      brimo.LFER_tau2  = 0.0;
      brimo.LFER_tau3  = 0.0;
      brimo.LFER_tau4  = 0.0;
      brimo.LFER_tau_stray     = 0.0;
      brimo.LFER_sphase        = 0.0;
      brimo.LFER_global_offset = 0.0;
    }

    brimo.lpf_fgauss = GetIMODouble( brimo.boloid, ":NoiseAndSyst:TimeResp:LowPassFilter:par1", imo_ptr);
    brimo.lpf_fc     = GetIMODouble( brimo.boloid, ":NoiseAndSyst:TimeResp:LowPassFilter:par2", imo_ptr);
    brimo.lpf_ffact  = GetIMODouble( brimo.boloid, ":NoiseAndSyst:TimeResp:LowPassFilter:par3", imo_ptr);

    // noise spectrum keywords
    // if noise spectrum DMC VECT doesn't exist, values will be 0.0
    PIOSTRING specname, comment;
    double *photonoise, gain;
    brimo.oof_slope = 0.0;
    brimo.oof_fknee = 0.0;
    brimo.electronic_whitenoise = 0.0;
    brimo.photonic_whitenoise   = 0.0;
    sprintf( specname, "/data/dmc/MISS03/DATA/detnoise_simu/%s_WhiteNoisePlusOneOverF_Watts_Spectrum_Jun15", brimo.pixname);
    if (PIOCheckObject( specname, NULL) == 0) {
      assert( PIOReadKeywordObject( (void *) &gain,            comment, (char *)"GAIN",  (char *)"PIODOUBLE", specname, NULL ) == 0);
      assert( PIOReadKeywordObject( (void *) &brimo.oof_slope, comment, (char *)"SLOPE", (char *)"PIODOUBLE", specname, NULL ) == 0);
      assert( PIOReadKeywordObject( (void *) &brimo.oof_fknee, comment, (char *)"FKNEE", (char *)"PIODOUBLE", specname, NULL ) == 0);
      assert( PIOReadKeywordObject( (void *) &brimo.electronic_whitenoise, comment, (char *)"elec_whitenoise", (char *)"PIODOUBLE", specname, NULL ) == 0);
      assert( PIOReadVECTObject( (void **) &photonoise, specname, "PIODOUBLE", "", NULL) > 0);
      brimo.photonic_whitenoise = photonoise[0] * gain;
    }
    else {
      // added on 9 Dec. 2016 for 545-857GHz, for reference
      sprintf( specname, "/data/dmc/MISS03/DATA/detnoise_simu/%s_MeanSpectrumFit_deconv_byRing_v64", brimo.pixname);
      if (PIOCheckObject( specname, NULL) == 0) {
        assert( PIOReadKeywordObject( (void *) &gain,            comment, (char *)"GAIN",  (char *)"PIODOUBLE", specname, NULL ) == 0);
        assert( PIOReadKeywordObject( (void *) &brimo.oof_slope, comment, (char *)"SLOPE", (char *)"PIODOUBLE", specname, NULL ) == 0);
        assert( PIOReadKeywordObject( (void *) &brimo.oof_fknee, comment, (char *)"FKNEE", (char *)"PIODOUBLE", specname, NULL ) == 0);
        assert( PIOReadVECTObject( (void **) &photonoise, specname, "PIODOUBLE", "", NULL) > 0);
        brimo.photonic_whitenoise   = photonoise[0] * gain * get_KCMB2MJyPerSr( brimo.pixname);
        brimo.electronic_whitenoise = brimo.photonic_whitenoise;
      }
    }

    printBRIMOtotable( &brimo);
    if (brimo_filename[0] != 0) {
      assert( writeBRIMO( brimo_filename, &brimo) == 0);
    }
  }
  
  PIOCloseIMO( &imo_ptr);

  exit( 0);
}
