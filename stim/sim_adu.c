#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>   
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <math.h>
#include <fftw3.h> 

#include "sim_adu.h"
#include "stim_tools.h"
#include "no_dmc_data_access.h"
#include "no_dmc_metadata.h"

#include "binary_file_toi.h"


// debug: set to a valid rank number to save 4K lines value per sample in a TOI for this rank
#define SAVE_4KLINES (-1)


long getinl( double *inl, double x) {

  long v=0;
  int i;

  for (i=0; i<ADCBITS; i++) {
    v += 1<<(ADCBITS-1-i);
    if (inl[v]>x) v -= 1<<(ADCBITS-1-i);
  }
  return( v);
}


int sim_adu( stimParameters *stimPar, PIOFLOAT *signal) {

  stim_parContent *Param = &stimPar->Param;
  stim_Data       *data  = &stimPar->data;
  BRIMO           *brimo = &stimPar->brimo;

  long   i, j;
  int    sphase    = brimo->sphase;
  double KCMB2DSN  = stimPar->calib / (brimo->DSN2V * brimo->g0);

  int BeginRing  = stimPar->first_full_ring;
  int EndRing    = stimPar->last_full_ring;

  if (data->raw_begin_ring == 0) {
    // raw data not already read, will start at BeginRing
    data->raw_begin_ring = BeginRing;
    if (data->raw_begin_ring < 240) {
      data->raw_begin_ring = 240;
    }
  }
  int ring_count = EndRing - data->raw_begin_ring + 1;

  struct timeval t0;
  gettimeofday( &t0, NULL);

  STIM_TRACE( 1, "converting signal from KCMB to ADU (x%g), adding ADCNL and 4K lines for rings %d to %d", KCMB2DSN, BeginRing, EndRing);

  //////////////////////////////////////////////////////////////////////////////
  // read DNL: 65536 ADC steps (double)

  PIODOUBLE inl[INLSIZE];
  if (strcmp( Param->sim_inl_name, "0") != 0) {
    memcpy( inl, data->sim_inl, INLSIZE * sizeof( PIODOUBLE));

    if (strstr( Param->sim_inl_name, "_offsets") != NULL) {
      // <inl> contains offsets to add to the correction INL to produce the simulation INL
      PIODOUBLE offsets[INLSIZE];
      float sigma = inl[0];
      for (i = 0; i < INLSIZE; i++) offsets[i] = 0.0;

      if (strstr( Param->sim_inl_name, "SIMINL2803_offsets") != NULL) {
        STIM_TRACE( 1, "creating INL realisation (SIMINL2803) with sigma=%f", sigma);
        // stack cumulative offsets
        for (i = 30000; i < 35000; i++) {
          if (inl[i] != 0.0) {
            inl[i] += RANDOM( stimPar, 0) * sigma;
            STIM_TRACE( 2, "adding %f to corr_INL at pos=%ld", inl[i], i);
            for (j = -16; j < 16; j++) {
              // spread the offset to avoid decreasing INL
              for (int k = i+j; k < INLSIZE; k++) {
                offsets[k] += inl[i] / 32;
              }
            }
          }
        }
      } else if (strstr( Param->sim_inl_name, "SIMINL1406_offsets") != NULL) {
        STIM_TRACE( 1, "creating INL realisation (SIMINL1406) with sigma=%f", sigma);
        for (j = 0; j <= 15; j++) {
          inl[10+j] += RANDOM( stimPar, 0) * sigma;
        }
        for (i = 1; i < INLSIZE; i++) {
          for (j = 0; j <= 15; j++) {
            if (((i-1) & (1<<j)) != 0) { // i-1 because HFI INLs have 2^15 at 32769 instead of 32768
              offsets[i] += inl[10+j];
            }
          }
        }
      } else {
        STIM_TRACE( 0, "don't know what to do with %s", Param->sim_inl_name);
        exit( -1);
      }
      // add the offsets to the correction INL to produce the simulation INL
      memcpy( inl, data->corr_inl, INLSIZE * sizeof( PIODOUBLE));
      for (i = 0; i < INLSIZE; i++) {
        inl[i] += offsets[i];
      }
    }
  } else {
    // sim_inl_name==0, produce a linear INL
    for (i = 0; i < INLSIZE; i++) {
      inl[i] = i;
    }
  }

  // assert that INL is never decreasing
  for (i = 0; i < INLSIZE-1; i++) {
    if (inl[i+1] < inl[i]) {
      inl[i+1] = inl[i];
    }
  }
  
  if (Param->flag_save_sim_inl) {
    PIOSTRING save_name;
    snprintf( save_name, PIOSTRINGMAXLEN, "%s_%s_%03d", Param->save_sim_inl, Param->bolometer, Param->random_seed);
    assert( BFT_writeObject( save_name, "float64", 0, INLSIZE, inl) >= 0
            && "sim_adu.c: error while writing <save_sim_inl> to disk.");
  }

  //////////////////////////////////////////////////////////////////////////////
  // read raw constant: 80 fast sample values (double) per ring
  // for interval BeginRing-EndRing

  if (data->raw_cst == NULL) {
    data->raw_cst = malloc( 80 * ring_count * sizeof( PIODOUBLE));
    assert( data->raw_cst != NULL);
    if (strcmp( Param->rawcst_name, "0") != 0) {
      assert( noDMC_readObject_PIODOUBLE( Param->rawcst_name, 80 * data->raw_begin_ring, 80 * ring_count, data->raw_cst) >= 0);
    } else {
      for (i = 0; i < 80 * ring_count; i++) {
        data->raw_cst[i] = INLSIZE/2;
      }
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  // read 4K vect

  read_4k_vect( Param->fourk_name, data->raw_begin_ring, ring_count, data);

  //////////////////////////////////////////////////////////////////////////////
  // convert signal from KCMB to ADU

  for (i=0; i < stimPar->total_signal_length; i++) {
    signal[i] *= KCMB2DSN;
  }

  //////////////////////////////////////////////////////////////////////////////
  // optionally add per-bolometer DSN baseline

  if (strcmp( Param->add_baseline, "0") != 0) {
    if (data->dsn_baseline == NULL) {
      data->dsn_baseline = malloc( stimPar->total_signal_length * sizeof( PIOFLOAT));
      assert( data->dsn_baseline != NULL);
      assert( noDMC_readObject_PIOFLOAT( Param->add_baseline,
                                         stimPar->signal_first_sample_number,
                                         stimPar->total_signal_length,
                                         data->dsn_baseline) >= 0);
    }
    for (i=0; i < stimPar->total_signal_length; i++) {
      signal[i] += data->dsn_baseline[i];
      assert( !isnan( signal[i]));
    }
    if (!Param->stay_in_memory) {
      free( data->dsn_baseline);
      data->dsn_baseline = NULL;
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  // iterate on rings

  #pragma omp parallel for
  for (int ring = BeginRing; ring <= EndRing; ring++) {
    long begin_ring_index = BEGINRINGINDEX( ring);
    long ndata = ENDRINGINDEX( ring) - begin_ring_index + 1;
    int PARITY0 = PARITY( begin_ring_index);
    double ring_4k[80*9];
    int raw_ring = ring - data->raw_begin_ring; // ring number in raw data local arrays
    if ( raw_ring < 0) {
      raw_ring = 0;
    }

    STIM_TRACE( 2, "processing ring %d, n_harms=%d", ring, data->n_harmonics);

    PIOFLOAT *ring_signal = signal + begin_ring_index - stimPar->signal_first_sample_number;

    // prepare raw_cst for the ring
    PIODOUBLE ring_raw_cst[80];
    for (int i = 0; i < 80; i++) {
      ring_raw_cst[i] = data->raw_cst[80 * raw_ring + (i + sphase) % 80];
    }

    // compute 4K lines for the ring
    compute_4K_lines( raw_ring, begin_ring_index, data, ring_4k);

    for (int i = 0; i < ndata; i++) {
      double tmp = 0.0;
      if ((begin_ring_index + i) > stimPar->signal_first_sample_number + stimPar->total_signal_length) {
        // processing sample after signal array end, abort
        break;
      }
      if ((begin_ring_index + i) < stimPar->signal_first_sample_number) {
        // processing sample before signal array beginning
        // get 40 electronic white noise to keep RNG in sync, but discard result
        for (int j = 0; j < 40; j++) {
          whitenoise( stimPar, ring);
        }
        continue;
      }
      // modulate signal, add raw constant and 4K lines and ADC nonlinearity
      for (int j = 0; j < 40; j++) {
        tmp += getinl( inl, ring_raw_cst[PARITY( begin_ring_index + i) * 40 + j]
                          + ring_signal[i] * data->raw_gain[PARITY( begin_ring_index + i) * 40 + j]
                          + ring_4k[j + 40 * ((i+18 - PARITY0) % 18)]
                          + whitenoise( stimPar, ring) * Param->electronic_noise_adu); // 4 ADU electronic white noise
      }

      ring_signal[i] = tmp;
    }
  }

  
  if (!Param->stay_in_memory) {
    free( data->raw_cst);
    data->raw_cst = NULL;
    data->raw_begin_ring = 0;
    if (data->n_harmonics != 0) {
      free( data->fourk_harmonics);
      data->fourk_harmonics = NULL;
      free( data->fourk_amplitudes);
      data->fourk_amplitudes = NULL;
      data->n_harmonics = 0;
    }
  }

  STIM_TRACE( 1, "duration: %s", get_elapsed_time( &t0));

  return( 0);
}
