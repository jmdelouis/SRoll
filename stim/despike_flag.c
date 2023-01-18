
#define _XOPEN_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>   
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <math.h>

#include "no_dmc_data_access.h"

#include "despike_flag.h"
#include "stim_tools.h"

#include "binary_file_toi.h"
#include "ring_index.h"

#define PBRSIZE (10822l)
#define HPRSIZE (27664l)

float GetThreshold( char* pixname);

int despike_flag( stimParameters *stimPar, PIOFLOAT *signal) {

  stim_parContent *Param = &stimPar->Param;
  stim_Data       *data  = &stimPar->data;
  BRIMO           *brimo = &stimPar->brimo;

  long   i;
  int    ring;
  double pbr_signal[PBRSIZE];
  int    pbr_hitcount[PBRSIZE];
  double mean, sigma;
  float  noise_threshold = Param->do_despike_flag;
  double sigma_mean = 0.0;
  int    sigma_count = 0;

  if (Param->do_despike_flag != 1.0) {
    noise_threshold = Param->do_despike_flag; // custom value
  } else {
    noise_threshold = GetThreshold( brimo->pixname); // default value
  }

  int BeginRing = stimPar->first_full_ring;
  int EndRing   = stimPar->last_full_ring;

  STIM_TRACE( 1, "despike flagging with threshold = %g for rings %d-%d", noise_threshold, BeginRing, EndRing);

  // read Phase_dx11 object
  if (data->phase == NULL) {
    data->phase = malloc( stimPar->total_signal_length * sizeof( PIODOUBLE));
    assert(( data->phase != NULL) && "despike_flag: Not enough memory to allocate phase");
    assert( noDMC_readObject_PIODOUBLE( Param->phase,
                                        stimPar->signal_first_sample_number,
                                        stimPar->total_signal_length,
                                        data->phase) >= 0);
  }

  // loop on rings
  for (ring = BeginRing; ring <= EndRing; ring++) {
    long first_sample = BEGINRINGINDEX( ring) - stimPar->signal_first_sample_number;
    long ring_length  = ENDRINGINDEX( ring) - BEGINRINGINDEX( ring) + 1;
    int  bin; // temporary bin number in PBR

    double good_sum[HPRSIZE]; // sum of noise values below threshold
    long   good_count[HPRSIZE];
    double bad_sum[HPRSIZE];  // sum of noise values above threshold
    long   bad_count[HPRSIZE];
    double pixel_value[HPRSIZE]; // mean of "good" samples in the same pixel

    float  *ring_signal = signal            + first_sample;
    double *ring_phase  = data->phase       + first_sample;
    PIOINT *ring_hpridx = data->pixel_index + first_sample; // hpridx[i] == -1 when sample i is flagged

    // roll signal+noise to PBR without flag
    TOI2PBR( ring_signal,
             NULL,
             ring_phase,
             ring_length,
             PBRSIZE,
             pbr_signal,
             pbr_hitcount);

    // produce noise-only TOI by subtracting signal PBR from signal+noise TOI
    PIOFLOAT *noisetoi = malloc( ring_length * sizeof( PIOFLOAT));
    assert( noisetoi != NULL);
    for (i = 0; i < ring_length; i++) {
      bin = (int)(ring_phase[i] * PBRSIZE / (2 * M_PI));
      assert( pbr_hitcount[bin] != 0);
      noisetoi[i] = ring_signal[i] - pbr_signal[bin];
    }

    // get ring noise-only TOI stddev
    stddev( noisetoi, ring_length, &mean, &sigma);
    if (sigma != NAN) {
      sigma_mean = (sigma_mean * sigma_count + sigma) / (sigma_count + 1);
      sigma_count++;
    }

    // build the sum and count of noise values above and below the threshold, per HPR pixel
    for (i = 0; i < HPRSIZE; i++) {
      good_sum[i]    = 0.0;
      good_count[i]  = 0;
      bad_sum[i]     = 0.0;
      bad_count[i]   = 0;
      pixel_value[i] = 0.0;
    }
    long n_over = 0; // number of samples over noise threshold value, for debug
    for (i = 0; i < ring_length; i++) {
      if (ring_hpridx[i] >= 0) {
        if (!isnan( noisetoi[i])) {
          if (fabs( noisetoi[i]) > noise_threshold * sigma) {
            bad_sum[ring_hpridx[i]]   += noisetoi[i];
            bad_count[ring_hpridx[i]] += 1;
            n_over += 1;
          } else {
            good_sum[ring_hpridx[i]]    += noisetoi[i];
            pixel_value[ring_hpridx[i]] += ring_signal[i];
            good_count[ring_hpridx[i]]  += 1;
          }
        }
      }
    }

    // compute each pixel value (=mean of unflagged samples with noise below threshold)
    for (i = 0; i < HPRSIZE; i++) {
      if (good_count[i] != 0) {
        pixel_value[i] /= good_count[i];
      }
    }

    // replace samples which noise is above threshold with pixel value
    for (i = 0; i < ring_length; i++) {
      if (ring_hpridx[i] != -1) {
        if (fabs( noisetoi[i]) > noise_threshold * sigma) {
          if (good_count[ring_hpridx[i]] > 0) {
            ring_signal[i] = pixel_value[ring_hpridx[i]];
          }
        }
      }
    }

#if 0
    STIM_TRACE( 1, "ring=%d, mean=%g, sigma=%g, bad samples=%ld/%ld (%.2f%%)", ring, mean, sigma, n_over, ring_length, ((float)n_over/ring_length)*100);

    assert( BFT_writeObject( "/pscratch1/RD12_data/dmc/MISS03/DATA/sylvain2_VEC/despike_flag_testnodau_pbr_signal",
                             "float32", ring*PBRSIZE, PBRSIZE, pbr_signal) >= 0);

    assert( BFT_writeObject( "/pscratch1/RD12_data/dmc/MISS03/DATA/sylvain2_VEC/despike_flag_testnoadu_noisetoi",
                             "float32", BEGINRINGINDEX(ring), ring_length, noisetoi) >= 0);

    assert( BFT_writeObject( "/pscratch1/RD12_data/dmc/MISS03/DATA/sylvain2_VEC/despike_flag_testnoadu_outsignal",
                             "float32", BEGINRINGINDEX(ring), ring_length, ring_signal) >= 0);
#endif

    free( noisetoi);
  }

  if (!Param->stay_in_memory) {
    free( data->phase);
    data->phase = NULL;
  }

  STIM_TRACE( 1, "complete, mean noise standard deviation = %g", sigma_mean);

  return( 0);
}


////////////////////////////////////////////////////////////////////////////////

float GetThreshold( char* pixname) {
  // 100ghz,143psb,217psb default values measured by Luca Pagano on 29 May 2016 plus 143swb and 217swb measured on 22 July 2016
  if (!strcmp( pixname, "100-1a")) return( 2.98520130823);
  if (!strcmp( pixname, "100-1b")) return( 2.90248985749);
  if (!strcmp( pixname, "100-2a")) return( 2.93047045829);
  if (!strcmp( pixname, "100-2b")) return( 2.93255428637);
  if (!strcmp( pixname, "100-3a")) return( 2.99172835377);
  if (!strcmp( pixname, "100-3b")) return( 2.84079613444);
  if (!strcmp( pixname, "100-4a")) return( 3.08909941657);
  if (!strcmp( pixname, "100-4b")) return( 3.03979246286);
  if (!strcmp( pixname, "143-1a")) return( 3.03015823145);
  if (!strcmp( pixname, "143-1b")) return( 2.79999145914);
  if (!strcmp( pixname, "143-2a")) return( 2.74958023541);
  if (!strcmp( pixname, "143-2b")) return( 2.96763702737);
  if (!strcmp( pixname, "143-3a")) return( 2.83686321703);
  if (!strcmp( pixname, "143-3b")) return( 2.76409703608);
  if (!strcmp( pixname, "143-4a")) return( 2.7733293191);
  if (!strcmp( pixname, "143-4b")) return( 2.89880820027);
  if (!strcmp( pixname, "143-5" )) return( 2.89698908163);
  if (!strcmp( pixname, "143-6" )) return( 2.83217666713);
  if (!strcmp( pixname, "143-7" )) return( 2.82990696068);
  if (!strcmp( pixname, "217-1" )) return( 2.90299160001);
  if (!strcmp( pixname, "217-2" )) return( 2.97047732833);
  if (!strcmp( pixname, "217-3" )) return( 2.9556098268);
  if (!strcmp( pixname, "217-4" )) return( 2.91198556997);
  if (!strcmp( pixname, "217-5a")) return( 2.93068203776);
  if (!strcmp( pixname, "217-5b")) return( 2.81248631621);
  if (!strcmp( pixname, "217-6a")) return( 2.96102235022);
  if (!strcmp( pixname, "217-6b")) return( 2.83576829032);
  if (!strcmp( pixname, "217-7a")) return( 3.09584893285);
  if (!strcmp( pixname, "217-7b")) return( 2.78212581432);
  if (!strcmp( pixname, "217-8a")) return( 2.95380669502);
  if (!strcmp( pixname, "217-8b")) return( 2.79703576323);
  if (!strcmp( pixname, "353-1" )) return( 2.88024393761);
  if (!strcmp( pixname, "353-2" )) return( 2.89369400879);
  if (!strcmp( pixname, "353-3a")) return( 2.8927469582);
  if (!strcmp( pixname, "353-3b")) return( 2.84059437282);
  if (!strcmp( pixname, "353-4a")) return( 2.75473230984);
  if (!strcmp( pixname, "353-4b")) return( 2.96102816276);
  if (!strcmp( pixname, "353-5a")) return( 2.81750207098);
  if (!strcmp( pixname, "353-5b")) return( 3.02236882148);
  if (!strcmp( pixname, "353-6a")) return( 2.76555570713);
  if (!strcmp( pixname, "353-6b")) return( 2.73805574865);
  if (!strcmp( pixname, "353-7" )) return( 2.85609646689);
  if (!strcmp( pixname, "353-8" )) return( 2.8853212746);
  return( 100.0); // no despike_flagging
}
