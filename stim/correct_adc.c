
#define _XOPEN_SOURCE

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#include <complex.h> 
#include <fftw3.h> 

#include "correct_adc.h"
#include "stim_tools.h"
#include "no_dmc_data_access.h"

#include "binary_file_toi.h"


#define INVSIZE          (10000000l)
#define INVFOURIERSIZE (INVSIZE/2+1)
#define NP                   (4000l)

// debug: set to 1 to save 4K lines value per sample in a TOI
#define SAVE_4KLINES 0

int TRACE_CORRECTADC = 0;


int correct_adc( stimParameters *stimPar, PIOFLOAT *signal) {

  stim_parContent *Param = &stimPar->Param; // defined in sroll/stim/stim_param.h
  stim_Data       *data  = &stimPar->data;
  BRIMO           *brimo = &stimPar->brimo; // defined in sroll/brimo.h

  int BeginRing  = stimPar->first_full_ring;
  int EndRing    = stimPar->last_full_ring;
  int sphase     = brimo->sphase;

  long i, j;

  STIM_TRACE( 1, "correcting ADC nonlinearity and removing 4K lines for rings %d-%d", BeginRing, EndRing);

  if (data->raw_begin_ring == 0) {
    // raw data not already read, will start at BeginRing
    data->raw_begin_ring = BeginRing;
    if (data->raw_begin_ring < 240) {
      data->raw_begin_ring = 240;
    }
  }
  int ring_count = EndRing - data->raw_begin_ring + 1;

  //////////////////////////////////////////////////////////////////////////////
  // read 4K vect

  read_4k_vect( Param->fourk_name, data->raw_begin_ring, ring_count, data);

  //////////////////////////////////////////////////////////////////////////////
  // read RAW CST (dd0 in correctadc.pro)

  STIM_TRACE( 2, "read RAW CST: %s", Param->rawcst_name);

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
  // get RAW GAIN

  double *sigpart = data->raw_gain; // keep corradc.pro names...

  //////////////////////////////////////////////////////////////////////////////
  // read and invert INL

  STIM_TRACE( 2, "invert INL: %s", Param->corr_inl_name);

  // invert INL
  double *sss = (double *) malloc( sizeof( double) * INVSIZE);
  assert( sss != NULL);
  for (double dbli = 0; dbli < INVSIZE; dbli++) {
    sss[(int) dbli] = dbli / INVSIZE * INLSIZE - INLSIZE/2;
  }

  double         *inlinv =         fftw_alloc_real(    INVSIZE);
  double complex *inlinv_fourier = fftw_alloc_complex( INVFOURIERSIZE);
  assert( inlinv != NULL);
  assert( inlinv_fourier != NULL);
  memset( inlinv_fourier, 0, sizeof(double complex) * INVFOURIERSIZE);
  
  fftw_plan fft_inlinv_forward  = fftw_plan_dft_r2c_1d( INVSIZE, inlinv, inlinv_fourier, FFTW_ESTIMATE);
  fftw_plan fft_inlinv_backward = fftw_plan_dft_c2r_1d( INVSIZE, inlinv_fourier, inlinv, FFTW_ESTIMATE);

  for (i = 10 * INVSIZE / (INLSIZE/2); i < INVSIZE - 10 * INVSIZE / (INLSIZE/2) - 2; i++) {
    long icou = -10;
    while (data->corr_inl[(int) sss[i] + icou + INLSIZE/2] < (sss[i] + INLSIZE/2)) icou++;
    inlinv[i] = (int) sss[i] + icou - 1;
  }
  free( sss);

  // transform inverse INL to Fourier space into inlinv_fourier[]
  fftw_execute( fft_inlinv_forward);

  // build smoothing filter
  // smoothing_filter = fft(exp(-dist(n_elements(adu),1)^2d0/2d0/(4d0*n_elements(adu)/65536d0)^2d0),-1)
  double         *filter =         fftw_alloc_real(    INVSIZE);
  double complex *filter_fourier = fftw_alloc_complex( INVFOURIERSIZE);
  assert( filter != NULL);
  assert( filter_fourier != NULL);
  memset( filter_fourier, 0, sizeof(double complex) * INVFOURIERSIZE);

  fftw_plan fft_filter_forward = fftw_plan_dft_r2c_1d( INVSIZE, filter, filter_fourier, FFTW_ESTIMATE);

  for (i=0; i<=INVSIZE/2; i++) {
    filter[i] = exp(-((double)i*i)/2.0/(4.0*INVSIZE/INLSIZE)/(4.0*INVSIZE/INLSIZE));
    if (i != 0) {
      filter[INVSIZE-i] = filter[i];
    }
  }

  // transform smoothing filter to Fourier space into filter_fourier[]
  fftw_execute( fft_filter_forward);

  // convolve inverse INL with smoothing filter
  double filter_max = 0;
  for (i=0; i<INVFOURIERSIZE; i++) {
    if (0) {
      // this expression drives valgrind crazy, not the next one...
      inlinv_fourier[i] *= filter_fourier[i];
    } else {
      inlinv_fourier[i] =   creal(inlinv_fourier[i]) * creal(filter_fourier[i])
                          - cimag(inlinv_fourier[i]) * cimag(filter_fourier[i])
                          + I * ( creal(inlinv_fourier[i]) * cimag(filter_fourier[i])
                                + cimag(inlinv_fourier[i]) * creal(filter_fourier[i]));
    }
    if (creal( filter_fourier[i]) > filter_max) {
      filter_max = creal( filter_fourier[i]);
    }
  }

  // back to time domain and normalise
  fftw_execute( fft_inlinv_backward);
  for (i=0; i<INVSIZE; i++) {
    inlinv[i] /= filter_max * INVSIZE;
  }

  // FFTW cleanup
  fftw_destroy_plan( fft_inlinv_forward);
  fftw_destroy_plan( fft_filter_forward);
  fftw_destroy_plan( fft_inlinv_backward);
  fftw_free( inlinv_fourier);
  fftw_free( filter_fourier);
  fftw_free( filter);


  //////////////////////////////////////////////////////////////////////////////
  // process signal per ring

  double *ccp0 = (double *) malloc( sizeof( double) * NP);
  assert (ccp0 != NULL);
  double *ccp1 = (double *) malloc( sizeof( double) * NP);
  assert (ccp1 != NULL);
  double *ccm0 = (double *) malloc( sizeof( double) * NP);
  assert (ccm0 != NULL);
  double *ccm1 = (double *) malloc( sizeof( double) * NP);
  assert (ccm1 != NULL);
  double ring_4k[80*9];

  for (int ring = BeginRing; ring <= EndRing; ring++) {

    long begin_ring_index = BEGINRINGINDEX( ring);
    long ns = ENDRINGINDEX( ring) - begin_ring_index + 1;
    int PARITY0 = PARITY( begin_ring_index);
    float *ring_signal = signal + begin_ring_index - stimPar->signal_first_sample_number;
    double corrp[9][NP][2];
    double corrm[9][NP][2];

    STIM_TRACE( 2, "ring %d, sample count=%ld", ring, ns);
    if (stimPar->trace_level >= 2) {
      long mem_total, mem_used, mem_free;
      stim_FreeMem( &mem_total, &mem_used, &mem_free);
      STIM_TRACE( 2, "ring %d, mem total=%ldMB, used=%ldMB, free=%ldMB", ring, mem_total, mem_used, mem_free);
    }
  
    assert( ring_signal + ns <= signal + stimPar->total_signal_length);

    double *fourkpersample = NULL;
    if (SAVE_4KLINES) {
      fourkpersample = (double *) malloc( ns * sizeof( double));
      assert( fourkpersample != NULL);
    }

    double convf = (double) INVSIZE/INLSIZE;

    // apply a circular permutation of sphase on ring raw constant
    int raw_ring = ring - data->raw_begin_ring;
    if ( raw_ring < 0) {
      raw_ring = 0;
    }
    double dstat[80];
    for (i=0; i<80; i++) {
      dstat[i] = data->raw_cst[raw_ring * 80 + (i+sphase) % 80] - INLSIZE/2;
      assert( (!isnan(dstat[i])) && (!isinf(dstat[i])));
    }

    STIM_TRACE( 2, "ring %d, computing 4K lines", ring);
    compute_4K_lines( raw_ring, begin_ring_index, data, ring_4k);

    for (int isamp=0; isamp<=8; isamp++) {

      // produce a ramp signal (sig)
      // convert it to ADU (corrp[][][1] for parity 0, corrm[][][1] for parity 1)
      // and pass it through ADCNL (corrp[][][0] for parity 0, corrm[][][0] for parity 1)
      for (i=0; i<NP; i++) {
        double sig = ((double)i - NP/2) * 25.0;
        if (i == 0)    sig = -0.08 * INVSIZE;
        if (i == NP-1) sig = +0.08 * INVSIZE;

        double raww[80];
        double mraw[80];
        for (j=0; j<80; j++) {
          raww[j] = dstat[j] + ring_4k[isamp*80 + j] + sig * sigpart[j];
          mraw[j] = inlinv[(long)((raww[j] + INLSIZE/2) * convf)];
        }

        corrp[isamp][i][0] = 0.0;
        corrp[isamp][i][1] = 0.0;
        corrm[isamp][i][0] = 0.0;
        corrm[isamp][i][1] = 0.0;
        for (j=0; j<40; j++) {
          corrp[isamp][i][1] += raww[j];
          corrm[isamp][i][1] += raww[j+40];
         
          corrp[isamp][i][0] += mraw[j];
          corrm[isamp][i][0] += mraw[j+40];
        }
      }
    }

    for (int isamp=0; isamp<=8; isamp++) {

      long   *selecp = NULL;
      long   n_out = 0;
      long   selecp_len = 0;
      double total_raw4kall;

      selecp = (long *) malloc( sizeof( long) * (ns/18+2)*2+1);
      assert( selecp != NULL);

      // IDL: selecp = reform([transpose(dindgen(ns/18+2)*18)+isamp*2,transpose(dindgen(ns/18+2)*18+isamp*2+1)],(ns/18+2)*2) + pparity(0)
      for (i=0; i<ns/18+2; i++) {
        selecp[2*i]   = i * 18 + isamp * 2     + PARITY0;
        selecp[2*i+1] = i * 18 + isamp * 2 + 1 + PARITY0;
        if (selecp[2*i]   < ns) selecp_len = 2*i;
        if (selecp[2*i+1] < ns) selecp_len = 2*i+1;
      }
      selecp_len += 1; // length of array of sample values less than the total number of samples in ring
      if ((isamp == 8) && (PARITY0 == 1)) { // ? special case
        selecp[selecp_len++] = 0;
      }
      assert( selecp_len <= (ns/18+2)*2+1);


      // correct parity 0
      STIM_TRACE( 2, "ring %d, correcting parity 0 for isamp=%d", ring, isamp);

      // IDL: ccp = reform(corrp(isamp mod 9,*,*))
      for (i=0; i<NP; i++) {
        ccp0[i] = corrp[isamp][i][0];
        ccp1[i] = corrp[isamp][i][1];
      }

      // IDL: dd(selecp(we))
      double *dd_selp = NULL;  // parity 0 signal to correct
      dd_selp = (double *) malloc( sizeof( double) * (selecp_len / 2 + 1));
      assert( dd_selp != NULL);
      j = 0;
      for (i=0; i<selecp_len; i++) {
        if (PARITY( begin_ring_index + selecp[i]) == 0) {
          dd_selp[j++] = ring_signal[selecp[i]] - INLSIZE/2 * 40.0;
        }
      }
      n_out = j;
      assert( n_out <= selecp_len / 2 + 1);

      STIM_TRACE( 2, "ring %d, before raw4k", ring);
      // IDL: total(raw4kall(0:39,isamp))
      total_raw4kall = 0.0;
      for (j = 0; j < 40; j++) {
        total_raw4kall += ring_4k[isamp*80 + j];
      }

      STIM_TRACE( 2, "ring %d, before linear interpolation", ring);
      // IDL: ddc(selecp(we)) = interpol(ccp(*,1),ccp(*,0),dd(selecp(we)))
      double *ddc_selp = NULL;  // parity 0 corrected signal
      ddc_selp = (double *) malloc( sizeof( double) * n_out);
      assert( ddc_selp != NULL);
      linear_interpolate( ccp0, ccp1, dd_selp, ddc_selp, NP, n_out, TRACE_CORRECTADC);

      STIM_TRACE( 2, "ring %d, before signal correction", ring);
      // IDL: ddc = ddc - dd4k + 32768l*40
      j = 0;
      for (i=0; i<selecp_len; i++) {
        if (PARITY( begin_ring_index + selecp[i]) == 0) {
          if (ring_signal[selecp[i]] != 0.0) {
            ring_signal[selecp[i]] = trunc( ddc_selp[j] - total_raw4kall + INLSIZE/2 * 40);
          }
          j++;
          assert( !isnan( ring_signal[selecp[i]]));

          if (fourkpersample != NULL) {
            fourkpersample[selecp[i]] = total_raw4kall;
          }
        }
      }
      free( dd_selp);
      free( ddc_selp);


      // correct parity 1
      STIM_TRACE( 2, "ring %d, correcting parity 1 for isamp=%d", ring, isamp);

      // IDL: ccm = reform(corrm(isamp mod 9,*,*))
      for (i=0; i<NP; i++) {
        ccm0[i] = corrm[isamp][i][0];
        ccm1[i] = corrm[isamp][i][1];
      }

      // IDL: dd(selecp(wee))
      double *dd_selm = NULL;  // parity 1 signal to correct
      dd_selm = (double *) malloc( sizeof( double) * (selecp_len / 2 + 1));
      assert( dd_selm != NULL);
      j = 0;
      for (i=0; i<selecp_len; i++) {
        if (PARITY( begin_ring_index + selecp[i]) == 1) {
          dd_selm[j++] = ring_signal[selecp[i]] - INLSIZE/2 * 40.0;
        }
      }
      n_out = j;
      assert( n_out <= selecp_len / 2 + 1);

      // IDL: total(raw4kall(40:79,isamp))
      total_raw4kall = 0.0;
      for (j = 40; j < 80; j++) {
        total_raw4kall += ring_4k[isamp*80 + j];
      }

      // IDL: ddc(selm) = interpol(ccm(*,1),ccm(*,0),dd(selm))
      double *ddc_selm = NULL;
      ddc_selm = (double *) malloc( sizeof( double) * n_out); // parity 1 corrected signal
      assert( ddc_selm != NULL);
      linear_interpolate( ccm0, ccm1, dd_selm, ddc_selm, NP, n_out, TRACE_CORRECTADC);

      // IDL: ddc = ddc - dd4k + 32768l*40
      j = 0;
      for (i=0; i<selecp_len; i++) {
        if (PARITY( begin_ring_index + selecp[i]) == 1) {
          if (ring_signal[selecp[i]] != 0.0) {
            ring_signal[selecp[i]] = trunc( ddc_selm[j] - total_raw4kall + INLSIZE/2 * 40);
          }
          j++;
          assert( !isnan( ring_signal[selecp[i]]));

          if (fourkpersample != NULL) {
            fourkpersample[selecp[i]] = total_raw4kall;
          }
        }
      }
      free( dd_selm);
      free( ddc_selm);
      free( selecp);
    }

    if ((fourkpersample != NULL) && (stimPar->mpi_rank == 0)) {
      assert( BFT_writeObject( "/pscratch1/RD12_data/dmc/MISS03/DATA/sylvain2_VEC/143-1a_correctadc_4kpersample2",
                               "float64", begin_ring_index, ns, fourkpersample) >= 0);
      free( fourkpersample);
      fourkpersample = NULL;
    }
  }

  free( ccp0);
  free( ccp1);
  free( ccm0);
  free( ccm1);
  fftw_free( inlinv);
  fftw_cleanup();

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

  return( 0);
}
