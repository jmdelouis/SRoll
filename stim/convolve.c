
#define _GNU_SOURCE

#include <unistd.h>

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <fcntl.h>
#include <omp.h>

#include <complex.h> 
#include <fftw3.h> 

#include "convolve.h"
#include "taudeconv_lib.h"
#include "stim_tools.h"


// replaces read_imofile() from tau_deconv.cc
void fill_partaudeconv( stimParameters *stimPar, struct par_taudeconv *tauparams) {

  BRIMO *brimo = &stimPar->brimo;

  tauparams->fsampling = SAMPLING_FREQ;

  strcpy(  tauparams->method,     "LFER8");
  strcpy(  tauparams->n_method,   "NONE");
  strncpy( tauparams->lpf_method, stimPar->conv_LPF_FILTER, PIOSTRINGMAXLEN);
  strncpy( tauparams->r_method,   stimPar->conv_R_FILTER,   PIOSTRINGMAXLEN);

  // Check convolution flag 
  tauparams->convol = stimPar->conv_CONVOLVE;

  // READ TIME SHIFT PARAMETERS    
  tauparams->sphase_offset = 0.0;
  tauparams->global_offset = brimo->LFER_global_offset;
  tauparams->total_offset  = tauparams->sphase_offset + tauparams->global_offset;   

  // Read TF clip parameter
  tauparams->tfmodclip = 1.0e-5;

  // LFER8 parameters
  tauparams->tau1 = brimo->LFER_tau1;
  tauparams->tau2 = brimo->LFER_tau2;
  tauparams->tau3 = brimo->LFER_tau3;
  tauparams->tau4 = brimo->LFER_tau4;
  tauparams->tau5 = brimo->LFER_tau5;
  tauparams->tau6 = brimo->LFER_tau6;
  tauparams->tau7 = brimo->LFER_tau7;
  tauparams->tau8 = brimo->LFER_tau8;
  tauparams->A1   = brimo->LFER_A1;
  tauparams->A2   = brimo->LFER_A2;
  tauparams->A3   = brimo->LFER_A3;
  tauparams->A4   = brimo->LFER_A4;
  tauparams->A5   = brimo->LFER_A5;
  tauparams->A6   = brimo->LFER_A6;
  tauparams->A7   = brimo->LFER_A7;
  tauparams->A8   = brimo->LFER_A8;
  tauparams->tau_stray = brimo->LFER_tau_stray;
  tauparams->sphase    = brimo->LFER_sphase;

  // lowpass filter parameters
  if (strcmp( stimPar->conv_LPF_FILTER, "GAUSSCOSSQR") == 0) {
    tauparams->fgauss_lpf = brimo->lpf_fgauss;
    tauparams->fc_lpf     = brimo->lpf_fc;
    tauparams->ffact_lpf  = brimo->lpf_ffact;
  } else {
//    printf( "stimPar->conv_LPF_FILTER = %s\n", stimPar->conv_LPF_FILTER);
    assert( strcmp( stimPar->conv_LPF_FILTER, "NONE") == 0);
  }

  // regularisation filter parameters
  if (strcmp( stimPar->conv_R_FILTER, "COSINE") == 0) {
    tauparams->fwidth_rfilter = brimo->fwidth_rfilter;
  } else {
    assert( strcmp( stimPar->conv_R_FILTER, "NONE") == 0);
  }
}


// replace tau_deconv() and timeresp_deconv() functions from tau_deconv.cc
int convolve( stimParameters *stimPar, PIOFLOAT *signal) {

//  stim_parContent *Param = &stimPar->Param;
  struct par_taudeconv tauparams;

  long     i, j;
  long     conv_start = 0;
  long     conv_end = stimPar->total_signal_length; // avoid -Werror=maybe-uninitialized
  long     global_conv_start, global_conv_end;
  PIOFLOAT *signal_out = NULL;

  double *toideconv = NULL;
  double *timeresp = NULL;

  double         *fftw_time_data = NULL;
  double complex *fftw_fourier_data = NULL;
  fftw_plan       fft_forward = NULL, fft_backward = NULL;

  const long ndata              = 2 * LDECONV;
  const long fftw_complex_ndata = ndata / 2 + 1;

  struct timeval t0;
  gettimeofday( &t0, NULL);

  // initialise tau_deconv library parameters
  fill_partaudeconv( stimPar, &tauparams);

  if (stimPar->conv_CONVOLVE == CONVOLVE) {
    STIM_TRACE( 1, "time constant convolution of rings %d-%d, Rfilter=%s, LPfilter=%s", stimPar->BeginRing, stimPar->EndRing, stimPar->conv_R_FILTER, stimPar->conv_LPF_FILTER);
    global_conv_start = stimPar->conv_start - stimPar->signal_first_sample_number;
    global_conv_end   = stimPar->conv_end - stimPar->signal_first_sample_number;
  } else {
    STIM_TRACE( 1, "time constant deconvolution of rings %d-%d, Rfilter=%s, LPfilter=%s", stimPar->BeginRing, stimPar->EndRing, stimPar->conv_R_FILTER, stimPar->conv_LPF_FILTER);
    global_conv_start = stimPar->deconv_start - stimPar->signal_first_sample_number;
    global_conv_end   = stimPar->deconv_end - stimPar->signal_first_sample_number;
  }

  // initialise FFTW
  assert( fftw_init_threads() != 0);
  fftw_plan_with_nthreads( omp_get_max_threads());
  fftw_time_data    = fftw_alloc_real(    ndata);
  fftw_fourier_data = fftw_alloc_complex( fftw_complex_ndata);
  assert( fftw_time_data    != NULL);
  assert( fftw_fourier_data != NULL);
  // FFTW_MEASURE takes 2 minutes
  fft_forward  = fftw_plan_dft_r2c_1d( ndata, fftw_time_data, fftw_fourier_data, FFTW_ESTIMATE); //, FFTW_MEASURE);
  fft_backward = fftw_plan_dft_c2r_1d( ndata, fftw_fourier_data, fftw_time_data, FFTW_ESTIMATE); //, FFTW_MEASURE);

  // alloc temporary output buffer
  signal_out = (PIOFLOAT *) malloc( sizeof( PIOFLOAT) * stimPar->total_signal_length);
  assert( signal_out != NULL);

  // alloc and build time response (LFER8) array
  timeresp = (double *) malloc( sizeof( double) * ndata);
  assert( timeresp != NULL);
  construct_timeresponse_inv( tauparams, ndata, timeresp);

  // build lowpass filter and multiply time response with it
  if (strcmp( stimPar->conv_LPF_FILTER, "GAUSSCOSSQR") == 0) {
    double *lowpassfilter = (double *) malloc( sizeof( double) * ndata);
    assert( lowpassfilter != NULL);
    construct_lpf_filter( tauparams, ndata, lowpassfilter);
    array_complex_mult( ndata, timeresp, lowpassfilter);
    free( lowpassfilter);
  }

  // build regularisation filter and multiply time response with it
  if (strcmp( stimPar->conv_R_FILTER, "COSINE") == 0) {
    double *regulfilter = (double *) malloc( sizeof( double) * ndata);
    assert( regulfilter != NULL);
    construct_regularization_filter( tauparams, ndata, regulfilter);
    array_complex_mult( ndata, timeresp, regulfilter);
    free( regulfilter);
  }

  // check that time response inversion didn't produce nans (spoiler: it does, at least at 545GHz)
  for (i = 0; i < ndata; i++) {
    if (isnan( timeresp[i])) {
      timeresp[i] = 0;
    }
  }

  // debug save timeresp
  if (0 && (stimPar->mpi_rank == 0)) {
    PIOSTRING filename;
    sprintf( filename, "/wrk/symottet/srollex_DEC16v1/timeresp_%s", stimPar->brimo.pixname);
    float *toto = (float *) malloc( sizeof( float) * fftw_complex_ndata);
    assert( toto != NULL);
    for (i = 1; i < fftw_complex_ndata; i++) toto[i] = cabs( timeresp[i]);
    int timeresp_filedesc = open( filename, O_RDWR | O_CREAT, S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH);
    assert( timeresp_filedesc >= 0);
    assert( pwrite64( timeresp_filedesc, toto, sizeof( float) * fftw_complex_ndata, 0) == sizeof( float) * fftw_complex_ndata); 
    close( timeresp_filedesc);
    free( toto);
  }

  // convolve signal in overlapping 2*LDECONV chunks
  for (i = global_conv_start; i <= global_conv_end; i += LDECONV) {
    conv_start = i - HALFPERIOD;
    conv_end   = i + LDECONV + HALFPERIOD - 1; // included
    assert( conv_start >= 0);
    assert( conv_end < stimPar->total_signal_length );
    STIM_TRACE( 2, "processing samples %ld-%ld", stimPar->signal_first_sample_number + conv_start,
                                                 stimPar->signal_first_sample_number + conv_end);

    // copy input signal to FFT input buffer
    for (j = 0; j < ndata; j++) {
      fftw_time_data[j] = signal[conv_start + j];
    }

    // convert signal in fftw_time_data to frequency domain in fftw_fourier_data
    fftw_execute( fft_forward);

    // convert FFTW complex output to NAG format
    toideconv = fftw_time_data; // reuse now useless FFTW input double array
    for (j = 0; j <= ndata / 2; j++) {
      toideconv[j] = creal(fftw_fourier_data[j]);
    }
    for (j = ndata / 2 + 1; j < ndata ; j++) {
      toideconv[j] = cimag(fftw_fourier_data[ndata - j]);
    }

    // multiply fftw_fourier_data by time response, lowpass filter and regularisation filter
    array_complex_mult( ndata, toideconv, timeresp);   

    // convert NAG format back to FFTW complex array
    fftw_fourier_data[0] = toideconv[0];
    for (j = 1; j <= ndata / 2; j++) {
      fftw_fourier_data[j] = (toideconv[j] + I * toideconv[ndata - j]);
    }

    // convert fftw_fourier_data back to fftw_time_data in time domain
    // WARNING: fftw_plan_dft_c2r_1d() destroys fftw_fourier_data
    fftw_execute( fft_backward);

    // copy result of convolution to temporary signal_out buffer
    // and apply normalisation factor left out by FFTW (1/ndata)
    for (j = i; j < i + LDECONV; j++) {
      signal_out[j] = fftw_time_data[j - conv_start] / ndata;
    }
  }

  // copy convolution result to signal buffer
  for (i = global_conv_start; i <= global_conv_end; i++) {
    signal[i] = signal_out[i];
  }

  free( signal_out);
  free( timeresp);
  fftw_destroy_plan( fft_forward);
  fftw_destroy_plan( fft_backward);
  fftw_free( fftw_time_data);
  fftw_free( fftw_fourier_data);
  fftw_cleanup_threads();

  STIM_TRACE( 1, "duration: %s", get_elapsed_time( &t0));

  return( 0);
}
