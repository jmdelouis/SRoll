
#define _XOPEN_SOURCE

#include <assert.h>
#include <math.h>

#include "bl_demod.h"
#include "stim_tools.h"

#include "binary_file_toi.h"

/*==============================================================================

bl_demod(): signal demodulation using a baseline

inputs:
- signal: modulated signal in Volts (adu_to_volts output)

processing:
- signal is smoothed with a length of LSMOOTH (648000 samples, roughly one hour)
and demodulated using this smoothed signal and PARITY

==============================================================================*/


int bl_demod( stimParameters *stimPar, PIOFLOAT *signal) {

  long   i;
  double sliding_sum = 0.0;

  STIM_TRACE( 1, "smoothing modulated signal with length=%ld for rings %d-%d", LSMOOTH, stimPar->BeginRing, stimPar->EndRing);

  // check that we have enough left and right margin to fully smooth all the useful samples
  assert( stimPar->beginring_first_sample >= LSMOOTH/2);
  assert( stimPar->endring_last_sample < stimPar->total_signal_length - LSMOOTH/2);

  // array to temporarily store demodulated signal
  float *demodsig = (float *) malloc( stimPar->total_signal_length * sizeof( float));
  assert(demodsig != NULL);

  // initialise sliding average by computing it for first useful sample at LSMOOTH/2
  for (i = 0; i < LSMOOTH; i++) {
    sliding_sum += signal[i];
  }

  // move LSMOOTH sum through all sample interval
  for (i = LSMOOTH/2; i < stimPar->total_signal_length-LSMOOTH/2; i++) {

    // demodulate signal using smoothed baseline and PARITY
    demodsig[i] = (signal[i] - sliding_sum / LSMOOTH) * PARITY_MOD(i + stimPar->signal_first_sample_number);

    // move to next sample
    sliding_sum += signal[i + LSMOOTH/2] - signal[i - LSMOOTH/2];
  }

  // replace signal with demodulated signal
  for (i = LSMOOTH/2; i < stimPar->total_signal_length-LSMOOTH/2; i++) {
    signal[i] = demodsig[i];
  }

  // clear parts of signal which are now useless
  for (i = 0; i < LSMOOTH/2; i++) {
    signal[i] = 0.0;
    signal[stimPar->total_signal_length-i-1] = 0.0;
  }
  
  free( demodsig);
  return( 0);
}
