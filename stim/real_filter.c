
#define _XOPEN_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "real_filter.h"
#include "stim_tools.h"


#define FILTERLEN (3)
float FILTER[FILTERLEN] = {0.25, 0.5, 0.25};


int real_filter( stimParameters *stimPar, PIOFLOAT *signal) {

  long i, j;
  long sig_length = stimPar->total_signal_length;
  
  STIM_TRACE( 1, "real-space filtering for rings %d-%d", stimPar->BeginRing, stimPar->EndRing);
  assert( FILTERLEN % 2 == 1);

  float *filtsig = (float *) malloc( sig_length * sizeof( float));
  assert( filtsig != NULL);

  // filter signal in interval where filter is complete
  for (i = FILTERLEN/2; i < sig_length - FILTERLEN/2; i++) {
    filtsig[i] = 0.0;
    for (j = 0; j < FILTERLEN; j++) {
      filtsig[i] += FILTER[j] * signal[i - FILTERLEN/2 + j];
    }
  }

  // copy back to signal, simply propagating to values outside of filter definition interval
  for (i = 0; i < FILTERLEN/2; i++) {
    signal[i] = filtsig[FILTERLEN/2];
  }
  for (i = FILTERLEN/2; i < sig_length - FILTERLEN/2; i++) {
    signal[i] = filtsig[i];
  }
  for (i = sig_length - FILTERLEN/2; i < sig_length; i++) {
    signal[i] = filtsig[sig_length - FILTERLEN/2 - 2];
  }

  free( filtsig);

  return( 0);
}
