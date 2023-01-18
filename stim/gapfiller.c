
#define _XOPEN_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "no_dmc_data_access.h"
#include "binary_file_toi.h"

#include "gapfiller.h"


#define PBRSIZE (10822l)
// The factor to be used in recursiv function decomposeSignalToBinForRing()
#define PBR_FACTOR (10)

//#define _DEBUG_GAPFILLER_
int DEBUG_RING_INFO;


//-----------------------------------------------------------------------------

/**
 * Allow to decompose the input signal in 'bin' using a specified PBRSize.
 * This is a recusive function: it succevely try with wider range of bin until
 * we reach a currentPBRSize < 10
 *
 * Return non zero in case of unsuccess (ie. unable to decomposed the input signal
 * with the specified currentPBRSize for avoiding gap), otherwise return 0 in case
 * of success.
 */
int decomposeSignalToBinForRing(PIOFLOAT *signal,
    PIOFLAG *sigflag, PIODOUBLE *phase, long currentPBRSize,
    long ring_first_sample, long ring_last_sample,
    PIOFLOAT **pbr_signal, PIOINT **pbr_hitcount) {
    int ret = 0;
    //fprintf(stderr, "DEBUG: Attempt to bin with %ld\n", currentPBRSize);

    // Recursion stop condition
    if (currentPBRSize < 10) {
      fprintf(stderr, "DEBUG: STOP RECURSION!\n");
      return ( 1);
    }

    // allocate Phase Binned Ring object
    *pbr_signal   = (PIOFLOAT *) calloc( currentPBRSize, sizeof( PIOFLOAT));
    *pbr_hitcount = (PIOINT *)   calloc( currentPBRSize, sizeof( PIOINT));

    assert(*pbr_signal != NULL);
    assert(*pbr_hitcount != NULL);

    // bin signal in PBR
    for (long i = ring_first_sample; i <= ring_last_sample; i++) {
      if (sigflag[i] == 0) {
        int bin = phase[i] * currentPBRSize / (2 * M_PI);
        (*pbr_signal)[bin] += signal[i];
        (*pbr_hitcount)[bin] += 1;
      }
    }

    PIOFLOAT *sub_pbr_signal = NULL;
    PIOINT   *sub_pbr_hitcount = NULL;

    // Final step for computing PBR
    for (int bin = 0; bin < currentPBRSize; bin++) {
      if ((*pbr_hitcount)[bin] > 0) {
        (*pbr_signal)[bin] /= (*pbr_hitcount)[bin];
      } else { // Case where we need to decompose in wider bin in order to avoid gap zone
        int subPBRSize = currentPBRSize/PBR_FACTOR;
        if (sub_pbr_signal == NULL) { // We need to compute sub_pbr_signal (first time)
          // Here we make use of currentPBRSize/PBR_FACTOR to have wider bin
          // and therefore a better chance to avoid gap zone.
          ret = decomposeSignalToBinForRing(signal, sigflag, phase,
              subPBRSize, ring_first_sample, ring_last_sample, &sub_pbr_signal,
              &sub_pbr_hitcount);
          if (ret != 0) { // Case of error
            break;
          }
        }

        // We already have compute sub_pbr_signal without error, so we use it...

        // Find the subBin index corresponding to our current bin
        // Integer division: so we only keep the integer part.
        int subbin = (bin * subPBRSize) / currentPBRSize;

        // Use the corresponding value as a "good" replacement for missing value at "bin".
        (*pbr_signal)[bin] = sub_pbr_signal[subbin];
#ifdef _DEBUG_GAPFILLER_
        // Just for debug we store the subPBRSize used for filling the gap
        (*pbr_hitcount)[bin] = sub_pbr_hitcount[subbin];
#endif

        fprintf(stderr, "[INFO] Gapfiller(): pbr bin #%d set to a subbinned value (%d)!\n", bin, subPBRSize);
        DEBUG_RING_INFO = 1;
      }
    }

    // Free unused memory
    if (sub_pbr_signal != NULL) {
      free(sub_pbr_signal);
    }
    if (sub_pbr_hitcount != NULL) {
      free(sub_pbr_hitcount);
    }

  return ( ret);
}

//-----------------------------------------------------------------------------

int gapfiller( stimParameters *stimPar, PIOFLOAT *signal, PIOFLAG *sigflag) {
  //fprintf(stderr, "* %s()\n", __FUNCTION__);

  int ret = 0;

  stim_parContent *Param = &stimPar->Param; // defined in sroll/stim/stim_param.h
  //BRIMO           *brimo = &stimPar->brimo; // defined in sroll/brimo.h

  // Read phase object
  PIODOUBLE *phase = (PIODOUBLE *) malloc( sizeof( PIODOUBLE) * stimPar->total_signal_length);
  assert( phase != NULL);
  assert( noDMC_readObject_PIODOUBLE( Param->phase,
                                     stimPar->signal_first_sample_number,
                                     stimPar->total_signal_length,
                                     phase) >= 0);

#if 0
  // Read badring object
  //fprintf(stderr, "[DEBUG] Param->badring='%s'\n", Param->badring);
  PIOINT *badring = (PIOINT *) malloc( sizeof( PIOINT) * stimPar->total_signal_length);
  assert( badring != NULL);
  assert( noDMC_readObject_PIOINT( Param->badring,
                                   0,
                                   27005,
                                   badring) >= 0);
#endif

#ifdef _DEBUG_GAPFILLER_
  ///fprintf(stderr, "DEBUG: BEGINRINGINDEX(%d) = %ld\n", stimPar->BeginRing, BEGINRINGINDEX( stimPar->BeginRing));
  // Write the input signal
  char str[500] = "\0";
  sprintf(str, "/pscratch1/cmadsen/OUTPUT_FROM_STIM/TMP_DATA_INPUT%d", stimPar->BeginRing);
  assert( BFT_writeObject( str, "float32", stimPar->beginring_first_sample, stimPar->endring_last_sample - stimPar->beginring_first_sample, signal) >= 0);
  // Write also input flag...
  sprintf(str, "/pscratch1/cmadsen/OUTPUT_FROM_STIM/TMP_DATA_INPUT_FLAG%d", stimPar->BeginRing);
  assert( BFT_writeObject( str, "byte", stimPar->beginring_first_sample, stimPar->endring_last_sample - stimPar->beginring_first_sample, sigflag) >= 0);

  // Allocate mem for output
  PIOFLOAT *tmp_data_pbr = (PIOFLOAT *) calloc( PBRSIZE * (stimPar->EndRing - stimPar->BeginRing + 1), sizeof( PIOFLOAT));
  assert( tmp_data_pbr != NULL);
  PIOINT *tmp_data_pbr_hitcount = (PIOINT *) calloc( PBRSIZE * (stimPar->EndRing - stimPar->BeginRing + 1), sizeof( PIOINT));
  assert( tmp_data_pbr_hitcount != NULL);
#endif


  // Process signal per ring
  for (int ring = stimPar->BeginRing; ring <= stimPar->EndRing; ring++) {
#if 0
    // Check that the current ring is not a bad one (otherwise just skip to the next...)
    if (badring[ring] != 0) {
      fprintf(stderr, "[DEBUG] Gapfiller(): detect bad ring %d (just skip it)\n", ring);
      continue;
    }
#endif
    long ring_first_sample = BEGINRINGINDEX( ring) - stimPar->signal_first_sample_number;
    long ring_last_sample  = ENDRINGINDEX(   ring) - stimPar->signal_first_sample_number;

    PIOFLOAT *pbr_signal = NULL;
    PIOINT   *pbr_hitcount = NULL;

    //fprintf(stderr, "[DEBUG] Gapfiller(): processing ring %d\n", ring);
    DEBUG_RING_INFO = 0;
    ret = decomposeSignalToBinForRing(signal, sigflag, phase, PBRSIZE,
        ring_first_sample, ring_last_sample, &pbr_signal, &pbr_hitcount);
    if (ret != 0) {
      fprintf(stderr, "ERROR in Gapfiller(): Unable to obtain appropriate PBR (even with adaptative resolution) for ring %d!\n", ring);

      break; // Stop processing!
    }

    // DEBUG
    if (DEBUG_RING_INFO != 0) {
      fprintf(stderr, "[DEBUG] Ring for which we have used subbinned: %d (see previous log)\n", ring);
    }

    // Now that we have PBR, we can "correct" gap in the input signal
    // ie. everywhere an invalid flag is set we replace value by the one from the PBR
    for (long i = ring_first_sample; i <= ring_last_sample; i++) {
      if (sigflag[i] != 0) { // Indicate an invalid value in input signal
        int bin = phase[i] * PBRSIZE / (2 * M_PI);
        signal[i] = pbr_signal[bin];

#ifdef _DEBUG_GAPFILLER_
        tmp_data_pbr[(ring-stimPar->BeginRing) * PBRSIZE + bin] = pbr_signal[bin];
        tmp_data_pbr_hitcount[(ring-stimPar->BeginRing) * PBRSIZE + bin] = pbr_hitcount[bin];
#endif
      }
    }

    free( pbr_signal);
    free( pbr_hitcount);
  }

  free( phase);

#ifdef _DEBUG_GAPFILLER_
  // Write updated signal
  sprintf(str, "/pscratch1/cmadsen/OUTPUT_FROM_STIM/TMP_DATA_OUTPUT%d", stimPar->BeginRing);
  assert( BFT_writeObject( str, "float32", stimPar->beginring_first_sample, stimPar->endring_last_sample - stimPar->beginring_first_sample, signal) >= 0);
  // Write PBR
  sprintf(str, "/pscratch1/cmadsen/OUTPUT_FROM_STIM/TMP_DATA_PBR_OUTPUT%d", stimPar->BeginRing);
  assert( BFT_writeObject( str, "float32", 0, PBRSIZE * (stimPar->EndRing - stimPar->BeginRing + 1), tmp_data_pbr) >= 0);
  sprintf(str, "/pscratch1/cmadsen/OUTPUT_FROM_STIM/TMP_DATA_PBR_OUTPUT_HITCOUNT%d", stimPar->BeginRing);
  assert( BFT_writeObject( str, "int32", 0, PBRSIZE * (stimPar->EndRing - stimPar->BeginRing + 1), tmp_data_pbr_hitcount) >= 0);
  free( tmp_data_pbr);
  free( tmp_data_pbr_hitcount);
#endif

  return( ret);
}
