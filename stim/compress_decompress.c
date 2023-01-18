
#define _XOPEN_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "compress_decompress.h"
#include "stim_tools.h"


int compress_decompress( stimParameters *stimPar, PIOFLOAT *signal) {

  double CompStep = stimPar->brimo.compstep;

  long i, j, comp_chunk_start;
  double mean1, mean2;

  STIM_TRACE( 1, "DNLSIM compressing/decompressing rings %d-%d with compstep=%f", stimPar->BeginRing, stimPar->EndRing, CompStep);
  comp_chunk_start = (stimPar->signal_first_sample_number / 254 + 1) * 254 - stimPar->signal_first_sample_number;
  for (i = comp_chunk_start; i < stimPar->total_signal_length - 254; i += 254) {
    mean1 = 0.0;
    mean2 = 0.0;
    for (j = 0; j < 127; j++) {
      mean1 += signal[i+2*j];
      mean2 += signal[i+2*j+1];
    }
    mean1 /= 127;
    mean2 /= 127;
    // compress / decompress
    for (j = 0; j < 127; j++) {
      signal[i+2*j]   = round( (signal[i+2*j]   - mean1) / CompStep) * CompStep + mean1;
      signal[i+2*j+1] = round( (signal[i+2*j+1] - mean2) / CompStep) * CompStep + mean2;
    }
  }
  return( 0);
}


/******************************************************************************/
/* DESIRE implementation */
/******************************************************************************/

// Desire lib types
typedef uint32_t U32BIT;
typedef int32_t  S32BIT;

float MeanCalculate (U32BIT* inData, U32BIT nbSamples, U32BIT sampling) {
// compute the mean of the slice,
// after having removed the min and max values of the slice
// if their deviation is more than 5 times the mean deviation

  U32BIT glitchCut  = 5;

  float meanChannel = 0.0;
  float newMean     = 0.0;
  U32BIT iValue     = 0;
  U32BIT nbValues   = 0;

  U32BIT data;
  U32BIT minVal     = 0xFFFFFFFF;
  U32BIT maxVal     = 0;
  U32BIT meanDev    = 0;
  U32BIT maxDev     = 0;
  U32BIT dev        = 0;

  for (iValue=0; iValue<nbSamples; iValue++) {
    data = inData[iValue];
    if (data > 0) {
      meanChannel += (float)data;
      if (iValue > 0) {
        dev = (U32BIT) abs( (S32BIT)inData[iValue] - (S32BIT)inData[iValue-1]);
        meanDev += dev;
        if (dev>maxDev) maxDev = dev;
      }
      if (data < minVal) minVal = data;
      if (data > maxVal) maxVal = data;
      nbValues++;
    }
  }

  if (nbValues > 0) {
    meanChannel /= nbValues;
  }

  //enleve-t-on minVal?
  if (nbValues > 1) {
    //on enleve maxDev a meanDev:
    meanDev /= (nbValues-1);
    if (nbValues > 2) {
      meanDev = (U32BIT) abs( (S32BIT)((nbValues-1)*meanDev) - (S32BIT)maxDev) / (nbValues-2);
    }
    //on enleve minVal a meanChannel
    newMean = (float)((S32BIT)nbValues*meanChannel - (S32BIT)minVal) / (nbValues-1);
    if ((U32BIT)abs((S32BIT)newMean-(S32BIT)minVal) > (glitchCut*meanDev)) {
      //printf("MeanCalculate: removing min=%d old mean=%f new mean=%f\n",minVal,meanChannel,newMean);
      meanChannel = newMean;
      nbValues--;
    }
  }

  //enleve-t-on maxVal?
  if (nbValues > 1) {
    //CM: deja divise
    //meanDev /= (nbValues-1);
    if (nbValues > 2) {
      meanDev = (U32BIT)abs((S32BIT)(nbValues-1)*meanDev - (S32BIT)maxDev) / (nbValues-2);
    }
    newMean = (float)((S32BIT)nbValues*meanChannel - (S32BIT)maxVal) / (nbValues-1);
    if ((U32BIT)abs((S32BIT)maxVal-(S32BIT)newMean) > (glitchCut*meanDev)) {
      //printf("MeanCalculate: removing max=%d old mean=%f new mean=%f\n",maxVal,meanChannel,newMean);
      meanChannel = newMean;
    }
  }
  return (float)meanChannel;
} /* End of function MeanCalculate */


long ChannelCpSliceScaling( U32BIT* inData, U32BIT nbSamples, float offset, U32BIT sampling) {
  U32BIT iSamp=1;
  S32BIT scaled;
  long errors=0;

  while (iSamp<nbSamples) {
    // demodulate positive parity signal to negative parity signal
    scaled = 2 * (U32BIT) (offset+0.5) - inData[iSamp];
    if (scaled >= 0) {
      inData[iSamp]=(U32BIT)scaled;
    }
    else {
//      printf("ChannelCpSliceScaling: error sample=%d ofset=%f putting 0\n",inData[iSamp],offset);
      errors += 1;
      inData[iSamp]=0;
    }
    iSamp += 2*(iSamp%2)*(sampling-1)+2;
  }
  return errors;
}


float MeanCut( U32BIT* inData, U32BIT nbSamples, float sigma) {

  /* moyenne en retirant les paires d'echantillons qui fluctuent a plus de glitchCut sigma */
  /* si pas assez d'echantillon a la fin (<minSample) prend la moyenne brute */
  /*                                                     S.Plaszczynski 5 feb 07   */
  /*-----------------------------------------------------------------------------------------------------*/

  U32BIT glitchCut     = 10;
  float fracminSamples = 0.90;

  U32BIT  meanChannel  = 0;
  U32BIT newMean       = 0;
  float zeMean         = 0.0; 


  U32BIT data   = 0;
  U32BIT iValue = 0;
  
  float dev     = 0;
  float maxDev  = glitchCut * sigma;

  U32BIT nbValues    = 0;
  U32BIT nbCutValues = 0;

  if (nbSamples !=0) {
    zeMean= (float)inData[0];
  }
  
  //la boucle commence a 1
  for(iValue=1; iValue<nbSamples; iValue++) {
    data = inData[iValue];
    /*printf("%d data=%d\n",iValue,data);*/

    //protection si la moyenne n'a pas pu etre soustraite avant:
    if (data > 0) {
      meanChannel += data;
      dev = (float) abs( (S32BIT)inData[iValue] - (S32BIT)inData[iValue-1]);
      /*printf("%d deviation =%f\n",iValue,dev);*/
      if (dev < maxDev) {
        newMean += data;
        nbCutValues++;
      }
      nbValues++;
    }
  }

  if ((nbValues != 0) && ((float)nbCutValues < fracminSamples * nbValues)) {
    //printf("*******MeanCut: pas assez d'echantillons nval=%d nbcut=%d sigma=%f maxdev=%f\n",nbValues,nbCutValues,sigma,maxDev);
    zeMean = (float)meanChannel / (float)nbValues;
  }
  else {
    //printf("MeanCut: nombre d'echantillons %d\n",nbCutValues);
    if (nbCutValues !=0) {
      zeMean = (float)newMean / (float)nbCutValues;
    }
  }
  return zeMean;
}


int sampleRegeneration( int inData, int quantifStep, int mean, int ref, int parityData) {
  int outData;

//	DataOut=(double)moy+(double)quantif*(0.5+DataIn);
  outData = mean + quantifStep * inData;

  if (parityData%2) {
    outData = 2 * ref - outData;
  }

  return(outData);
}


int desire_codec( stimParameters *stimPar, PIOFLOAT *signal) {

  const unsigned csize=254;
//  float sig_Q = 19.0;
  float sig_Q = 2.5;
  long i, j, islice;

  int CompStep = stimPar->brimo.compstep;
  long comp_start = (stimPar->signal_first_sample_number / csize + 1) * csize - stimPar->signal_first_sample_number;
  long N = (stimPar->total_signal_length - comp_start) / csize * csize;
  float sigEst = CompStep * sig_Q;
  int firstPar = PARITY( stimPar->signal_first_sample_number + comp_start);

  STIM_TRACE( 1, "DESIRE compressing/decompressing rings %d-%d with sigma=%g, compstep=%d", stimPar->BeginRing, stimPar->EndRing, sigEst, CompStep);

  U32BIT *adu = malloc( N * sizeof( U32BIT));
  for (i=0; i<N; i++) {
    adu[i] = signal[comp_start + i];
  }

  long scaling_errors = 0;

  //on saute la fin
  //loop on CS
  for (islice=0; islice<N; islice+=csize) {
    U32BIT *p_inScData = &adu[islice];

    //demodulation: modifies in place
    float ref = MeanCalculate( p_inScData, csize, 1);
    scaling_errors += ChannelCpSliceScaling( p_inScData, csize, ref, 1);
    float channelMean = MeanCut( p_inScData, csize, sigEst);

    //COMPRESSION
    U32BIT quantifStep = sigEst / sig_Q;

    float sample, outData;
    S32BIT compData;
    S32BIT bin[csize];
    for (j=0; j<csize; j++) {
      sample = (float)adu[islice+j];
      outData = floorf( (sample - channelMean + quantifStep / 2.0) / quantifStep);
      compData = (S32BIT)outData;
      //sp: si sample=0 il y a erreur : mettre > 16bits
      if ((U32BIT)sample == 0) {
        compData=0xFFFF+1;
      }
      //could put a test on 22bits here
      bin[j] = compData;
    }

    //DECOMPRESSION--->PUS packets
    //headers
    U32BIT compMean = channelMean + 0.5;
    U32BIT compQuantStep = quantifStep + 0.5;
    U32BIT compRef = ref + 0.5;
    int parity = firstPar;
    S32BIT data2Recup;

    long ipos = islice;
    int dpData;
    for (j=0; j<csize; j++) {
      data2Recup = bin[j]; //lossless
      //recover REU modulated data
      dpData = sampleRegeneration( data2Recup, compQuantStep, compMean, compRef, parity);
      //recopie les donnees
      signal[comp_start + ipos++] = dpData;
      parity++;
    }
  }
  STIM_TRACE( 1, "DESIRE rings %d-%d, scaling errors = %ld (%.2f%%)", stimPar->BeginRing, stimPar->EndRing, scaling_errors, (float)scaling_errors/N*100 );
  return 0;
}
