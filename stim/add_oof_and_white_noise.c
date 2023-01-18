#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>   
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <math.h>

#include "add_oof_and_white_noise.h"
#include "stim_tools.h"


float GetRD12PhotonicWhiteNoise( char* pixname);
float GetRD12ElectronicWhiteNoise( char* pixname);
float GetWhiteNoiseAlignFactor( char* pixname);


int add_oof_and_white_noise( stimParameters *stimPar, PIOFLOAT *signal, int noise_type) {

  stim_parContent *Param = &stimPar->Param;
  BRIMO           *brimo = &stimPar->brimo;

  // process parameters
  float wn_level;
  switch (noise_type) {
    // for white noise level, use do_*_noise parameter value, or Luca's RD12 fitted values, or default DX11 value, in this order
    case PHOTONIC_NOISE:
      if (Param->do_photonic_noise != 1) {
        wn_level = Param->do_photonic_noise;
      } else {
        wn_level = GetRD12PhotonicWhiteNoise( brimo->pixname);
        if (wn_level == 1.0) {
          wn_level = brimo->photonic_whitenoise;
        }
      }
      break;
    case ELECTRONIC_NOISE:
      if (Param->do_electronic_noise != 1) {
        wn_level = Param->do_electronic_noise;
      } else {
        wn_level = GetRD12ElectronicWhiteNoise( brimo->pixname);
        if (wn_level == 1.0) {
          wn_level = brimo->electronic_whitenoise;
        }
      }
      break;
    default:
      assert( !"Unknown <noise_type> in add_oof_and_white_noise.c");
  }

  int BeginRing = stimPar->BeginRing;
  int EndRing   = stimPar->EndRing;
  // find begin ring to cover left sample margin
  while (BEGINRINGINDEX( BeginRing) > stimPar->signal_first_sample_number) {
    BeginRing -= 1;
  }
  // find end ring to cover right sample margin
  while (ENDRINGINDEX( EndRing) < stimPar->signal_first_sample_number + stimPar->total_signal_length) {
    EndRing += 1;
  }

  if (noise_type == PHOTONIC_NOISE) {
    STIM_TRACE( 1, "adding photonic white noise (%g) to rings %d-%d", wn_level, BeginRing, EndRing);
  } else {
    STIM_TRACE( 1, "adding electronic white noise (%g) to rings %d-%d", wn_level, BeginRing, EndRing);
  }

  // loop on rings (use a noise seed per ring for realisation reproducibility)
  #pragma omp parallel for
  for (int ring = BeginRing; ring <= EndRing; ring++) {
    long ring_length  = ENDRINGINDEX( ring) - BEGINRINGINDEX( ring) + 1;
    long first_sample = BEGINRINGINDEX( ring) - stimPar->signal_first_sample_number;
    float noise;

    // add white noise
    for (long i = 0; i < ring_length; i++) {
      // always get noise value, even for samples outside of signal array
      // to always have the same noise realisation per ring
      noise = whitenoise( stimPar, ring) * wn_level;
      if ((first_sample + i >= 0) && (first_sample + i < stimPar->total_signal_length)) {
        signal[first_sample + i] += noise;
      }
    }
  }

  return( 0);
}


////////////////////////////////////////////////////////////////////////////////

float GetRD12PhotonicWhiteNoise( char* pixname) {
  // 1.0 means read the default value form IMO/brimo
  if (!strcmp( pixname, "100-1a")) return( 0.00034571);
  if (!strcmp( pixname, "100-1b")) return( 0.00028840);
  if (!strcmp( pixname, "100-2a")) return( 0.00039086);
  if (!strcmp( pixname, "100-2b")) return( 0.00037242);
  if (!strcmp( pixname, "100-3a")) return( 0.00045172);
  if (!strcmp( pixname, "100-3b")) return( 0.00049599);
  if (!strcmp( pixname, "100-4a")) return( 0.00049484);
  if (!strcmp( pixname, "100-4b")) return( 0.00049415);
  if (!strcmp( pixname, "143-1a")) return( 0.00041450);
  if (!strcmp( pixname, "143-1b")) return( 0.00037917);
  if (!strcmp( pixname, "143-2a")) return( 0.00039217);
  if (!strcmp( pixname, "143-2b")) return( 0.00035802);
  if (!strcmp( pixname, "143-3a")) return( 0.00038717);
  if (!strcmp( pixname, "143-3b")) return( 0.00031075);
  if (!strcmp( pixname, "143-4a")) return( 0.00040356);
  if (!strcmp( pixname, "143-4b")) return( 0.00038379);
  if (!strcmp( pixname, "143-5" )) return( 0.00028218);
  if (!strcmp( pixname, "143-6" )) return( 0.00019734);
  if (!strcmp( pixname, "143-7" )) return( 0.00031734);
  if (!strcmp( pixname, "217-1" )) return( 0.00033870);
  if (!strcmp( pixname, "217-2" )) return( 0.00045801);
  if (!strcmp( pixname, "217-3" )) return( 0.00041398);
  if (!strcmp( pixname, "217-4" )) return( 0.00056304);
  if (!strcmp( pixname, "217-5a")) return( 0.00066226);
  if (!strcmp( pixname, "217-5b")) return( 0.00062616);
  if (!strcmp( pixname, "217-6a")) return( 0.00062089);
  if (!strcmp( pixname, "217-6b")) return( 0.00049907);
  if (!strcmp( pixname, "217-7a")) return( 0.00068609);
  if (!strcmp( pixname, "217-7b")) return( 0.00065394);
  if (!strcmp( pixname, "217-8a")) return( 0.00054159);
  if (!strcmp( pixname, "217-8b")) return( 0.00064071);
  if (!strcmp( pixname, "353-1" )) return( 0.00238147);
  if (!strcmp( pixname, "353-2" )) return( 0.00257177);
  if (!strcmp( pixname, "353-3a")) return( 0.00188733);
  if (!strcmp( pixname, "353-3b")) return( 0.00254130);
  if (!strcmp( pixname, "353-4a")) return( 0.00282905);
  if (!strcmp( pixname, "353-4b")) return( 0.00266302);
  if (!strcmp( pixname, "353-5a")) return( 0.00221111);
  if (!strcmp( pixname, "353-5b")) return( 0.00220439);
  if (!strcmp( pixname, "353-6a")) return( 0.00220347);
  if (!strcmp( pixname, "353-6b")) return( 0.00180626);
  if (!strcmp( pixname, "353-7" )) return( 0.00143505);
  if (!strcmp( pixname, "353-8" )) return( 0.00190904);

/*
import piolib, numpy, healpy, LSCtools
pixlist = ["545-1", "545-2", "545-4", "857-1", "857-2", "857-3", "857-4"]
for pixname in pixlist:
  specname = "/data/dmc/MISS03/DATA/detnoise_simu/%s_MeanSpectrumFit_deconv_byRing_v64" % pixname
  white = piolib.read( specname)[0]
  gain, comment = piolib.ReadKeywordObject("GAIN", "PIODOUBLE", specname)
  print '  if (!strcmp( pixname, "%s" )) return( %f);' % (pixname, white * gain / 2.0)
*/

  // 4 Jan 2017 WN fit by SM
  if (!strcmp( pixname, "545-1" )) return( 0.75309);
  if (!strcmp( pixname, "545-2" )) return( 0.61891);
  if (!strcmp( pixname, "545-4" )) return( 0.73764);
  if (!strcmp( pixname, "857-1" )) return( 1.00856);
  if (!strcmp( pixname, "857-2" )) return( 0.98798);
  if (!strcmp( pixname, "857-3" )) return( 0.95622);
  if (!strcmp( pixname, "857-4" )) return( 1.32123);
  return( 1.0);
}


////////////////////////////////////////////////////////////////////////////////

float GetRD12ElectronicWhiteNoise( char* pixname) {
  // 1.0 means read the default value form IMO/brimo
  if (!strcmp( pixname, "100-1a")) return( 0.00084235);
  if (!strcmp( pixname, "100-1b")) return( 0.00073877);
  if (!strcmp( pixname, "100-2a")) return( 0.00057298);
  if (!strcmp( pixname, "100-2b")) return( 0.00060772);
  if (!strcmp( pixname, "100-3a")) return( 0.00065136);
  if (!strcmp( pixname, "100-3b")) return( 0.00075391);
  if (!strcmp( pixname, "100-4a")) return( 0.00056035);
  if (!strcmp( pixname, "100-4b")) return( 0.00073395);
  if (!strcmp( pixname, "143-1a")) return( 0.00051517);
  if (!strcmp( pixname, "143-1b")) return( 0.00060097);
  if (!strcmp( pixname, "143-2a")) return( 0.00054112);
  if (!strcmp( pixname, "143-2b")) return( 0.00054861);
  if (!strcmp( pixname, "143-3a")) return( 0.00057377);
  if (!strcmp( pixname, "143-3b")) return( 0.00061866);
  if (!strcmp( pixname, "143-4a")) return( 0.00058646);
  if (!strcmp( pixname, "143-4b")) return( 0.00066164);
  if (!strcmp( pixname, "143-5" )) return( 0.00044676);
  if (!strcmp( pixname, "143-6" )) return( 0.00050403);
  if (!strcmp( pixname, "143-7" )) return( 0.00044566);
  if (!strcmp( pixname, "217-1" )) return( 0.00085068);
  if (!strcmp( pixname, "217-2" )) return( 0.00081754);
  if (!strcmp( pixname, "217-3" )) return( 0.00079584);
  if (!strcmp( pixname, "217-4" )) return( 0.00075406);
  if (!strcmp( pixname, "217-5a")) return( 0.00086706);
  if (!strcmp( pixname, "217-5b")) return( 0.00090026);
  if (!strcmp( pixname, "217-6a")) return( 0.00082782);
  if (!strcmp( pixname, "217-6b")) return( 0.00091021);
  if (!strcmp( pixname, "217-7a")) return( 0.00078561);
  if (!strcmp( pixname, "217-7b")) return( 0.00080007);
  if (!strcmp( pixname, "217-8a")) return( 0.00085736);
  if (!strcmp( pixname, "217-8b")) return( 0.00086918);
  if (!strcmp( pixname, "353-1" )) return( 0.00132837);
  if (!strcmp( pixname, "353-2" )) return( 0.00088594);
  if (!strcmp( pixname, "353-3a")) return( 0.00354900);
  if (!strcmp( pixname, "353-3b")) return( 0.00317568);
  if (!strcmp( pixname, "353-4a")) return( 0.00326704);
  if (!strcmp( pixname, "353-4b")) return( 0.00337240);
  if (!strcmp( pixname, "353-5a")) return( 0.00284789);
  if (!strcmp( pixname, "353-5b")) return( 0.00294562);
  if (!strcmp( pixname, "353-6a")) return( 0.00468056);
  if (!strcmp( pixname, "353-6b")) return( 0.00528149);
  if (!strcmp( pixname, "353-7" )) return( 0.00294794);
  if (!strcmp( pixname, "353-8" )) return( 0.00294136);

  // 4 Jan 2017 WN fit by SM
  if (!strcmp( pixname, "545-1" )) return( 0.69662);
  if (!strcmp( pixname, "545-2" )) return( 0.67799);
  if (!strcmp( pixname, "545-4" )) return( 0.71091);
  if (!strcmp( pixname, "857-1" )) return( 0.91262);
  if (!strcmp( pixname, "857-2" )) return( 0.96747);
  if (!strcmp( pixname, "857-3" )) return( 0.95245);
  if (!strcmp( pixname, "857-4" )) return( 1.23404);
  return( 1.0);
}
