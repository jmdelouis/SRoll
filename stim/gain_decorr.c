
#define _XOPEN_SOURCE

#include <assert.h>
#include <math.h>

#include "no_dmc_data_access.h"

#include "gain_decorr.h"
#include "stim_tools.h"


//out_sig = ((in_sig - glitchtails) * g0 * (1.0 + (in_sig - glitchtails) / v0)) - coeffT90 * thermal_baseline

int gain_decorr( stimParameters *stimPar, PIOFLOAT *signal, PIOFLOAT *glitchtails_toi) {

  stim_parContent *Param = &stimPar->Param;
  stim_Data       *data  = &stimPar->data;
  BRIMO           *brimo = &stimPar->brimo;

  float g0 = brimo->g0;
  float v0 = brimo->v0;
  float coeffT90 = brimo->coeffT90;

  STIM_TRACE( 1, "converting signal from Volts to Watts and subtracting thermal baseline (x%g) for rings %d-%d", coeffT90, stimPar->BeginRing, stimPar->EndRing);

  if (Param->use_bolometer_nonlinearity == 0) {
    v0 = 1e20;
  }


  if (strcmp( Param->thermal_baseline, "0") != 0) {
    if (data->thermal_baseline == NULL) {
      data->thermal_baseline = malloc( stimPar->total_signal_length * sizeof( PIOFLOAT));
      assert( data->thermal_baseline != NULL);
      assert( noDMC_readObject_PIOFLOAT( Param->thermal_baseline,
                                         stimPar->signal_first_sample_number,
                                         stimPar->total_signal_length,
                                         data->thermal_baseline) >= 0);
    }
  }

  #pragma omp parallel for
  for (long i = 0; i < stimPar->total_signal_length; i++) {
    //subtract glitch tails TOI
    if (glitchtails_toi != NULL) {
      signal[i] -= glitchtails_toi[i];
    }
    // convert from Volts to Watts
    signal[i] = signal[i] * g0 * (1.0 + signal[i] / v0);
    // subtract thermal baseline converted to Watts
    if (data->thermal_baseline != NULL) {
      signal[i] -= data->thermal_baseline[i] * coeffT90;
    }
    assert( !isnan( signal[i]));
  }

  if (!Param->stay_in_memory) {
    if (data->thermal_baseline != NULL) {
      free( data->thermal_baseline);
      data->thermal_baseline = NULL;
    }
  }
  
  return( 0);
}
