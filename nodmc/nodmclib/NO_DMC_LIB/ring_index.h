
#ifndef _RINGINDEX_H
#define _RINGINDEX_H

#include "no_dmc_piolib_type_def.h"

#define MISSION_BEGIN_RING   (240)
#define MISSION_END_RING   (27005)
#define LAST_RING_NUMBER   (27020)

// return the first sample number of <ring_number>
PIOLONG BEGINRINGINDEX( PIOINT ring_number);

// return the last sample number of <ring_number>
PIOLONG ENDRINGINDEX( PIOINT ring_number);

// return the number of samples in <ring_number>
PIOLONG GetRingLength( PIOINT ring_number);

// return 0 if sample number is even, else 1
PIOINT PARITY( PIOLONG sample_number);

// return 1 if sample number is even, else -1, used for signal modulation / demodulation
PIOINT PARITY_MOD( PIOLONG sample_number);

// array containing the values from MISS03/Sa_HFI_C_Bolo/BEGINRINGINDEX
const PIOLONG _BEGINRINGINDEX[LAST_RING_NUMBER + 2];

#endif
