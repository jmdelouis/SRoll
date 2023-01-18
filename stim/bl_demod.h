#ifndef _BLDEMOD_H_
#define _BLDEMOD_H_


#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "stim_param.h"

// LSMOOTH must be even
#define LSMOOTH (648000l)

int bl_demod( stimParameters *stimPar, PIOFLOAT *signal);


#endif
