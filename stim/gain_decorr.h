#ifndef _GAINDECORR_H_
#define _GAINDECORR_H_


#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "stim_param.h"


int gain_decorr( stimParameters *stimPar, PIOFLOAT *signal, PIOFLOAT *glitchtails_toi);


#endif
