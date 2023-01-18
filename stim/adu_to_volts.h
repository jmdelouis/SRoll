#ifndef _DEMODADU_H_
#define _DEMODADU_H_


#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "stim_param.h"


int adu_to_volts( stimParameters *stimPar, PIOFLOAT *signal, PIOFLOAT *signal_3ptdemod);


#endif
