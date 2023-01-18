#ifndef _ADDOOFNOISE_H_
#define _ADDOOFNOISE_H_


#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "stim_param.h"

#define NOISEKERNEL (10800.)

int add_oof_noise( stimParameters *stimPar, PIOFLOAT *signal);
void add_oof_noise_setdata(double the_loc_amp);


#endif
