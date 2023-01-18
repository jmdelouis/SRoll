#ifndef _ADDOOFANDWHITENOISE_H_
#define _ADDOOFANDWHITENOISE_H_


#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "stim_param.h"

#define PHOTONIC_NOISE   0
#define ELECTRONIC_NOISE 1

int add_oof_and_white_noise( stimParameters *stimPar, PIOFLOAT *signal, int noise_type);


#endif
