#ifndef _CONVOLVE_H_
#define _CONVOLVE_H_


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "stim_param.h"

// get M_PI from <math.h>
#define _USE_MATH_DEFINES

// 524288 = 2^19
#define LDECONV    (524288l)
#define HALFPERIOD (LDECONV/2)

#define SPIN_FREQ  (0.016)

// use in stimPar.conv_CONVOLVE to switch between convolution and deconvolution
#define CONVOLVE   (1)
#define DECONVOLVE (0)

int convolve( stimParameters *stimPar, PIOFLOAT *signal);

#endif
