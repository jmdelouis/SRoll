#ifndef _COMPRESS_H_
#define _COMPRESS_H_


#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "stim_param.h"


int compress_decompress( stimParameters *stimPar, PIOFLOAT *signal);

int desire_codec( stimParameters *stimPar, PIOFLOAT *signal);


#endif
