#ifndef _STIM_H_
#define _STIM_H_

#include <stdio.h>
#include <stdlib.h>

#include "no_dmc_metadata.h"

#ifndef RINGSIZE
#define RINGSIZE (27664l)
#endif

int init_stim( stimParameters *stimPar,   // pointer to an allocated structure to recieve parameter file, brimo and input data
               PIOSTRING param_filename); // parameter file name

int free_stim( stimParameters *stimPar); // pointer to an allocated structure containing data pointers to free

int stim( stimParameters *stimPar, // stim parameter file and brimo
          int begin_ring,          // first ring number to process
          int rings_to_process,    // number of rings to process
          PIOFLOAT *hpr_out);      // if not NULL, will contain the produced TOI projected to a HPR

#endif
