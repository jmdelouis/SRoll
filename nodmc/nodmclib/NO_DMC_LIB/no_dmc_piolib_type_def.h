/******************************************************************************
 * "no_dmc_piolib_type_def.h"
 *
 *
 * This file is a partial copy of type definitions from:
 * http://cvs.planck.fr/cvs/Level2/Lib_pkg/HL2_DMC/PioLib/HL2_PIOLIB/PIOLib.h?revision=1.175
 *****************************************************************************/

#ifndef ISDMC // don't duplicate Lib_pkg/HL2_DMC/PioLib/HL2_PIOLIB/PIOLib.h definitions if compiled inside the DMC

#ifndef _NO_DMC_PIOLIB_TYPE_DEF_H_
#define _NO_DMC_PIOLIB_TYPE_DEF_H_

#include <stdint.h>
#include <inttypes.h>
#include <limits.h>

#define PIOSTRINGMAXLEN (256)

typedef unsigned char PIOBYTE;

typedef unsigned char PIOFLAG; /* 0 or 1 only */

typedef short PIOSHORT;

typedef int32_t PIOINT;
#define PIOINT_FMT "%"PRId32
#define PIOINT_MIN INT32_MIN
#define PIOINT_MAX INT32_MAX

typedef int64_t PIOLONG;
#define PIOLONG_FMT "%"PRId64
#define PIOLONG_MIN INT64_MIN
#define PIOLONG_MAX INT64_MAX

typedef float PIOFLOAT;

typedef double PIODOUBLE;

typedef char PIOSTRING[PIOSTRINGMAXLEN];

typedef struct {
  float real;
  float imaginary;
} PIOCOMPLEX;

#endif // ifndef _NO_DMC_PIOLIB_TYPE_DEF_H_

#endif // ifndef ISDMC
