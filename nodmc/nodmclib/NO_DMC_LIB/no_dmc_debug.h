/******************************************************************************
 * "no_dmc_debug.h"
 * 
 * author:  Christan Madsen
 * date:    2015-01-30
 * version: 1.0
 *****************************************************************************/

#ifndef _NO_DMC_DEBUG_H_
#define _NO_DMC_DEBUG_H_


#ifndef DEBUG
#define DEBUG 0
#else
#include <stdio.h>
#endif

#define debug_print(fmt, ...) \
  do { if (DEBUG) fprintf(stderr, "[DBG] " fmt, ##__VA_ARGS__); } while (0)
/*  do { if (DEBUG) fprintf(stderr, "[DEBUG %s:%d:%s()] " fmt, __FILE__, \
      __LINE__, __func__, ##__VA_ARGS__); } while (0)
*/

/* same as previous but without the prefix [DBG] */
#define debug_printNP(fmt, ...) \
  do { if (DEBUG) fprintf(stderr, fmt, ##__VA_ARGS__); } while (0)

#endif
