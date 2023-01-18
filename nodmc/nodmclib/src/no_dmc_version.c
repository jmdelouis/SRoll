/******************************************************************************
 * "no_dmc_version.c"
 * 
 * author:  Christan Madsen
 * date:    2015-04-02 (initial)
 * version: STABLE
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include "no_dmc_version.h"


/*
 * Allow to get the no dmc library version as a string.
 * Note that it is the caller responsability to free the memory returned by this funtion!
 *
 * @return the no dmc library version as a string, NULL otherwise.
 */
char *noDMC_getVersion() {
  const int MAX_SIZE_VER = 50;
  
  char *ver = malloc(MAX_SIZE_VER * sizeof(char));

  if (ver == NULL) {
    perror("Error");
    return NULL;
  }

  snprintf(ver, MAX_SIZE_VER, "%d.%d", __NODMCLIB_MAJOR__, __NODMCLIB_MINOR__);

  return ver;
}
