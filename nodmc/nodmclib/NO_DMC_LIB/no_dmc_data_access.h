/******************************************************************************
 * "no_dmc_data_access.h"
 * 
 * author:  Christan Madsen
 * date:    2014-12-22
 * version: STABLE
 *****************************************************************************/

#ifndef _NO_DMC_DATA_ACCESS_H_
#define _NO_DMC_DATA_ACCESS_H_

#include "no_dmc_piolib_type_def.h"


int noDMC_readObject_PIOBYTE(const char *object_name, PIOLONG offset, PIOLONG nbsample, PIOBYTE *data);
int noDMC_readObject_PIOFLAG(const char *object_name, PIOLONG offset, PIOLONG nbsample, PIOFLAG *data);
int noDMC_readObject_PIOSHORT(const char *object_name, PIOLONG offset, PIOLONG nbsample, PIOSHORT *data);
int noDMC_readObject_PIOINT(const char *object_name, PIOLONG offset, PIOLONG nbsample, PIOINT *data);
int noDMC_readObject_PIOLONG(const char *object_name, PIOLONG offset, PIOLONG nbsample, PIOLONG *data);
int noDMC_readObject_PIOFLOAT(const char *object_name, PIOLONG offset, PIOLONG nbsample, PIOFLOAT *data);
int noDMC_readObject_PIODOUBLE(const char *object_name, PIOLONG offset, PIOLONG nbsample, PIODOUBLE *data);
int noDMC_readObject_PIOSTRING(const char *object_name, PIOLONG offset, PIOLONG nbsample, PIOSTRING *data);
int noDMC_readObject_PIOCOMPLEX(const char *object_name, PIOLONG offset, PIOLONG nbsample, PIOCOMPLEX *data);

int noDMC_readObject_PIOINTorPIOFLOAT(const char *object_name, PIOLONG offset, PIOLONG nbsample, PIOFLOAT *data);

// TOSEE: in fact we can know the type of data since it is stored in metadata: can we use this info for a better (smarter) interface?

/*
  Allow to read flags ("Written") corresponding to the specified object name.

  @param object_name the object name for which we want to retrieve the flag info.
  @param offset the offset from which we want to start the reading.
  @param nbsample the number of flags to be read.
  @param flags the pointer to memory allocated by user for storing read values.
 */
int noDMC_readFlagWritten(const char *object_name, PIOLONG offset, PIOLONG nbsample, PIOFLAG *flags);

#endif
