###############################################################################
# "cno_dmc_data_access.pxd"
# This is the no_dmc_lib wrapper for python.
#
# Author: C. Madsen
# Date: 2015-01-26
###############################################################################


cdef extern from "no_dmc_data_access.h":
  bint noDMC_readObject_PIOBYTE(char *object_name, long offset, long nbsample, char *data)
  bint noDMC_readObject_PIOFLAG(char *object_name, long offset, long nbsample, unsigned char *data)
  bint noDMC_readObject_PIOSHORT(char *object_name, long offset, long nbsample, short *data)
  bint noDMC_readObject_PIOINT(char *object_name, long offset, long nbsample, int *data)
  bint noDMC_readObject_PIOLONG(char *object_name, long offset, long nbsample, long *data)
  bint noDMC_readObject_PIOFLOAT(char *object_name, long offset, long nbsample, float *data)
  bint noDMC_readObject_PIODOUBLE(char *object_name, long offset, long nbsample, double *data)
  bint noDMC_readObject_PIOSTRING(char *object_name, long offset, long nbsample, char *data[256])

  bint noDMC_readFlagWritten(char *object_name, long offset, long nbsample, unsigned char *flags)


cdef extern from "no_dmc_metadata.h":
  ctypedef struct metadata:
    char *PIOtype
    char *Datatype
    long BeginIndex
    long EndIndex
    int BeginRing
    int EndRing
    char *Author
    char *Date

  bint getMetadataFor(char *object_name, metadata *metadata_p)
