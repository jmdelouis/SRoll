
#ifndef _BINARY_FILE_TOI_H_
#define _BINARY_FILE_TOI_H_

#include "no_dmc_data_access.h"

typedef struct {
  PIOLONG first_sample;
  PIOLONG last_sample;
} written_chunk_struct;

// read BinaryFileToi data, returns 0 if ok, failed assert if nok. data is malloced and freed by the caller
int BFT_readObject( const char *object_name,
                    const char *data_type,
                    long first_sample,
                    long nbsample,
                    void *data);

// write BinaryFileToi data, returns nbsample (>=0) if ok, failed assert if nok
int BFT_writeObject( const char *object_name,
                     const char *data_type,
                     long first_sample,
                     long nbsample,
                     void *data);

#endif
