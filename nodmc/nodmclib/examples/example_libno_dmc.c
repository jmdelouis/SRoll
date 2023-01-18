/******************************************************************************
 * "example_libno_dmc.c"
 *
 * This program aims to be a simple example for using the NO_DMC_LIB.
 * It demonstrate the most usual function for a specified object:
 * - how to access the metadata
 * - how to retrieve data from object
 * - how to retrieve flag (written) information from object
 *
 * Author: C. Madsen
 * Date: 2015-03-06
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include "no_dmc_metadata.h"
#include "no_dmc_data_access.h"
#include "no_dmc_piolib_type_def.h"
#include "no_dmc_util.h"


#define TEST_OK 0
#define TEST_ERROR -1


//-----------------------------------------------------------------------------
// Check 'res' value: if non zero print error message and exit, unless it was expected!
//-----------------------------------------------------------------------------
void check_res(int res, int expected_res) {
  bool success = (res == expected_res);

  fprintf(stderr, ">> Result = %d   (expected %d)   %s\n\n", res, expected_res, success?"SUCCESS":"FAILED");

  if (!success) {
    fprintf(stderr, "TEST FAILED detected in last call!\n");
    exit(res);
  }
}


//-----------------------------------------------------------------------------
// This function demonstrate how to read metadata.
//-----------------------------------------------------------------------------
int read_no_dmc_metadata(char *obj) {
  metadata m_test; // the metadata used to store info

  // Call metadata read function for test obj
  int ret = getMetadataFor(obj, &m_test);
  if (ret != 0) {
    return ret;
  }

  // Print retrieved info
  fprintf(stderr, "METADATA read for DMC object '%s'\n", obj);
  fprintf(stderr, "  noDMCmetadata_date    = '%s'\n", m_test.noDMCmetadata_date);
  fprintf(stderr, "  noDMCmetadata_version = '%s'\n", m_test.noDMCmetadata_version);
  fprintf(stderr, "  PIOtype               = '%s'\n", m_test.PIOtype);
  fprintf(stderr, "  Datatype              = '%s'\n", m_test.Datatype);
  fprintf(stderr, "  BeginIndex            = "PIOLONG_FMT"\n", m_test.BeginIndex);
  fprintf(stderr, "  EndIndex              = "PIOLONG_FMT"\n", m_test.EndIndex);
  fprintf(stderr, "  BeginRing             = %d\n", m_test.BeginRing);
  fprintf(stderr, "  EndRing               = %d\n", m_test.EndRing);
  fprintf(stderr, "  Author                = '%s'\n", m_test.Author);
  fprintf(stderr, "  Date                  = '%s'\n", m_test.Date);
  fprintf(stderr, "  Keyword_list          = '%s'\n", m_test.Keyword_list);
  fprintf(stderr, "  Backendname           = '%s'\n", m_test.Backendname);
  fprintf(stderr, "  Iooffset              = "PIOLONG_FMT"\n", m_test.Iooffset);
  fprintf(stderr, "  Flagchunksize         = "PIOLONG_FMT"\n", m_test.Flagchunksize);

  // Free mem
  freeMetadata(&m_test);

  return 0;
}


//-----------------------------------------------------------------------------
// Demonstrate capacity to read data for PIODOUBLE type.
//-----------------------------------------------------------------------------
int read_data_as_PIODOUBLE(char *obj, PIOLONG offset, PIOLONG nbSample) {
  int ret;
  void *data;
  PIOSHORT sizeofType = sizeof(PIODOUBLE);

  fprintf(stderr, "  READ data from offset = "PIOLONG_FMT" for nbSample = "PIOLONG_FMT"\n", offset, nbSample);

  // Allocate memory for data to be read
  data = calloc(nbSample, sizeofType);
  if (data == NULL) {
    fprintf(stderr, "Error: unable to allocate memory for data\n");
    return -1;
  }

  // Retrieve data using the new library
  // Select the appropriate function according to data type to be read
  ret = noDMC_readObject_PIODOUBLE(obj, offset, nbSample, data);

  // Check return
  if (ret != 0) {
    free(data);
    return ret;
  }

  // Dump retrieve data
  PIOLONG i;
  for (i = 0; i < nbSample; i++) {
    fprintf(stderr, "  val #"PIOLONG_FMT" = %g\n", i, ((PIODOUBLE *)data)[i]);
  }

  // free mem
  free(data);

  return 0;
}


//-----------------------------------------------------------------------------
// Demonstrate capacity to read Flag Written for a specified object.
//-----------------------------------------------------------------------------
int read_Flag(char *obj, PIOLONG offset, PIOLONG nbSample) {
  int ret;
  PIOFLAG *flags;

  fprintf(stderr, "  READ FLAG data from offset = "PIOLONG_FMT" for nbSample = "PIOLONG_FMT"\n", offset, nbSample);

  // Allocate memory for data to be read
  flags = calloc(nbSample, sizeof(PIOFLAG));
  if (flags == NULL) {
    fprintf(stderr, "Error: unable to allocate memory for flags\n");
    return -1;
  } 

  // Retrieve flags using the new library
  ret = noDMC_readFlagWritten(obj, offset, nbSample, flags);

  // Check return
  if (ret != 0) {
    free(flags);
    return ret;
  }

  // For debug purpose: dump retrieve flags with the no_dmc_lib
  PIOLONG i;
  for (i = 0; i < nbSample; i++) {
    fprintf(stderr, "  flag #"PIOLONG_FMT" = '%s'\n", i, flags[i] == 0 ? "False":"True");
  }

  // free mem
  free(flags);

  return 0;
}


//*****************************************************************************
// MAIN
//*****************************************************************************
int main() {
  int res;

  fprintf(stderr, "**************************************************\n");
  fprintf(stderr, "** This is the example program for 'no dmc lib' **\n");
  fprintf(stderr, "**************************************************\n");

  char *obj = "/data/dmc/MISS03/DATA/RingInfo/SpinPhase"; // The DMC object used for the example
  PIOLONG offset = 10; // Offset is the index from which we start to read
  PIOLONG nbSample = 12; // Number of samples to be read

  fprintf(stderr, "\nData used for this example:\n");
  fprintf(stderr, "    >> DMC Object fullpath = %s\n", obj);
  fprintf(stderr, "    >> offset = "PIOLONG_FMT"\n", offset);
  fprintf(stderr, "    >> number of samples to be read = "PIOLONG_FMT"\n", nbSample);
  fprintf(stderr, "\n\n");

  // read_no_dmc_metadata()
  fprintf(stderr, "===============================================\n");
  fprintf(stderr, "read_no_dmc_metadata()\n");
  fprintf(stderr, "-----------------------------------------------\n");
  res = read_no_dmc_metadata(obj);
  fprintf(stderr, "-----------------------------------------------\n");
  check_res(res, TEST_OK); // Always expect to succeed!

  // read_data_as_PIODOUBLE()
  fprintf(stderr, "===============================================\n");
  fprintf(stderr, "read_data_as_PIODOUBLE()\n");
  fprintf(stderr, "-----------------------------------------------\n");
  res = read_data_as_PIODOUBLE(obj, offset, nbSample);
  fprintf(stderr, "-----------------------------------------------\n");
  check_res(res, TEST_OK);

  // read_Flag()
  fprintf(stderr, "===============================================\n");
  fprintf(stderr, "read_Flag()\n");
  fprintf(stderr, "-----------------------------------------------\n");
  res = read_Flag(obj, offset, nbSample);
  fprintf(stderr, "-----------------------------------------------\n");
  check_res(res, TEST_OK);

  fprintf(stderr, "\n\nALL SUCCESS !\n");

  return 0;
}
