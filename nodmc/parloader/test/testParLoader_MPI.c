/******************************************************************************
 * "testParLoader_MPI.c"
 *
 * This is the test program for the 'parLoader' facility using MPI.
 *
 * Author: C. Madsen
 * Date: 2015-06-09
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <mpi.h>


#include "test_parLoader.h"


//-----------------------------------------------------------------------------

// Check 'res' value and if non zero print error message and exit, unless it was expected!
void check_res(int res, int expected_res, bool verbose) {
  bool success = (res == expected_res);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  if (verbose) {
    fprintf(stderr, "rank#%03d   Result: %s\n", rank, success?"SUCCESS":"FAILED");
  } else {
    fprintf(stderr, " --> %s\n", success?"SUCCESS":"FAILED");
  }

  if (!success) {
    fprintf(stderr, "rank#%03d TEST FAILED detected in last call!\n", rank);
    exit(res);
  }
}

//-----------------------------------------------------------------------------


//*****************************************************************************
// MAIN
//*****************************************************************************
int main(int argc, char *argv[]) {
  int res;
  int rank,size;

  char filename[] = "samples/Param_OK_143_RD12.txt";
  test_parContent par;

  // Start MPI
  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);


  if (rank == 0) {
    fprintf(stderr, "******************************************************************\n");
    fprintf(stderr, "** This is the test program for 'parLoader' facility using MPI. **\n");
    fprintf(stderr, "******************************************************************\n");


    fprintf(stderr, "===============================================\n");
    fprintf(stderr, "Testing param loading from: \"%s\"\n", filename);
    fprintf(stderr, "-----------------------------------------------\n");

    fprintf(stderr, "Total rank = %d\n", size);
    fprintf(stderr, "-----------------------------------------------\n\n");
  }


  MPI_Barrier(MPI_COMM_WORLD); // Synchro for improve display stderr...

  fprintf(stderr, "rank#%03d PHASE 1 - Loading param file...\n", rank);

  res = test_readParam(&par, filename);
  check_res(res, 0, true); // Check success!


  fprintf(stderr, "rank#%03d PHASE 2 - Other tests...\n", rank);

  // Print some value
  fprintf(stderr, "rank#%03d   par.Theo_Dust_I = \"%s\"\n", rank, par.Theo_Dust_I); // PIOSTRING

  // Test: Allocate memory
  char *tmpstr = malloc(1000000*sizeof(char));
  if (tmpstr == NULL) {
    fprintf(stderr, "rank#%03d   . Alloc memory: ERROR\n", rank);
    return 1;
  } else {
    fprintf(stderr, "rank#%03d   . Alloc memory: OK\n", rank);
  }
  
  fprintf(stderr, "rank#%03d   . Printing some values...\n", rank);
  for (PIOLONG i=0; i < par.n_NEP; i++) {
    fprintf(stderr, "rank#%03d     par.NEP["PIOLONG_FMT"]      = %.12g\n", rank, i, par.NEP[i]); // PIODOUBLE
  }

  // Just for synchro before end!
  MPI_Barrier(MPI_COMM_WORLD);

  fprintf(stderr, "rank#%03d   . Printing some MORE values (after barrier)...\n", rank);
  for (PIOLONG i=0; i < par.n_fsl; i++) {
    fprintf(stderr, "rank#%03d     par.fsl["PIOLONG_FMT"]      = \"%s\"\n", rank, i, par.fsl[i]); // PIOSTRING
  }

  // Free memory after barrier to ensure worst condition (ie; all ranks have their memory allocated!)
  free(tmpstr);
  tmpstr = NULL;

  fprintf(stderr, "rank#%03d   . Printing AGAIN some MORE values (after free)...\n", rank);
  for (PIOLONG i=0; i < par.n_bolomask; i++) {
    fprintf(stderr, "rank#%03d     par.bolomask["PIOLONG_FMT"]      = "PIOINT_FMT"\n", rank, i, par.bolomask[i]); // PIOLONG
  }

  if (rank == 0) {
    fprintf(stderr, "\n\nALL SUCCESS !\n");
  }

  MPI_Finalize();

  return 0;
}
