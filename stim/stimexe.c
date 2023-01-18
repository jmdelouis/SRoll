
#define _XOPEN_SOURCE 500

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <sys/time.h>
#include <time.h>
#include <omp.h>

#include "no_dmc_piolib_type_def.h"

#include "stim_param.h"
#include "stim_parLoader.h"
#include "binary_file_toi.h"
#include "stim_tools.h"

#include "stim.h"


int main( int argc,char *argv[]) {

  int             mpi_rank, mpi_size;
  stimParameters  stimPar;
  stim_parContent *Param;
  PIOFLOAT        *hpr_out = NULL;
  long            vmem, phymem;

  // initialise MPI
  MPI_Init( &argc, &argv);  
  MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size( MPI_COMM_WORLD, &mpi_size);

  // read parameter file and brimo only once with rank 0
  if (mpi_rank == 0) {
    char rpath[0124];
    realpath( argv[0], rpath);
    time_t now = time( NULL);
    fprintf( stderr, "%s: --------------------------\n", __FILE__ );
    fprintf( stderr, "%s: Starting %s\n",                __FILE__, rpath);
    fprintf( stderr, "%s: at %s",                        __FILE__, ctime( &now));
    fprintf( stderr, "%s: with %d MPI ranks and %d OMP threads\n", __FILE__, mpi_size, omp_get_max_threads());
    fprintf( stderr, "%s: --------------------------\n", __FILE__ );
  }

  // initialise stim parameters
  init_stim( &stimPar, argv[1]);
  Param = &stimPar.Param;
  Param->stay_in_memory = 1;

  // load balancing per sample count, starting from the end where very long rings are
  int BeginRing[mpi_size];
  int EndRing[mpi_size];
  for (int irank = mpi_size-1; irank >= 0 ; irank--) {
    if (irank == mpi_size-1) {
      EndRing[irank] = Param->EndRing;
    } else {
      EndRing[irank] = BeginRing[irank+1] - 1;
    }
    if (irank == 0) {
      BeginRing[irank] = Param->BeginRing;
    } else {
      long samples_to_process = (ENDRINGINDEX( EndRing[irank]) - BEGINRINGINDEX( Param->BeginRing)) / (irank+1);
      int temp_beg_ring = EndRing[irank] - 1;
      while (ENDRINGINDEX( EndRing[irank]) - BEGINRINGINDEX( temp_beg_ring) < samples_to_process) {
        temp_beg_ring--;
      }
      BeginRing[irank] = temp_beg_ring+1;
    }
  }

#if 1
  if (mpi_rank == 0) {
    int max_samples = 0;
    for (int irank = 0; irank < mpi_size; irank++) {
      printf( "  rank#%d/%d begin=%d end=%d, (%d rings, %de6 samples)\n",
              irank, mpi_size, BeginRing[irank], EndRing[irank],
              EndRing[irank] - BeginRing[irank] + 1,
              (int)((ENDRINGINDEX( EndRing[irank]) - BEGINRINGINDEX( BeginRing[irank])) / 1e6));
      if (max_samples < ENDRINGINDEX( EndRing[irank]) - BEGINRINGINDEX( BeginRing[irank])) {
        max_samples = ENDRINGINDEX( EndRing[irank]) - BEGINRINGINDEX( BeginRing[irank]);
      }
    }
    printf( "max sample count per MPI rank = %.2fe6\n", max_samples / 1e6);
  }
#endif

  int rank_BeginRing = BeginRing[mpi_rank];
  int rank_EndRing   = EndRing[mpi_rank];

  // overwrite rings_per_read parameter, read all needed rings by the rank
  Param->rings_per_read = rank_EndRing - rank_BeginRing + 1;

  hpr_out = (PIOFLOAT *) malloc( Param->rings_per_read * RINGSIZE * sizeof(PIOFLOAT));
  assert( hpr_out != NULL);

  // call simulation for this mpi rank
  for (int iter = 0; iter < Param->iterations; iter++) { 
    for (long i=0; i<Param->rings_per_read * RINGSIZE; i++) {
      hpr_out[i] = 0.0;
    }
    stim( &stimPar, rank_BeginRing, Param->rings_per_read, hpr_out);
    if (mpi_rank == 0) {
      stim_GetProcMem( &vmem, &phymem);
      printf( "%s **** %s after realisation %d, VMem=%dMB PhyMem=%dMB\n", stimPar.msg_prefix, __FILE__, Param->random_seed, (int)(vmem/1e6), (int)(phymem/1e6));
    }
    Param->random_seed++;
    MPI_Barrier( MPI_COMM_WORLD);
  }

  if (hpr_out != NULL) {
    free( hpr_out);
  }

  if (mpi_rank == 0) {
    time_t now = time( NULL);
    fprintf( stderr, "\n");
    fprintf( stderr, "%s: --------------------------\n", __FILE__ );
    fprintf( stderr, "%s: Finished sucessfully at %s",   __FILE__, ctime( &now));
    fprintf( stderr, "%s: --------------------------\n", __FILE__ );
  }
  
  MPI_Finalize();
  exit( 0);
}
