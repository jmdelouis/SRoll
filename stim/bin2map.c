
#define _XOPEN_SOURCE 500

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <mpi.h>
#include <math.h>
#include <sys/time.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/param.h>
#include <fcntl.h>

#include "chealpix.h"

#include "no_dmc_metadata.h"
#include "no_dmc_data_access.h"
#include "no_dmc_piolib_type_def.h"
#include "no_dmc_util.h"
#include "no_dmc_debug.h"
#include "no_dmc_version.h"
#include "ring_index.h"
#include "brimo.h"
#include "stim_tools.h"


#define HPRSIZE (27664l)
#define NSIDE   (2048l)
#define MAPSIZE (12*NSIDE*NSIDE)
#define HPRIDX  "/e2e_common_TOI/%s_HPRIDX_ABER_TotalFlag_dx11"
#define HITHPR  "/PBR_JMD/%s_REP6_hit"
#define PTGPHI  "/PBR_JMD/%s_REP6_ptg"
#define PTGTHE  "/PBR_JMD/%s_REP6_ptg_TUPLE_1"
#define DIPHPR  "/HPRREPSM_PBR/%s_v64_MAY16_diporb"
//#define DIPHPR  "/PBR_JMD/%s_dipHFI17_quad_hprbin"
#define BADRING "/calROIs/%s_bad_rings_and_offsets_rd12"
#define STRLEN  (260l)
#define MAXREAD (1e8)


typedef struct {
  char toiname[STRLEN];
  int  flag_toiname;
  char hprname[STRLEN];
  int  flag_hprname;
  char sigmap[STRLEN];
  int  flag_sigmap;
  char hitmap[STRLEN];
  int  flag_hitmap;
  int  beginring;
  int  endring;
  int  subdiphpr;
  int  badrings;
} bin2map_parContent;

int mpi_rank=0;
int mpi_size=0;

char *survey_name[] = {"_full", "_hm1", "_hm2", "_odd", "_even"};
int survey_beg[] = {240, 240, 13145, 240, 240};
int survey_end[] = {26050, 13144, 26050, 26050, 26050};
int n_surveys = 5;


////////////////////////////////////////////////////////////////////////////////

void usage( char *msg) {

  if (mpi_rank == 0) {
    fprintf( stderr,"%s", msg);
    fprintf( stderr, "\n\n");
    fprintf( stderr, "bin2map: just bin (flat average per pixel) a signal TOI or HPR in a nside=2048 map using RD12 pointing and flagging\n");
    fprintf( stderr, "         you must provide one input (TOI or HPR) and one output (HPR or MAP)\n");
    fprintf( stderr, "         if no beginring and endring are given, maps for full mission, half-missions and odd and even rings will be produced\n");
    fprintf( stderr, "usage: bin2map --toiname=<toiname> --hprname=<hprname> --sigmap=<sigmap> --hitmap=<hitmap> [--subdiphpr] [--beginring=<beginring>] [--endring=<endring>]\n");
    fprintf( stderr, "       - <toiname>: name of input float32 TOI in HFI DMC or flat bianry format\n");
    fprintf( stderr, "       - <hprname>: name of input float32 HPR in HFI DMC or flat bianry format\n");
    fprintf( stderr, "       - --subdiphpr: if present, total dipole HPR will be subtracted from input signal HPR\n");
    fprintf( stderr, "       - --badrings: if present, bad rings from '{pixname}_bad_rings_and_offsets_rd12' won't be projected (but processed in the HPR)\n");
    fprintf( stderr, "       - <sigmap>: name of output float32 map containing binned signal, mandatory\n");
    fprintf( stderr, "       - <hitmap>: name of output int32 map containing hitcount, optional\n");
    fprintf( stderr, "       - <beginring>: first ring number to process\n");
    fprintf( stderr, "       - <endring>: last ring number to process\n");
    fprintf( stderr, "example: mpirun -np 8 ./bin2map --toiname=/redtruck/SimuData/e2e_common_TOI/143-1a_dustapr17_2019_stimTOI_000.float32.bin --hprname=/redtruck/SimuData/e2e_common_TOI/143-1a_dustapr17_2019_000_hprbin --sigmap=/redtruck/SimuData/e2e_common_TOI/143-1a_dustapr17_2019_MAP --badrings --beginring=240 --endring=26050\n");
    fprintf( stderr, "\n");
  }
  MPI_Finalize();
  exit( -1);
}


////////////////////////////////////////////////////////////////////////////////

void readParameters( bin2map_parContent *param, int argc, char *argv[]) {

  param->flag_toiname = 0;
  param->flag_hprname = 0;
  param->flag_sigmap = 0;
  param->flag_hitmap = 0;
  param->beginring = 0;
  param->endring = 0;
  param->subdiphpr = 0;
  param->badrings = 0;

  for (int argi=1; argi<argc; argi++) {
    if (startswith( argv[argi], "--toiname=")) {
      param->flag_toiname = 1;
      strncpy( param->toiname, argv[argi] + strlen( "--toiname="), STRLEN);
      if (mpi_rank==0) fprintf( stderr, "parameter toiname=%s\n", param->toiname);
    }
    else if (startswith( argv[argi], "--hprname=")) {
      param->flag_hprname = 1;
      strncpy( param->hprname, argv[argi] + strlen( "--hprname="), STRLEN);
      if (mpi_rank==0) fprintf( stderr, "parameter hprname=%s\n", param->hprname);
    }
    else if (startswith( argv[argi], "--sigmap=")) {
      param->flag_sigmap = 1;
      strncpy( param->sigmap, argv[argi] + strlen( "--sigmap="), STRLEN);
      if (mpi_rank==0) fprintf( stderr, "parameter sigmap=%s\n", param->sigmap);
    }
    else if (startswith( argv[argi], "--hitmap=")) {
      param->flag_hitmap = 1;
      strncpy( param->hitmap, argv[argi] + strlen( "--hitmap="), STRLEN);
      if (mpi_rank==0) fprintf( stderr, "parameter hitmap=%s\n", param->hitmap);
    }
    else if (startswith( argv[argi], "--beginring=")) {
      param->beginring = strtol( argv[argi] + strlen( "--beginring="), NULL, 0);
      if (mpi_rank==0) fprintf( stderr, "parameter beginring=%d\n", param->beginring);
    }
    else if (startswith( argv[argi], "--endring=")) {
      param->endring = strtol( argv[argi] + strlen( "--endring="), NULL, 0);
      if (mpi_rank==0) fprintf( stderr, "parameter endring=%d\n", param->endring);
    }
    else if (strcmp( argv[argi], "--subdiphpr") == 0) {
      param->subdiphpr = 1;
      if (mpi_rank==0) fprintf( stderr, "parameter subdiphpr=1\n");
    }
    else if (strcmp( argv[argi], "--badrings") == 0) {
      param->badrings = 1;
      if (mpi_rank==0) fprintf( stderr, "parameter badrings=1\n");
    }
    else {
      fprintf( stderr, "Unknown parameter: %s\n\n", argv[argi]);
      usage( NULL);
      exit( -1);
    }
  }
  if ((param->flag_toiname == 0) && (param->flag_hprname == 0))
    usage( "Error: input <toiname> or <hprname> is mandatory");
  if ((param->flag_hprname == 0) && (param->flag_sigmap == 0))
    usage( "Error: output <hprname> or <sigmap> is mandatory");
  if (!param->beginring != !param->endring)
    usage( "Error: you must give both <beginring> and <endring> or none");
}


////////////////////////////////////////////////////////////////////////////////

int main( int argc, char *argv[]) {

  long i, ipix;
  int  dmc_err;
  DETNAME pixname;
  char tempstr[STRLEN];
  char objname[STRLEN];
  char *dbpath = NULL;

  bin2map_parContent param;

  // initialise MPI
  MPI_Init( &argc, &argv);
  MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size( MPI_COMM_WORLD, &mpi_size);

  if (mpi_rank == 0) {
    char rpath[PATH_MAX];
    fprintf( stderr, "starting %s\n" , realpath( argv[0], rpath));
  }

  // get DBPATH from environment variable
  dbpath = getenv("NODMCDATA");
  if (dbpath == NULL) {
    dbpath = getenv("DMCDATA");
    if (dbpath == NULL) {
      if (mpi_rank == 0) {
        fprintf( stderr, "ERROR: you must have either NODMCDATA or DMCDATA evironment variables set.\n");
        fprintf( stderr, "(source ./srollex_setenv.sh).\n");
      }
      MPI_Finalize();
      exit( -1);
    }
  }

  readParameters( &param, argc, argv);

  // find a pixel name in input object
  char *inputobject = NULL;
  if (param.flag_toiname) {
    inputobject = param.toiname;
  }
  else {
    inputobject = param.hprname;
  }
  strncpy( pixname, get_pixname( inputobject), DETNAMELEN);

  if (param.endring != 0) {
    n_surveys = 1;
    survey_beg[0] = param.beginring;
    survey_end[0] = param.endring;
  }
  else {
    param.beginring = 240;
    param.endring = 26050;
  }

  // get begin and end ring per mpi rank
  int ring_per_rank = (param.endring - param.beginring + 1) / mpi_size;
  int ring_remainder = (param.endring - param.beginring + 1) % mpi_size;
  int begring_offset = ring_remainder;
  int endring_offset = ring_remainder;
  if (mpi_rank < ring_remainder) begring_offset = mpi_rank;
  if (mpi_rank+1 < ring_remainder) endring_offset = mpi_rank+1;
  int rank_begring = param.beginring + mpi_rank * ring_per_rank + begring_offset;
  int rank_endring = param.beginring + (mpi_rank + 1) * ring_per_rank + endring_offset - 1;
  int nrings  = rank_endring - rank_begring + 1;
  fprintf( stderr, "rank %d/%d, rank_begring=%d, rank_endring=%d\n", mpi_rank, mpi_size, rank_begring, rank_endring);

  if (param.subdiphpr) {
    fprintf( stderr, "using dipole HPR: %s\n", DIPHPR);
  }

  // allocate signal and hitcount HPRs memory
  float *sighpr = malloc( nrings * HPRSIZE * sizeof( float));
  assert( sighpr != NULL);
  float *hithpr = malloc( nrings * HPRSIZE * sizeof( float));
  assert( hithpr != NULL);

//
// TOI to HPR
//

  if (param.flag_toiname) {
    float *sigtoi = malloc( MAXREAD * sizeof( float));
    assert( sigtoi != NULL);
    int *pixel_index = malloc( MAXREAD * sizeof( int));
    assert( pixel_index != NULL);
    int toi_begring = rank_begring;

    sprintf( tempstr, HPRIDX, pixname);
    sprintf( objname, "%s%s", dbpath, tempstr);

    while (toi_begring <= rank_endring) {
      int toi_endring = toi_begring;
      long read_first_sample = BEGINRINGINDEX( toi_begring);
      while ( (toi_endring < rank_endring) && (ENDRINGINDEX( toi_endring+1) - read_first_sample < MAXREAD)) {
        toi_endring += 1;
      }
      long read_last_sample = ENDRINGINDEX( toi_endring);
      long read_sample_count = read_last_sample - read_first_sample + 1;
      fprintf( stderr, "rank %d/%d, TOI2HPR rings %d-%d, sample_count=%ld\n", mpi_rank, mpi_size, toi_begring, toi_endring, read_sample_count);
      assert( read_sample_count < MAXREAD);

      dmc_err = noDMC_readObject_PIOFLOAT( param.toiname, read_first_sample, read_sample_count, sigtoi);
      if (dmc_err < 0) {
        fprintf( stderr, "Error %d when reading %s rings %d-%d\n", dmc_err, param.toiname, toi_begring, toi_endring);
        exit( -1);
      }
      dmc_err = noDMC_readObject_PIOINT( objname, read_first_sample, read_sample_count, pixel_index);
      if (dmc_err < 0) {
        fprintf( stderr, "Error %d when reading %s rings %d-%d\n", dmc_err, objname, toi_begring, toi_endring);
        exit( -1);
      }

      // project TOI to HPR
      for (int hpr_ring = toi_begring; hpr_ring <= toi_endring; hpr_ring += 1) {
        long beginringindex = BEGINRINGINDEX( hpr_ring);
        long endringindex   = ENDRINGINDEX(   hpr_ring);
        long hpr_begin      = (hpr_ring - rank_begring) * HPRSIZE;

        for (i = 0; i < HPRSIZE; i++) {
          sighpr[hpr_begin + i] = 0.0;
          hithpr[hpr_begin + i] = 0;
        }

        for (i = beginringindex - read_first_sample; i < endringindex + 1 - read_first_sample; i++) {
          if ((pixel_index[i] >= 0) && (pixel_index[i]) < HPRSIZE) {
            sighpr[hpr_begin + pixel_index[i]] += sigtoi[i];
            hithpr[hpr_begin + pixel_index[i]] += 1;
          }
        }

        for (i = 0; i < HPRSIZE; i++) {
          if (hithpr[hpr_begin + i] > 0) {
            sighpr[hpr_begin + i] /= hithpr[hpr_begin + i];
          }
          if ( isnan( sighpr[hpr_begin + i])) {
            printf("nan at offset %ld in ring %d\n", i, hpr_ring);
            exit(-1);
          }
        }
      }

      toi_begring = toi_endring + 1;
    }
    free( sigtoi);
    free( pixel_index);

    if (param.flag_hprname) {
      // write HPR to disk
      for (i = 0; i < mpi_size; i++) {
        if (mpi_rank == i) {
          writebin( param.hprname, sighpr, (long)nrings * HPRSIZE * sizeof( float),
                                           (long)rank_begring * HPRSIZE * sizeof( float));
#if 1
          // write hitcount HPR to disk for comparison with PBR_JMD/{pixname}_REP6_hit
          sprintf( objname, "%s%s", param.hprname, "_hit");
          writebin( objname, hithpr, (long)nrings * HPRSIZE * sizeof( float),
                                     (long)rank_begring * HPRSIZE * sizeof( float));
#endif
        }
        MPI_Barrier( MPI_COMM_WORLD);
      }
    }
  }

  else {
    // read signal HPR
    dmc_err = noDMC_readObject_PIOFLOAT( param.hprname, rank_begring * HPRSIZE, nrings * HPRSIZE, sighpr);
    if (dmc_err < 0) {
      fprintf( stderr, "Error %d when reading %s\n", dmc_err, param.hprname);
      exit( -1);
    }
    // read hitcount HPR
    sprintf( tempstr, HITHPR, pixname);
    sprintf( objname, "%s%s", dbpath, tempstr);
    dmc_err = noDMC_readObject_PIOFLOAT( objname, rank_begring * HPRSIZE, nrings * HPRSIZE, hithpr);
    if (dmc_err < 0) {
      fprintf( stderr, "Error %d when reading %s\n", dmc_err, objname);
      exit( -1);
    }
  }

//
// HPR to MAP
//

  if (param.flag_sigmap) {
    // allocate signal and hit maps
    double *mapsig = malloc( MAPSIZE * sizeof( double));
    assert( mapsig != NULL);
    int *maphit = malloc( MAPSIZE * sizeof( int));
    assert( maphit != NULL);

    // read phi pointing HPR
    double *ptgphi = malloc( nrings * HPRSIZE * sizeof( double));
    assert( ptgphi != NULL);
    sprintf( tempstr, PTGPHI, pixname);
    sprintf( objname, "%s%s", dbpath, tempstr);
    dmc_err = noDMC_readObject_PIODOUBLE( objname, rank_begring * HPRSIZE, nrings * HPRSIZE, ptgphi);
    if (dmc_err < 0) {
      fprintf( stderr, "Error %d when reading %s\n", dmc_err, objname);
      exit( -1);
    }

    // read theta pointing HPR
    double *ptgthe = malloc( nrings * HPRSIZE * sizeof( double));
    assert( ptgthe != NULL);
    sprintf( tempstr, PTGTHE, pixname);
    sprintf( objname, "%s%s", dbpath, tempstr);
    dmc_err = noDMC_readObject_PIODOUBLE( objname, rank_begring * HPRSIZE, nrings * HPRSIZE, ptgthe);
    if (dmc_err < 0) {
      fprintf( stderr, "Error %d when reading %s\n", dmc_err, objname);
      exit( -1);
    }

    if (param.subdiphpr) {
      // read total dipole HPR
      float *diphpr = malloc( nrings * HPRSIZE * sizeof( float));
      assert( diphpr != NULL);
      sprintf( tempstr, DIPHPR, pixname);
      sprintf( objname, "%s%s", dbpath, tempstr);
      dmc_err = noDMC_readObject_PIOFLOAT( objname, rank_begring * HPRSIZE, nrings * HPRSIZE, diphpr);
      if (dmc_err < 0) {
        fprintf( stderr, "Error %d when reading %s\n", dmc_err, objname);
        exit( -1);
      }
      for (i=0; i<nrings * HPRSIZE; i++) {
        sighpr[i] -= diphpr[i];
      }
      free( diphpr);
    }

    int *badring = NULL;
    if (param.badrings) {
      // read badrings list
      badring = malloc( nrings * sizeof( int));
      assert( badring != NULL);
      sprintf( tempstr, BADRING, pixname);
      sprintf( objname, "%s%s", dbpath, tempstr);
      dmc_err = noDMC_readObject_PIOINT( objname, rank_begring, nrings, badring);
      if (dmc_err < 0) {
        fprintf( stderr, "Error %d when reading %s\n", dmc_err, objname);
        exit( -1);
      }
    }

    for (int s=0; s<n_surveys; s++) {
      char sname[10];
      if (n_surveys == 1) {
        strncpy( sname, "", 10);
      }
      else {
        strncpy( sname, survey_name[s], 10);
      }
      fprintf( stderr, "rank %d/%d, HPR2MAP rings %d-%d for survey %s\n", mpi_rank, mpi_size, rank_begring, rank_endring, sname);

      for (i=0; i<MAPSIZE; i++) {
        mapsig[i] = 0.0;
        maphit[i] = 0;
      }

      // sum HPR in mapsig
      for (int ring = rank_begring; ring <= rank_endring; ring++) {
        if (badring[ring-rank_begring] != 0) {
          continue;
        }
        if ((ring < survey_beg[s]) || (ring > survey_end[s])) {
          continue;
        }
        if ((strcmp( "_odd", sname) == 0) && (ring%2==0)) {
          continue;
        }
        if ((strcmp( "_even", sname) == 0) && (ring%2==1)) {
          continue;
        }
        for (i=0; i<HPRSIZE; i++) {
          int offset = (ring - rank_begring) * HPRSIZE;
          if (hithpr[i + offset] != 0.0) {
            ang2pix_ring( NSIDE, ptgthe[i + offset], ptgphi[i + offset], &ipix);
            mapsig[ipix] += sighpr[i + offset] * hithpr[i + offset];
            maphit[ipix] += hithpr[i + offset];
          }
        }
      }

      if (mpi_rank==0) {
        // sum all rank maps in rank 0
        MPI_Status mpistat;
        double *ranksig = malloc( MAPSIZE * sizeof( double));
        assert( ranksig != NULL);
        int *rankhit = malloc( MAPSIZE * sizeof( int));
        assert( rankhit != NULL);
        for (int rank=1; rank<mpi_size; rank++) {
          fprintf( stderr, "adding rank %d MAP\n", rank);
          MPI_Recv( ranksig, MAPSIZE * sizeof( double), MPI_BYTE, rank, 1031, MPI_COMM_WORLD, &mpistat);
          MPI_Recv( rankhit, MAPSIZE * sizeof( int),    MPI_BYTE, rank, 1032, MPI_COMM_WORLD, &mpistat);
          for (i=0; i<MAPSIZE; i++) {
            mapsig[i] += ranksig[i];
            maphit[i] += rankhit[i];
          }
        }
        free( ranksig);
        free( rankhit);

        // average all pixels
        float *mapfloat = malloc( MAPSIZE * sizeof( float));
        assert( mapfloat != NULL);
        for (i=0; i<MAPSIZE; i++) {
          if (maphit[i] != 0) {
            mapfloat[i] = mapsig[i] / maphit[i];
          }
          else {
            mapfloat[i] = NAN;
          }
        }

        // write map(s) to disk
        char mapname[STRLEN];
        sprintf( mapname, "%s%s", param.sigmap, sname);
        writebin( mapname, mapfloat, MAPSIZE * sizeof( float), 0);
    // TODO: requires #DEFINE_FITSIO in chealpix.*, and fitsio.*
    //    sprintf( tempstr, "%s.fits", param.sigmap);
    //    write_healpix_map( mapfloat, NSIDE, tempstr, 0, "GALACTIC");
        if (param.flag_hitmap) {
          sprintf( mapname, "%s%s", param.hitmap, sname);
          writebin( mapname, maphit, MAPSIZE * sizeof( int), 0);
        }
        free( mapfloat);
      }

      else {
        // send signal and hit maps to rank 0
        MPI_Send( mapsig, MAPSIZE * sizeof( double), MPI_BYTE, 0, 1031, MPI_COMM_WORLD);
        MPI_Send( maphit, MAPSIZE * sizeof( int),    MPI_BYTE, 0, 1032, MPI_COMM_WORLD);
      }

    }

    free( ptgthe);
    free( ptgphi);
    free( mapsig);
    free( maphit);
  }

  free( sighpr);
  free( hithpr);
  MPI_Finalize();
  exit( 0);
}
