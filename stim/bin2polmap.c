
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


#define STRLEN  (260l)
#define HPRSIZE (27664l)
#define NSIDE   (2048l)
#define MAPSIZE (12*NSIDE*NSIDE)
#define MAPSIZE128 (12*128*128)
#define HM1LASTRING (13144)

typedef double MAPTYPE;

enum SURVEY {FULL, HM1, HM2, ODD, EVEN};
char *survey_name[] = {"_full", "_hm1", "_hm2", "_odd", "_even"};
int survey_beg[] = {240, 240, 13145, 240, 240};
int survey_end[] = {26050, HM1LASTRING, 26050, 26050, 26050};
int n_surveys = 5;

typedef struct {
  char dsname[STRLEN];
  char hprname[STRLEN];
  char mapname[STRLEN];
  char hitname[STRLEN];
  int  beginring;
  int  endring;
  int  subdiphpr;
  int  nobandpasscorr;
  int  nopol;
  int  removepsmoffset;
  enum SURVEY survey;
} detset2map_parContent;

int mpi_rank=0;
int mpi_size=0;
MPI_Status mpistat;
char *dbpath = NULL;
char *srollhost = NULL;

enum DATSET {DEC16, DUST17};
enum DATSET DATASET = DEC16;

// functions defined at the end of this file to avoid clutter

void usage( char *msg);
void readParameters( detset2map_parContent *param, int argc, char *argv[]);
void loadbal( int mpi_rank,  int mpi_size, int global_first, int global_last, int *rank_first, int *rank_last); 
MAPTYPE solvemap( MAPTYPE *x, MAPTYPE *y, MAPTYPE *z,
                  MAPTYPE II, MAPTYPE IQ, MAPTYPE IU,
                  MAPTYPE QQ, MAPTYPE QU, MAPTYPE UU);
double get_crosspol( char *pixname);
double get_DX11calib( char *pixname);
double get_DX11nep( char *pixname);
double get_RD12calib( char *pixname);
double get_RD12nep( char *pixname);
double get_FFP10PSM_offset( char *pixname);
double get_freefree_coef( char *pixname);
double get_dust_coef( char *pixname);
double get_co_coef( char *pixname);
void   get_DEC16_bandpass_correction( char *pixname, double *dust_coeff, double *freefree_coeff, double *co_coeff);
void   get_julia_bandpass_correction( char *pixname, double *dust_coeff, double *freefree_coeff, double *co_coeff);


////////////////////////////////////////////////////////////////////////////////

int main( int argc, char *argv[]) {

  long i, ibol, ipix, ipix128;
  int  dmc_err;
  char pixname[10];
  char *boloid;
  char objname[STRLEN];
  char rpath[PATH_MAX];

  int  n_bolo;
  char **bolo_list;

  detset2map_parContent param;
  
  // initialise MPI
  MPI_Init( &argc, &argv);  
  MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size( MPI_COMM_WORLD, &mpi_size);

  realpath( argv[0], rpath);
  if (mpi_rank == 0) {
    fprintf( stderr, "starting %s\n" , rpath);
  }

  // get DBPATH from environment variable
  srollhost = getenv("SROLLHOST");
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
  
  // read command line parameters
  readParameters( &param, argc, argv);

  // get list of bolometer in detset
  getDetset( param.dsname, &bolo_list, &n_bolo);
  
  // get begin and end ring per mpi rank
  if (mpi_rank == 0) fprintf( stderr, "\n");
  int rank_begring, rank_endring;
  loadbal( mpi_rank,  mpi_size, param.beginring, param.endring, &rank_begring, &rank_endring);
  int nrings  = rank_endring - rank_begring + 1;
  for (i=0; i<mpi_size; i++) {
    if (i == mpi_rank) {
      fprintf( stderr, "rank %d/%d, rank_begring=%d, rank_endring=%d (%d rings)\n", mpi_rank, mpi_size, rank_begring, rank_endring, nrings);
    }
    MPI_Barrier( MPI_COMM_WORLD);
  }

  // allocate input HPRs memory
  float *sighpr = malloc( nrings * HPRSIZE * sizeof( float));
  assert( sighpr != NULL);
  float *hithpr = malloc( nrings * HPRSIZE * sizeof( float));
  assert( hithpr != NULL);
  double *ptgphi = malloc( nrings * HPRSIZE * sizeof( double));
  assert( ptgphi != NULL);
  double *ptgthe = malloc( nrings * HPRSIZE * sizeof( double));
  assert( ptgthe != NULL);
  double *ptgpsi = malloc( nrings * HPRSIZE * sizeof( double));
  assert( ptgpsi != NULL);

  float *diphpr = NULL;
  if (param.subdiphpr) {
    diphpr = malloc( nrings * HPRSIZE * sizeof( float));
    assert( diphpr != NULL);
  }

  int *badrings = malloc( 26051 * sizeof( int));
  assert( badrings != NULL);

  // allocate maps
  MAPTYPE *mapSI=NULL, *mapII=NULL;
  mapSI = malloc( MAPSIZE * sizeof( MAPTYPE));
  assert( mapSI != NULL);
  mapII = malloc( MAPSIZE * sizeof( MAPTYPE));
  assert( mapII != NULL);

  int *mapHIT=NULL;
  if (param.hitname[0] != 0) {
    mapHIT = malloc( MAPSIZE * sizeof( int));
    assert( mapHIT != NULL);
  }

  MAPTYPE *mapSQ=NULL, *mapSU=NULL, *mapIQ=NULL, *mapIU=NULL, *mapQQ=NULL, *mapQU=NULL, *mapUU=NULL;
  if (param.nopol == 0) {
    mapSQ = malloc( MAPSIZE * sizeof( MAPTYPE));
    assert( mapSQ != NULL);
    mapSU = malloc( MAPSIZE * sizeof( MAPTYPE));
    assert( mapSU != NULL);
    mapIQ = malloc( MAPSIZE * sizeof( MAPTYPE));
    assert( mapIQ != NULL);
    mapIU = malloc( MAPSIZE * sizeof( MAPTYPE));
    assert( mapIU != NULL);
    mapQQ = malloc( MAPSIZE * sizeof( MAPTYPE));
    assert( mapQQ != NULL);
    mapQU = malloc( MAPSIZE * sizeof( MAPTYPE));
    assert( mapQU != NULL);
    mapUU = malloc( MAPSIZE * sizeof( MAPTYPE));
    assert( mapUU != NULL);
  }

  for (i=0; i<MAPSIZE; i++) {
    mapSI[i] = 0.0;
    mapII[i] = 0.0;
    if (param.hitname[0] != 0) {
      mapHIT[i] = 0;
    }
    if (param.nopol == 0) {
      mapSQ[i] = 0.0;
      mapSU[i] = 0.0;
      mapIQ[i] = 0.0;
      mapIU[i] = 0.0;
      mapQQ[i] = 0.0;
      mapQU[i] = 0.0;
      mapUU[i] = 0.0;
    }
  }

  float *comap = NULL;
  float *dustmapI = NULL;
  float *dustmapQ = NULL;
  float *dustmapU = NULL;
  float *freefree = NULL;

  if ( strstr( param.hprname, "dustapr17") != NULL) {
    DATASET = DUST17;
    if (mpi_rank == 0) fprintf( stderr, "Using DUST17 bandpasses and PSM monopoles\n");
  }
  else {
    DATASET = DEC16;
    if (mpi_rank == 0) fprintf( stderr, "Using DEC16 bandpasses and PSM monopoles\n");
  }

  for (ibol=0; ibol<n_bolo; ibol++) {
    strcpy( pixname, bolo_list[ibol]);
    double crosspol       = get_crosspol( pixname);
    double calib          = get_DX11calib( pixname); // DEC16v1 sims use DX11 calib and RD12 nep...
    double nep            = get_RD12nep( pixname);
    double offset         = get_FFP10PSM_offset( pixname);
    double dust_coeff     = get_dust_coef( pixname);
    double freefree_coeff = get_freefree_coef( pixname);
    double co_coeff       = get_co_coef( pixname);
    boloid = toBoloID( pixname);

    if (param.nopol == 1) {
      crosspol = 1.0;
    }

    if (param.nobandpasscorr == 1) {
      dust_coeff     = 0.0;
      freefree_coeff = 0.0;
      co_coeff       = 0.0;
    }

    if (param.removepsmoffset == 0) {
      offset = 0.0;
    }

    double etha   = (1.0 - crosspol) / (1.0 + crosspol);
    double weight = (calib / nep) * (calib / nep);

    if (mpi_rank == 0) {
      fprintf( stderr, "\nprocessing %s\n", pixname);
      fprintf( stderr, "  crosspol   = %g\n", crosspol);
      fprintf( stderr, "  calib      = %.12g\n", calib);
      fprintf( stderr, "  nep        = %.12g\n", nep);
      fprintf( stderr, "  PSM offset = %.12g\n", offset);
      fprintf( stderr, "  co coef    = %.12g\n", co_coeff);
      fprintf( stderr, "  dust coef  = %.12g\n", dust_coeff);
      fprintf( stderr, "  freefree   = %.12g\n", freefree_coeff);
    }

    // read bandpass correction maps
    if ((co_coeff != 0) && (comap == NULL)) {
      comap = malloc( MAPSIZE128 * sizeof( float));
      assert ( comap != NULL); 
      sprintf( objname, "%s/MAP_JMD_128/COMAP", dbpath);
      if (mpi_rank == 0) fprintf( stderr, "reading: %s\n", objname);
      dmc_err = noDMC_readObject_PIOFLOAT( objname, 0, MAPSIZE128, comap);
      if (dmc_err < 0) {
        fprintf( stderr, "Error %d when reading %s\n", dmc_err, objname);
        exit( -1);
      }
    }
  
    if ((dust_coeff != 0) && (dustmapI == NULL)) {
      dustmapI = malloc( MAPSIZE128 * sizeof( float));
      assert ( dustmapI != NULL); 
      sprintf( objname, "%s/MAP_JMD_128/Dust_I", dbpath);
      if (mpi_rank == 0) fprintf( stderr, "reading: %s\n", objname);
      dmc_err = noDMC_readObject_PIOFLOAT( objname, 0, MAPSIZE128, dustmapI);
      if (dmc_err < 0) {
        fprintf( stderr, "Error %d when reading %s\n", dmc_err, objname);
        exit( -1);
      }

      dustmapQ = malloc( MAPSIZE128 * sizeof( float));
      assert ( dustmapQ != NULL); 
      sprintf( objname, "%s/MAP_JMD_128/Dust_Q", dbpath);
      if (mpi_rank == 0) fprintf( stderr, "reading: %s\n", objname);
      dmc_err = noDMC_readObject_PIOFLOAT( objname, 0, MAPSIZE128, dustmapQ);
      if (dmc_err < 0) {
        fprintf( stderr, "Error %d when reading %s\n", dmc_err, objname);
        exit( -1);
      }

      dustmapU = malloc( MAPSIZE128 * sizeof( float));
      assert ( dustmapU != NULL); 
      sprintf( objname, "%s/MAP_JMD_128/Dust_U", dbpath);
      if (mpi_rank == 0) fprintf( stderr, "reading: %s\n", objname);
      dmc_err = noDMC_readObject_PIOFLOAT( objname, 0, MAPSIZE128, dustmapU);
      if (dmc_err < 0) {
        fprintf( stderr, "Error %d when reading %s\n", dmc_err, objname);
        exit( -1);
      }
    }
  
    if ((freefree_coeff != 0) && (freefree == NULL)) {
      freefree = malloc( MAPSIZE128 * sizeof( float));
      assert ( freefree != NULL); 
      sprintf( objname, "%s/MAP_JMD_128/freefree_NORM_I", dbpath);
      if (mpi_rank == 0) fprintf( stderr, "reading: %s\n", objname);
      dmc_err = noDMC_readObject_PIOFLOAT( objname, 0, MAPSIZE128, freefree);
      if (dmc_err < 0) {
        fprintf( stderr, "Error %d when reading %s\n", dmc_err, objname);
        exit( -1);
      }
    }
  
  
    // read signal HPR
    char *hprname = str_replace ( param.hprname, "{pixname}", pixname);
    if (mpi_rank == 0) fprintf( stderr, "reading input HPR: %s\n", hprname);

    dmc_err = noDMC_readObject_PIOFLOAT( hprname, rank_begring * HPRSIZE, nrings * HPRSIZE, sighpr);
    if (dmc_err < 0) {
      fprintf( stderr, "Error %d when reading %s\n", dmc_err, param.hprname);
      exit( -1);
    }

    // read bad rings list
//    sprintf( objname, "%s/calROIs/%s_bad_rings_and_offsets_rd12", dbpath, pixname);
    sprintf( objname, "%s/calROIs/%s_discarded_rings_dx11", dbpath, boloid);
    if (mpi_rank == 0) fprintf( stderr, "reading: %s\n", objname);
    dmc_err = noDMC_readObject_PIOINT( objname, 0, 26051, badrings);
    if (dmc_err < 0) {
      if (mpi_rank == 0) fprintf( stderr, "!!! Error %d when reading %s, won't use badrings\n", dmc_err, objname);
      exit( -1);
      for (i=0; i<26051; i++) {
        badrings[i] = 0;
      }
    }

    // read hitcount HPR
    sprintf( objname, "%s/PBR_JMD/%s_REP6_hit", dbpath, pixname);
    if (mpi_rank == 0) fprintf( stderr, "reading: %s\n", objname);
    dmc_err = noDMC_readObject_PIOFLOAT( objname, rank_begring * HPRSIZE, nrings * HPRSIZE, hithpr);
    if (dmc_err < 0) {
      fprintf( stderr, "Error %d when reading %s\n", dmc_err, objname);
      exit( -1);
    }

    // read phi pointing HPR
    sprintf( objname, "%s/PBR_JMD/%s_REP6_ptg", dbpath, pixname);
    if (mpi_rank == 0) fprintf( stderr, "reading: %s\n", objname);
    dmc_err = noDMC_readObject_PIODOUBLE( objname, rank_begring * HPRSIZE, nrings * HPRSIZE, ptgphi);
    if (dmc_err < 0) {
      fprintf( stderr, "Error %d when reading %s\n", dmc_err, objname);
      exit( -1);
    }
  
    // read theta pointing HPR
    sprintf( objname, "%s/PBR_JMD/%s_REP6_ptg_TUPLE_1", dbpath, pixname);
    if (mpi_rank == 0) fprintf( stderr, "reading: %s\n", objname);
    dmc_err = noDMC_readObject_PIODOUBLE( objname, rank_begring * HPRSIZE, nrings * HPRSIZE, ptgthe);
    if (dmc_err < 0) {
      fprintf( stderr, "Error %d when reading %s\n", dmc_err, objname);
      exit( -1);
    }

    // read psi pointing HPR
    sprintf( objname, "%s/PBR_JMD/%s_REP6_ptg_TUPLE_2", dbpath, pixname);
    if (mpi_rank == 0) fprintf( stderr, "reading: %s\n", objname);
    dmc_err = noDMC_readObject_PIODOUBLE( objname, rank_begring * HPRSIZE, nrings * HPRSIZE, ptgpsi);
    if (dmc_err < 0) {
      fprintf( stderr, "Error %d when reading %s\n", dmc_err, objname);
      exit( -1);
    }

    // read and subtract dipole HPR
    if (param.subdiphpr) {
      sprintf( objname, "%s/HPRREPSM_PBR/%s_v64_MAY16_diporb", dbpath, pixname);
      if (mpi_rank == 0) fprintf( stderr, "reading: %s\n", objname);
      dmc_err = noDMC_readObject_PIOFLOAT( objname, rank_begring * HPRSIZE, nrings * HPRSIZE, diphpr);
      if (dmc_err < 0) {
        fprintf( stderr, "Error %d when reading %s\n", dmc_err, objname);
        exit( -1);
      }
      for (i=0; i<nrings * HPRSIZE; i++) {
        // input signal at 545-857 is in MJYy/Sr, while diphpr is in KCMB
        sighpr[i] -= diphpr[i] * get_KCMB2MJyPerSr( pixname);
      }
    }

    for (int ring = rank_begring; ring <= rank_endring; ring++) {
      // don't process excluded rings
      if ((param.survey == HM1)  && (ring >  HM1LASTRING)) continue;
      if ((param.survey == HM2)  && (ring <= HM1LASTRING)) continue;
      if ((param.survey == ODD)  && (ring%2 == 0)) continue;
      if ((param.survey == EVEN) && (ring%2 == 1)) continue;
      if (badrings[ring] != 0) continue;

      for (i = (ring - rank_begring) * HPRSIZE; i < (ring+1 - rank_begring) * HPRSIZE; i++) {
        if (hithpr[i] == 0) continue; // empty pixel
        ang2pix_ring( NSIDE, ptgthe[i], ptgphi[i], &ipix);
        ang2pix_ring( 128,   ptgthe[i], ptgphi[i], &ipix128);
        
        // sroll.c hpix.co and hpix.si are PIOFLOAT
        double ethacos = etha * (float)cos( 2 * ptgpsi[i]);
        double ethasin = etha * (float)sin( 2 * ptgpsi[i]);

        // bandpass leakage correction
        if (dust_coeff != 0) {
          sighpr[i] -= dust_coeff * (dustmapI[ipix128] + ethacos * dustmapQ[ipix128] + ethasin * dustmapU[ipix128]);
        }
        if (co_coeff != 0) {
          sighpr[i] -= co_coeff * comap[ipix128];
        }
        if (freefree_coeff != 0) {
          sighpr[i] -= freefree_coeff * freefree[ipix128];
        }

        // sroll.c hpix.w is PIOFLOAT
        float wh = weight * hithpr[i];
  
        // temperature
        mapSI[ipix] += wh * (sighpr[i] - offset);
        mapII[ipix] += wh;

        // hitcount
        if (param.hitname[0] != 0) {
          mapHIT[ipix] += hithpr[i];
        }

        // polarisation
        if (param.nopol == 0) {
          mapSQ[ipix] += wh * (sighpr[i] - offset) * ethacos;
          mapSU[ipix] += wh * (sighpr[i] - offset) * ethasin;
          mapIQ[ipix] += wh * ethacos;
          mapIU[ipix] += wh * ethasin;
          mapQQ[ipix] += wh * ethacos * ethacos;
          mapQU[ipix] += wh * ethacos * ethasin;
          mapUU[ipix] += wh * ethasin * ethasin;
        }
      }
    } // for (ring = rank_begring)
  } // while (ibol++ < n_bolo)


// sum hitcount maps on rank 0 and write to disk
  if (param.hitname[0] != 0) {
    if (mpi_rank==0) {
      int *buf = malloc( MAPSIZE * sizeof( int));
      assert( buf != NULL);
      for (int send_rank=1; send_rank<mpi_size; send_rank++) {
        MPI_Recv( buf, MAPSIZE * sizeof( int), MPI_BYTE, send_rank, 1039, MPI_COMM_WORLD, &mpistat);
        for (i=0; i<MAPSIZE; i++) mapHIT[i] += buf[i];
      }
      free( buf);
      writebin( param.hitname, mapHIT, MAPSIZE * sizeof( int), 0);
    }
    else {
      MPI_Send( mapHIT, MAPSIZE * sizeof( int), MPI_BYTE, 0, 1039, MPI_COMM_WORLD);
    }
  }


// gather pixel interval per rank
  if (mpi_rank == 0) fprintf( stderr, "\n");
  int first_pix, last_pix, pix_count;
  for (int recv_rank=0; recv_rank<mpi_size; recv_rank++) {
    loadbal( recv_rank,  mpi_size, 0, MAPSIZE-1, &first_pix, &last_pix);
    pix_count = last_pix - first_pix + 1;
    if (recv_rank == mpi_rank) {
      fprintf( stderr, "rank %d gathering data for pixels %d-%d\n", recv_rank, first_pix, last_pix);
      MAPTYPE *buf = malloc( pix_count * sizeof( MAPTYPE));
      assert( buf != NULL);
      for (int send_rank=0; send_rank<mpi_size; send_rank++) {
        if (send_rank == mpi_rank) continue;
        MPI_Recv( buf, pix_count * sizeof( MAPTYPE), MPI_BYTE, send_rank, 1030, MPI_COMM_WORLD, &mpistat);
        for (i=0; i<pix_count; i++) mapSI[first_pix + i] += buf[i];
        MPI_Recv( buf, pix_count * sizeof( MAPTYPE), MPI_BYTE, send_rank, 1033, MPI_COMM_WORLD, &mpistat);
        for (i=0; i<pix_count; i++) mapII[first_pix + i] += buf[i];
        if (param.nopol == 0) {
          MPI_Recv( buf, pix_count * sizeof( MAPTYPE), MPI_BYTE, send_rank, 1031, MPI_COMM_WORLD, &mpistat);
          for (i=0; i<pix_count; i++) mapSQ[first_pix + i] += buf[i];
          MPI_Recv( buf, pix_count * sizeof( MAPTYPE), MPI_BYTE, send_rank, 1032, MPI_COMM_WORLD, &mpistat);
          for (i=0; i<pix_count; i++) mapSU[first_pix + i] += buf[i];
          MPI_Recv( buf, pix_count * sizeof( MAPTYPE), MPI_BYTE, send_rank, 1034, MPI_COMM_WORLD, &mpistat);
          for (i=0; i<pix_count; i++) mapIQ[first_pix + i] += buf[i];
          MPI_Recv( buf, pix_count * sizeof( MAPTYPE), MPI_BYTE, send_rank, 1035, MPI_COMM_WORLD, &mpistat);
          for (i=0; i<pix_count; i++) mapIU[first_pix + i] += buf[i];
          MPI_Recv( buf, pix_count * sizeof( MAPTYPE), MPI_BYTE, send_rank, 1036, MPI_COMM_WORLD, &mpistat);
          for (i=0; i<pix_count; i++) mapQQ[first_pix + i] += buf[i];
          MPI_Recv( buf, pix_count * sizeof( MAPTYPE), MPI_BYTE, send_rank, 1037, MPI_COMM_WORLD, &mpistat);
          for (i=0; i<pix_count; i++) mapQU[first_pix + i] += buf[i];
          MPI_Recv( buf, pix_count * sizeof( MAPTYPE), MPI_BYTE, send_rank, 1038, MPI_COMM_WORLD, &mpistat);
          for (i=0; i<pix_count; i++) mapUU[first_pix + i] += buf[i];
        }
      }
      free( buf);
    }
    else {
      MPI_Send( &mapSI[first_pix], pix_count * sizeof( MAPTYPE), MPI_BYTE, recv_rank, 1030, MPI_COMM_WORLD);
      MPI_Send( &mapII[first_pix], pix_count * sizeof( MAPTYPE), MPI_BYTE, recv_rank, 1033, MPI_COMM_WORLD);
      if (param.nopol == 0) {
        MPI_Send( &mapSQ[first_pix], pix_count * sizeof( MAPTYPE), MPI_BYTE, recv_rank, 1031, MPI_COMM_WORLD);
        MPI_Send( &mapSU[first_pix], pix_count * sizeof( MAPTYPE), MPI_BYTE, recv_rank, 1032, MPI_COMM_WORLD);
        MPI_Send( &mapIQ[first_pix], pix_count * sizeof( MAPTYPE), MPI_BYTE, recv_rank, 1034, MPI_COMM_WORLD);
        MPI_Send( &mapIU[first_pix], pix_count * sizeof( MAPTYPE), MPI_BYTE, recv_rank, 1035, MPI_COMM_WORLD);
        MPI_Send( &mapQQ[first_pix], pix_count * sizeof( MAPTYPE), MPI_BYTE, recv_rank, 1036, MPI_COMM_WORLD);
        MPI_Send( &mapQU[first_pix], pix_count * sizeof( MAPTYPE), MPI_BYTE, recv_rank, 1037, MPI_COMM_WORLD);
        MPI_Send( &mapUU[first_pix], pix_count * sizeof( MAPTYPE), MPI_BYTE, recv_rank, 1038, MPI_COMM_WORLD);
      }
    }
  }

// solve for current mpi_rank pixel interval

  int ninv_count = 0; // number of pixels with non-inversible matrix
  if (mpi_rank == 0) fprintf( stderr, "\nsolving MAP\n\n");
  loadbal( mpi_rank,  mpi_size, 0, MAPSIZE-1, &first_pix, &last_pix);
  pix_count = last_pix -first_pix + 1;
  for (i=first_pix; i<=last_pix; i++) {
    if (mapII[i] != 0.0) {
      if (param.nopol == 0) {
        solvemap( &mapSI[i], &mapSQ[i], &mapSU[i],
                   mapII[i],  mapIQ[i],  mapIU[i],
                   mapQQ[i],  mapQU[i],  mapUU[i]);
        if (!isfinite( mapSI[i])) {
          ninv_count++;
          mapSI[i] = HEALPIX_NULLVAL;
          mapSQ[i] = HEALPIX_NULLVAL;
          mapSU[i] = HEALPIX_NULLVAL;
        }
      }
      else {
        mapSI[i] /= mapII[i];
      }
    }
    else {
      mapSI[i] = HEALPIX_NULLVAL;
      if (param.nopol == 0) {
        mapSQ[i] = HEALPIX_NULLVAL;
        mapSU[i] = HEALPIX_NULLVAL;
      }
    }
  }
  if (ninv_count != 0) {
    fprintf( stderr, "rank %d, non inversible pixel count = %d\n", mpi_rank, ninv_count);
  }

// gather maps on rank 0 and write to disk
  
  if (mpi_rank==0) {
    // sum all rank maps in rank 0
    for (int send_rank=1; send_rank<mpi_size; send_rank++) {
      loadbal( send_rank,  mpi_size, 0, MAPSIZE-1, &first_pix, &last_pix);
      fprintf( stderr, "recieving pixels %d-%d from rank %d\n", first_pix, last_pix, send_rank);
      pix_count = last_pix - first_pix + 1;
      MPI_Recv( &mapSI[first_pix], pix_count * sizeof( MAPTYPE), MPI_BYTE, send_rank, 1040, MPI_COMM_WORLD, &mpistat);
      MPI_Recv( &mapII[first_pix], pix_count * sizeof( MAPTYPE), MPI_BYTE, send_rank, 1043, MPI_COMM_WORLD, &mpistat);
      if (param.nopol == 0) {
        MPI_Recv( &mapSQ[first_pix], pix_count * sizeof( MAPTYPE), MPI_BYTE, send_rank, 1041, MPI_COMM_WORLD, &mpistat);
        MPI_Recv( &mapSU[first_pix], pix_count * sizeof( MAPTYPE), MPI_BYTE, send_rank, 1042, MPI_COMM_WORLD, &mpistat);
        MPI_Recv( &mapIQ[first_pix], pix_count * sizeof( MAPTYPE), MPI_BYTE, send_rank, 1044, MPI_COMM_WORLD, &mpistat);
        MPI_Recv( &mapIU[first_pix], pix_count * sizeof( MAPTYPE), MPI_BYTE, send_rank, 1045, MPI_COMM_WORLD, &mpistat);
        MPI_Recv( &mapQQ[first_pix], pix_count * sizeof( MAPTYPE), MPI_BYTE, send_rank, 1046, MPI_COMM_WORLD, &mpistat);
        MPI_Recv( &mapQU[first_pix], pix_count * sizeof( MAPTYPE), MPI_BYTE, send_rank, 1047, MPI_COMM_WORLD, &mpistat);
        MPI_Recv( &mapUU[first_pix], pix_count * sizeof( MAPTYPE), MPI_BYTE, send_rank, 1048, MPI_COMM_WORLD, &mpistat);
      }
    }

    // write map(s) to disk
    float *floatmap = malloc( MAPSIZE * sizeof( float));
    assert( floatmap != NULL);
    char mapname[STRLEN];
    if (mpi_rank == 0) fprintf( stderr, "\n");
  
    sprintf( mapname, "%s_I", param.mapname);
    for (i=0; i<MAPSIZE; i++) floatmap[i] = mapSI[i];
    writebin( mapname, floatmap, MAPSIZE * sizeof( float), 0);

    sprintf( mapname, "%s_II", param.mapname);
    for (i=0; i<MAPSIZE; i++) floatmap[i] = mapII[i];
    writebin( mapname, floatmap, MAPSIZE * sizeof( float), 0);

    if (param.nopol == 0) {
      sprintf( mapname, "%s_Q", param.mapname);
      for (i=0; i<MAPSIZE; i++) floatmap[i] = mapSQ[i];
      writebin( mapname, floatmap, MAPSIZE * sizeof( float), 0);
  
      sprintf( mapname, "%s_U", param.mapname);
      for (i=0; i<MAPSIZE; i++) floatmap[i] = mapSU[i];
      writebin( mapname, floatmap, MAPSIZE * sizeof( float), 0);
  
      sprintf( mapname, "%s_IQ", param.mapname);
      for (i=0; i<MAPSIZE; i++) floatmap[i] = mapIQ[i];
      writebin( mapname, floatmap, MAPSIZE * sizeof( float), 0);
  
      sprintf( mapname, "%s_IU", param.mapname);
      for (i=0; i<MAPSIZE; i++) floatmap[i] = mapIU[i];
      writebin( mapname, floatmap, MAPSIZE * sizeof( float), 0);
  
      sprintf( mapname, "%s_QQ", param.mapname);
      for (i=0; i<MAPSIZE; i++) floatmap[i] = mapQQ[i];
      writebin( mapname, floatmap, MAPSIZE * sizeof( float), 0);
  
      sprintf( mapname, "%s_QU", param.mapname);
      for (i=0; i<MAPSIZE; i++) floatmap[i] = mapQU[i];
      writebin( mapname, floatmap, MAPSIZE * sizeof( float), 0);
  
      sprintf( mapname, "%s_UU", param.mapname);
      for (i=0; i<MAPSIZE; i++) floatmap[i] = mapUU[i];
      writebin( mapname, floatmap, MAPSIZE * sizeof( float), 0);
    }

    free( floatmap);
  }

  else {
    // send current rank pixel interval and hit maps to rank 0
    loadbal( mpi_rank,  mpi_size, 0, MAPSIZE-1, &first_pix, &last_pix);
    pix_count = last_pix - first_pix + 1;
    MPI_Send( &mapSI[first_pix], pix_count * sizeof( MAPTYPE), MPI_BYTE, 0, 1040, MPI_COMM_WORLD);
    MPI_Send( &mapII[first_pix], pix_count * sizeof( MAPTYPE), MPI_BYTE, 0, 1043, MPI_COMM_WORLD);
    if (param.nopol == 0) {
      MPI_Send( &mapSQ[first_pix], pix_count * sizeof( MAPTYPE), MPI_BYTE, 0, 1041, MPI_COMM_WORLD);
      MPI_Send( &mapSU[first_pix], pix_count * sizeof( MAPTYPE), MPI_BYTE, 0, 1042, MPI_COMM_WORLD);
      MPI_Send( &mapIQ[first_pix], pix_count * sizeof( MAPTYPE), MPI_BYTE, 0, 1044, MPI_COMM_WORLD);
      MPI_Send( &mapIU[first_pix], pix_count * sizeof( MAPTYPE), MPI_BYTE, 0, 1045, MPI_COMM_WORLD);
      MPI_Send( &mapQQ[first_pix], pix_count * sizeof( MAPTYPE), MPI_BYTE, 0, 1046, MPI_COMM_WORLD);
      MPI_Send( &mapQU[first_pix], pix_count * sizeof( MAPTYPE), MPI_BYTE, 0, 1047, MPI_COMM_WORLD);
      MPI_Send( &mapUU[first_pix], pix_count * sizeof( MAPTYPE), MPI_BYTE, 0, 1048, MPI_COMM_WORLD);
    }
  }

  // free mem and bye
  free( sighpr);
  free( hithpr);
  free( ptgthe);
  free( ptgphi);
  free( ptgpsi);
  if (param.subdiphpr) {
    free( diphpr);
  }
  free( mapSI);
  free( mapII);
  if (param.nopol == 0) {
    free( mapIQ);
    free( mapIU);
    free( mapQQ);
    free( mapQU);
    free( mapUU);
    free( mapSQ);
    free( mapSU);
  }

  if (comap != NULL) {
    free( comap);
  }
  if (dustmapI != NULL) {
    free( dustmapI);
    free( dustmapQ);
    free( dustmapU);
  }
  if (freefree != NULL) {
    free( freefree);
  }

  if (mpi_rank == 0) {
    fprintf( stderr, "\n%s Finished Successfully\n\n" , rpath);
  }

  MPI_Finalize();
  exit( 0);
}


////////////////////////////////////////////////////////////////////////////////

void usage( char *msg) {

  if (mpi_rank == 0) {
    fprintf( stderr, msg);
    fprintf( stderr, "\n\n");
    fprintf( stderr, "bin2polmap: multi-bolometer projection of HPR to MAP\n");
    fprintf( stderr, "usage: bin2polmap --detset=<detset> --hprname=<hprname> --mapname=<mapname> [--survey=<survey>] [--subdiphpr] [--beginring=<beginring>] [--endring=<endring>]\n");
    fprintf( stderr, "       - <detset>: name of detset (e.g. 100GHz, 143ds1, 353PSB), case insensitive, mandatory\n");
    fprintf( stderr, "       - <hprname>: name of input float32 HPR in HFI DMC or flat bianry format, must include the {pixname} or {boloid} strings, mandatory\n");
    fprintf( stderr, "       - <mapname>: name of output float32 map containing binned signal, mandatory\n");
    fprintf( stderr, "       - <hitname>: optional, name of output int32 hitcount map\n");
    fprintf( stderr, "       - <survey>: one of [hm1|hm2|odd|even], only these rings will be used for projection\n");
    fprintf( stderr, "       - --subdiphpr: if present, total dipole HPR will be subtracted from input signal HPR\n");
    fprintf( stderr, "       - --nobandpasscorr: if present, no bandpass correction will be applied\n");
    fprintf( stderr, "       - --nopol: if present, maps will be produced without polarisation\n");
    fprintf( stderr, "       - --removepsmoffset: if present, the PSM offsets per bolometer will be subtracted from the HPR signal\n");
    fprintf( stderr, "       - <beginring>: first ring number to process, 240 if omitted\n");
    fprintf( stderr, "       - <endring>: last ring number to process, 26050 if omitted\n");
    fprintf( stderr, "\n");
  }
  MPI_Finalize();
  exit( -1);
}

////////////////////////////////////////////////////////////////////////////////

void readParameters( detset2map_parContent *param, int argc, char *argv[]) {

  char tempstr[STRLEN];
  param->dsname[0]  = 0;
  param->hprname[0] = 0;
  param->mapname[0] = 0;
  param->hitname[0] = 0;
  param->survey     = FULL;
  param->beginring  = 240;
  param->endring    = 26050;
  param->subdiphpr  = 0;
  param->nopol      = 0;
  param->nobandpasscorr  = 0;
  param->removepsmoffset = 0;


  for (int argi=1; argi<argc; argi++) {
    if (startswith( argv[argi], "--detset=")) {
      strncpy( param->dsname, argv[argi] + strlen( "--detset="), STRLEN);
      if (mpi_rank==0) fprintf( stderr, "parameter detset=%s\n", param->dsname);
    }
    else if (startswith( argv[argi], "--hprname=")) {
      strncpy( param->hprname, argv[argi] + strlen( "--hprname="), STRLEN);
      if (mpi_rank==0) fprintf( stderr, "parameter hprname=%s\n", param->hprname);
    }
    else if (startswith( argv[argi], "--mapname=")) {
      strncpy( param->mapname, argv[argi] + strlen( "--mapname="), STRLEN);
      if (mpi_rank==0) fprintf( stderr, "parameter mapname=%s\n", param->mapname);
    }
    else if (startswith( argv[argi], "--hitname=")) {
      strncpy( param->hitname, argv[argi] + strlen( "--hitname="), STRLEN);
      if (mpi_rank==0) fprintf( stderr, "parameter hitname=%s\n", param->hitname);
    }
    else if (startswith( argv[argi], "--survey=")) {
      strncpy( tempstr, argv[argi] + strlen( "--survey="), STRLEN);
      lowercase( tempstr);
      if (strcmp( tempstr, "hm1") == 0) {
        param->survey = HM1;
        if (param->endring > HM1LASTRING) param->endring = HM1LASTRING;
      }
      else if (strcmp( tempstr, "hm2") == 0) {
        param->survey = HM2;
        if (param->beginring <= HM1LASTRING) param->beginring = HM1LASTRING+1;
      }
      else if (strcmp( tempstr, "odd") == 0) param->survey = ODD;
      else if (strcmp( tempstr, "even") == 0) param->survey = EVEN;
      else usage( "Error: unknown <survey> value");
      if (mpi_rank==0) fprintf( stderr, "parameter survey=%s\n", tempstr);
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
    else if (strcmp( argv[argi], "--nobandpasscorr") == 0) {
      param->nobandpasscorr = 1;
      if (mpi_rank==0) fprintf( stderr, "parameter nobandpasscorr=1\n");
    }
    else if (strcmp( argv[argi], "--nopol") == 0) {
      param->nopol = 1;
      if (mpi_rank==0) fprintf( stderr, "parameter nopol=1\n");
    }
    else if (strcmp( argv[argi], "--removepsmoffset") == 0) {
      param->removepsmoffset = 1;
      if (mpi_rank==0) fprintf( stderr, "parameter removepsmoffset=1\n");
    }
    else {
      fprintf( stderr, "Unknown parameter: %s\n\n", argv[argi]);
      usage( NULL);
      exit( -1);
    }
  }
  if (param->dsname[0] == 0)
    usage( "Error: input <detset> is mandatory");
  if (param->hprname[0] == 0)
    usage( "Error: input <hprname> is mandatory");
  if (param->mapname[0] == 0)
    usage( "Error: input <mapname> is mandatory");
}


////////////////////////////////////////////////////////////////////////

void loadbal( int mpi_rank,  int mpi_size, int global_first, int global_last, int *rank_first, int *rank_last) {
  int task_per_rank  = (global_last - global_first + 1) / mpi_size;
  int task_remainder = (global_last - global_first + 1) % mpi_size;
  int rankfirst_offset = task_remainder;
  int ranklast_offset  = task_remainder;
  if (mpi_rank < task_remainder) rankfirst_offset = mpi_rank;
  if (mpi_rank+1 < task_remainder) ranklast_offset = mpi_rank+1;
  *rank_first = global_first + mpi_rank * task_per_rank + rankfirst_offset;
  *rank_last  = global_first + (mpi_rank + 1) * task_per_rank + ranklast_offset - 1;
}


////////////////////////////////////////////////////////////////////////

MAPTYPE solvemap( MAPTYPE *x, MAPTYPE *y, MAPTYPE *z,
                  MAPTYPE II, MAPTYPE IQ, MAPTYPE IU,
                  MAPTYPE QQ, MAPTYPE QU, MAPTYPE UU)
{
  //double determinant =    +A(0,0)*(A(1,1)*A(2,2)-A(2,1)*A(1,2))
  //                     -A(0,1)*(A(1,0)*A(2,2)-A(1,2)*A(2,0))
  //                    +A(0,2)*(A(1,0)*A(2,1)-A(1,1)*A(2,0));
  MAPTYPE determinant = II*(QQ*UU-QU*QU)
                       -IQ*(IQ*UU-QU*IU)
                       +IU*(IQ*QU-QQ*IU);
  
  MAPTYPE invdet = 1/determinant;
  MAPTYPE result_0_0 =  (QQ*UU - QU*QU)*invdet; // (A(1,1)*A(2,2)-A(2,1)*A(1,2))*invdet;
  MAPTYPE result_1_0 = -(IQ*UU - IU*QU)*invdet; //-(A(0,1)*A(2,2)-A(0,2)*A(2,1))*invdet;
  MAPTYPE result_2_0 =  (IQ*QU - IU*QQ)*invdet; // (A(0,1)*A(1,2)-A(0,2)*A(1,1))*invdet;
  MAPTYPE result_0_1 = -(IQ*UU - QU*IU)*invdet; //-(A(1,0)*A(2,2)-A(1,2)*A(2,0))*invdet;
  MAPTYPE result_1_1 =  (II*UU - IU*IU)*invdet; // (A(0,0)*A(2,2)-A(0,2)*A(2,0))*invdet;
  MAPTYPE result_2_1 = -(II*QU - IQ*IU)*invdet; //-(A(0,0)*A(1,2)-A(1,0)*A(0,2))*invdet;
  MAPTYPE result_0_2 =  (IQ*QU - IU*QQ)*invdet; // (A(1,0)*A(2,1)-A(2,0)*A(1,1))*invdet;
  MAPTYPE result_1_2 = -(II*QU - IU*IQ)*invdet; //-(A(0,0)*A(2,1)-A(2,0)*A(0,1))*invdet;
  MAPTYPE result_2_2 =  (II*QQ - IQ*IQ)*invdet; // (A(0,0)*A(1,1)-A(1,0)*A(0,1))*invdet;

  MAPTYPE ox = *x*result_0_0 + *y*result_1_0 + *z*result_2_0;
  MAPTYPE oy = *x*result_0_1 + *y*result_1_1 + *z*result_2_1;
  MAPTYPE oz = *x*result_0_2 + *y*result_1_2 + *z*result_2_2;

  *x=ox;
  *y=oy;
  *z=oz;

  return( determinant);
}


////////////////////////////////////////////////////////////////////////
// RD12 parameter values

double get_crosspol( char *pixname) {
  if (!strcmp( pixname, "100-1a")) return( 0.0272);
  if (!strcmp( pixname, "100-1b")) return( 0.0293);
  if (!strcmp( pixname, "100-2a")) return( 0.0195);
  if (!strcmp( pixname, "100-2b")) return( 0.0513);
  if (!strcmp( pixname, "100-3a")) return( 0.0521);
  if (!strcmp( pixname, "100-3b")) return( 0.0339);
  if (!strcmp( pixname, "100-4a")) return( 0.0219);
  if (!strcmp( pixname, "100-4b")) return( 0.0402);
  if (!strcmp( pixname, "143-1a")) return( 0.0915);
  if (!strcmp( pixname, "143-1b")) return( 0.0835);
  if (!strcmp( pixname, "143-2a")) return( 0.067);
  if (!strcmp( pixname, "143-2b")) return( 0.0572);
  if (!strcmp( pixname, "143-3a")) return( 0.0875);
  if (!strcmp( pixname, "143-3b")) return( 0.053);
  if (!strcmp( pixname, "143-4a")) return( 0.0357);
  if (!strcmp( pixname, "143-4b")) return( 0.0371);
  if (!strcmp( pixname, "143-5"))  return( 0.882);
  if (!strcmp( pixname, "143-6"))  return( 0.92);
  if (!strcmp( pixname, "143-7"))  return( 0.974);
  if (!strcmp( pixname, "217-1"))  return( 0.926);
  if (!strcmp( pixname, "217-2"))  return( 0.961);
  if (!strcmp( pixname, "217-3"))  return( 0.924);
  if (!strcmp( pixname, "217-4"))  return( 0.916);
  if (!strcmp( pixname, "217-5a")) return( 0.0256);
  if (!strcmp( pixname, "217-5b")) return( 0.0246);
  if (!strcmp( pixname, "217-6a")) return( 0.026);
  if (!strcmp( pixname, "217-6b")) return( 0.0236);
  if (!strcmp( pixname, "217-7a")) return( 0.0307);
  if (!strcmp( pixname, "217-7b")) return( 0.0327);
  if (!strcmp( pixname, "217-8a")) return( 0.0298);
  if (!strcmp( pixname, "217-8b")) return( 0.0303);
  if (!strcmp( pixname, "353-1"))  return( 0.938);
  if (!strcmp( pixname, "353-2"))  return( 0.91);
  if (!strcmp( pixname, "353-3a")) return( 0.0583);
  if (!strcmp( pixname, "353-3b")) return( 0.0413);
  if (!strcmp( pixname, "353-4a")) return( 0.0683);
  if (!strcmp( pixname, "353-4b")) return( 0.0439);
  if (!strcmp( pixname, "353-5a")) return( 0.0829);
  if (!strcmp( pixname, "353-5b")) return( 0.0661);
  if (!strcmp( pixname, "353-6a")) return( 0.0664);
  if (!strcmp( pixname, "353-6b")) return( 0.0595);
  if (!strcmp( pixname, "353-7"))  return( 0.852);
  if (!strcmp( pixname, "353-8"))  return( 0.855);
  if (!strcmp( pixname, "545-1"))  return( 0.913);
  if (!strcmp( pixname, "545-2"))  return( 0.894);
  if (!strcmp( pixname, "545-4"))  return( 0.89);
  if (!strcmp( pixname, "857-1"))  return( 0.887);
  if (!strcmp( pixname, "857-2"))  return( 0.883);
  if (!strcmp( pixname, "857-3"))  return( 0.856);
  if (!strcmp( pixname, "857-4"))  return( 0.897);
  fprintf( stderr, "get_crosspol(): Unkown pixname (%s)", pixname);
  exit( 1);
}

double get_RD12nep( char *pixname) {
  if (!strcmp( pixname, "100-1a")) return( 2.47672e-16);
  if (!strcmp( pixname, "100-1b")) return( 2.5673e-16);
  if (!strcmp( pixname, "100-2a")) return( 1.84486e-16);
  if (!strcmp( pixname, "100-2b")) return( 2.68522e-16);
  if (!strcmp( pixname, "100-3a")) return( 1.45183e-16);
  if (!strcmp( pixname, "100-3b")) return( 1.46603e-16);
  if (!strcmp( pixname, "100-4a")) return( 2.26421e-16);
  if (!strcmp( pixname, "100-4b")) return( 2.45552e-16);
  if (!strcmp( pixname, "143-1a")) return( 1.43495e-16);
  if (!strcmp( pixname, "143-1b")) return( 1.91453e-16);
  if (!strcmp( pixname, "143-2a")) return( 1.34663e-16);
  if (!strcmp( pixname, "143-2b")) return( 1.47483e-16);
  if (!strcmp( pixname, "143-3a")) return( 1.49003e-16);
  if (!strcmp( pixname, "143-3b")) return( 1.29915e-16);
  if (!strcmp( pixname, "143-4a")) return( 1.46059e-16);
  if (!strcmp( pixname, "143-4b")) return( 1.49668e-16);
  if (!strcmp( pixname, "143-5"))  return( 1.81481e-16);
  if (!strcmp( pixname, "143-6"))  return( 1.65717e-16);
  if (!strcmp( pixname, "143-7"))  return( 1.57075e-16);
  if (!strcmp( pixname, "217-1"))  return( 1.52802e-16);
  if (!strcmp( pixname, "217-2"))  return( 1.52707e-16);
  if (!strcmp( pixname, "217-3"))  return( 1.52137e-16);
  if (!strcmp( pixname, "217-4"))  return( 1.39316e-16);
  if (!strcmp( pixname, "217-5a")) return( 1.67331e-16);
  if (!strcmp( pixname, "217-5b")) return( 1.46629e-16);
  if (!strcmp( pixname, "217-6a")) return( 1.65907e-16);
  if (!strcmp( pixname, "217-6b")) return( 1.59924e-16);
  if (!strcmp( pixname, "217-7a")) return( 1.54986e-16);
  if (!strcmp( pixname, "217-7b")) return( 1.46629e-16);
  if (!strcmp( pixname, "217-8a")) return( 1.69801e-16);
  if (!strcmp( pixname, "217-8b")) return( 1.65907e-16);
  if (!strcmp( pixname, "353-1"))  return( 1.39791e-16);
  if (!strcmp( pixname, "353-2"))  return( 1.5584e-16);
  if (!strcmp( pixname, "353-3a")) return( 1.99715e-16);
  if (!strcmp( pixname, "353-3b")) return( 1.63248e-16);
  if (!strcmp( pixname, "353-4a")) return( 1.40076e-16);
  if (!strcmp( pixname, "353-4b")) return( 1.45774e-16);
  if (!strcmp( pixname, "353-5a")) return( 1.5603e-16);
  if (!strcmp( pixname, "353-5b")) return( 1.5717e-16);
  if (!strcmp( pixname, "353-6a")) return( 1.62393e-16);
  if (!strcmp( pixname, "353-6b")) return( 1.40551e-16);
  if (!strcmp( pixname, "353-7")) return(  1.41975e-16);
  if (!strcmp( pixname, "353-8")) return(  1.46629e-16);
  if (!strcmp( pixname, "545-1")) return(  2.06268e-16);
  if (!strcmp( pixname, "545-2")) return(  1.78632e-16);
  if (!strcmp( pixname, "545-4")) return(  1.63153e-16);
  if (!strcmp( pixname, "857-1")) return(  2.17664e-16);
  if (!strcmp( pixname, "857-2")) return(  2.36372e-16);
  if (!strcmp( pixname, "857-3")) return(  2.08547e-16);
  if (!strcmp( pixname, "857-4")) return(  2.01804e-16);
  fprintf( stderr, "get_RD12nep(): Unkown pixname (%s)", pixname);
  exit( 1);
}

double get_DX11calib( char *pixname) {
  if (!strcmp( pixname, "100-1a")) return( 9.98929e-14);
  if (!strcmp( pixname, "100-1b")) return( 1.22488e-13);
  if (!strcmp( pixname, "100-2a")) return( 1.51546e-13);
  if (!strcmp( pixname, "100-2b")) return( 1.58035e-13);
  if (!strcmp( pixname, "100-3a")) return( 1.38029e-13);
  if (!strcmp( pixname, "100-3b")) return( 1.14552e-13);
  if (!strcmp( pixname, "100-4a")) return( 1.46123e-13);
  if (!strcmp( pixname, "100-4b")) return( 1.16753e-13);
  if (!strcmp( pixname, "143-1a")) return( 1.87259e-13);
  if (!strcmp( pixname, "143-1b")) return( 1.60741e-13);
  if (!strcmp( pixname, "143-2a")) return( 1.77637e-13);
  if (!strcmp( pixname, "143-2b")) return( 1.81202e-13);
  if (!strcmp( pixname, "143-3a")) return( 1.78495e-13);
  if (!strcmp( pixname, "143-3b")) return( 1.60738e-13);
  if (!strcmp( pixname, "143-4a")) return( 1.65044e-13);
  if (!strcmp( pixname, "143-4b")) return( 1.54786e-13);
  if (!strcmp( pixname, "143-5"))  return( 2.6395e-13);
  if (!strcmp( pixname, "143-6"))  return( 2.36995e-13);
  if (!strcmp( pixname, "143-7"))  return( 2.56619e-13);
  if (!strcmp( pixname, "217-1"))  return( 1.67389e-13);
  if (!strcmp( pixname, "217-2"))  return( 1.61614e-13);
  if (!strcmp( pixname, "217-3"))  return( 1.70059e-13);
  if (!strcmp( pixname, "217-4"))  return( 1.65354e-13);
  if (!strcmp( pixname, "217-5a")) return( 1.13427e-13);
  if (!strcmp( pixname, "217-5b")) return( 1.14376e-13);
  if (!strcmp( pixname, "217-6a")) return( 1.14989e-13);
  if (!strcmp( pixname, "217-6b")) return( 1.16227e-13);
  if (!strcmp( pixname, "217-7a")) return( 1.23498e-13);
  if (!strcmp( pixname, "217-7b")) return( 1.18642e-13);
  if (!strcmp( pixname, "217-8a")) return( 1.1918e-13);
  if (!strcmp( pixname, "217-8b")) return( 1.13096e-13);
  if (!strcmp( pixname, "353-1"))  return( 5.70346e-14);
  if (!strcmp( pixname, "353-2"))  return( 6.20394e-14);
  if (!strcmp( pixname, "353-3a")) return( 3.48625e-14);
  if (!strcmp( pixname, "353-3b")) return( 3.37349e-14);
  if (!strcmp( pixname, "353-4a")) return( 2.86904e-14);
  if (!strcmp( pixname, "353-4b")) return( 2.83174e-14);
  if (!strcmp( pixname, "353-5a")) return( 3.27368e-14);
  if (!strcmp( pixname, "353-5b")) return( 3.23329e-14);
  if (!strcmp( pixname, "353-6a")) return( 2.34751e-14);
  if (!strcmp( pixname, "353-6b")) return( 2.12574e-14);
  if (!strcmp( pixname, "353-7"))  return( 4.70167e-14);
  if (!strcmp( pixname, "353-8"))  return( 4.42892e-14);
  if (!strcmp( pixname, "545-1"))  return( 3.30836e-16);
  if (!strcmp( pixname, "545-2"))  return( 3.17769e-16);
  if (!strcmp( pixname, "545-4"))  return( 2.70598e-16);
  if (!strcmp( pixname, "857-1"))  return( 3.43811e-16);
  if (!strcmp( pixname, "857-2"))  return( 3.78339e-16);
  if (!strcmp( pixname, "857-3"))  return( 3.25099e-16);
  if (!strcmp( pixname, "857-4"))  return( 2.28563e-16);
  fprintf( stderr, "get_DX11calib(): Unkown pixname (%s)", pixname);
  exit( 1);
}

double get_DX11nep( char *pixname) {
  if (!strcmp( pixname, "100-1a")) return( 3.5671e-16);
  if (!strcmp( pixname, "100-1b")) return( 3.62753e-16);
  if (!strcmp( pixname, "100-2a")) return( 2.69412e-16);
  if (!strcmp( pixname, "100-2b")) return( 3.9512e-16);
  if (!strcmp( pixname, "100-3a")) return( 2.06693e-16);
  if (!strcmp( pixname, "100-3b")) return( 2.08439e-16);
  if (!strcmp( pixname, "100-4a")) return( 3.24477e-16);
  if (!strcmp( pixname, "100-4b")) return( 3.53218e-16);
  if (!strcmp( pixname, "143-1a")) return( 2.00246e-16);
  if (!strcmp( pixname, "143-1b")) return( 2.62832e-16);
  if (!strcmp( pixname, "143-2a")) return( 1.8695e-16);
  if (!strcmp( pixname, "143-2b")) return( 2.04678e-16);
  if (!strcmp( pixname, "143-3a")) return( 2.08707e-16);
  if (!strcmp( pixname, "143-3b")) return( 1.78489e-16);
  if (!strcmp( pixname, "143-4a")) return( 2.03469e-16);
  if (!strcmp( pixname, "143-4b")) return( 2.07096e-16);
  if (!strcmp( pixname, "143-5"))  return( 2.52759e-16);
  if (!strcmp( pixname, "143-6"))  return( 2.30867e-16);
  if (!strcmp( pixname, "143-7"))  return( 2.1878e-16);
  if (!strcmp( pixname, "217-1"))  return( 2.14214e-16);
  if (!strcmp( pixname, "217-2"))  return( 2.13811e-16);
  if (!strcmp( pixname, "217-3"))  return( 2.12468e-16);
  if (!strcmp( pixname, "217-4"))  return( 1.93128e-16);
  if (!strcmp( pixname, "217-5a")) return( 2.33016e-16);
  if (!strcmp( pixname, "217-5b")) return( 2.05618e-16);
  if (!strcmp( pixname, "217-6a")) return( 2.31136e-16);
  if (!strcmp( pixname, "217-6b")) return( 2.21332e-16);
  if (!strcmp( pixname, "217-7a")) return( 2.13408e-16);
  if (!strcmp( pixname, "217-7b")) return( 2.02126e-16);
  if (!strcmp( pixname, "217-8a")) return( 2.34091e-16);
  if (!strcmp( pixname, "217-8b")) return( 2.3033e-16);
  if (!strcmp( pixname, "353-1"))  return( 1.95277e-16);
  if (!strcmp( pixname, "353-2"))  return( 2.169e-16);
  if (!strcmp( pixname, "353-3a")) return( 2.76396e-16);
  if (!strcmp( pixname, "353-3b")) return( 2.28181e-16);
  if (!strcmp( pixname, "353-4a")) return( 1.95143e-16);
  if (!strcmp( pixname, "353-4b")) return( 2.02798e-16);
  if (!strcmp( pixname, "353-5a")) return( 2.17034e-16);
  if (!strcmp( pixname, "353-5b")) return( 2.18646e-16);
  if (!strcmp( pixname, "353-6a")) return( 2.25092e-16);
  if (!strcmp( pixname, "353-6b")) return( 1.95411e-16);
  if (!strcmp( pixname, "353-7"))  return( 1.97829e-16);
  if (!strcmp( pixname, "353-8"))  return( 2.03201e-16);
  if (!strcmp( pixname, "545-1"))  return( 3.00571e-16);
  if (!strcmp( pixname, "545-2"))  return( 2.55311e-16);
  if (!strcmp( pixname, "545-4"))  return( 2.41477e-16);
  if (!strcmp( pixname, "857-1"))  return( 3.06749e-16);
  if (!strcmp( pixname, "857-2"))  return( 3.21925e-16);
  if (!strcmp( pixname, "857-3"))  return( 2.87006e-16);
  if (!strcmp( pixname, "857-4"))  return( 2.80291e-16);
  fprintf( stderr, "get_DX11nep(): Unkown pixname (%s)", pixname);
  exit( 1);
}

double get_FFP10PSM_offset( char *pixname) {

/*
# mean of destriping offsets from one sroll run per freq with rstep=10, gainstep=1, no Theo_noPS
# (magique4:/wrk/symottet/srollex_JAN17v1/run_stim/psmproj*)
import numpy
from IMO_4_27 import *
for freq in ["100", "143", "217", "353"]:
  detset = DETSETS["%sghz"%freq]
  for pixname in detset:
    offsets = numpy.fromfile( "/pscratch1/RD12_data/dmc/MISS03/DATA/stimDEC16v1_VEC/psmprojrstep10_%sghz_%s_offset_OFF" % (freq, pixname))
    offsets[offsets==0] = numpy.nan
    offsets[offsets<-5000] = numpy.nan
    print( '  if (!strcmp( pixname, "%s")) return( %g);' % (pixname, numpy.nanmean( offsets)))
*/

  if (DATASET == DEC16) {
    if (!strcmp( pixname, "100-1a")) return( -4.03765e-07);
    if (!strcmp( pixname, "100-1b")) return( -2.72012e-07);
    if (!strcmp( pixname, "100-4a")) return(  3.2449e-07);
    if (!strcmp( pixname, "100-4b")) return(  3.04028e-07);
    if (!strcmp( pixname, "100-2a")) return(  1.11234e-07);
    if (!strcmp( pixname, "100-2b")) return(  2.07847e-08);
    if (!strcmp( pixname, "100-3a")) return(  2.5425e-07);
    if (!strcmp( pixname, "100-3b")) return( -3.12475e-07);
    if (!strcmp( pixname, "143-1a")) return( -6.01675e-07);
    if (!strcmp( pixname, "143-1b")) return( -1.74636e-07);
    if (!strcmp( pixname, "143-3a")) return( -1.49019e-06);
    if (!strcmp( pixname, "143-3b")) return( -6.39797e-07);
    if (!strcmp( pixname, "143-2a")) return( -5.07831e-07);
    if (!strcmp( pixname, "143-2b")) return( -7.02835e-08);
    if (!strcmp( pixname, "143-4a")) return(  2.85206e-07);
    if (!strcmp( pixname, "143-4b")) return( -2.93546e-07);
    if (!strcmp( pixname, "143-5"))  return(  1.27254e-06);
    if (!strcmp( pixname, "143-6"))  return(  3.27276e-07);
    if (!strcmp( pixname, "143-7"))  return(  1.49325e-06);
    if (!strcmp( pixname, "217-5a")) return( -1.87701e-06);
    if (!strcmp( pixname, "217-5b")) return( -1.60209e-06);
    if (!strcmp( pixname, "217-7a")) return( -1.40219e-06);
    if (!strcmp( pixname, "217-7b")) return( -2.368e-06);
    if (!strcmp( pixname, "217-6a")) return( -1.5309e-06);
    if (!strcmp( pixname, "217-6b")) return( -1.56646e-06);
    if (!strcmp( pixname, "217-8a")) return( -1.61019e-06);
    if (!strcmp( pixname, "217-8b")) return( -1.47054e-06);
    if (!strcmp( pixname, "217-1"))  return(  2.96939e-06);
    if (!strcmp( pixname, "217-2"))  return(  4.17944e-06);
    if (!strcmp( pixname, "217-3"))  return(  3.56126e-06);
    if (!strcmp( pixname, "217-4"))  return(  2.62157e-06);
    if (!strcmp( pixname, "353-3a")) return( -1.40848e-05);
    if (!strcmp( pixname, "353-3b")) return( -1.26097e-05);
    if (!strcmp( pixname, "353-5a")) return( -2.42422e-05);
    if (!strcmp( pixname, "353-5b")) return( -2.28664e-05);
    if (!strcmp( pixname, "353-4a")) return(  2.06008e-05);
    if (!strcmp( pixname, "353-4b")) return(  2.16748e-05);
    if (!strcmp( pixname, "353-6a")) return( -8.42712e-06);
    if (!strcmp( pixname, "353-6b")) return( -5.81768e-05);
    if (!strcmp( pixname, "353-1"))  return( -3.94025e-06);
    if (!strcmp( pixname, "353-2"))  return(  3.42585e-06);
    if (!strcmp( pixname, "353-7"))  return(  3.67164e-05);
    if (!strcmp( pixname, "353-8"))  return(  6.16436e-05);
  }
  else if (DATASET == DUST17) {
    if (!strcmp( pixname, "100-1a")) return( -2.19821e-07);
    if (!strcmp( pixname, "100-1b")) return( -1.81266e-07);
    if (!strcmp( pixname, "100-4a")) return(  2.05915e-07);
    if (!strcmp( pixname, "100-4b")) return(  1.89824e-07);
    if (!strcmp( pixname, "100-2a")) return(  7.22098e-08);
    if (!strcmp( pixname, "100-2b")) return(  5.93942e-09);
    if (!strcmp( pixname, "100-3a")) return(  1.70552e-07);
    if (!strcmp( pixname, "100-3b")) return( -2.26818e-07);
    if (!strcmp( pixname, "143-1a")) return( -3.7996e-07);
    if (!strcmp( pixname, "143-1b")) return( -1.19775e-07);
    if (!strcmp( pixname, "143-3a")) return( -9.51305e-07);
    if (!strcmp( pixname, "143-3b")) return( -4.23685e-07);
    if (!strcmp( pixname, "143-2a")) return( -3.60709e-07);
    if (!strcmp( pixname, "143-2b")) return( -6.7952e-08);
    if (!strcmp( pixname, "143-4a")) return(  1.88314e-07);
    if (!strcmp( pixname, "143-4b")) return( -2.20707e-07);
    if (!strcmp( pixname, "143-5"))  return(  8.61606e-07);
    if (!strcmp( pixname, "143-6"))  return(  2.01754e-07);
    if (!strcmp( pixname, "143-7"))  return(  1.01037e-06);
    if (!strcmp( pixname, "217-5a")) return( -1.09928e-06);
    if (!strcmp( pixname, "217-5b")) return( -9.45949e-07);
    if (!strcmp( pixname, "217-7a")) return( -8.20351e-07);
    if (!strcmp( pixname, "217-7b")) return( -1.51385e-06);
    if (!strcmp( pixname, "217-6a")) return( -8.94687e-07);
    if (!strcmp( pixname, "217-6b")) return( -1.04652e-06);
    if (!strcmp( pixname, "217-8a")) return( -1.04274e-06);
    if (!strcmp( pixname, "217-8b")) return( -8.69178e-07);
    if (!strcmp( pixname, "217-1"))  return(  1.82494e-06);
    if (!strcmp( pixname, "217-2"))  return(  2.22355e-06);
    if (!strcmp( pixname, "217-3"))  return(  2.31029e-06);
    if (!strcmp( pixname, "217-4"))  return(  1.81177e-06);
    if (!strcmp( pixname, "353-3a")) return( -4.19052e-06);
    if (!strcmp( pixname, "353-3b")) return( -8.34711e-06);
    if (!strcmp( pixname, "353-5a")) return( -1.12343e-05);
    if (!strcmp( pixname, "353-5b")) return( -1.41981e-05);
    if (!strcmp( pixname, "353-4a")) return(  1.11967e-05);
    if (!strcmp( pixname, "353-4b")) return(  9.6841e-06);
    if (!strcmp( pixname, "353-6a")) return( -4.64339e-06);
    if (!strcmp( pixname, "353-6b")) return( -2.81353e-05);
    if (!strcmp( pixname, "353-1"))  return( -2.15907e-06);
    if (!strcmp( pixname, "353-2"))  return(  1.90524e-06);
    if (!strcmp( pixname, "353-7"))  return(  2.01955e-05);
    if (!strcmp( pixname, "353-8"))  return(  2.97797e-05);
  }
  fprintf( stderr, "get_FFP10PSM_offset(): Unkown pixname (%s)", pixname);
  exit( 1);
}

////////////////////////////////////////////////////////////////////////
// bandpass correction coeff from the mean of 15 DEC16v1 iterations

double get_dust_coef( char *pixname) {
  if (DATASET == DEC16) {
    if (!strcmp( pixname, "100-1a")) return( -0.000498288);
    if (!strcmp( pixname, "100-1b")) return( -0.000221353);
    if (!strcmp( pixname, "100-4a")) return(  0.000285105);
    if (!strcmp( pixname, "100-4b")) return(  0.000308142);
    if (!strcmp( pixname, "100-2a")) return(  9.78852e-05);
    if (!strcmp( pixname, "100-2b")) return(  4.10159e-05);
    if (!strcmp( pixname, "100-3a")) return(  0.000177145);
    if (!strcmp( pixname, "100-3b")) return( -0.000189653);
    if (!strcmp( pixname, "143-1a")) return( -0.000516647);
    if (!strcmp( pixname, "143-1b")) return( -9.60593e-05);
    if (!strcmp( pixname, "143-3a")) return( -0.00127784);
    if (!strcmp( pixname, "143-3b")) return( -0.000496937);
    if (!strcmp( pixname, "143-2a")) return( -0.000331263);
    if (!strcmp( pixname, "143-2b")) return(  2.79656e-05);
    if (!strcmp( pixname, "143-4a")) return(  0.000275303);
    if (!strcmp( pixname, "143-4b")) return( -0.000159278);
    if (!strcmp( pixname, "143-5"))  return(  0.00101935);
    if (!strcmp( pixname, "143-6"))  return(  0.000352608);
    if (!strcmp( pixname, "143-7"))  return(  0.0012028);
    if (!strcmp( pixname, "217-5a")) return( -0.00201468);
    if (!strcmp( pixname, "217-5b")) return( -0.00157757);
    if (!strcmp( pixname, "217-7a")) return( -0.00144063);
    if (!strcmp( pixname, "217-7b")) return( -0.00218591);
    if (!strcmp( pixname, "217-6a")) return( -0.00165795);
    if (!strcmp( pixname, "217-6b")) return( -0.00128935);
    if (!strcmp( pixname, "217-8a")) return( -0.00140382);
    if (!strcmp( pixname, "217-8b")) return( -0.00150728);
    if (!strcmp( pixname, "217-1"))  return(  0.00288672);
    if (!strcmp( pixname, "217-2"))  return(  0.00496127);
    if (!strcmp( pixname, "217-3"))  return(  0.00315886);
    if (!strcmp( pixname, "217-4"))  return(  0.00207034);
    if (!strcmp( pixname, "353-3a")) return( -0.0240008);
    if (!strcmp( pixname, "353-3b")) return( -0.0106631);
    if (!strcmp( pixname, "353-5a")) return( -0.0320224);
    if (!strcmp( pixname, "353-5b")) return( -0.0225553);
    if (!strcmp( pixname, "353-4a")) return(  0.0246175);
    if (!strcmp( pixname, "353-4b")) return(  0.0303393);
    if (!strcmp( pixname, "353-6a")) return( -0.0094937);
    if (!strcmp( pixname, "353-6b")) return( -0.0753309);
    if (!strcmp( pixname, "353-1"))  return( -0.0043237);
    if (!strcmp( pixname, "353-2"))  return(  0.00351919);
    if (!strcmp( pixname, "353-7"))  return(  0.0409882);
    if (!strcmp( pixname, "353-8"))  return(  0.0789257);
  }
  else if (DATASET == DUST17) {
    if (!strcmp( pixname, "100-1a")) return( -0.000442072);
    if (!strcmp( pixname, "100-1b")) return( -0.000216695);
    if (!strcmp( pixname, "100-4a")) return(  0.000298929);
    if (!strcmp( pixname, "100-4b")) return(  0.00024499);
    if (!strcmp( pixname, "100-2a")) return(  0.000104405);
    if (!strcmp( pixname, "100-2b")) return(  3.75277e-05);
    if (!strcmp( pixname, "100-3a")) return(  0.000159839);
    if (!strcmp( pixname, "100-3b")) return( -0.000186923);
    if (!strcmp( pixname, "143-1a")) return( -0.000479022);
    if (!strcmp( pixname, "143-1b")) return( -9.98968e-05);
    if (!strcmp( pixname, "143-3a")) return( -0.00118185);
    if (!strcmp( pixname, "143-3b")) return( -0.000444634);
    if (!strcmp( pixname, "143-2a")) return( -0.000306351);
    if (!strcmp( pixname, "143-2b")) return(  2.91494e-05);
    if (!strcmp( pixname, "143-4a")) return(  0.000237514);
    if (!strcmp( pixname, "143-4b")) return( -0.000142312);
    if (!strcmp( pixname, "143-5"))  return(  0.000937954);
    if (!strcmp( pixname, "143-6"))  return(  0.00033593);
    if (!strcmp( pixname, "143-7"))  return(  0.00111352);
    if (!strcmp( pixname, "217-5a")) return( -0.00187856);
    if (!strcmp( pixname, "217-5b")) return( -0.00144737);
    if (!strcmp( pixname, "217-7a")) return( -0.00134755);
    if (!strcmp( pixname, "217-7b")) return( -0.00203825);
    if (!strcmp( pixname, "217-6a")) return( -0.00154367);
    if (!strcmp( pixname, "217-6b")) return( -0.00118794);
    if (!strcmp( pixname, "217-8a")) return( -0.00126437);
    if (!strcmp( pixname, "217-8b")) return( -0.00141228);
    if (!strcmp( pixname, "217-1"))  return(  0.00267427);
    if (!strcmp( pixname, "217-2"))  return(  0.0046019);
    if (!strcmp( pixname, "217-3"))  return(  0.00292752);
    if (!strcmp( pixname, "217-4"))  return(  0.00191631);
    if (!strcmp( pixname, "353-3a")) return( -0.0219551);
    if (!strcmp( pixname, "353-3b")) return( -0.00910228);
    if (!strcmp( pixname, "353-5a")) return( -0.0288053);
    if (!strcmp( pixname, "353-5b")) return( -0.0206396);
    if (!strcmp( pixname, "353-4a")) return(  0.0220797);
    if (!strcmp( pixname, "353-4b")) return(  0.0278787);
    if (!strcmp( pixname, "353-6a")) return( -0.00920003);
    if (!strcmp( pixname, "353-6b")) return( -0.0686782);
    if (!strcmp( pixname, "353-1"))  return( -0.00395162);
    if (!strcmp( pixname, "353-2"))  return(  0.0033748);
    if (!strcmp( pixname, "353-7"))  return(  0.0374276);
    if (!strcmp( pixname, "353-8"))  return(  0.0715714);
  }
  return( 0.0);
}

double get_co_coef( char *pixname) {
  if (DATASET == DEC16) {
    if (!strcmp( pixname, "100-1a")) return( -2.74154);
    if (!strcmp( pixname, "100-1b")) return( -1.41745);
    if (!strcmp( pixname, "100-4a")) return(  4.15712);
    if (!strcmp( pixname, "100-4b")) return(  1.44879);
    if (!strcmp( pixname, "100-2a")) return(  0.38489);
    if (!strcmp( pixname, "100-2b")) return( -1.81044);
    if (!strcmp( pixname, "100-3a")) return(  2.01032);
    if (!strcmp( pixname, "100-3b")) return( -2.03168);
    if (!strcmp( pixname, "217-5a")) return( -0.271675);
    if (!strcmp( pixname, "217-5b")) return( -0.503364);
    if (!strcmp( pixname, "217-7a")) return(  0.552799);
    if (!strcmp( pixname, "217-7b")) return( -0.466202);
    if (!strcmp( pixname, "217-6a")) return( -2.89154);
    if (!strcmp( pixname, "217-6b")) return( -1.99609);
    if (!strcmp( pixname, "217-8a")) return( -0.298847);
    if (!strcmp( pixname, "217-8b")) return( -1.56593);
    if (!strcmp( pixname, "217-1"))  return(  3.27367);
    if (!strcmp( pixname, "217-2"))  return( -1.37482);
    if (!strcmp( pixname, "217-3"))  return(  3.77113);
    if (!strcmp( pixname, "217-4"))  return(  1.77087);
    if (!strcmp( pixname, "353-3a")) return(  3.59778);
    if (!strcmp( pixname, "353-3b")) return(  7.75699);
    if (!strcmp( pixname, "353-5a")) return( -4.88332);
    if (!strcmp( pixname, "353-5b")) return( -3.49354);
    if (!strcmp( pixname, "353-4a")) return(  0.971087);
    if (!strcmp( pixname, "353-4b")) return( -8.46291);
    if (!strcmp( pixname, "353-6a")) return( -7.0019);
    if (!strcmp( pixname, "353-6b")) return( -1.81139);
    if (!strcmp( pixname, "353-1"))  return( -0.858896);
    if (!strcmp( pixname, "353-2"))  return(  1.41224);
    if (!strcmp( pixname, "353-7"))  return(  7.99562);
    if (!strcmp( pixname, "353-8"))  return(  4.77823);
  }
  else if (DATASET == DUST17) {
    if (!strcmp( pixname, "100-1a")) return( -2.71355);
    if (!strcmp( pixname, "100-1b")) return( -1.41029);
    if (!strcmp( pixname, "100-4a")) return(  4.1554);
    if (!strcmp( pixname, "100-4b")) return(  1.42854);
    if (!strcmp( pixname, "100-2a")) return(  0.35614);
    if (!strcmp( pixname, "100-2b")) return( -1.81706);
    if (!strcmp( pixname, "100-3a")) return(  2.05138);
    if (!strcmp( pixname, "100-3b")) return( -2.05056);
    if (!strcmp( pixname, "217-5a")) return( -0.191966);
    if (!strcmp( pixname, "217-5b")) return( -0.479777);
    if (!strcmp( pixname, "217-7a")) return(  0.627107);
    if (!strcmp( pixname, "217-7b")) return( -0.445939);
    if (!strcmp( pixname, "217-6a")) return( -2.83853);
    if (!strcmp( pixname, "217-6b")) return( -1.99354);
    if (!strcmp( pixname, "217-8a")) return( -0.304437);
    if (!strcmp( pixname, "217-8b")) return( -1.54099);
    if (!strcmp( pixname, "217-1"))  return(  3.19979);
    if (!strcmp( pixname, "217-2"))  return( -1.46217);
    if (!strcmp( pixname, "217-3"))  return(  3.70334);
    if (!strcmp( pixname, "217-4"))  return(  1.72711);
    if (!strcmp( pixname, "353-3a")) return(  4.09424);
    if (!strcmp( pixname, "353-3b")) return(  7.73603);
    if (!strcmp( pixname, "353-5a")) return( -4.36325);
    if (!strcmp( pixname, "353-5b")) return( -3.14017);
    if (!strcmp( pixname, "353-4a")) return(  0.809547);
    if (!strcmp( pixname, "353-4b")) return( -9.07817);
    if (!strcmp( pixname, "353-6a")) return( -6.50679);
    if (!strcmp( pixname, "353-6b")) return( -0.731458);
    if (!strcmp( pixname, "353-1"))  return( -0.822329);
    if (!strcmp( pixname, "353-2"))  return(  1.3064);
    if (!strcmp( pixname, "353-7"))  return(  7.36087);
    if (!strcmp( pixname, "353-8"))  return(  3.33509);
  }
  return( 0.0);
}

double get_freefree_coef( char *pixname) {
  if (DATASET == DEC16) {
    if (!strcmp( pixname, "100-1a")) return(  0.00600425);
    if (!strcmp( pixname, "100-1b")) return(  0.00152695);
    if (!strcmp( pixname, "100-4a")) return( -0.00415516);
    if (!strcmp( pixname, "100-4b")) return( -0.00386588);
    if (!strcmp( pixname, "100-2a")) return( -0.000272462);
    if (!strcmp( pixname, "100-2b")) return(  0.000555563);
    if (!strcmp( pixname, "100-3a")) return( -0.00233568);
    if (!strcmp( pixname, "100-3b")) return(  0.00254241);
  }
  else if (DATASET == DUST17) {
    if (!strcmp( pixname, "100-1a")) return(  0.00406279);
    if (!strcmp( pixname, "100-1b")) return(  0.00228952);
    if (!strcmp( pixname, "100-4a")) return( -0.00367969);
    if (!strcmp( pixname, "100-4b")) return( -0.00384853);
    if (!strcmp( pixname, "100-2a")) return( -0.00180881);
    if (!strcmp( pixname, "100-2b")) return(  0.00292988);
    if (!strcmp( pixname, "100-3a")) return( -0.00223639);
    if (!strcmp( pixname, "100-3b")) return(  0.00229123);
  }
  return( 0.0);
}

// deprecated, replaced by mean of 15 iterations
void get_DEC16_bandpass_correction( char *pixname, double *dust_coeff, double *freefree_coeff, double *co_coeff) {

  int dmc_err;
  int bidx = -1;
  int nbolo = -1;
  int ndust = 0, nco = 0, nfree = 0; // 1 if dust, co, freefree map factors in X2 vects, else 0
  int gainstep = 0; // sroll GAINSTEP param
  int npixbeam = 0; // number of Theo_noPS param per bolometer

// sroll.c:
// fprintf(stderr,"Write MAT  %lld\n",(long long) PIOWriteVECT(saveg,x3+newnr[nbolo],0,(nbolo*(GAINSTEP+npixbeam)+nmatdust+nmatco+nfreefree)*sizeof(PIODOUBLE)));

  char X2name[STRLEN];
  if (!strncmp( pixname, "100", 3)) {
    sprintf( X2name, "%s/stimDEC16v1_VEC/DEC16v1_100ghz_100-1a_offset_0_X2", dbpath);
    nbolo = 8;
    gainstep = 128;
    npixbeam = 4;
    ndust = 1;
    nco = 1;
    nfree = 1;
  }
  if (!strncmp( pixname, "143", 3)) {
    sprintf( X2name, "%s/stimDEC16v1_VEC/DEC16v1_143ghz_143-1a_offset_0_X2", dbpath);
    nbolo = 11;
    gainstep = 128;
    npixbeam = 4;
    ndust = 1;
    nco = 0;
    nfree = 0;
  }
  if (!strncmp( pixname, "217", 3)) {
    sprintf( X2name, "%s/stimDEC16v1_VEC/DEC16v1_217ghz_217-5a_offset_0_X2", dbpath);
    nbolo = 12;
    gainstep = 128;
    npixbeam = 4;
    ndust = 1;
    nco = 1;
    nfree = 0;
  }
  if (!strncmp( pixname, "353", 3)) {
    sprintf( X2name, "%s/stimDEC16v1_VEC/DEC16v1_353ghz_353-3a_offset_0_X2", dbpath);
    nbolo = 12;
    gainstep = 32;
    npixbeam = 7;
    ndust = 1;
    nco = 1;
    nfree = 0;
  }
  assert( nbolo != -1);

  // index of bolometer in parameter file and output vects
  if (!strcmp( pixname, "100-1a")) bidx = 0;
  if (!strcmp( pixname, "100-1b")) bidx = 1;
  if (!strcmp( pixname, "100-2a")) bidx = 4;
  if (!strcmp( pixname, "100-2b")) bidx = 5;
  if (!strcmp( pixname, "100-3a")) bidx = 6;
  if (!strcmp( pixname, "100-3b")) bidx = 7;
  if (!strcmp( pixname, "100-4a")) bidx = 2;
  if (!strcmp( pixname, "100-4b")) bidx = 3;
  if (!strcmp( pixname, "143-1a")) bidx = 0;
  if (!strcmp( pixname, "143-1b")) bidx = 1;
  if (!strcmp( pixname, "143-2a")) bidx = 4;
  if (!strcmp( pixname, "143-2b")) bidx = 5;
  if (!strcmp( pixname, "143-3a")) bidx = 2;
  if (!strcmp( pixname, "143-3b")) bidx = 3;
  if (!strcmp( pixname, "143-4a")) bidx = 6;
  if (!strcmp( pixname, "143-4b")) bidx = 7;
  if (!strcmp( pixname, "143-5"))  bidx = 8;
  if (!strcmp( pixname, "143-6"))  bidx = 9;
  if (!strcmp( pixname, "143-7"))  bidx = 10;
  if (!strcmp( pixname, "217-5a")) bidx = 0;
  if (!strcmp( pixname, "217-5b")) bidx = 1;
  if (!strcmp( pixname, "217-6a")) bidx = 4;
  if (!strcmp( pixname, "217-6b")) bidx = 5;
  if (!strcmp( pixname, "217-7a")) bidx = 2;
  if (!strcmp( pixname, "217-7b")) bidx = 3;
  if (!strcmp( pixname, "217-8a")) bidx = 6;
  if (!strcmp( pixname, "217-8b")) bidx = 7;
  if (!strcmp( pixname, "217-1"))  bidx = 8;
  if (!strcmp( pixname, "217-2"))  bidx = 9;
  if (!strcmp( pixname, "217-3"))  bidx = 10;
  if (!strcmp( pixname, "217-4"))  bidx = 11;
  if (!strcmp( pixname, "353-3a")) bidx = 0;
  if (!strcmp( pixname, "353-3b")) bidx = 1;
  if (!strcmp( pixname, "353-4a")) bidx = 4;
  if (!strcmp( pixname, "353-4b")) bidx = 5;
  if (!strcmp( pixname, "353-5a")) bidx = 2;
  if (!strcmp( pixname, "353-5b")) bidx = 3;
  if (!strcmp( pixname, "353-6a")) bidx = 6;
  if (!strcmp( pixname, "353-6b")) bidx = 7;
  if (!strcmp( pixname, "353-1"))  bidx = 8;
  if (!strcmp( pixname, "353-2"))  bidx = 9;
  if (!strcmp( pixname, "353-7"))  bidx = 10;
  if (!strcmp( pixname, "353-8"))  bidx = 11;

  int X2len = nbolo * (gainstep + npixbeam + ndust + nco + nfree);
  double *X2 = malloc( X2len * sizeof( double));
  assert( X2 != NULL);
  
  dmc_err = noDMC_readObject_PIODOUBLE( X2name, 0, X2len, X2);
  if (dmc_err < 0) {
    fprintf( stderr, "Error %d when reading %s\n", dmc_err, X2name);
    exit( -1);
  }

  int coefoffset = nbolo * (gainstep + npixbeam);
  
  if (nco > 0) {
    *co_coeff = X2[coefoffset + bidx];
  }
  else {
    *co_coeff = 0.0;
  }

  *dust_coeff = X2[coefoffset + nbolo * nco + bidx];

  if (nfree > 0) {
    *freefree_coeff = X2[coefoffset + nbolo * (nco + ndust) + bidx];
  }
  else {
    *freefree_coeff = 0.0;
  }
  
  free( X2);
}


void get_julia_bandpass_correction( char *pixname, double *dust_coeff, double *freefree_coeff, double *co_coeff) {
  if (!strncmp( pixname, "100", 3)) {
    *co_coeff       = 1.7e-3;
    *dust_coeff     = 2.0e-3;
    *freefree_coeff = 2.5e-3;
    return;
  }
  if (!strncmp( pixname, "143", 3)) {
    *co_coeff       = 0.0;
    *dust_coeff     = 9.0e-4;
    *freefree_coeff = 0.0;
    return;
  }
  if (!strncmp( pixname, "217", 3)) {
    *co_coeff       = 1.0e-3;
    *dust_coeff     = 1.0e-4;
    *freefree_coeff = 0.0;
    return;
  }
  if (!strncmp( pixname, "353", 3)) {
    *co_coeff       = 0.022;
    *dust_coeff     = 0.0036;
    *freefree_coeff = 0.0;
    return;
  }
  fprintf( stderr, "get_julia_bandpass_correction(): unknown frequency for %s\n", pixname);
  exit( -1);
}
