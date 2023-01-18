// this is for srand48(), drand48(), realpath()
#define _XOPEN_SOURCE 500

#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <inttypes.h>
#include <math.h>

#include "HPR_Cleaner.h"
#include "spline.h"       
#include <stdlib.h>
#include "mpi.h"

#include "no_dmc_metadata.h"
#include "no_dmc_data_access.h"
#include "no_dmc_piolib_type_def.h"
#include "no_dmc_util.h"
#include "no_dmc_debug.h"
#include "no_dmc_version.h"

#include "hprclean_param.h"
#include "hprclean_parLoader.h"

int main(int argc,char *argv[])
{
  int ir,i;

  MPI_Init(&argc, &argv);

  int rank, size;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank);
  MPI_Comm_size( MPI_COMM_WORLD, &size);

  /* read params */
  hprclean_parContent Param;

  int res = hprclean_readParam(&Param, argv[1] );
  if (res != 0) {
    fprintf(stderr, "Unable to parse the parameter file.\n");
    exit(res);
  }

  if (rank >= Param.n_Signal) {
    fprintf( stderr, "No parameters for %d MPI rank, exiting\n", rank);
    exit( 0);
  }

  PIOLONG ridx[30000];

  PIOLONG nsa = noDMC_readObject_PIOLONG( Param.BEGINRINGINDEX, 0, 30000, ridx);
  if (nsa<0) {
    fprintf(stderr,"error reading %s\n", Param.BEGINRINGINDEX);
    exit(-1);
  }

  FILE *fout = NULL;
  if (Param.flag_hpr == _PAR_TRUE) {
    fout=fopen(Param.hpr[rank],"w");
  }
  FILE *fout2=fopen(Param.hpr_corr[rank],"w");
  
  for (ir=Param.BeginRing;ir<=Param.EndRing;ir++) {
    int ndata=ridx[ir+1]-ridx[ir];
    
    PIOINT   *hidx = (PIOINT *) malloc(sizeof(PIOINT)*ndata);
    PIOFLOAT *data = (PIOFLOAT *) malloc(sizeof(PIOFLOAT)*ndata);
    PIOFLOAT *zodi = (PIOFLOAT *) malloc(sizeof(PIOFLOAT)*ndata);

    nsa = noDMC_readObject_PIOINT(Param.HIDX[rank],ridx[ir],ndata,hidx);
    if (nsa<0) {
      fprintf(stderr,"error reading %s ring %d\n",Param.HIDX[rank], ir);
      exit( -1);
    }
    nsa = noDMC_readObject_PIOFLOAT(Param.Signal[rank],ridx[ir],ndata,data);
    if (nsa<0) {
      fprintf(stderr,"error reading %s ring %d\n",Param.Signal[rank], ir);
      exit( -1);
    }
    if (Param.flag_zodi == _PAR_TRUE) {
      nsa = noDMC_readObject_PIOFLOAT(Param.zodi[rank],ridx[ir],ndata,zodi);
      if (nsa<0) {
        fprintf(stderr,"error reading %s ring %d\n",Param.zodi[rank], ir);
        exit(-1);
      }
    }
    else {
      // if input zodi TOI omitted, just add zeros to <data>...
      memset( zodi, 0, sizeof(PIOFLOAT) * ndata);
    }

    float hpr[RGSIZE];
    float hpr_corr[RGSIZE];

    //============================================================
    //   TRAVAIL EN INTERNE EN DOUBLE

    double dhpr[RGSIZE];
    double dnhpr[RGSIZE];

    memset( dhpr,     0, RGSIZE * sizeof(double)); // output HPR classic in float64
    memset( dnhpr,    0, RGSIZE * sizeof(double)); // number of hits per dhpr pixel
    memset( hpr_corr, 0, RGSIZE * sizeof(float));  // output HPR with spline

    int npoint=0;
    for (i=0;i<ndata;i++) {
      if (hidx[i]>=0) {
        int iii=1;
        while (isfinite(zodi[i])==0) {
          // propagate zodi values to fill NaN
          fprintf(stderr,"NaN detected in %s Ring=%d offset %d\n",Param.zodi[rank],ir,i+iii);
          if (i>iii-1) zodi[i]=zodi[i-iii];
          else zodi[i]=zodi[i+iii];
          iii++;
        }
        data[i]=(data[i]/Param.calib[rank]-zodi[i]); // To avoid numerical problem
        if (isfinite(data[i])==0) {
          fprintf(stderr,"DATA PROBLEM %s %f %f\n",Param.Signal[rank],data[i],zodi[i]);
          exit(-1);
        }
        npoint++;
      }
    }

    // print HPR_Cleaner message only when ring number is divisible by 1K
    char *msg = NULL;
    char message[128];
    sprintf(message,"rank %d, ring %d, ", rank, ir);
    if ((ir%1000 == 0) || (ir==Param.BeginRing) || (ir==Param.EndRing)) {
      msg = message;
    }
    // clean HPR from splined 1/f noise
    int test = HPR_Cleaner( hidx, data, hpr_corr, ndata, msg);
    for (i=0;i<RGSIZE;i++) {
      hpr_corr[i] = hpr_corr[i] * Param.calib[rank];
    }
    // write splined HPR to disk
    fseek(fout2, (ir)*RGSIZE*sizeof(float),SEEK_SET);
    fwrite(hpr_corr, RGSIZE*sizeof(float),1,fout2);

    if (Param.flag_hpr == _PAR_TRUE) {
      if ((npoint>100) && (test==0)) {
        // produce legacy HPR
        for (i=0;i<ndata;i++) {
          if (hidx[i]>=0) {
            dhpr[hidx[i]]+=data[i];
            dnhpr[hidx[i]]+=1;
          }
        }
        for (i=0;i<RGSIZE;i++) {
          if (dnhpr[i]>0) dhpr[i]=dhpr[i]/dnhpr[i];
        }      
      }

      for (i=0;i<RGSIZE;i++) {
        hpr[i] = dhpr[i] * Param.calib[rank];
      }
      // write legacy HPR to disk
      fseek(fout, (ir)*RGSIZE*sizeof(float),SEEK_SET);
      fwrite(hpr, RGSIZE*sizeof(float),1,fout);
    }

    free(hidx);
    free(data);
    free(zodi);
  }

  if (Param.flag_hpr == _PAR_TRUE) {
    fclose(fout);
  }
  fclose(fout2);
  
  MPI_Barrier(MPI_COMM_WORLD);
  if (rank==0) {
    fprintf(stderr, "%s: --------------------------\n", __FILE__ );
    fprintf(stderr, "%s: Finished sucessfully      \n", __FILE__ );
    fprintf(stderr, "%s: --------------------------\n", __FILE__ );
  }

  MPI_Finalize();
}
