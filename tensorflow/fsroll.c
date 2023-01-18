#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "spline.h"
#include "chealpix.h"
#include <math.h>

#ifdef DOMPI
#include "mpi.h"
#endif

#if 0
#define PDEBUG {fprintf(stderr,"%s %d %d\n",__FILE__,New->rank,__LINE__);}
#else
#define PDEBUG {}
#endif

#define FUNCTYPE_UNDEF 0
#define FUNCTYPE_SPLINE1D 1

void ang2pix_ring_table(long nside, double *theta, double *phi, long *ipix,long ndata)
{
  int i;
  for (i=0;i<ndata;i++) {
    ang2pix_ring(nside, *(theta+i), *(phi+i), ipix+i);
  }
}

void pix2ang_ring_table(long nside, long *ipix,double *theta, double *phi, long ndata)
{
  int i;
  for (i=0;i<ndata;i++) {
    pix2ang_ring(nside,*(ipix+i), theta+i, phi+i);
  }
}


typedef struct {
  double **vals;
  int **index;
} Function ;

typedef struct {
  double **AddData;
  double **Data;
  double **calib;
  double **Hdata;
  double **Hdata2;
  double **allval;
  double **tmp;
  double **ltmp;
  double ***tmp_Data;
  double *PixWeight;
  double *wfunc;
  double **correction;
  int *ncorrection;
  int ***icorrection;
  float ***iwcorrection;
  int **bidx;
  int **ridx;
  int *ncorr;
  long *ndata;
  long *nring;
  int *relativ;
  long nbolo;
  long ntmp;
  long nfunc;
  long npixel;
  long ntotring;
  long ntotfit;
  int rank;
  int size;
  int docalib;
  int doadd;
  int *maskinv;
  void **func_buffer;
  double ***func_data;
  int ***func_index;
  unsigned char *func_type;
  int *NFUNFC_IDX;
  int *NFUNFC_MIDX;
  int Nside;
  double A00;
} sroll_buffer;

long init_buffer(int Nside,int npixel,int nbolo,int ntmp,int nfunc,int ntotfunc,int ntotMfunc,int docalib,int doadd,int docorrection,int *incorr,int rank, int size)
{
  sroll_buffer *New = (sroll_buffer *) malloc(sizeof(sroll_buffer));
  long ptr=(long) New;
int i,j;
#if DOMPI
  int argc=1;
  char **argv;
  argv= (char **) malloc(sizeof(char *)*3);
  argv[0]=(char *) malloc(10);
  argv[1]=(char *) malloc(10);
  argv[2]=(char *) malloc(10);
  strcpy(argv[0],"SROll4");
  strcpy(argv[1],"");
  strcpy(argv[2],"");
  MPI_Init(&argc, &argv);
#endif
  New->Nside=Nside;
  New->rank=rank;
  New->size=size;
  New->nbolo=nbolo;
  New->ntmp=ntmp;
  New->nfunc=nfunc;
  New->wfunc=NULL;
  New->ncorr= (int *) malloc(sizeof(int)*(nbolo+1));
  memcpy(New->ncorr,incorr,(nbolo+1)*sizeof(int));
  if (nfunc>0) {
    New->func_buffer = (void **) malloc(sizeof(void *)*nfunc);
    New->func_data  = (double ***) malloc(sizeof(double ***)*npixel);
    New->func_index = (int ***) malloc(sizeof(int ***)*npixel);
    for (i=0;i<npixel;i++) {
      New->func_data[i]  = (double **) malloc(sizeof(double *)*ntotfunc);
      New->func_index[i] = (int **) malloc(sizeof(int *)*ntotfunc);
    }
    New->func_type = (unsigned char *) malloc(sizeof(unsigned char)*nfunc);
    for (i=0;i<nfunc;i++) New->func_type[i]=FUNCTYPE_UNDEF;
  }
  New->NFUNFC_IDX = (int *) malloc(sizeof(int)*(nfunc+1));
  New->NFUNFC_MIDX = (int *) malloc(sizeof(int)*(nfunc+1));
  for (i=0;i<nfunc+1;i++) {
    New->NFUNFC_MIDX[i]=0;
    New->NFUNFC_IDX[i]=0;
  }

  if (doadd==1) {
    New->AddData    = (double **) malloc(sizeof(double *)*npixel);
  }
  New->doadd      = doadd;
  New->Data       = (double **) malloc(sizeof(double *)*npixel);
  New->Hdata      = (double **) malloc(sizeof(double *)*npixel);
  New->Hdata2     = (double **) malloc(sizeof(double *)*npixel);
  New->tmp_Data   = (double ***) malloc(sizeof(double **)*npixel);
  New->allval     = (double **) malloc(sizeof(double *)*npixel);
  New->tmp        = (double **) malloc(sizeof(double *)*npixel);
  New->ltmp       = (double **) malloc(sizeof(double *)*npixel);
  if (docalib==1) New->calib = (double **) malloc(sizeof(double *)*npixel);

  if (docorrection>0) {
    New->icorrection = (int ***) malloc(sizeof(int **)*npixel);
    New->iwcorrection = (float ***) malloc(sizeof(float **)*npixel);
    New->ncorrection = (int *) malloc(sizeof(int)*docorrection);
    New->correction = (double **) malloc(sizeof(double *)*docorrection);
    for (i=0;i<npixel;i++) {
      New->icorrection[i] = (int **) malloc(sizeof(int *)*docorrection);
      New->iwcorrection[i] = (float **) malloc(sizeof(float *)*docorrection);
      for (j=0;j<docorrection;j++) {
	New->icorrection[i][j] = NULL;
	New->iwcorrection[i][j] = NULL;
      }
    }
  }
  for (i=0;i<npixel;i++) {
    New->tmp_Data[i] = (double **) malloc(sizeof(double *)*ntmp*nbolo);
    New->allval[i] = (double *) malloc(sizeof(double)*(ntmp+ntotMfunc)*nbolo);
  }
  New->bidx = (int **) malloc(sizeof(int *)*npixel);
  New->ridx = (int **) malloc(sizeof(int *)*npixel);
  New->relativ = (int *) malloc(sizeof(int)*ntmp);
  memset(New->relativ,0,sizeof(int)*ntmp);
  New->ndata = (long *) malloc(sizeof(long)*npixel);
  New->nring = (long *) malloc(sizeof(long)*nbolo);
  memset(New->ndata,0,sizeof(long)*npixel);
  memset(New->nring,0,sizeof(long)*nbolo);
  New->ntotring=0;
  New->ntmp=ntmp;
  New->nbolo=nbolo;
  New->npixel=npixel;
  New->docalib=docalib;
  New->PixWeight=(double *) malloc(sizeof(double)*npixel);
  New->maskinv=(int *) malloc(sizeof(int)*npixel);
  return(ptr);
}

void count_data_per_pixel(long buffer,int *ipidx,int ndata)
{
  sroll_buffer *New = (sroll_buffer *) buffer;
  int i,j,pidx;
  
  for (i=0;i<ndata;i++) {
    pidx=ipidx[i];
    if (pidx>=New->npixel) {
      fprintf(stderr,"Nb Pixel %d , found pix idx %d: PROBLEM\n",(int) pidx,(int) New->npixel);
      exit(0);
    }
    New->ndata[pidx]+=1;
  }
}

void InitSpline1D(long buffer,int ifunc,int nspline,double min,double max)
{

  sroll_buffer *New = (sroll_buffer *) buffer;
  int istart,iend,i,j;

  PDEBUG;
  BSpline *bspline = bspline_alloc (3,nspline-2);
  bspline_init_uniform (bspline, min, max);

  New->func_type[ifunc] = FUNCTYPE_SPLINE1D;
  New->func_buffer[ifunc] = (void *) bspline;
  New->NFUNFC_IDX[ifunc+1]=New->NFUNFC_IDX[ifunc]+4;
  New->NFUNFC_MIDX[ifunc+1]=New->NFUNFC_MIDX[ifunc]+nspline;
  PDEBUG;
}

void allocbolo(long buffer)
{
  sroll_buffer *New = (sroll_buffer *) buffer;
  int i,j,k,pidx;

  PDEBUG;  
  for (pidx=0;pidx<New->npixel;pidx++) {
    New->Data[pidx]            = (double *) malloc(sizeof(double)*New->ndata[pidx]);
    if (New->doadd==1) New->AddData[pidx] = (double *) malloc(sizeof(double)*New->ndata[pidx]);
    if (New->docalib==1) New->calib[pidx] = (double *) malloc(sizeof(double)*New->ndata[pidx]);
    New->Hdata[pidx]           = (double *) malloc(sizeof(double)*New->ndata[pidx]);
    New->Hdata2[pidx]          = (double *) malloc(sizeof(double)*New->ndata[pidx]);
    for (j=0;j<New->ntmp*New->nbolo;j++) {
      New->tmp_Data[pidx][j]   = (double *) malloc(sizeof(double)*New->ndata[pidx]);
    }
    for (j=0;j<New->ncorr[New->nbolo];j++) {
      New->icorrection[pidx][j] = (int *) malloc(sizeof(int)*New->ndata[pidx]);
      New->iwcorrection[pidx][j] = (float *) malloc(sizeof(float)*New->ndata[pidx]);
      for (k=0;k<New->ndata[pidx];k++) New->icorrection[pidx][j][k]= -1.0;
      for (k=0;k<New->ndata[pidx];k++) New->iwcorrection[pidx][j][k]= 1.0;
    }
    if (New->nfunc>0) {
      for (j=0;j<New->NFUNFC_IDX[New->nfunc];j++) {
	New->func_data[pidx][j]   = (double *) malloc(sizeof(double)*New->ndata[pidx]);
	New->func_index[pidx][j]  = (int *) malloc(sizeof(int)*New->ndata[pidx]);
      }
    }
    New->ridx[pidx]            = (int *) malloc(sizeof(int)*New->ndata[pidx]);
    New->bidx[pidx]            = (int *) malloc(sizeof(int)*New->ndata[pidx]);
  }
    
  memset(New->ndata,0,New->npixel*sizeof(long));
  PDEBUG;
}


void setbolo(long buffer,double *iadddata,double *idata,double *ihdata,
	     double *ihdata2,double *tmp_data,double *tmp_func,
	     double *icalib,int *ipidx,
	     int *iridx,int *icorrection,int ib,int ndata)
 {
  sroll_buffer *New = (sroll_buffer *) buffer;
  int i,j,k;
  int pidx;
  double *histo;
  PDEBUG;
  if (New->nfunc>0) {
    histo=(double *) malloc(sizeof(double)*New->NFUNFC_MIDX[New->nfunc]*New->nbolo);
    memset(histo,0,sizeof(double)*New->NFUNFC_MIDX[New->nfunc]*New->nbolo);
  }
  
  for (i=0;i<ndata;i++) {
    pidx=ipidx[i];
    New->AddData[pidx][New->ndata[pidx]]= iadddata[i];
    New->Data[pidx][New->ndata[pidx]]   = idata[i];
    if (New->docalib==1) New->calib[pidx][New->ndata[pidx]] = icalib[i];
    New->Hdata[pidx][New->ndata[pidx]]  = ihdata[i];
    New->Hdata2[pidx][New->ndata[pidx]] = ihdata2[i];
    
    for (j=0;j<New->ntmp;j++) {
      New->tmp_Data[pidx][j][New->ndata[pidx]] = tmp_data[i+ndata*j];
    }
    for (j=0;j<New->nfunc;j++) {
      if (New->func_type[j] == FUNCTYPE_SPLINE1D) {
	BSpline *bspline = (BSpline *) New->func_buffer[j];
	int istart,iend;
	bspline_value(bspline,tmp_func[i+ndata*j],&istart,&iend);
	double *vals = bspline->vals;
	if (istart==-1) {
	  fprintf(stderr,"Wrong domain %d %d %g\n",i,j,tmp_func[i+ndata*j]);
	  exit(0);
	}
	for (k=istart;k<=iend;k++) {
	  New->func_data[pidx][New->NFUNFC_IDX[j]+k-istart][New->ndata[pidx]]=vals[k];
	  New->func_index[pidx][New->NFUNFC_IDX[j]+k-istart][New->ndata[pidx]]=ib*New->NFUNFC_MIDX[New->nfunc]+New->NFUNFC_MIDX[j]+k;
	  histo[ib*New->NFUNFC_MIDX[New->nfunc]+New->NFUNFC_MIDX[j]+k]+=vals[k];
	}
      }
   }
    New->ridx[pidx][New->ndata[pidx]] = iridx[i];
    New->bidx[pidx][New->ndata[pidx]] = ib;
    for (j=New->ncorr[New->bidx[pidx][New->ndata[pidx]]];j<New->ncorr[New->bidx[pidx][New->ndata[pidx]]+1];j++) {
      New->icorrection[pidx][j][New->ndata[pidx]] = icorrection[i+ndata*(j-New->ncorr[New->bidx[pidx][New->ndata[pidx]]])];
      New->iwcorrection[pidx][j][New->ndata[pidx]] = 1.0;
    }
    New->ndata[pidx]+=1;
  }
#if 0
  if (New->nfunc>0) {
    fprintf(stderr,"HISTO ");
    for (i=0;i<New->NFUNFC_MIDX[New->nfunc]*New->nbolo;i++) fprintf(stderr,"%g ",histo[i]);
    fprintf(stderr,"\n");
    free(histo);
  }
#endif


  New->nring[ib]=0;
  for (i=0;i<ndata;i++) {
    if (New->nring[ib]<=iridx[i]) New->nring[ib]=iridx[i]+1;
  }
  PDEBUG;
}

void alloccorrection(long buffer,int ncorr,int icorr)
{
  int j;
  sroll_buffer *New = (sroll_buffer *) buffer;
  New->correction[icorr] = (double *) malloc(ncorr*sizeof(double));
  memset(New->correction[icorr],0,ncorr*sizeof(double));
  New->ncorrection[icorr] = ncorr;
}

void setcorrection(long buffer,double *correction,int icorr,int itemp)
{
  int j;
  sroll_buffer *New = (sroll_buffer *) buffer;
  
  if (itemp==-1) {
    memcpy(New->correction[icorr],correction,New->ncorrection[icorr]*sizeof(double));
  }
  else {
    int i,j,k;
    for (i=0;i<New->npixel;i++) {
      for (j=0;j<New->ndata[i];j++) {
	if (New->icorrection[i][icorr][j]!=-1) {
	  New->tmp_Data[i][itemp][j]=New->iwcorrection[i][icorr][j]*correction[New->icorrection[i][icorr][j]];
	}
      }
    }
  }
}

void setrelativ(long buffer,int *relativ)
{
  int j;
  sroll_buffer *New = (sroll_buffer *) buffer;
  memcpy(New->relativ,relativ,New->ntmp*sizeof(int));
  
  PDEBUG;
  if (New->rank==0) {
    fprintf(stderr,"\n=======================================\n");
    for (j=0;j<New->ntmp;j++) {
      if (New->relativ[j]==1) {
	fprintf(stderr,"Relativ %d\n" ,j);
      }
    }
    fprintf(stderr,"=======================================\n\n");
  }
  PDEBUG;
}

void nocalib(long buffer,int *icalib)
{
  sroll_buffer *New = (sroll_buffer *) buffer;
  int k;

  if (New->wfunc==NULL) New->wfunc= (double *) malloc(sizeof(double)*(New->ntmp+New->NFUNFC_MIDX[New->nfunc])*New->nbolo);
  
  for (k=0;k<New->ntmp*New->nbolo;k++) {
    New->wfunc[k]=1.0;
  }
  for (k=0;k<New->nbolo;k++) {
    New->wfunc[k*New->ntmp+icalib[k]]=0.0;
  }
  
  for (k=New->ntmp*New->nbolo;k<(New->ntmp+New->NFUNFC_MIDX[New->nfunc])*New->nbolo;k++) {
    New->wfunc[k]=1.0;
  }
}
void docalib(long buffer)
{
  sroll_buffer *New = (sroll_buffer *) buffer;
  int k;

  if (New->wfunc==NULL) New->wfunc= (double *) malloc(sizeof(double)*(New->ntmp+New->NFUNFC_MIDX[New->nfunc])*New->nbolo);

  for (k=0;k<New->ntmp*New->nbolo;k++) {
    New->wfunc[k]=1.0;
  }
  for (k=New->ntmp*New->nbolo;k<(New->ntmp+New->NFUNFC_MIDX[New->nfunc])*New->nbolo;k++) {
    New->wfunc[k]=0.0;
  }
  
}

void setbolo_gain(long buffer,double *gain,int cleancalib,int docorr)
{
  sroll_buffer *New = (sroll_buffer *) buffer;
  int i,j,k,l;

  for (i=0;i<New->npixel;i++) {
    for (j=0;j<New->ndata[i];j++) {
      double val=gain[New->bidx[i][j]]*New->Data[i][j];
      if (New->docalib==1&&cleancalib==1) val-=New->calib[i][j];
      if (docorr==1) {
	for (k=New->ncorr[New->bidx[i][j]];k<New->ncorr[New->bidx[i][j]+1];k++) 
	  {
	    if (New->icorrection[i][k][j]!=-1) {
	      val-=New->correction[k][New->icorrection[i][k][j]];
	    }
	  }
      }
      if (New->doadd==1) val+=New->AddData[i][j];
      New->tmp[i][j]=val;
    }
  }
  PDEBUG;
}

void initgrad(long buffer)
{
  sroll_buffer *New = (sroll_buffer *) buffer;
  int ib,i,j;

  PDEBUG;
  for (i=0;i<New->npixel;i++) {
    New->tmp[i]  = (double *) malloc(sizeof(double)*New->ndata[i]);
    New->ltmp[i]  = (double *) malloc(sizeof(double)*New->ndata[i]);
  }

  New->ntotring=0;
  for (i=0;i<New->nbolo;i++) {
    New->ntotring+=New->nring[i];
  }
  
  int *rg_pos = (int *) malloc(New->nbolo*sizeof(int));
  rg_pos[0]=0;
  for (i=1;i<New->nbolo;i++) {
    rg_pos[i]=rg_pos[i-1]+New->nring[i-1];
  }
  
  for (i=0;i<New->npixel;i++) {
    for (j=0;j<New->ndata[i];j++) {
      New->ridx[i][j]+=rg_pos[New->bidx[i][j]];
    }
  }
  free(rg_pos);

  New->ntotfit=New->ntotring+(New->ntmp+New->NFUNFC_MIDX[New->nfunc])*New->nbolo;
  if (New->rank==0) {
    fprintf(stderr,"nfit : %d\n",(int) New->ntotfit);
    fprintf(stderr,"ntot : %d\n",(int) New->ntotring);
    fprintf(stderr,"ntmp : %d\n",(int) New->ntmp);
    fprintf(stderr,"ntfunc : %d\n",(int) New->NFUNFC_MIDX[New->nfunc]);
  }
  for (i=0;i<New->npixel;i++) {
    double ohmap=0;
    for (j=0;j<New->ndata[i];j++) {
      ohmap += New->Hdata2[i][j];
    }
    if (ohmap>0) {
      New->PixWeight[i]=ohmap;
      for (j=0;j<New->ndata[i];j++) {
      	New->Hdata2[i][j]/=ohmap;
      }
      if (New->ndata[i]>1) New->maskinv[i]=1;
      else {
	New->PixWeight[i]=0;
	New->maskinv[i]=0;
      }
    }
    else {
      New->PixWeight[i]=0;
      New->maskinv[i]=0;
    }
  }
  PDEBUG;
}

void calc_A0(long buffer,double *A0)
{
  sroll_buffer *New = (sroll_buffer *) buffer;
  int i,j,k;

  PDEBUG;
  double A00=0;

  for (i=0;i<New->npixel;i++)  {
    if (New->maskinv[i]) {
      double omap=0;
      double ohmap=0;
     
      for (j=0;j<New->ndata[i];j++) {
	if (New->ridx[i][j]==0) {
	  New->ltmp[i][j]=1.0;
	  omap  += New->Hdata2[i][j];
	  ohmap += New->Hdata2[i][j];
	}
	else {
	  New->ltmp[i][j]=0.0;
	}
      }

      double vmap=0;
      for (j=0;j<New->ndata[i];j++) {
	vmap+=ohmap*New->ltmp[i][j]-omap;
      }
      for (j=0;j<New->ndata[i];j++) {
	if (New->ridx[i][j]==0) {
	  double destrip=ohmap*New->ltmp[i][j]-omap;
	  double val=ohmap*destrip;
	  A00+=New->PixWeight[i]*(val-New->Hdata2[i][j]*vmap);
	}
      }
    }
  }
  *A0=A00;
  New->A00=A00;
  PDEBUG;
}

void projdata(long buffer,double A0,double *y,int nocalib) 
{
  sroll_buffer *New = (sroll_buffer *) buffer;
  int ib,i,j,k;
  PDEBUG;
  memset(y,0,New->ntotfit*sizeof(double));
  
  for (i=0;i<New->npixel;i++) {
    if (New->maskinv[i]) {
      memset(New->allval[i],0,sizeof(double)*(New->ntmp+New->NFUNFC_MIDX[New->nfunc])*New->nbolo);
    }
  }
  
  for (i=0;i<New->npixel;i++) {
    if (New->maskinv[i])  {
      double omap=0;
      double ohmap=0;
      for (j=0;j<New->ndata[i];j++) {
	omap  += New->Hdata2[i][j]*New->tmp[i][j];
	ohmap += New->Hdata2[i][j];
	for (k=0;k<New->ntmp;k++) {
	  New->allval[i][k+New->bidx[i][j]*New->ntmp]+=New->tmp_Data[i][k][j]*New->Hdata2[i][j];
	}
	for (k=0;k<New->NFUNFC_IDX[New->nfunc];k++) {
	  New->allval[i][New->ntmp*New->nbolo+New->func_index[i][k][j]]+=New->func_data[i][k][j]*New->Hdata2[i][j];
	}
      }
    
      double vmap=0;
      for (j=0;j<New->ndata[i];j++) {
	double destrip=ohmap*New->tmp[i][j]-omap;
	vmap+=destrip;
      }
      for (k=0;k<(New->ntmp+New->NFUNFC_MIDX[New->nfunc])*New->nbolo;k++) {
	y[New->ntotring+k]-=New->wfunc[k]*New->PixWeight[i]*New->allval[i][k]*vmap;
      }

      for (j=0;j<New->ndata[i];j++) {
	double destrip=ohmap*New->tmp[i][j]-omap;
	double val=ohmap*destrip;
	y[New->ridx[i][j]]+=New->PixWeight[i]*(val-New->Hdata2[i][j]*vmap);
	for (k=0;k<New->ntmp;k++) {
	  y[New->ntotring+k+New->bidx[i][j]*New->ntmp]+=New->wfunc[k+New->bidx[i][j]*New->ntmp]*New->PixWeight[i]*(val*New->tmp_Data[i][k][j]);
	}
	for (k=0;k<New->NFUNFC_IDX[New->nfunc];k++) {
	  y[New->ntotring+New->nbolo*New->ntmp+New->func_index[i][k][j]]+=New->wfunc[New->nbolo*New->ntmp+New->func_index[i][k][j]]
	    *New->PixWeight[i]*(val*New->func_data[i][k][j]);
	}
      }
    }
  }
  PDEBUG;
}
   



void getndata(long buffer,int *ndata) 
{
 sroll_buffer *New = (sroll_buffer *) buffer;
 int ib,i,j,k,n=0;
  
 for (i=0;i<New->npixel;i++) {
   for (j=0;j<New->ndata[i];j++) {
     if (New->Hdata2[i][j]>0) n++;
   }
 }
 *ndata=n;
}

void getidxdata(long buffer,int *idxdata,float *wdata,int icorr,int icov) 
{
 sroll_buffer *New = (sroll_buffer *) buffer;
 int ib,i,j,k,ndata=0;
 

 PDEBUG;
 for (i=0;i<New->npixel;i++) {
   for (j=0;j<New->ndata[i];j++) {
     if (New->Hdata2[i][j]>0) {
       if (New->icorrection[i][icorr][j]!=-1) {
	 idxdata[ndata]=New->icorrection[i][icorr][j];
	 wdata[ndata]=New->iwcorrection[i][icorr][j];
       }
       else {
	 idxdata[ndata]=0;
	 wdata[ndata]=0.0;
       }
       ndata+=1;
     }
   }
 }

 
}

void getdataidx(long buffer,float *h2,int *pidx) 
{
 sroll_buffer *New = (sroll_buffer *) buffer;
 int i,j,k,ndata=0;
 
 
 for (i=0;i<New->npixel;i++) {
   for (j=0;j<New->ndata[i];j++) {
     if (New->Hdata2[i][j]>0) {
       if (New->ndata[i]>1) 
	 h2[ndata]=New->Hdata2[i][j];
       else
	 h2[ndata]=0;
       pidx[ndata]=i;
       ndata+=1;
     }
   }
 }
}

void getdata(long buffer,double *x,double *data,double *init_data,int icorr,int dotemplate,int tmpidx) 
{
 sroll_buffer *New = (sroll_buffer *) buffer;
 int ib,i,j,k,ndata=0;
 

 PDEBUG;
 memset(init_data,0,New->ncorrection[icorr]*sizeof(double));
 double *ninit_data = (double *) malloc(sizeof(double)*New->ncorrection[icorr]);
 memset(ninit_data,0,New->ncorrection[icorr]*sizeof(double));
   
 for (i=0;i<New->npixel;i++) {
   double avvmap=0;
   double navvmap=0;
   for (j=0;j<New->ndata[i];j++) {
     if (New->Hdata2[i][j]>0) {
#if 1
       double val=x[New->ridx[i][j]];
       for (k=0;k<New->ntmp;k++) if (k!=tmpidx) {
	 val+=New->tmp_Data[i][k][j]*x[New->ntotring+k+New->bidx[i][j]*New->ntmp];
       }
       for (k=0;k<New->NFUNFC_IDX[New->nfunc];k++) {
	 val+=x[New->ntotring+New->ntmp*New->nbolo+New->func_index[i][k][j]]*New->func_data[i][k][j];
       }
#endif
       for (k=New->ncorr[New->bidx[i][j]];k<New->ncorr[New->bidx[i][j]+1];k++)
	 if (New->icorrection[i][k][j]!=-1) {
	   val+=New->iwcorrection[i][k][j]*New->correction[k][New->icorrection[i][k][j]];
	 }
       avvmap+=New->Hdata2[i][j]*(New->tmp[i][j]-val);
       navvmap+=New->Hdata2[i][j];
     }
   }

   if (navvmap>0) avvmap/=navvmap;
   
   for (j=0;j<New->ndata[i];j++) {
     if (New->Hdata2[i][j]>0) {
#if 1
       double val=x[New->ridx[i][j]];
       for (k=0;k<New->ntmp;k++)  if (k!=tmpidx) {
	 val+=New->tmp_Data[i][k][j]*x[New->ntotring+k+New->bidx[i][j]*New->ntmp];
       }
       for (k=0;k<New->NFUNFC_IDX[New->nfunc];k++) {
	 val+=x[New->ntotring+New->ntmp*New->nbolo+New->func_index[i][k][j]]*New->func_data[i][k][j];
       }
#endif
       
       if (dotemplate==1) {
	   data[ndata]=New->Data[i][j];
	   if (New->doadd==1) data[ndata]+=New->AddData[i][j];
       }
       else {
	   data[ndata]=New->tmp[i][j] - val;
       }
       
       if (New->icorrection[i][icorr][j]!=-1) {
	 if (dotemplate==1) {
	   if (New->doadd==1) init_data[New->icorrection[i][icorr][j]]+=New->Hdata2[i][j]*(New->Data[i][j]+New->AddData[i][j]);
	   else init_data[New->icorrection[i][icorr][j]]+=New->Hdata2[i][j]*New->Data[i][j];
	 }
	 else {
	   init_data[New->icorrection[i][icorr][j]]+=New->Hdata2[i][j]*(data[ndata]);
	 }
	 ninit_data[New->icorrection[i][icorr][j]]+=New->Hdata2[i][j];
       }
       ndata+=1;
     }
   }
 }
 
 for (i=0;i<New->ncorrection[icorr];i++) {
   if (ninit_data[i]>0) init_data[i]/=ninit_data[i];
 }
 free(ninit_data);

 PDEBUG;
}
//  void domap(long buffer,double *x,int *validx,float *omap,float *ohmap,int *fbolo) 
  
void domap(long buffer,double *x,float *imap,float *hmap,int *fbolo,int *fring) 
{
 sroll_buffer *New = (sroll_buffer *) buffer;
 int ib,i,j,k;
 

 PDEBUG;
#if 0
 float *imap = (float *) malloc(sizeof(float)*New->npixel);
 float *hmap = (float *) malloc(sizeof(float)*New->npixel);
#endif

 for (i=0;i<New->npixel;i++) {
   imap[i]=0;
   hmap[i]=0;
   for (j=0;j<New->ndata[i];j++) {
     if (fbolo[New->bidx[i][j]]==1&&fring[New->ridx[i][j]]==1) {
       double val=x[New->ridx[i][j]];
       for (k=0;k<New->ntmp;k++) {
	 val+=New->tmp_Data[i][k][j]*x[New->ntotring+k+New->bidx[i][j]*New->ntmp];
       }
       for (k=0;k<New->NFUNFC_IDX[New->nfunc];k++) {
	 val+=x[New->ntotring+New->ntmp*New->nbolo+New->func_index[i][k][j]]*New->func_data[i][k][j];
       }
     
       imap[i] += New->Hdata[i][j]*(New->tmp[i][j]-val);
       hmap[i] += New->Hdata[i][j];
     }
   }
 }
 PDEBUG;
#if 0
 if (New->rank==0) {
   memset(omap,0,sizeof(float)*New->Nside*New->Nside*12);
   memset(ohmap,0,sizeof(float)*New->Nside*New->Nside*12);
   for (i=0;i<New->npixel;i++) omap[validx[i]]=imap[i];
   for (i=0;i<New->npixel;i++) ohmap[validx[i]]=hmap[i];
 }
#ifdef DOMPI
 if (New->rank==0) {
   int rk;
   for (rk=1;rk<New->size;rk++) {
     MPI_Status statu;
     fprintf(stderr,"%d\n",rk);
     long nvalue;
     MPI_Recv(&nvalue,sizeof(long), MPI_BYTE, rk,100, MPI_COMM_WORLD,&statu);
     int *lval_pidx = (int *) malloc(sizeof(int)*nvalue);
     float *lval = (float *) malloc(sizeof(float)*nvalue);
     MPI_Recv(lval_pidx,sizeof(int)*nvalue, MPI_BYTE, rk,101, MPI_COMM_WORLD,&statu);
     MPI_Recv(lval,sizeof(float)*nvalue, MPI_BYTE, rk,102, MPI_COMM_WORLD,&statu);
     for (i=0;i<nvalue;i++) omap[lval_pidx[i]]=lval[i];
     MPI_Recv(lval,sizeof(float)*nvalue, MPI_BYTE, rk,103, MPI_COMM_WORLD,&statu);
     for (i=0;i<nvalue;i++) ohmap[lval_pidx[i]]=lval[i];
     free(lval);
     free(lval_pidx);
   }
 }
 else {
   MPI_Send(&((long) (New->npixel)), sizeof(long), MPI_BYTE, 0, 100, MPI_COMM_WORLD);
   MPI_Send(validx, sizeof(int)*New->npixel, MPI_BYTE, 0, 101, MPI_COMM_WORLD);
   MPI_Send(imap, sizeof(float)*New->npixel, MPI_BYTE, 0, 102, MPI_COMM_WORLD);
   MPI_Send(hmap, sizeof(float)*New->npixel, MPI_BYTE, 0, 103, MPI_COMM_WORLD);
 }
#endif
 free(imap);
 free(hmap);
#endif
}


void proj(long buffer,double *x,double *y) 
{
 sroll_buffer *New = (sroll_buffer *) buffer;
 int i,j,k;

 PDEBUG;
 memset(y,0,New->ntotfit*sizeof(double));

 if (New->rank==0) {

   double normoffset=0;
   for (i=0;i<New->ntotring;i++) {
     normoffset+=x[i];
   }

   for (i=0;i<New->ntotring;i++) {
     y[i]=New->A00*normoffset;
   }

   for (j=0;j<New->ntmp;j++) {
     if (New->relativ[j]==1) {
       double val=0;

       for (i=0;i<New->nbolo;i++) {
	 val+=x[New->ntotring+i*New->ntmp+j];
       }
       for (i=0;i<New->nbolo;i++) {
	 y[New->ntotring+i*New->ntmp+j]=New->A00*val;
       }
     }
   }
   for (j=0;j<New->nfunc;j++) {
     for (k=0;k<New->nbolo;k++) {
       double val=0;

       for (i=k*New->NFUNFC_MIDX[New->nfunc]+New->NFUNFC_MIDX[j];i<k*New->NFUNFC_MIDX[New->nfunc]+New->NFUNFC_MIDX[j+1];i++) {
	 val+=x[New->ntotring+New->nbolo*New->ntmp+i];
       }

       for (i=k*New->NFUNFC_MIDX[New->nfunc]+New->NFUNFC_MIDX[j];i<k*New->NFUNFC_MIDX[New->nfunc]+New->NFUNFC_MIDX[j+1];i++) {
	 y[New->ntotring+New->nbolo*New->ntmp+i]+=New->A00*val;
       }

     }
   }
 }

 for (i=0;i<New->npixel;i++) {
   if (New->maskinv[i]) {
     double omap=0;
     double ohmap=0;
     
     for (j=0;j<New->ndata[i];j++) {
       double val=x[New->ridx[i][j]];
       for (k=0;k<New->ntmp;k++) {
	 val+=New->tmp_Data[i][k][j]*x[New->ntotring+k+New->bidx[i][j]*New->ntmp];
       }
       for (k=0;k<New->NFUNFC_IDX[New->nfunc];k++) {
	 val+=x[New->ntotring+New->ntmp*New->nbolo+New->func_index[i][k][j]]*New->func_data[i][k][j];
       }

       New->ltmp[i][j]=val;
       omap  += New->Hdata2[i][j]*val;
       ohmap += New->Hdata2[i][j];
     }
     
     double vmap=0;
     for (j=0;j<New->ndata[i];j++) {
       vmap+=ohmap*New->ltmp[i][j]-omap;
     }

     
     for (k=0;k<(New->ntmp+New->NFUNFC_MIDX[New->nfunc])*New->nbolo;k++) {
       y[New->ntotring+k]-=New->wfunc[k]*New->PixWeight[i]*New->allval[i][k]*vmap;
      }
     
     for (j=0;j<New->ndata[i];j++) {
       double destrip=ohmap*New->ltmp[i][j]-omap;
       double val=ohmap*destrip;
       y[New->ridx[i][j]]+=New->PixWeight[i]*(val-New->Hdata2[i][j]*vmap);
       for (k=0;k<New->ntmp;k++) {
	 y[New->ntotring+k+New->bidx[i][j]*New->ntmp]+=New->wfunc[k+New->bidx[i][j]*New->ntmp]*
	   New->PixWeight[i]*val*New->tmp_Data[i][k][j];
       }
       for (k=0;k<New->NFUNFC_IDX[New->nfunc];k++) {
	 y[New->ntotring+New->ntmp*New->nbolo+New->func_index[i][k][j]]+=New->wfunc[New->nbolo*New->ntmp+New->func_index[i][k][j]]
	    *New->PixWeight[i]*val*New->func_data[i][k][j];
       }
     }
   }
 }
 PDEBUG;
}
