#include <stdio.h>
#include <math.h>
#include "spline.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <inttypes.h>
#include "HPR_Cleaner.h"


void Proj(double *Ax0,double *x0,double *him,int *hidx,double *sval,int *istart,int ndata,int nspline)
{

  int i,j;
  double im[RGSIZE];
  memset(im,0,sizeof(double)*RGSIZE);
  memset(Ax0,0,sizeof(double)*nspline);

  double tot=0;
  for (i=0;i<nspline;i++) {
    tot+=x0[i];
  }
  for (i=0;i<nspline;i++) {
    Ax0[i]=tot;
  }
  
  for (i=0;i<ndata;i++) {
    if (hidx[i]>=0) {
      double tmp=0;
      for (j=0;j<4;j++) {
        tmp += sval[j+i*4] * x0[istart[i]+j];
      }
      im[hidx[i]]+=tmp;
    }
  }

  for (i=0;i<ndata;i++) {
    if (hidx[i]>=0) {
      double tmp=0;
      for (j=0;j<4;j++) {
        tmp += sval[j+i*4] * x0[istart[i]+j];
      }
    
      for (j=0;j<4;j++) {
        Ax0[istart[i]+j] += sval[j+i*4] * (him[hidx[i]]*tmp-im[hidx[i]]);
      }
    }
  }
  
}


int HPR_Cleaner( int *hidx,     // in: one ring of {pixname}_HPRIDX_ABER_TotalFlag_dx11 TOI
                 float *data,   // in: one ring of TOI to project to Splined HPR (SHPR)
                 float *ohpr,   // out: one produced SHPR ring of length RGSIZE, allocated by the caller
                 int ndata,     // number of samples in <hidx> and <data>
                 char *message) // prefix of message to display at the end of function, set to NULL to display nothing
{
  int i,j,k;
  BSpline *bspline  = NULL;
  double im[RGSIZE];
  double him[RGSIZE];
  double hpr[RGSIZE];
  int minidx=-1; // first index of valid data in <data>
  int maxidx=-1; // last index of valid data in <data>
  int nvalid=0;
  
  memset(im,0,sizeof(double)*RGSIZE);
  memset(him,0,sizeof(double)*RGSIZE);
  memset(ohpr,0,sizeof(float)*RGSIZE);
  memset(hpr,0,sizeof(double)*RGSIZE);

  for (i=0;i<ndata;i++) {
    if (hidx[i]>=0) {
      im[hidx[i]] += data[i];
      him[hidx[i]] += 1;
      if (minidx==-1) minidx = i;
      maxidx = i;
      nvalid++;
    }
  }

  // number of circles of data in the TOI
  int nspline=(int)((maxidx-minidx)/10800)+1;
  if (nspline<4) {
    memset(ohpr,0,sizeof(float)*RGSIZE);
    return(-1);
  }

  int rg_start,rg_end;
  double *sval = (double *) malloc(4*sizeof(double)*ndata);
  int *istart = (int *) malloc(sizeof(int)*ndata);
  double *vec = (double *) malloc(sizeof(double)*nspline);
  double *x0 = (double *) malloc(sizeof(double)*nspline);
  double *Ap0 = (double *) malloc(sizeof(double)*nspline);
  double *r0 = (double *) malloc(sizeof(double)*nspline);
  double *p0 = (double *) malloc(sizeof(double)*nspline);

  memset(x0,0,sizeof(double)*nspline);
  memset(p0,0,sizeof(double)*nspline);
  memset(Ap0,0,sizeof(double)*nspline);
  memset(r0,0,sizeof(double)*nspline);
  
  memset(vec,0,sizeof(double)*nspline);
  
  bspline = bspline_alloc (3,nspline-2);
  bspline_init_uniform (bspline, 0.0, nvalid-1.0);
  
  nvalid = 0;
  for (i=0;i<ndata;i++) {
    if (hidx[i]>=0) {
      bspline_value(bspline,(double)(nvalid),&rg_start,&rg_end);
      double *vals = bspline->vals;
      for (k=rg_start;k<=rg_end;k++) sval[k-rg_start+4*i] = vals[k];
      istart[i] = rg_start;
      nvalid++;
    }
  }
  
  for (i=0;i<ndata;i++) {
    if (hidx[i]>=0) {
      for (j=0;j<4;j++) {
        vec[istart[i]+j] += sval[j+i*4] * (him[hidx[i]]*data[i]-im[hidx[i]]);
      }
    }
  }

  Proj(r0,x0,him,hidx,sval,istart,ndata,nspline);
  double delta0 = 0;
  
  for (i=0;i<nspline;i++) {
    double tmp = vec[i]-r0[i];
    r0[i] = tmp;
    delta0 += tmp*tmp;
  }
  memcpy(p0,r0,nspline*sizeof(double));
  double deltaold = delta0;

  int itt=0;
  
  while(delta0>1E-24&&itt<nspline) {
    Proj(Ap0,p0,him,hidx,sval,istart,ndata,nspline);
    
    double denum_alpha=0;
    for (i=0;i<nspline;i++) {
      denum_alpha += p0[i]*Ap0[i];
    }
    double num_alpha=0;
    for (i=0;i<nspline;i++) {
      num_alpha += r0[i]*r0[i];
    }
    double alpha = num_alpha/denum_alpha;

    for (i=0;i<nspline;i++) x0[i] += alpha*p0[i];
    delta0=0;
    for (i=0;i<nspline;i++) {
      r0[i] -= alpha*Ap0[i];
    }
    double num_beta=0;
    for (i=0;i<nspline;i++) {
      num_beta += r0[i]*r0[i];
    }
    delta0 = num_beta;

    double beta = num_beta/num_alpha;
    
    for (i=0;i<nspline;i++) p0[i] = r0[i]+beta*p0[i];

    itt++;
  }
  if (message != NULL) {
    fprintf(stderr,"%s %d [%d,%d] %g %g\n",message,nspline,minidx,maxidx,deltaold,delta0);
  }

  for (i=0;i<ndata;i++) {
    if (hidx[i]>=0) {
      double tmp=0;
      for (j=0;j<4;j++) {
        tmp += sval[j+i*4] * x0[istart[i]+j];
      }
      
      hpr[hidx[i]] += data[i]-tmp;
    }
  }

  for (i=0;i<RGSIZE;i++) {
    if (him[i]>0) ohpr[i] = hpr[i]/him[i];
  }

  free(sval);
  free(istart);
  free(vec);
  free(x0);
  free(Ap0);
  free(r0);
  free(p0);
  bspline_free (bspline);

  return(0);
}
