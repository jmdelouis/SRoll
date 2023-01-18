
// this is for srand48(), drand48(), realpath()
#define _XOPEN_SOURCE 500

#define USEDII

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <mpi.h>
#include <float.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>   
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <fftw3.h>
#include <assert.h>
#include <errno.h>
#include <omp.h>

#include "no_dmc_metadata.h"
#include "no_dmc_data_access.h"
#include "no_dmc_piolib_type_def.h"
#include "no_dmc_util.h"
#include "no_dmc_debug.h"
#include "no_dmc_version.h"
#include "chealpix.h"

#include "sroll_param.h"

#include "stim_parLoader.h"
#include "stim.h"
#include "stim_tools.h"

#define ADURGSTEP (270)
#define ADUSTEP (3200)
#define NADUSTEP (ADUSTEP*ADURGSTEP)


enum MAPRINGS_values {
  FULL        = 1,
  HM12        = 2,
  S12345      = 4,
  YEAR12      = 8,
  FULLODDEVEN = 16,
  HM12ODDEVEN = 32
};

enum MAPRINGS_maps {
  MAPFULL  = 1,
  HM1      = 2,
  HM2      = 4,
  S1       = 8,
  S2       = 16,
  S3       = 32,
  S4       = 64,
  S5       = 128,
  YEAR1    = 256,
  YEAR2    = 512,
  FULLODD  = 1024,
  FULLEVEN = 2048,
  HM1ODD   = 4096,
  HM1EVEN  = 8192,
  HM2ODD   = 16384,
  HM2EVEN  = 32768
};


int getx12(int idx)
{
  int res=(idx%2)
    +((idx/4)%2)*2
    +((idx/16)%2)*4
    +((idx/64)%2)*8
    +((idx/256)%2)*16
    +((idx/1024)%2)*32
    +((idx/4096)%2)*64
    +((idx/16384)%2)*128
    +((idx/65536)%2)*256
    +((idx/262144)%2)*512
    +((idx/1048576)%2)*1024
    +((idx/4194304)%2)*2048;
  return(res);
}

int gety12(int idx)
{
  int res=((idx/2)%2)
    +((idx/8)%2)*2
    +((idx/32)%2)*4
    +((idx/128)%2)*8
    +((idx/512)%2)*16
    +((idx/2048)%2)*32
    +((idx/8192)%2)*64
    +((idx/32768)%2)*128
    +((idx/131072)%2)*256
    +((idx/524288)%2)*512
    +((idx/2097152)%2)*1024
    +((idx/8388608)%2)*2048;
  return(res);
}

int set12(int x,int y)
{
  int res=((x)%2)
    +((x/2)%2)*4
    +((x/4)%2)*16
    +((x/8)%2)*64
    +((x/16)%2)*256
    +((x/32)%2)*1024
    +((x/64)%2)*4096
    +((x/128)%2)*16384
    +((x/256)%2)*65536
    +((x/512)%2)*262144
    +((x/1024)%2)*1048576
    +((x/2048)%2)*4194304
    +((y)%2)*2
    +((y/2)%2)*8
    +((y/4)%2)*32
    +((y/8)%2)*128
    +((y/16)%2)*512
    +((y/32)%2)*2048
    +((y/64)%2)*8192
    +((y/128)%2)*32768
    +((y/256)%2)*131072
    +((y/512)%2)*524288
    +((y/1024)%2)*2097152
    +((y/2048)%2)*8388608;
  return(res);
}

ssize_t pwrite(int fildes, const void *buf, size_t nbyte,
              off_t offset);

int isPowerOfTwo (unsigned int x)
{
  return ((x != 0) && ((x & (~x + 1)) == x));
}

#define RINGSIZE (27664l)
void GetProcMem(long *vmem,long *phymem)
{
    FILE *fp;
    char name[256];
    int pid;
    char comm[512];
    char state;
    int ppid,pgrp,session,tty_nr,tpgid;
    unsigned long flags,minflt,cminflt,majflt,cmajflt,utime;
    long cutime,cstime,priority,nice,xxx,itrealvalue,starttime,vsize,rss,stime;

    sprintf(name,"/proc/%d/stat",(int) getpid());

    fp=fopen(name,"r");
    fscanf(fp,"%d %s %c %d %d %d %d %d %lu %lu %lu %lu %lu %lu %lu %ld %ld %ld %ld %ld %ld %lu %lu %ld",
           &pid,comm,&state,&ppid,&pgrp,&session,&tty_nr,&tpgid,&flags,&minflt,&cminflt,&majflt,
           &cmajflt,&utime,&stime,&cutime,&cstime,&priority,&nice,&xxx,&itrealvalue,&starttime,&vsize,
           &rss);
    fclose(fp);

    *vmem=vsize;
    *phymem=rss*4096;
}

float GetLoadAvg()
{
    char Line[256];
    float a,b,c;

    FILE *fp=fopen("/proc/loadavg","r");
    fgets(Line,256,fp);
    fclose(fp);
    sscanf(Line,"%f %f %f",&a,&b,&c);
    return(a);
}


#include "sroll_parLoader.h"

#ifdef NOGCITER
#define NOGCITER
#endif

#if 0
#define DORGG
#endif

#ifdef TRUECP
#define TRUECP
#endif

#if 0
#define TESPT fprintf(stderr,"***DBG*** %s:%d\n",__FILE__,__LINE__)
#else
#define TESPT {}
#endif

#define UNSEENPIX (-1.6375000E+30)

#ifdef FITTHETA
#define FITTHETA
#endif


#define DONSIDE 
#ifdef DONSIDE 
#endif
#define NSIDE (8)
#define ORIENT (4)

#define DOMAP
#ifdef DOMAP
#endif

//#define SOLDIPX (-0.00049639860)
//#define SOLDIPY (-0.0022099761)   
//#define SOLDIPZ (0.0024331235)

#define SOLDIPX (-0.00022735409)
#define SOLDIPY (-0.0022259792)   
#define SOLDIPZ (0.0025066439)
#define TINY 1.0e-303;


int ludcmp(double *a,double *d,int *indx,int n)
{
  int i,imax=0,j,k;
  double big,dum,sum,temp;
  double *vv = (double *) malloc(sizeof(double)*n);

  *d=1.0;
  for (i=0;i<n;i++) {
    big=0.0;
    for (j=0;j<n;j++)
      if ((temp=fabs(a[n*i+j])) > big) big=temp;
    if (big == 0.0) {
      fprintf(stderr,"Singular matrix in routine LUDCMP\n");
      return(-1);
    }
    vv[i]=1.0/big;
  }
  for (j=0;j<n;j++) {
    for (i=0;i<j;i++) {
      sum=a[n*i+j];
      for (k=0;k<i;k++) sum -= a[n*i+k]*a[n*k+j];
      a[n*i+j]=sum;
    }
    big=0.0;
    for (i=j;i<n;i++) {
      sum=a[n*i+j];
      for (k=0;k<j;k++)
        sum -= a[n*i+k]*a[n*k+j];
      a[n*i+j]=sum;
      if ( (dum=vv[i]*fabs(sum)) >= big) {
        big=dum;
        imax=i;
      }
    }
    if (j != imax) {
      for (k=0;k<n;k++) {
        dum=a[n*imax+k];
        a[n*imax+k]=a[n*j+k];
        a[n*j+k]=dum;
      }
      *d = -(*d);
      vv[imax]=vv[j];
    }
    indx[j]=imax;
    if (a[n*j+j] == 0.0) a[n*j+j]=TINY;
    if (j != n) {
      dum=1.0/(a[n*j+j]);
      for (i=j+1;i<n;i++) a[n*i+j] *= dum;
    }
  }
  free(vv);
  return(0);
}

#undef TINY

void lubksb(double *a,double *b,int *indx,int n)
{
  int i,ii=-1,ip,j;
  double sum;

  for (i=0;i<n;i++) {
    ip=indx[i];
    sum=b[ip];
    b[ip]=b[i];
    if (ii>-1)
      for (j=ii;j<=i-1;j++) sum -= a[n*i+j]*b[j];
    else if (sum!=0) ii=i;
    b[i]=sum;
  }

  for (i=n-1;i>=0;i--) {
    sum=b[i];
    for (j=i+1;j<n;j++) sum -= a[n*i+j]*b[j];
    b[i]=sum/a[n*i+i];
  }
}

int lusol(double *a,double *b,int n)
{
  double d;
  int *indx= (int *) malloc(n*sizeof(int));
  if (ludcmp(a,&d,indx,n)!=0) return(-1);
  lubksb(a,b,indx,n);
  free(indx);
  return(0);
}

void invert(double *a,double *y,int n)
{
  double d;
  double *col= (double *) malloc(n*sizeof(double));
  double *l_a= (double *) malloc(n*n*sizeof(double));
  memcpy(l_a,a,n*n*sizeof(double));
  int i,j;
  int *indx= (int *) malloc(n*sizeof(int));
  ludcmp(l_a,&d,indx,n);

  for (j=0;j<n;j++) {
    memset(col,0,n*sizeof(double));
    col[j]=1.0;
    lubksb(l_a,col,indx,n);
    for (i=0;i<n;i++) y[j+i*n]=col[i];
  }
      
  free(l_a);
  free(col);
  free(indx);
}

double norm(double *a,int n)
{
  int i,j;
  double maxd=0;
  double *col= (double *) malloc(n*sizeof(double));
  memset(col,0,n*sizeof(double));
  
  for (j=0;j<n;j++) {
    for (i=0;i<n;i++) col[j]+=fabs(a[i+j*n]);
    if (maxd<col[j]) maxd=col[j];
  }
  
  free(col);
  return(maxd);
}

double Dcond(double *a,int n)
{
  double *b=(double *) malloc(n*n*sizeof(double));
  
  invert(a,b,n);
  double cond=norm(a,n)*norm(b,n);
  free(b);
  return(cond);
}


int invert_3_3(double *mat,double *res)
{
  double a=mat[0];
  double b=mat[1];
  double c=mat[2];
  double d=mat[3];
  double e=mat[4];
  double f=mat[5];
  double g=mat[6];
  double h=mat[7];
  double i=mat[8];

  double determinant = a*(e*i-f*h)-d*(b*i-h*c)+g*(b*f-e*c);
  if (determinant==0) return(-1);
  double invdet = 1/determinant;

  res[0] =  (e*i-f*h)*invdet;
  res[3] = -(d*i-g*f)*invdet;
  res[6] =  (d*h-g*e)*invdet;
  res[1] = -(b*i-h*c)*invdet;
  res[4] =  (a*i-g*c)*invdet;
  res[7] = -(a*h-b*g)*invdet;
  res[2] =  (b*f-c*e)*invdet;
  res[5] = -(a*f-c*d)*invdet;
  res[8] =  (a*e-b*d)*invdet;
  return(0);
}


double donorm_3_3(double *mat)
{
  double a=mat[0];
  double b=mat[1];
  double c=mat[2];
  double d=mat[3];
  double e=mat[4];
  double f=mat[5];
  double g=mat[6];
  double h=mat[7];
  double i=mat[8];
  
  double res1= fabs(a)+fabs(b)+fabs(c);
  double res2= fabs(d)+fabs(e)+fabs(f);
  double res3= fabs(g)+fabs(h)+fabs(i);
  
  if (res2>res1) {
    if (res3>res2) return(res3);
    else return(res2);
  }
  if (res3>res1) return(res3);
  else return(res1);
}

double cond_3_3_thres(double a,double b,double c,
                      double d,double e,double f,
                      double g,double h,double i)
{
  double mat[9];
  double res[9];
  int ii,j,k;
  
  mat[0]=a;
  mat[1]=b;
  mat[2]=c;
  mat[3]=d;
  mat[4]=e;
  mat[5]=f;
  mat[6]=g;
  mat[7]=h;
  mat[8]=i;

  if (invert_3_3(mat,res)==-1) return(1E30);

  double cond=donorm_3_3(mat)*donorm_3_3(res);

  if (cond<1000) {
    double err=0;
    for (k=1;k<3;k++) {
      for (ii=0;ii<3;ii++) {
        for (j=0;j<3;j++) {
          err += res[3*ii+k]*res[3*k+j]*mat[3*ii+j];
        }
      }
    }

    cond=sqrt(err/2);
  }

  return(cond);
}

double legendre(double x,int n)
{
  int i;
  double P0=1;
  double P1=x;
  double P2=0;
  if (n==0) return(P0);
  if (n==1) return(P1);
  for (i=2;i<=n;i++) {
    P2=(2*(i-1)+1)*x*P1-(i-1)*P0;
    P0=P1;
    P1=P2;
  }
  return(P2);
}

double solvemap2(double a,double b,double c,double d,
                 double *v0,double *v1)
{

  double det=1./(a*d-b*c);
  double x0=det*(d*(*v0)-b*(*v1));
  double x1=det*(-c*(*v0)+a*(*v1));
  
  *v0=x0;
  *v1=x1;
  return(det);
}

double solvemap(double *x,double *y,double *z,
                double II, double IQ, double IU,
                double QQ, double QU, double UU)
{
  //double determinant =    +A(0,0)*(A(1,1)*A(2,2)-A(2,1)*A(1,2))
  //                     -A(0,1)*(A(1,0)*A(2,2)-A(1,2)*A(2,0))
  //                    +A(0,2)*(A(1,0)*A(2,1)-A(1,1)*A(2,0));
  double determinant = II*(QQ*UU-QU*QU)
    -IQ*(IQ*UU-QU*IU)
    +IU*(IQ*QU-QQ*IU);
  
  double invdet = 1/determinant;
 double result_0_0 =  (QQ*UU - QU*QU)*invdet; // (A(1,1)*A(2,2)-A(2,1)*A(1,2))*invdet;
 double result_1_0 = -(IQ*UU - IU*QU)*invdet; //-(A(0,1)*A(2,2)-A(0,2)*A(2,1))*invdet;
 double result_2_0 =  (IQ*QU - IU*QQ)*invdet; // (A(0,1)*A(1,2)-A(0,2)*A(1,1))*invdet;
 double result_0_1 = -(IQ*UU - QU*IU)*invdet; //-(A(1,0)*A(2,2)-A(1,2)*A(2,0))*invdet;
 double result_1_1 =  (II*UU - IU*IU)*invdet; // (A(0,0)*A(2,2)-A(0,2)*A(2,0))*invdet;
 double result_2_1 = -(II*QU - IQ*IU)*invdet; //-(A(0,0)*A(1,2)-A(1,0)*A(0,2))*invdet;
 double result_0_2 =  (IQ*QU - IU*QQ)*invdet; // (A(1,0)*A(2,1)-A(2,0)*A(1,1))*invdet;
 double result_1_2 = -(II*QU - IU*IQ)*invdet; //-(A(0,0)*A(2,1)-A(2,0)*A(0,1))*invdet;
 double result_2_2 =  (II*QQ - IQ*IQ)*invdet; // (A(0,0)*A(1,1)-A(1,0)*A(0,1))*invdet;

 double ox = *x*result_0_0 + *y*result_1_0 + *z*result_2_0;
 double oy = *x*result_0_1 + *y*result_1_1 + *z*result_2_1;
 double oz = *x*result_0_2 + *y*result_1_2 + *z*result_2_2;

 *x=ox;
 *y=oy;
 *z=oz;

 return(determinant);
}

// maximum number of stim iterations in one sroll run
#ifndef MAXSIMU
#define MAXSIMU (1)
#endif

// * maximum number of Theo_noPS HPR templates per bolometer
// * from main():
//     npixbeam = Param->n_Theo_noPS / nbolo;
//     assert( npixbeam <= MAXTHEOHPR);
// * can be set at compilation time with '-DMAXTHEOHPR=xxx'
// * values used in RD12 are 4 at 100ghz-217ghz, 7 at 353ghz and 9 at 545ghz-857ghz
// * need to remember why it was set to 11 for RD12...
#ifndef MAXTHEOHPR
#define MAXTHEOHPR (7)
#endif


typedef struct {
  PIODOUBLE sig;
  PIODOUBLE listp[MAXSIMU];
  PIODOUBLE listofpix[MAXTHEOHPR];
  PIOFLOAT comap;
  PIOFLOAT dustmap;
  PIOFLOAT freefree;
  PIOINT ipix;
  PIOINT rg;
  PIOINT hrg;
  PIOINT adu;
  PIOINT gi;
  //PIOFLOAT ggi;
  //PIOFLOAT xrg;
  PIODOUBLE dip;
  //PIOFLOAT soldip;
  //PIODOUBLE thsig;
  PIOFLOAT fsl;
  PIOFLOAT w;
  PIOFLOAT co;
  PIOFLOAT si; // pixel pour 32 orientation et 12*32*32 pixel
  PIOBYTE  surv;
  PIOBYTE  ib;
  PIODOUBLE hit;
  //PIODOUBLE sig1;
  //PIODOUBLE sri;
  PIODOUBLE wp;
  PIODOUBLE vi,vq,vu;
  PIODOUBLE model;
} hpix;

#define _PIOMALLOC malloc
#define _PIOFREE free

double gcmat_mpi(double *mat,double *vec,int n,int nn,long begr,long edr)
{
  int itermax = 1000;
  double eps = 1.e-14;
  int iter;
  int i,j;
  int rank,mpi_size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&mpi_size);
  int *tab_begr =_PIOMALLOC(sizeof(int)*mpi_size);
  int *tab_edr  =_PIOMALLOC(sizeof(int)*mpi_size);

  tab_begr[rank]=begr;
  tab_edr[rank]=edr;

  for (i=0;i<mpi_size;i++) {
    MPI_Bcast(tab_begr+i, sizeof(int), MPI_BYTE, i, MPI_COMM_WORLD);
    MPI_Bcast(tab_edr+i , sizeof(int), MPI_BYTE, i, MPI_COMM_WORLD);
  }

  double delta0, delta_new, alpha, delta_old, beta;
  

  double *x = (double *) calloc (n, sizeof (double));
  double *b = (double *) calloc (n, sizeof (double));
  double *d = (double *) calloc (n, sizeof (double));
  double *q = (double *) calloc (n, sizeof (double));
  double *r = (double *) calloc (n, sizeof (double));
  double *s = (double *) calloc (n, sizeof (double));
  double *hit = (double *) calloc (n, sizeof (double));
  
  iter = 0;
  // Compute second member
  for (i=0;i<mpi_size;i++) {
    memcpy(b+tab_begr[i],vec+tab_begr[i],(tab_edr[i]-tab_begr[i]+1)*sizeof(double));
    MPI_Bcast(b+tab_begr[i], sizeof(double)*(tab_edr[i]-tab_begr[i]+1), MPI_BYTE, i, MPI_COMM_WORLD);
  }
  // starting point x = solution
  for (i=0; i<n; i++) x[i] = 0.; 

  // Compute Ax
  //healspline_fill_s_with_PtP_s (hs_rec, Ndata, vec, x, q);
  for (i=begr;i<= edr;i++) {
    double tmp=0;
    for (j=0;j<n;j++) tmp=tmp+x[j]*mat[j+(i-begr)*nn];
    q[i]=tmp;
  }

  for (i=begr; i <= edr; i++)
    {
      r[i] = b[i] - q[i];
      d[i] = r[i];
    }

  // Preconditionnement
  for (i=begr; i <=edr; i++) {
    hit[i]=mat[i+(i-begr)*nn];
    if (hit[i]==0) hit[i]=1;
  }
  for (i=begr; i<= edr; i++)
    s[i] = r[i] / hit[i];
  
  double delta_new_tmp = 0.0;
  for (i = begr; i <= edr; i++) delta_new_tmp += r[i] * s[i];
  delta_new=0;
  for (i=0;i<mpi_size;i++) {
    double tmp=delta_new_tmp;
    MPI_Bcast(&tmp, sizeof(double), MPI_BYTE, i, MPI_COMM_WORLD);
    delta_new+=tmp;
  }
  delta0 = delta_new;
  

  if (rank==0) fprintf (stderr, "iter = %d - delta0 = %lg - delta_new = %lg\n", iter, delta0, delta_new);

  while (iter < itermax && delta_new > eps*eps*delta0)
    {
      // q <= Ad
      for (i=begr; i <= edr; i++)
        q[i] = 0.0;
   
      //healspline_fill_s_with_PtP_s (hs_rec, Ndata, vec, d, q);
      for (i=0;i<mpi_size;i++) {
        MPI_Bcast(d+tab_begr[i], sizeof(double)*(tab_edr[i]-tab_begr[i]+1), MPI_BYTE, i, MPI_COMM_WORLD);
      }
      for (i=begr;i<= edr;i++) {
        double tmp=0;
        for (j=0;j<n;j++) tmp=tmp+d[j]*mat[j+(i-begr)*nn];
        q[i]=tmp;
      }
      
      //
      double dtq = 0.0;
      double dtq_tmp = 0.0;
      for (i=begr; i <= edr; i++)
        dtq_tmp += d[i] * q[i];
      for (i=0;i<mpi_size;i++) {
        double tmp=dtq_tmp;
        MPI_Bcast(&tmp, sizeof(double), MPI_BYTE, i, MPI_COMM_WORLD);
        dtq+=tmp;
      }
      alpha = delta_new / dtq;
      for (i=begr; i<= edr; i++)
        x[i] += alpha * d[i];

      if (iter % 100 == 0)
        {
          if (rank==0) fprintf (stderr,"gcmat_mpi() iter = %d - delta0 = %lg - delta_new = %lg\n", iter, delta0, delta_new);
          for (i=0; i < n; i++)
            q[i] = 0.0;

          for (i=0;i<mpi_size;i++) {
            MPI_Bcast(x+tab_begr[i], sizeof(double)*(tab_edr[i]-tab_begr[i]+1), MPI_BYTE, i, MPI_COMM_WORLD);
          }

          //healspline_fill_s_with_PtP_s (hs_rec, Ndata, vec, x, q);
          for (i=begr;i<= edr;i++) {
            double tmp=0;
            for (j=0;j<n;j++) tmp=tmp+x[j]*mat[j+(i-begr)*nn];
            q[i]=tmp;
          }
          for (i=begr; i <= edr; i++) {
            r[i] = b[i] - q[i];
          }
        }
      else
        {
          for (i=begr; i<= edr; i++)
            r[i] -= alpha * q[i];
        }
      for (i=begr; i<= edr; i++)
        s[i] = r[i] / hit[i];

      delta_old = delta_new;
      delta_new = 0.0;

      delta_new_tmp = 0.0;
      for (i = begr; i <= edr; i++) delta_new_tmp += r[i] * s[i];
      delta_new=0;
      for (i=0;i<mpi_size;i++) {
        double tmp=delta_new_tmp;
        MPI_Bcast(&tmp, sizeof(double), MPI_BYTE, i, MPI_COMM_WORLD);
        delta_new+=tmp;
      }
      beta = delta_new / delta_old;

      for (i=begr; i<= edr; i++)
        d[i] = s[i] + beta * d[i];
      iter ++;
    }
  if (rank==0) fprintf (stderr,"gcmat_mpi2() iter = %d - delta0 = %lg - delta_new = %lg\n", 
          iter, delta0, delta_new);
  if (rank==0) fprintf (stderr,"CG in iter = %d (max=%d)\n", iter, itermax);
  for (i=0;i<mpi_size;i++) {
    memcpy(vec+tab_begr[i],x+tab_begr[i],(tab_edr[i]-tab_begr[i]+1)*sizeof(double));
    MPI_Bcast(vec+tab_begr[i], sizeof(double)*(tab_edr[i]-tab_begr[i]+1), MPI_BYTE, i, MPI_COMM_WORLD);
  } 

  _PIOFREE(tab_begr);
  _PIOFREE(tab_edr);

  return(delta_new/delta0);
}


double gcmat(double *mat,double *vec,int n,int nn)
{
  int itermax = 1000;
  double eps = 1.e-6;
  int iter;
  int i,j;
  
  long begr=0;
  long edr =n-1;
  double delta0, delta_new, alpha, delta_old, beta;
  double *x = (double *) calloc (n, sizeof (double));
  double *b = (double *) calloc (n, sizeof (double));
  double *d = (double *) calloc (n, sizeof (double));
  double *q = (double *) calloc (n, sizeof (double));
  double *r = (double *) calloc (n, sizeof (double));
  double *s = (double *) calloc (n, sizeof (double));
  double *hit = (double *) calloc (n, sizeof (double));
  
  memcpy(b,vec,n*sizeof(double));

  iter = 0;
  // starting point x = solution
  for (i=0; i<n; i++) x[i] = 0.; 

  // Compute Ax
  //healspline_fill_s_with_PtP_s (hs_rec, Ndata, vec, x, q);
  for (i=begr;i<= edr;i++) {
    double tmp=0;
    for (j=0;j<n;j++) tmp=tmp+x[j]*mat[j+(i-begr)*nn];
    q[i]=tmp;
  }

  for (i=begr; i <= edr; i++)
    {
      r[i] = b[i] - q[i];
      d[i] = r[i];
    }

  // Preconditionnement
  for (i=begr; i <=edr; i++) {
    hit[i]=mat[i+(i-begr)*nn];
    if (hit[i]==0) hit[i]=1;
  }
  for (i=begr; i<= edr; i++) {
    s[i] = r[i] / hit[i];
    d[i] /= hit[i];
  }
  
  delta_new = 0.0;
  for (i = begr; i <= edr; i++) delta_new += r[i] * s[i];
  delta0 = delta_new;

  fprintf (stderr, "gcmat() iter = %d - delta0 = %lg - delta_new = %lg\n", iter, delta0, delta_new);

  while (iter < itermax && delta_new > eps*eps*delta0)
    {

      // q <= Ad
      for (i=begr; i <= edr; i++)
        q[i] = 0.0;
   
      //healspline_fill_s_with_PtP_s (hs_rec, Ndata, vec, d, q);
      for (i=begr;i<= edr;i++) {
        double tmp=0;
        for (j=0;j<n;j++) tmp=tmp+d[j]*mat[j+(i-begr)*nn];
        q[i]=tmp;
      }
      
      //
      double dtq = 0.0;
      for (i=begr; i <= edr; i++) {
        dtq += d[i] * q[i];
      }

      alpha = delta_new / dtq;
      for (i=begr; i<= edr; i++)
        x[i] += alpha * d[i];

      fprintf (stderr,"gcmat2() iter = %d - delta0 = %lg - delta_new = %lg\n", iter, delta0, delta_new);
      if (iter % 100 == 0 && iter != 0)
        {
          for (i=0; i < n; i++)
            q[i] = 0.0;

          //healspline_fill_s_with_PtP_s (hs_rec, Ndata, vec, x, q);
          for (i=begr;i<= edr;i++) {
            double tmp=0;
            for (j=0;j<n;j++) tmp=tmp+x[j]*mat[j+(i-begr)*nn];
            q[i]=tmp;
          }
          for (i=begr; i <= edr; i++) {
            r[i] = b[i] - q[i];
          }
        }
      else
        {
          for (i=begr; i<= edr; i++)
            r[i] -= alpha * q[i];
        }
      for (i=begr; i<= edr; i++)
        s[i] = r[i] / hit[i];

      delta_old = delta_new;
      delta_new = 0.0;
      for (i = begr; i <= edr; i++) delta_new += r[i] * s[i];
      beta = delta_new / delta_old;

      for (i=begr; i<= edr; i++)
        d[i] = s[i] + beta * d[i];
      iter ++;
    }
  fprintf (stderr,"gcmat3() iter = %d - delta0 = %lg - delta_new = %lg\n", 
          iter, delta0, delta_new);
  fprintf (stderr,"CG in iter = %d (max=%d)\n", iter, itermax);

  memcpy(vec,x,n*sizeof(double));
  return(delta_new/delta0);
}
           
void matmul(double *in1,double *in2,double *out,int n)
{
  int i,j,k;

  for (i=0;i<n;i++) {
    for (j=0;j<n;j++) {
      double tmp=0;
      for (k=0;k<n;k++) {
        tmp+=in1[i+n*k]*in2[k+n*j];
      }
      out[i+j*n]=tmp;
    }
  }
}

void euler_matrix_new(double a1,double a2,double a3,int itype,double *mat)
{
  double c1 = cos(a1);
  double s1 = sin(a1);
  double c2 = cos(a2);
  double s2 = sin(a2);
  double c3 = cos(a3);
  double s3 = sin(a3);
  
  double ze = 0.0;
  double un = 1.0;

  double m1[9], m2[9], m3[9];

  //-- form full rotation matrix --
  //-- FOR REFERENCE:
  //
  // itype : 1
  // m1 = [[ c1,-s1,  0],[ s1, c1,  0],[  0,  0,  1]] ; around   z
  // m2 = [[ c2,  0, s2],[  0,  1,  0],[-s2,  0, c2]] ; around   y
  // m3 = [[  1,  0,  0],[  0, c3,-s3],[  0, s3, c3]] ; around   x
  //
  // itype : 2
  // m1 = [[ c1,-s1,  0],[ s1, c1,  0],[  0,  0,  1]] ; around   z
  // m2 = [[  1,  0,  0],[  0, c2,-s2],[  0, s2, c2]] ; around   x
  // m3 = [[ c3,-s3,  0],[ s3, c3,  0],[  0,  0,  1]] ; around   z
  //
  // itype : 3
  // m1 = [[ c1,-s1,  0],[ s1, c1,  0],[  0,  0,  1]] ; around   z
  // m2 = [[ c2,  0, s2],[  0,  1,  0],[-s2,  0, c2]] ; around   y
  // m3 = [[ c3,-s3,  0],[ s3, c3,  0],[  0,  0,  1]] ; around   z
  //
  //
  // matrix = m1 ( m2 m3 )

  //m1(:,1) = (/ c1,-s1, ze /)
  //m1(:,2) = (/ s1, c1, ze /)
  // m1(:,3) = (/ ze, ze, un /)
  m1[0]=c1;  m1[3]=-s1;  m1[6]=ze;
  m1[1]=s1;  m1[4]= c1;  m1[7]=ze;
  m1[2]=ze;  m1[5]= ze;  m1[8]=un;

  switch(itype) {

  case 1: // ZYX

    //m2(:,1) = (/ c2, ze, s2 /)
    //m2(:,2) = (/ ze, un, ze /)
    //m2(:,3) = (/-s2, ze, c2 /)
    m2[0]=c2;  m2[3]=ze;  m2[6]=s2;
    m2[1]=ze;  m2[4]=un;  m2[7]=ze;
    m2[2]=-s2; m2[5]=ze;  m2[8]=c2;

    //m3(:,1) = (/ un, ze, ze /)
    //m3(:,2) = (/ ze, c3,-s3 /)
    //m3(:,3) = (/ ze, s3, c3 /)
    m3[0]=un;  m3[3]=ze;  m3[6]=ze;
    m3[1]=ze;  m3[4]=c3;  m3[7]=-s3;
    m3[2]=ze;  m3[5]=s3;  m3[8]=c3;
    break;

  case 2: // ZXZ

    //m2(:,1) = (/ un, ze, ze /)
    //m2(:,2) = (/ ze, c2,-s2 /)
    //m2(:,3) = (/ ze, s2, c2 /)
    m2[0]=un;  m2[1]=ze;  m2[2]=ze;
    m2[3]=ze;  m2[4]=c2;  m2[5]=-s2;
    m2[6]=ze;  m2[7]=s2;  m2[8]=c2;


    //m3(:,1) = (/ c3,-s3, ze /)
    //m3(:,2) = (/ s3, c3, ze /)
    //m3(:,3) = (/ ze, ze, un /)
    m3[0]=c3;  m3[1]=-s3; m3[2]=ze;
    m3[3]=s3;  m3[4]=c3;  m3[5]=ze;
    m3[6]=ze;  m3[7]=ze;  m3[8]=un;
    break;

  case 3: // ZYZ

    //m2(:,1) = (/ c2, ze, s2 /)
    //m2(:,2) = (/ ze, un, ze /)
    //m2(:,3) = (/-s2, ze, c2 /)
    m2[0]=c2;  m2[1]=ze;  m2[2]=s2;
    m2[3]=ze;  m2[4]=un;  m2[5]=ze;
    m2[6]=-s2; m2[7]=ze;  m2[8]=c2;

    //m3(:,1) = (/ c3,-s3, ze /)
    //m3(:,2) = (/ s3, c3, ze /)
    //m3(:,3) = (/ ze, ze, un /)
    m3[0]=c3;  m3[1]=-s3; m3[2]=ze;
    m3[3]=s3;  m3[4]=c3;  m3[5]=ze;
    m3[6]=ze;  m3[7]=ze;  m3[8]=un;
    break;

  default:
    assert( 0 && "sroll.c:euler_matrix_new() error, invalid itype value.");

  }
  double res[9];
  //matrix = matmul(m1,matmul(m2,m3))
  matmul(m2,m3,res,3);
  matmul(m1,res,mat,3);
}

//==================================================================================
//
//    DFPMIN MINIMISATION
//
//==================================================================================
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
int ncom;       /* defined in DLINMIN */
double *pcom=0,*xicom=0,(*nrfunc)(double *);
void (*nrdfun)(double *,double *);


#define ITMAX 10000
#define EPS 1.0e-16
#define FREEALL free_vector(xi,1,n);free_vector(h,1,n);free_vector(g,1,n);
#define TOL 2.0e-4

int ncom=0;     /* defining declarations */

#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define NRMAX(a,b) ((a) > (b) ? (a) : (b))
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

void mnbrak(double *ax,double *bx,double *cx,double *fa,double *fb,double *fc,double (*func)(double))
{
  double ulim,u,r,q,fu,dum;

  *fa=(*func)(*ax);
  *fb=(*func)(*bx);
  if (*fb > *fa) {
    SHFT(dum,*ax,*bx,dum)
    SHFT(dum,*fb,*fa,dum)
  }
  *cx=(*bx)+GOLD*(*bx-*ax);
  *fc=(*func)(*cx);
  while (*fb > *fc) {
    r=(*bx-*ax)*(*fb-*fc);
    q=(*bx-*cx)*(*fb-*fa);
    u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
      (2.0*SIGN(NRMAX(fabs(q-r),TINY),q-r));
    ulim=(*bx)+GLIMIT*(*cx-*bx);
    if ((*bx-u)*(u-*cx) > 0.0) {
      fu=(*func)(u);
      if (fu < *fc) {
        *ax=(*bx);
        *bx=u;
        *fa=(*fb);
        *fb=fu;
        return;
      } else if (fu > *fb) {
        *cx=u;
        *fc=fu;
        return;
      }
      u=(*cx)+GOLD*(*cx-*bx);
      fu=(*func)(u);
    } else if ((*cx-u)*(u-ulim) > 0.0) {
      fu=(*func)(u);
      if (fu < *fc) {
        SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
        SHFT(*fb,*fc,fu,(*func)(u))
      }
    } else if ((u-ulim)*(ulim-*cx) >= 0.0) {
      u=ulim;
      fu=(*func)(u);
    } else {
      u=(*cx)+GOLD*(*cx-*bx);
      fu=(*func)(u);
    }
    SHFT(*ax,*bx,*cx,u)
    SHFT(*fa,*fb,*fc,fu)
  }
}

#undef GOLD
#undef GLIMIT
#undef TINY
#undef NRMAX
#undef SHFT

#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

double brent(double ax,double bx,double cx,double (*f)(double ),double tol,double *xmin)
{
  int iter;
  double a,b,d=0.0,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  double e=0.0;

  a=((ax < cx) ? ax : cx);
  b=((ax > cx) ? ax : cx);
  x=w=v=bx;
  fw=fv=fx=(*f)(x);
  for (iter=0;iter<ITMAX;iter++) {
    xm=0.5*(a+b);
    tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
    if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
      *xmin=x;
      return fx;
    }
    if (fabs(e) > tol1) {
      r=(x-w)*(fx-fv);
      q=(x-v)*(fx-fw);
      p=(x-v)*q-(x-w)*r;
      q=2.0*(q-r);
      if (q > 0.0) p = -p;
      q=fabs(q);
      etemp=e;
      e=d;
      if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
        d=CGOLD*(e=(x >= xm ? a-x : b-x));
      else {
        d=p/q;
        u=x+d;
        if (u-a < tol2 || b-u < tol2)
          d=SIGN(tol1,xm-x);
      }
    } else {
      d=CGOLD*(e=(x >= xm ? a-x : b-x));
    }
    u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
    fu=(*f)(u);
    if (fu <= fx) {
      if (u >= x) a=x; else b=x;
      SHFT(v,w,x,u)
      SHFT(fv,fw,fx,fu)
    } else {
      if (u < x) a=u; else b=u;
      if (fu <= fw || w == x) {
        v=w;
        w=u;
        fv=fw;
        fw=fu;
      } else if (fu <= fv || v == x || v == w) {
        v=u;
        fv=fu;
      }
    }
  }
  fprintf(stderr,"Too many iterations in BRENT");
  *xmin=x;
  return fx;
}

#undef CGOLD
#undef ZEPS

#define NRANSI
#define ZEPS 1.0e-10
#define MOV3(a,b,c, d,e,f) (a)=(d);(b)=(e);(c)=(f);

double dbrent(double ax, double bx, double cx, double (*f)(double),
        double (*df)(double), double tol, double *xmin)
{
  int iter,ok1,ok2;
  double a,b,d=0.0,d1,d2,du,dv,dw,dx,e=0.0;
  double fu,fv,fw,fx,olde,tol1,tol2,u,u1,u2,v,w,x,xm;
  
  a=(ax < cx ? ax : cx);
  b=(ax > cx ? ax : cx);
  x=w=v=bx;
  fw=fv=fx=(*f)(x);
  dw=dv=dx=(*df)(x);
  for (iter=0;iter<ITMAX;iter++) {
    xm=0.5*(a+b);
    tol1=tol*fabs(x)+ZEPS;
    tol2=2.0*tol1;
    if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
      *xmin=x;
      return fx;
    }
    if (fabs(e) > tol1) {
      d1=2.0*(b-a);
      d2=d1;
      if (dw != dx) d1=(w-x)*dx/(dx-dw);
      if (dv != dx) d2=(v-x)*dx/(dx-dv);
      u1=x+d1;
      u2=x+d2;
      ok1 = (a-u1)*(u1-b) > 0.0 && dx*d1 <= 0.0;
      ok2 = (a-u2)*(u2-b) > 0.0 && dx*d2 <= 0.0;
      olde=e;
      e=d;
      if (ok1 || ok2) {
        if (ok1 && ok2)
          d=(fabs(d1) < fabs(d2) ? d1 : d2);
        else if (ok1)
          d=d1;
        else
          d=d2;
        if (fabs(d) <= fabs(0.5*olde)) {
          u=x+d;
          if (u-a < tol2 || b-u < tol2)
            d=SIGN(tol1,xm-x);
        } else {
          d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
        }
      } else {
        d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
      }
    } else {
      d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
    }
    if (fabs(d) >= tol1) {
      u=x+d;
      fu=(*f)(u);
    } else {
      u=x+SIGN(tol1,d);
      fu=(*f)(u);
      if (fu > fx) {
        *xmin=x;
        return fx;
      }
    }
    du=(*df)(u);
    if (fu <= fx) {
      if (u >= x) a=x; else b=x;
      MOV3(v,fv,dv, w,fw,dw)
      MOV3(w,fw,dw, x,fx,dx)
      MOV3(x,fx,dx, u,fu,du)
    } else {
      if (u < x) a=u; else b=u;
      if (fu <= fw || w == x) {
        MOV3(v,fv,dv, w,fw,dw)
        MOV3(w,fw,dw, u,fu,du)
      } else if (fu < fv || v == x || v == w) {
        MOV3(v,fv,dv, u,fu,du)
          }
    }
  }
  fprintf(stderr,"Too many iterations in routine dbrent");
  return 0.0;
}

double f1dim(double x)
{
  int j;
  double f,*xt;

  xt= (double *) malloc(ncom*sizeof(double));
  for (j=0;j<ncom;j++) xt[j]=pcom[j]+x*xicom[j];
  f=(*nrfunc)(xt);
  free(xt);
  return f;
}

double df1dim(double x)
{

  int j;
  double df1=0.0;
  double *xt,*df;

  xt= (double *) malloc(ncom*sizeof(double));
  df= (double *) malloc(ncom*sizeof(double));
  for (j=0;j<ncom;j++) xt[j]=pcom[j]+x*xicom[j];
  (*nrdfun)(xt,df);
  for (j=0;j<ncom;j++) df1 += df[j]*xicom[j];
  free(df);
  free(xt);
  return df1;
}



void dlinmin(double *p,double *xi,int n,double *fret,
            double (*func)(double *p),
            void (*dfunc)(double *p,double *xi))
{
  int j;
  double xx,fx,fb,fa,bx,ax;
  double xmin = DBL_MIN;
  void mnbrak();

  ncom=n;
  double *pcom=(double *) malloc(sizeof(double)*n);
  double *xicom=(double *) malloc(sizeof(double)*n);
  nrfunc=func;
  nrdfun=dfunc;
  for (j=0;j<n;j++) {
    pcom[j]=p[j];
    xicom[j]=xi[j];
  }
  ax=0.0;
  xx=1.0;
  mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,f1dim);
  *fret=dbrent(ax,xx,bx,f1dim,df1dim,TOL,&xmin);
  for (j=0;j<n;j++) {
    xi[j] *= xmin;
    p[j] += xi[j];
  }
  free(xicom);free(pcom);
}
#define NRANSI
#define TOL 2.0e-4

int gainoff=0;
PIOLONG nbolo;
PIOLONG GAINSTEP;
PIOLONG GAINSTEPADU;

void linmin(double *p, double *xi, int n, double *fret, 
            double (*func)(double *))
{
  int j;
  double xx,xmin,fx,fb,fa,bx,ax;
  
  ncom=n;
  pcom= (double *)malloc(n*sizeof(double));
  xicom= (double *)malloc(n*sizeof(double));
  nrfunc=func;
  for (j=0;j<n;j++) {
    pcom[j]=p[j];
    xicom[j]=xi[j];
  }
  ax=0.0;
  xx=1.0;
  // peut-etre ici que j'evite de me taper les >>1%
  double total=1;
  while (total>1E-3) {
    total=0;
    for (j=0;j<GAINSTEP*nbolo;j++) {
      total=(p[j+gainoff]+xx*xi[j+gainoff]-1)*(p[j+gainoff]+xx*xi[j+gainoff]-1);
    }
    total/=((double) n);
    if (total>1E-3) xx/=2;
  }
  
  mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,f1dim);
  *fret=brent(ax,xx,bx,f1dim,TOL,&xmin);

  // autre endroit ou j'evite de me taper les >>1%
  total=1;
  while (total>1E-3) {
    total=0;
    for (j=0;j<GAINSTEP*nbolo;j++) {
      total=(xmin*xi[j+gainoff])*(xmin*xi[j+gainoff]);
    }
    total/=((double) n);
    if (total>1E-3) xmin/=2;
  }

  for (j=0;j<n;j++) {
    xi[j] *= xmin;
    p[j] += xi[j];
  }
  free(xicom);
  free(pcom);
}
#undef TOL
#undef NRANSI

#undef TOL
PIOLONG *newnr;

void frprmn(double *p,int n,double ftol,int *iter,double *fret,
            double (*func)(double *p),
            void  (*dfunc)(double *p,double *xi))
{
  int j,its;
  double gg,gam,fp,dgg;
  double *g = (double *) malloc(sizeof(double)*n);
  double *h = (double *) malloc(sizeof(double)*n);
  double *xi= (double *) malloc(sizeof(double)*n);

  int rank_size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank_size);
  int rank=rank_size;

  memset(h,0,sizeof(double)*n);
  fp=(*func)(p);
  (*dfunc)(p,xi);
  for (j=0;j<n;j++) {
    g[j] = -xi[j];
    xi[j]=h[j]=g[j];
  }
  for (its=0;its<ITMAX;its++) {
    *iter=its;
    linmin(p,xi,n,fret,func);
    if (2.0*fabs(*fret-fp) <= ftol*(fabs(*fret)+fabs(fp)+EPS)||*fret>fp) {
      free(g);free(h);free(xi);
      if (rank==0) fprintf(stderr,"ABS OK fret[%ld] = %lg %lg %lg %lg %lg %lg\n",(long) its,fp,*fret,p[0],p[1],p[2],p[3]);
      return;
    }
    if (rank==0) fprintf(stderr,"fret[%ld] = %lg %lg %lg %lg %lg %lg : %lg\n",(long) its,fp,*fret,
                         p[0]-1.01,p[1]-1.01,p[2]-1.01,p[3]-1.01,
                          (p[0]-1.01)*(p[0]-1.01)+
                          (p[1]-1.01)*(p[1]-1.01)+
                          (p[2]-1.01)*(p[2]-1.01)+
                          (p[3]-1.01)*(p[3]-1.01));

    fp=(*func)(p);
    (*dfunc)(p,xi);
    dgg=gg=0.0;
    for (j=0;j<n;j++) {
      gg += g[j]*g[j];
      dgg += xi[j]*xi[j];       
      //dgg += (xi[j]+g[j])*xi[j];
    }
    if (gg == 0.0) {
      free(g);free(h);free(xi);
      if (rank==0) fprintf(stderr,"gg==0 fret[%ld] = %lg %lg %lg %lg %lg\n",(long) its,*fret,p[newnr[nbolo]-1],
                           p[newnr[nbolo]],p[newnr[nbolo]+1],p[newnr[nbolo]+2]);
      return;
    }
    gam=dgg/gg;
    for (j=0;j<n;j++) {
      g[j] = -xi[j];
      xi[j]=h[j]=g[j]+gam*h[j];
    }
  }
  free(g);free(h);free(xi);
  //fprintf(stderr,"Too many iterations in FRPRMN\n");
}

#undef ITMAX
#undef EPS
#undef FREEALL


#define ITMAX 2000
#define EPS 1.0e-10

void dfpmin(double *p,int n,double ftol,int *iter,double *fret,double (*func)(),void (*dfunc)())
{
  int rank_size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank_size);
  int rank=rank_size;
  int j,i,its;
  double fp,fae,fad,fac;
  double *xi,*g,*dg,*hdg,*vector();
  double **hessin,**matrix();
  double *p0=(double *)malloc(sizeof(double)*n);
  for (i=0;i<n;i++) p0[i]=2;

  hessin=(double **) malloc(sizeof(double *)*n);
  for (i=0;i<n;i++) hessin[i]=(double *) malloc(sizeof(double)*n);
  
  xi=(double *) malloc(sizeof(double)*n);
  g=(double *) malloc(sizeof(double)*n);
  dg=(double *) malloc(sizeof(double)*n);
  hdg=(double *) malloc(sizeof(double)*n);
  fp=(*func)(p);
  (*dfunc)(p,g);

  for (i=0;i<n;i++) {
    for (j=0;j<n;j++) hessin[i][j]=0.0;
    hessin[i][i]=1.0;
    xi[i] = -g[i];
  }
  for (its=0;its<ITMAX;its++) {
    *iter=its;
    linmin(p,xi,n,fret,func);
    double diffp=0;
    for (i=0;i<n;i++) diffp+=(p0[i]-p[i])*(p0[i]-p[i]);
    diffp=sqrt(diffp/n);

    if (2.0*fabs(*fret-fp) <= ftol*(fabs(*fret)+fabs(fp)+EPS)||*fret<1E-10||diffp<1E-10 ||fabs(*fret-fp)<1E-12) {
      if (rank==0) fprintf(stderr,"END ABS fret[%ld] = %.10lg %.10lg %lg %lg %lg %lg : %lg %lg %lg\n",(long) its,fp,*fret,
                           p[gainoff+0],p[gainoff+1],p[gainoff+2],p[gainoff+3],ftol*(fabs(*fret)+fabs(fp)+EPS),fabs(*fret-fp),diffp);
      free(hdg);
      free(dg);
      free(g);
      free(xi);
      free(p0);
      for (i=0;i<n;i++) free(hessin[i]);
      free(hessin);
      return;
    }
    memcpy(p0,p,n*sizeof(double));

    if (rank==0) fprintf(stderr,"fret[%ld] %lg =\t%.10lg\t%.10lg\t%lg\t[ %lg %lg %lg %lg ]\n",(long) its,diffp,fp,*fret,fp-*fret,
                         p[gainoff+0],p[gainoff+1],p[gainoff+2],p[gainoff+3]);
    fp=(*fret);
    for (i=0;i<n;i++) dg[i]=g[i];
    *fret=(*func)(p);
    (*dfunc)(p,g);
    for (i=0;i<n;i++) dg[i]=g[i]-dg[i];
    for (i=0;i<n;i++) {
      hdg[i]=0.0;
      for (j=0;j<n;j++) hdg[i] += hessin[i][j]*dg[j];
    }
    fac=fae=0.0;
    for (i=0;i<n;i++) {
      fac += dg[i]*xi[i];
      fae += dg[i]*hdg[i];
    }
    fac=1.0/fac;
    fad=1.0/fae;
    for (i=0;i<n;i++) dg[i]=fac*xi[i]-fad*hdg[i];
    for (i=0;i<n;i++)
      for (j=0;j<n;j++)
        hessin[i][j] += fac*xi[i]*xi[j]
          -fad*hdg[i]*hdg[j]+fae*dg[i]*dg[j];
    for (i=0;i<n;i++) {
      xi[i]=0.0;
      for (j=0;j<n;j++) xi[i] -= hessin[i][j]*g[j];
    }
  }
  fprintf(stderr,"Too many iterations in DFPMIN");
}

#undef ITMAX
#undef EPS

// TODO TODO
long nnbpix;
hpix **loc_hpix;
PIOLONG *loc_nhpix;

double *SSI;
double *SSQ;
double *SSU;
double *SSI2;
double *SSQ2;
double *SSU2;
double *II;
double *IQ;
double *IU;
double *QQ;
double *UU;
double *QU;

double *dthetai;
double *dthetaq;
double *dthetau;

double *dii;
double *dqq;
double *duu;

//double *ddegi;
//double *ddegq;
//double *ddegu;

double *dcoi;
double *dcoq;
double *dcou;

double *dfri;
double *dfrq;
double *dfru;

double *ddusti;
double *ddustq;
double *ddustu;

double *dpixi;
double *dpixq;
double *dpixu;

double *cdip;  
double *cco;   
double *ccfree; 
double *ctheta;
double *cdust; 
double *cpix;
double *ctmp;  
double *g;
PIOBYTE *flgpix;

sroll_parContent* Param;

PIOLONG globalBeginRing;
PIOLONG globalEndRing;
PIOLONG globalRangeRing;
// This new struct allow to store the ring distribution (begin/end) among each rank
typedef struct {
  PIOLONG *BeginRing;
  PIOLONG *EndRing;
} rankInfo;

rankInfo globalRankInfo;


PIOINT **rgord;

double *eta_dest;
double *dpsico;
double *dpsisi;

PIOLONG nadu=0;


PIOBYTE **flg_rg;
PIOLONG DEGREE=0;
long npixbeam;
PIOLONG nmatco;
PIOLONG nfreefree;
PIOLONG nmatdust;
long ittt,itbogo;
double *qmat;
double *umat;
long nmatres;
int testzero=-1;

double normaoff=0;
double NEP_tab[20];


#if 0
double funcxi3(double *p)
{ 
  double res=0;
  int i;
  long k,l1,ib,j,m,rrk;
  
  int rank_size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank_size);
  int rank=rank_size;
  int size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  int mpi_size=size;

  for (k=0;k<nnbpix;k++)  {
    long ndata = loc_nhpix[k];
    hpix *htmp = loc_hpix[k];

    if (ndata>0&&flgpix[k]>0) {
      memset(dii ,0,sizeof(double)*nbolo*nadu);
      memset(dqq ,0,sizeof(double)*nbolo*nadu);
      memset(duu ,0,sizeof(double)*nbolo*nadu);

      for (l1=0;l1<ndata;l1++) {
        long ri1=htmp[l1].rg-globalBeginRing;
        long iri1=rgord[htmp[l1].ib][ri1]+newnr[htmp[l1].ib];
#ifdef DORGG
        double g1=g[iri1]; 
#else
        double g1=g[htmp[l1].gi+htmp[l1].ib*GAINSTEP]*htmp[l1].ggi;
#endif
        double CO1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].co
                                          -dpsisi[htmp[l1].ib]*htmp[l1].si);
        double SI1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].si
                                          +dpsisi[htmp[l1].ib]*htmp[l1].co);
        if (flg_rg[htmp[l1].ib][ri1]!=0) {
          
          dii[htmp[l1].adu+htmp[l1].ib*nadu ]+=htmp[l1].w*htmp[l1].dip;
          dqq[htmp[l1].adu+htmp[l1].ib*nadu ]+=htmp[l1].w*htmp[l1].dip*CO1;
          duu[htmp[l1].adu+htmp[l1].ib*nadu ]+=htmp[l1].w*htmp[l1].dip*SI1;
          
        }
      }
    
      for (ib=0;ib<nbolo;ib++) {
        for (j=0;j<nadu;j++) 
          solvemap(dii+j+ib*nadu,dqq+j+ib*nadu,duu+j+ib*nadu,II[k],IQ[k],IU[k],QQ[k],QU[k],UU[k]);
      }
   
      for (l1=0;l1<ndata;l1++) {
        long ri1=htmp[l1].rg-globalBeginRing;
        long iri1=rgord[htmp[l1].ib][ri1]+newnr[htmp[l1].ib];
#ifdef DORGG
        double g1=g[iri1]; 
#else
        double g1=g[htmp[l1].gi+htmp[l1].ib*GAINSTEP]*htmp[l1].ggi;
#endif
        double CO1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].co
                                          -dpsisi[htmp[l1].ib]*htmp[l1].si);
        double SI1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].si
                                          +dpsisi[htmp[l1].ib]*htmp[l1].co);
        if (flg_rg[htmp[l1].ib][ri1]!=0) {
          double ww=1;//htmp[l1].w*htmp[l1].wp/(htmp[l1].w+htmp[l1].wp);
          double tmp=htmp[l1].dip-(SSI[k]+CO1*SSQ[k]+SI1*SSU[k]); 
             

          double divi=NEP_tab[htmp[l1].ib]*htmp[l1].hit*p[newnr[nbolo]+htmp[l1].ib*nadu+htmp[l1].adu];
          divi=divi*divi;
          long l3;

          for (l3=0;l3<ndata;l3++) {
            long ri3=htmp[l3].rg-globalBeginRing;
            if (flg_rg[htmp[l3].ib][ri3]!=0) {
              double sigma=NEP_tab[htmp[l3].ib]*htmp[l3].hit*p[newnr[nbolo]+htmp[l3].ib*nadu+htmp[l3].adu]
                *(htmp[l3].vi+CO1*htmp[l3].vq+SI1*htmp[l3].vu);
              divi+=sigma*sigma;
            }
          }

          for (l3=0;l3<ndata;l3++) {
            long ri3=htmp[l3].rg-globalBeginRing;
            if (flg_rg[htmp[l3].ib][ri3]!=0) {
              long ir3=rgord[htmp[l3].ib][ri3]+newnr[htmp[l3].ib];
              tmp+=p[ir3]*(htmp[l3].vi+CO1*htmp[l3].vq+SI1*htmp[l3].vu);
            }
          }
          tmp-=p[iri1];

          double tmptest=0;
          for (ib=0;ib<nbolo;ib++) 
            for (j=0;j<nadu;j++) {
              tmptest+=p[newnr[nbolo]+ib*nadu+j]*(dii[j+ib*nadu]+CO1*dqq[j+ib*nadu]+SI1*duu[j+ib*nadu]);
            }
          
          tmp+=tmptest-p[newnr[nbolo]+htmp[l1].ib*nadu+htmp[l1].adu]*htmp[l1].thsig*htmp[l1].ggi;

          if (testzero==0) {
            fprintf(stderr,"testzero %lg %lg %lg %lg %lg %lg  %lg \t:: %lg %lg %lg \t:: %lg %lg %lg \n",
                    p[newnr[nbolo]+htmp[l1].ib*nadu+htmp[l1].adu]*htmp[l1].thsig*htmp[l1].ggi,htmp[l1].dip,htmp[l1].dip/(1+1E-3),htmp[l1].sig,
                    p[newnr[nbolo]+htmp[l1].ib*nadu+htmp[l1].adu]*htmp[l1].thsig*htmp[l1].ggi-htmp[l1].dip,htmp[l1].dip/(1+1E-3)-htmp[l1].sig,
                    p[newnr[nbolo]+htmp[l1].ib*nadu+htmp[l1].adu]/(1+1E-3)-1.,
                    tmptest,(SSI[k]+CO1*SSQ[k]+SI1*SSU[k]),
                    tmptest-(SSI[k]+CO1*SSQ[k]+SI1*SSU[k]),
                    p[newnr[nbolo]+htmp[l1].ib*nadu+htmp[l1].adu]*htmp[l1].thsig*htmp[l1].ggi-tmptest,
                    htmp[l1].dip-(SSI[k]+CO1*SSQ[k]+SI1*SSU[k]),
                    p[newnr[nbolo]+htmp[l1].ib*nadu+htmp[l1].adu]*htmp[l1].thsig*htmp[l1].ggi-tmptest -
                    (htmp[l1].dip-(SSI[k]+CO1*SSQ[k]+SI1*SSU[k])));
            
            testzero++;     
          }

          //if (testzero==0) {
          //  fprintf(stderr,"%lg %lg %lg %lg\n",htmp[l1].dip-(SSI[k]+CO1*SSQ[k]+SI1*SSU[k]),tmp,
          //        htmp[l1].dip,p[newnr[nbolo]+htmp[l1].ib*nadu+htmp[l1].adu]*htmp[l1].thsig*htmp[l1].ggi);
          //  testzero++;
          //}
          
          res+=ww*tmp*tmp/divi;
        }
      }
    }
  } 

  if (rank==0) {
    MPI_Status statu;
    double lb;
    for (rrk=1;rrk<mpi_size;rrk++) {
      MPI_Recv(&lb,sizeof(double), MPI_BYTE, rrk,1053, MPI_COMM_WORLD,&statu);
      res+=lb;
    }
  }
  else {
    MPI_Send(&res, sizeof(double), MPI_BYTE, 0, 1053, MPI_COMM_WORLD);
  }
  MPI_Bcast(&res,sizeof(double), MPI_BYTE, 0, MPI_COMM_WORLD);

  /// Average offset equal to zero
  double tmpo=0;
  for (i=0;i<newnr[nbolo];i++) tmpo+=p[i];
  res+=normaoff*tmpo*tmpo;

  return(res);

}


void dfuncxi3(double *p,double *xi)
{
  
  double *ccdip = (double *) malloc(sizeof(double)*nadu*nbolo);

  memset(xi,0,sizeof(double)*nmatres);
  int i;
  long k,l1,ib,j,m,rrk;
  
  int rank_size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank_size);
  int rank=rank_size;
  int size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  int mpi_size=size;

  int testzero=0;
  for (k=0;k<nnbpix;k++)  {
    long ndata = loc_nhpix[k];
    hpix *htmp = loc_hpix[k];

    if (ndata>0&&flgpix[k]>0) {
      memset(dii ,0,sizeof(double)*nbolo*nadu);
      memset(dqq ,0,sizeof(double)*nbolo*nadu);
      memset(duu ,0,sizeof(double)*nbolo*nadu);

      for (l1=0;l1<ndata;l1++) {
        long ri1=htmp[l1].rg-globalBeginRing;
        long iri1=rgord[htmp[l1].ib][ri1]+newnr[htmp[l1].ib];
#ifdef DORGG
        double g1=g[iri1]; 
#else
        double g1=g[htmp[l1].gi+htmp[l1].ib*GAINSTEP]*htmp[l1].ggi;
#endif
        double CO1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].co
                                          -dpsisi[htmp[l1].ib]*htmp[l1].si);
        double SI1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].si
                                          +dpsisi[htmp[l1].ib]*htmp[l1].co);
        if (flg_rg[htmp[l1].ib][ri1]!=0) {
          
          dii[htmp[l1].adu+htmp[l1].ib*nadu ]+=htmp[l1].w*htmp[l1].thsig*htmp[l1].ggi;
          dqq[htmp[l1].adu+htmp[l1].ib*nadu ]+=htmp[l1].w*htmp[l1].thsig*CO1*htmp[l1].ggi;
          duu[htmp[l1].adu+htmp[l1].ib*nadu ]+=htmp[l1].w*htmp[l1].thsig*SI1*htmp[l1].ggi;
          
        }
      }
    
      for (ib=0;ib<nbolo;ib++) {
        for (j=0;j<nadu;j++) 
          solvemap(dii+j+ib*nadu,dqq+j+ib*nadu,duu+j+ib*nadu,II[k],IQ[k],IU[k],QQ[k],QU[k],UU[k]);
      }

      memset(ccdip,0,sizeof(double)*nadu*nbolo);
   
      for (l1=0;l1<ndata;l1++) {
        long ri1=htmp[l1].rg-globalBeginRing;
        long iri1=rgord[htmp[l1].ib][ri1]+newnr[htmp[l1].ib];
#ifdef DORGG
        double g1=g[iri1]; 
#else
        double g1=g[htmp[l1].gi+htmp[l1].ib*GAINSTEP]*htmp[l1].ggi;
#endif
        double CO1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].co
                                          -dpsisi[htmp[l1].ib]*htmp[l1].si);
        double SI1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].si
                                          +dpsisi[htmp[l1].ib]*htmp[l1].co);
        if (flg_rg[htmp[l1].ib][ri1]!=0) {
          double ww=1;//htmp[l1].w*htmp[l1].wp/(htmp[l1].w+htmp[l1].wp);
          double tmp=htmp[l1].dip-(SSI[k]+CO1*SSQ[k]+SI1*SSU[k]); 
        
          double divi=NEP_tab[htmp[l1].ib]*htmp[l1].hit*p[newnr[nbolo]+htmp[l1].ib*nadu+htmp[l1].adu];
          divi=divi*divi;
          long l3;

          for (l3=0;l3<ndata;l3++) {
            long ri3=htmp[l3].rg-globalBeginRing;
            if (flg_rg[htmp[l3].ib][ri3]!=0) {
              double sigma=NEP_tab[htmp[l3].ib]*htmp[l3].hit*p[newnr[nbolo]+htmp[l3].ib*nadu+htmp[l3].adu]*(htmp[l3].vi+CO1*htmp[l3].vq+SI1*htmp[l3].vu);
              ccdip[htmp[l3].ib*nadu+htmp[l3].adu]+=NEP_tab[htmp[l3].ib]*htmp[l3].hit*(htmp[l3].vi+CO1*htmp[l3].vq+SI1*htmp[l3].vu);
              divi+=sigma*sigma;
            }
          }

          for (l3=0;l3<ndata;l3++) {
            long ri3=htmp[l3].rg-globalBeginRing;
            if (flg_rg[htmp[l3].ib][ri3]!=0) {
              long ir3=rgord[htmp[l3].ib][ri3]+newnr[htmp[l3].ib];
              tmp+=p[ir3]*(htmp[l3].vi+CO1*htmp[l3].vq+SI1*htmp[l3].vu);
            }
          }
          tmp-=p[iri1];
        
          for (ib=0;ib<nbolo;ib++) 
            for (j=0;j<nadu;j++) {
              tmp+=p[newnr[nbolo]+ib*nadu+j]*(dii[j+ib*nadu]+CO1*dqq[j+ib*nadu]+SI1*duu[j+ib*nadu]);
            }
          tmp-=p[newnr[nbolo]+htmp[l1].ib*nadu+htmp[l1].adu]*htmp[l1].thsig*htmp[l1].ggi;

          for (l3=0;l3<ndata;l3++) {
            long ri3=htmp[l3].rg-globalBeginRing;
            if (flg_rg[htmp[l3].ib][ri3]!=0) {
              long ir3=rgord[htmp[l3].ib][ri3]+newnr[htmp[l3].ib];
              xi[ir3]+=2*(htmp[l3].vi+CO1*htmp[l3].vq+SI1*htmp[l3].vu)*tmp/divi;
            }
          }
          xi[iri1]-=2*ww*tmp/divi;

          for (ib=0;ib<nbolo;ib++) 
            for (j=0;j<nadu;j++) {
              xi[newnr[nbolo]+ib*nadu+j]+=
                ww*(2*(dii[j+ib*nadu]+CO1*dqq[j+ib*nadu]+SI1*duu[j+ib*nadu])*tmp*divi
                    -2*p[newnr[nbolo]+ib*nadu+j]*ccdip[ib*nadu+j]*ccdip[ib*nadu+j]*tmp*tmp)/divi/divi;
            }
          xi[newnr[nbolo]+htmp[l1].ib*nadu+htmp[l1].adu]
            -=ww*(2*htmp[l1].thsig*htmp[l1].ggi*tmp*divi
               +2*p[newnr[nbolo]+htmp[l1].ib*nadu+htmp[l1].adu]*tmp*tmp)/divi/divi;
        }
      }
    }
  } 

  if (rank==0) {
    MPI_Status statu;
    double *lb = (double *) malloc(sizeof(double)*nmatres);
    for (rrk=1;rrk<mpi_size;rrk++) {
      MPI_Recv(lb,sizeof(double)*nmatres, MPI_BYTE, rrk,1053, MPI_COMM_WORLD,&statu);
      for (i=0;i<nmatres;i++) xi[i]+=lb[i];
    }
    free(lb);
  }
  else {
    MPI_Send(xi, sizeof(double)*nmatres, MPI_BYTE, 0, 1053, MPI_COMM_WORLD);
  }
  
  MPI_Bcast(xi,sizeof(double)*nmatres, MPI_BYTE, 0, MPI_COMM_WORLD);

  free(ccdip);

  /// Average offset equal to zero
  double tmpo=0;
  for (i=0;i<newnr[nbolo];i++) tmpo+=p[i];
  for (i=0;i<newnr[nbolo];i++) xi[i]+=2*normaoff*tmpo;

  return;
}
#endif

#if 0
void Force_Dip(double *p)
{
 
  double res=0;
  int i;
  long k,l1,ib,j,m,rrk;
  
  int rank_size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank_size);
  int rank=rank_size;
  int size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  int mpi_size=size;

  // 2 steps : 
  // - Compute maps
  // - Compute detector maps
  // - Align detector on dipole


  for (k=0;k<nnbpix;k++)  {
    long ndata = loc_nhpix[k];
    hpix *htmp = loc_hpix[k];

    double si=0,sq=0,su=0;
    double sii=0,siq=0,siu=0,sqq=0,suu=0,squ=0;

    if (ndata>0&&flgpix[k]>0) {
      for (l1=0;l1<ndata;l1++) {
        long idg=rgord[htmp[l1].ib][ri]+newnr[htmp[l1].ib];
        ri=rgord[ib][ri]+newnr[htmp[l1].ib];
        
        double gg2=p[newnr[nbolo]+htmp[l1].adu+htmp[l1].ib*GAINSTEP];
        
        double sig_corr = (htmp[l1].sig*gg2-htmp[l1].fsl-htmp[l1].dip-gainadu[htmp[l1].adu+htmp[l1].ib*GAINSTEPADU]);
        
        sig_corr-=p[ri];
        
        if (nmatco!=0) {
        sig_corr -=htmp[l1].comap*p[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP)+htmp[l1].ib];
        }
        if (nmatdust!=0) {
          sig_corr -=htmp[l1].dustmap*p[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP)+nmatco+htmp[l1].ib];
        }
        if (nfreefree!=0) {
          sig_corr -=htmp[l1].freefree*p[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP)+nmatco+nmatdust+htmp[l1].ib];
        }
        
        for (m=0;m<npixbeam;m++)  {
          sig_corr -= htmp[l1].listofpix[m]*p[newnr[nbolo]+nbolo*(GAINSTEP)+m*nbolo+htmp[l1].ib];
        }
        
        double CO1=eta[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].co
                                     -dpsisi[htmp[l1].ib]*htmp[l1].si);
        double SI1=eta[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].si
                                     +dpsisi[htmp[l1].ib]*htmp[l1].co);
        
        si +=htmp[l1].w*sig_corr;
        sq +=htmp[l1].w*CO1*sig_corr;
        su +=htmp[l1].w*SI1*sig_corr;
        
        sii +=htmp[l1].w;
        siq +=htmp[l1].w*CO1;
        siu +=htmp[l1].w*SI1;
        sqq +=htmp[l1].w*CO1*CO1;
        squ +=htmp[l1].w*CO1*SI1;
        suu +=htmp[l1].w*SI1*SI1;
      }

      
    }
  } 

  if (rank==0) {
    MPI_Status statu;
    double lb;
    for (rrk=1;rrk<mpi_size;rrk++) {
      MPI_Recv(&lb,sizeof(double), MPI_BYTE, rrk,1053, MPI_COMM_WORLD,&statu);
      res+=lb;
    }
  }
  else {
    MPI_Send(&res, sizeof(double), MPI_BYTE, 0, 1053, MPI_COMM_WORLD);
  }
  MPI_Bcast(&res,sizeof(double), MPI_BYTE, 0, MPI_COMM_WORLD);
}
#endif
double funcxi2(double *p)
{
 
  double res=0;
  int i;
  long k,l1,ib,j,m,rrk;
  
  int rank_size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank_size);
  int rank=rank_size;
  int size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  int mpi_size=size;

  for (k=0;k<nnbpix;k++)  {
    long ndata = loc_nhpix[k];
    hpix *htmp = loc_hpix[k];

    if (ndata>0&&flgpix[k]>0) {
      memset(dii ,0,sizeof(double)*nbolo*GAINSTEP);
      memset(dqq ,0,sizeof(double)*nbolo*GAINSTEP);
      memset(duu ,0,sizeof(double)*nbolo*GAINSTEP);

      for (l1=0;l1<ndata;l1++) {
        long ri1=htmp[l1].rg-globalBeginRing;

        double CO1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].co
                                          -dpsisi[htmp[l1].ib]*htmp[l1].si);
        double SI1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].si
                                          +dpsisi[htmp[l1].ib]*htmp[l1].co);
        if (flg_rg[htmp[l1].ib][ri1]!=0) {
          
          dii[htmp[l1].adu+htmp[l1].ib*GAINSTEP ]+=htmp[l1].w*htmp[l1].sig;
          dqq[htmp[l1].adu+htmp[l1].ib*GAINSTEP ]+=htmp[l1].w*htmp[l1].sig*CO1;
          duu[htmp[l1].adu+htmp[l1].ib*GAINSTEP ]+=htmp[l1].w*htmp[l1].sig*SI1;

        }
      }
    
      for (ib=0;ib<nbolo;ib++) {
        for (j=0;j<GAINSTEP;j++) 
          solvemap(dii+j+ib*GAINSTEP,dqq+j+ib*GAINSTEP,duu+j+ib*GAINSTEP,II[k],IQ[k],IU[k],QQ[k],QU[k],UU[k]);
      }
   
      for (l1=0;l1<ndata;l1++) {
        long ri1=htmp[l1].rg-globalBeginRing;
        long iri1=rgord[htmp[l1].ib][ri1]+newnr[htmp[l1].ib];

        double CO1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].co
                                          -dpsisi[htmp[l1].ib]*htmp[l1].si);
        double SI1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].si
                                          +dpsisi[htmp[l1].ib]*htmp[l1].co);
        if (flg_rg[htmp[l1].ib][ri1]!=0) {
          double ww=1;//htmp[l1].w*htmp[l1].wp/(htmp[l1].w+htmp[l1].wp);
          //double ddip=htmp[l1].dip;
          double tmp=htmp[l1].dip+htmp[l1].fsl-(SSI2[k]+CO1*SSQ2[k]+SI1*SSU2[k]); 
             
          double divi=NEP_tab[htmp[l1].ib]*htmp[l1].hit*p[newnr[nbolo]+htmp[l1].ib*GAINSTEP+htmp[l1].adu];
          divi=divi*divi;
          long l3;

          for (l3=0;l3<ndata;l3++) {
            long ri3=htmp[l3].rg-globalBeginRing;
            if (flg_rg[htmp[l3].ib][ri3]!=0) {
              double sigma=NEP_tab[htmp[l3].ib]*htmp[l3].hit*p[newnr[nbolo]+htmp[l3].ib*GAINSTEP+htmp[l3].adu]
                *(htmp[l3].vi+CO1*htmp[l3].vq+SI1*htmp[l3].vu);
              divi+=sigma*sigma;
            }
          }

          if (divi==0) {
            ww=0;divi=1;
          }
 
          for (l3=0;l3<ndata;l3++) {
            long ri3=htmp[l3].rg-globalBeginRing;
            if (flg_rg[htmp[l3].ib][ri3]!=0) {
              long ir3=rgord[htmp[l3].ib][ri3]+newnr[htmp[l3].ib];
              tmp-=p[ir3]*(htmp[l3].vi+CO1*htmp[l3].vq+SI1*htmp[l3].vu);
            }
          }
          tmp+=p[iri1];

          double tmptest=0;
          for (ib=0;ib<nbolo;ib++) 
            for (j=0;j<GAINSTEP;j++) {
              tmptest+=p[newnr[nbolo]+ib*GAINSTEP+j]*(dii[j+ib*GAINSTEP]+CO1*dqq[j+ib*GAINSTEP]+SI1*duu[j+ib*GAINSTEP]);
            }
          
          tmp+=tmptest-p[newnr[nbolo]+htmp[l1].ib*GAINSTEP+htmp[l1].adu]*htmp[l1].sig;

#if 1
          double tmpco=0;
          double tmpdust=0;
          double tmpfree=0;

          if (nmatco>0) {
            for (ib=0;ib<nbolo;ib++) tmpco-=p[newnr[nbolo]+nbolo*(DEGREE+GAINSTEP+npixbeam)+ib]*(dcoi[ib+k*nbolo]+CO1*dcoq[ib+k*nbolo]+SI1*dcou[ib+k*nbolo]);
            tmp+=p[newnr[nbolo]+nbolo*(DEGREE+GAINSTEP+npixbeam)+htmp[l1].ib]*htmp[l1].comap;
          }
          tmp+=tmpco;
          if (nmatdust>0) {
            for (ib=0;ib<nbolo;ib++) tmpdust-=p[newnr[nbolo]+nbolo*(DEGREE+GAINSTEP+npixbeam)+nmatco+ib]*(ddusti[ib+k*nbolo]+CO1*ddustq[ib+k*nbolo]+SI1*ddustu[ib+k*nbolo]);
            tmp+=p[newnr[nbolo]+nbolo*(DEGREE+GAINSTEP+npixbeam)+nmatco+htmp[l1].ib]*htmp[l1].dustmap;
          }
          tmp+=tmpdust;
          if (nfreefree>0) {
            for (ib=0;ib<nbolo;ib++) tmpfree-=p[newnr[nbolo]+nbolo*(DEGREE+GAINSTEP+npixbeam)+nmatco+nmatdust+ib]*(dfri[ib+k*nbolo]+CO1*dfrq[ib+k*nbolo]+SI1*dfru[ib+k*nbolo]);
            tmp+=p[newnr[nbolo]+nbolo*(DEGREE+GAINSTEP+npixbeam)+nmatco+nmatdust+htmp[l1].ib]*htmp[l1].freefree;
          }
          tmp+=tmpfree;

          for (m=0;m<npixbeam;m++) {
            for (ib=0;ib<nbolo;ib++) tmp-=p[newnr[nbolo]+(DEGREE+GAINSTEP)*nbolo+m*nbolo+ib]
                                       *(dpixi[ib+m*nbolo+k*nbolo*npixbeam]
                                         +CO1*dpixq[ib+m*nbolo+k*nbolo*npixbeam]
                                         +SI1*dpixu[ib+m*nbolo+k*nbolo*npixbeam]);

            tmp+=p[newnr[nbolo]+(DEGREE+GAINSTEP)*nbolo+m*nbolo+htmp[l1].ib]*htmp[l1].listofpix[m];
          }
#endif
          res+=ww*tmp*tmp/divi;
        }
      }
    }
  } 

  if (rank==0) {
    MPI_Status statu;
    double lb;
    for (rrk=1;rrk<mpi_size;rrk++) {
      MPI_Recv(&lb,sizeof(double), MPI_BYTE, rrk,1053, MPI_COMM_WORLD,&statu);
      res+=lb;
    }
  }
  else {
    MPI_Send(&res, sizeof(double), MPI_BYTE, 0, 1053, MPI_COMM_WORLD);
  }
  MPI_Bcast(&res,sizeof(double), MPI_BYTE, 0, MPI_COMM_WORLD);

  /// Average offset equal to zero
  double tmpo=0;
  for (i=0;i<newnr[nbolo];i++) tmpo+=p[i];
  res+=normaoff*tmpo*tmpo;

#if 1
  // DO CO average =0
  if (nmatco>0) {
    tmpo=0;
    for (ib=0;ib<nbolo;ib++) tmpo+=p[newnr[nbolo]+nbolo*(DEGREE+nadu+npixbeam)+ib];
    res+=normaoff*tmpo*tmpo;
    //if (rank==0) fprintf(stderr,"normaoff %lg %lg \n",normaoff,normaoff*tmpo*tmpo);
  }
  // DO DUST average =0
  if (nmatdust>0) {
    tmpo=0;
    for (ib=0;ib<nbolo;ib++) tmpo+=p[newnr[nbolo]+nbolo*(DEGREE+nadu+npixbeam)+nmatco+ib];
    res+=normaoff*tmpo*tmpo;
    //if (rank==0) fprintf(stderr,"normaoff %lg %lg \n",normaoff,normaoff*tmpo*tmpo);
  }
  // DO FREE FREE average =0
  if (nfreefree>0) {
    tmpo=0;
    for (ib=0;ib<nbolo;ib++) tmpo+=p[newnr[nbolo]+nbolo*(DEGREE+nadu+npixbeam)+nmatco+nmatdust+ib];
    res+=normaoff*tmpo*tmpo;
    //if (rank==0) fprintf(stderr,"normaoff %lg %lg \n",normaoff,normaoff*tmpo*tmpo);
  }
#endif
  return(res);

}



void dfuncxi2_light(double *p,double *xi)
{
 
  int i;
  long k,l1,ib,j,m,rrk;
  
  int rank_size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank_size);
  int rank=rank_size;
  int size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  int mpi_size=size;
  
  double *ccdip = (double *) malloc(sizeof(double)*GAINSTEP*nbolo);

  memset(xi,0,sizeof(double)*(newnr[nbolo]+nbolo*npixbeam+nbolo*(GAINSTEP)+nmatco+nmatdust+nfreefree));

  for (k=0;k<nnbpix;k++)  {
    long ndata = loc_nhpix[k];
    hpix *htmp = loc_hpix[k];

    if (ndata>0&&flgpix[k]>0) {
      memset(dii ,0,sizeof(double)*nbolo*GAINSTEP);
      memset(dqq ,0,sizeof(double)*nbolo*GAINSTEP);
      memset(duu ,0,sizeof(double)*nbolo*GAINSTEP);

      for (l1=0;l1<ndata;l1++) {
        long ri1=htmp[l1].rg-globalBeginRing;

        double CO1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].co
                                          -dpsisi[htmp[l1].ib]*htmp[l1].si);
        double SI1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].si
                                          +dpsisi[htmp[l1].ib]*htmp[l1].co);
        if (flg_rg[htmp[l1].ib][ri1]!=0) {
          
          dii[htmp[l1].adu+htmp[l1].ib*GAINSTEP ]+=htmp[l1].w*htmp[l1].sig;
          dqq[htmp[l1].adu+htmp[l1].ib*GAINSTEP ]+=htmp[l1].w*htmp[l1].sig*CO1;
          duu[htmp[l1].adu+htmp[l1].ib*GAINSTEP ]+=htmp[l1].w*htmp[l1].sig*SI1;
          
        }
      }
    
      for (ib=0;ib<nbolo;ib++) {
        for (j=0;j<GAINSTEP;j++) 
          solvemap(dii+j+ib*GAINSTEP,dqq+j+ib*GAINSTEP,duu+j+ib*GAINSTEP,II[k],IQ[k],IU[k],QQ[k],QU[k],UU[k]);
      }
   
      memset(ccdip,0,sizeof(double)*GAINSTEP*nbolo);

      for (l1=0;l1<ndata;l1++) {
        long ri1=htmp[l1].rg-globalBeginRing;
        long iri1=rgord[htmp[l1].ib][ri1]+newnr[htmp[l1].ib];

        double CO1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].co
                                          -dpsisi[htmp[l1].ib]*htmp[l1].si);
        double SI1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].si
                                          +dpsisi[htmp[l1].ib]*htmp[l1].co);
        if (flg_rg[htmp[l1].ib][ri1]!=0) {
          double ww=1;//htmp[l1].w*htmp[l1].wp/(htmp[l1].w+htmp[l1].wp);
          double tmp=htmp[l1].dip+htmp[l1].fsl-(SSI2[k]+CO1*SSQ2[k]+SI1*SSU2[k]); 
             
          double divi=NEP_tab[htmp[l1].ib]*htmp[l1].hit*p[newnr[nbolo]+htmp[l1].ib*GAINSTEP+htmp[l1].adu];
          divi=divi*divi;
          long l3;

          for (l3=0;l3<ndata;l3++) {
            long ri3=htmp[l3].rg-globalBeginRing;
            if (flg_rg[htmp[l3].ib][ri3]!=0) {
              double sigma=NEP_tab[htmp[l3].ib]*htmp[l3].hit*p[newnr[nbolo]+htmp[l3].ib*GAINSTEP+htmp[l3].adu]
                *(htmp[l3].vi+CO1*htmp[l3].vq+SI1*htmp[l3].vu);
              ccdip[htmp[l3].ib*GAINSTEP+htmp[l3].adu]+=NEP_tab[htmp[l3].ib]*htmp[l3].hit*(htmp[l3].vi+CO1*htmp[l3].vq+SI1*htmp[l3].vu);
              divi+=sigma*sigma;
            }
          }

          if (divi==0) {
            ww=0;divi=1;
          }

          for (l3=0;l3<ndata;l3++) {
            long ri3=htmp[l3].rg-globalBeginRing;
            if (flg_rg[htmp[l3].ib][ri3]!=0) {
              long ir3=rgord[htmp[l3].ib][ri3]+newnr[htmp[l3].ib];
              tmp-=p[ir3]*(htmp[l3].vi+CO1*htmp[l3].vq+SI1*htmp[l3].vu);
            }
          }
          tmp+=p[iri1];

          double tmptest=0;
          for (ib=0;ib<nbolo;ib++) 
            for (j=0;j<GAINSTEP;j++) {
              tmptest+=p[newnr[nbolo]+ib*GAINSTEP+j]*(dii[j+ib*GAINSTEP]+CO1*dqq[j+ib*GAINSTEP]+SI1*duu[j+ib*GAINSTEP]);
            }
          
          tmp+=tmptest-p[newnr[nbolo]+htmp[l1].ib*GAINSTEP+htmp[l1].adu]*htmp[l1].sig;


          double tmpco=0;
          double tmpdust=0;
          double tmpfree=0;

          if (nmatco>0) {
            for (ib=0;ib<nbolo;ib++) tmpco-=p[newnr[nbolo]+nbolo*(DEGREE+GAINSTEP+npixbeam)+ib]*(dcoi[ib+k*nbolo]+CO1*dcoq[ib+k*nbolo]+SI1*dcou[ib+k*nbolo]);
            tmp+=p[newnr[nbolo]+nbolo*(DEGREE+GAINSTEP+npixbeam)+htmp[l1].ib]*htmp[l1].comap;
          }
          tmp+=tmpco;
          if (nmatdust>0) {
            for (ib=0;ib<nbolo;ib++) tmpdust-=p[newnr[nbolo]+nbolo*(DEGREE+GAINSTEP+npixbeam)+nmatco+ib]*(ddusti[ib+k*nbolo]+CO1*ddustq[ib+k*nbolo]+SI1*ddustu[ib+k*nbolo]);
            tmp+=p[newnr[nbolo]+nbolo*(DEGREE+GAINSTEP+npixbeam)+nmatco+htmp[l1].ib]*htmp[l1].dustmap;
          }
          tmp+=tmpdust;
          if (nfreefree>0) {
            for (ib=0;ib<nbolo;ib++) tmpfree-=p[newnr[nbolo]+nbolo*(DEGREE+GAINSTEP+npixbeam)+nmatco+nmatdust+ib]*(dfri[ib+k*nbolo]+CO1*dfrq[ib+k*nbolo]+SI1*dfru[ib+k*nbolo]);
            tmp+=p[newnr[nbolo]+nbolo*(DEGREE+GAINSTEP+npixbeam)+nmatco+nmatdust+htmp[l1].ib]*htmp[l1].freefree;
          }
          tmp+=tmpfree;

          for (m=0;m<npixbeam;m++) {
            for (ib=0;ib<nbolo;ib++) tmp-=p[newnr[nbolo]+(DEGREE+GAINSTEP)*nbolo+m*nbolo+ib]
                                       *(dpixi[ib+m*nbolo+k*nbolo*npixbeam]
                                         +CO1*dpixq[ib+m*nbolo+k*nbolo*npixbeam]
                                         +SI1*dpixu[ib+m*nbolo+k*nbolo*npixbeam]);

            tmp+=p[newnr[nbolo]+(DEGREE+GAINSTEP)*nbolo+m*nbolo+htmp[l1].ib]*htmp[l1].listofpix[m];
          }

          /*========================================================================================
            AND now derivate the gain only
            ========================================================================================*/

          for (ib=0;ib<nbolo;ib++) 
            for (j=0;j<GAINSTEP;j++) {
              xi[newnr[nbolo]+ib*GAINSTEP+j]+=
                ww*(2*(dii[j+ib*GAINSTEP]+CO1*dqq[j+ib*GAINSTEP]+SI1*duu[j+ib*GAINSTEP])*tmp*divi
                    -2*p[newnr[nbolo]+ib*GAINSTEP+j]*ccdip[ib*GAINSTEP+j]*ccdip[ib*GAINSTEP+j]*tmp*tmp)/divi/divi;
            }
          xi[newnr[nbolo]+htmp[l1].ib*GAINSTEP+htmp[l1].adu]
            -=ww*(2*htmp[l1].sig*tmp*divi
                  +2*p[newnr[nbolo]+htmp[l1].ib*GAINSTEP+htmp[l1].adu]*ccdip[htmp[l1].ib*GAINSTEP+htmp[l1].adu]
                  *ccdip[htmp[l1].ib*GAINSTEP+htmp[l1].adu]*tmp*tmp)/divi/divi;

          
        }
      }
    }
  } 

  if (rank==0) {
    MPI_Status statu;
    double *lb = (double *) malloc(sizeof(double)*(newnr[nbolo]+nbolo*npixbeam+nbolo*(GAINSTEP)+nmatco+nmatdust+nfreefree));
    for (rrk=1;rrk<mpi_size;rrk++) {
      MPI_Recv(lb,sizeof(double)*(newnr[nbolo]+nbolo*npixbeam+nbolo*(GAINSTEP)+nmatco+nmatdust+nfreefree), MPI_BYTE, rrk,1053, MPI_COMM_WORLD,&statu);
      for (i=0;i<(newnr[nbolo]+nbolo*npixbeam+nbolo*(GAINSTEP)+nmatco+nmatdust+nfreefree);i++) xi[i]+=lb[i];
    }
    free(lb);
  }
  else {
    MPI_Send(xi, sizeof(double)*(newnr[nbolo]+nbolo*npixbeam+nbolo*(GAINSTEP)+nmatco+nmatdust+nfreefree), MPI_BYTE, 0, 1053, MPI_COMM_WORLD);
  }
  
  MPI_Bcast(xi,sizeof(double)*(newnr[nbolo]+nbolo*npixbeam+nbolo*(GAINSTEP)+nmatco+nmatdust+nfreefree), MPI_BYTE, 0, MPI_COMM_WORLD);

  
  free(ccdip);

}

void dfuncxi2(double *p,double *xi)
{
  
  double *ccdip = (double *) malloc(sizeof(double)*GAINSTEP*nbolo);

  memset(xi,0,sizeof(double)*(newnr[nbolo]+nbolo*npixbeam+nbolo*(GAINSTEP)+nmatco+nmatdust+nfreefree));
  int i;
  long k,l1,ib,j,m,rrk;
  
  int rank_size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank_size);
  int rank=rank_size;
  int size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  int mpi_size=size;

  for (k=0;k<nnbpix;k++)  {
    long ndata = loc_nhpix[k];
    hpix *htmp = loc_hpix[k];

    if (ndata>0&&flgpix[k]>0) {
      memset(dii ,0,sizeof(double)*nbolo*GAINSTEP);
      memset(dqq ,0,sizeof(double)*nbolo*GAINSTEP);
      memset(duu ,0,sizeof(double)*nbolo*GAINSTEP);

      for (l1=0;l1<ndata;l1++) {
        long ri1=htmp[l1].rg-globalBeginRing;
        double CO1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].co
                                          -dpsisi[htmp[l1].ib]*htmp[l1].si);
        double SI1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].si
                                          +dpsisi[htmp[l1].ib]*htmp[l1].co);
        if (flg_rg[htmp[l1].ib][ri1]!=0) {
          
          dii[htmp[l1].adu+htmp[l1].ib*GAINSTEP ]+=htmp[l1].w*htmp[l1].sig;
          dqq[htmp[l1].adu+htmp[l1].ib*GAINSTEP ]+=htmp[l1].w*htmp[l1].sig*CO1;
          duu[htmp[l1].adu+htmp[l1].ib*GAINSTEP ]+=htmp[l1].w*htmp[l1].sig*SI1;
          
        }
      }
    
      for (ib=0;ib<nbolo;ib++) {
        for (j=0;j<GAINSTEP;j++) 
          solvemap(dii+j+ib*GAINSTEP,dqq+j+ib*GAINSTEP,duu+j+ib*GAINSTEP,II[k],IQ[k],IU[k],QQ[k],QU[k],UU[k]);
      }

      memset(ccdip,0,sizeof(double)*GAINSTEP*nbolo);
   
      for (l1=0;l1<ndata;l1++) {
        long ri1=htmp[l1].rg-globalBeginRing;
        long iri1=rgord[htmp[l1].ib][ri1]+newnr[htmp[l1].ib];

        double CO1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].co
                                          -dpsisi[htmp[l1].ib]*htmp[l1].si);
        double SI1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].si
                                          +dpsisi[htmp[l1].ib]*htmp[l1].co);
        if (flg_rg[htmp[l1].ib][ri1]!=0) {
          double ww=1;//htmp[l1].w*htmp[l1].wp/(htmp[l1].w+htmp[l1].wp);
          double tmp=htmp[l1].dip+htmp[l1].fsl-(SSI2[k]+CO1*SSQ2[k]+SI1*SSU2[k]); 
        
          double divi=NEP_tab[htmp[l1].ib]*htmp[l1].hit*p[newnr[nbolo]+htmp[l1].ib*GAINSTEP+htmp[l1].adu];
          divi=divi*divi;
          long l3;

          for (l3=0;l3<ndata;l3++) {
            long ri3=htmp[l3].rg-globalBeginRing;
            if (flg_rg[htmp[l3].ib][ri3]!=0) {
              double sigma=NEP_tab[htmp[l3].ib]*htmp[l3].hit*p[newnr[nbolo]+htmp[l3].ib*GAINSTEP+htmp[l3].adu]
                *(htmp[l3].vi+CO1*htmp[l3].vq+SI1*htmp[l3].vu);
              ccdip[htmp[l3].ib*GAINSTEP+htmp[l3].adu]+=NEP_tab[htmp[l3].ib]*htmp[l3].hit*(htmp[l3].vi+CO1*htmp[l3].vq+SI1*htmp[l3].vu);
              divi+=sigma*sigma;
            }
          }

          for (l3=0;l3<ndata;l3++) {
            long ri3=htmp[l3].rg-globalBeginRing;
            if (flg_rg[htmp[l3].ib][ri3]!=0) {
              long ir3=rgord[htmp[l3].ib][ri3]+newnr[htmp[l3].ib];
              tmp-=p[ir3]*(htmp[l3].vi+CO1*htmp[l3].vq+SI1*htmp[l3].vu);
            }
          }
          tmp+=p[iri1];
        
          for (ib=0;ib<nbolo;ib++) 
            for (j=0;j<GAINSTEP;j++) {
              tmp+=p[newnr[nbolo]+ib*GAINSTEP+j]*(dii[j+ib*GAINSTEP]+CO1*dqq[j+ib*GAINSTEP]+SI1*duu[j+ib*GAINSTEP]);
            }
          tmp-=p[newnr[nbolo]+htmp[l1].ib*GAINSTEP+htmp[l1].adu]*htmp[l1].sig;

          
          if (nmatco>0) {
            for (ib=0;ib<nbolo;ib++) tmp-=p[newnr[nbolo]+nbolo*(DEGREE+GAINSTEP+npixbeam)+ib]*(dcoi[ib+k*nbolo]+CO1*dcoq[ib+k*nbolo]+SI1*dcou[ib+k*nbolo]);
            tmp+=p[newnr[nbolo]+nbolo*(DEGREE+GAINSTEP+npixbeam)+htmp[l1].ib]*htmp[l1].comap;
          }

          if (nmatdust>0) {
            for (ib=0;ib<nbolo;ib++) tmp-=p[newnr[nbolo]+nbolo*(DEGREE+GAINSTEP+npixbeam)+nmatco+ib]*(ddusti[ib+k*nbolo]+CO1*ddustq[ib+k*nbolo]+SI1*ddustu[ib+k*nbolo]);
            tmp+=p[newnr[nbolo]+nbolo*(DEGREE+GAINSTEP+npixbeam)+nmatco+htmp[l1].ib]*htmp[l1].dustmap;
          }
          if (nfreefree>0) {
            for (ib=0;ib<nbolo;ib++) tmp-=p[newnr[nbolo]+nbolo*(DEGREE+GAINSTEP+npixbeam)+nmatco+nmatdust+ib]*(dfri[ib+k*nbolo]+CO1*dfrq[ib+k*nbolo]+SI1*dfru[ib+k*nbolo]);
            tmp+=p[newnr[nbolo]+nbolo*(DEGREE+GAINSTEP+npixbeam)+nmatco+nmatdust+htmp[l1].ib]*htmp[l1].freefree;
          }
          for (m=0;m<npixbeam;m++) {
            for (ib=0;ib<nbolo;ib++) tmp-=p[newnr[nbolo]+(DEGREE+GAINSTEP)*nbolo+m*nbolo+ib]
                                       *(dpixi[ib+m*nbolo+k*nbolo*npixbeam]
                                         +CO1*dpixq[ib+m*nbolo+k*nbolo*npixbeam]
                                         +SI1*dpixu[ib+m*nbolo+k*nbolo*npixbeam]);

            tmp+=p[newnr[nbolo]+(DEGREE+GAINSTEP)*nbolo+m*nbolo+htmp[l1].ib]*htmp[l1].listofpix[m];
          }


          for (l3=0;l3<ndata;l3++) {
            long ri3=htmp[l3].rg-globalBeginRing;
            if (flg_rg[htmp[l3].ib][ri3]!=0) {
              long ir3=rgord[htmp[l3].ib][ri3]+newnr[htmp[l3].ib];
              xi[ir3]-=2*(htmp[l3].vi+CO1*htmp[l3].vq+SI1*htmp[l3].vu)*tmp/divi;
            }
          }
          xi[iri1]+=2*ww*tmp/divi;


          for (ib=0;ib<nbolo;ib++) 
            for (j=0;j<GAINSTEP;j++) {
              xi[newnr[nbolo]+ib*GAINSTEP+j]+=
                ww*(2*(dii[j+ib*GAINSTEP]+CO1*dqq[j+ib*GAINSTEP]+SI1*duu[j+ib*GAINSTEP])*tmp*divi
                    -2*p[newnr[nbolo]+ib*GAINSTEP+j]*ccdip[ib*GAINSTEP+j]*ccdip[ib*GAINSTEP+j]*tmp*tmp)/divi/divi;
            }
          xi[newnr[nbolo]+htmp[l1].ib*GAINSTEP+htmp[l1].adu]
            -=ww*(2*htmp[l1].sig*tmp*divi
               +2*p[newnr[nbolo]+htmp[l1].ib*GAINSTEP+htmp[l1].adu]*tmp*tmp)/divi/divi;

          if (nmatco>0) {
            for (ib=0;ib<nbolo;ib++) xi[newnr[nbolo]+nbolo*(DEGREE+GAINSTEP+npixbeam)+ib]
                                       -=2*ww*(dcoi[ib+k*nbolo]+CO1*dcoq[ib+k*nbolo]+SI1*dcou[ib+k*nbolo])*tmp/divi;
            xi[newnr[nbolo]+nbolo*(DEGREE+GAINSTEP+npixbeam)+htmp[l1].ib]
              +=2*ww*htmp[l1].comap*tmp/divi;
          }
            
          if (nmatdust>0) {
            for (ib=0;ib<nbolo;ib++) xi[newnr[nbolo]+nbolo*(DEGREE+GAINSTEP+npixbeam)+nmatco+ib]
                                       -=2*ww*(ddusti[ib+k*nbolo]+CO1*ddustq[ib+k*nbolo]+SI1*ddustu[ib+k*nbolo])*tmp/divi;
            xi[newnr[nbolo]+nbolo*(DEGREE+GAINSTEP+npixbeam)+nmatco+htmp[l1].ib]
              +=2*ww*htmp[l1].dustmap*tmp/divi;
          }

          if (nfreefree>0) {
            for (ib=0;ib<nbolo;ib++) xi[newnr[nbolo]+nbolo*(DEGREE+GAINSTEP+npixbeam)+nmatco+nmatdust+ib]
                                       -=2*ww*(dfri[ib+k*nbolo]+CO1*dfrq[ib+k*nbolo]+SI1*dfru[ib+k*nbolo])*tmp/divi;
            xi[newnr[nbolo]+nbolo*(DEGREE+GAINSTEP+npixbeam)+nmatco+nmatdust+htmp[l1].ib]
              +=2*ww*htmp[l1].freefree*tmp/divi;
          }

          for (m=0;m<npixbeam;m++) {
            for (ib=0;ib<nbolo;ib++) xi[newnr[nbolo]+(DEGREE+GAINSTEP)*nbolo+m*nbolo+ib]
                                       -=2*ww*(dpixi[ib+m*nbolo+k*nbolo*npixbeam]
                                         +CO1*dpixq[ib+m*nbolo+k*nbolo*npixbeam]
                                         +SI1*dpixu[ib+m*nbolo+k*nbolo*npixbeam])*tmp/divi;

            xi[newnr[nbolo]+(DEGREE+GAINSTEP)*nbolo+m*nbolo+htmp[l1].ib]
              +=2*ww*htmp[l1].listofpix[m]*tmp/divi;
          }
        }
      }
    }
  } 

  if (rank==0) {
    MPI_Status statu;
    double *lb = (double *) malloc(sizeof(double)*(newnr[nbolo]+nbolo*npixbeam+nbolo*(GAINSTEP)+nmatco+nmatdust+nfreefree));
    for (rrk=1;rrk<mpi_size;rrk++) {
      MPI_Recv(lb,sizeof(double)*(newnr[nbolo]+nbolo*npixbeam+nbolo*(GAINSTEP)+nmatco+nmatdust+nfreefree), MPI_BYTE, rrk,1053, MPI_COMM_WORLD,&statu);
      for (i=0;i<nmatres;i++) xi[i]+=lb[i];
    }
    free(lb);
  }
  else {
    MPI_Send(xi, sizeof(double)*(newnr[nbolo]+nbolo*npixbeam+nbolo*(GAINSTEP)+nmatco+nmatdust+nfreefree), MPI_BYTE, 0, 1053, MPI_COMM_WORLD);
  }
  
  MPI_Bcast(xi,sizeof(double)*(newnr[nbolo]+nbolo*npixbeam+nbolo*(GAINSTEP)+nmatco+nmatdust+nfreefree), MPI_BYTE, 0, MPI_COMM_WORLD);

  /// Average offset equal to zero
  double tmpo=0;
  for (i=0;i<newnr[nbolo];i++) tmpo+=p[i];
  for (i=0;i<newnr[nbolo];i++) xi[i]+=2*normaoff*tmpo;
  // DO CO average =0
  if (nmatco>0) {
    tmpo=0;
    for (ib=0;ib<nbolo;ib++) tmpo+=p[newnr[nbolo]+nbolo*(GAINSTEP+npixbeam)+ib];
    for (ib=0;ib<nbolo;ib++) xi[newnr[nbolo]+nbolo*(GAINSTEP+npixbeam)+ib]+=2*normaoff*tmpo;
  }
  // DO DUST average =0
  if (nmatdust>0) {
    tmpo=0;
    for (ib=0;ib<nbolo;ib++) tmpo+=p[newnr[nbolo]+nbolo*(GAINSTEP+npixbeam)+ib+nmatco];
    for (ib=0;ib<nbolo;ib++) xi[newnr[nbolo]+nbolo*(GAINSTEP+npixbeam)+ib+nmatco]+=2*normaoff*tmpo;
  }
  // DO FREE FREE average =0
  if (nfreefree>0) {
    tmpo=0;
    for (ib=0;ib<nbolo;ib++) tmpo+=p[newnr[nbolo]+nbolo*(GAINSTEP+npixbeam)+ib+nmatco+nmatdust];
    for (ib=0;ib<nbolo;ib++) xi[newnr[nbolo]+nbolo*(GAINSTEP+npixbeam)+ib+nmatco+nmatdust]+=2*normaoff*tmpo;
  }

  free(ccdip);
  return;
}

long REMHDIP;
long DOFITTHETA=0;
double *x2;
double *x2old;
double *x2init;
double *b2;
double *d2;
double *q2;
double *r2;
double *s2;
double *hit2;
PIOLONG docutrg;

long nmatpix;
#define FITGAIN
#ifdef FITGAIN
  int nitbogo=4;
#else
  int nitbogo=1;
#endif
PIOLONG CUTRG;
PIOLONG *newnr2;
PIOINT **rgord2;
//int *the_stat_pix;
double *xadu;
double delta0;

double **cache_xi2;
double **cache_gain_xi2;
long n_cache_xi2=0;


void minimize(double *ix2,double *gain)
{
  MPI_Status statu;
  long i,rrk,j,k,l1,l2,ib;
  int itermax = 500;
  int iter;
  double  delta_new, alpha=0, delta_old, beta;


  int rank;
  int size;
  int mpi_size;
  int rank_size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank_size);
  rank=rank_size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  mpi_size=size;

  double mintest=1E30;
  k=-1;
  for (i=0;i<n_cache_xi2;i++) {
    double tmptest=0;
    for (j=0;j<nbolo*GAINSTEP;j++) tmptest+=(gain[j]-cache_gain_xi2[i][j])*(gain[j]-cache_gain_xi2[i][j]);
    if (tmptest<mintest) {
      mintest=tmptest;
      k=i;
    }
  }
  if (k!=-1) {
    if (mintest==0) {
      if (rank==0) fprintf(stderr,"FIND CACHE GAIN %lg\n",mintest);
      memcpy(ix2,cache_xi2[k],sizeof(double)*(newnr[nbolo]+nbolo*npixbeam+nmatco+nmatdust+nfreefree));
      return;
    }
  }

  //memset(ix2,0,sizeof(double)*(newnr[nbolo]+nbolo*npixbeam+nmatco+nmatdust+nfreefree));

  if (rank==0) {
    fprintf(stderr,"==============================\n\nITERATION %ld\n\n==============================\n",itbogo);
    fprintf(stderr,"GAIN ");
    for (i=0;i<nbolo;i++) fprintf(stderr,"%lg ",gain[i*GAINSTEP]);
    fprintf(stderr,"\n");
    if (GAINSTEP>1) {
      fprintf(stderr,"...\nGAIN ");
      for (i=0;i<nbolo;i++) fprintf(stderr,"%lg ",gain[i*GAINSTEP+(GAINSTEP-1)]);
      fprintf(stderr,"\n");
    }
  }
  if (itbogo==0) delta0=0;
  MPI_Barrier(MPI_COMM_WORLD);
  
  struct timeval tp1,tp2;
  gettimeofday(&tp1,NULL);

  nmatres=newnr[nbolo]+nbolo*npixbeam+nmatco+nmatdust+nfreefree;
  
  iter = 0;
  memset(b2  ,0,nmatres*sizeof (double));
  memset(d2  ,0,nmatres*sizeof (double));
  memset(q2  ,0,nmatres*sizeof (double));
  memset(r2  ,0,nmatres*sizeof (double));
  memset(s2  ,0,nmatres*sizeof (double));
  memset(hit2,0,nmatres*sizeof (double));

  //==========================================
  //=  Compute second member
  //=
  //=
  long l,m;

  ////// BUILD B2
  //GetProcMem(&vmem,&phymem);
  //if (rank==0) fprintf(stderr,"Rank: %ld Line=%d MEM %.1lf[%.1lf]MB\n",
  //              (long) rank, __LINE__,
  //              (double) vmem/1024./1024.,
  //              (double) phymem/1024./1024.);
  //int ptest=0;
  for (k=0;k<nnbpix;k++)  {
    //long imat=the_stat_pix[k];
    long ndata = loc_nhpix[k];
    hpix *htmp = loc_hpix[k];
    II[k]=0;
    IQ[k]=0;
    IU[k]=0;
    QQ[k]=0;
    UU[k]=0;
    QU[k]=0;
      

    for (l1=0;l1<ndata;l1++) {
      long ri1=htmp[l1].rg-globalBeginRing;
      
      double CO1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].co
                                        -dpsisi[htmp[l1].ib]*htmp[l1].si);
      double SI1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].si
                                        +dpsisi[htmp[l1].ib]*htmp[l1].co);
      if (flg_rg[htmp[l1].ib][ri1]!=0) {
        II[k]+=htmp[l1].w;
        IQ[k]+=htmp[l1].w*CO1;
        IU[k]+=htmp[l1].w*SI1;
        QQ[k]+=htmp[l1].w*CO1*CO1;
        QU[k]+=htmp[l1].w*SI1*CO1;
        UU[k]+=htmp[l1].w*SI1*SI1;
      }
    }

    if (ndata>0&&flgpix[k]>0) {

      double SI=0;
      double SQ=0;
      double SU=0;

      memset(dcoi+k*nbolo,0,sizeof(double)*nbolo);
      memset(dcoq+k*nbolo,0,sizeof(double)*nbolo);
      memset(dcou+k*nbolo,0,sizeof(double)*nbolo);
      memset(dfri+k*nbolo,0,sizeof(double)*nbolo);
      memset(dfrq+k*nbolo,0,sizeof(double)*nbolo);
      memset(dfru+k*nbolo,0,sizeof(double)*nbolo);
      memset(dthetai+k*nbolo,0,sizeof(double)*nbolo);
      memset(dthetaq+k*nbolo,0,sizeof(double)*nbolo);
      memset(dthetau+k*nbolo,0,sizeof(double)*nbolo);
      memset(ddusti+k*nbolo,0,sizeof(double)*nbolo);
      memset(ddustq+k*nbolo,0,sizeof(double)*nbolo);
      memset(ddustu+k*nbolo,0,sizeof(double)*nbolo);
      memset(dpixi+k*nbolo*npixbeam,0,sizeof(double)*nbolo*npixbeam);
      memset(dpixq+k*nbolo*npixbeam,0,sizeof(double)*nbolo*npixbeam);
      memset(dpixu+k*nbolo*npixbeam,0,sizeof(double)*nbolo*npixbeam);
      

      for (l1=0;l1<ndata;l1++) {
        long ri1=htmp[l1].rg-globalBeginRing;
        
        double g1=gain[htmp[l1].adu+htmp[l1].ib*GAINSTEP];
        
        double CO1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].co
                                          -dpsisi[htmp[l1].ib]*htmp[l1].si);
        double SI1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].si
                                          +dpsisi[htmp[l1].ib]*htmp[l1].co);
        if (flg_rg[htmp[l1].ib][ri1]!=0) {
          
          htmp[l1].wp=0;
          //htmp[l1].thsig=0;
          
          SI+=htmp[l1].w*(htmp[l1].sig*g1-htmp[l1].fsl-htmp[l1].dip);
          SQ+=htmp[l1].w*(htmp[l1].sig*g1-htmp[l1].fsl-htmp[l1].dip)*CO1;
          SU+=htmp[l1].w*(htmp[l1].sig*g1-htmp[l1].fsl-htmp[l1].dip)*SI1;
          
          if (nmatco>0) {
            dcoi[htmp[l1].ib+k*nbolo] += htmp[l1].w*htmp[l1].comap;
            dcoq[htmp[l1].ib+k*nbolo] += htmp[l1].w*htmp[l1].comap*CO1;
            dcou[htmp[l1].ib+k*nbolo] += htmp[l1].w*htmp[l1].comap*SI1;
          }
          if (nmatdust>0) {
            ddusti[htmp[l1].ib+k*nbolo] += htmp[l1].w*htmp[l1].dustmap;
            ddustq[htmp[l1].ib+k*nbolo] += htmp[l1].w*htmp[l1].dustmap*CO1;
            ddustu[htmp[l1].ib+k*nbolo] += htmp[l1].w*htmp[l1].dustmap*SI1;
          }
          if (nfreefree>0) {
            dfri[htmp[l1].ib+k*nbolo] += htmp[l1].w*htmp[l1].freefree;
            dfrq[htmp[l1].ib+k*nbolo] += htmp[l1].w*htmp[l1].freefree*CO1;
            dfru[htmp[l1].ib+k*nbolo] += htmp[l1].w*htmp[l1].freefree*SI1;
          }
          
          if (ittt>0) {

            for (m=0;m<npixbeam;m++)  {
              dpixi[nbolo*m+htmp[l1].ib+k*nbolo*npixbeam] += htmp[l1].w*htmp[l1].listofpix[m];
              dpixq[nbolo*m+htmp[l1].ib+k*nbolo*npixbeam] += htmp[l1].w*htmp[l1].listofpix[m]*CO1;
              dpixu[nbolo*m+htmp[l1].ib+k*nbolo*npixbeam] += htmp[l1].w*htmp[l1].listofpix[m]*SI1;
            }
          } 
        }
      }

      solvemap(&SI,&SQ,&SU,II[k],IQ[k],IU[k],QQ[k],QU[k],UU[k]);
      SSI[k]=SI;
      SSQ[k]=SQ;
      SSU[k]=SU;
      
      for (ib=0;ib<nbolo;ib++) {
        
        if (nmatco>0) {
          solvemap(dcoi+ib+k*nbolo,dcoq+ib+k*nbolo,dcou+ib+k*nbolo,II[k],IQ[k],IU[k],QQ[k],QU[k],UU[k]);
        }
        if (nmatdust>0) {
          solvemap(ddusti+ib+k*nbolo,ddustq+ib+k*nbolo,ddustu+ib+k*nbolo,II[k],IQ[k],IU[k],QQ[k],QU[k],UU[k]);
        }
        if (nfreefree>0) {
          solvemap(dfri+ib+k*nbolo,dfrq+ib+k*nbolo,dfru+ib+k*nbolo,II[k],IQ[k],IU[k],QQ[k],QU[k],UU[k]);
        }
        
        //for (m=0;m<DEGREE;m++) {
        //  solvemap(ddegi+ib+m*nbolo+k*nbolo*DEGREE,
        //         ddegq+ib+m*nbolo+k*nbolo*DEGREE,
        //         ddegu+ib+m*nbolo+k*nbolo*DEGREE,II[k],IQ[k],IU[k],QQ[k],QU[k],UU[k]);
        //}
        
        if (ittt>0) {
          for (m=0;m<npixbeam;m++)  {
            solvemap(dpixi+ib+m*nbolo+k*nbolo*npixbeam,
                     dpixq+ib+m*nbolo+k*nbolo*npixbeam,
                     dpixu+ib+m*nbolo+k*nbolo*npixbeam,II[k],IQ[k],IU[k],QQ[k],QU[k],UU[k]);
          }
        } 
      }
      for (l1=0;l1<ndata;l1++) {
        long ri1=htmp[l1].rg-globalBeginRing;
        double CO1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].co
                                          -dpsisi[htmp[l1].ib]*htmp[l1].si);
        double SI1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].si
                                          +dpsisi[htmp[l1].ib]*htmp[l1].co);
        if (flg_rg[htmp[l1].ib][ri1]!=0) {
          double li=htmp[l1].w,lco=CO1*htmp[l1].w,lsi=SI1*htmp[l1].w;
          solvemap(&li,&lco,&lsi,II[k],IQ[k],IU[k],QQ[k],QU[k],UU[k]);
          htmp[l1].vi=li;
          htmp[l1].vq=lco;
          htmp[l1].vu=lsi;
        }
      }
      for (l1=0;l1<ndata;l1++) {
        long ri1=htmp[l1].rg-globalBeginRing;
        
        double g1=gain[htmp[l1].adu+htmp[l1].ib*GAINSTEP];
        
        double CO1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].co
                                          -dpsisi[htmp[l1].ib]*htmp[l1].si);
        double SI1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].si
                                          +dpsisi[htmp[l1].ib]*htmp[l1].co);
        if (flg_rg[htmp[l1].ib][ri1]!=0) {
          double divi=NEP_tab[htmp[l1].ib]*htmp[l1].hit*g1;
          divi=divi*divi;
          
          long l3;
          
          for (l3=0;l3<ndata;l3++) {
            long ri3=htmp[l3].rg-globalBeginRing;
            if (flg_rg[htmp[l3].ib][ri3]!=0) {
              double g3=gain[htmp[l3].adu+htmp[l3].ib*GAINSTEP];
              double sigma=NEP_tab[htmp[l3].ib]*htmp[l3].hit*(htmp[l3].vi+CO1*htmp[l3].vq+SI1*htmp[l3].vu)*g3;
              divi+=sigma*sigma;
            }
          }
          
          
          if (divi==0) htmp[l1].wp=0;
          else htmp[l1].wp=1/divi;
          
          normaoff+=NEP_tab[htmp[l1].ib]*htmp[l1].hit;
          
        }
      }
      
      for (l1=0;l1<ndata;l1++) {
        long ri1=htmp[l1].rg-globalBeginRing;
        long iri1=rgord[htmp[l1].ib][ri1]+newnr[htmp[l1].ib];
        
        double g1=gain[htmp[l1].adu+htmp[l1].ib*GAINSTEP];
        
        double CO1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].co
                                          -dpsisi[htmp[l1].ib]*htmp[l1].si);
        double SI1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].si
                                          +dpsisi[htmp[l1].ib]*htmp[l1].co);
        if (flg_rg[htmp[l1].ib][ri1]!=0) {
          double ww=htmp[l1].wp;
          double tmp=((htmp[l1].sig*g1-htmp[l1].fsl-htmp[l1].dip)-(SI+CO1*SQ+SI1*SU));
          long l3;
          
          for (l3=0;l3<ndata;l3++) {
            long ri3=htmp[l3].rg-globalBeginRing;
            if (flg_rg[htmp[l3].ib][ri3]!=0) {
              long ir3=rgord[htmp[l3].ib][ri3]+newnr[htmp[l3].ib];
              ctmp[ir3]=-(htmp[l3].vi+CO1*htmp[l3].vq+SI1*htmp[l3].vu);
            }
          }
          ctmp[iri1]+=1;
          
          for (ib=0;ib<nbolo;ib++) cco[ib]=-(dcoi[ib+k*nbolo]+CO1*dcoq[ib+k*nbolo]+SI1*dcou[ib+k*nbolo]);
          cco[htmp[l1].ib]+=htmp[l1].comap;
          
          for (ib=0;ib<nbolo;ib++) cdust[ib]=-(ddusti[ib+k*nbolo]+CO1*ddustq[ib+k*nbolo]+SI1*ddustu[ib+k*nbolo]);
          cdust[htmp[l1].ib]+=htmp[l1].dustmap;
          
          for (ib=0;ib<nbolo;ib++) ccfree[ib]=-(dfri[ib+k*nbolo]+CO1*dfrq[ib+k*nbolo]+SI1*dfru[ib+k*nbolo]);
          ccfree[htmp[l1].ib]+=htmp[l1].freefree;

          for (m=0;m<npixbeam;m++) {
            for (ib=0;ib<nbolo;ib++) cpix[ib+m*nbolo]=-(dpixi[ib+m*nbolo+k*nbolo*npixbeam]
                                                        +CO1*dpixq[ib+m*nbolo+k*nbolo*npixbeam]
                                                        +SI1*dpixu[ib+m*nbolo+k*nbolo*npixbeam]);
            cpix[htmp[l1].ib+m*nbolo]+=htmp[l1].listofpix[m];
          }
          long ir;
          /////////////////  OFFSET
          
          for (l2=0;l2<ndata;l2++) {
            long ri2=htmp[l2].rg-globalBeginRing;
            if (flg_rg[htmp[l2].ib][ri2]!=0) {
              long ir=rgord[htmp[l2].ib][ri2]+newnr[htmp[l2].ib]; 
              b2[ir]+=ww*tmp*ctmp[ir];
              hit2[ir]+=ww*ctmp[ir]*ctmp[ir];
            }
          }
          
          
          /////////////////  DIPOLE FIT
          
          for (ir=0;ir<nbolo;ir++) {
            
            ///////////// CO
            if (nmatco>0) {
              b2[newnr[nbolo]+nbolo*(npixbeam)+ir]+=ww*tmp*cco[ir];
              hit2[newnr[nbolo]+nbolo*(npixbeam)+ir]+=ww*cco[ir]*cco[ir];
            }
            
            ///////////// DUST
            
            if (nmatdust>0) {
              b2[newnr[nbolo]+nbolo*(npixbeam)+nmatco+ir]+=ww*tmp*cdust[ir];
              hit2[newnr[nbolo]+nbolo*(npixbeam)+nmatco+ir]+=ww*cdust[ir]*cdust[ir];
            }
            
            ///////////// FREEFREE
            
            if (nfreefree>0) {
              b2[newnr[nbolo]+nbolo*(npixbeam)+nmatco+nmatdust+ir]+=ww*tmp*ccfree[ir];
              hit2[newnr[nbolo]+nbolo*(npixbeam)+nmatco+nmatdust+ir]+=ww*ccfree[ir]*ccfree[ir];
            }
            
            ////////// PIXBEAM
            for (j=0;j<npixbeam;j++) {
              b2[newnr[nbolo]+j*nbolo+ir]+=ww*tmp*cpix[ir+j*nbolo];
              hit2[newnr[nbolo]+j*nbolo+ir]+=ww*cpix[ir+j*nbolo]*cpix[ir+j*nbolo];
            }
          }
        }
      }
    }
  }


  if (rank==0) {
    double *lb = (double *) malloc(sizeof(double)*(nmatres));
    for (rrk=1;rrk<mpi_size;rrk++) {
      MPI_Recv(lb,sizeof(double)*(nmatres), MPI_BYTE, rrk,1031, MPI_COMM_WORLD,&statu);
      for (l=0;l<nmatres;l++) b2[l]+=lb[l];
    }
    free(lb);
  }
  else MPI_Send(b2, sizeof(double)*(nmatres), MPI_BYTE, 0, 1031, MPI_COMM_WORLD);
  
  if (rank==0) {
    double *lb = (double *) malloc(sizeof(double)*(nmatres));
    for (rrk=1;rrk<mpi_size;rrk++) {
      MPI_Recv(lb,sizeof(double)*(nmatres), MPI_BYTE, rrk,1033, MPI_COMM_WORLD,&statu);

      for (l=0;l<nmatres;l++) hit2[l]+=lb[l];
    }
    free(lb);
  }
  else {
    MPI_Send(hit2, sizeof(double)*(nmatres), MPI_BYTE, 0, 1033, MPI_COMM_WORLD);
  }
  

  //==========================================================
  // Compute Ax
  //healspline_fill_s_with_PtP_s (hs_rec, Ndata, vec, x, q);
  // + 
  // Preconditionnement
  //
  //
  
  for (l=0;l<nmatres;l++) q2[l]=0;
  if (rank==0) {

    double soff=0;
    for (i=0;i<newnr[nbolo];i++) soff+=hit2[0]*ix2[i];
    for (i=0;i<newnr[nbolo];i++) q2[i]=soff;
    if (Param->flag_AVGR0==(_PAR_TRUE)) {
      soff=0;
      for (i=0;i<nbolo;i++) soff+=hit2[newnr[nbolo]]*
                              ix2[newnr[nbolo]+i];
      for (i=0;i<nbolo;i++) q2[newnr[nbolo]+i]=soff; 
      for (i=0;i<nbolo;i++) b2[newnr[nbolo]+i]+=
                              hit2[newnr[nbolo]]*Param->AVGR0;
    }
    if (nmatco>0) {
      soff=0;
      for (i=0;i<nbolo;i++) soff+=hit2[newnr[nbolo]+nbolo*(npixbeam)]
                              *ix2[newnr[nbolo]+nbolo*(npixbeam)+i];
      for (i=0;i<nbolo;i++) q2[newnr[nbolo]+nbolo*(npixbeam)+i]=soff;
      for (i=0;i<nbolo;i++) b2[newnr[nbolo]+nbolo*(npixbeam)+i]+=
                              hit2[newnr[nbolo]+nbolo*(npixbeam)]*Param->AVGCO;
    }
    if (nmatdust>0) {
      soff=0;
      for (i=0;i<nbolo;i++) 
        soff+=hit2[newnr[nbolo]+nbolo*(npixbeam)+nmatco]*
          ix2[newnr[nbolo]+nbolo*(npixbeam)+nmatco+i];
      for (i=0;i<nbolo;i++) q2[newnr[nbolo]+nbolo*(npixbeam)+nmatco+i]=soff;
      for (i=0;i<nbolo;i++) b2[newnr[nbolo]+nbolo*(npixbeam)+nmatco+i]+=
                              hit2[newnr[nbolo]+nbolo*(npixbeam)+nmatco]*Param->AVGDUST;
    }
    if (nfreefree>0) {
      soff=0;
      for (i=0;i<nbolo;i++) 
        soff+=hit2[newnr[nbolo]+nbolo*(npixbeam)+nmatco+nmatdust]*
          ix2[newnr[nbolo]+nbolo*(npixbeam)+nmatco+nmatdust+i];
      for (i=0;i<nbolo;i++) q2[newnr[nbolo]+nbolo*(npixbeam)+nmatco+nmatdust+i]=soff;
    }
    
  }

  //ptest=0;
  for (k=0;k<nnbpix;k++) if (flgpix[k]>0) {
    
    long ndata = loc_nhpix[k];
    if (ndata>0) {
#ifdef TIMING
      gettimeofday(&tp1,NULL);
#endif
      hpix *htmp = loc_hpix[k];

      double vali=0,valq=0,valu=0;

      long l3;
      for (l3=0;l3<ndata;l3++) {
        long ri3=htmp[l3].rg-globalBeginRing;
        if (flg_rg[htmp[l3].ib][ri3]!=0) {
          long ir3=rgord[htmp[l3].ib][ri3]+newnr[htmp[l3].ib];
          vali+=htmp[l3].vi*ix2[ir3];
          valq+=htmp[l3].vq*ix2[ir3];
          valu+=htmp[l3].vu*ix2[ir3];
        }
      }
      
      for (ib=0;ib<nbolo;ib++) {
        if (nmatco>0) {
          vali+=dcoi[ib+k*nbolo]  *ix2[newnr[nbolo]+nbolo*(npixbeam)+ib];
          valq+=dcoq[ib+k*nbolo]  *ix2[newnr[nbolo]+nbolo*(npixbeam)+ib];
          valu+=dcou[ib+k*nbolo]  *ix2[newnr[nbolo]+nbolo*(npixbeam)+ib];
        }           
        if (nmatdust>0) {
          vali+=ddusti[ib+k*nbolo]  *ix2[newnr[nbolo]+nbolo*(npixbeam)+nmatco+ib];
          valq+=ddustq[ib+k*nbolo]  *ix2[newnr[nbolo]+nbolo*(npixbeam)+nmatco+ib];
          valu+=ddustu[ib+k*nbolo]  *ix2[newnr[nbolo]+nbolo*(npixbeam)+nmatco+ib];
        }           
        if (nfreefree>0) {
          vali+=dfri[ib+k*nbolo]  *ix2[newnr[nbolo]+nbolo*(npixbeam)+nmatco+nmatdust+ib];
          valq+=dfrq[ib+k*nbolo]  *ix2[newnr[nbolo]+nbolo*(npixbeam)+nmatco+nmatdust+ib];
          valu+=dfru[ib+k*nbolo]  *ix2[newnr[nbolo]+nbolo*(npixbeam)+nmatco+nmatdust+ib];
        }           
        
        for (m=0;m<npixbeam;m++) {
          vali+=ix2[newnr[nbolo]+ib+m*nbolo]*dpixi[ib+m*nbolo+k*nbolo*npixbeam];
          valq+=ix2[newnr[nbolo]+ib+m*nbolo]*dpixq[ib+m*nbolo+k*nbolo*npixbeam];
          valu+=ix2[newnr[nbolo]+ib+m*nbolo]*dpixu[ib+m*nbolo+k*nbolo*npixbeam];
        }  
        
      }
      double qri=0;
      double qrq=0;
      double qru=0;
      
      for (l1=0;l1<ndata;l1++) {
        long ri1=htmp[l1].rg-globalBeginRing;
        long iri1=rgord[htmp[l1].ib][ri1]+newnr[htmp[l1].ib];
        
        double CO1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].co
                                          -dpsisi[htmp[l1].ib]*htmp[l1].si);
        double SI1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].si
                                          +dpsisi[htmp[l1].ib]*htmp[l1].co);
        if (flg_rg[htmp[l1].ib][ri1]!=0) {
          double ww=htmp[l1].wp;
          double val2=ix2[iri1]-(vali+CO1*valq+SI1*valu);
          
          if (nmatco>0)   
            val2+=ix2[newnr[nbolo]+nbolo*(npixbeam)+htmp[l1].ib]*htmp[l1].comap;
          if (nmatdust>0) 
            val2+=ix2[newnr[nbolo]+nbolo*(npixbeam)+nmatco+htmp[l1].ib]*htmp[l1].dustmap;
          if (nfreefree>0) 
            val2+=ix2[newnr[nbolo]+nbolo*(npixbeam)+nmatco+nmatdust+htmp[l1].ib]*htmp[l1].freefree;
          
          for (m=0;m<npixbeam;m++) {
            val2+=ix2[newnr[nbolo]+htmp[l1].ib+m*nbolo]*htmp[l1].listofpix[m];
          }
          
          
          qri-=ww*val2;
          qrq-=ww*val2*CO1;
          qru-=ww*val2*SI1;
          
          
          q2[iri1]+=ww*val2;
          
          ///////////// CO
          if (nmatco>0) {
            q2[newnr[nbolo]+nbolo*(npixbeam)+htmp[l1].ib]+=ww*val2*htmp[l1].comap;
          }
          
          ///////////// DUST
          if (nmatdust>0) {
            q2[newnr[nbolo]+nbolo*(npixbeam)+nmatco+htmp[l1].ib]+=ww*val2*htmp[l1].dustmap;
          }
          
          ///////////// FREEFREE
          if (nfreefree>0) {
            q2[newnr[nbolo]+nbolo*(npixbeam)+nmatco+nmatdust+htmp[l1].ib]+=ww*val2*htmp[l1].freefree;
          }
          
          ////////// SYSTE
          for (j=0;j<npixbeam;j++) {
            q2[newnr[nbolo]+j*nbolo+htmp[l1].ib]+=ww*val2*htmp[l1].listofpix[j];  
          }
          
          
        }
      }
    

      for (l2=0;l2<ndata;l2++) {
        long ri2=htmp[l2].rg-globalBeginRing;
        if (flg_rg[htmp[l2].ib][ri2]!=0) {
          long ir=rgord[htmp[l2].ib][ri2]+newnr[htmp[l2].ib]; 
          q2[ir]+=qri*htmp[l2].vi+qrq*htmp[l2].vq+qru*htmp[l2].vu;
        }
      }
    
        
      for (ib=0;ib<nbolo;ib++) {
        
        ///////////// CO
        if (nmatco>0) {
          q2[newnr[nbolo]+nbolo*(npixbeam)+ib]+=qri*dcoi[ib+k*nbolo]+qrq*dcoq[ib+k*nbolo]+qru*dcou[ib+k*nbolo];
        }
        ///////////// DUST
        if (nmatdust>0) {
          q2[newnr[nbolo]+nbolo*(npixbeam)+nmatco+ib]+=qri*ddusti[ib+k*nbolo]+qrq*ddustq[ib+k*nbolo]+qru*ddustu[ib+k*nbolo];
        }
        ///////////// FREEFREE
        if (nfreefree>0) {
          q2[newnr[nbolo]+nbolo*(npixbeam)+nmatco+nmatdust+ib]+=qri*dfri[ib+k*nbolo]+qrq*dfrq[ib+k*nbolo]+qru*dfru[ib+k*nbolo];
        }
        
        for (j=0;j<npixbeam;j++) {
          q2[newnr[nbolo]+j*nbolo+ib]+=qri*dpixi[ib+j*nbolo+k*nbolo*npixbeam]+qrq*dpixq[ib+j*nbolo+k*nbolo*npixbeam]+qru*dpixu[ib+j*nbolo+k*nbolo*npixbeam];  
        }
      }
    }
      
#ifdef TIMING
      gettimeofday(&tp2,NULL);
      double dt=(double)(tp2.tv_sec-tp1.tv_sec)+(1E-6)*(tp2.tv_usec-tp1.tv_usec);
      dthit[ndata]+=dt;
      ndthit[ndata]+=1;
#endif
  }

  if (rank==0) {
    double *lb = (double *) malloc(sizeof(double)*(nmatres));
    for (rrk=1;rrk<mpi_size;rrk++) {
      MPI_Recv(lb,sizeof(double)*(nmatres), MPI_BYTE, rrk,1032, MPI_COMM_WORLD,&statu);
      for (l=0;l<nmatres;l++) q2[l]+=lb[l];
    }
    free(lb);
  }
  else {
    MPI_Send(q2, sizeof(double)*(nmatres), MPI_BYTE, 0, 1032, MPI_COMM_WORLD);
  }

  if (rank==0) {
    for (i=0; i < nmatres; i++)
      {
        r2[i] = b2[i] - q2[i];
        d2[i] = r2[i] / hit2[i];
      }
  }

  
  double delta_new_tmp = 0.0;
  if (rank==0) {
    for (i = 0; i < nmatres; i++) {
      delta_new_tmp += b2[i] ;
      if (isnan(b2[i])) {
        fprintf(stderr,"NAN B2 PBS %ld\n",(long) i);
      }
      if (isnan(d2[i])) {
        fprintf(stderr,"NAN D2 PBS %ld %lg\n",(long) i,hit2[i]);
      }
    }
    //fprintf(stderr,"B2 %lg\n",delta_new_tmp);
  
    delta_new_tmp = 0.0;
    for (i = 0; i < nmatres; i++) {
      delta_new_tmp += q2[i] ;
      if (isnan(q2[i])) {
        fprintf(stderr,"NAN Q2 PBS %ld\n",(long) i);
      }
    }
    //fprintf(stderr,"Q2 %lg\n",delta_new_tmp);
    delta_new_tmp = 0.0;
    for (i = 0; i < nmatres; i++) {
      delta_new_tmp += q2[i]-b2[i] ;
      //fprintf(stderr,"B2 Q2 B2-Q2 [%ld]: %lg\t%lg\t%lg\n",(long) i,b2[i],q2[i],q2[i]-b2[i]);
    }
  }
  
  delta_new_tmp = 0.0;
  if (rank==0) for (i=0; i < nmatres; i++) {
    delta_new_tmp += r2[i] * d2[i];
  }
  
  delta_new=0;
  for (i=0;i<mpi_size;i++) {
    double tmp=delta_new_tmp;
    MPI_Bcast(&tmp, sizeof(double), MPI_BYTE, i, MPI_COMM_WORLD);
    delta_new+=tmp;
  }
  if (itbogo==0) delta0 = delta_new;
  if (rank==0) fprintf (stderr, "min() iter = %d - delta0 = %lg - delta_new = %lg\n", iter, delta0, delta_new);
  //int testwrit=0;

  MPI_Bcast(&normaoff, sizeof(double), MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&delta0, sizeof(double), MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&delta_new, sizeof(double), MPI_BYTE, 0, MPI_COMM_WORLD);


  while ((iter < itermax)  && ((delta_new) > delta0*1E-24) && ((delta_new) > 1E-20)) //Param->XI2STOP))
    {
      // q <= Ad
      if (rank==0) {
        memcpy(x2old,ix2,nmatres*sizeof(double));
        //testwrit=1;
     }
     //else {
     //  //testwrit=0;
     //}

      //healspline_fill_s_with_PtP_s (hs_rec, Ndata, vec, d, q);
      //for (i=0;i<mpi_size;i++) {
      //MPI_Bcast(d+tab_begr[i]*2, sizeof(double)*(tab_edr[i]-tab_begr[i]+1)*2, MPI_BYTE, i, MPI_COMM_WORLD);
      //}
      MPI_Bcast(d2, sizeof(double)*(nmatres), MPI_BYTE, 0, MPI_COMM_WORLD);
      
      // ===========================================================================
      // PROJECTION DE d dans q
      //
      
      for (l=0;l<nmatres;l++) q2[l]=0;
      if (rank==0) {
        double soff=0;
        for (i=0;i<newnr[nbolo];i++) soff+=hit2[0]*d2[i];
        for (i=0;i<newnr[nbolo];i++) q2[i]=soff;
        if (Param->flag_AVGR0==_PAR_TRUE) {
          soff=0;
          for (i=0;i<nbolo;i++) soff+=hit2[newnr[nbolo]]*
                                  d2[newnr[nbolo]+i];
          for (i=0;i<nbolo;i++) q2[newnr[nbolo]+i]=soff; 
        }
        if (nmatco>0) {
          soff=0;
          for (i=0;i<nbolo;i++) soff+=hit2[newnr[nbolo]+nbolo*(npixbeam)]
                                  *d2[newnr[nbolo]+nbolo*(npixbeam)+i];
          for (i=0;i<nbolo;i++) q2[newnr[nbolo]+nbolo*(npixbeam)+i]=soff;
        }
        if (nfreefree>0) {
          soff=0;
          for (i=0;i<nbolo;i++) soff+=hit2[newnr[nbolo]+nbolo*(npixbeam)+nmatco+nmatdust]
                                  *d2[newnr[nbolo]+nbolo*(npixbeam)+nmatco+nmatdust+i];
          for (i=0;i<nbolo;i++) q2[newnr[nbolo]+nbolo*(npixbeam)+i+nmatco+nmatdust]=soff;
        }
        if (nmatdust>0) {
          soff=0;
          for (i=0;i<nbolo;i++) 
            soff+=hit2[newnr[nbolo]+nbolo*(npixbeam)+nmatco]*
              d2[newnr[nbolo]+nbolo*(npixbeam)+nmatco+i];
          for (i=0;i<nbolo;i++) q2[newnr[nbolo]+nbolo*(npixbeam)+nmatco+i]=soff;
        }
      }
      for (k=0;k<nnbpix;k++) if (flgpix[k]>0) {
        
        //long imat=the_stat_pix[k];
        long ndata = loc_nhpix[k];
        if (ndata>0) {
          hpix *htmp = loc_hpix[k];

          double vali=0,valq=0,valu=0;

          long l3;
          for (l3=0;l3<ndata;l3++) {
            long ri3=htmp[l3].rg-globalBeginRing;
            if (flg_rg[htmp[l3].ib][ri3]!=0) {
              long ir3=rgord[htmp[l3].ib][ri3]+newnr[htmp[l3].ib];
              vali+=htmp[l3].vi*d2[ir3];
              valq+=htmp[l3].vq*d2[ir3];
              valu+=htmp[l3].vu*d2[ir3];
            }
          }
          for (ib=0;ib<nbolo;ib++) {
            if (nmatco>0) {
              vali+=dcoi[ib+k*nbolo]  *d2[newnr[nbolo]+nbolo*(npixbeam)+ib];
              valq+=dcoq[ib+k*nbolo]  *d2[newnr[nbolo]+nbolo*(npixbeam)+ib];
              valu+=dcou[ib+k*nbolo]  *d2[newnr[nbolo]+nbolo*(npixbeam)+ib];
            }       
            if (nmatdust>0) {
              vali+=ddusti[ib+k*nbolo]  *d2[newnr[nbolo]+nbolo*(npixbeam)+nmatco+ib];
              valq+=ddustq[ib+k*nbolo]  *d2[newnr[nbolo]+nbolo*(npixbeam)+nmatco+ib];
              valu+=ddustu[ib+k*nbolo]  *d2[newnr[nbolo]+nbolo*(npixbeam)+nmatco+ib];
            }       
            if (nfreefree>0) {
              vali+=dfri[ib+k*nbolo]  *d2[newnr[nbolo]+nbolo*(npixbeam)+nmatco+nmatdust+ib];
              valq+=dfrq[ib+k*nbolo]  *d2[newnr[nbolo]+nbolo*(npixbeam)+nmatco+nmatdust+ib];
              valu+=dfru[ib+k*nbolo]  *d2[newnr[nbolo]+nbolo*(npixbeam)+nmatco+nmatdust+ib];
            }       
            for (m=0;m<npixbeam;m++) {
              vali+=d2[newnr[nbolo]+ib+m*nbolo]*dpixi[ib+m*nbolo+k*nbolo*npixbeam];
              valq+=d2[newnr[nbolo]+ib+m*nbolo]*dpixq[ib+m*nbolo+k*nbolo*npixbeam];
              valu+=d2[newnr[nbolo]+ib+m*nbolo]*dpixu[ib+m*nbolo+k*nbolo*npixbeam];
            }   
          }

          double qri=0;
          double qrq=0;
          double qru=0;

          for (l1=0;l1<ndata;l1++) {
            long ri1=htmp[l1].rg-globalBeginRing;
            long iri1=rgord[htmp[l1].ib][ri1]+newnr[htmp[l1].ib];

            double CO1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].co
                                              -dpsisi[htmp[l1].ib]*htmp[l1].si);
            double SI1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].si
                                              +dpsisi[htmp[l1].ib]*htmp[l1].co);
            if (flg_rg[htmp[l1].ib][ri1]!=0) {
              double ww=htmp[l1].wp;
              double val2=d2[iri1]-(vali+CO1*valq+SI1*valu);

              if (nmatco>0)   
                val2+=d2[newnr[nbolo]+nbolo*(npixbeam)+htmp[l1].ib]*htmp[l1].comap;
              if (nmatdust>0) 
                val2+=d2[newnr[nbolo]+nbolo*(npixbeam)+nmatco+htmp[l1].ib]*htmp[l1].dustmap;
              if (nfreefree>0) 
                val2+=d2[newnr[nbolo]+nbolo*(npixbeam)+nmatco+nmatdust+htmp[l1].ib]*htmp[l1].freefree;
                
              
              for (m=0;m<npixbeam;m++) {
                val2+=d2[newnr[nbolo]+htmp[l1].ib+m*nbolo]*htmp[l1].listofpix[m];
              }
              
                
              qri-=ww*val2;
              qrq-=ww*val2*CO1;
              qru-=ww*val2*SI1;
                
              q2[iri1]+=ww*val2;

              ///////////// CO
              if (nmatco>0) {
                q2[newnr[nbolo]+nbolo*(npixbeam)+htmp[l1].ib]+=ww*val2*htmp[l1].comap;
              }
                
              ///////////// DUST
              if (nmatdust>0) {
                q2[newnr[nbolo]+nbolo*(npixbeam)+nmatco+htmp[l1].ib]+=ww*val2*htmp[l1].dustmap;
              }
                
              ///////////// FREEFREE
              if (nfreefree>0) {
                q2[newnr[nbolo]+nbolo*(npixbeam)+nmatco+nmatdust+htmp[l1].ib]+=ww*val2*htmp[l1].freefree;
              }
            
              for (j=0;j<npixbeam;j++) {
                q2[newnr[nbolo]+j*nbolo+htmp[l1].ib]+=ww*val2*htmp[l1].listofpix[j];  
              }
                
            }
          }


          for (l2=0;l2<ndata;l2++) {
            long ri2=htmp[l2].rg-globalBeginRing;
            if (flg_rg[htmp[l2].ib][ri2]!=0) {
              long ir=rgord[htmp[l2].ib][ri2]+newnr[htmp[l2].ib]; 
              q2[ir]+=qri*htmp[l2].vi+qrq*htmp[l2].vq+qru*htmp[l2].vu;
            }
          }
            
          for (ib=0;ib<nbolo;ib++) {

            ///////////// CO
            if (nmatco>0) {
              q2[newnr[nbolo]+nbolo*(npixbeam)+ib]+=qri*dcoi[ib+k*nbolo]+qrq*dcoq[ib+k*nbolo]+qru*dcou[ib+k*nbolo];
            }
            ///////////// DUST
            if (nmatdust>0) {
              q2[newnr[nbolo]+nbolo*(npixbeam)+nmatco+ib]+=qri*ddusti[ib+k*nbolo]+qrq*ddustq[ib+k*nbolo]+qru*ddustu[ib+k*nbolo];
            }
            ///////////// FREEFREE
            if (nfreefree>0) {
              q2[newnr[nbolo]+nbolo*(npixbeam)+nmatco+nmatdust+ib]+=qri*dfri[ib+k*nbolo]+qrq*dfrq[ib+k*nbolo]+qru*dfru[ib+k*nbolo];
            }

            for (j=0;j<npixbeam;j++) {
              q2[newnr[nbolo]+j*nbolo+ib]+=qri*dpixi[ib+j*nbolo+k*nbolo*npixbeam]+qrq*dpixq[ib+j*nbolo+k*nbolo*npixbeam]+qru*dpixu[ib+j*nbolo+k*nbolo*npixbeam];  
            }
          }
        }
      }
      if (rank==0) {
        double *lb = (double *) malloc(sizeof(double)*(nmatres));
        for (rrk=1;rrk<mpi_size;rrk++) {
          MPI_Recv(lb,sizeof(double)*(nmatres), MPI_BYTE, rrk,1034, MPI_COMM_WORLD,&statu);
          for (l=0;l<nmatres;l++) q2[l]+=lb[l];
        }
        free(lb);
      }
      else {
        MPI_Send(q2, sizeof(double)*(nmatres), MPI_BYTE, 0, 1034, MPI_COMM_WORLD);
      }
      double dtq = 0.0;
      if (rank==0) {
        for (i=0; i < nmatres; i++) dtq += d2[i] * q2[i];
        alpha = delta_new / dtq;
        for (i=0; i < nmatres ; i++) ix2[i] += alpha * d2[i];
      }


      gettimeofday(&tp2,NULL);
      //if (rank==0) fprintf (stderr,"iter = %d - delta0 = %lg - delta_new = %lg %ld %3lfs\n", iter, delta0, delta_new,
      //                    (long) testwrit,(double)(tp2.tv_sec-tp1.tv_sec)+(1E-6)*(tp2.tv_usec-tp1.tv_usec));

      if (rank==0) fprintf (stderr,".");
      gettimeofday(&tp1,NULL);
          
      if (iter % 100 == 0 && iter !=0)
        {
          // Use the best case
          memcpy(ix2,x2old,nmatres*sizeof(double));
          //for (i=0;i<mpi_size;i++) {
          //  MPI_Bcast(x+tab_begr[i]*2, sizeof(double)*(tab_edr[i]-tab_begr[i]+1)*2, MPI_BYTE, i, MPI_COMM_WORLD);
          //}
          MPI_Bcast(ix2, sizeof(double)*(nmatres), MPI_BYTE, 0, MPI_COMM_WORLD);

          
      // ===========================================================================
      // PROJECTION DE x dans q
      //
          
          for (l=0;l<nmatres;l++) q2[l]=0;
          if (rank==0) {
            double soff=0;
            for (i=0;i<newnr[nbolo];i++) soff+=hit2[0]*ix2[i];
            for (i=0;i<newnr[nbolo];i++) q2[i]=soff;
            if (Param->flag_AVGR0==_PAR_TRUE) {
              soff=0;
              for (i=0;i<nbolo;i++) soff+=hit2[newnr[nbolo]]*
                                      ix2[newnr[nbolo]+i];
              for (i=0;i<nbolo;i++) q2[newnr[nbolo]+i]=soff; 
            }
            if (nmatco>0) {
              soff=0;
              for (i=0;i<nbolo;i++) soff+=hit2[newnr[nbolo]+nbolo*(npixbeam)]
                                      *ix2[newnr[nbolo]+nbolo*(npixbeam)+i];
              for (i=0;i<nbolo;i++) q2[newnr[nbolo]+nbolo*(npixbeam)+i]=soff;
            }
            if (nmatdust>0) {
              soff=0;
              for (i=0;i<nbolo;i++) 
                soff+=hit2[newnr[nbolo]+nbolo*(npixbeam)+nmatco]*
                  ix2[newnr[nbolo]+nbolo*(npixbeam)+nmatco+i];
              for (i=0;i<nbolo;i++) q2[newnr[nbolo]+nbolo*(npixbeam)+nmatco+i]=soff;
            }
            if (nfreefree>0) {
              soff=0;
              for (i=0;i<nbolo;i++) 
                soff+=hit2[newnr[nbolo]+nbolo*(npixbeam)+nmatco+nmatdust]*
                  ix2[newnr[nbolo]+nbolo*(npixbeam)+nmatco+nmatdust+i];
              for (i=0;i<nbolo;i++) q2[newnr[nbolo]+nbolo*(npixbeam)+nmatco+nmatdust+i]=soff;
            }
            
            
          }

          for (k=0;k<nnbpix;k++) if (flgpix[k]>0) {
        
            //long imat=the_stat_pix[k];
            long ndata = loc_nhpix[k];
            if (ndata>0) {
              hpix *htmp = loc_hpix[k];

              double vali=0,valq=0,valu=0;

              long l3;
              for (l3=0;l3<ndata;l3++) {
                long ri3=htmp[l3].rg-globalBeginRing;
                if (flg_rg[htmp[l3].ib][ri3]!=0) {
                  long ir3=rgord[htmp[l3].ib][ri3]+newnr[htmp[l3].ib];
                  vali+=htmp[l3].vi*ix2[ir3];
                  valq+=htmp[l3].vq*ix2[ir3];
                  valu+=htmp[l3].vu*ix2[ir3];
                }
              }
                
              for (ib=0;ib<nbolo;ib++) {

                if (nmatco>0) {
                  vali+=dcoi[ib+k*nbolo]  *ix2[newnr[nbolo]+nbolo*(npixbeam)+ib];
                  valq+=dcoq[ib+k*nbolo]  *ix2[newnr[nbolo]+nbolo*(npixbeam)+ib];
                  valu+=dcou[ib+k*nbolo]  *ix2[newnr[nbolo]+nbolo*(npixbeam)+ib];
                }           
                if (nmatdust>0) {
                  vali+=ddusti[ib+k*nbolo]  *ix2[newnr[nbolo]+nbolo*(npixbeam)+nmatco+ib];
                  valq+=ddustq[ib+k*nbolo]  *ix2[newnr[nbolo]+nbolo*(npixbeam)+nmatco+ib];
                  valu+=ddustu[ib+k*nbolo]  *ix2[newnr[nbolo]+nbolo*(npixbeam)+nmatco+ib];
                }           
                if (nfreefree>0) {
                  vali+=dfri[ib+k*nbolo]  *ix2[newnr[nbolo]+nbolo*(npixbeam)+nmatco+nmatdust+ib];
                  valq+=dfrq[ib+k*nbolo]  *ix2[newnr[nbolo]+nbolo*(npixbeam)+nmatco+nmatdust+ib];
                  valu+=dfru[ib+k*nbolo]  *ix2[newnr[nbolo]+nbolo*(npixbeam)+nmatco+nmatdust+ib];
                }           
                for (m=0;m<npixbeam;m++) {
                  vali+=ix2[newnr[nbolo]+ib+m*nbolo]*dpixi[ib+m*nbolo+k*nbolo*npixbeam];
                  valq+=ix2[newnr[nbolo]+ib+m*nbolo]*dpixq[ib+m*nbolo+k*nbolo*npixbeam];
                  valu+=ix2[newnr[nbolo]+ib+m*nbolo]*dpixu[ib+m*nbolo+k*nbolo*npixbeam];
                }   
              }

              double qri=0;
              double qrq=0;
              double qru=0;

              for (l1=0;l1<ndata;l1++) {
                long ri1=htmp[l1].rg-globalBeginRing;
                long iri1=rgord[htmp[l1].ib][ri1]+newnr[htmp[l1].ib];

                double CO1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].co
                                                  -dpsisi[htmp[l1].ib]*htmp[l1].si);
                double SI1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].si
                                                  +dpsisi[htmp[l1].ib]*htmp[l1].co);
                if (flg_rg[htmp[l1].ib][ri1]!=0) {
                  double ww=htmp[l1].wp;
                  double val2=ix2[iri1]-(vali+CO1*valq+SI1*valu);

                  if (nmatco>0)   
                    val2+=ix2[newnr[nbolo]+nbolo*(npixbeam)+htmp[l1].ib]*htmp[l1].comap;
                  if (nmatdust>0) 
                    val2+=ix2[newnr[nbolo]+nbolo*(npixbeam)+nmatco+htmp[l1].ib]*htmp[l1].dustmap;
                  if (nfreefree>0) 
                    val2+=ix2[newnr[nbolo]+nbolo*(npixbeam)+nmatco+nfreefree+htmp[l1].ib]*htmp[l1].freefree;

                  for (m=0;m<npixbeam;m++) {
                    val2+=ix2[newnr[nbolo]+htmp[l1].ib+m*nbolo]*htmp[l1].listofpix[m];
                  }
                    
                  qri-=ww*val2;
                  qrq-=ww*val2*CO1;
                  qru-=ww*val2*SI1;
                    
                  q2[iri1]+=ww*val2;

                  ///////////// CO
                  if (nmatco>0) {
                    q2[newnr[nbolo]+nbolo*(npixbeam)+htmp[l1].ib]+=ww*val2*htmp[l1].comap;
                  }
                    
                  ///////////// DUST
                  if (nmatdust>0) {
                    q2[newnr[nbolo]+nbolo*(npixbeam)+nmatco+htmp[l1].ib]+=ww*val2*htmp[l1].dustmap;
                  }

                  ///////////// FREEFREE
                  if (nfreefree>0) {
                    q2[newnr[nbolo]+nbolo*(npixbeam)+nmatco+nfreefree+htmp[l1].ib]+=ww*val2*htmp[l1].freefree;
                  }
                    
                  for (j=0;j<npixbeam;j++) {
                    q2[newnr[nbolo]+j*nbolo+htmp[l1].ib]+=ww*val2*htmp[l1].listofpix[j];  
                  }
                    
                }
              }

              for (l2=0;l2<ndata;l2++) {
                long ri2=htmp[l2].rg-globalBeginRing;
                if (flg_rg[htmp[l2].ib][ri2]!=0) {
                  long ir=rgord[htmp[l2].ib][ri2]+newnr[htmp[l2].ib]; 
                  q2[ir]+=qri*htmp[l2].vi+qrq*htmp[l2].vq+qru*htmp[l2].vu;
                }
              }
                
              for (ib=0;ib<nbolo;ib++) {

                ///////////// CO
                if (nmatco>0) {
                  q2[newnr[nbolo]+nbolo*(npixbeam)+ib]+=qri*dcoi[ib+k*nbolo]+qrq*dcoq[ib+k*nbolo]+qru*dcou[ib+k*nbolo];
                }
                ///////////// DUST
                if (nmatdust>0) {
                  q2[newnr[nbolo]+nbolo*(npixbeam)+nmatco+ib]+=qri*ddusti[ib+k*nbolo]+qrq*ddustq[ib+k*nbolo]+qru*ddustu[ib+k*nbolo];
                }

                ///////////// FREEFREE
                if (nfreefree>0) {
                  q2[newnr[nbolo]+nbolo*(npixbeam)+nmatco+nmatdust+ib]+=qri*dfri[ib+k*nbolo]+qrq*dfrq[ib+k*nbolo]+qru*dfru[ib+k*nbolo];
                }
                  
                for (j=0;j<npixbeam;j++) {
                  q2[newnr[nbolo]+j*nbolo+ib]+=qri*dpixi[ib+j*nbolo+k*nbolo*npixbeam]+qrq*dpixq[ib+j*nbolo+k*nbolo*npixbeam]+qru*dpixu[ib+j*nbolo+k*nbolo*npixbeam];  
                }
              }
            }
          }
          if (rank==0) {
            double *lb = (double *) malloc(sizeof(double)*(nmatres));
            for (rrk=1;rrk<mpi_size;rrk++) {
              MPI_Recv(lb,sizeof(double)*(nmatres), MPI_BYTE, rrk,1034, MPI_COMM_WORLD,&statu);
              for (l=0;l<nmatres;l++) q2[l]+=lb[l];
            }

            free(lb);
          }
          else {
            MPI_Send(q2, sizeof(double)*(nmatres), MPI_BYTE, 0, 1034, MPI_COMM_WORLD);
          }
          
          if (rank==0) for (i=0; i < nmatres; i++) {
            r2[i] = b2[i] - q2[i];
          }
          
        }
      else
        {
          if (rank==0) for (i=0; i < nmatres; i++) r2[i] -= alpha * q2[i];
        }
      
      if (rank==0) for (i=0; i < nmatres; i++) s2[i] = r2[i] / hit2[i];

      delta_old = delta_new;
      if (rank==0) {
        delta_new=0;
        for (i=0; i < nmatres ; i++) delta_new += r2[i] * s2[i];
        beta = delta_new / delta_old;
        for (i=0; i < nmatres ; i++) d2[i] = s2[i] + beta * d2[i];
      }
      iter ++;
      MPI_Bcast(&delta_new, sizeof(double), MPI_BYTE, 0, MPI_COMM_WORLD);
      
    }
  
  if (rank==0) fprintf (stderr,"\nmin2() iter = %d - delta0 = %lg - delta_new = %lg\n", 
                        iter, delta0, delta_new);
  if (rank==0) fprintf (stderr,"CG in iter = %d (max=%d)\n", iter, itermax);
  
  if (rank==0) memcpy(ix2  ,x2old,nmatres*sizeof (double));
  MPI_Bcast(ix2, sizeof(double)*(nmatres), MPI_BYTE, 0, MPI_COMM_WORLD);
  itbogo++;
#if 0
  if (mintest!=0) {
    if (n_cache_xi2==0) {
      cache_xi2 = (double **) malloc(sizeof(double *));
      cache_gain_xi2 = (double **) malloc(sizeof(double *));
    }
    else {
      cache_xi2 = (double **) realloc(cache_xi2,sizeof(double *)*(n_cache_xi2+1));
      cache_gain_xi2 = (double **) realloc(cache_gain_xi2,sizeof(double *)*(n_cache_xi2+1));
    }
    cache_xi2[n_cache_xi2]= (double *) malloc(sizeof(double)*(newnr[nbolo]+nbolo*npixbeam+nmatco+nmatdust+nfreefree));
    cache_gain_xi2[n_cache_xi2]= (double *) malloc(sizeof(double)*(nbolo*GAINSTEP));

    memcpy(cache_xi2[n_cache_xi2],ix2,sizeof(double)*(newnr[nbolo]+nbolo*npixbeam+nmatco+nmatdust+nfreefree));
    memcpy(cache_gain_xi2[n_cache_xi2],gain,sizeof(double)*(nbolo*GAINSTEP));
    n_cache_xi2++;
  }
#endif
}


void minimize_gain_gi(double *ix2,double *gain,double *gaingi)
{
  MPI_Status statu;
  long i,rrk,j,k,l1,l2,ib;
  int itermax = 500;
  int iter;
  double  delta_new, delta_old, beta;
  double  alpha=1.0; // get rid of gcc "maybe-uninitialized" warning depending on optimsation level

  PIOLONG GAINSTEP2;
  int rank;
  int size;
  int mpi_size;
  int rank_size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank_size);
  rank=rank_size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  mpi_size=size;

  GAINSTEP2=GAINSTEP;

  nmatres=newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+nfreefree;

  if (rank==0) {
    fprintf(stderr,"\nmatres=%ld, nnbpix=%ld\n\n", nmatres, nnbpix);

    fprintf(stderr,"==============================\n\nminimize_gain_gi(): ITERATION %ld\n\n==============================\n",itbogo);
    fprintf(stderr,"GAIN ");
    for (i=0;i<nbolo;i++) fprintf(stderr,"%lg ",gaingi[i*GAINSTEP]);
    fprintf(stderr,"\n");
    if (GAINSTEP>1) {
      fprintf(stderr,"...\nGAIN ");
      for (i=0;i<nbolo;i++) fprintf(stderr,"%lg ",gaingi[i*GAINSTEP+(GAINSTEP-1)]);
      fprintf(stderr,"\n");
    }
  }
  if (itbogo==0) delta0=0;
  MPI_Barrier(MPI_COMM_WORLD);
  
  struct timeval tp1,tp2;
  gettimeofday(&tp1,NULL);
  
  iter = 0;
  memset(b2  ,0,nmatres*sizeof (double));
  memset(d2  ,0,nmatres*sizeof (double));
  memset(q2  ,0,nmatres*sizeof (double));
  memset(r2  ,0,nmatres*sizeof (double));
  memset(s2  ,0,nmatres*sizeof (double));
  memset(hit2,0,nmatres*sizeof (double));

  //==========================================
  //=  Compute second member
  //=
  //=
  long l,m;

  ////// BUILD B2
  //GetProcMem(&vmem,&phymem);
  //if (rank==0) fprintf(stderr,"Rank: %ld Line=%d MEM %.1lf[%.1lf]MB\n",
  //                  (long) rank, __LINE__,
  //                  (double) vmem/1024./1024.,
  //                  (double) phymem/1024./1024.);

  for (k=0;k<nnbpix;k++)  {
    //long imat=the_stat_pix[k];
    long ndata = loc_nhpix[k];
    hpix *htmp = loc_hpix[k];
    II[k]=0;
    IQ[k]=0;
    IU[k]=0;
    QQ[k]=0;
    UU[k]=0;
    QU[k]=0;
      

    for (l1=0;l1<ndata;l1++) {
      long ri1=htmp[l1].rg-globalBeginRing;
      
      double CO1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].co
                                        -dpsisi[htmp[l1].ib]*htmp[l1].si);
      double SI1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].si
                                        +dpsisi[htmp[l1].ib]*htmp[l1].co);
      if (flg_rg[htmp[l1].ib][ri1]!=0) {
        II[k]+=htmp[l1].w;
        IQ[k]+=htmp[l1].w*CO1;
        IU[k]+=htmp[l1].w*SI1;
        QQ[k]+=htmp[l1].w*CO1*CO1;
        QU[k]+=htmp[l1].w*SI1*CO1;
        UU[k]+=htmp[l1].w*SI1*SI1;
      }
    }

    if (ndata>0&&flgpix[k]>0) {

      double SI=0;
      double SQ=0;
      double SU=0;

#ifdef USEDII
      memset(dii ,0,sizeof(double)*nbolo*GAINSTEP);
      memset(dqq ,0,sizeof(double)*nbolo*GAINSTEP);
      memset(duu ,0,sizeof(double)*nbolo*GAINSTEP);
#endif
      memset(dcoi+k*nbolo,0,sizeof(double)*nbolo);
      memset(dcoq+k*nbolo,0,sizeof(double)*nbolo);
      memset(dcou+k*nbolo,0,sizeof(double)*nbolo);
      memset(dfri+k*nbolo,0,sizeof(double)*nbolo);
      memset(dfrq+k*nbolo,0,sizeof(double)*nbolo);
      memset(dfru+k*nbolo,0,sizeof(double)*nbolo);
      memset(dthetai+k*nbolo,0,sizeof(double)*nbolo);
      memset(dthetaq+k*nbolo,0,sizeof(double)*nbolo);
      memset(dthetau+k*nbolo,0,sizeof(double)*nbolo);
      memset(ddusti+k*nbolo,0,sizeof(double)*nbolo);
      memset(ddustq+k*nbolo,0,sizeof(double)*nbolo);
      memset(ddustu+k*nbolo,0,sizeof(double)*nbolo);
      memset(dpixi+k*nbolo*npixbeam,0,sizeof(double)*nbolo*npixbeam);
      memset(dpixq+k*nbolo*npixbeam,0,sizeof(double)*nbolo*npixbeam);
      memset(dpixu+k*nbolo*npixbeam,0,sizeof(double)*nbolo*npixbeam);
      

      for (l1=0;l1<ndata;l1++) {
        long ri1=htmp[l1].rg-globalBeginRing;
        
        double g1=gaingi[htmp[l1].gi+htmp[l1].ib*GAINSTEP];
        
        double CO1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].co
                                          -dpsisi[htmp[l1].ib]*htmp[l1].si);
        double SI1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].si
                                          +dpsisi[htmp[l1].ib]*htmp[l1].co);
        if (flg_rg[htmp[l1].ib][ri1]!=0) {
          
          htmp[l1].wp=0;
          //htmp[l1].thsig=0;
          
          if (Param->REMHDIP==0) {
            SI+=htmp[l1].w*(htmp[l1].sig*g1-htmp[l1].fsl-htmp[l1].dip-gain[htmp[l1].adu+htmp[l1].ib*GAINSTEPADU]);
            SQ+=htmp[l1].w*(htmp[l1].sig*g1-htmp[l1].fsl-htmp[l1].dip-gain[htmp[l1].adu+htmp[l1].ib*GAINSTEPADU])*CO1;
            SU+=htmp[l1].w*(htmp[l1].sig*g1-htmp[l1].fsl-htmp[l1].dip-gain[htmp[l1].adu+htmp[l1].ib*GAINSTEPADU])*SI1;
          }
          else {
            SI+=htmp[l1].w*(htmp[l1].sig*g1-htmp[l1].fsl-htmp[l1].freefree-gain[htmp[l1].adu+htmp[l1].ib*GAINSTEPADU]);
            SQ+=htmp[l1].w*(htmp[l1].sig*g1-htmp[l1].fsl-htmp[l1].freefree-gain[htmp[l1].adu+htmp[l1].ib*GAINSTEPADU])*CO1;
            SU+=htmp[l1].w*(htmp[l1].sig*g1-htmp[l1].fsl-htmp[l1].freefree-gain[htmp[l1].adu+htmp[l1].ib*GAINSTEPADU])*SI1;            
          }
          
#ifdef USEDII
          dii[htmp[l1].gi+htmp[l1].ib*GAINSTEP ]+=htmp[l1].w*htmp[l1].dip;
          dqq[htmp[l1].gi+htmp[l1].ib*GAINSTEP ]+=htmp[l1].w*htmp[l1].dip*CO1;
          duu[htmp[l1].gi+htmp[l1].ib*GAINSTEP ]+=htmp[l1].w*htmp[l1].dip*SI1;
          
#endif
            
          if (nmatco>0) {
            dcoi[htmp[l1].ib+k*nbolo] += htmp[l1].w*htmp[l1].comap;
            dcoq[htmp[l1].ib+k*nbolo] += htmp[l1].w*htmp[l1].comap*CO1;
            dcou[htmp[l1].ib+k*nbolo] += htmp[l1].w*htmp[l1].comap*SI1;
          }
          if (nmatdust>0) {
            ddusti[htmp[l1].ib+k*nbolo] += htmp[l1].w*htmp[l1].dustmap;
            ddustq[htmp[l1].ib+k*nbolo] += htmp[l1].w*htmp[l1].dustmap*CO1;
            ddustu[htmp[l1].ib+k*nbolo] += htmp[l1].w*htmp[l1].dustmap*SI1;
          }
          if (nfreefree>0) {
            dfri[htmp[l1].ib+k*nbolo] += htmp[l1].w*htmp[l1].freefree;
            dfrq[htmp[l1].ib+k*nbolo] += htmp[l1].w*htmp[l1].freefree*CO1;
            dfru[htmp[l1].ib+k*nbolo] += htmp[l1].w*htmp[l1].freefree*SI1;
          }
          
          if (ittt>0) {
            for (m=0;m<npixbeam;m++)  {
              dpixi[nbolo*m+htmp[l1].ib+k*nbolo*npixbeam] += htmp[l1].w*htmp[l1].listofpix[m];
              dpixq[nbolo*m+htmp[l1].ib+k*nbolo*npixbeam] += htmp[l1].w*htmp[l1].listofpix[m]*CO1;
              dpixu[nbolo*m+htmp[l1].ib+k*nbolo*npixbeam] += htmp[l1].w*htmp[l1].listofpix[m]*SI1;
            }
          } 
        }
      }

      solvemap(&SI,&SQ,&SU,II[k],IQ[k],IU[k],QQ[k],QU[k],UU[k]);
      SSI[k]=SI;
      SSQ[k]=SQ;
      SSU[k]=SU;
      
      for (ib=0;ib<nbolo;ib++) {
        
#ifdef USEDII
        for (j=0;j<GAINSTEP;j++) {
          solvemap(dii+j+ib*GAINSTEP,dqq+j+ib*GAINSTEP,duu+j+ib*GAINSTEP,II[k],IQ[k],IU[k],QQ[k],QU[k],UU[k]);

        }
#endif

        if (nmatco>0) {
          solvemap(dcoi+ib+k*nbolo,dcoq+ib+k*nbolo,dcou+ib+k*nbolo,II[k],IQ[k],IU[k],QQ[k],QU[k],UU[k]);
        }
        if (nmatdust>0) {
          solvemap(ddusti+ib+k*nbolo,ddustq+ib+k*nbolo,ddustu+ib+k*nbolo,II[k],IQ[k],IU[k],QQ[k],QU[k],UU[k]);
        }
        if (nfreefree>0) {
          solvemap(dfri+ib+k*nbolo,dfrq+ib+k*nbolo,dfru+ib+k*nbolo,II[k],IQ[k],IU[k],QQ[k],QU[k],UU[k]);
        }
        
        //for (m=0;m<DEGREE;m++) {
        //  solvemap(ddegi+ib+m*nbolo+k*nbolo*DEGREE,
        //           ddegq+ib+m*nbolo+k*nbolo*DEGREE,
        //           ddegu+ib+m*nbolo+k*nbolo*DEGREE,II[k],IQ[k],IU[k],QQ[k],QU[k],UU[k]);
        //}
        
        if (ittt>0) {
          for (m=0;m<npixbeam;m++)  {
            solvemap(dpixi+ib+m*nbolo+k*nbolo*npixbeam,
                     dpixq+ib+m*nbolo+k*nbolo*npixbeam,
                     dpixu+ib+m*nbolo+k*nbolo*npixbeam,II[k],IQ[k],IU[k],QQ[k],QU[k],UU[k]);
          }
        } 
      }
      for (l1=0;l1<ndata;l1++) {
        long ri1=htmp[l1].rg-globalBeginRing;

        double CO1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].co
                                          -dpsisi[htmp[l1].ib]*htmp[l1].si);
        double SI1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].si
                                          +dpsisi[htmp[l1].ib]*htmp[l1].co);
        if (flg_rg[htmp[l1].ib][ri1]!=0) {
          double li=htmp[l1].w,lco=CO1*htmp[l1].w,lsi=SI1*htmp[l1].w;
          solvemap(&li,&lco,&lsi,II[k],IQ[k],IU[k],QQ[k],QU[k],UU[k]);
          htmp[l1].vi=li;
          htmp[l1].vq=lco;
          htmp[l1].vu=lsi;
#ifndef USEDII
          li=htmp[l1].w*htmp[l1].dip;
          lco=CO1*htmp[l1].w*htmp[l1].dip;
          lsi=SI1*htmp[l1].w*htmp[l1].dip;

          solvemap(&li,&lco,&lsi,II[k],IQ[k],IU[k],QQ[k],QU[k],UU[k]);
          htmp[l1].lvi=li;
          htmp[l1].lvq=lco;
          htmp[l1].lvu=lsi;
#endif
        }
      }
      for (l1=0;l1<ndata;l1++) {
        long ri1=htmp[l1].rg-globalBeginRing;
        
        double g1=gaingi[htmp[l1].gi+htmp[l1].ib*GAINSTEP];
        
        if (flg_rg[htmp[l1].ib][ri1]!=0) {
          double divi=NEP_tab[htmp[l1].ib]*htmp[l1].hit*g1;
          divi=divi*divi;
          
          if (divi==0) htmp[l1].wp=0;
          else htmp[l1].wp=1/divi;
          
          if (itbogo==0) normaoff+=NEP_tab[htmp[l1].ib]*htmp[l1].hit;
          
        }
      }
      
      for (l1=0;l1<ndata;l1++) {
        long ri1=htmp[l1].rg-globalBeginRing;
        long iri1=rgord[htmp[l1].ib][ri1]+newnr[htmp[l1].ib];
        
        double g1=gaingi[htmp[l1].gi+htmp[l1].ib*GAINSTEP];
        
        double CO1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].co
                                          -dpsisi[htmp[l1].ib]*htmp[l1].si);
        double SI1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].si
                                          +dpsisi[htmp[l1].ib]*htmp[l1].co);
        if (flg_rg[htmp[l1].ib][ri1]!=0) {
          double ww=htmp[l1].wp; 
          double tmp;
          if (REMHDIP==0) {
            tmp=((htmp[l1].sig*g1-htmp[l1].fsl-htmp[l1].dip-gain[htmp[l1].adu+htmp[l1].ib*GAINSTEPADU])-(SI+CO1*SQ+SI1*SU));
          }
          else {
            tmp=((htmp[l1].sig*g1-htmp[l1].fsl-htmp[l1].freefree-gain[htmp[l1].adu+htmp[l1].ib*GAINSTEPADU])-(SI+CO1*SQ+SI1*SU));
          }

          long l3;
#ifndef USEDII
          memset(cdip,0,GAINSTEP*nbolo*sizeof(double));
#endif
          
          for (l3=0;l3<ndata;l3++) {
            long ri3=htmp[l3].rg-globalBeginRing;
            if (flg_rg[htmp[l3].ib][ri3]!=0) {
              long ir3=rgord[htmp[l3].ib][ri3]+newnr[htmp[l3].ib];
              ctmp[ir3]=-(htmp[l3].vi+CO1*htmp[l3].vq+SI1*htmp[l3].vu);
#ifndef USEDII
              cdip[htmp[l3].ib*GAINSTEP+htmp[l3].gi]-=(htmp[l3].lvi+CO1*htmp[l3].lvq+SI1*htmp[l3].lvu);
#endif
            }
          }
          ctmp[iri1]+=1;
             
#ifdef USEDII
          for (ib=0;ib<nbolo;ib++) 
            for (j=0;j<GAINSTEP;j++) {
              cdip[ib*GAINSTEP+j]=-(dii[j+ib*GAINSTEP]+CO1*dqq[j+ib*GAINSTEP]+SI1*duu[j+ib*GAINSTEP]);
            }
#endif
          cdip[htmp[l1].ib*GAINSTEP+htmp[l1].gi]+=htmp[l1].dip;

          for (ib=0;ib<nbolo;ib++) cco[ib]=-(dcoi[ib+k*nbolo]+CO1*dcoq[ib+k*nbolo]+SI1*dcou[ib+k*nbolo]);
          cco[htmp[l1].ib]+=htmp[l1].comap;
          
          for (ib=0;ib<nbolo;ib++) cdust[ib]=-(ddusti[ib+k*nbolo]+CO1*ddustq[ib+k*nbolo]+SI1*ddustu[ib+k*nbolo]);
          cdust[htmp[l1].ib]+=htmp[l1].dustmap;
          
          for (ib=0;ib<nbolo;ib++) ccfree[ib]=-(dfri[ib+k*nbolo]+CO1*dfrq[ib+k*nbolo]+SI1*dfru[ib+k*nbolo]);
          ccfree[htmp[l1].ib]+=htmp[l1].freefree;
          
          for (m=0;m<npixbeam;m++) {
            for (ib=0;ib<nbolo;ib++) cpix[ib+m*nbolo]=-(dpixi[ib+m*nbolo+k*nbolo*npixbeam]
                                                        +CO1*dpixq[ib+m*nbolo+k*nbolo*npixbeam]
                                                        +SI1*dpixu[ib+m*nbolo+k*nbolo*npixbeam]);
            cpix[htmp[l1].ib+m*nbolo]+=htmp[l1].listofpix[m];
          }
          long ir;
          /////////////////  OFFSET
          
          for (l2=0;l2<ndata;l2++) {
            long ri2=htmp[l2].rg-globalBeginRing;
            if (flg_rg[htmp[l2].ib][ri2]!=0) {
              long ir=rgord[htmp[l2].ib][ri2]+newnr[htmp[l2].ib]; 
              b2[ir]+=ww*tmp*ctmp[ir];
              hit2[ir]+=ww*ctmp[ir]*ctmp[ir];
#ifndef USEDII
              b2[newnr[nbolo]+htmp[l2].ib*GAINSTEP+htmp[l2].gi]+=ww*tmp*cdip[htmp[l2].ib*GAINSTEP+htmp[l2].gi];
              hit2[newnr[nbolo]+htmp[l2].ib*GAINSTEP+htmp[l2].gi]+=ww*cdip[htmp[l2].ib*GAINSTEP+htmp[l2].gi]
                *cdip[htmp[l2].ib*GAINSTEP+htmp[l2].gi];
#endif
            }
          }
          
          
          /////////////////  DIPOLE FIT
          
          for (ir=0;ir<nbolo;ir++) {
            
#ifdef USEDII
            
            for (j=0;j<GAINSTEP;j++) {
              b2[newnr[nbolo]+ir*GAINSTEP+j]+=ww*tmp*cdip[ir*GAINSTEP+j];
              hit2[newnr[nbolo]+ir*GAINSTEP+j]+=ww*cdip[ir*GAINSTEP+j]*cdip[ir*GAINSTEP+j];
            }
            
#endif
              
            ///////////// CO
            if (nmatco>0) {
              b2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+ir]+=ww*tmp*cco[ir];
              hit2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+ir]+=ww*cco[ir]*cco[ir];
            }
            
            ///////////// DUST
            
            if (nmatdust>0) {
              b2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+ir]+=ww*tmp*cdust[ir];
              hit2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+ir]+=ww*cdust[ir]*cdust[ir];
            }
            
            ///////////// FREEFREE
            
            if (nfreefree>0) {
              b2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+ir]+=ww*tmp*ccfree[ir];
              hit2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+ir]+=ww*ccfree[ir]*ccfree[ir];
            }
            
            ////////// PIXBEAM
            for (j=0;j<npixbeam;j++) {
              b2[newnr[nbolo]+nbolo*GAINSTEP2+j*nbolo+ir]+=ww*tmp*cpix[ir+j*nbolo];
              hit2[newnr[nbolo]+nbolo*GAINSTEP2+j*nbolo+ir]+=ww*cpix[ir+j*nbolo]*cpix[ir+j*nbolo];
            }
          }
        }
      }
    }
  }


  if (rank==0) {
    double *lb = (double *) malloc(sizeof(double)*(nmatres));
    for (rrk=1;rrk<mpi_size;rrk++) {
      MPI_Recv(lb,sizeof(double)*(nmatres), MPI_BYTE, rrk,1031, MPI_COMM_WORLD,&statu);
      for (l=0;l<nmatres;l++) b2[l]+=lb[l];
    }
    free(lb);
  }
  else MPI_Send(b2, sizeof(double)*(nmatres), MPI_BYTE, 0, 1031, MPI_COMM_WORLD);
  
  if (rank==0) {
    double *lb = (double *) malloc(sizeof(double)*(nmatres));
    for (rrk=1;rrk<mpi_size;rrk++) {
      MPI_Recv(lb,sizeof(double)*(nmatres), MPI_BYTE, rrk,1033, MPI_COMM_WORLD,&statu);

      for (l=0;l<nmatres;l++) hit2[l]+=lb[l];
    }
    free(lb);
  }
  else {
    MPI_Send(hit2, sizeof(double)*(nmatres), MPI_BYTE, 0, 1033, MPI_COMM_WORLD);
  }
  

  //==========================================================
  // Compute Ax
  //healspline_fill_s_with_PtP_s (hs_rec, Ndata, vec, x, q);
  // + 
  // Preconditionnement
  //
  //
  
  for (l=0;l<nmatres;l++) q2[l]=0;
  if (rank==0) {

    double soff=0;
    for (i=0;i<newnr[nbolo];i++) soff+=hit2[0]*ix2[i];
    for (i=0;i<newnr[nbolo];i++) q2[i]=soff;

    if (Param->REMHDIP==1) {
      soff=0;
      for (i=newnr[nbolo];i<newnr[nbolo]+GAINSTEP*nbolo;i++) soff+=hit2[newnr[nbolo]]*ix2[i];
      for (i=newnr[nbolo];i<newnr[nbolo]+GAINSTEP*nbolo;i++) q2[i]=soff;
    }

    if (Param->flag_AVGR0==_PAR_TRUE) {
      soff=0;
      for (i=0;i<nbolo;i++) soff+=hit2[newnr[nbolo]+nbolo*(GAINSTEP2)]*
                              ix2[newnr[nbolo]+nbolo*(GAINSTEP2)+i];
      for (i=0;i<nbolo;i++) q2[newnr[nbolo]+nbolo*(GAINSTEP2)+i]=soff;
      for (i=0;i<nbolo;i++) b2[newnr[nbolo]+nbolo*(GAINSTEP2)+i]+=
                              hit2[newnr[nbolo]+nbolo*(GAINSTEP2)]*Param->AVGR0;
      
    }

    if (nmatco>0) {
      soff=0;
      for (i=0;i<nbolo;i++) soff+=hit2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)]
                              *ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+i];
      for (i=0;i<nbolo;i++) q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+i]=soff;
      for (i=0;i<nbolo;i++) b2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+i]+=
                              hit2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)]*Param->AVGCO;
    }
    if (nmatdust>0) {
      soff=0;
      for (i=0;i<nbolo;i++) 
        soff+=hit2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco]*
          ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+i];
      for (i=0;i<nbolo;i++) q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+i]=soff;
      for (i=0;i<nbolo;i++) b2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+i]+=
                              hit2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco]*Param->AVGDUST;
    }
    if (nfreefree>0) {
      soff=0;
      for (i=0;i<nbolo;i++) 
        soff+=hit2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust]*
          ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+i];
      for (i=0;i<nbolo;i++) q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+i]=soff;
    }
    if (DOFITTHETA==1) {
      soff=0;
      for (i=0;i<nbolo;i++) soff+=hit2[newnr[nbolo]+nbolo*(GAINSTEP2+npixbeam-1)]*
                              ix2[newnr[nbolo]+nbolo*(GAINSTEP2+npixbeam-1)+i];
      for (i=0;i<nbolo;i++) q2[newnr[nbolo]+nbolo*(GAINSTEP2+npixbeam-1)+i]=soff;
    }
  }

  for (k=0;k<nnbpix;k++) if (flgpix[k]>0) {
    
    long ndata = loc_nhpix[k];
    if (ndata>0) {
#ifdef TIMING
      gettimeofday(&tp1,NULL);
#endif
      hpix *htmp = loc_hpix[k];

      double vali=0,valq=0,valu=0;

#ifdef USEDII
      memset(dii ,0,sizeof(double)*nbolo*GAINSTEP);
      memset(dqq ,0,sizeof(double)*nbolo*GAINSTEP);
      memset(duu ,0,sizeof(double)*nbolo*GAINSTEP);


      for (l1=0;l1<ndata;l1++) {
        long ri1=htmp[l1].rg-globalBeginRing;
        double CO1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].co
                                          -dpsisi[htmp[l1].ib]*htmp[l1].si);
        double SI1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].si
                                          +dpsisi[htmp[l1].ib]*htmp[l1].co);
        if (flg_rg[htmp[l1].ib][ri1]!=0) {
          
          dii[htmp[l1].gi+htmp[l1].ib*GAINSTEP ]+=htmp[l1].w*htmp[l1].dip;
          dqq[htmp[l1].gi+htmp[l1].ib*GAINSTEP ]+=htmp[l1].w*htmp[l1].dip*CO1;
          duu[htmp[l1].gi+htmp[l1].ib*GAINSTEP ]+=htmp[l1].w*htmp[l1].dip*SI1;
        }
      }

      for (ib=0;ib<nbolo;ib++) {
        for (j=0;j<GAINSTEP;j++) {
          solvemap(dii+j+ib*GAINSTEP,dqq+j+ib*GAINSTEP,duu+j+ib*GAINSTEP,II[k],IQ[k],IU[k],QQ[k],QU[k],UU[k]);
        }
      } 
#endif
        
        
      long l3;
      for (l3=0;l3<ndata;l3++) {
        long ri3=htmp[l3].rg-globalBeginRing;
        if (flg_rg[htmp[l3].ib][ri3]!=0) {
          long ir3=rgord[htmp[l3].ib][ri3]+newnr[htmp[l3].ib];
          vali+=htmp[l3].vi*ix2[ir3];
          valq+=htmp[l3].vq*ix2[ir3];
          valu+=htmp[l3].vu*ix2[ir3];
#ifndef USEDII
          vali+=htmp[l3].lvi*ix2[newnr[nbolo]+htmp[l3].ib*GAINSTEP+htmp[l3].gi];
          valq+=htmp[l3].lvq*ix2[newnr[nbolo]+htmp[l3].ib*GAINSTEP+htmp[l3].gi];
          valu+=htmp[l3].lvu*ix2[newnr[nbolo]+htmp[l3].ib*GAINSTEP+htmp[l3].gi];
#endif
        }
      }
      //fprintf(stderr,"I0 vali %lg %lg %lg\n",vali,valq,valu);
      
      for (ib=0;ib<nbolo;ib++) {
#ifdef USEDII
        for (j=0;j<GAINSTEP;j++) {
          vali+=dii[j+ib*GAINSTEP]  *ix2[newnr[nbolo]+ib*GAINSTEP+j];
          valq+=dqq[j+ib*GAINSTEP]  *ix2[newnr[nbolo]+ib*GAINSTEP+j];
          valu+=duu[j+ib*GAINSTEP]  *ix2[newnr[nbolo]+ib*GAINSTEP+j];
        }
        //fprintf(stderr,"I1 vali %lg %lg %lg\n",vali,valq,valu);
#endif

        if (nmatco>0) {
          vali+=dcoi[ib+k*nbolo]  *ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+ib];
          valq+=dcoq[ib+k*nbolo]  *ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+ib];
          valu+=dcou[ib+k*nbolo]  *ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+ib];
        }            
        if (nmatdust>0) {
          vali+=ddusti[ib+k*nbolo]  *ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+ib];
          valq+=ddustq[ib+k*nbolo]  *ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+ib];
          valu+=ddustu[ib+k*nbolo]  *ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+ib];
        }            
        if (nfreefree>0) {
          vali+=dfri[ib+k*nbolo]  *ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+ib];
          valq+=dfrq[ib+k*nbolo]  *ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+ib];
          valu+=dfru[ib+k*nbolo]  *ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+ib];
        }
        //fprintf(stderr,"I3 vali %lg %lg %lg : %lg %lg %lg %lg\n",vali,valq,valu,
        //dfri[ib+k*nbolo],dfrq[ib+k*nbolo],dfru[ib+k*nbolo],
        //ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+ib]);            
        
        for (m=0;m<npixbeam;m++) {
          vali+=ix2[newnr[nbolo]+nbolo*GAINSTEP2+ib+m*nbolo]*dpixi[ib+m*nbolo+k*nbolo*npixbeam];
          valq+=ix2[newnr[nbolo]+nbolo*GAINSTEP2+ib+m*nbolo]*dpixq[ib+m*nbolo+k*nbolo*npixbeam];
          valu+=ix2[newnr[nbolo]+nbolo*GAINSTEP2+ib+m*nbolo]*dpixu[ib+m*nbolo+k*nbolo*npixbeam];
        }  
        
        //fprintf(stderr,"I4 vali %lg %lg %lg\n",vali,valq,valu);
      }
      double qri=0;
      double qrq=0;
      double qru=0;
      
      for (l1=0;l1<ndata;l1++) {
        long ri1=htmp[l1].rg-globalBeginRing;
        long iri1=rgord[htmp[l1].ib][ri1]+newnr[htmp[l1].ib];
        
        double CO1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].co
                                          -dpsisi[htmp[l1].ib]*htmp[l1].si);
        double SI1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].si
                                          +dpsisi[htmp[l1].ib]*htmp[l1].co);
        if (flg_rg[htmp[l1].ib][ri1]!=0) {
          double ww=htmp[l1].wp;
          double val2=ix2[iri1]-(vali+CO1*valq+SI1*valu);
          
          //fprintf(stderr,"I1 val2 %lg %lg %lg\n",val2,ix2[iri1],vali+CO1*valq+SI1*valu);
          val2+=ix2[newnr[nbolo]+htmp[l1].ib*GAINSTEP+htmp[l1].gi]*htmp[l1].dip;

          //fprintf(stderr,"I2 val2 %lg\n",val2);
          if (nmatco>0)   
            val2+=ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+htmp[l1].ib]*htmp[l1].comap;
          if (nmatdust>0) 
            val2+=ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+htmp[l1].ib]*htmp[l1].dustmap;
          if (nfreefree>0) 
            val2+=ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+htmp[l1].ib]*htmp[l1].freefree;
          
          //fprintf(stderr,"I3 val2 %lg\n",val2);
          for (m=0;m<npixbeam;m++) {
            val2+=ix2[newnr[nbolo]+nbolo*GAINSTEP2+htmp[l1].ib+m*nbolo]*htmp[l1].listofpix[m];
          }
          
          
          qri-=ww*val2;
          qrq-=ww*val2*CO1;
          qru-=ww*val2*SI1;
          
          //fprintf(stderr,"val2 I4 %lg\n",val2);
          q2[iri1]+=ww*val2;
          
          q2[newnr[nbolo]+htmp[l1].ib*GAINSTEP+htmp[l1].gi]+=ww*val2*htmp[l1].dip;

          ///////////// CO
          if (nmatco>0) {
            q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+htmp[l1].ib]+=ww*val2*htmp[l1].comap;
          }
          
          ///////////// DUST
          if (nmatdust>0) {
            q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+htmp[l1].ib]+=ww*val2*htmp[l1].dustmap;
          }
          
          ///////////// FREEFREE
          if (nfreefree>0) {
            q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+htmp[l1].ib]+=ww*val2*htmp[l1].freefree;
          }
          
          ////////// SYSTE
          for (j=0;j<npixbeam;j++) {
            q2[newnr[nbolo]+nbolo*GAINSTEP2+j*nbolo+htmp[l1].ib]+=ww*val2*htmp[l1].listofpix[j];  
          }
          
          
        }
      }
    

      for (l2=0;l2<ndata;l2++) {
        long ri2=htmp[l2].rg-globalBeginRing;
        if (flg_rg[htmp[l2].ib][ri2]!=0) {
          long ir=rgord[htmp[l2].ib][ri2]+newnr[htmp[l2].ib]; 
          q2[ir]+=qri*htmp[l2].vi+qrq*htmp[l2].vq+qru*htmp[l2].vu;
#ifndef USEDII
          q2[newnr[nbolo]+htmp[l2].ib*GAINSTEP+htmp[l2].gi]+=qri*htmp[l2].lvi+qrq*htmp[l2].lvq+qru*htmp[l2].lvu;
#endif
        }
      }
    
        
      for (ib=0;ib<nbolo;ib++) {
        
#ifdef USEDII
        for (j=0;j<GAINSTEP;j++) {
          q2[newnr[nbolo]+GAINSTEP*ib+j]+=qri*dii[j+ib*GAINSTEP]+qrq*dqq[j+ib*GAINSTEP]+qru*duu[j+ib*GAINSTEP];
        }
#endif

        ///////////// CO
        if (nmatco>0) {
          q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+ib]+=qri*dcoi[ib+k*nbolo]+qrq*dcoq[ib+k*nbolo]+qru*dcou[ib+k*nbolo];
        }
        ///////////// DUST
        if (nmatdust>0) {
          q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+ib]+=qri*ddusti[ib+k*nbolo]+qrq*ddustq[ib+k*nbolo]+qru*ddustu[ib+k*nbolo];
        }
        ///////////// FREEFREE
        if (nfreefree>0) {
          q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+ib]+=qri*dfri[ib+k*nbolo]+qrq*dfrq[ib+k*nbolo]+qru*dfru[ib+k*nbolo];
        }
        
        for (j=0;j<npixbeam;j++) {
          q2[newnr[nbolo]+nbolo*GAINSTEP2+j*nbolo+ib]+=qri*dpixi[ib+j*nbolo+k*nbolo*npixbeam]+qrq*dpixq[ib+j*nbolo+k*nbolo*npixbeam]+qru*dpixu[ib+j*nbolo+k*nbolo*npixbeam];  
        }
      }
    }
      
#ifdef TIMING
      gettimeofday(&tp2,NULL);
      double dt=(double)(tp2.tv_sec-tp1.tv_sec)+(1E-6)*(tp2.tv_usec-tp1.tv_usec);
      dthit[ndata]+=dt;
      ndthit[ndata]+=1;
#endif
    }

  if (rank==0) {
    double *lb = (double *) malloc(sizeof(double)*(nmatres));
    for (rrk=1;rrk<mpi_size;rrk++) {
      MPI_Recv(lb,sizeof(double)*(nmatres), MPI_BYTE, rrk,1032, MPI_COMM_WORLD,&statu);
      for (l=0;l<nmatres;l++) q2[l]+=lb[l];
    }
    free(lb);
  }
  else {
    MPI_Send(q2, sizeof(double)*(nmatres), MPI_BYTE, 0, 1032, MPI_COMM_WORLD);
  }

  if (rank==0) fprintf(stderr,"QQ2 %lg\n",q2[0]);

  if (rank==0) {
    for (i=0; i < nmatres; i++)
      {
        r2[i] = b2[i] - q2[i];
        d2[i] = r2[i] / hit2[i];
      }
  }
  
  double delta_new_tmp = 0.0;
  if (rank==0) {
    for (i = 0; i < nmatres; i++) {
      delta_new_tmp += b2[i] ;
      if (isnan(b2[i])) {
        fprintf(stderr,"NAN B2 PBS %ld\n",(long) i);
      }
      if (isnan(d2[i])) {
        fprintf(stderr,"NAN D2 PBS %ld %lg\n",(long) i,hit2[i]);
      }
    }
    //fprintf(stderr,"B2 %lg\n",delta_new_tmp);
  
    delta_new_tmp = 0.0;
    for (i = 0; i < nmatres; i++) {
      delta_new_tmp += q2[i] ;
      if (isnan(q2[i])) {
        fprintf(stderr,"NAN Q2 PBS %ld\n",(long) i);
      }
    }
    //fprintf(stderr,"Q2 %lg\n",delta_new_tmp);
    delta_new_tmp = 0.0;
    for (i = 0; i < nmatres; i++) {
      delta_new_tmp += q2[i]-b2[i] ;
      //fprintf(stderr,"B2 Q2 B2-Q2 [%ld]: %lg\t%lg\t%lg\n",(long) i,b2[i],q2[i],q2[i]-b2[i]);
    }
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
  delta_new_tmp = 0.0;
  if (rank==0) for (i=0; i < nmatres; i++) {
    delta_new_tmp += r2[i] * d2[i];
  }
  
  delta_new=0;
  for (i=0;i<mpi_size;i++) {
    double tmp=delta_new_tmp;
    MPI_Bcast(&tmp, sizeof(double), MPI_BYTE, i, MPI_COMM_WORLD);
    delta_new+=tmp;
  }
  if (itbogo==0) delta0 = delta_new;
  if (rank==0) fprintf (stderr, "iter = %d - delta0 = %lg - delta_new = %lg\n", iter, delta0, delta_new);
  int testwrit=0;

  if (itbogo==0) MPI_Bcast(&normaoff, sizeof(double), MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&delta0, sizeof(double), MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&delta_new, sizeof(double), MPI_BYTE, 0, MPI_COMM_WORLD);


  while ((iter < itermax)  && ((delta_new) > delta0*1E-24) && ((delta_new) > 1E-20)) //Param->XI2STOP))
    {
      // q <= Ad
      //if (rank==0&&mindelta>delta_new) {
      if (rank==0) {
        memcpy(x2old,ix2,nmatres*sizeof(double));
        testwrit=1;
      }
      else testwrit=0;

      //healspline_fill_s_with_PtP_s (hs_rec, Ndata, vec, d, q);
      //for (i=0;i<mpi_size;i++) {
      //MPI_Bcast(d+tab_begr[i]*2, sizeof(double)*(tab_edr[i]-tab_begr[i]+1)*2, MPI_BYTE, i, MPI_COMM_WORLD);
      //}
      MPI_Bcast(d2, sizeof(double)*(nmatres), MPI_BYTE, 0, MPI_COMM_WORLD);
      
      // ===========================================================================
      // PROJECTION DE d dans q
      //
      
      for (l=0;l<nmatres;l++) q2[l]=0;
      if (rank==0) {
        double soff=0;
        for (i=0;i<newnr[nbolo];i++) soff+=hit2[0]*d2[i];
        for (i=0;i<newnr[nbolo];i++) q2[i]=soff;
        
        if (Param->REMHDIP==1) {
          soff=0;
          for (i=newnr[nbolo];i<newnr[nbolo]+GAINSTEP*nbolo;i++) soff+=hit2[newnr[nbolo]]*d2[i];
          for (i=newnr[nbolo];i<newnr[nbolo]+GAINSTEP*nbolo;i++) q2[i]=soff;
        }

        if (Param->flag_AVGR0==_PAR_TRUE) {
          soff=0;
          for (i=0;i<nbolo;i++) soff+=hit2[newnr[nbolo]+nbolo*(GAINSTEP2)]*
                                  d2[newnr[nbolo]+nbolo*(GAINSTEP2)+i];
          for (i=0;i<nbolo;i++) q2[newnr[nbolo]+nbolo*(GAINSTEP2)+i]=soff; 
        }
        if (nmatco>0) {
          soff=0;
          for (i=0;i<nbolo;i++) soff+=hit2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)]
                                  *d2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+i];
          for (i=0;i<nbolo;i++) q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+i]=soff;
        }
        if (nfreefree>0) {
          soff=0;
          for (i=0;i<nbolo;i++) soff+=hit2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust]
                                  *d2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+i];
          for (i=0;i<nbolo;i++) q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+i+nmatco+nmatdust]=soff;
        }
        if (nmatdust>0) {
          soff=0;
          for (i=0;i<nbolo;i++) 
            soff+=hit2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco]*
              d2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+i];
          for (i=0;i<nbolo;i++) q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+i]=soff;
        }
        if (DOFITTHETA==1) {
          soff=0;
          for (i=0;i<nbolo;i++) soff+=hit2[newnr[nbolo]+nbolo*(GAINSTEP2+npixbeam-1)]*
                                  d2[newnr[nbolo]+nbolo*(GAINSTEP2+npixbeam-1)+i];
          for (i=0;i<nbolo;i++) q2[newnr[nbolo]+nbolo*(GAINSTEP2+npixbeam-1)+i]=soff;
        }
      }
      for (k=0;k<nnbpix;k++) if (flgpix[k]>0) {
        
        //long imat=the_stat_pix[k];
        long ndata = loc_nhpix[k];
        if (ndata>0) {
          hpix *htmp = loc_hpix[k];

          double vali=0,valq=0,valu=0;

#ifdef USEDII
          memset(dii ,0,sizeof(double)*nbolo*GAINSTEP);
          memset(dqq ,0,sizeof(double)*nbolo*GAINSTEP);
          memset(duu ,0,sizeof(double)*nbolo*GAINSTEP);


          for (l1=0;l1<ndata;l1++) {
            long ri1=htmp[l1].rg-globalBeginRing;
            double CO1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].co
                                              -dpsisi[htmp[l1].ib]*htmp[l1].si);
            double SI1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].si
                                              +dpsisi[htmp[l1].ib]*htmp[l1].co);
            if (flg_rg[htmp[l1].ib][ri1]!=0) {
              
              dii[htmp[l1].gi+htmp[l1].ib*GAINSTEP ]+=htmp[l1].w*htmp[l1].dip;
              dqq[htmp[l1].gi+htmp[l1].ib*GAINSTEP ]+=htmp[l1].w*htmp[l1].dip*CO1;
              duu[htmp[l1].gi+htmp[l1].ib*GAINSTEP ]+=htmp[l1].w*htmp[l1].dip*SI1;
            }
          }
          for (ib=0;ib<nbolo;ib++) {
            for (j=0;j<GAINSTEP;j++) {
              solvemap(dii+j+ib*GAINSTEP,dqq+j+ib*GAINSTEP,duu+j+ib*GAINSTEP,II[k],IQ[k],IU[k],QQ[k],QU[k],UU[k]);
            }
          } 
#endif
          long l3;
          for (l3=0;l3<ndata;l3++) {
            long ri3=htmp[l3].rg-globalBeginRing;
            if (flg_rg[htmp[l3].ib][ri3]!=0) {
              long ir3=rgord[htmp[l3].ib][ri3]+newnr[htmp[l3].ib];
              vali+=htmp[l3].vi*d2[ir3];
              valq+=htmp[l3].vq*d2[ir3];
              valu+=htmp[l3].vu*d2[ir3];
#ifndef USEDII
              vali+=htmp[l3].lvi*d2[newnr[nbolo]+htmp[l3].ib*GAINSTEP+htmp[l3].gi];
              valq+=htmp[l3].lvq*d2[newnr[nbolo]+htmp[l3].ib*GAINSTEP+htmp[l3].gi];
              valu+=htmp[l3].lvu*d2[newnr[nbolo]+htmp[l3].ib*GAINSTEP+htmp[l3].gi];
#endif
            }
          }
          for (ib=0;ib<nbolo;ib++) {
            
            
#ifdef USEDII
            for (j=0;j<GAINSTEP;j++) {
              vali+=dii[j+ib*GAINSTEP]  *d2[newnr[nbolo]+ib*GAINSTEP+j];
              valq+=dqq[j+ib*GAINSTEP]  *d2[newnr[nbolo]+ib*GAINSTEP+j];
              valu+=duu[j+ib*GAINSTEP]  *d2[newnr[nbolo]+ib*GAINSTEP+j];
            }
#endif

            if (nmatco>0) {
              vali+=dcoi[ib+k*nbolo]  *d2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+ib];
              valq+=dcoq[ib+k*nbolo]  *d2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+ib];
              valu+=dcou[ib+k*nbolo]  *d2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+ib];
            }            
            if (nmatdust>0) {
              vali+=ddusti[ib+k*nbolo]  *d2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+ib];
              valq+=ddustq[ib+k*nbolo]  *d2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+ib];
              valu+=ddustu[ib+k*nbolo]  *d2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+ib];
            }            
            if (nfreefree>0) {
              vali+=dfri[ib+k*nbolo]  *d2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+ib];
              valq+=dfrq[ib+k*nbolo]  *d2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+ib];
              valu+=dfru[ib+k*nbolo]  *d2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+ib];
            }            
            for (m=0;m<npixbeam;m++) {
              vali+=d2[newnr[nbolo]+nbolo*GAINSTEP2+ib+m*nbolo]*dpixi[ib+m*nbolo+k*nbolo*npixbeam];
              valq+=d2[newnr[nbolo]+nbolo*GAINSTEP2+ib+m*nbolo]*dpixq[ib+m*nbolo+k*nbolo*npixbeam];
              valu+=d2[newnr[nbolo]+nbolo*GAINSTEP2+ib+m*nbolo]*dpixu[ib+m*nbolo+k*nbolo*npixbeam];
            }   
          }

          double qri=0;
          double qrq=0;
          double qru=0;

          for (l1=0;l1<ndata;l1++) {
            long ri1=htmp[l1].rg-globalBeginRing;
            long iri1=rgord[htmp[l1].ib][ri1]+newnr[htmp[l1].ib];

            double CO1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].co
                                              -dpsisi[htmp[l1].ib]*htmp[l1].si);
            double SI1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].si
                                              +dpsisi[htmp[l1].ib]*htmp[l1].co);
            if (flg_rg[htmp[l1].ib][ri1]!=0) {
              double ww=htmp[l1].wp;
              double val2=d2[iri1]-(vali+CO1*valq+SI1*valu);

              val2+=d2[newnr[nbolo]+htmp[l1].ib*GAINSTEP+htmp[l1].gi]*htmp[l1].dip;
              if (nmatco>0)   
                val2+=d2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+htmp[l1].ib]*htmp[l1].comap;
              if (nmatdust>0) 
                val2+=d2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+htmp[l1].ib]*htmp[l1].dustmap;
              if (nfreefree>0) 
                val2+=d2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+htmp[l1].ib]*htmp[l1].freefree;
                
              for (m=0;m<npixbeam;m++) {
                val2+=d2[newnr[nbolo]+nbolo*GAINSTEP2+htmp[l1].ib+m*nbolo]*htmp[l1].listofpix[m];
              }
              
                
              qri-=ww*val2;
              qrq-=ww*val2*CO1;
              qru-=ww*val2*SI1;
                
              q2[iri1]+=ww*val2;

              q2[newnr[nbolo]+htmp[l1].ib*GAINSTEP+htmp[l1].gi]+=ww*val2*htmp[l1].dip;

              ///////////// CO
              if (nmatco>0) {
                q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+htmp[l1].ib]+=ww*val2*htmp[l1].comap;
              }
                
              ///////////// DUST
              if (nmatdust>0) {
                q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+htmp[l1].ib]+=ww*val2*htmp[l1].dustmap;
              }
                
              ///////////// FREEFREE
              if (nfreefree>0) {
                q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+htmp[l1].ib]+=ww*val2*htmp[l1].freefree;
              }
            
                
              for (j=0;j<npixbeam;j++) {
                q2[newnr[nbolo]+nbolo*GAINSTEP2+j*nbolo+htmp[l1].ib]+=ww*val2*htmp[l1].listofpix[j];  
              }
                
            }
          }


          for (l2=0;l2<ndata;l2++) {
            long ri2=htmp[l2].rg-globalBeginRing;
            if (flg_rg[htmp[l2].ib][ri2]!=0) {
              long ir=rgord[htmp[l2].ib][ri2]+newnr[htmp[l2].ib]; 
              q2[ir]+=qri*htmp[l2].vi+qrq*htmp[l2].vq+qru*htmp[l2].vu;
#ifndef USEDII
              q2[newnr[nbolo]+htmp[l2].ib*GAINSTEP+htmp[l2].gi]+=qri*htmp[l2].lvi+qrq*htmp[l2].lvq+qru*htmp[l2].lvu;
#endif
            }
          }
            
          for (ib=0;ib<nbolo;ib++) {


#ifdef USEDII
            for (j=0;j<GAINSTEP;j++) {
              q2[newnr[nbolo]+GAINSTEP*ib+j]+=qri*dii[j+ib*GAINSTEP]+qrq*dqq[j+ib*GAINSTEP]+qru*duu[j+ib*GAINSTEP];
            }
#endif

            ///////////// CO
            if (nmatco>0) {
              q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+ib]+=qri*dcoi[ib+k*nbolo]+qrq*dcoq[ib+k*nbolo]+qru*dcou[ib+k*nbolo];
            }
            ///////////// DUST
            if (nmatdust>0) {
              q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+ib]+=qri*ddusti[ib+k*nbolo]+qrq*ddustq[ib+k*nbolo]+qru*ddustu[ib+k*nbolo];
            }
            ///////////// FREEFREE
            if (nfreefree>0) {
              q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+ib]+=qri*dfri[ib+k*nbolo]+qrq*dfrq[ib+k*nbolo]+qru*dfru[ib+k*nbolo];
            }

            for (j=0;j<npixbeam;j++) {
              q2[newnr[nbolo]+nbolo*GAINSTEP2+j*nbolo+ib]+=qri*dpixi[ib+j*nbolo+k*nbolo*npixbeam]+qrq*dpixq[ib+j*nbolo+k*nbolo*npixbeam]+qru*dpixu[ib+j*nbolo+k*nbolo*npixbeam];  
            }
          }
        }
      }
      if (rank==0) {
        double *lb = (double *) malloc(sizeof(double)*(nmatres));
        for (rrk=1;rrk<mpi_size;rrk++) {
          MPI_Recv(lb,sizeof(double)*(nmatres), MPI_BYTE, rrk,1034, MPI_COMM_WORLD,&statu);
          for (l=0;l<nmatres;l++) q2[l]+=lb[l];
        }
        free(lb);
      }
      else {
        MPI_Send(q2, sizeof(double)*(nmatres), MPI_BYTE, 0, 1034, MPI_COMM_WORLD);
      }
      double dtq = 0.0;
      if (rank==0) {
        for (i=0; i < nmatres; i++) dtq += d2[i] * q2[i];
        alpha = delta_new / dtq;
        for (i=0; i < nmatres ; i++) ix2[i] += alpha * d2[i];
      }


      gettimeofday(&tp2,NULL);
      if (rank==0) fprintf (stderr,"iter = %ld-%d - delta_new = %.3lg %ld %.3lfs\n", itbogo, iter, delta_new,
                          (long) testwrit,(double)(tp2.tv_sec-tp1.tv_sec)+(1E-6)*(tp2.tv_usec-tp1.tv_usec));

      //if (rank==0) fprintf (stderr,".");
      gettimeofday(&tp1,NULL);
          
      if (iter % 100 == 0 && iter !=0)
        {
          // Use the best case
          memcpy(ix2,x2old,nmatres*sizeof(double));
          //for (i=0;i<mpi_size;i++) {
          //  MPI_Bcast(x+tab_begr[i]*2, sizeof(double)*(tab_edr[i]-tab_begr[i]+1)*2, MPI_BYTE, i, MPI_COMM_WORLD);
          //}
          MPI_Bcast(ix2, sizeof(double)*(nmatres), MPI_BYTE, 0, MPI_COMM_WORLD);

          
      // ===========================================================================
      // PROJECTION DE x dans q
      //
          
          for (l=0;l<nmatres;l++) q2[l]=0;
          if (rank==0) {
            double soff=0;
            for (i=0;i<newnr[nbolo];i++) soff+=hit2[0]*ix2[i];
            for (i=0;i<newnr[nbolo];i++) q2[i]=soff;

            if (Param->REMHDIP==1) {
              soff=0;
              for (i=newnr[nbolo];i<newnr[nbolo]+GAINSTEP*nbolo;i++) soff+=hit2[newnr[nbolo]]*ix2[i];
              for (i=newnr[nbolo];i<newnr[nbolo]+GAINSTEP*nbolo;i++) q2[i]=soff;
            }

            if (Param->flag_AVGR0==_PAR_TRUE) {
              soff=0;
              for (i=0;i<nbolo;i++) soff+=hit2[newnr[nbolo]+nbolo*(GAINSTEP2)]*
                                      ix2[newnr[nbolo]+nbolo*(GAINSTEP2)+i];
              for (i=0;i<nbolo;i++) q2[newnr[nbolo]+nbolo*(GAINSTEP2)+i]=soff; 
            }
            if (nmatco>0) {
              soff=0;
              for (i=0;i<nbolo;i++) soff+=hit2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)]
                                      *ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+i];
              for (i=0;i<nbolo;i++) q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+i]=soff;
            }
            if (nmatdust>0) {
              soff=0;
              for (i=0;i<nbolo;i++) 
                soff+=hit2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco]*
                  ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+i];
              for (i=0;i<nbolo;i++) q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+i]=soff;
            }
            if (nfreefree>0) {
              soff=0;
              for (i=0;i<nbolo;i++) 
                soff+=hit2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust]*
                  ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+i];
              for (i=0;i<nbolo;i++) q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+i]=soff;
            }
            
            if (DOFITTHETA==1) {
              soff=0;
              for (i=0;i<nbolo;i++) soff+=hit2[newnr[nbolo]+nbolo*(GAINSTEP2+npixbeam-1)]*
                                      ix2[newnr[nbolo]+nbolo*(GAINSTEP2+npixbeam-1)+i];
              for (i=0;i<nbolo;i++) q2[newnr[nbolo]+nbolo*(GAINSTEP2+npixbeam-1)+i]=soff;
            }
            
          }

          for (k=0;k<nnbpix;k++) if (flgpix[k]>0) {
        
            //long imat=the_stat_pix[k];
            long ndata = loc_nhpix[k];
            if (ndata>0) {
              hpix *htmp = loc_hpix[k];

              double vali=0,valq=0,valu=0;

#ifdef USEDII
              memset(dii ,0,sizeof(double)*nbolo*GAINSTEP);
              memset(dqq ,0,sizeof(double)*nbolo*GAINSTEP);
              memset(duu ,0,sizeof(double)*nbolo*GAINSTEP);

              for (l1=0;l1<ndata;l1++) {
                long ri1=htmp[l1].rg-globalBeginRing;
                double CO1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].co
                                                  -dpsisi[htmp[l1].ib]*htmp[l1].si);
                double SI1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].si
                                                  +dpsisi[htmp[l1].ib]*htmp[l1].co);
                if (flg_rg[htmp[l1].ib][ri1]!=0) {
                  
                  dii[htmp[l1].gi+htmp[l1].ib*GAINSTEP ]+=htmp[l1].w*htmp[l1].dip;
                  dqq[htmp[l1].gi+htmp[l1].ib*GAINSTEP ]+=htmp[l1].w*htmp[l1].dip*CO1;
                  duu[htmp[l1].gi+htmp[l1].ib*GAINSTEP ]+=htmp[l1].w*htmp[l1].dip*SI1;
                }
              }
              for (ib=0;ib<nbolo;ib++) {
                for (j=0;j<GAINSTEP;j++) {
                  solvemap(dii+j+ib*GAINSTEP,dqq+j+ib*GAINSTEP,duu+j+ib*GAINSTEP,II[k],IQ[k],IU[k],QQ[k],QU[k],UU[k]);
                }
              } 
#endif
              long l3;
              for (l3=0;l3<ndata;l3++) {
                long ri3=htmp[l3].rg-globalBeginRing;
                if (flg_rg[htmp[l3].ib][ri3]!=0) {
                  long ir3=rgord[htmp[l3].ib][ri3]+newnr[htmp[l3].ib];
                  vali+=htmp[l3].vi*ix2[ir3];
                  valq+=htmp[l3].vq*ix2[ir3];
                  valu+=htmp[l3].vu*ix2[ir3];
#ifndef USEDII
                  vali+=htmp[l3].lvi*ix2[newnr[nbolo]+htmp[l3].ib*GAINSTEP+htmp[l3].gi];
                  valq+=htmp[l3].lvq*ix2[newnr[nbolo]+htmp[l3].ib*GAINSTEP+htmp[l3].gi];
                  valu+=htmp[l3].lvu*ix2[newnr[nbolo]+htmp[l3].ib*GAINSTEP+htmp[l3].gi];
#endif
                }
              }
                
              for (ib=0;ib<nbolo;ib++) {


#ifdef USEDII
                for (j=0;j<GAINSTEP;j++) {
                  vali+=dii[j+ib*GAINSTEP]  *ix2[newnr[nbolo]+ib*GAINSTEP+j];
                  valq+=dqq[j+ib*GAINSTEP]  *ix2[newnr[nbolo]+ib*GAINSTEP+j];
                  valu+=duu[j+ib*GAINSTEP]  *ix2[newnr[nbolo]+ib*GAINSTEP+j];
                }
#endif

                if (nmatco>0) {
                  vali+=dcoi[ib+k*nbolo]  *ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+ib];
                  valq+=dcoq[ib+k*nbolo]  *ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+ib];
                  valu+=dcou[ib+k*nbolo]  *ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+ib];
                }            
                if (nmatdust>0) {
                  vali+=ddusti[ib+k*nbolo]  *ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+ib];
                  valq+=ddustq[ib+k*nbolo]  *ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+ib];
                  valu+=ddustu[ib+k*nbolo]  *ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+ib];
                }            
                if (nfreefree>0) {
                  vali+=dfri[ib+k*nbolo]  *ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+ib];
                  valq+=dfrq[ib+k*nbolo]  *ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+ib];
                  valu+=dfru[ib+k*nbolo]  *ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+ib];
                }            
                for (m=0;m<npixbeam;m++) {
                  vali+=ix2[newnr[nbolo]+nbolo*GAINSTEP2+ib+m*nbolo]*dpixi[ib+m*nbolo+k*nbolo*npixbeam];
                  valq+=ix2[newnr[nbolo]+nbolo*GAINSTEP2+ib+m*nbolo]*dpixq[ib+m*nbolo+k*nbolo*npixbeam];
                  valu+=ix2[newnr[nbolo]+nbolo*GAINSTEP2+ib+m*nbolo]*dpixu[ib+m*nbolo+k*nbolo*npixbeam];
                }   
              }

              double qri=0;
              double qrq=0;
              double qru=0;

              for (l1=0;l1<ndata;l1++) {
                long ri1=htmp[l1].rg-globalBeginRing;
                long iri1=rgord[htmp[l1].ib][ri1]+newnr[htmp[l1].ib];

                double CO1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].co
                                                  -dpsisi[htmp[l1].ib]*htmp[l1].si);
                double SI1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].si
                                                  +dpsisi[htmp[l1].ib]*htmp[l1].co);
                if (flg_rg[htmp[l1].ib][ri1]!=0) {
                  double ww=htmp[l1].wp;
                  double val2=ix2[iri1]-(vali+CO1*valq+SI1*valu);

                  val2+=ix2[newnr[nbolo]+htmp[l1].ib*GAINSTEP+htmp[l1].gi]*htmp[l1].dip;
                  if (nmatco>0)   
                    val2+=ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+htmp[l1].ib]*htmp[l1].comap;
                  if (nmatdust>0) 
                    val2+=ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+htmp[l1].ib]*htmp[l1].dustmap;
                  if (nfreefree>0) 
                    val2+=ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nfreefree+htmp[l1].ib]*htmp[l1].freefree;

                  for (m=0;m<npixbeam;m++) {
                    val2+=ix2[newnr[nbolo]+nbolo*GAINSTEP2+htmp[l1].ib+m*nbolo]*htmp[l1].listofpix[m];
                  }
                    
                  qri-=ww*val2;
                  qrq-=ww*val2*CO1;
                  qru-=ww*val2*SI1;
                    
                  q2[iri1]+=ww*val2;

                  q2[newnr[nbolo]+htmp[l1].ib*GAINSTEP+htmp[l1].gi]+=ww*val2*htmp[l1].dip;

                  ///////////// CO
                  if (nmatco>0) {
                    q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+htmp[l1].ib]+=ww*val2*htmp[l1].comap;
                  }
                    
                  ///////////// DUST
                  if (nmatdust>0) {
                    q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+htmp[l1].ib]+=ww*val2*htmp[l1].dustmap;
                  }

                  ///////////// FREEFREE
                  if (nfreefree>0) {
                    q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nfreefree+htmp[l1].ib]+=ww*val2*htmp[l1].freefree;
                  }
                    
                  for (j=0;j<npixbeam;j++) {
                    q2[newnr[nbolo]+nbolo*GAINSTEP2+j*nbolo+htmp[l1].ib]+=ww*val2*htmp[l1].listofpix[j];  
                  }
                    
                }
              }

              for (l2=0;l2<ndata;l2++) {
                long ri2=htmp[l2].rg-globalBeginRing;
                if (flg_rg[htmp[l2].ib][ri2]!=0) {
                  long ir=rgord[htmp[l2].ib][ri2]+newnr[htmp[l2].ib]; 
                  q2[ir]+=qri*htmp[l2].vi+qrq*htmp[l2].vq+qru*htmp[l2].vu;
#ifndef USEDII
                  q2[newnr[nbolo]+htmp[l2].ib*GAINSTEP+htmp[l2].gi]+=qri*htmp[l2].lvi+qrq*htmp[l2].lvq+qru*htmp[l2].lvu;
#endif
                }
              }
                
              for (ib=0;ib<nbolo;ib++) {


#ifdef USEDII
                for (j=0;j<GAINSTEP;j++) {
                  q2[newnr[nbolo]+GAINSTEP*ib+j]+=qri*dii[j+ib*GAINSTEP]+qrq*dqq[j+ib*GAINSTEP]+qru*duu[j+ib*GAINSTEP];
                }
#endif

                ///////////// CO
                if (nmatco>0) {
                  q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+ib]+=qri*dcoi[ib+k*nbolo]+qrq*dcoq[ib+k*nbolo]+qru*dcou[ib+k*nbolo];
                }
                ///////////// DUST
                if (nmatdust>0) {
                  q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+ib]+=qri*ddusti[ib+k*nbolo]+qrq*ddustq[ib+k*nbolo]+qru*ddustu[ib+k*nbolo];
                }

                ///////////// FREEFREE
                if (nfreefree>0) {
                  q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+ib]+=qri*dfri[ib+k*nbolo]+qrq*dfrq[ib+k*nbolo]+qru*dfru[ib+k*nbolo];
                }
                  
                for (j=0;j<npixbeam;j++) {
                  q2[newnr[nbolo]+nbolo*GAINSTEP2+j*nbolo+ib]+=qri*dpixi[ib+j*nbolo+k*nbolo*npixbeam]+qrq*dpixq[ib+j*nbolo+k*nbolo*npixbeam]+qru*dpixu[ib+j*nbolo+k*nbolo*npixbeam];  
                }
              }
            }
          }
          if (rank==0) {
            double *lb = (double *) malloc(sizeof(double)*(nmatres));
            for (rrk=1;rrk<mpi_size;rrk++) {
              MPI_Recv(lb,sizeof(double)*(nmatres), MPI_BYTE, rrk,1034, MPI_COMM_WORLD,&statu);
              for (l=0;l<nmatres;l++) q2[l]+=lb[l];
            }

            free(lb);
          }
          else {
            MPI_Send(q2, sizeof(double)*(nmatres), MPI_BYTE, 0, 1034, MPI_COMM_WORLD);
          }
          
          if (rank==0) for (i=0; i < nmatres; i++) {
            r2[i] = b2[i] - q2[i];
          }
          
        }
      else
        {
          if (rank==0) for (i=0; i < nmatres; i++) r2[i] -= alpha * q2[i];
        }
      
      if (rank==0) for (i=0; i < nmatres; i++) s2[i] = r2[i] / hit2[i];

      delta_old = delta_new;
      if (rank==0) {
        delta_new=0;
        for (i=0; i < nmatres ; i++) delta_new += r2[i] * s2[i];
        beta = delta_new / delta_old;
        for (i=0; i < nmatres ; i++) d2[i] = s2[i] + beta * d2[i];
      }
      iter ++;
      MPI_Bcast(&delta_new, sizeof(double), MPI_BYTE, 0, MPI_COMM_WORLD);
      
    }
  
  if (rank==0) fprintf (stderr,"\niter = %d - delta0 = %lg - delta_new = %lg\n", 
                        iter, delta0, delta_new);
  if (rank==0) fprintf (stderr,"CG in iter = %d (max=%d)\n", iter, itermax);
  
  if (rank==0) memcpy(ix2  ,x2old,nmatres*sizeof (double));
  MPI_Bcast(ix2, sizeof(double)*(nmatres), MPI_BYTE, 0, MPI_COMM_WORLD);
  itbogo++;
}

#ifdef UPDATE_DIP
#define UPDATE_DIP
#endif

void minimize_gain_nopol(double *ix2,double *gain,double *gaingi)
{
  MPI_Status statu;
  long i,rrk,j,k,l1,l2,ib;
  int itermax = 500;
  int iter;
  double  delta_new, delta_old, beta;
  double  alpha=1.0; // get rid of gcc "maybe-uninitialized" warning depending on optimsation level

  PIOLONG GAINSTEP2;
  int rank;
  int size;
  int mpi_size;
  int rank_size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank_size);
  rank=rank_size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  mpi_size=size;

  GAINSTEP2=GAINSTEP;

  nmatres=newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+nfreefree;

  if (rank==0) {
    fprintf(stderr,"==============================\n\nminimize_gain_nopol(): ITERATION NOPOL %ld %ld \n\n==============================\n",itbogo,(long) GAINSTEP);
    fprintf(stderr,"GAIN ");
    for (i=0;i<nbolo;i++) fprintf(stderr,"%lg ",gaingi[i*GAINSTEP]);
    fprintf(stderr,"\n");
    if (GAINSTEP>1) {
      fprintf(stderr,"...\nGAIN ");
      for (i=0;i<nbolo;i++) fprintf(stderr,"%lg ",gaingi[i*GAINSTEP+(GAINSTEP-1)]);
      fprintf(stderr,"\n");
    }
  }
  if (itbogo==0) delta0=0;
  MPI_Barrier(MPI_COMM_WORLD);
  
  struct timeval tp1,tp2;
  gettimeofday(&tp1,NULL);

  
  iter = 0;
  memset(b2  ,0,nmatres*sizeof (double));
  memset(d2  ,0,nmatres*sizeof (double));
  memset(q2  ,0,nmatres*sizeof (double));
  memset(r2  ,0,nmatres*sizeof (double));
  memset(s2  ,0,nmatres*sizeof (double));
  memset(hit2,0,nmatres*sizeof (double));

  //==========================================
  //=  Compute second member
  //=
  //=
  long l,m;

  ////// BUILD B2
  //GetProcMem(&vmem,&phymem);
  //if (rank==0) fprintf(stderr,"Rank: %ld Line=%d MEM %.1lf[%.1lf]MB\n",
  //                  (long) rank, __LINE__,
  //                  (double) vmem/1024./1024.,
  //                  (double) phymem/1024./1024.);


  for (k=0;k<nnbpix;k++)  {
    //long imat=the_stat_pix[k];
    long ndata = loc_nhpix[k];
    hpix *htmp = loc_hpix[k];

    II[k]=0;
      

    for (l1=0;l1<ndata;l1++) {
      long ri1=htmp[l1].rg-globalBeginRing;
      if (flg_rg[htmp[l1].ib][ri1]!=0) {
        II[k]+=htmp[l1].w;
      }
    }

    if (ndata>0&&flgpix[k]>0) {

      double SI=0;

#ifdef USEDII
      memset(dii ,0,sizeof(double)*nbolo*GAINSTEP);
#endif
      memset(dcoi+k*nbolo,0,sizeof(double)*nbolo);
      memset(dfri+k*nbolo,0,sizeof(double)*nbolo);
      memset(dthetai+k*nbolo,0,sizeof(double)*nbolo);
      memset(ddusti+k*nbolo,0,sizeof(double)*nbolo);
      memset(dpixi+k*nbolo*npixbeam,0,sizeof(double)*nbolo*npixbeam);
      

      for (l1=0;l1<ndata;l1++) {
        long ri1=htmp[l1].rg-globalBeginRing;
        
        double g1=gaingi[htmp[l1].gi+htmp[l1].ib*GAINSTEP];
        
        if (flg_rg[htmp[l1].ib][ri1]!=0) {
          
          htmp[l1].wp=0;
          //htmp[l1].thsig=0;
          
          if (REMHDIP==0) SI+=htmp[l1].w*(htmp[l1].sig*g1-htmp[l1].fsl-htmp[l1].dip-gain[htmp[l1].adu+htmp[l1].ib*GAINSTEPADU]);
          else SI+=htmp[l1].w*(htmp[l1].sig*g1-htmp[l1].fsl-htmp[l1].freefree-gain[htmp[l1].adu+htmp[l1].ib*GAINSTEPADU]);
          
#ifdef USEDII
          dii[htmp[l1].gi+htmp[l1].ib*GAINSTEP ]+=htmp[l1].w*htmp[l1].dip;
          
#endif
            
          if (nmatco>0) {
            dcoi[htmp[l1].ib+k*nbolo] += htmp[l1].w*htmp[l1].comap;
          }
          if (nmatdust>0) {
            ddusti[htmp[l1].ib+k*nbolo] += htmp[l1].w*htmp[l1].dustmap;
          }
          if (nfreefree>0) {
            dfri[htmp[l1].ib+k*nbolo] += htmp[l1].w*htmp[l1].freefree;
          }
          
          if (ittt>0) {
            for (m=0;m<npixbeam;m++)  {
              dpixi[nbolo*m+htmp[l1].ib+k*nbolo*npixbeam] += htmp[l1].w*htmp[l1].listofpix[m];
            }
          } 
        }
      }
      SI/=II[k];
      SSI[k]=SI;
      
      for (ib=0;ib<nbolo;ib++) {
        
#ifdef USEDII
        for (j=0;j<GAINSTEP;j++) {
          dii[j+ib*GAINSTEP]=dii[j+ib*GAINSTEP]/II[k];
        }
#endif

        if (nmatco>0) {
          dcoi[ib+k*nbolo]=dcoi[ib+k*nbolo]/II[k];
        }
        if (nmatdust>0) {
          ddusti[ib+k*nbolo]=ddusti[ib+k*nbolo]/II[k];
        }
        if (nfreefree>0) {
          dfri[ib+k*nbolo]=dfri[ib+k*nbolo]/II[k];
        }
        
        if (ittt>0) {
          for (m=0;m<npixbeam;m++)  {
            dpixi[ib+m*nbolo+k*nbolo*npixbeam]=dpixi[ib+m*nbolo+k*nbolo*npixbeam]/II[k];
          }
        } 
      }

      for (l1=0;l1<ndata;l1++) {
        long ri1=htmp[l1].rg-globalBeginRing;
        if (flg_rg[htmp[l1].ib][ri1]!=0) {
          htmp[l1].vi=htmp[l1].w/II[k];
#ifndef USEDII
          htmp[l1].lvi=htmp[l1].w*htmp[l1].dip/II[k];
#endif
        }
      }
      for (l1=0;l1<ndata;l1++) {
        long ri1=htmp[l1].rg-globalBeginRing;
        
        double g1=gaingi[htmp[l1].gi+htmp[l1].ib*GAINSTEP];
        
        if (flg_rg[htmp[l1].ib][ri1]!=0) {
          double divi=NEP_tab[htmp[l1].ib]*htmp[l1].hit*g1;
          divi=divi*divi;
          
          
          if (divi==0) htmp[l1].wp=0;
          else htmp[l1].wp=1/divi;
          
          if (itbogo==0) normaoff+=NEP_tab[htmp[l1].ib]*htmp[l1].hit;
          
        }
      }


      
      for (l1=0;l1<ndata;l1++) {
        long ri1=htmp[l1].rg-globalBeginRing;
        long iri1=rgord[htmp[l1].ib][ri1]+newnr[htmp[l1].ib];
        
        double g1=gaingi[htmp[l1].gi+htmp[l1].ib*GAINSTEP];
        
        if (flg_rg[htmp[l1].ib][ri1]!=0) {
          double ww=htmp[l1].wp;
          double tmp=((htmp[l1].sig*g1-htmp[l1].fsl-gain[htmp[l1].adu+htmp[l1].ib*GAINSTEPADU])-(SI));
          if (REMHDIP==0) tmp-=htmp[l1].dip;
          else tmp-=htmp[l1].freefree;
            
          long l3;
#ifndef USEDII
          memset(cdip,0,GAINSTEP*nbolo*sizeof(double));
#endif
          
          for (l3=0;l3<ndata;l3++) {
            long ri3=htmp[l3].rg-globalBeginRing;
            if (flg_rg[htmp[l3].ib][ri3]!=0) {
              long ir3=rgord[htmp[l3].ib][ri3]+newnr[htmp[l3].ib];
              ctmp[ir3]=-htmp[l3].vi;
#ifndef USEDII
              cdip[htmp[l3].ib*GAINSTEP+htmp[l3].gi]-=htmp[l3].lvi;
#endif
            }
          }
          ctmp[iri1]+=1;
             
#ifdef USEDII
          for (ib=0;ib<nbolo;ib++) 
            for (j=0;j<GAINSTEP;j++) {
              cdip[ib*GAINSTEP+j]=-dii[j+ib*GAINSTEP];
            }
#endif
          cdip[htmp[l1].ib*GAINSTEP+htmp[l1].gi]+=htmp[l1].dip;

          for (ib=0;ib<nbolo;ib++) cco[ib]=-dcoi[ib+k*nbolo];
          cco[htmp[l1].ib]+=htmp[l1].comap;
          
          for (ib=0;ib<nbolo;ib++) cdust[ib]=-ddusti[ib+k*nbolo];
          cdust[htmp[l1].ib]+=htmp[l1].dustmap;
          
          for (ib=0;ib<nbolo;ib++) ccfree[ib]=-dfri[ib+k*nbolo];
          ccfree[htmp[l1].ib]+=htmp[l1].freefree;
          
          for (m=0;m<npixbeam;m++) {
            for (ib=0;ib<nbolo;ib++) cpix[ib+m*nbolo]=-dpixi[ib+m*nbolo+k*nbolo*npixbeam];
            cpix[htmp[l1].ib+m*nbolo]+=htmp[l1].listofpix[m];
          }

          long ir;
          /////////////////  OFFSET
          
          for (l2=0;l2<ndata;l2++) {
            long ri2=htmp[l2].rg-globalBeginRing;
            if (flg_rg[htmp[l2].ib][ri2]!=0) {
              long ir=rgord[htmp[l2].ib][ri2]+newnr[htmp[l2].ib]; 
              b2[ir]+=ww*tmp*ctmp[ir];
              hit2[ir]+=ww*ctmp[ir]*ctmp[ir];
#ifndef USEDII
              b2[newnr[nbolo]+htmp[l2].ib*GAINSTEP+htmp[l2].gi]+=ww*tmp*cdip[htmp[l2].ib*GAINSTEP+htmp[l2].gi];
              hit2[newnr[nbolo]+htmp[l2].ib*GAINSTEP+htmp[l2].gi]+=ww*cdip[htmp[l2].ib*GAINSTEP+htmp[l2].gi]
                *cdip[htmp[l2].ib*GAINSTEP+htmp[l2].gi];
#endif
            }
          }
          
          
          /////////////////  DIPOLE FIT
          
          for (ir=0;ir<nbolo;ir++) {
            
#ifdef USEDII
            
            for (j=0;j<GAINSTEP;j++) {
              b2[newnr[nbolo]+ir*GAINSTEP+j]+=ww*tmp*cdip[ir*GAINSTEP+j];
              hit2[newnr[nbolo]+ir*GAINSTEP+j]+=ww*cdip[ir*GAINSTEP+j]*cdip[ir*GAINSTEP+j];
            }
            
#endif
              
            ///////////// CO
            if (nmatco>0) {
              b2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+ir]+=ww*tmp*cco[ir];
              hit2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+ir]+=ww*cco[ir]*cco[ir];
            }
            
            ///////////// DUST
            
            if (nmatdust>0) {
              b2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+ir]+=ww*tmp*cdust[ir];
              hit2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+ir]+=ww*cdust[ir]*cdust[ir];
            }
            
            ///////////// FREEFREE
            
            if (nfreefree>0) {
              b2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+ir]+=ww*tmp*ccfree[ir];
              hit2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+ir]+=ww*ccfree[ir]*ccfree[ir];
            }
            
            ////////// PIXBEAM
            for (j=0;j<npixbeam;j++) {
              b2[newnr[nbolo]+nbolo*GAINSTEP2+j*nbolo+ir]+=ww*tmp*cpix[ir+j*nbolo];
              hit2[newnr[nbolo]+nbolo*GAINSTEP2+j*nbolo+ir]+=ww*cpix[ir+j*nbolo]*cpix[ir+j*nbolo];
            }
          }
        }
      }
    }
  }

  if (rank==0) {
    double *lb = (double *) malloc(sizeof(double)*(nmatres));
    for (rrk=1;rrk<mpi_size;rrk++) {
      MPI_Recv(lb,sizeof(double)*(nmatres), MPI_BYTE, rrk,1031, MPI_COMM_WORLD,&statu);
      for (l=0;l<nmatres;l++) b2[l]+=lb[l];
    }
    free(lb);
  }
  else MPI_Send(b2, sizeof(double)*(nmatres), MPI_BYTE, 0, 1031, MPI_COMM_WORLD);
  
  if (rank==0) {
    double *lb = (double *) malloc(sizeof(double)*(nmatres));
    for (rrk=1;rrk<mpi_size;rrk++) {
      MPI_Recv(lb,sizeof(double)*(nmatres), MPI_BYTE, rrk,1033, MPI_COMM_WORLD,&statu);

      for (l=0;l<nmatres;l++) hit2[l]+=lb[l];
    }
    free(lb);
  }
  else {
    MPI_Send(hit2, sizeof(double)*(nmatres), MPI_BYTE, 0, 1033, MPI_COMM_WORLD);
  }
  

  //==========================================================
  // Compute Ax
  //healspline_fill_s_with_PtP_s (hs_rec, Ndata, vec, x, q);
  // + 
  // Preconditionnement
  //
  //
  
  for (l=0;l<nmatres;l++) q2[l]=0;
  if (rank==0) {

    double soff=0;
    for (i=0;i<newnr[nbolo];i++) soff+=hit2[0]*ix2[i];
    for (i=0;i<newnr[nbolo];i++) q2[i]=soff;
    if (REMHDIP==1) {
      soff=0;
      for (i=newnr[nbolo];i<newnr[nbolo]+GAINSTEP*nbolo;i++) soff+=hit2[newnr[nbolo]]*ix2[i];
      for (i=newnr[nbolo];i<newnr[nbolo]+GAINSTEP*nbolo;i++) q2[i]=soff;
    }
    if (Param->flag_AVGR0==_PAR_TRUE) {
      soff=0;
      for (i=0;i<nbolo;i++) soff+=hit2[newnr[nbolo]+nbolo*(GAINSTEP2)]*
                              ix2[newnr[nbolo]+nbolo*(GAINSTEP2)+i];
      for (i=0;i<nbolo;i++) q2[newnr[nbolo]+nbolo*(GAINSTEP2)+i]=soff;
      for (i=0;i<nbolo;i++) b2[newnr[nbolo]+nbolo*(GAINSTEP2)+i]+=
                              hit2[newnr[nbolo]+nbolo*(GAINSTEP2)]*Param->AVGR0;
      
    }

    if (nmatco>0) {
      soff=0;
      for (i=0;i<nbolo;i++) soff+=hit2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)]
                              *ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+i];
      for (i=0;i<nbolo;i++) q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+i]=soff;
      for (i=0;i<nbolo;i++) b2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+i]+=
                              hit2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)]*Param->AVGCO;
    }
    if (nmatdust>0) {
      soff=0;
      for (i=0;i<nbolo;i++) 
        soff+=hit2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco]*
          ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+i];
      for (i=0;i<nbolo;i++) q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+i]=soff;
      for (i=0;i<nbolo;i++) b2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+i]+=
                              hit2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco]*Param->AVGDUST;
    }
    if (nfreefree>0) {
      soff=0;
      for (i=0;i<nbolo;i++) 
        soff+=hit2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust]*
          ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+i];
      for (i=0;i<nbolo;i++) q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+i]=soff;
    }
    
  }

  for (k=0;k<nnbpix;k++) if (flgpix[k]>0) {
    
    long ndata = loc_nhpix[k];
    if (ndata>0) {
#ifdef TIMING
      gettimeofday(&tp1,NULL);
#endif
      hpix *htmp = loc_hpix[k];

      double vali=0;

#ifdef USEDII
      memset(dii ,0,sizeof(double)*nbolo*GAINSTEP);

      for (l1=0;l1<ndata;l1++) {
        long ri1=htmp[l1].rg-globalBeginRing;

        if (flg_rg[htmp[l1].ib][ri1]!=0) {
          
          dii[htmp[l1].gi+htmp[l1].ib*GAINSTEP ]+=htmp[l1].w*htmp[l1].dip;
          
        }
      }
      for (ib=0;ib<nbolo;ib++) {
        for (j=0;j<GAINSTEP;j++) {
          dii[j+ib*GAINSTEP]=dii[j+ib*GAINSTEP]/II[k];
        }
      } 
#endif
        
      long l3;
      for (l3=0;l3<ndata;l3++) {
        long ri3=htmp[l3].rg-globalBeginRing;
        if (flg_rg[htmp[l3].ib][ri3]!=0) {
          long ir3=rgord[htmp[l3].ib][ri3]+newnr[htmp[l3].ib];
          vali+=htmp[l3].vi*ix2[ir3];
#ifndef USEDII
          vali+=htmp[l3].lvi*ix2[newnr[nbolo]+htmp[l3].ib*GAINSTEP+htmp[l3].gi];
#endif
        }
      }
      //fprintf(stderr,"I0 vali %lg %lg %lg\n",vali,valq,valu);
      
      for (ib=0;ib<nbolo;ib++) {
#ifdef USEDII
        for (j=0;j<GAINSTEP;j++) {
          vali+=dii[j+ib*GAINSTEP]  *ix2[newnr[nbolo]+ib*GAINSTEP+j];
        }
        //fprintf(stderr,"I1 vali %lg %lg %lg\n",vali,valq,valu);
#endif

        if (nmatco>0) {
          vali+=dcoi[ib+k*nbolo]  *ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+ib];
        }            
        if (nmatdust>0) {
          vali+=ddusti[ib+k*nbolo]  *ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+ib];
        }            
        if (nfreefree>0) {
          vali+=dfri[ib+k*nbolo]  *ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+ib];
        }
        //fprintf(stderr,"I3 vali %lg %lg %lg : %lg %lg %lg %lg\n",vali,valq,valu,
        //dfri[ib+k*nbolo],dfrq[ib+k*nbolo],dfru[ib+k*nbolo],
        //ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+ib]);            
        
        for (m=0;m<npixbeam;m++) {
          vali+=ix2[newnr[nbolo]+nbolo*GAINSTEP2+ib+m*nbolo]*dpixi[ib+m*nbolo+k*nbolo*npixbeam];
        }  
        
        //fprintf(stderr,"I4 vali %lg %lg %lg\n",vali,valq,valu);
      }
      double qri=0;
      
      for (l1=0;l1<ndata;l1++) {
        long ri1=htmp[l1].rg-globalBeginRing;
        long iri1=rgord[htmp[l1].ib][ri1]+newnr[htmp[l1].ib];

        if (flg_rg[htmp[l1].ib][ri1]!=0) {
          double ww=htmp[l1].wp;
          double val2=ix2[iri1]-(vali);
          
          val2+=ix2[newnr[nbolo]+htmp[l1].ib*GAINSTEP+htmp[l1].gi]*htmp[l1].dip;

          //fprintf(stderr,"I2 val2 %lg\n",val2);
          if (nmatco>0)   
            val2+=ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+htmp[l1].ib]*htmp[l1].comap;
          if (nmatdust>0) 
            val2+=ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+htmp[l1].ib]*htmp[l1].dustmap;
          if (nfreefree>0) 
            val2+=ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+htmp[l1].ib]*htmp[l1].freefree;
          
          //fprintf(stderr,"I3 val2 %lg\n",val2);
          for (m=0;m<npixbeam;m++) {
            val2+=ix2[newnr[nbolo]+nbolo*GAINSTEP2+htmp[l1].ib+m*nbolo]*htmp[l1].listofpix[m];
          }
          
          
          qri-=ww*val2;
          
          q2[iri1]+=ww*val2;
          
          q2[newnr[nbolo]+htmp[l1].ib*GAINSTEP+htmp[l1].gi]+=ww*val2*htmp[l1].dip;

          ///////////// CO
          if (nmatco>0) {
            q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+htmp[l1].ib]+=ww*val2*htmp[l1].comap;
          }
          
          ///////////// DUST
          if (nmatdust>0) {
            q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+htmp[l1].ib]+=ww*val2*htmp[l1].dustmap;
          }
          
          ///////////// FREEFREE
          if (nfreefree>0) {
            q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+htmp[l1].ib]+=ww*val2*htmp[l1].freefree;
          }
          
          ////////// SYSTE
          for (j=0;j<npixbeam;j++) {
            q2[newnr[nbolo]+nbolo*GAINSTEP2+j*nbolo+htmp[l1].ib]+=ww*val2*htmp[l1].listofpix[j];  
          }
          
          
        }
      }
    

      for (l2=0;l2<ndata;l2++) {
        long ri2=htmp[l2].rg-globalBeginRing;
        if (flg_rg[htmp[l2].ib][ri2]!=0) {
          long ir=rgord[htmp[l2].ib][ri2]+newnr[htmp[l2].ib]; 
          q2[ir]+=qri*htmp[l2].vi;
#ifndef USEDII
          q2[newnr[nbolo]+htmp[l2].ib*GAINSTEP+htmp[l2].gi]+=qri*htmp[l2].lvi;
#endif
        }
      }
    
        
      for (ib=0;ib<nbolo;ib++) {
        
#ifdef USEDII
        for (j=0;j<GAINSTEP;j++) {
          q2[newnr[nbolo]+GAINSTEP*ib+j]+=qri*dii[j+ib*GAINSTEP];
        }
#endif

        ///////////// CO
        if (nmatco>0) {
          q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+ib]+=qri*dcoi[ib+k*nbolo];
        }
        ///////////// DUST
        if (nmatdust>0) {
          q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+ib]+=qri*ddusti[ib+k*nbolo];
        }
        ///////////// FREEFREE
        if (nfreefree>0) {
          q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+ib]+=qri*dfri[ib+k*nbolo];
        }
        
        for (j=0;j<npixbeam;j++) {
          q2[newnr[nbolo]+nbolo*GAINSTEP2+j*nbolo+ib]+=qri*dpixi[ib+j*nbolo+k*nbolo*npixbeam];  
        }
      }
    }
      
#ifdef TIMING
      gettimeofday(&tp2,NULL);
      double dt=(double)(tp2.tv_sec-tp1.tv_sec)+(1E-6)*(tp2.tv_usec-tp1.tv_usec);
      dthit[ndata]+=dt;
      ndthit[ndata]+=1;
#endif
  }

  if (rank==0) {
    double *lb = (double *) malloc(sizeof(double)*(nmatres));
    for (rrk=1;rrk<mpi_size;rrk++) {
      MPI_Recv(lb,sizeof(double)*(nmatres), MPI_BYTE, rrk,1032, MPI_COMM_WORLD,&statu);
      for (l=0;l<nmatres;l++) q2[l]+=lb[l];
    }
    free(lb);
  }
  else {
    MPI_Send(q2, sizeof(double)*(nmatres), MPI_BYTE, 0, 1032, MPI_COMM_WORLD);
  }

  if (rank==0) fprintf(stderr,"QQ2 %lg\n",q2[0]);

  if (rank==0) {
    for (i=0; i < nmatres; i++)
      {
        r2[i] = b2[i] - q2[i];
        d2[i] = r2[i] / hit2[i];
      }
  }

  double delta_new_tmp = 0.0;
  if (rank==0) {
    for (i = 0; i < nmatres; i++) {
      delta_new_tmp += b2[i] ;
      if (isnan(b2[i])) {
        fprintf(stderr,"NAN B2 PBS %ld\n",(long) i);
      }
      if (isnan(d2[i])) {
        fprintf(stderr,"NAN D2 PBS %ld %lg\n",(long) i,hit2[i]);
      }
    }
    //fprintf(stderr,"B2 %lg\n",delta_new_tmp);
  
    delta_new_tmp = 0.0;
    for (i = 0; i < nmatres; i++) {
      delta_new_tmp += q2[i] ;
      if (isnan(q2[i])) {
        fprintf(stderr,"NAN Q2 PBS %ld\n",(long) i);
      }
    }
    //fprintf(stderr,"Q2 %lg\n",delta_new_tmp);
    delta_new_tmp = 0.0;
    for (i = 0; i < nmatres; i++) {
      delta_new_tmp += q2[i]-b2[i] ;
      //fprintf(stderr,"B2 Q2 B2-Q2 [%ld]: %lg\t%lg\t%lg\n",(long) i,b2[i],q2[i],q2[i]-b2[i]);
    }
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
  delta_new_tmp = 0.0;
  if (rank==0) for (i=0; i < nmatres; i++) {
    delta_new_tmp += r2[i] * d2[i];
  }
  
  delta_new=0;
  for (i=0;i<mpi_size;i++) {
    double tmp=delta_new_tmp;
    MPI_Bcast(&tmp, sizeof(double), MPI_BYTE, i, MPI_COMM_WORLD);
    delta_new+=tmp;
  }

  if (itbogo==0) delta0 = delta_new;
  if (rank==0) fprintf (stderr, "iter = %d - delta0 = %lg - delta_new = %lg\n", iter, delta0, delta_new);
  int testwrit=0;

  if (itbogo==0) MPI_Bcast(&normaoff, sizeof(double), MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&delta0, sizeof(double), MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&delta_new, sizeof(double), MPI_BYTE, 0, MPI_COMM_WORLD);


  while ((iter < itermax)  && ((delta_new) > delta0*1E-24) && ((delta_new) > 1E-20)) //Param->XI2STOP))
    {
      // q <= Ad
      //if (rank==0&&mindelta>delta_new) {
      if (rank==0) {
        memcpy(x2old,ix2,nmatres*sizeof(double));
        testwrit=1;
      }
    //  else testwrit=0;

      //healspline_fill_s_with_PtP_s (hs_rec, Ndata, vec, d, q);
      //for (i=0;i<mpi_size;i++) {
      //MPI_Bcast(d+tab_begr[i]*2, sizeof(double)*(tab_edr[i]-tab_begr[i]+1)*2, MPI_BYTE, i, MPI_COMM_WORLD);
      //}
      MPI_Bcast(d2, sizeof(double)*(nmatres), MPI_BYTE, 0, MPI_COMM_WORLD);
      
      // ===========================================================================
      // PROJECTION DE d dans q
      //
      
      for (l=0;l<nmatres;l++) q2[l]=0;
      if (rank==0) {
        double soff=0;
        for (i=0;i<newnr[nbolo];i++) soff+=hit2[0]*d2[i];
        for (i=0;i<newnr[nbolo];i++) q2[i]=soff;
        if (REMHDIP==1) {
          soff=0;
          for (i=newnr[nbolo];i<newnr[nbolo]+GAINSTEP*nbolo;i++) soff+=hit2[newnr[nbolo]]*d2[i];
          for (i=newnr[nbolo];i<newnr[nbolo]+GAINSTEP*nbolo;i++) q2[i]=soff;
        }
        
        if (Param->flag_AVGR0==_PAR_TRUE) {
          soff=0;
          for (i=0;i<nbolo;i++) soff+=hit2[newnr[nbolo]+nbolo*(GAINSTEP2)]*
                                  d2[newnr[nbolo]+nbolo*(GAINSTEP2)+i];
          for (i=0;i<nbolo;i++) q2[newnr[nbolo]+nbolo*(GAINSTEP2)+i]=soff; 
        }
        if (nmatco>0) {
          soff=0;
          for (i=0;i<nbolo;i++) soff+=hit2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)]
                                  *d2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+i];
          for (i=0;i<nbolo;i++) q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+i]=soff;
        }
        if (nfreefree>0) {
          soff=0;
          for (i=0;i<nbolo;i++) soff+=hit2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust]
                                  *d2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+i];
          for (i=0;i<nbolo;i++) q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+i+nmatco+nmatdust]=soff;
        }
        if (nmatdust>0) {
          soff=0;
          for (i=0;i<nbolo;i++) 
            soff+=hit2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco]*
              d2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+i];
          for (i=0;i<nbolo;i++) q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+i]=soff;
        }
      }
      for (k=0;k<nnbpix;k++) if (flgpix[k]>0) {
        
        //long imat=the_stat_pix[k];
        long ndata = loc_nhpix[k];
        if (ndata>0) {
          hpix *htmp = loc_hpix[k];

          double vali=0;

#ifdef USEDII
          memset(dii ,0,sizeof(double)*nbolo*GAINSTEP);

          for (l1=0;l1<ndata;l1++) {
            long ri1=htmp[l1].rg-globalBeginRing;

            if (flg_rg[htmp[l1].ib][ri1]!=0) {
              
              dii[htmp[l1].gi+htmp[l1].ib*GAINSTEP ]+=htmp[l1].w*htmp[l1].dip;
            }
          }
          for (ib=0;ib<nbolo;ib++) {
            for (j=0;j<GAINSTEP;j++) {
              dii[j+ib*GAINSTEP]=dii[j+ib*GAINSTEP]/II[k];
            }
          } 
#endif
          long l3;
          for (l3=0;l3<ndata;l3++) {
            long ri3=htmp[l3].rg-globalBeginRing;
            if (flg_rg[htmp[l3].ib][ri3]!=0) {
              long ir3=rgord[htmp[l3].ib][ri3]+newnr[htmp[l3].ib];
              vali+=htmp[l3].vi*d2[ir3];
#ifndef USEDII
              vali+=htmp[l3].lvi*d2[newnr[nbolo]+htmp[l3].ib*GAINSTEP+htmp[l3].gi];
#endif
            }
          }
          for (ib=0;ib<nbolo;ib++) {
            
            
#ifdef USEDII
            for (j=0;j<GAINSTEP;j++) {
              vali+=dii[j+ib*GAINSTEP]  *d2[newnr[nbolo]+ib*GAINSTEP+j];
            }
#endif

            if (nmatco>0) {
              vali+=dcoi[ib+k*nbolo]  *d2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+ib];
            }            
            if (nmatdust>0) {
              vali+=ddusti[ib+k*nbolo]  *d2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+ib];
            }            
            if (nfreefree>0) {
              vali+=dfri[ib+k*nbolo]  *d2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+ib];
            }            
            for (m=0;m<npixbeam;m++) {
              vali+=d2[newnr[nbolo]+nbolo*GAINSTEP2+ib+m*nbolo]*dpixi[ib+m*nbolo+k*nbolo*npixbeam];
            }   
          }

          double qri=0;

          for (l1=0;l1<ndata;l1++) {
            long ri1=htmp[l1].rg-globalBeginRing;
            long iri1=rgord[htmp[l1].ib][ri1]+newnr[htmp[l1].ib];

            if (flg_rg[htmp[l1].ib][ri1]!=0) {
              double ww=htmp[l1].wp;
              double val2=d2[iri1]-(vali);

              val2+=d2[newnr[nbolo]+htmp[l1].ib*GAINSTEP+htmp[l1].gi]*htmp[l1].dip;
              if (nmatco>0)   
                val2+=d2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+htmp[l1].ib]*htmp[l1].comap;
              if (nmatdust>0) 
                val2+=d2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+htmp[l1].ib]*htmp[l1].dustmap;
              if (nfreefree>0) 
                val2+=d2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+htmp[l1].ib]*htmp[l1].freefree;
                
              for (m=0;m<npixbeam;m++) {
                val2+=d2[newnr[nbolo]+nbolo*GAINSTEP2+htmp[l1].ib+m*nbolo]*htmp[l1].listofpix[m];
              }
              
                
              qri-=ww*val2;
                
              q2[iri1]+=ww*val2;

              q2[newnr[nbolo]+htmp[l1].ib*GAINSTEP+htmp[l1].gi]+=ww*val2*htmp[l1].dip;

              ///////////// CO
              if (nmatco>0) {
                q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+htmp[l1].ib]+=ww*val2*htmp[l1].comap;
              }
                
              ///////////// DUST
              if (nmatdust>0) {
                q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+htmp[l1].ib]+=ww*val2*htmp[l1].dustmap;
              }
                
              ///////////// FREEFREE
              if (nfreefree>0) {
                q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+htmp[l1].ib]+=ww*val2*htmp[l1].freefree;
              }
            
                
              for (j=0;j<npixbeam;j++) {
                q2[newnr[nbolo]+nbolo*GAINSTEP2+j*nbolo+htmp[l1].ib]+=ww*val2*htmp[l1].listofpix[j];  
              }
                
            }
          }


          for (l2=0;l2<ndata;l2++) {
            long ri2=htmp[l2].rg-globalBeginRing;
            if (flg_rg[htmp[l2].ib][ri2]!=0) {
              long ir=rgord[htmp[l2].ib][ri2]+newnr[htmp[l2].ib]; 
              q2[ir]+=qri*htmp[l2].vi;
#ifndef USEDII
              q2[newnr[nbolo]+htmp[l2].ib*GAINSTEP+htmp[l2].gi]+=qri*htmp[l2].lvi;
#endif
            }
          }
            
          for (ib=0;ib<nbolo;ib++) {


#ifdef USEDII
            for (j=0;j<GAINSTEP;j++) {
              q2[newnr[nbolo]+GAINSTEP*ib+j]+=qri*dii[j+ib*GAINSTEP];
            }
#endif

            ///////////// CO
            if (nmatco>0) {
              q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+ib]+=qri*dcoi[ib+k*nbolo];
            }
            ///////////// DUST
            if (nmatdust>0) {
              q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+ib]+=qri*ddusti[ib+k*nbolo];
            }
            ///////////// FREEFREE
            if (nfreefree>0) {
              q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+ib]+=qri*dfri[ib+k*nbolo];
            }

            for (j=0;j<npixbeam;j++) {
              q2[newnr[nbolo]+nbolo*GAINSTEP2+j*nbolo+ib]+=qri*dpixi[ib+j*nbolo+k*nbolo*npixbeam];  
            }
          }
        }
      }
      if (rank==0) {
        double *lb = (double *) malloc(sizeof(double)*(nmatres));
        for (rrk=1;rrk<mpi_size;rrk++) {
          MPI_Recv(lb,sizeof(double)*(nmatres), MPI_BYTE, rrk,1034, MPI_COMM_WORLD,&statu);
          for (l=0;l<nmatres;l++) q2[l]+=lb[l];
        }
        free(lb);
      }
      else {
        MPI_Send(q2, sizeof(double)*(nmatres), MPI_BYTE, 0, 1034, MPI_COMM_WORLD);
      }
      double dtq = 0.0;
      if (rank==0) {
        for (i=0; i < nmatres; i++) dtq += d2[i] * q2[i];
        alpha = delta_new / dtq;
        for (i=0; i < nmatres ; i++) ix2[i] += alpha * d2[i];
      }


      gettimeofday(&tp2,NULL);
      if (rank==0) fprintf (stderr,"iter = %ld-%d - delta_new = %lg %ld %3lfs\n", itbogo, iter, delta_new,
                          (long) testwrit,(double)(tp2.tv_sec-tp1.tv_sec)+(1E-6)*(tp2.tv_usec-tp1.tv_usec));

      //if (rank==0) fprintf (stderr,".");
      gettimeofday(&tp1,NULL);
          
      if (iter % 100 == 0 && iter !=0)
        {
          // Use the best case
          memcpy(ix2,x2old,nmatres*sizeof(double));
          //for (i=0;i<mpi_size;i++) {
          //  MPI_Bcast(x+tab_begr[i]*2, sizeof(double)*(tab_edr[i]-tab_begr[i]+1)*2, MPI_BYTE, i, MPI_COMM_WORLD);
          //}
          MPI_Bcast(ix2, sizeof(double)*(nmatres), MPI_BYTE, 0, MPI_COMM_WORLD);

          
      // ===========================================================================
      // PROJECTION DE x dans q
      //
          
          for (l=0;l<nmatres;l++) q2[l]=0;
          if (rank==0) {
            double soff=0;
            for (i=0;i<newnr[nbolo];i++) soff+=hit2[0]*ix2[i];
            for (i=0;i<newnr[nbolo];i++) q2[i]=soff;
            if (REMHDIP==1) {
              soff=0;
              for (i=newnr[nbolo];i<newnr[nbolo]+GAINSTEP*nbolo;i++) soff+=hit2[newnr[nbolo]]*ix2[i];
              for (i=newnr[nbolo];i<newnr[nbolo]+GAINSTEP*nbolo;i++) q2[i]=soff;
            }
            

            if (Param->flag_AVGR0==_PAR_TRUE) {
              soff=0;
              for (i=0;i<nbolo;i++) soff+=hit2[newnr[nbolo]+nbolo*(GAINSTEP2)]*
                                      ix2[newnr[nbolo]+nbolo*(GAINSTEP2)+i];
              for (i=0;i<nbolo;i++) q2[newnr[nbolo]+nbolo*(GAINSTEP2)+i]=soff; 
            }
            if (nmatco>0) {
              soff=0;
              for (i=0;i<nbolo;i++) soff+=hit2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)]
                                      *ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+i];
              for (i=0;i<nbolo;i++) q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+i]=soff;
            }
            if (nmatdust>0) {
              soff=0;
              for (i=0;i<nbolo;i++) 
                soff+=hit2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco]*
                  ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+i];
              for (i=0;i<nbolo;i++) q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+i]=soff;
            }
            if (nfreefree>0) {
              soff=0;
              for (i=0;i<nbolo;i++) 
                soff+=hit2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust]*
                  ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+i];
              for (i=0;i<nbolo;i++) q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+i]=soff;
            }
            
            
          }

          for (k=0;k<nnbpix;k++) if (flgpix[k]>0) {
        
            //long imat=the_stat_pix[k];
            long ndata = loc_nhpix[k];
            if (ndata>0) {
              hpix *htmp = loc_hpix[k];

              double vali=0;

#ifdef USEDII
              memset(dii ,0,sizeof(double)*nbolo*GAINSTEP);


              for (l1=0;l1<ndata;l1++) {
                long ri1=htmp[l1].rg-globalBeginRing;

                if (flg_rg[htmp[l1].ib][ri1]!=0) {
                  
                  dii[htmp[l1].gi+htmp[l1].ib*GAINSTEP ]+=htmp[l1].w*htmp[l1].dip;
                }
              }
              for (ib=0;ib<nbolo;ib++) {
                for (j=0;j<GAINSTEP;j++) {
                  dii[j+ib*GAINSTEP]=dii[j+ib*GAINSTEP]/II[k];
                }
              } 
#endif
              long l3;
              for (l3=0;l3<ndata;l3++) {
                long ri3=htmp[l3].rg-globalBeginRing;
                if (flg_rg[htmp[l3].ib][ri3]!=0) {
                  long ir3=rgord[htmp[l3].ib][ri3]+newnr[htmp[l3].ib];
                  vali+=htmp[l3].vi*ix2[ir3];
#ifndef USEDII
                  vali+=htmp[l3].lvi*ix2[newnr[nbolo]+htmp[l3].ib*GAINSTEP+htmp[l3].gi];
#endif
                }
              }
                
              for (ib=0;ib<nbolo;ib++) {


#ifdef USEDII
                for (j=0;j<GAINSTEP;j++) {
                  vali+=dii[j+ib*GAINSTEP]  *ix2[newnr[nbolo]+ib*GAINSTEP+j];
                }
#endif

                if (nmatco>0) {
                  vali+=dcoi[ib+k*nbolo]  *ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+ib];
                }            
                if (nmatdust>0) {
                  vali+=ddusti[ib+k*nbolo]  *ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+ib];
                }            
                if (nfreefree>0) {
                  vali+=dfri[ib+k*nbolo]  *ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+ib];
                }            
                for (m=0;m<npixbeam;m++) {
                  vali+=ix2[newnr[nbolo]+nbolo*GAINSTEP2+ib+m*nbolo]*dpixi[ib+m*nbolo+k*nbolo*npixbeam];
                }   
              }

              double qri=0;

              for (l1=0;l1<ndata;l1++) {
                long ri1=htmp[l1].rg-globalBeginRing;
                long iri1=rgord[htmp[l1].ib][ri1]+newnr[htmp[l1].ib];

                if (flg_rg[htmp[l1].ib][ri1]!=0) {
                  double ww=htmp[l1].wp;
                  double val2=ix2[iri1]-(vali);

                  val2+=ix2[newnr[nbolo]+htmp[l1].ib*GAINSTEP+htmp[l1].gi]*htmp[l1].dip;

                  if (nmatco>0)   
                    val2+=ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+htmp[l1].ib]*htmp[l1].comap;
                  if (nmatdust>0) 
                    val2+=ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+htmp[l1].ib]*htmp[l1].dustmap;
                  if (nfreefree>0) 
                    val2+=ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nfreefree+htmp[l1].ib]*htmp[l1].freefree;

                  for (m=0;m<npixbeam;m++) {
                    val2+=ix2[newnr[nbolo]+nbolo*GAINSTEP2+htmp[l1].ib+m*nbolo]*htmp[l1].listofpix[m];
                  }
                    
                  qri-=ww*val2;
                    
                  q2[iri1]+=ww*val2;

                  q2[newnr[nbolo]+htmp[l1].ib*GAINSTEP+htmp[l1].gi]+=ww*val2*htmp[l1].dip;

                  ///////////// CO
                  if (nmatco>0) {
                    q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+htmp[l1].ib]+=ww*val2*htmp[l1].comap;
                  }
                    
                  ///////////// DUST
                  if (nmatdust>0) {
                    q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+htmp[l1].ib]+=ww*val2*htmp[l1].dustmap;
                  }

                  ///////////// FREEFREE
                  if (nfreefree>0) {
                    q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nfreefree+htmp[l1].ib]+=ww*val2*htmp[l1].freefree;
                  }
                    
                  for (j=0;j<npixbeam;j++) {
                    q2[newnr[nbolo]+nbolo*GAINSTEP2+j*nbolo+htmp[l1].ib]+=ww*val2*htmp[l1].listofpix[j];  
                  }
                    
                }
              }

              for (l2=0;l2<ndata;l2++) {
                long ri2=htmp[l2].rg-globalBeginRing;
                if (flg_rg[htmp[l2].ib][ri2]!=0) {
                  long ir=rgord[htmp[l2].ib][ri2]+newnr[htmp[l2].ib]; 
                  q2[ir]+=qri*htmp[l2].vi;
#ifndef USEDII
                  q2[newnr[nbolo]+htmp[l2].ib*GAINSTEP+htmp[l2].gi]+=qri*htmp[l2].lvi;
#endif
                }
              }
                
              for (ib=0;ib<nbolo;ib++) {


#ifdef USEDII
                for (j=0;j<GAINSTEP;j++) {
                  q2[newnr[nbolo]+GAINSTEP*ib+j]+=qri*dii[j+ib*GAINSTEP];
                }
#endif

                ///////////// CO
                if (nmatco>0) {
                  q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+ib]+=qri*dcoi[ib+k*nbolo];
                }
                ///////////// DUST
                if (nmatdust>0) {
                  q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+ib]+=qri*ddusti[ib+k*nbolo];
                }

                ///////////// FREEFREE
                if (nfreefree>0) {
                  q2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+ib]+=qri*dfri[ib+k*nbolo];
                }
                  
                for (j=0;j<npixbeam;j++) {
                  q2[newnr[nbolo]+nbolo*GAINSTEP2+j*nbolo+ib]+=qri*dpixi[ib+j*nbolo+k*nbolo*npixbeam];  
                }
              }
            }
          }
          if (rank==0) {
            double *lb = (double *) malloc(sizeof(double)*(nmatres));
            for (rrk=1;rrk<mpi_size;rrk++) {
              MPI_Recv(lb,sizeof(double)*(nmatres), MPI_BYTE, rrk,1034, MPI_COMM_WORLD,&statu);
              for (l=0;l<nmatres;l++) q2[l]+=lb[l];
            }

            free(lb);
          }
          else {
            MPI_Send(q2, sizeof(double)*(nmatres), MPI_BYTE, 0, 1034, MPI_COMM_WORLD);
          }
          
          if (rank==0) for (i=0; i < nmatres; i++) {
            r2[i] = b2[i] - q2[i];
          }
          
        }
      else
        {
          if (rank==0) for (i=0; i < nmatres; i++) r2[i] -= alpha * q2[i];
        }
      
      if (rank==0) for (i=0; i < nmatres; i++) s2[i] = r2[i] / hit2[i];

      delta_old = delta_new;
      if (rank==0) {
        delta_new=0;
        for (i=0; i < nmatres ; i++) delta_new += r2[i] * s2[i];
        beta = delta_new / delta_old;
        for (i=0; i < nmatres ; i++) d2[i] = s2[i] + beta * d2[i];
      }
      iter ++;
      MPI_Bcast(&delta_new, sizeof(double), MPI_BYTE, 0, MPI_COMM_WORLD);
      
    }
  
  if (rank==0) fprintf (stderr,"\niter = %d - delta0 = %lg - delta_new = %lg\n", 
                        iter, delta0, delta_new);
  if (rank==0) fprintf (stderr,"CG in iter = %d (max=%d)\n", iter, itermax);
  
  if (rank==0) memcpy(ix2  ,x2old,nmatres*sizeof (double));
  MPI_Bcast(ix2, sizeof(double)*(nmatres), MPI_BYTE, 0, MPI_COMM_WORLD);
  itbogo++;
}


void minimize_adu_nopol(double *in_off,double *ix2,double *gain,double *gaingi)
{
  MPI_Status statu;
  long i,rrk,j,k,l1,l2,ib;
  int itermax = 500;
  int iter;
  double  delta_new, delta_old, beta;
  double  alpha=1.0; // get rid of gcc "maybe-uninitialized" warning depending on optimsation level

  PIOLONG GAINSTEP2;
  int rank;
  int size;
  int mpi_size;
  int rank_size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank_size);
  rank=rank_size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  mpi_size=size;

  GAINSTEP2=GAINSTEPADU;

  nmatres=nbolo*(npixbeam+GAINSTEPADU)+nmatco+nmatdust+nfreefree;

  if (rank==0) {
    fprintf(stderr,"==============================\n\nITERATION %ld\n\n==============================\n",itbogo);
    fprintf(stderr,"GAIN ");
    for (i=0;i<nbolo;i++) fprintf(stderr,"%lg ",gaingi[i*GAINSTEP]);
    fprintf(stderr,"\n");
    if (GAINSTEP>1) {
      fprintf(stderr,"...\nGAIN ");
      for (i=0;i<nbolo;i++) fprintf(stderr,"%lg ",gaingi[i*GAINSTEP+(GAINSTEP-1)]);
      fprintf(stderr,"\n");
    }
  }
  if (itbogo==0) delta0=0;
  MPI_Barrier(MPI_COMM_WORLD);
  
  struct timeval tp1,tp2;
  gettimeofday(&tp1,NULL);

  
  iter = 0;
  memset(b2  ,0,nmatres*sizeof (double));
  memset(d2  ,0,nmatres*sizeof (double));
  memset(q2  ,0,nmatres*sizeof (double));
  memset(r2  ,0,nmatres*sizeof (double));
  memset(s2  ,0,nmatres*sizeof (double));
  memset(hit2,0,nmatres*sizeof (double));

  //==========================================
  //=  Compute second member
  //=
  //=
  long l,m;

  for (k=0;k<nnbpix;k++)  {
    //long imat=the_stat_pix[k];
    long ndata = loc_nhpix[k];
    hpix *htmp = loc_hpix[k];

    II[k]=0;
      
    for (l1=0;l1<ndata;l1++) {
      long ri1=htmp[l1].rg-globalBeginRing;
      if (flg_rg[htmp[l1].ib][ri1]!=0) {
        II[k]+=htmp[l1].w;
      }
    }

    if (ndata>0&&flgpix[k]>0) {

      double SI=0;

      memset(dcoi+k*nbolo,0,sizeof(double)*nbolo);
      memset(dfri+k*nbolo,0,sizeof(double)*nbolo);
      memset(dthetai+k*nbolo,0,sizeof(double)*nbolo);
      memset(ddusti+k*nbolo,0,sizeof(double)*nbolo);
      memset(dpixi+k*nbolo*npixbeam,0,sizeof(double)*nbolo*npixbeam);
      

      for (l1=0;l1<ndata;l1++) {
        long ri1=htmp[l1].rg-globalBeginRing;
        long iri1=rgord[htmp[l1].ib][ri1]+newnr[htmp[l1].ib];
        
        double g1=gaingi[htmp[l1].gi+htmp[l1].ib*GAINSTEP];
        
        if (flg_rg[htmp[l1].ib][ri1]!=0) {
          
          htmp[l1].wp=0;
          //htmp[l1].thsig=0;
          
          if (REMHDIP==0) SI+=htmp[l1].w*(htmp[l1].sig*g1-htmp[l1].fsl-htmp[l1].dip-gain[htmp[l1].adu+htmp[l1].ib*GAINSTEPADU]-in_off[iri1]);
          else SI+=htmp[l1].w*(htmp[l1].sig*g1-htmp[l1].fsl-htmp[l1].freefree-gain[htmp[l1].adu+htmp[l1].ib*GAINSTEPADU]-in_off[iri1]);
            
          if (nmatco>0) {
            dcoi[htmp[l1].ib+k*nbolo]   += htmp[l1].w*htmp[l1].comap;
          }
          if (nmatdust>0) {
            ddusti[htmp[l1].ib+k*nbolo] += htmp[l1].w*htmp[l1].dustmap;
          }
          if (nfreefree>0) {
            dfri[htmp[l1].ib+k*nbolo]   += htmp[l1].w*htmp[l1].freefree;
          }
          
          if (ittt>0) {
            for (m=0;m<npixbeam;m++)  {
              dpixi[nbolo*m+htmp[l1].ib+k*nbolo*npixbeam] += htmp[l1].w*htmp[l1].listofpix[m];
            }
          } 
        }
      }
      SI/=II[k];
      SSI[k]=SI;
      
      for (ib=0;ib<nbolo;ib++) {

        if (nmatco>0) {
          dcoi[ib+k*nbolo]=dcoi[ib+k*nbolo]/II[k];
        }
        if (nmatdust>0) {
          ddusti[ib+k*nbolo]=ddusti[ib+k*nbolo]/II[k];
        }
        if (nfreefree>0) {
          dfri[ib+k*nbolo]=dfri[ib+k*nbolo]/II[k];
        }
        
        if (ittt>0) {
          for (m=0;m<npixbeam;m++)  {
            dpixi[ib+m*nbolo+k*nbolo*npixbeam]=dpixi[ib+m*nbolo+k*nbolo*npixbeam]/II[k];
          }
        } 
      }

      for (l1=0;l1<ndata;l1++) {
        long ri1=htmp[l1].rg-globalBeginRing;
        if (flg_rg[htmp[l1].ib][ri1]!=0) {
          htmp[l1].vi=htmp[l1].w/II[k];
        }
      }
      for (l1=0;l1<ndata;l1++) {
        long ri1=htmp[l1].rg-globalBeginRing;
        
        double g1=gaingi[htmp[l1].gi+htmp[l1].ib*GAINSTEP];
        
        if (flg_rg[htmp[l1].ib][ri1]!=0) {
          double divi=NEP_tab[htmp[l1].ib]*htmp[l1].hit*g1;
          divi=divi*divi;          
          
          if (divi==0) htmp[l1].wp=0;
          else htmp[l1].wp=1/divi;
          
          if (itbogo==0) normaoff+=NEP_tab[htmp[l1].ib]*htmp[l1].hit;
        }
      }

      
      for (l1=0;l1<ndata;l1++) {
        long ri1=htmp[l1].rg-globalBeginRing;
        long iri1=rgord[htmp[l1].ib][ri1]+newnr[htmp[l1].ib];
        
        double g1=gaingi[htmp[l1].gi+htmp[l1].ib*GAINSTEP];
        
        if (flg_rg[htmp[l1].ib][ri1]!=0) {
          double ww=htmp[l1].wp;
          double tmp=((htmp[l1].sig*g1-htmp[l1].fsl-gain[htmp[l1].adu+htmp[l1].ib*GAINSTEPADU]-in_off[iri1])-(SI));
          if (REMHDIP==0) tmp-=htmp[l1].dip;
          else tmp-=htmp[l1].freefree;
            
          long l3;
          memset(ctmp,0,sizeof(double)*GAINSTEPADU*nbolo);
          
          for (l3=0;l3<ndata;l3++) {
            long ri3=htmp[l3].rg-globalBeginRing;
            if (flg_rg[htmp[l3].ib][ri3]!=0) {
              ctmp[htmp[l3].adu+GAINSTEPADU*htmp[l3].ib]-=htmp[l3].vi;
            }
          }
          ctmp[htmp[l1].adu+GAINSTEPADU*htmp[l1].ib]+=1;

          for (ib=0;ib<nbolo;ib++) cco[ib]=-dcoi[ib+k*nbolo];
          cco[htmp[l1].ib]+=htmp[l1].comap;
          
          for (ib=0;ib<nbolo;ib++) cdust[ib]=-ddusti[ib+k*nbolo];
          cdust[htmp[l1].ib]+=htmp[l1].dustmap;
          
          for (ib=0;ib<nbolo;ib++) ccfree[ib]=-dfri[ib+k*nbolo];
          ccfree[htmp[l1].ib]+=htmp[l1].freefree;
          
          for (m=0;m<npixbeam;m++) {
            for (ib=0;ib<nbolo;ib++) cpix[ib+m*nbolo]=-dpixi[ib+m*nbolo+k*nbolo*npixbeam];
            cpix[htmp[l1].ib+m*nbolo]+=htmp[l1].listofpix[m];
          }

          long ir;
          /////////////////  OFFSET
          
          for (l2=0;l2<ndata;l2++) {
            long ri2=htmp[l2].rg-globalBeginRing;
            if (flg_rg[htmp[l2].ib][ri2]!=0) {
              long iadu=htmp[l2].adu+GAINSTEPADU*htmp[l2].ib;
              b2[iadu]+=ww*tmp*ctmp[iadu];
              hit2[iadu]+=ww*ctmp[iadu]*ctmp[iadu];
            }
          }
          
          /////////////////  DIPOLE FIT
          
          for (ir=0;ir<nbolo;ir++) {
              
            ///////////// CO
            if (nmatco>0) {
              b2[nbolo*(npixbeam+GAINSTEP2)+ir]+=ww*tmp*cco[ir];
              hit2[nbolo*(npixbeam+GAINSTEP2)+ir]+=ww*cco[ir]*cco[ir];
            }
            
            ///////////// DUST
            
            if (nmatdust>0) {
              b2[nbolo*(npixbeam+GAINSTEP2)+nmatco+ir]+=ww*tmp*cdust[ir];
              hit2[nbolo*(npixbeam+GAINSTEP2)+nmatco+ir]+=ww*cdust[ir]*cdust[ir];
            }
            
            ///////////// FREEFREE
            
            if (nfreefree>0) {
              b2[nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+ir]+=ww*tmp*ccfree[ir];
              hit2[nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+ir]+=ww*ccfree[ir]*ccfree[ir];
            }
            
            ////////// PIXBEAM
            for (j=0;j<npixbeam;j++) {
              b2[nbolo*GAINSTEP2+j*nbolo+ir]+=ww*tmp*cpix[ir+j*nbolo];
              hit2[nbolo*GAINSTEP2+j*nbolo+ir]+=ww*cpix[ir+j*nbolo]*cpix[ir+j*nbolo];
            }
          }
        }
      }
    }
  }



  if (rank==0) {
    double *lb = (double *) malloc(sizeof(double)*(nmatres));
    for (rrk=1;rrk<mpi_size;rrk++) {
      MPI_Recv(lb,sizeof(double)*(nmatres), MPI_BYTE, rrk,1031, MPI_COMM_WORLD,&statu);
      for (l=0;l<nmatres;l++) b2[l]+=lb[l];
    }
    free(lb);
  }
  else MPI_Send(b2, sizeof(double)*(nmatres), MPI_BYTE, 0, 1031, MPI_COMM_WORLD);
  
  if (rank==0) {
    double *lb = (double *) malloc(sizeof(double)*(nmatres));
    for (rrk=1;rrk<mpi_size;rrk++) {
      MPI_Recv(lb,sizeof(double)*(nmatres), MPI_BYTE, rrk,1033, MPI_COMM_WORLD,&statu);

      for (l=0;l<nmatres;l++) hit2[l]+=lb[l];
    }
    free(lb);
  }
  else {
    MPI_Send(hit2, sizeof(double)*(nmatres), MPI_BYTE, 0, 1033, MPI_COMM_WORLD);
  }
  

  //==========================================================
  // Compute Ax
  //healspline_fill_s_with_PtP_s (hs_rec, Ndata, vec, x, q);
  // + 
  // Preconditionnement
  //
  //
  
  for (l=0;l<nmatres;l++) q2[l]=0;
  if (rank==0) {

    double soff=0;
    for (i=0;i<GAINSTEPADU*nbolo;i++) soff+=hit2[0]*ix2[i];
    for (i=0;i<GAINSTEPADU*nbolo;i++) q2[i]=soff;

    if (REMHDIP==1) {
      soff=0;
      for (i=newnr[nbolo];i<newnr[nbolo]+GAINSTEP*nbolo;i++) soff+=hit2[newnr[nbolo]]*ix2[i];
      for (i=newnr[nbolo];i<newnr[nbolo]+GAINSTEP*nbolo;i++) q2[i]=soff;
    }
    if (Param->flag_AVGR0==_PAR_TRUE) {
      soff=0;
      for (i=0;i<nbolo;i++) soff+=hit2[nbolo*(GAINSTEP2)]*
                              ix2[nbolo*(GAINSTEP2)+i];
      for (i=0;i<nbolo;i++) q2[nbolo*(GAINSTEP2)+i]=soff;
      for (i=0;i<nbolo;i++) b2[nbolo*(GAINSTEP2)+i]+=
                              hit2[nbolo*(GAINSTEP2)]*Param->AVGR0;
      
    }

    if (nmatco>0) {
      soff=0;
      for (i=0;i<nbolo;i++) soff+=hit2[nbolo*(npixbeam+GAINSTEP2)]
                              *ix2[nbolo*(npixbeam+GAINSTEP2)+i];
      for (i=0;i<nbolo;i++) q2[nbolo*(npixbeam+GAINSTEP2)+i]=soff;
      for (i=0;i<nbolo;i++) b2[nbolo*(npixbeam+GAINSTEP2)+i]+=
                              hit2[nbolo*(npixbeam+GAINSTEP2)]*Param->AVGCO;
    }
    if (nmatdust>0) {
      soff=0;
      for (i=0;i<nbolo;i++) 
        soff+=hit2[nbolo*(npixbeam+GAINSTEP2)+nmatco]*
          ix2[nbolo*(npixbeam+GAINSTEP2)+nmatco+i];
      for (i=0;i<nbolo;i++) q2[nbolo*(npixbeam+GAINSTEP2)+nmatco+i]=soff;
      for (i=0;i<nbolo;i++) b2[nbolo*(npixbeam+GAINSTEP2)+nmatco+i]+=
                              hit2[nbolo*(npixbeam+GAINSTEP2)+nmatco]*Param->AVGDUST;
    }
    if (nfreefree>0) {
      soff=0;
      for (i=0;i<nbolo;i++) 
        soff+=hit2[nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust]*
          ix2[nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+i];
      for (i=0;i<nbolo;i++) q2[nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+i]=soff;
    }
    
  }

  for (k=0;k<nnbpix;k++) if (flgpix[k]>0) {
    
    long ndata = loc_nhpix[k];
    if (ndata>0) {
#ifdef TIMING
      gettimeofday(&tp1,NULL);
#endif
      hpix *htmp = loc_hpix[k];

      double vali=0;
        
      long l3;
      for (l3=0;l3<ndata;l3++) {
        long ri3=htmp[l3].rg-globalBeginRing;
        if (flg_rg[htmp[l3].ib][ri3]!=0) {
          vali+=htmp[l3].vi*(ix2[htmp[l3].ib*GAINSTEPADU+htmp[l3].adu]);
        }
      }
      //fprintf(stderr,"I0 vali %lg %lg %lg\n",vali,valq,valu);
      
      for (ib=0;ib<nbolo;ib++) {

        if (nmatco>0) {
          vali+=dcoi[ib+k*nbolo]  *ix2[nbolo*(npixbeam+GAINSTEP2)+ib];
        }            
        if (nmatdust>0) {
          vali+=ddusti[ib+k*nbolo]  *ix2[nbolo*(npixbeam+GAINSTEP2)+nmatco+ib];
        }            
        if (nfreefree>0) {
          vali+=dfri[ib+k*nbolo]  *ix2[nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+ib];
        }
        //fprintf(stderr,"I3 vali %lg %lg %lg : %lg %lg %lg %lg\n",vali,valq,valu,
        //dfri[ib+k*nbolo],dfrq[ib+k*nbolo],dfru[ib+k*nbolo],
        //ix2[nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+ib]);            
        
        for (m=0;m<npixbeam;m++) {
          vali+=ix2[nbolo*GAINSTEP2+ib+m*nbolo]*dpixi[ib+m*nbolo+k*nbolo*npixbeam];
        }  
        
        //fprintf(stderr,"I4 vali %lg %lg %lg\n",vali,valq,valu);
      }
      double qri=0;
      
      for (l1=0;l1<ndata;l1++) {
        long ri1=htmp[l1].rg-globalBeginRing;

        if (flg_rg[htmp[l1].ib][ri1]!=0) {
          double ww=htmp[l1].wp;
          double val2=ix2[htmp[l1].ib*GAINSTEPADU+htmp[l1].adu]-(vali);

          //fprintf(stderr,"I2 val2 %lg\n",val2);
          if (nmatco>0)   

            val2+=ix2[nbolo*(npixbeam+GAINSTEP2)+htmp[l1].ib]*htmp[l1].comap;
          if (nmatdust>0) 
            val2+=ix2[nbolo*(npixbeam+GAINSTEP2)+nmatco+htmp[l1].ib]*htmp[l1].dustmap;
          if (nfreefree>0) 
            val2+=ix2[nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+htmp[l1].ib]*htmp[l1].freefree;
          
          //fprintf(stderr,"I3 val2 %lg\n",val2);
          for (m=0;m<npixbeam;m++) {
            val2+=ix2[nbolo*GAINSTEP2+htmp[l1].ib+m*nbolo]*htmp[l1].listofpix[m];
          }
          
          qri-=ww*val2;          

          memset(ctmp,0,sizeof(double)*GAINSTEPADU*nbolo);
          for (l3=0;l3<ndata;l3++) {
            long ri3=htmp[l3].rg-globalBeginRing;
            if (flg_rg[htmp[l3].ib][ri3]!=0) {
              ctmp[htmp[l3].adu+GAINSTEPADU*htmp[l3].ib]-=htmp[l3].vi;
            }
          }
          ctmp[htmp[l1].adu+GAINSTEPADU*htmp[l1].ib]+=1;

          for (l2=0;l2<ndata;l2++) {
            long ri2=htmp[l2].rg-globalBeginRing;
            if (flg_rg[htmp[l2].ib][ri2]!=0) {
              long iadu=htmp[l2].adu+GAINSTEPADU*htmp[l2].ib;
              q2[iadu]+=ww*val2*ctmp[iadu];
            }
          }

          ///////////// CO
          if (nmatco>0) {
            q2[nbolo*(npixbeam+GAINSTEP2)+htmp[l1].ib]+=ww*val2*htmp[l1].comap;
          }
          
          ///////////// DUST
          if (nmatdust>0) {
            q2[nbolo*(npixbeam+GAINSTEP2)+nmatco+htmp[l1].ib]+=ww*val2*htmp[l1].dustmap;
          }
          
          ///////////// FREEFREE
          if (nfreefree>0) {
            q2[nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+htmp[l1].ib]+=ww*val2*htmp[l1].freefree;
          }
          
          ////////// SYSTE
          for (j=0;j<npixbeam;j++) {
            q2[nbolo*GAINSTEP2+j*nbolo+htmp[l1].ib]+=ww*val2*htmp[l1].listofpix[j];  
          }
          
          
        }
      }
    
        
      for (ib=0;ib<nbolo;ib++) {

        ///////////// CO
        if (nmatco>0) {
          q2[nbolo*(npixbeam+GAINSTEP2)+ib]+=qri*dcoi[ib+k*nbolo];
        }
        ///////////// DUST
        if (nmatdust>0) {
          q2[nbolo*(npixbeam+GAINSTEP2)+nmatco+ib]+=qri*ddusti[ib+k*nbolo];
        }
        ///////////// FREEFREE
        if (nfreefree>0) {
          q2[nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+ib]+=qri*dfri[ib+k*nbolo];
        }
        
        for (j=0;j<npixbeam;j++) {
          q2[nbolo*GAINSTEP2+j*nbolo+ib]+=qri*dpixi[ib+j*nbolo+k*nbolo*npixbeam];  
        }
      }
    }
      
#ifdef TIMING
      gettimeofday(&tp2,NULL);
      double dt=(double)(tp2.tv_sec-tp1.tv_sec)+(1E-6)*(tp2.tv_usec-tp1.tv_usec);
      dthit[ndata]+=dt;
      ndthit[ndata]+=1;
#endif
  }

  if (rank==0) {
    double *lb = (double *) malloc(sizeof(double)*(nmatres));
    for (rrk=1;rrk<mpi_size;rrk++) {
      MPI_Recv(lb,sizeof(double)*(nmatres), MPI_BYTE, rrk,1032, MPI_COMM_WORLD,&statu);
      for (l=0;l<nmatres;l++) q2[l]+=lb[l];
    }
    free(lb);
  }
  else {
    MPI_Send(q2, sizeof(double)*(nmatres), MPI_BYTE, 0, 1032, MPI_COMM_WORLD);
  }

  if (rank==0) fprintf(stderr,"QQ2 %lg\n",q2[0]);

  if (rank==0) {
    for (i=0; i < nmatres; i++)
      {
        r2[i] = b2[i] - q2[i];
        d2[i] = r2[i] / hit2[i];
      }
  }

  
  double delta_new_tmp = 0.0;
  if (rank==0) {
    for (i = 0; i < nmatres; i++) {
      delta_new_tmp += b2[i] ;
      if (isnan(b2[i])) {
        fprintf(stderr,"NAN B2 PBS %ld\n",(long) i);
      }
      if (isnan(d2[i])) {
        fprintf(stderr,"NAN D2 PBS %ld %lg\n",(long) i,hit2[i]);
      }
    }
    fprintf(stderr,"B2 %lg\n",delta_new_tmp);
  
    delta_new_tmp = 0.0;
    for (i = 0; i < nmatres; i++) {
      delta_new_tmp += q2[i] ;
      if (isnan(q2[i])) {
        fprintf(stderr,"NAN Q2 PBS %ld\n",(long) i);
      }
    }
    fprintf(stderr,"Q2 %lg\n",delta_new_tmp);
    delta_new_tmp = 0.0;
    for (i = 0; i < nmatres; i++) {
      delta_new_tmp += q2[i]-b2[i] ;
      //fprintf(stderr,"B2 Q2 B2-Q2 [%ld]: %lg\t%lg\t%lg\n",(long) i,b2[i],q2[i],q2[i]-b2[i]);
    }
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
  delta_new_tmp = 0.0;
  if (rank==0) for (i=0; i < nmatres; i++) {
    delta_new_tmp += r2[i] * d2[i];
  }
  
  delta_new=0;
  for (i=0;i<mpi_size;i++) {
    double tmp=delta_new_tmp;
    MPI_Bcast(&tmp, sizeof(double), MPI_BYTE, i, MPI_COMM_WORLD);
    delta_new+=tmp;
  }
  if (itbogo==0) delta0 = delta_new;
  if (rank==0) fprintf (stderr, "min_adu_nop() iter = %d - delta0 = %lg - delta_new = %lg\n", iter, delta0, delta_new);
  int testwrit=0;

  if (itbogo==0) MPI_Bcast(&normaoff, sizeof(double), MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&delta0, sizeof(double), MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&delta_new, sizeof(double), MPI_BYTE, 0, MPI_COMM_WORLD);
  if (rank==0) {
    memcpy(x2old,ix2,nmatres*sizeof(double));
    testwrit=1;
  }


  while ((iter < itermax)  && ((delta_new) > delta0*1E-24) && ((delta_new) > 1E-20)) //Param->XI2STOP))
    {
      // q <= Ad
      //if (rank==0&&mindelta>delta_new) {
      if (rank==0) {
        memcpy(x2old,ix2,nmatres*sizeof(double));
        testwrit=1;
      }
    //  else testwrit=0;

      //healspline_fill_s_with_PtP_s (hs_rec, Ndata, vec, d, q);
      //for (i=0;i<mpi_size;i++) {
      //MPI_Bcast(d+tab_begr[i]*2, sizeof(double)*(tab_edr[i]-tab_begr[i]+1)*2, MPI_BYTE, i, MPI_COMM_WORLD);
      //}
      MPI_Bcast(d2, sizeof(double)*(nmatres), MPI_BYTE, 0, MPI_COMM_WORLD);
      
      // ===========================================================================
      // PROJECTION DE d dans q
      //
      
      for (l=0;l<nmatres;l++) q2[l]=0;
      if (rank==0) {
        double soff=0;
        for (i=0;i<GAINSTEPADU*nbolo;i++) soff+=hit2[0]*d2[i];
        for (i=0;i<GAINSTEPADU*nbolo;i++) q2[i]=soff;

        if (REMHDIP==1) {
          soff=0;
          for (i=newnr[nbolo];i<newnr[nbolo]+GAINSTEP*nbolo;i++) soff+=hit2[newnr[nbolo]]*d2[i];
          for (i=newnr[nbolo];i<newnr[nbolo]+GAINSTEP*nbolo;i++) q2[i]=soff;
        }

        if (Param->flag_AVGR0==_PAR_TRUE) {
          soff=0;
          for (i=0;i<nbolo;i++) soff+=hit2[nbolo*(GAINSTEP2)]*
                                  d2[nbolo*(GAINSTEP2)+i];
          for (i=0;i<nbolo;i++) q2[nbolo*(GAINSTEP2)+i]=soff; 
        }
        if (nmatco>0) {
          soff=0;
          for (i=0;i<nbolo;i++) soff+=hit2[nbolo*(npixbeam+GAINSTEP2)]
                                  *d2[nbolo*(npixbeam+GAINSTEP2)+i];
          for (i=0;i<nbolo;i++) q2[nbolo*(npixbeam+GAINSTEP2)+i]=soff;
        }
        if (nfreefree>0) {
          soff=0;
          for (i=0;i<nbolo;i++) soff+=hit2[nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust]
                                  *d2[nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+i];
          for (i=0;i<nbolo;i++) q2[nbolo*(npixbeam+GAINSTEP2)+i+nmatco+nmatdust]=soff;
        }
        if (nmatdust>0) {
          soff=0;
          for (i=0;i<nbolo;i++) 
            soff+=hit2[nbolo*(npixbeam+GAINSTEP2)+nmatco]*
              d2[nbolo*(npixbeam+GAINSTEP2)+nmatco+i];
          for (i=0;i<nbolo;i++) q2[nbolo*(npixbeam+GAINSTEP2)+nmatco+i]=soff;
        }
      }
      for (k=0;k<nnbpix;k++) if (flgpix[k]>0) {
        
        //long imat=the_stat_pix[k];
        long ndata = loc_nhpix[k];
        if (ndata>0) {
          hpix *htmp = loc_hpix[k];

          double vali=0;

          long l3;
          for (l3=0;l3<ndata;l3++) {
            long ri3=htmp[l3].rg-globalBeginRing;
            if (flg_rg[htmp[l3].ib][ri3]!=0) {
              vali+=htmp[l3].vi*(d2[htmp[l3].ib*GAINSTEPADU+htmp[l3].adu]);
            }
          }
          for (ib=0;ib<nbolo;ib++) {
            
            if (nmatco>0) {
              vali+=dcoi[ib+k*nbolo]  *d2[nbolo*(npixbeam+GAINSTEP2)+ib];
            }            
            if (nmatdust>0) {
              vali+=ddusti[ib+k*nbolo]  *d2[nbolo*(npixbeam+GAINSTEP2)+nmatco+ib];
            }            
            if (nfreefree>0) {
              vali+=dfri[ib+k*nbolo]  *d2[nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+ib];
            }            
            for (m=0;m<npixbeam;m++) {
              vali+=d2[nbolo*GAINSTEP2+ib+m*nbolo]*dpixi[ib+m*nbolo+k*nbolo*npixbeam];
            }   
          }

          double qri=0;

          for (l1=0;l1<ndata;l1++) {
            long ri1=htmp[l1].rg-globalBeginRing;

            if (flg_rg[htmp[l1].ib][ri1]!=0) {
              double ww=htmp[l1].wp;
              double val2=d2[htmp[l1].ib*GAINSTEPADU+htmp[l1].adu]-(vali);

              if (nmatco>0)   
                val2+=d2[nbolo*(npixbeam+GAINSTEP2)+htmp[l1].ib]*htmp[l1].comap;
              if (nmatdust>0) 
                val2+=d2[nbolo*(npixbeam+GAINSTEP2)+nmatco+htmp[l1].ib]*htmp[l1].dustmap;
              if (nfreefree>0) 
                val2+=d2[nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+htmp[l1].ib]*htmp[l1].freefree;
                
              for (m=0;m<npixbeam;m++) {
                val2+=d2[nbolo*GAINSTEP2+htmp[l1].ib+m*nbolo]*htmp[l1].listofpix[m];
              }
              
                
              qri-=ww*val2;

              memset(ctmp,0,sizeof(double)*GAINSTEPADU*nbolo);
              for (l3=0;l3<ndata;l3++) {
                long ri3=htmp[l3].rg-globalBeginRing;
                if (flg_rg[htmp[l3].ib][ri3]!=0) {
                  ctmp[htmp[l3].adu+GAINSTEPADU*htmp[l3].ib]-=htmp[l3].vi;
                }
              }
              ctmp[htmp[l1].adu+GAINSTEPADU*htmp[l1].ib]+=1;

              for (l2=0;l2<ndata;l2++) {
                long ri2=htmp[l2].rg-globalBeginRing;
                if (flg_rg[htmp[l2].ib][ri2]!=0) {
                  long iadu=htmp[l2].adu+GAINSTEPADU*htmp[l2].ib;
                  q2[iadu]+=ww*val2*ctmp[iadu];
                }
              }

              ///////////// CO
              if (nmatco>0) {
                q2[nbolo*(npixbeam+GAINSTEP2)+htmp[l1].ib]+=ww*val2*htmp[l1].comap;
              }
                
              ///////////// DUST
              if (nmatdust>0) {
                q2[nbolo*(npixbeam+GAINSTEP2)+nmatco+htmp[l1].ib]+=ww*val2*htmp[l1].dustmap;
              }
                
              ///////////// FREEFREE
              if (nfreefree>0) {
                q2[nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+htmp[l1].ib]+=ww*val2*htmp[l1].freefree;
              }
            
                
              for (j=0;j<npixbeam;j++) {
                q2[nbolo*GAINSTEP2+j*nbolo+htmp[l1].ib]+=ww*val2*htmp[l1].listofpix[j];  
              }
                
            }
          }

          for (ib=0;ib<nbolo;ib++) {

            ///////////// CO
            if (nmatco>0) {
              q2[nbolo*(npixbeam+GAINSTEP2)+ib]+=qri*dcoi[ib+k*nbolo];
            }
            ///////////// DUST
            if (nmatdust>0) {
              q2[nbolo*(npixbeam+GAINSTEP2)+nmatco+ib]+=qri*ddusti[ib+k*nbolo];
            }
            ///////////// FREEFREE
            if (nfreefree>0) {
              q2[nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+ib]+=qri*dfri[ib+k*nbolo];
            }

            for (j=0;j<npixbeam;j++) {
              q2[nbolo*GAINSTEP2+j*nbolo+ib]+=qri*dpixi[ib+j*nbolo+k*nbolo*npixbeam];  
            }
          }
        }
      }
      if (rank==0) {
        double *lb = (double *) malloc(sizeof(double)*(nmatres));
        for (rrk=1;rrk<mpi_size;rrk++) {
          MPI_Recv(lb,sizeof(double)*(nmatres), MPI_BYTE, rrk,1034, MPI_COMM_WORLD,&statu);
          for (l=0;l<nmatres;l++) q2[l]+=lb[l];
        }
        free(lb);
      }
      else {
        MPI_Send(q2, sizeof(double)*(nmatres), MPI_BYTE, 0, 1034, MPI_COMM_WORLD);
      }
      double dtq = 0.0;
      if (rank==0) {
        for (i=0; i < nmatres; i++) dtq += d2[i] * q2[i];
        alpha = delta_new / dtq;
        for (i=0; i < nmatres ; i++) ix2[i] += alpha * d2[i];
      }


      gettimeofday(&tp2,NULL);
      if (rank==0&&iter%20==0) fprintf (stderr,"min_adu_nop2() iter = %d - delta0 = %lg - delta_new = %lg %ld %3lfs\n", iter, delta0, delta_new,
                          (long) testwrit,(double)(tp2.tv_sec-tp1.tv_sec)+(1E-6)*(tp2.tv_usec-tp1.tv_usec));

      //if (rank==0) fprintf (stderr,".");
      gettimeofday(&tp1,NULL);
          
      if (iter % 100 == 0 && iter !=0)
        {
          // Use the best case
          memcpy(ix2,x2old,nmatres*sizeof(double));
          //for (i=0;i<mpi_size;i++) {
          //  MPI_Bcast(x+tab_begr[i]*2, sizeof(double)*(tab_edr[i]-tab_begr[i]+1)*2, MPI_BYTE, i, MPI_COMM_WORLD);
          //}
          MPI_Bcast(ix2, sizeof(double)*(nmatres), MPI_BYTE, 0, MPI_COMM_WORLD);

          
      // ===========================================================================
      // PROJECTION DE x dans q
      //
          
          for (l=0;l<nmatres;l++) q2[l]=0;
          if (rank==0) {
            double soff=0;
            for (i=0;i<GAINSTEPADU*nbolo;i++) soff+=hit2[0]*ix2[i];
            for (i=0;i<GAINSTEPADU*nbolo;i++) q2[i]=soff;

            if (REMHDIP==1) {
              soff=0;
              for (i=newnr[nbolo];i<newnr[nbolo]+GAINSTEP*nbolo;i++) soff+=hit2[newnr[nbolo]]*ix2[i];
              for (i=newnr[nbolo];i<newnr[nbolo]+GAINSTEP*nbolo;i++) q2[i]=soff;
            }

            if (Param->flag_AVGR0==_PAR_TRUE) {
              soff=0;
              for (i=0;i<nbolo;i++) soff+=hit2[nbolo*(GAINSTEP2)]*
                                      ix2[nbolo*(GAINSTEP2)+i];
              for (i=0;i<nbolo;i++) q2[nbolo*(GAINSTEP2)+i]=soff; 
            }
            if (nmatco>0) {
              soff=0;
              for (i=0;i<nbolo;i++) soff+=hit2[nbolo*(npixbeam+GAINSTEP2)]
                                      *ix2[nbolo*(npixbeam+GAINSTEP2)+i];
              for (i=0;i<nbolo;i++) q2[nbolo*(npixbeam+GAINSTEP2)+i]=soff;
            }
            if (nmatdust>0) {
              soff=0;
              for (i=0;i<nbolo;i++) 
                soff+=hit2[nbolo*(npixbeam+GAINSTEP2)+nmatco]*
                  ix2[nbolo*(npixbeam+GAINSTEP2)+nmatco+i];
              for (i=0;i<nbolo;i++) q2[nbolo*(npixbeam+GAINSTEP2)+nmatco+i]=soff;
            }
            if (nfreefree>0) {
              soff=0;
              for (i=0;i<nbolo;i++) 
                soff+=hit2[nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust]*
                  ix2[nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+i];
              for (i=0;i<nbolo;i++) q2[nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+i]=soff;
            }
            
            
          }

          for (k=0;k<nnbpix;k++) if (flgpix[k]>0) {
        
            //long imat=the_stat_pix[k];
            long ndata = loc_nhpix[k];
            if (ndata>0) {
              hpix *htmp = loc_hpix[k];

              double vali=0;

              long l3;
              for (l3=0;l3<ndata;l3++) {
                long ri3=htmp[l3].rg-globalBeginRing;
                if (flg_rg[htmp[l3].ib][ri3]!=0) {
                  vali+=htmp[l3].vi*ix2[htmp[l3].ib*GAINSTEPADU+htmp[l3].adu];
                }
              }
                
              for (ib=0;ib<nbolo;ib++) {

                if (nmatco>0) {
                  vali+=dcoi[ib+k*nbolo]  *ix2[nbolo*(npixbeam+GAINSTEP2)+ib];
                }            
                if (nmatdust>0) {
                  vali+=ddusti[ib+k*nbolo]  *ix2[nbolo*(npixbeam+GAINSTEP2)+nmatco+ib];
                }            
                if (nfreefree>0) {
                  vali+=dfri[ib+k*nbolo]  *ix2[nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+ib];
                }            
                for (m=0;m<npixbeam;m++) {
                  vali+=ix2[nbolo*GAINSTEP2+ib+m*nbolo]*dpixi[ib+m*nbolo+k*nbolo*npixbeam];
                }   
              }

              double qri=0;

              for (l1=0;l1<ndata;l1++) {
                long ri1=htmp[l1].rg-globalBeginRing;

                if (flg_rg[htmp[l1].ib][ri1]!=0) {
                  double ww=htmp[l1].wp;
                  double val2=ix2[htmp[l1].ib*GAINSTEPADU+htmp[l1].adu]-(vali);

                  if (nmatco>0)   
                    val2+=ix2[nbolo*(npixbeam+GAINSTEP2)+htmp[l1].ib]*htmp[l1].comap;
                  if (nmatdust>0) 
                    val2+=ix2[nbolo*(npixbeam+GAINSTEP2)+nmatco+htmp[l1].ib]*htmp[l1].dustmap;
                  if (nfreefree>0) 
                    val2+=ix2[nbolo*(npixbeam+GAINSTEP2)+nmatco+nfreefree+htmp[l1].ib]*htmp[l1].freefree;

                  for (m=0;m<npixbeam;m++) {
                    val2+=ix2[nbolo*GAINSTEP2+htmp[l1].ib+m*nbolo]*htmp[l1].listofpix[m];
                  }
                    
                  qri-=ww*val2;

                  memset(ctmp,0,sizeof(double)*GAINSTEPADU*nbolo);
                  for (l3=0;l3<ndata;l3++) {
                    long ri3=htmp[l3].rg-globalBeginRing;
                    if (flg_rg[htmp[l3].ib][ri3]!=0) {
                      ctmp[htmp[l3].adu+GAINSTEPADU*htmp[l3].ib]-=htmp[l3].vi;
                    }
                  }
                  ctmp[htmp[l1].adu+GAINSTEPADU*htmp[l1].ib]+=1;
                  
                  for (l2=0;l2<ndata;l2++) {
                    long ri2=htmp[l2].rg-globalBeginRing;
                    if (flg_rg[htmp[l2].ib][ri2]!=0) {
                      long iadu=htmp[l2].adu+GAINSTEPADU*htmp[l2].ib;
                      q2[iadu]+=ww*val2*ctmp[iadu];
                    }
                  }


                  ///////////// CO
                  if (nmatco>0) {
                    q2[nbolo*(npixbeam+GAINSTEP2)+htmp[l1].ib]+=ww*val2*htmp[l1].comap;
                  }
                    
                  ///////////// DUST
                  if (nmatdust>0) {
                    q2[nbolo*(npixbeam+GAINSTEP2)+nmatco+htmp[l1].ib]+=ww*val2*htmp[l1].dustmap;
                  }

                  ///////////// FREEFREE
                  if (nfreefree>0) {
                    q2[nbolo*(npixbeam+GAINSTEP2)+nmatco+nfreefree+htmp[l1].ib]+=ww*val2*htmp[l1].freefree;
                  }
                    
                  for (j=0;j<npixbeam;j++) {
                    q2[nbolo*GAINSTEP2+j*nbolo+htmp[l1].ib]+=ww*val2*htmp[l1].listofpix[j];  
                  }
                    
                }
              }
                
              for (ib=0;ib<nbolo;ib++) {

                ///////////// CO
                if (nmatco>0) {
                  q2[nbolo*(npixbeam+GAINSTEP2)+ib]+=qri*dcoi[ib+k*nbolo];
                }
                ///////////// DUST
                if (nmatdust>0) {
                  q2[nbolo*(npixbeam+GAINSTEP2)+nmatco+ib]+=qri*ddusti[ib+k*nbolo];
                }

                ///////////// FREEFREE
                if (nfreefree>0) {
                  q2[nbolo*(npixbeam+GAINSTEP2)+nmatco+nmatdust+ib]+=qri*dfri[ib+k*nbolo];
                }
                  
                for (j=0;j<npixbeam;j++) {
                  q2[nbolo*GAINSTEP2+j*nbolo+ib]+=qri*dpixi[ib+j*nbolo+k*nbolo*npixbeam];  
                }
              }
            }
          }
          if (rank==0) {
            double *lb = (double *) malloc(sizeof(double)*(nmatres));
            for (rrk=1;rrk<mpi_size;rrk++) {
              MPI_Recv(lb,sizeof(double)*(nmatres), MPI_BYTE, rrk,1034, MPI_COMM_WORLD,&statu);
              for (l=0;l<nmatres;l++) q2[l]+=lb[l];
            }

            free(lb);
          }
          else {
            MPI_Send(q2, sizeof(double)*(nmatres), MPI_BYTE, 0, 1034, MPI_COMM_WORLD);
          }
          
          if (rank==0) for (i=0; i < nmatres; i++) {
            r2[i] = b2[i] - q2[i];
          }
          
        }
      else
        {
          if (rank==0) for (i=0; i < nmatres; i++) r2[i] -= alpha * q2[i];
        }
      
      if (rank==0) for (i=0; i < nmatres; i++) s2[i] = r2[i] / hit2[i];

      delta_old = delta_new;
      if (rank==0) {
        delta_new=0;
        for (i=0; i < nmatres ; i++) delta_new += r2[i] * s2[i];
        beta = delta_new / delta_old;
        for (i=0; i < nmatres ; i++) d2[i] = s2[i] + beta * d2[i];
      }
      iter ++;
      MPI_Bcast(&delta_new, sizeof(double), MPI_BYTE, 0, MPI_COMM_WORLD);
      
    }
  
  if (rank==0) fprintf (stderr,"\nmin_adu_nop3() iter = %d - delta0 = %lg - delta_new = %lg\n", 
                        iter, delta0, delta_new);
  if (rank==0) fprintf (stderr,"CG in iter = %d (max=%d)\n", iter, itermax);
  
  if (rank==0) memcpy(ix2  ,x2old,nmatres*sizeof (double));
  MPI_Bcast(ix2, sizeof(double)*(nmatres), MPI_BYTE, 0, MPI_COMM_WORLD);
  itbogo++;
}



void minimize_optimize(double *ix2,double *x3,double *gain)
{
  MPI_Status statu;
  long i,rrk,j,k,l1,l2;
  int itermax = 200;
  int iter;
  double delta_new, alpha=0, delta_old, beta;


  int rank;
  int size;
  int mpi_size;
  int rank_size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank_size);
  rank=rank_size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  mpi_size=size;

  if (rank==0) {
    fprintf(stderr,"==============================\n\nOPTIMIZE %ld\n\n==============================\n",itbogo);
    for (j=0;j<GAINSTEP;j++) {
      fprintf(stderr,"GAIN ");
      for (i=0;i<nbolo;i++) fprintf(stderr,"%lg ",gain[i*GAINSTEP+j]);
      fprintf(stderr,"\n");
    } 
  }
  if (itbogo==0) delta0=0;
  MPI_Barrier(MPI_COMM_WORLD);
  
  struct timeval tp1,tp2;
  gettimeofday(&tp1,NULL);

  nmatres=newnr2[nbolo];
  
  iter = 0;
  memset(b2  ,0,nmatres*sizeof (double));
  memset(d2  ,0,nmatres*sizeof (double));
  memset(q2  ,0,nmatres*sizeof (double));
  memset(r2  ,0,nmatres*sizeof (double));
  memset(s2  ,0,nmatres*sizeof (double));
  memset(hit2,0,nmatres*sizeof (double));

  //==========================================
  //=  Compute second member
  //=
  //=
  long l,m;

  ////// BUILD B2

  for (k=0;k<nnbpix;k++)  {
    //long imat=the_stat_pix[k];
    long ndata = loc_nhpix[k];
    hpix *htmp = loc_hpix[k];
    II[k]=0;
    IQ[k]=0;
    IU[k]=0;
    QQ[k]=0;
    UU[k]=0;
    QU[k]=0;

    for (l1=0;l1<ndata;l1++) {
      long ri1=htmp[l1].hrg-globalBeginRing*CUTRG;
      double CO1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].co
                                        -dpsisi[htmp[l1].ib]*htmp[l1].si);
      double SI1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].si
                                        +dpsisi[htmp[l1].ib]*htmp[l1].co);
      if (flg_rg[htmp[l1].ib][ri1/CUTRG]!=0) {
        II[k]+=htmp[l1].w;
        IQ[k]+=htmp[l1].w*CO1;
        IU[k]+=htmp[l1].w*SI1;
        QQ[k]+=htmp[l1].w*CO1*CO1;
        QU[k]+=htmp[l1].w*SI1*CO1;
        UU[k]+=htmp[l1].w*SI1*SI1;
      }
    }
  
    if (ndata>0&&flgpix[k]>0) {

      double SI=0;
      double SQ=0;
      double SU=0;
      for (l1=0;l1<ndata;l1++) {
        long ri1=htmp[l1].rg-globalBeginRing;
        //long iri1=rgord[htmp[l1].ib][ri1]+newnr[htmp[l1].ib];
        
        double g1=gain[htmp[l1].adu+htmp[l1].ib*GAINSTEP];
        
        ri1=htmp[l1].hrg-globalBeginRing*CUTRG;
        //iri1=rgord2[htmp[l1].ib][ri1]+newnr2[htmp[l1].ib];
        double CO1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].co
                                          -dpsisi[htmp[l1].ib]*htmp[l1].si);
        double SI1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].si
                                          +dpsisi[htmp[l1].ib]*htmp[l1].co);
        if (flg_rg[htmp[l1].ib][ri1/CUTRG]!=0) {
          
          double tmp=0;
          
          if (nmatco>0) {
            tmp+=ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP)+htmp[l1].ib]*htmp[l1].comap;
          }
          if (nmatdust>0) {
            tmp+=ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP)+nmatco+htmp[l1].ib]*htmp[l1].dustmap;
          }
          if (nfreefree>0) {
            tmp+=ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP)+nmatco+nmatdust+htmp[l1].ib]*htmp[l1].freefree;
          }
          
          if (ittt>0) {
            for (m=0;m<npixbeam;m++) {
              tmp+= htmp[l1].listofpix[m]*
                ix2[newnr[nbolo]+nbolo*GAINSTEP+htmp[l1].ib+m*nbolo];
            }
          }
          
          double msig=((htmp[l1].sig*g1-htmp[l1].fsl-htmp[l1].dip)-tmp);
          
          SI+=htmp[l1].w*msig;
          SQ+=htmp[l1].w*msig*CO1;
          SU+=htmp[l1].w*msig*SI1;
        }
      }
    
    
      solvemap(&SI,&SQ,&SU,II[k],IQ[k],IU[k],QQ[k],QU[k],UU[k]);
      SSI[k]=SI;
      SSQ[k]=SQ;
      SSU[k]=SU;

      for (l1=0;l1<ndata;l1++) {
        long ri1=htmp[l1].hrg-globalBeginRing*CUTRG;
        double CO1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].co
                                          -dpsisi[htmp[l1].ib]*htmp[l1].si);
        double SI1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].si
                                          +dpsisi[htmp[l1].ib]*htmp[l1].co);
        if (flg_rg[htmp[l1].ib][ri1/CUTRG]!=0) {
          double li=htmp[l1].w,lco=CO1*htmp[l1].w,lsi=SI1*htmp[l1].w;
          solvemap(&li,&lco,&lsi,II[k],IQ[k],IU[k],QQ[k],QU[k],UU[k]);
          htmp[l1].vi=li;
          htmp[l1].vq=lco;
          htmp[l1].vu=lsi;
        }
      }
    
  
      for (l1=0;l1<ndata;l1++) {
        long ri1=htmp[l1].rg-globalBeginRing;
        long iri1=rgord[htmp[l1].ib][ri1]+newnr[htmp[l1].ib];
        
        double g1=gain[htmp[l1].adu+htmp[l1].ib*GAINSTEP];
        
        ri1=htmp[l1].hrg-globalBeginRing*CUTRG;
        iri1=rgord2[htmp[l1].ib][ri1]+newnr2[htmp[l1].ib];
        double CO1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].co
                                          -dpsisi[htmp[l1].ib]*htmp[l1].si);
        double SI1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].si
                                          +dpsisi[htmp[l1].ib]*htmp[l1].co);
        if (flg_rg[htmp[l1].ib][ri1/CUTRG]!=0) {
          
          double tmp=0;
          
          if (nmatco>0) {
            tmp+=ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP)+htmp[l1].ib]*htmp[l1].comap;
          }
          if (nmatdust>0) {
            tmp+=ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP)+nmatco+htmp[l1].ib]*htmp[l1].dustmap;
          }
          if (nfreefree>0) {
            tmp+=ix2[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP)+nmatco+nmatdust+htmp[l1].ib]*htmp[l1].freefree;
          }
          
          if (ittt>0) {
            for (m=0;m<npixbeam;m++) {
              tmp+= htmp[l1].listofpix[m]*
                ix2[newnr[nbolo]+nbolo*GAINSTEP+htmp[l1].ib+m*nbolo];
            }
          }
          
          double msig=((htmp[l1].sig*g1-htmp[l1].fsl-htmp[l1].dip)-tmp);
          
          double ww=htmp[l1].wp;
          tmp=(msig-(SSI[k]+CO1*SSQ[k]+SI1*SSU[k]));
          long l3;
          for (l3=0;l3<ndata;l3++) {
            long ri3=htmp[l3].hrg-globalBeginRing*CUTRG;
            if (flg_rg[htmp[l3].ib][ri3/CUTRG]!=0) {
              long ir3=rgord2[htmp[l3].ib][ri3]+newnr2[htmp[l3].ib];
              ctmp[ir3]=-(htmp[l3].vi+CO1*htmp[l3].vq+SI1*htmp[l3].vu);
            }
          }
          ctmp[iri1]+=1;

          /////////////////  OFFSET
          for (l2=0;l2<ndata;l2++) {
            long ri2=htmp[l2].hrg-globalBeginRing*CUTRG;
            if (flg_rg[htmp[l2].ib][ri2/CUTRG]!=0) {
              long ir=rgord2[htmp[l2].ib][ri2]+newnr2[htmp[l2].ib]; 
              b2[ir]+=ww*tmp*ctmp[ir];
              hit2[ir]+=ww*ctmp[ir]*ctmp[ir];
            }
          }
        }
      } 
    }
  }

  

  if (rank==0) {
    double *lb = (double *) malloc(sizeof(double)*(nmatres));
    for (rrk=1;rrk<mpi_size;rrk++) {
      MPI_Recv(lb,sizeof(double)*(nmatres), MPI_BYTE, rrk,1031, MPI_COMM_WORLD,&statu);
      for (l=0;l<nmatres;l++) b2[l]+=lb[l];
    }
    free(lb);
  }
  else MPI_Send(b2, sizeof(double)*(nmatres), MPI_BYTE, 0, 1031, MPI_COMM_WORLD);
  
  if (rank==0) {
    double *lb = (double *) malloc(sizeof(double)*(nmatres));
    for (rrk=1;rrk<mpi_size;rrk++) {
      MPI_Recv(lb,sizeof(double)*(nmatres), MPI_BYTE, rrk,1033, MPI_COMM_WORLD,&statu);

      for (l=0;l<nmatres;l++) hit2[l]+=lb[l];
    }
    free(lb);
  }
  else {
    MPI_Send(hit2, sizeof(double)*(nmatres), MPI_BYTE, 0, 1033, MPI_COMM_WORLD);
  }
  
  // starting point x = solution
  memset(x3  ,0,nmatres*sizeof (double));

  //==========================================================
  // Compute Ax
  //healspline_fill_s_with_PtP_s (hs_rec, Ndata, vec, x, q);
  // + 
  // Preconditionnement
  //
  //
  
  for (l=0;l<nmatres;l++) q2[l]=0;
  if (rank==0) {

    double soff=0;
    for (i=0;i<newnr2[nbolo];i++) soff+=hit2[0]*x3[i];
    for (i=0;i<newnr2[nbolo];i++) q2[i]=soff;
  }

#ifdef TIMING
#define TIMING
  long maxhit=0;
  for (k=0;k<nnbpix;k++) if (flgpix[k]>0) {
    long ndata = loc_nhpix[k];
    if (ndata>maxhit) maxhit=ndata;
  }
  if (rank==0) fprintf(stderr,"rank %ld %ld\n",(long) rank,(long) maxhit);
  double *dthit = (double *) malloc(maxhit*sizeof(double));
  double *ndthit = (double *) malloc(maxhit*sizeof(double));
  memset(dthit,0,maxhit*sizeof(double));
  memset(ndthit,0,maxhit*sizeof(double));
#endif
  for (k=0;k<nnbpix;k++) if (flgpix[k]>0) {
    
    //long imat=the_stat_pix[k];
    long ndata = loc_nhpix[k];
    if (ndata>0) {
#ifdef TIMING
      gettimeofday(&tp1,NULL);
#endif
      hpix *htmp = loc_hpix[k];

#ifdef TRUECP
#else
      double vali=0,valq=0,valu=0;
      long l3;
      for (l3=0;l3<ndata;l3++) {
        long ri3=htmp[l3].hrg-globalBeginRing*CUTRG;
        if (flg_rg[htmp[l3].ib][ri3/CUTRG]!=0) {
          long ir3=rgord2[htmp[l3].ib][ri3]+newnr2[htmp[l3].ib];
          vali+=htmp[l3].vi*x3[ir3];
          valq+=htmp[l3].vq*x3[ir3];
          valu+=htmp[l3].vu*x3[ir3];
        }
      }
      double qri=0;
      double qrq=0;
      double qru=0;
#endif
      for (l1=0;l1<ndata;l1++) {
        long ri1=htmp[l1].hrg-globalBeginRing*CUTRG;
        long iri1=rgord2[htmp[l1].ib][ri1]+newnr2[htmp[l1].ib];
        double CO1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].co
                                          -dpsisi[htmp[l1].ib]*htmp[l1].si);
        double SI1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].si
                                          +dpsisi[htmp[l1].ib]*htmp[l1].co);
        if (flg_rg[htmp[l1].ib][ri1/CUTRG]!=0) {
          double ww=htmp[l1].wp;
          double val2=x3[iri1]-(vali+CO1*valq+SI1*valu);
          
          qri-=ww*val2;
          qrq-=ww*val2*CO1;
          qru-=ww*val2*SI1;
          
          q2[iri1]+=ww*val2;
        }
      }
#ifdef TRUECP
#else

      for (l2=0;l2<ndata;l2++) {
        long ri2=htmp[l2].hrg-globalBeginRing*CUTRG;
        if (flg_rg[htmp[l2].ib][ri2/CUTRG]!=0) {
          long ir=rgord2[htmp[l2].ib][ri2]+newnr2[htmp[l2].ib]; 
          q2[ir]+=qri*htmp[l2].vi+qrq*htmp[l2].vq+qru*htmp[l2].vu;
        }
      }
#endif
      
#ifdef TIMING
      gettimeofday(&tp2,NULL);
      double dt=(double)(tp2.tv_sec-tp1.tv_sec)+(1E-6)*(tp2.tv_usec-tp1.tv_usec);
      dthit[ndata]+=dt;
      ndthit[ndata]+=1;
#endif
    }
  }

  if (rank==0) {
    double *lb = (double *) malloc(sizeof(double)*(nmatres));
    for (rrk=1;rrk<mpi_size;rrk++) {
      MPI_Recv(lb,sizeof(double)*(nmatres), MPI_BYTE, rrk,1032, MPI_COMM_WORLD,&statu);
      for (l=0;l<nmatres;l++) q2[l]+=lb[l];
    }
    free(lb);
  }
  else {
    MPI_Send(q2, sizeof(double)*(nmatres), MPI_BYTE, 0, 1032, MPI_COMM_WORLD);
  }

  if (rank==0) {
    for (i=0; i < nmatres; i++)
      {
        r2[i] = b2[i] - q2[i];
        d2[i] = r2[i] / hit2[i];
      }
  }

  
  double delta_new_tmp = 0.0;
  if (rank==0) {
    for (i = 0; i < nmatres; i++) {
      delta_new_tmp += b2[i] ;
      if (isnan(b2[i])) {
        fprintf(stderr,"NAN B2 PBS %ld\n",(long) i);
      }
      if (isnan(d2[i])) {
        fprintf(stderr,"NAN D2 PBS %ld %lg\n",(long) i,hit2[i]);
      }
    }
    //fprintf(stderr,"B2 %lg\n",delta_new_tmp);
  
    delta_new_tmp = 0.0;
    for (i = 0; i < nmatres; i++) {
      delta_new_tmp += q2[i] ;
      if (isnan(q2[i])) {
        fprintf(stderr,"NAN Q2 PBS %ld\n",(long) i);
      }
    }
    //fprintf(stderr,"Q2 %lg\n",delta_new_tmp);
    delta_new_tmp = 0.0;
    for (i = 0; i < nmatres; i++) {
      delta_new_tmp += q2[i]-b2[i] ;
      //fprintf(stderr,"B2 Q2 B2-Q2 [%ld]: %lg\t%lg\t%lg\n",(long) i,b2[i],q2[i],q2[i]-b2[i]);
    }
  }
  
  delta_new_tmp = 0.0;
  if (rank==0) for (i=0; i < nmatres; i++) {
    delta_new_tmp += r2[i] * d2[i];
  }
  
  delta_new=0;
  for (i=0;i<mpi_size;i++) {
    double tmp=delta_new_tmp;
    MPI_Bcast(&tmp, sizeof(double), MPI_BYTE, i, MPI_COMM_WORLD);
    delta_new+=tmp;
  }
  if (itbogo==0) delta0 = delta_new;
  if (rank==0) fprintf (stderr, "min_opt() iter = %d - delta0 = %lg - delta_new = %lg\n", iter, delta0, delta_new);
  int testwrit=0;

  MPI_Bcast(&delta0, sizeof(double), MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&delta_new, sizeof(double), MPI_BYTE, 0, MPI_COMM_WORLD);

  while ((iter < itermax)  && ((delta_new) > delta0*1E-24)) //Param->XI2STOP))
    {
      // q <= Ad
      if (rank==0) {
        memcpy(x2old,x3,nmatres*sizeof(double));
        testwrit=1;
      }
        //else testwrit=0;

      //healspline_fill_s_with_PtP_s (hs_rec, Ndata, vec, d, q);
      //for (i=0;i<mpi_size;i++) {
      //MPI_Bcast(d+tab_begr[i]*2, sizeof(double)*(tab_edr[i]-tab_begr[i]+1)*2, MPI_BYTE, i, MPI_COMM_WORLD);
      //}
      MPI_Bcast(d2, sizeof(double)*(nmatres), MPI_BYTE, 0, MPI_COMM_WORLD);
      
      // ===========================================================================
      // PROJECTION DE d dans q
      //
      
      for (l=0;l<nmatres;l++) q2[l]=0;
      if (rank==0) {
        double soff=0;
        for (i=0;i<newnr2[nbolo];i++) soff+=hit2[0]*d2[i];
        for (i=0;i<newnr2[nbolo];i++) q2[i]=soff;
      }

      for (k=0;k<nnbpix;k++) if (flgpix[k]>0) {
        
        //long imat=the_stat_pix[k];
        long ndata = loc_nhpix[k];
        if (ndata>0) {
          hpix *htmp = loc_hpix[k];

#ifdef TRUECP
#else
          double vali=0,valq=0,valu=0;
          long l3;
          for (l3=0;l3<ndata;l3++) {
            long ri3=htmp[l3].hrg-globalBeginRing*CUTRG;
            if (flg_rg[htmp[l3].ib][ri3/CUTRG]!=0) {
              long ir3=rgord2[htmp[l3].ib][ri3]+newnr2[htmp[l3].ib];
              vali+=htmp[l3].vi*d2[ir3];
              valq+=htmp[l3].vq*d2[ir3];
              valu+=htmp[l3].vu*d2[ir3];
            }
          }

          double qri=0;
          double qrq=0;
          double qru=0;
#endif
          for (l1=0;l1<ndata;l1++) {
            long ri1=htmp[l1].hrg-globalBeginRing*CUTRG;
            long iri1=rgord2[htmp[l1].ib][ri1]+newnr2[htmp[l1].ib];
            double CO1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].co
                                              -dpsisi[htmp[l1].ib]*htmp[l1].si);
            double SI1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].si
                                              +dpsisi[htmp[l1].ib]*htmp[l1].co);
            if (flg_rg[htmp[l1].ib][ri1/CUTRG]!=0) {
              double ww=htmp[l1].wp;
              double val2=d2[iri1]-(vali+CO1*valq+SI1*valu);
              qri-=ww*val2;
              qrq-=ww*val2*CO1;
              qru-=ww*val2*SI1;
              
              q2[iri1]+=ww*val2;
            }
          }
#ifdef TRUECP
#else
          for (l2=0;l2<ndata;l2++) {
            long ri2=htmp[l2].hrg-globalBeginRing*CUTRG;
            if (flg_rg[htmp[l2].ib][ri2/CUTRG]!=0) {
              long ir=rgord2[htmp[l2].ib][ri2]+newnr2[htmp[l2].ib]; 
              q2[ir]+=qri*htmp[l2].vi+qrq*htmp[l2].vq+qru*htmp[l2].vu;
            }
          }
        }
#endif
      }
      if (rank==0) {
        double *lb = (double *) malloc(sizeof(double)*(nmatres));
        for (rrk=1;rrk<mpi_size;rrk++) {
          MPI_Recv(lb,sizeof(double)*(nmatres), MPI_BYTE, rrk,1034, MPI_COMM_WORLD,&statu);
          for (l=0;l<nmatres;l++) q2[l]+=lb[l];
        }
        free(lb);
      }
      else {
        MPI_Send(q2, sizeof(double)*(nmatres), MPI_BYTE, 0, 1034, MPI_COMM_WORLD);
      }
      double dtq = 0.0;
      if (rank==0) {
        for (i=0; i < nmatres; i++) dtq += d2[i] * q2[i];
        alpha = delta_new / dtq;
        for (i=0; i < nmatres ; i++) x3[i] += alpha * d2[i];
#if 0
        char tmpfile[256];
        sprintf(tmpfile,"/data/delouis/TESTXI2/ITT_%d",(int) iter);
        FILE* fp=fopen(tmpfile,"w");
        fwrite(x3,nmatres*sizeof(double),1,fp);
        fclose(fp);
#endif
      }


      gettimeofday(&tp2,NULL);
      if (rank==0) fprintf (stderr,"min_opt() iter = %d - delta0 = %lg - delta_new = %lg %ld %3lfs\n", iter, delta0, delta_new,
                    (long) testwrit,(double)(tp2.tv_sec-tp1.tv_sec)+(1E-6)*(tp2.tv_usec-tp1.tv_usec));

      //if (rank==0) fprintf (stderr,".");
      gettimeofday(&tp1,NULL);
          
      if (iter % 100 == 0 && iter !=0)
        {
          // Use the best case
          memcpy(x3,x2old,nmatres*sizeof(double));
          //for (i=0;i<mpi_size;i++) {
          //  MPI_Bcast(x+tab_begr[i]*2, sizeof(double)*(tab_edr[i]-tab_begr[i]+1)*2, MPI_BYTE, i, MPI_COMM_WORLD);
          //}
          MPI_Bcast(x3, sizeof(double)*(nmatres), MPI_BYTE, 0, MPI_COMM_WORLD);

          
      // ===========================================================================
      // PROJECTION DE x dans q
      //
          
          for (l=0;l<nmatres;l++) q2[l]=0;
          double soff=0;
          for (i=0;i<newnr2[nbolo];i++) soff+=hit2[0]*x3[i];
          for (i=0;i<newnr2[nbolo];i++) q2[i]=soff;

          for (k=0;k<nnbpix;k++) if (flgpix[k]>0) {
        
            //long imat=the_stat_pix[k];
            long ndata = loc_nhpix[k];
            if (ndata>0) {
              hpix *htmp = loc_hpix[k];

#ifdef TRUECP
#else
              double vali=0,valq=0,valu=0;
              long l3;
              for (l3=0;l3<ndata;l3++) {
                long ri3=htmp[l3].hrg-globalBeginRing*CUTRG;
                if (flg_rg[htmp[l3].ib][ri3/CUTRG]!=0) {
                  long ir3=rgord2[htmp[l3].ib][ri3]+newnr2[htmp[l3].ib];
                  vali+=htmp[l3].vi*x3[ir3];
                  valq+=htmp[l3].vq*x3[ir3];
                  valu+=htmp[l3].vu*x3[ir3];
                }
              }
              double qri=0;
              double qrq=0;
              double qru=0;
#endif
              for (l1=0;l1<ndata;l1++) {
                long ri1=htmp[l1].hrg-globalBeginRing*CUTRG;
                long iri1=rgord2[htmp[l1].ib][ri1]+newnr2[htmp[l1].ib];
                double CO1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].co
                                                  -dpsisi[htmp[l1].ib]*htmp[l1].si);
                double SI1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].si
                                                  +dpsisi[htmp[l1].ib]*htmp[l1].co);
                if (flg_rg[htmp[l1].ib][ri1/CUTRG]!=0) {
                  double ww=htmp[l1].wp;
                  double val2=x3[iri1]-(vali+CO1*valq+SI1*valu);
                  qri-=ww*val2;
                  qrq-=ww*val2*CO1;
                  qru-=ww*val2*SI1;
                  
                  q2[iri1]+=ww*val2;
                }
              }
#ifdef TRUECP
#else
              for (l2=0;l2<ndata;l2++) {
                long ri2=htmp[l2].hrg-globalBeginRing*CUTRG;
                if (flg_rg[htmp[l2].ib][ri2/CUTRG]!=0) {
                  long ir=rgord2[htmp[l2].ib][ri2]+newnr2[htmp[l2].ib]; 
                  q2[ir]+=qri*htmp[l2].vi+qrq*htmp[l2].vq+qru*htmp[l2].vu;
                }
              }
#endif
              
            }
          }
          if (rank==0) {
            double *lb = (double *) malloc(sizeof(double)*(nmatres));
            for (rrk=1;rrk<mpi_size;rrk++) {
              MPI_Recv(lb,sizeof(double)*(nmatres), MPI_BYTE, rrk,1034, MPI_COMM_WORLD,&statu);
              for (l=0;l<nmatres;l++) q2[l]+=lb[l];
            }

            free(lb);
          }
          else {
            MPI_Send(q2, sizeof(double)*(nmatres), MPI_BYTE, 0, 1034, MPI_COMM_WORLD);
          }
          
          if (rank==0) for (i=0; i < nmatres; i++) {
            r2[i] = b2[i] - q2[i];
          }
          
        }
      else
        {
          if (rank==0) for (i=0; i < nmatres; i++) r2[i] -= alpha * q2[i];
        }
      
      if (rank==0) for (i=0; i < nmatres; i++) s2[i] = r2[i] / hit2[i];

      delta_old = delta_new;
      if (rank==0) {
        delta_new=0;
        for (i=0; i < nmatres ; i++) delta_new += r2[i] * s2[i];
        beta = delta_new / delta_old;
        for (i=0; i < nmatres ; i++) d2[i] = s2[i] + beta * d2[i];
      }
      iter ++;
      MPI_Bcast(&delta_new, sizeof(double), MPI_BYTE, 0, MPI_COMM_WORLD);
      
    }
  
  if (rank==0) fprintf (stderr,"\nmin_opt() iter = %d - delta0 = %lg - delta_new = %lg\n", 
                        iter, delta0, delta_new);
  if (rank==0) fprintf (stderr,"CG in iter = %d (max=%d)\n", iter, itermax);
  
  if (rank==0) memcpy(x3  ,x2old,nmatres*sizeof (double));
  MPI_Bcast(x3, sizeof(double)*(nmatres), MPI_BYTE, 0, MPI_COMM_WORLD);
  itbogo++;
}



double funcxigain(double *p)
{
  long i;
  int rank,size;
  //int mpi_size;

  int rank_size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank_size);
  rank=rank_size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  //mpi_size=size;

  double *x3= (double *) malloc(sizeof(double)*(newnr[nbolo]+nbolo*npixbeam+nbolo*(GAINSTEP)+nmatco+nmatdust+nfreefree));

  minimize(x2,p);

  for (i=0;i<newnr[nbolo];i++)   x3[i]=x2[i];
  for (i=0;i<nbolo*GAINSTEP;i++) {
    x3[newnr[nbolo]+i]=p[i];
  }
    
  for (i=0;i<nbolo*npixbeam+nmatco+nmatdust+nfreefree;i++) {
    x3[newnr[nbolo]+nbolo*(GAINSTEP)+i]=x2[newnr[nbolo]+i];
  }
  double xi2_res=funcxi2(x3);
  if (rank==0) fprintf(stderr,"XI2 %.10lg\n",xi2_res);
  free(x3);
  return(xi2_res);
}

void dfuncxigain(double *p,double *xi)
{
  long i;
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  double *x3= (double *) malloc(sizeof(double)*(newnr[nbolo]+nbolo*npixbeam+nbolo*(GAINSTEP)+nmatco+nmatdust+nfreefree));
  double *xi2= (double *) malloc(sizeof(double)*(newnr[nbolo]+nbolo*npixbeam+nbolo*(GAINSTEP)+nmatco+nmatdust+nfreefree));
   minimize(x2,p);
  for (i=0;i<newnr[nbolo];i++)   x3[i]=x2[i];
  for (i=0;i<nbolo*GAINSTEP;i++) x3[newnr[nbolo]+i]=p[i];
  for (i=0;i<nbolo*npixbeam+nmatco+nmatdust+nfreefree;i++) {
    x3[newnr[nbolo]+nbolo*(GAINSTEP)+i]=x2[newnr[nbolo]+i];
    if (rank==0) fprintf(stderr,"X2[%d] = %lg\n",(int) i,x2[newnr[nbolo]+i]);
  }
  

  dfuncxi2_light(x3,xi2);

  //double xi2_res=funcxi2(x3);
  for (i=0;i<nbolo*GAINSTEP;i++) {
    //x3[newnr[nbolo]+i]=p[i]+1E-7;
    //double xi2_res_plus=funcxi2(x3);
    //xi[i]=(xi2_res_plus-xi2_res)/1E-7;
    //x3[newnr[nbolo]+i]=p[i];
    //if (rank==0) fprintf(stderr,"DXI2[%d] =%lg %lg\n",(int) i,xi2[newnr[nbolo]+i],xi[i]);
    xi[i]=xi2[newnr[nbolo]+i]; // ?????? Je ne sais pas encore ....
  }

  //double total=0;
  //for (i=0;i<nbolo*GAINSTEP;i++) {
  //  total+=xi[i]*xi[i];
  //}
  //if (total>.0009*nbolo*GAINSTEP) {
  //  total=0.03;
  //  for (i=0;i<nbolo*GAINSTEP;i++) xi[i]/=0.03;
  //}


  free(xi2);
  free(x3);
}


int PIOMergeMAP(const char *path)
{
  //int rank;
  int mpi_size;

  int rank_size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank_size);
  //rank=rank_size;
  MPI_Comm_size(MPI_COMM_WORLD,&rank_size);
  mpi_size=rank_size;

  int err=0;
  double *map = (double *) malloc(sizeof(double)*12*2048*2048);
  
  char thepath[1024];
  
  fprintf(stderr,"Write %s\n",path);
  sprintf(thepath,"%s_dir_%d",path,mpi_size);
          

  int rrk;
  long nnn=0;


  for (rrk=0;rrk<mpi_size;rrk++) {
    struct stat buf;
    sprintf(thepath,"%s_dir_%d/proc_%d",path,mpi_size,rrk);
    
    stat(thepath,&buf);
    long sizefile = buf.st_size;
    unsigned char *value = (unsigned char *) malloc(sizefile);
    
    int fp=open(thepath,O_WRONLY|O_CREAT,0664);
    err=read(fp,value,sizefile);
    close(fp);
    
    memcpy(map+nnn,value,sizefile);
    fprintf(stderr,"nnn %lg\n",map[nnn]);
    nnn+=sizefile/sizeof(double);
  }
  
  int fp=open(path,O_WRONLY|O_CREAT,0664);
  err=write(fp,map,12*2048*2048*sizeof(double));
  close(fp);
  free(map);

  return(err);
}


int PIOWriteMAP(const char *path, double *value_in_double,int beg,int end)
{
  int rank,mpi_size;
  int map_size = 12*2048*2048;

  int rank_size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank_size);
  rank=rank_size;
  MPI_Comm_size(MPI_COMM_WORLD,&rank_size);
  mpi_size=rank_size;

  ssize_t err=0;
  // Convert array from double to float for smaller file to be written
  float *value = (float *)malloc(sizeof(float)*(end-beg+1));
  for (int i = 0; i < (end-beg+1); ++i) {
    value[i] = (float)value_in_double[i];
  }

  if (rank==0) {
    PIOSTRING fitspath;
    sprintf( fitspath, "%s.fits", path);
    MPI_Status statu;
    float *map=(float *)malloc(sizeof(float)*map_size);
    int i,rrk;
    for (i=beg;i<=end;i++) map[i]=value[i-beg];
    for (rrk=1;rrk<mpi_size;rrk++) {
      int l_beg,l_end;
      MPI_Recv(&l_beg,sizeof(int), MPI_BYTE, rrk,2031, MPI_COMM_WORLD,&statu);
      MPI_Recv(&l_end,sizeof(int), MPI_BYTE, rrk,2032, MPI_COMM_WORLD,&statu);
      MPI_Recv(map+l_beg,(l_end-l_beg+1)*sizeof(float), MPI_BYTE, rrk,2033, MPI_COMM_WORLD,&statu);
    }
    if (remove( fitspath) == 0) {
      fprintf( stderr, "removed existing MAP: %s\n", fitspath);
    }
    fprintf( stderr, "WRITE MAP: %s\n", fitspath);
    write_healpix_map( map, 2048, fitspath, 0, "G");
    free(map);
  }

  else {
    MPI_Send(&beg, sizeof(int), MPI_BYTE, 0, 2031, MPI_COMM_WORLD);
    MPI_Send(&end, sizeof(int), MPI_BYTE, 0, 2032, MPI_COMM_WORLD);
    MPI_Send(value, (end-beg+1)*sizeof(float), MPI_BYTE, 0, 2033, MPI_COMM_WORLD);
  }

  free(value);
  return(err);
}

int PIOWriteVECT(const char *path,void *value,int off,int size)
{
  int fp=open(path,O_WRONLY|O_CREAT,0664);
  int err=pwrite(fp,value,size,off);
  close(fp);
  return(err);
}



int main(int argc,char *argv[])  {

  PIOSTRING Command;
  PIOSTRING hostname;
  int rank,size,mpi_size;
  long i,j;
  PIOLONG rg;
  int stim_first_seed = 0;
  time_t now;

  //PIODOUBLE  dip[]={-0.000233797238224 , -0.00222070369388 , 0.00250271842609};

  //TESPT;
  MPI_Init(&argc, &argv);

  int rank_size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank_size);
  rank=rank_size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  mpi_size=size;

  if (rank==0) {
    char rpath[PATH_MAX];
    realpath( argv[0], rpath);
    now = time( NULL);
    fprintf( stderr, "%s: --------------------------\n", __FILE__ );
    fprintf( stderr, "%s: Starting %s\n",                __FILE__, rpath);
    fprintf( stderr, "%s: at %s",                        __FILE__, ctime( &now));
    fprintf( stderr, "%s: with %d MPI ranks and %d OMP threads\n", __FILE__, mpi_size, omp_get_max_threads());
    fprintf( stderr, "%s: --------------------------\n", __FILE__ );

    /* Ensure mpi_rank is a power of two (required by some algo...) */
    if (!isPowerOfTwo(mpi_size)) {
      fprintf(stderr,"INVALID mpi_size(%d), shall be a power of two!\n", mpi_size);
      return 1;
    }
    fprintf( stderr, "\nstarting sroll with MAXSIMU=%d, MAXTHEOHPR=%d\n", MAXSIMU, MAXTHEOHPR);
  }

  /* read params */
  sroll_parContent par;
  GetHostname( hostname);

  int res = sroll_readParam(&par, argv[1] );
  if (res != 0) {
    fprintf(stderr, "Unable to parse the parameter file.\n");
    exit(res);
  }

  //TESPT;
  Param = &par;

  /*-------------------------------------------------------------------------*/
  /* parameters consistency checks                                           */
  /* and default values / legacy behavior for optional parameters            */
  /*-------------------------------------------------------------------------*/

  // global number of bolometers
  nbolo = Param->n_Ptg_noPS;

  // get the list of pixnames from Ptg_noPS objects
  DETNAME pixnames[nbolo];
  for (i=0; i<nbolo; i++) {
    strncpy( pixnames[i], get_pixname( Param->Ptg_noPS[i]), DETNAMELEN);
  }

  // bolometers frequency
  int freq = strtol( pixnames[0], NULL, 10);
  assert( freq != 0);

  assert( Param->n_OUT_NOPOL == Param->n_MAP);

  if (Param->flag_KCMBIN == 0) {
    Param->KCMBIN = 0;
  }

  // default value for REMHDIP is 1 for 545GHz-857GHz
  if ((freq >= 545) && (Param->flag_REMHDIP == 0)) {
    Param->REMHDIP = 1;
    Param->flag_REMHDIP = 1;
  }

  if (Param->n_bolomask == 0) {
    // if bolomask is empty, produce a map with all bolometers
    assert( Param->n_MAP == 1);
    Param->n_bolomask = nbolo;
    Param->bolomask = malloc( nbolo * sizeof( PIOINT));
    assert( Param->bolomask != NULL);
    for (i=0; i<nbolo; i++) {
      Param->bolomask[i] = 1;
    }
  }

  if (Param->flag_MAPRINGS == 0) {
    // if MAPRINGS not given, set its default value according to DOMAXVRAIE
    Param->flag_MAPRINGS = 1;
    Param->n_MAPRINGS = Param->n_MAP;
    Param->MAPRINGS = malloc( Param->n_MAPRINGS * sizeof( PIOLONG));
    assert( Param->MAPRINGS != NULL);
    for (i=0; i<Param->n_MAPRINGS; i++) {
      if (Param->DOMAXVRAIE == 0) {
        Param->MAPRINGS[i] = FULL + HM12 + S12345 + YEAR12 + FULLODDEVEN;
      }
      else {
        Param->MAPRINGS[i] = FULL + HM12;
      }
    }
  }
  assert( Param->n_MAPRINGS == Param->n_MAP);

  // manage addHPR lists and default values
  if (Param->flag_addHPR_name == 1) {
    assert( Param->n_addHPR_name % nbolo == 0);
    // if addHPR_factor is not given, set it to 1.0
    if (Param->flag_addHPR_factor == 0) {
      Param->addHPR_factor = malloc( sizeof( PIOFLOAT));
      if (Param->addHPR_factor == NULL) {
        fprintf( stderr, "ERROR: not enough memory to allocate Param->addHPR_factor\n");
        exit(-1);
      }
      Param->addHPR_factor[0] = 1.0;
      Param->n_addHPR_factor = 1;
      Param->flag_addHPR_factor = 1;
    }
    if ((Param->n_addHPR_factor == 1) && (Param->n_addHPR_name > 1)) {
      // if only one addHPR_factor value is given, use it for all addHPR_name objects
      assert( realloc( Param->addHPR_factor, Param->n_addHPR_name * sizeof( PIOFLOAT)) != NULL);
      for (i=1; i<Param->n_addHPR_name; i++) {
        Param->addHPR_factor[i] = Param->addHPR_factor[0];
      }
      Param->n_addHPR_factor = Param->n_addHPR_name;
    }
    // if addHPR_watts is not given, set it to 0
    if (Param->flag_addHPR_watts == 0) {
      Param->addHPR_watts = malloc( sizeof( PIOINT));
      if (Param->addHPR_watts == NULL) {
        fprintf( stderr, "ERROR: not enough memory to allocate Param->addHPR_watts\n");
        exit(-1);
      }
      Param->addHPR_watts[0] = 0;
      Param->n_addHPR_watts = 1;
      Param->flag_addHPR_watts = 1;
    }
    if ((Param->n_addHPR_watts == 1) && (Param->n_addHPR_name > 1)) {
      // if only one addHPR_watts value is given, use it for all addHPR_name objects
      assert( realloc( Param->addHPR_watts, Param->n_addHPR_name * sizeof( PIOINT)) != NULL);
      for (i=1; i<Param->n_addHPR_name; i++) {
        Param->addHPR_watts[i] = Param->addHPR_watts[0];
      }
      Param->n_addHPR_watts = Param->n_addHPR_name;
    }
  }


  /*-------------------------------------------------------------------------*/
  /* MPI: Ring dispatching between available ranks                           */
  /*-------------------------------------------------------------------------*/
  globalBeginRing = Param->BeginRing;
  globalEndRing   = Param->EndRing;

  /* Compute some frequently used values */
  globalRangeRing = globalEndRing - globalBeginRing + 1;
  globalRankInfo.BeginRing = (PIOLONG *)malloc(sizeof(PIOLONG) * mpi_size);
  if (globalRankInfo.BeginRing == NULL) {
    perror("Error");
    return 1;
  }
  globalRankInfo.EndRing   = (PIOLONG *)malloc(sizeof(PIOLONG) * mpi_size);
  if (globalRankInfo.EndRing == NULL) {
    perror("Error");
    return 1;
  }

  if (Param->flag_stim_paramfiles == 1) {
    /* Compute load balancing per sample count, starting from the end where very long rings are */
    for (int irank = mpi_size-1; irank >= 0 ; irank--) {
      if (irank == mpi_size-1) {
        globalRankInfo.EndRing[irank] = globalEndRing;
      } else {
        globalRankInfo.EndRing[irank] = globalRankInfo.BeginRing[irank+1] - 1;
      }
      if (irank == 0) {
        globalRankInfo.BeginRing[irank] = globalBeginRing;
      } else {
        long samples_to_process = (ENDRINGINDEX( globalRankInfo.EndRing[irank]) - BEGINRINGINDEX( globalBeginRing)) / (irank+1);
        int temp_beg_ring = globalRankInfo.EndRing[irank]-1;
        while (ENDRINGINDEX( globalRankInfo.EndRing[irank]) - BEGINRINGINDEX( temp_beg_ring) < samples_to_process) {
          temp_beg_ring--;
        }
        globalRankInfo.BeginRing[irank] = temp_beg_ring+1;
      }
    }
#if 0
    /* Compute load balancing per sample count*/
    for (int irank = 0; irank < mpi_size; irank++) {
      if (irank > 0) {
        globalRankInfo.BeginRing[irank] = globalRankInfo.EndRing[irank-1] + 1;
      } else {
        globalRankInfo.BeginRing[irank] = globalBeginRing;
      }
      if (irank == mpi_size - 1) {
        globalRankInfo.EndRing[irank] = globalEndRing;
      } else {
        long samples_to_process = (ENDRINGINDEX( globalEndRing) - BEGINRINGINDEX( globalRankInfo.BeginRing[irank])) / (mpi_size - irank);
        int temp_end_ring = globalRankInfo.BeginRing[irank];
        while (ENDRINGINDEX( temp_end_ring) - BEGINRINGINDEX( globalRankInfo.BeginRing[irank]) < samples_to_process) {
          temp_end_ring++;
        }
        globalRankInfo.EndRing[irank] = temp_end_ring;
      }
    }
#endif
  } else {
    /* Compute load balancing per ring count*/
    // number of ranks to get an extra ring to proceed
    PIOLONG balancing_correction = 0;
    PIOLONG rings_per_rank = globalRangeRing / mpi_size;
    /* Check border case: when more proc than ring to process */
    if (rings_per_rank == 0) {
      if (rank==0) {
        fprintf(stderr,"ERROR: too few data to be processed ("PIOLONG_FMT" rings) regarding the available ranks (%d)\n",
            globalRangeRing, mpi_size);
      }
      return 1;
    } else { // Take all available procs
      balancing_correction = globalRangeRing - (rings_per_rank * mpi_size);
    }
  
    if (rank==0) {
      fprintf( stderr, "LoadBalancing globalRangeRing      = "PIOLONG_FMT" \n", globalRangeRing);
      fprintf( stderr, "LoadBalancing rings_per_rank       = "PIOLONG_FMT" \n", rings_per_rank);
      fprintf( stderr, "LoadBalancing balancing_correction = "PIOLONG_FMT" \n", balancing_correction);
    }
  
    PIOLONG previous = globalBeginRing;
    for (int irank = 0; irank < mpi_size; irank++) {
      PIOLONG sup = (balancing_correction>0) ? 1:0;
      PIOLONG range = rings_per_rank + sup;
      globalRankInfo.BeginRing[irank] = previous;
      globalRankInfo.EndRing[irank]   = previous+range-1; // -1 since this is an index of ring!
      previous += range;
      if (balancing_correction != 0) {
        balancing_correction--;
      }
    }
  }

  /* Display ring dispatching between ranks */
  for (int irank = 0; irank < mpi_size; irank++) {
    if (rank == irank) {
      fprintf( stderr, "  rank#%d/%d (%s, %.2fGB) begin=%ld end=%ld, (%ld rings, %de6 samples)\n",
              irank, mpi_size, hostname, GetFreeMemGB(),
              globalRankInfo.BeginRing[irank], globalRankInfo.EndRing[irank],
              globalRankInfo.EndRing[irank] - globalRankInfo.BeginRing[irank] + 1,
              (int)((ENDRINGINDEX( globalRankInfo.EndRing[irank]) - BEGINRINGINDEX(globalRankInfo.BeginRing[irank])) / 1e6));
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  /*--- End : MPI: Ring dispatching --- */


  PIOINT Nside=Param->Nside;

  if (Param->flag_stim_paramfiles == 1) {
    assert( Param->n_stim_paramfiles == nbolo);
  }

#ifndef DOMAP
  PIOLONG LMAX=Param->LMAX*(Param->LMAX+1)+Param->LMAX;
#endif
  MPI_Barrier(MPI_COMM_WORLD);
  CUTRG=Param->CUTRG;

  PIOSTRING *mapout = (PIOSTRING *) malloc(sizeof(PIOSTRING)*Param->n_MAP);

  REMHDIP = Param->REMHDIP;
  Param->REMDIP = Param->REMDIP;

  double avvnep=0;
  for (i=0;i<nbolo;i++) avvnep+=Param->NEP[i]/Param->Calibration[i];
  for (i=0;i<nbolo;i++) NEP_tab[i]=nbolo*Param->NEP[i]/Param->Calibration[i]/avvnep;

  MPI_Barrier(MPI_COMM_WORLD);

  PIOINT  *badring;
#ifdef GAIN_RATIO
  PIODOUBLE  **gain_ratio = (PIODOUBLE **) malloc(nbolo*sizeof(PIODOUBLE *));
  PIOLONG     **gain_ratio_off = (PIOLONG **) malloc(nbolo*sizeof(PIOLONG *));
#endif
  
  srand48((long)(rank+Param->SEED[0]));

  PIOLONG *l_nhpix  = (PIOLONG *) malloc(sizeof(PIOLONG)*12*Nside*Nside);
  hpix **l_hpix = (hpix **) malloc(sizeof(hpix *)*12*Nside*Nside);
  memset(l_nhpix,0,sizeof(PIOLONG)*12*Nside*Nside);

  assert( Param->n_Theo_noPS % nbolo == 0);
  npixbeam = Param->n_Theo_noPS / nbolo;
  assert( npixbeam <= MAXTHEOHPR);
  
  /*======================================================
    =
    =      read data
    =
    =*/

  long step=Param->RSTEP;


  PIOFLOAT **theo = (PIOFLOAT **) malloc(sizeof(PIOFLOAT *)*npixbeam);

  long vmem,phymem;

  PIOLONG ib;
  double *eta=(double *) malloc(sizeof(double)*nbolo); // (1-crosspol)/(1+crosspol)
  eta_dest=(double *) malloc(sizeof(double)*nbolo);
  dpsico=(double *) malloc(sizeof(double)*nbolo);
  dpsisi=(double *) malloc(sizeof(double)*nbolo);

  PIOFLOAT *FREEFREE=NULL;
  if (Param->flag_Theo_FREEFREE==_PAR_TRUE) {
    FREEFREE = (PIOFLOAT *) malloc(sizeof(PIOFLOAT)*12*128*128);
    PIOLONG nsa  = noDMC_readObject_PIOFLOAT(Param->Theo_FREEFREE,0,12*128*128,FREEFREE);
    if (nsa<0) {
      fprintf(stderr, "Impossible to read Theo_FREEFREE: %s %d\n",Param->Theo_FREEFREE,(int) nsa);
      
      exit ( -1);
    }
  } 

  PIOLONG nadu3=Param->NADU;
  GAINSTEPADU=Param->NADU;
  long DODISTOR=0;

  if (GAINSTEPADU>1) {
    DODISTOR=Param->NADU;
    GAINSTEPADU=Param->GAINSTEP;
    nadu3=GAINSTEPADU;
    Param->GAINSTEP=1;
  }
  if (rank==0) fprintf(stderr,"DODISTOR : %d\n",(int) DODISTOR);
  if (rank==0) fprintf(stderr,"REMDIP : %d\n",(int) Param->REMDIP);
  
  PIOFLOAT *comapI=NULL;
  if (Param->flag_Theo_CO==_PAR_TRUE) {
    comapI = (PIOFLOAT *) malloc(sizeof(PIOFLOAT)*12*128*128);
    PIOLONG nsa  = noDMC_readObject_PIOFLOAT(Param->Theo_CO,0,12*128*128,comapI);
    if (nsa<0) {
      fprintf(stderr, "Impossible to read Theo_CO: %s %d\n",Param->Theo_CO,(int) nsa);
      
      exit ( -1);
    }
  } 

  PIOFLOAT *dustmapI=NULL;
  PIOFLOAT *dustmapQ=NULL;
  PIOFLOAT *dustmapU=NULL;
  if (Param->flag_Theo_Dust_I==_PAR_TRUE) {
    dustmapI = (PIOFLOAT *) malloc(sizeof(PIOFLOAT)*12*128*128);
    PIOLONG nsa  = noDMC_readObject_PIOFLOAT(Param->Theo_Dust_I,0,12*128*128,dustmapI);
    if (nsa<0) {
      fprintf(stderr, "Impossible to read Theo_Dust_I: %s %d\n",Param->Theo_Dust_I,(int) nsa);
      
      exit ( -1);
    }
    dustmapQ = (PIOFLOAT *) malloc(sizeof(PIOFLOAT)*12*128*128);
    nsa  = noDMC_readObject_PIOFLOAT(Param->Theo_Dust_Q,0,12*128*128,dustmapQ);
    if (nsa<0) {
      fprintf(stderr, "Impossible to read Theo_Dust_Q: %s %d\n",Param->Theo_Dust_Q,(int) nsa);
      
      exit ( -1);
    }
    dustmapU = (PIOFLOAT *) malloc(sizeof(PIOFLOAT)*12*128*128);
    nsa  = noDMC_readObject_PIOFLOAT(Param->Theo_Dust_U,0,12*128*128,dustmapU);
    if (nsa<0) {
      fprintf(stderr, "Impossible to read Theo_Dust_U:%s %d\n",Param->Theo_Dust_U,(int) nsa);
      
      exit ( -1);
    }
  } 

  PIOFLOAT *simdustmapQ=NULL;
  PIOFLOAT *simdustmapU=NULL;

  double *histo_adu=(double *) malloc(NADUSTEP*sizeof(double)*nbolo);
  int *invadu = (int *) malloc(NADUSTEP*sizeof(int)*nbolo);
  xadu = (double *) malloc(GAINSTEPADU*sizeof(double)*nbolo);
  double *yadu = (double *) malloc(GAINSTEPADU*sizeof(double)*nbolo);
  double *nxadu = (double *) malloc(GAINSTEPADU*sizeof(double)*nbolo);

  long nnoisesim = 1048576;
  fftw_complex  *in_fft = NULL;
  fftw_complex  *out_fft = NULL;
  fftw_plan     p1 = NULL, p2 = NULL;
  
  if (Param->flag_stim_paramfiles == 1) {
    // stim_paramfiles parameter overrides TESTPOL
    Param->TESTPOL = -1;
  }

  if (Param->TESTPOL==0) {
    assert( (Param->n_Calibration == 5 * nbolo)
            && "Error: simulations inside sroll requires 5 calibration values per bolo for ADCNL residuals");
    // prepare FFTW for adding noise from noise Fourier Transform
    in_fft  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nnoisesim);
    out_fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nnoisesim);
    p1 = fftw_plan_dft_1d( nnoisesim, in_fft, out_fft, FFTW_FORWARD, FFTW_ESTIMATE );
    p2 = fftw_plan_dft_1d( nnoisesim, in_fft, out_fft, FFTW_BACKWARD, FFTW_ESTIMATE );
  }
  
  //long stepadu=(globalEndRing-globalBeginRing+1)/Param->REMDIP;
  int number_of_iterations = 1;

  for (ib=0;ib<nbolo;ib++) {

    double *tfnoise=NULL;
    if (Param->TESTPOL==0) {
      tfnoise=(double *)malloc(sizeof(double)*nnoisesim);
      int fp=open(Param->Signal_noPS[2*nbolo+ib],O_RDONLY,0664);
      int err=read(fp,tfnoise,sizeof(double)*nnoisesim);
      close(fp);
      if (err<0) {
        fprintf(stderr,"Do not succeed to read Signal_noPS[%ld]: %s %d %ld\n",2*nbolo+ib,Param->Signal_noPS[2*nbolo+ib],err,(long) nnoisesim);
        exit(-1);
      }
      double reg=sqrt(1.21)/((double)nnoisesim); //Take into account effect of the correlation
      for (i=0;i<nnoisesim;i++) tfnoise[i]*=reg;
    }

    eta[ib]=(1.-Param->CrossPol[ib])/(1.+Param->CrossPol[ib]);
    if (Param->D_NOPOL) eta_dest[ib]=0;
    else eta_dest[ib]=eta[ib];

    dpsico[ib]=1;
    dpsisi[ib]=0;
    double sxi= (Param->Calibration[ib]/Param->NEP[ib])*(Param->Calibration[ib]/Param->NEP[ib]);
    if (rank==0) fprintf(stderr,"SXI %d %lg\n",(int) ib,sxi);
    sprintf(Command,"begin=%lld;end=%lld",
            (long long) (globalRankInfo.BeginRing[rank]),
            (long long) (globalRankInfo.EndRing[rank]));

    int ring_count = globalRankInfo.EndRing[rank]-globalRankInfo.BeginRing[rank]+1;
    badring = (PIOINT *) malloc(sizeof(PIOINT)*ring_count);
    if (badring == NULL) {
      fprintf(stderr, "**** ERROR **** Memory allocation error !!!!\n");
      fprintf(stderr, "** rank = %d  globalRankInfo.BeginRing[rank]="PIOLONG_FMT"  globalRankInfo.EndRing[rank]="PIOLONG_FMT"\n", rank, globalRankInfo.BeginRing[rank], globalRankInfo.EndRing[rank]);
    }

    assert( check_pixname( Param->Badring[ib], pixnames[ib]));
    int tperr = noDMC_readObject_PIOINT(Param->Badring[ib],globalRankInfo.BeginRing[rank],ring_count,badring);
    if (tperr<0) {
      fprintf(stderr, "Impossible to read Badring[%ld]: %s %d\n",ib,Param->Badring[ib],tperr);
      exit ( -1);
    }


    int iter;
    PIOFLOAT *stim_hpr[MAXSIMU];
    for (iter = 0; iter < number_of_iterations; iter++) { 
      stim_hpr[iter] = NULL;
    }

    if (Param->flag_stim_paramfiles == 1) {
      
      if (rank==0) {
        GetProcMem(&vmem,&phymem);
        fprintf(stderr,"\nbefore stim bolo loop: used VMEM %.1lf[PHYS %.1lf]MB\n",
            (double) vmem/1024./1024., (double) phymem/1024./1024.);
      }

      stimParameters stimPar;
      // read stim parameter file for this bolometer and broadcast it to all MPI ranks
      init_stim( &stimPar, Param->stim_paramfiles[ib]);
      number_of_iterations = stimPar.Param.iterations;
      if (number_of_iterations > MAXSIMU) {
        fprintf(stderr, "Number of simulation iterations %d > MAXSIMU %d\n", number_of_iterations, MAXSIMU);
        exit ( -1);
      }
      if ((number_of_iterations > 1) && (stimPar.Param.stay_in_memory == -1)) {
        stimPar.Param.stay_in_memory = 1;
      } else {
        stimPar.Param.stay_in_memory = 0;
      }
      stim_first_seed = stimPar.Param.random_seed;

      // iterate on stim
      for (iter = 0; iter < number_of_iterations; iter++) { 
        stim_hpr[iter] = malloc( ring_count * RINGSIZE * sizeof( PIOFLOAT));
        assert( stim_hpr[iter] != NULL);
        stim( &stimPar, globalRankInfo.BeginRing[rank], ring_count, stim_hpr[iter]);
        stimPar.Param.random_seed++;
      }
      free_stim( &stimPar);

      if (rank==0) {
        GetProcMem(&vmem,&phymem);
        fprintf(stderr,"\nafter stim bolo loop: used VMEM %.1lf[PHYS %.1lf]MB\n",
            (double) vmem/1024./1024., (double) phymem/1024./1024.);
      }

    }

    for (rg=globalRankInfo.BeginRing[rank];rg<=globalRankInfo.EndRing[rank];rg+=step) {

      if (badring[rg-globalRankInfo.BeginRing[rank]]==0) {

        PIOFLOAT *h,*fsl; 
        PIODOUBLE *ph,*th,*psi,*phase;
        PIOFLOAT *diporb;
        PIOBYTE surv=-1;
        if (rg>=240&&rg<=5720)    surv=1;
        if (rg>=5721&&rg<=11194)  surv=2;
        if (rg>=11195&&rg<=16691) surv=3;
        if (rg>=16692&&rg<=21720) surv=4;
        if (rg>=21721)            surv=5;
        if (Param->EndRing == 27005) {
          // DX11 / RD12ll half-mission
          if (rg>13471)           surv+=10;
        }
        else {
          // RD12RC* half-mission
          if (rg>13144)           surv+=10;
        }
          
        PIOFLOAT *y[MAXSIMU];
        // get rid of gcc "maybe-uninitialized" warning depending on optimsation level
        for (iter = 0; iter < number_of_iterations; iter++) { 
          y[iter] = NULL;
        }
          
        sprintf(Command,"ring=%lld",(long long)(rg));
        assert( check_pixname( Param->Hit_noPS[ib], pixnames[ib]));
        h = (PIOFLOAT *) malloc(sizeof(PIOFLOAT)*RINGSIZE);
        assert( (h != NULL) && "Not enough memory to allocate hitcount HPR array");
        tperr = noDMC_readObject_PIOFLOAT(Param->Hit_noPS[ib],rg*RINGSIZE,RINGSIZE,h);
        if (tperr<0) {
          fprintf(stderr, "Impossible to read Hit_noPS[%ld]: %s %d\n",ib,Param->Hit_noPS[ib],tperr);          exit ( -1);
        }
  
        if (Param->flag_fsl == _PAR_TRUE) {
          assert( check_pixname( Param->fsl[ib], pixnames[ib]));
          fsl = (PIOFLOAT *) malloc(sizeof(PIOFLOAT)*RINGSIZE);
          assert( (fsl != NULL) && "Not enough memory to allocate FSL HPR array");
          tperr = noDMC_readObject_PIOFLOAT(Param->fsl[ib],rg*RINGSIZE,RINGSIZE,fsl);
          if (tperr<0) {
            fprintf(stderr, "Impossible to read fsl[%ld]: %s %d\n",ib,Param->fsl[ib],tperr);
            exit ( -1);
          }
        } else {
          fsl = NULL;
        }
  
        if (Param->flag_stim_paramfiles == 0) {
          assert( check_pixname( Param->Signal_noPS[ib], pixnames[ib]));
          y[0] = (PIOFLOAT *) malloc(sizeof(PIOFLOAT)*RINGSIZE);
          assert( (y[0] != NULL) && "Not enough memory to allocate signal HPR array");
          tperr = noDMC_readObject_PIOFLOAT(Param->Signal_noPS[ib],rg*RINGSIZE,RINGSIZE,y[0]);
          if (tperr<0) {
            fprintf(stderr, "Impossible to read Signal_noPS[%ld]: %s %d\n",ib,Param->Signal_noPS[ib],tperr);
            exit ( -1);
          }
        } else {
          for (iter = 0; iter < number_of_iterations; iter++) {
            y[iter] = stim_hpr[iter] + RINGSIZE * (rg - globalRankInfo.BeginRing[rank]);
          }
        }

  
        /*========================================================================
          =
          =                  AJOUT TOISIM
          =*/
  
        assert( check_pixname( Param->DipOrb_noPS[ib], pixnames[ib]));
        diporb = (PIOFLOAT *) malloc(sizeof(PIOFLOAT)*RINGSIZE);
        assert( (diporb != NULL) && "Not enough memory to allocate total dipole HPR array");
        tperr = noDMC_readObject_PIOFLOAT(Param->DipOrb_noPS[ib],rg*RINGSIZE,RINGSIZE,diporb);
        if (tperr<0) {
          fprintf(stderr, "Impossible to read DipOrb_noPS[%ld]: %s %d\n",ib,Param->DipOrb_noPS[ib],tperr);
          exit ( -1);
        }

        if (Param->KCMBIN == 1) {
          // input HPR are in KCMB, convert them to watts
          for (i=0; i<RINGSIZE; i++) {
            y[0][i] *= Param->Calibration[ib];
          }
        }

        if (Param->ADDDIP == 1) {
          // for projection only, signal must contain the total dipole, add it if it's not here
          for (i=0; i<RINGSIZE; i++) {
            y[0][i] = y[0][i] + diporb[i] * Param->Calibration[ib];
          }
        }



        if (Param->flag_addHPR_name == 1) {
          PIOFLOAT *addhpr = (PIOFLOAT *) malloc( sizeof( PIOFLOAT) * RINGSIZE);
          for (int i = 0; i < Param->n_addHPR_name / nbolo; i++) {
            assert( check_pixname( Param->addHPR_name[i*nbolo+ib], pixnames[ib]));
            tperr = noDMC_readObject_PIOFLOAT( Param->addHPR_name[i*nbolo+ib], rg*RINGSIZE, RINGSIZE, addhpr);
            if (tperr<0) {
              fprintf(stderr, "Error %d while reading addHPR_name[%ld]: %s\n", tperr, i*nbolo+ib, Param->addHPR_name[i*nbolo+ib]);
              exit ( -1);
            }
            for (j=0; j<RINGSIZE; j++) {
              addhpr[j] *= Param->addHPR_factor[i*nbolo+ib];
              if (!Param->addHPR_watts[i*nbolo+ib]) {
                // convert addhpr from KCMB to Watts
                addhpr[j] *= Param->Calibration[ib];
              }
              y[0][j] += addhpr[j];
            }
          }
          free( addhpr);
        }



        PIOSTRING tpname;
        for (i=0;i<npixbeam;i++) {
          assert( check_pixname( Param->Theo_noPS[i*nbolo+ib], pixnames[ib]));
          theo[i] = (PIOFLOAT *) malloc(sizeof(PIOFLOAT)*RINGSIZE);
          tperr = noDMC_readObject_PIOFLOAT(Param->Theo_noPS[i*nbolo+ib],rg*RINGSIZE,RINGSIZE,theo[i]);
          if (tperr<0) {
            fprintf(stderr, "Impossible to read Theo_noPS[%ld]: %s %d\n",i*nbolo+ib,Param->Theo_noPS[i*nbolo+ib],tperr);
            exit ( -1);
          }
        }
  
        PIOFLOAT *ADU;
        if (Param->flag_ADU==_PAR_TRUE) {
          assert( check_pixname( Param->ADU[ib], pixnames[ib]));
          ADU = (PIOFLOAT *) malloc(sizeof(PIOFLOAT)*RINGSIZE);
          tperr = noDMC_readObject_PIOFLOAT(Param->ADU[ib],rg*RINGSIZE,RINGSIZE,ADU);
          if (tperr<0) {
            fprintf(stderr, "Impossible to read ADU[%ld]: %s %d\n",ib,Param->ADU[ib],tperr);
            exit ( -1);
          }
        }
        else {
          ADU = (PIOFLOAT *) _PIOMALLOC(sizeof(PIOFLOAT)*RINGSIZE);
          for (i=0;i<RINGSIZE;i++) ADU[i]=0;
        }

        //====================================================================================
        //  E2E 
        //

        if (Param->TESTPOL==0) {
          // inline production of a colored noise TOI and projection in HPR
          assert( check_pixname( Param->Signal_noPS[nbolo+ib], pixnames[ib]));
          long ntoisample=(ENDRINGINDEX(rg)-BEGINRINGINDEX(rg)+1);
          PIOINT *posnoise=(PIOINT *) malloc(sizeof(PIOINT)*ntoisample);
          assert( (posnoise != NULL) && "Not enough memory to allocate HPR index TOI array");
          tperr = noDMC_readObject_PIOINT(Param->Signal_noPS[nbolo+ib],BEGINRINGINDEX(rg),ntoisample,posnoise);

          PIODOUBLE *tpin   = (PIODOUBLE *) in_fft;
          PIODOUBLE *tpout  = (PIODOUBLE *) out_fft;
          // build white noise TOI
          for (j=0;j<nnoisesim;j++) {
            tpin[2*j]=sqrt(-2*log( drand48()))*cos(2*M_PI*drand48());
            tpin[2*j+1]=0;
          }
          fftw_execute(p1);

          // convolve (multiply in frequency domain) white noise TOI with noise Fourier transform
          for (j=0;j<nnoisesim;j++) {
            tpin[2*j]=tpout[2*j]*tfnoise[j];
            tpin[2*j+1]=tpout[2*j+1]*tfnoise[j];
          }
          fftw_execute(p2);

          // projection of colored noise TOI into HPR
          for (int iter = 0; iter < number_of_iterations; iter++) {
            for (i=0;i<ntoisample;i++) {
              if (posnoise[i]>=0&&posnoise[i]<RINGSIZE) {
                if (h[posnoise[i]]>0) {
                  y[iter][posnoise[i]]+=tpout[2*i]/h[posnoise[i]];
                }
              }
            }
          }
          free(posnoise);

          // ADC nonlinearity simulation
          for (int iter = 0; iter < number_of_iterations; iter++) {
            for (i=0;i<RINGSIZE;i++) {
              double adux=ADU[i]*1E-3-Param->Calibration[nbolo+ib*4+1];
              y[iter][i] /= Param->Calibration[nbolo+ib*4]+
                            Param->Calibration[nbolo+ib*4+2]*
                            (adux)*exp(-(adux)*(adux)*Param->Calibration[nbolo+ib*4+3]);
            }
          }

#if 0
          PIOSTRING thepath;
          sprintf(thepath,"/redtruck/delouis/SROLL_OUT/test_%d",rank);
          FILE *fp=fopen(thepath,"w");
          fwrite(y,RINGSIZE*sizeof(PIOFLOAT),1,fp);
          fclose(fp);
          sprintf(thepath,"/redtruck/delouis/SROLL_OUT/testn_%d",rank);
          fp=fopen(thepath,"w");
          fwrite(tpout,nnoisesim*sizeof(PIODOUBLE),1,fp);
          fclose(fp);
          MPI_Barrier(MPI_COMM_WORLD);
          exit(0);
#endif
        }

        //
        //====================================================================================

        if (Param->TESTPOL==4) {
          for (int iter = 0; iter < number_of_iterations; iter++) {
            for (i=0;i<RINGSIZE;i++) {
              double adux=ADU[i]*1E-3-Param->Calibration[nbolo+ib*4+1];
              y[iter][i] *= Param->Calibration[nbolo+ib*4]+
                            Param->Calibration[nbolo+ib*4+2]*
                            (adux)*exp(-(adux)*(adux)*Param->Calibration[nbolo+ib*4+3]);
            }
          }
        }

        assert( check_pixname( Param->Ptg_noPS[ib], pixnames[ib]));
        ph = (PIODOUBLE *) malloc(sizeof(PIODOUBLE)*RINGSIZE);
        tperr = noDMC_readObject_PIODOUBLE(Param->Ptg_noPS[ib],rg*RINGSIZE,RINGSIZE,ph);
        if (tperr<0) {
          fprintf( stderr, "Impossible to read Ptg_noPS[%ld]: %s %d\n", ib, Param->Ptg_noPS[ib], tperr);
          exit ( -1);
        }

        sprintf( tpname,"%s_TUPLE_1",Param->Ptg_noPS[ib]);
        th = (PIODOUBLE *) malloc(sizeof(PIODOUBLE)*RINGSIZE);
        tperr = noDMC_readObject_PIODOUBLE(tpname,rg*RINGSIZE,RINGSIZE,th);
        if (tperr<0) {
          fprintf(stderr, "Impossible to read Ptg_noPS[%ld]: %s %d\n", ib, tpname, tperr);
          exit ( -1);
        }

        sprintf(tpname,"%s_TUPLE_2",Param->Ptg_noPS[ib]);
        psi = (PIODOUBLE *) malloc(sizeof(PIODOUBLE)*RINGSIZE);
        tperr = noDMC_readObject_PIODOUBLE(tpname,rg*RINGSIZE,RINGSIZE,psi);
        if (tperr<0) {
          fprintf(stderr, "Impossible to read Ptg_noPS[%ld]: %s %d\n", ib, tpname, tperr);
          exit ( -1);
        }

        phase=(PIODOUBLE *) malloc(sizeof(PIODOUBLE)*RINGSIZE);
        memset(phase,0,sizeof(PIODOUBLE)*RINGSIZE);
        double mmat[9],rres[3];
        double v0[3],v1[3],v2[3];
        int k,l;

        memset(mmat,0,9*sizeof(double));
        memset(rres,0,3*sizeof(double));


        for (i=0;i<RINGSIZE;i++) if (h[i]>0) {
          double vec[3];
          ang2vec(th[i],ph[i],vec);
          for (k=0;k<3;k++) for (l=0;l<3;l++) mmat[k+3*l]+=vec[k]*vec[l];
          for (k=0;k<3;k++) rres[k]+=vec[k];
        }

        if (lusol(mmat,rres,3)==0) {
          for (k=0;k<3;k++) v2[k]=rres[k];
          double tmp=0;
          for (k=0;k<3;k++) tmp+=v2[k]*v2[k];
          for (k=0;k<3;k++) v2[k]/=sqrt(tmp);
          v1[0]=0;v1[1]=1;v1[2]=-v2[1]/v2[2];
          tmp=0;
          for (k=0;k<3;k++) tmp+=v1[k]*v1[k];
          for (k=0;k<3;k++) v1[k]/=sqrt(tmp);
          mmat[0]=v1[0];mmat[1]=v1[1];mmat[2]=v2[0];mmat[3]=v2[1];
          rres[0]=v1[2];rres[1]=v2[2];
          if (lusol(mmat,rres,2)!=0) {
            fprintf(stderr,"%s %d %d\n",__FILE__,__LINE__,rank);
          }
          v0[0]=rres[0];v0[1]=rres[1];v0[2]=-1;
          tmp=0;
          for (k=0;k<3;k++) tmp+=v0[k]*v0[k];
          for (k=0;k<3;k++) v0[k]/=sqrt(tmp);
          double navr=0,avr=0,avr2=0;
          //fprintf(stderr,"PHASE v0 %ld %lf %lf %lf\n",(long) rg,v0[0],v0[1],v0[2]);
          //fprintf(stderr,"PHASE v1 %ld %lf %lf %lf\n",(long) rg,v1[0],v1[1],v1[2]);
          //fprintf(stderr,"PHASE v2 %ld %lf %lf %lf\n",(long) rg,v2[0],v2[1],v2[2]);

          for (i=0;i<RINGSIZE;i++) if (h[i]>0) {
            double vec[3];
            ang2vec(th[i],ph[i],vec);
            double xx=vec[0]*v0[0]+vec[1]*v0[1]+vec[2]*v0[2];
            double yy=vec[0]*v1[0]+vec[1]*v1[1]+vec[2]*v1[2];
            phase[i]=atan2(yy,xx);
            tmp=sqrt(xx*xx+yy*yy);
            navr+=1;
            avr+=tmp;
            avr2+=tmp*tmp;
          }

          long nrg_htmp = 0;
          for (i=0;i<RINGSIZE;i++) phase[i]=fmod(phase[i]+M_PI+rg*M_PI/30.,2*M_PI);

          //long off_adu=(32000/Param->REMDIP)*((rg-globalBeginRing)/stepadu);
          //if (off_adu>=32000) off_adu=(32000/Param->REMDIP)*(Param->REMDIP-1);

          for (i=0;i<RINGSIZE;i++) if (h[i]>0) {
            long ipix;
            long ipix128;
            nrg_htmp++;
            double vecpix[3];
            ang2pix_ring( Nside, th[i], ph[i], &ipix);
            ang2pix_ring( 128, th[i], ph[i], &ipix128);
            int testmap=0;
            if (dustmapI!=NULL) {
              if (dustmapI[ipix128]>-1E10) testmap=1;
              else testmap=0;
            }
            else testmap=1;

            if (testmap) {
              ang2vec(th[i], ph[i],vecpix);
              hpix *tp_hpix;

              if (l_nhpix[ipix]==0) {
                l_hpix[ipix] = (hpix *) malloc(sizeof(hpix));
              }
              else {
                l_hpix[ipix] = (hpix *) realloc(l_hpix[ipix],(l_nhpix[ipix]+1)*sizeof(hpix));
              }
              tp_hpix=l_hpix[ipix]+l_nhpix[ipix];
              for (int iter = 0; iter < number_of_iterations; iter++) {
                tp_hpix->listp[iter]=(y[iter][i]-Param->Monop[ib])/Param->Calibration[ib];
              }
              tp_hpix->sig=tp_hpix->listp[0];

              //double soldip=SOLDIPX*vecpix[0]+SOLDIPY*vecpix[1]+SOLDIPZ*vecpix[2];
              //tp_hpix->soldip = diporb[i];
              tp_hpix->dip    = diporb[i];
              //tp_hpix->thsig  = diporb[i];
              tp_hpix->hit    = 1/sqrt(h[i]);
              if (fsl != NULL) {
                tp_hpix->fsl = fsl[i]*Param->FSLCOEF[ib];
              } else {
                tp_hpix->fsl = 0.0;
              }



              tp_hpix->adu    = (PIOINT) (floor(ADU[i]));
              if (tp_hpix->adu+16000<0) tp_hpix->adu=-16000;
              if (tp_hpix->adu+16000>31999) tp_hpix->adu=15999;
              tp_hpix->adu=(rg-240)/100+ADURGSTEP*(tp_hpix->adu/10+1600);

#if 0
              tp_hpix->adu=(tp_hpix->adu+16000)/Param->REMDIP+off_adu;
              if (tp_hpix->adu<=0||tp_hpix->adu>31999) {
                fprintf(stderr,"ADUPROB %lg %lg %lg\n",ADU[i],(double) off_adu,diporb[i]);
              }
#endif


              tp_hpix->co= cos(2*psi[i]);
              tp_hpix->si= sin(2*psi[i]);
              for (j=0;j<npixbeam;j++) tp_hpix->listofpix[j]=theo[j][i]; //-diporb[i];

              tp_hpix->w=h[i]*sxi;
              if (comapI!=NULL) tp_hpix->comap=(comapI[ipix128]);
              else tp_hpix->comap=0;
              if (FREEFREE!=NULL) tp_hpix->freefree=(FREEFREE[ipix128]);
              else tp_hpix->freefree=0;

              if (dustmapI!=NULL) tp_hpix->dustmap=(dustmapI[ipix128]
                                                    +eta[ib]*tp_hpix->co*dustmapQ[ipix128]
                                                    +eta[ib]*tp_hpix->si*dustmapU[ipix128]);
              else tp_hpix->dustmap=0;

#if 1
              if (dustmapI!=NULL&&Param->OUT_NOPOL[0]%2==1) {
                if (nbolo==1) {
                  tp_hpix->dustmap =tp_hpix->co*comapI[ipix128]+tp_hpix->si*dustmapI[ipix128];
                  tp_hpix->comap   =-tp_hpix->si*comapI[ipix128]+tp_hpix->co*dustmapI[ipix128];
                }
                if (dustmapI!=NULL&&(Param->OUT_NOPOL[0]/2)%2==1) {
                  tp_hpix->listofpix[npixbeam]=tp_hpix->co*dustmapQ[ipix128]+tp_hpix->si*dustmapU[ipix128];
                }
              }

#endif

              if (simdustmapQ!=NULL) tp_hpix->sig+=(eta[ib]*tp_hpix->co*simdustmapQ[ipix]
                                                    +eta[ib]*tp_hpix->si*simdustmapU[ipix])*0.0007;
              tp_hpix->surv = surv;
              tp_hpix->ib = ib;
              tp_hpix->rg = rg;
#if GAIN_RATIO
              if (gain_ratio_off[ib][rg]>10000000) {
                tp_hpix->gi = gain_ratio_off[ib][rg];
                //tp_hpix->ggi = gain_ratio[ib][rg];
              }
              else {
              }
#endif
#define NBPH (4096)
              tp_hpix->gi = rg;
              //tp_hpix->ggi = 1.;
              tp_hpix->hrg = (phase[i]*NBPH)/(2*M_PI);
              //tp_hpix->xrg = (double) (rg-240)/25000.;
              tp_hpix->ipix = ipix;
              l_nhpix[ipix]+=1;
              if (isnan(tp_hpix->sig)||isnan(tp_hpix->w)||isnan(tp_hpix->dip)) {
                fprintf(stderr,"NAN NAN NAN %lf %lf %lf %ld %ld\n",tp_hpix->dip,tp_hpix->sig,tp_hpix->w,(long)rg,(long)i);
              }
              if (Param->TESTPOL==3) {
//                double xxx=((tp_hpix->adu-off_adu)-16000/Param->REMDIP)*Param->REMDIP/1000.;
//                tp_hpix->sig  = (0.99+0.001*tp_hpix->ib)*(tp_hpix->dip+tp_hpix->fsl) + 1E-5*exp(-xxx*xxx);
              }

            }
          }
          if ((rank==0) && (rg==globalRankInfo.BeginRing[rank]) ) {
            GetProcMem(&vmem,&phymem);
            if (Param->TESTPOL==0) {
              fprintf(stderr,"Com %s Rank: %ld[%d] MEM %.1lf[%.1lf]MB Nd=%ld %lg %lg [%lg,%lg,%lg,%lg]\n",
                  Command, (long) rank, getpid(),
                  (double) vmem/1024./1024.,
                  (double) phymem/1024./1024.,
                  (long) nrg_htmp,
                  avr/navr,sqrt(avr2/navr-(avr/navr)*(avr/navr)),
                               Param->Calibration[nbolo+ib*4],Param->Calibration[nbolo+ib*4+1],
                               Param->Calibration[nbolo+ib*4+2],Param->Calibration[nbolo+ib*4+3]);
            }
            else {
              fprintf(stderr,"Com %s Rank: %ld[%d] MEM %.1lf[%.1lf]MB Nd=%ld %lg %lg\n",
                  Command, (long) rank, getpid(),
                  (double) vmem/1024./1024.,
                  (double) phymem/1024./1024.,
                  (long) nrg_htmp,
                  avr/navr,sqrt(avr2/navr-(avr/navr)*(avr/navr)));
            }
          }
        }
  
        free(h);
        if (Param->flag_stim_paramfiles == 0) {
          free(y[0]);
        }
        free(fsl);
        free(diporb); 
        for (i=0;i<npixbeam;i++) free(theo[i]);
        free(phase);
        free(th);
        free(ph);
        free(psi);
        free(ADU);
      } // if not badring
    } // end ring loop
    
    if (Param->TESTPOL==0) {
      free(tfnoise);
    }

    if (Param->flag_stim_paramfiles == 1) {
      for (int iter = 0; iter < number_of_iterations; iter++) {
        free( stim_hpr[iter]);
      }
    }

    if (rank==0) {
      GetProcMem(&vmem,&phymem);
      fprintf(stderr,"\nafter sroll bolometer: used VMEM %.1lf[PHYS %.1lf]MB\n",
          (double) vmem/1024./1024., (double) phymem/1024./1024.);
    }

  } // end bolometer loop
  
  fftw_free(in_fft);
  fftw_free(out_fft);

  if (comapI!=NULL) {
    free(comapI);
  }
  if (FREEFREE!=NULL) {
    free(FREEFREE);
  }
  if (dustmapI!=NULL) {
    free(dustmapI);
    free(dustmapQ);
    free(dustmapU);
  }

  PrintFreeMemOnNodes( rank, mpi_size, "before pixel balancing");

  if ((Param->OUT_NOPOL[0]/2)%2==1) npixbeam++;
  /*======================================================
    =
    =      compute load balancing between proc
    =
    =*/
  double *stat_pix = (double *) malloc( 12*Nside*Nside*sizeof(double));
  memset(stat_pix,0,12*Nside*Nside*sizeof(double));

  long *begpix = (long *) malloc(sizeof(long)*mpi_size);
  long *edpix = (long *) malloc(sizeof(long)*mpi_size);
  long maxitt = 0;
  long rrk;

  //int *inv_stat_pix;
  {
    long l_nr=12*Nside*Nside;
    for (i=0;i<l_nr;i++) {
      stat_pix[i]=l_nhpix[i];
    }
   

    if (rank==0) {
      MPI_Status statu;
      double *l_stat_pix = (double *) malloc(12*Nside*Nside*sizeof(double));
      for (rrk=1;rrk<mpi_size;rrk++) {
        MPI_Recv(l_stat_pix,sizeof(double)*12*Nside*Nside, MPI_BYTE, rrk,31, MPI_COMM_WORLD,&statu);
        for (j=0;j<12*Nside*Nside;j++) stat_pix[j]+=l_stat_pix[j];
      }
      free(l_stat_pix);
      nmatpix=0;
      for (j=0;j<12*Nside*Nside;j++) {
        if (stat_pix[j]>0) {
          nmatpix++;
        }
      }
      stat_pix[0]=stat_pix[0];
      maxitt=stat_pix[0];
      for (j=1;j<12*Nside*Nside;j++) {
        if (maxitt<stat_pix[j]) maxitt=stat_pix[j];
        stat_pix[j]=stat_pix[j]+stat_pix[j-1];
        
      }
      double step_pix=stat_pix[12*Nside*Nside-1]/mpi_size;
      fprintf(stderr,"step_pix=%lf\n",step_pix);
      int k=1;
      begpix[0]=0;
      for (j=1;j<12*Nside*Nside;j++) {
        if (k*step_pix<stat_pix[j]) {
          edpix[k-1]=j-1;
          begpix[k]=j;
          k++;
        }
      }
      edpix[mpi_size-1]=12*Nside*Nside-1;
      for (j=0;j<mpi_size;j++) {
        fprintf(stderr,"rank=%ld begpix=%ld endpix=%ld npix=%ld %lg\n",j,begpix[j],edpix[j],edpix[j]-begpix[j],stat_pix[edpix[j]]-stat_pix[begpix[j]]);
      }
    }
    else  {
      MPI_Send(stat_pix, sizeof(double)*12*Nside*Nside, MPI_BYTE, 0, 31, MPI_COMM_WORLD);
    }
    MPI_Bcast(begpix,sizeof(long)*mpi_size, MPI_BYTE, 0, MPI_COMM_WORLD);
    MPI_Bcast(edpix,sizeof(long)*mpi_size, MPI_BYTE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&maxitt,sizeof(long), MPI_BYTE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nmatpix,sizeof(long),MPI_BYTE, 0, MPI_COMM_WORLD);
  }
  free(stat_pix);

  GetProcMem(&vmem,&phymem);
  if (rank==0) fprintf(stderr,"Rank: %ld[%d] MEM %.1lf[%.1lf]MB line=%d\n",
              (long) rank, getpid(),
              (double) vmem/1024./1024.,
              (double) phymem/1024./1024.,__LINE__);
  
  /*======================================================
    =
    =      order data per pixel and per proc
    =
    =*/
  
  nnbpix = edpix[rank]-begpix[rank]+1;

  loc_hpix = (hpix **) malloc(sizeof(hpix *)*nnbpix);
  loc_nhpix = (PIOLONG *) malloc(sizeof(PIOLONG)*nnbpix);
  memset(loc_nhpix,0,sizeof(PIOLONG)*nnbpix);

  MPI_Barrier(MPI_COMM_WORLD);

#if 0
  long ndata=0;
  for (rrk=0;rrk<mpi_size;rrk++) {
    MPI_Status statu;
    PIOLONG l_nsamp;
    if (rrk==rank) {
      int rk2;
      memcpy(loc_nhpix,l_nhpix+begpix[rank],sizeof(long)*(nnbpix));
      for (j=0;j<nnbpix;j++) loc_hpix[j]=l_hpix[j+begpix[rank]];
      for (rk2=0;rk2<mpi_size;rk2++) {
        hpix *tbs;
        long ntbs=0,k;
        if (rk2!=rrk) {
          MPI_Recv(&ntbs,sizeof(long), MPI_BYTE, rk2,450, MPI_COMM_WORLD,&statu);

          GetProcMem(&vmem,&phymem);
          if (rank==0) fprintf(stderr,"Rank: %ld[%d] MEM %.1lf[%.1lf]MB line=%d %ld\n",
                  (long) rank, rk2,
                  (double) vmem/1024./1024.,
                  (double) phymem/1024./1024.,__LINE__,ntbs);
          if (ntbs>0) {
            tbs=(hpix *) malloc(sizeof(hpix)*(ntbs));
            MPI_Recv(tbs,sizeof(hpix)*ntbs, MPI_BYTE, rk2,451, MPI_COMM_WORLD,&statu);
            
            for (k=0;k<ntbs;k++) {
              j=tbs[k].ipix-begpix[rrk];
              if (loc_nhpix[j]==0) {
                loc_hpix[j] = (hpix *) malloc(sizeof(hpix));
              }
              else {
                loc_hpix[j] = (hpix *) realloc(loc_hpix[j],(1+loc_nhpix[j])*sizeof(hpix));
              }
              memcpy(loc_hpix[j]+loc_nhpix[j],tbs+k,sizeof(hpix));
              loc_nhpix[j] +=1;
            }

            free(tbs);
          }
        }
      }
    }
    else {
      hpix *tbs;
      long ntbs=0;
      for (j=begpix[rrk];j<=edpix[rrk];j++) if (l_nhpix[j]>0) {
        if (ntbs==0) {
          tbs=(hpix *) malloc(sizeof(hpix)*(l_nhpix[j]));
        }
        else {
          tbs=(hpix *) realloc(tbs,sizeof(hpix)*(ntbs+l_nhpix[j]));
        }
          
        memcpy(tbs+ntbs,l_hpix[j], sizeof(hpix)*(l_nhpix[j]));
        free(l_hpix[j]);
        ntbs+=l_nhpix[j];
      }
      MPI_Send(&ntbs, sizeof(long), MPI_BYTE, rrk, 450, MPI_COMM_WORLD);
      if (ntbs>0) {
        MPI_Send(tbs, sizeof(hpix)*(ntbs), MPI_BYTE, rrk, 451, MPI_COMM_WORLD);
        free(tbs);
      }
    }
  }
#else

  //========================================================================
  // USE iterative approach to avoid slow down while exchanging data
  long nstep=log((double)(mpi_size))/log(2.0)+0.001;
  long exch=1;
  int is,js;
  for (is=0;is<nstep;is++) {
    MPI_Status statu;
    exch*=2;
    for (js=0;js<mpi_size;js+=exch) {
      long rk0=exch*((long)(rank/exch));
      int k;
      for (k=0;k<exch;k++) {
        if (rank%exch==k) {
          hpix *tbs;
          long ntbs=0,kk;
          //fprintf (stderr,"%d <- %d\n",rank,rk0+(k+exch/2)%(exch));
          MPI_Recv(&ntbs,sizeof(long), MPI_BYTE,rk0+(k+exch/2)%(exch),450, MPI_COMM_WORLD,&statu);
          
          if (rank==0) {
            GetProcMem( &vmem, &phymem);
            fprintf(stderr,"Rank: %d(%s, free=%.2fGB) [%d/%d] used MEM virt:%.1lf[phys:%.1lf]MB line=%d ntbs:%ldMB\n",
                  rank, hostname, GetFreeMemGB(), (int) (rk0+(k+exch/2)%(exch)),(int) exch,
                  (double) vmem/1024./1024.,
                  (double) phymem/1024./1024.,__LINE__,sizeof(hpix)*ntbs/1024/1024);
          }
          if (ntbs>0) {
            tbs=(hpix *) malloc(sizeof(hpix)*(ntbs));
            if (tbs == NULL) {
              GetProcMem( &vmem, &phymem);
              fprintf(stderr,"Rank: %d(%s, free=%.2fGB) [%d/%d] used MEM virt:%.1lf[phys:%.1lf]MB line=%d ntbs:%ldMB\n",
                    rank, hostname, GetFreeMemGB(), (int) (rk0+(k+exch/2)%(exch)),(int) exch,
                    (double) vmem/1024./1024.,
                    (double) phymem/1024./1024.,__LINE__,sizeof(hpix)*ntbs/1024/1024);
              sleep(10); // wait for other ranks to fail allocation and display error message
              return 1;
            }
            MPI_Recv(tbs,sizeof(hpix)*ntbs, MPI_BYTE, rk0+(k+exch/2)%(exch),451, MPI_COMM_WORLD,&statu);
            
            for (kk=0;kk<ntbs;kk++) {
              j=tbs[kk].ipix;
              if (l_nhpix[j]==0) {
                l_hpix[j] = (hpix *) malloc(sizeof(hpix));
              }
              else {
                l_hpix[j] = (hpix *) realloc(l_hpix[j],(1+l_nhpix[j])*sizeof(hpix));
              }
              memcpy(l_hpix[j]+l_nhpix[j],tbs+kk,sizeof(hpix));
              l_nhpix[j] +=1;
            }
            
            free(tbs);
          }
        }
        if ((rank+exch/2)%exch==k) {
          hpix *tbs=NULL;
          long ntbs=0;
          for (j=begpix[js+k];j<=edpix[js+k];j++) if (l_nhpix[j]>0) {
            if (ntbs==0) {
              tbs=(hpix *) malloc(sizeof(hpix)*(l_nhpix[j]));
              if (tbs == NULL) {
                GetProcMem( &vmem, &phymem);
                fprintf(stderr,"Rank: %d(%s, free=%.2fGB) [%d/%d] used MEM virt:%.1lf[phys:%.1lf]MB line=%d ntbs:%ldMB\n",
                      rank, hostname, GetFreeMemGB(), (int) (rk0+(k+exch/2)%(exch)),(int) exch,
                      (double) vmem/1024./1024.,
                      (double) phymem/1024./1024.,__LINE__,sizeof(hpix)*ntbs/1024/1024);
                sleep(10); // wait for other ranks to fail allocation and display error message
                return 1;
              }
            }
            else {
              tbs=(hpix *) realloc(tbs,sizeof(hpix)*(ntbs+l_nhpix[j]));
              if (tbs == NULL) {
                GetProcMem( &vmem, &phymem);
                fprintf(stderr,"Rank: %d(%s, free=%.2fGB) [%d/%d] used MEM virt:%.1lf[phys:%.1lf]MB line=%d ntbs:%ldMB\n",
                      rank, hostname, GetFreeMemGB(), (int) (rk0+(k+exch/2)%(exch)),(int) exch,
                      (double) vmem/1024./1024.,
                      (double) phymem/1024./1024.,__LINE__,sizeof(hpix)*ntbs/1024/1024);
                sleep(10); // wait for other ranks to fail allocation and display error message
                return 1;
              }
            }
            
            memcpy(tbs+ntbs,l_hpix[j], sizeof(hpix)*(l_nhpix[j]));
            free(l_hpix[j]);
            ntbs+=l_nhpix[j];
          }
          //fprintf (stderr,"%d -> %d\n",rank,rk0+k);
          MPI_Send(&ntbs, sizeof(long), MPI_BYTE, rk0+k, 450, MPI_COMM_WORLD);
          if (ntbs>0) {
            MPI_Send(tbs, sizeof(hpix)*(ntbs), MPI_BYTE, rk0+k, 451, MPI_COMM_WORLD);
            free(tbs);
          }
        }
      }
    }
  }

  memcpy(loc_nhpix,l_nhpix+begpix[rank],sizeof(long)*(nnbpix));
  for (j=0;j<nnbpix;j++) loc_hpix[j]=l_hpix[j+begpix[rank]];
#endif

  free(l_hpix);
  free(l_nhpix);


  PIOINT *mask;
  if (Param->flag_Mask==_PAR_TRUE) {
    PIOSTRING commask;
    sprintf(commask,"begin=%lld;end=%lld",(long long) begpix[rank],(long long) edpix[rank]);
    mask = (PIOINT *) malloc(sizeof(PIOINT)*(edpix[rank]-begpix[rank]+1));
    long resmask=(long) noDMC_readObject_PIOINT(Param->Mask,begpix[rank],edpix[rank]-begpix[rank]+1,mask);
    if (rank==0) fprintf(stderr,"Mask %ld\n",(long) resmask);
  }
  else {
    mask = (PIOINT *) _PIOMALLOC(sizeof(PIOINT)*nnbpix);
    memset(mask,1,sizeof(PIOINT)*nnbpix);
  }
  
#ifndef DOMAP
  PIOFLOAT **alm_map = (PIOFLOAT **) malloc(LMAX*sizeof(PIOFLOAT *));
  {
    PIOSTRING commask;
    PIOFLOAT *tmpalm;
    sprintf(commask,"begin=%lld;end=%lld",(long long) begpix[rank],(long long) edpix[rank]);
    for (i=0;i<LMAX;i++) {
      fprintf(stderr,"ALM %ld\n",(long) PIOReadMAPObject((void **) &tmpalm,Param->ALMMAP[i],"PIOFLOAT",commask,NULL));
      alm_map[i]=tmpalm;
    }
  }
#endif

  GetProcMem(&vmem,&phymem);
  if (rank%64==0) fprintf(stderr,"Rank: %ld[%d] MEM %.1lf[%.1lf]MB line=%d\n",
              (long) rank, getpid(),
              (double) vmem/1024./1024.,
              (double) phymem/1024./1024.,__LINE__);
  
  PrintFreeMemOnNodes( rank, mpi_size, "before matrix computation");

  /*======================================================
    =
    =     Compute matrix size
    =
    =*/
  flg_rg = (PIOBYTE **) malloc(nbolo*sizeof(PIOBYTE *));
  for (ib=0;ib<nbolo;ib++) {
    flg_rg[ib]= (PIOBYTE *) malloc(sizeof(PIOBYTE)*
                                       (globalRangeRing));

    memset(flg_rg[ib],0,sizeof(PIOBYTE)*(globalRangeRing));
  }
  //PIOBYTE **flg_rg2 = (PIOBYTE **) malloc(nbolo*sizeof(PIOBYTE *));
  //for (ib=0;ib<nbolo;ib++) {
  //  flg_rg2[ib]= (PIOBYTE *) malloc(sizeof(PIOBYTE)*
  //                                   (globalRangeRing)*CUTRG);
  //
  //  memset(flg_rg2[ib],0,sizeof(PIOBYTE)*(globalRangeRing)*CUTRG);
  //}

  PIOLONG k;
  long i0;
  flgpix = (PIOBYTE *) malloc(sizeof(PIOBYTE)*nnbpix);
  memset(flgpix,0,sizeof(PIOBYTE)*nnbpix);

  MPI_Status statu;

  long l_nmatpix=0;
  for (k=0;k<nnbpix;k++) if (mask[k]==1) {
    long ndata = loc_nhpix[k];
    if (ndata>2) {
      hpix *htmp = loc_hpix[k];
      double II=0,QQ=0,UU=0,QU=0,IQ=0,IU=0;
      for (i0=0;i0<ndata;i0++) {
        double CO0=eta[(int) htmp[i0].ib]*(dpsico[(int) htmp[i0].ib]*htmp[i0].co
                                     -dpsisi[(int) htmp[i0].ib]*htmp[i0].si);
        double SI0=eta[(int) htmp[i0].ib]*(dpsico[(int) htmp[i0].ib]*htmp[i0].si
                                     +dpsisi[(int) htmp[i0].ib]*htmp[i0].co);
        II+= htmp[i0].w;
        QQ+= htmp[i0].w*CO0*CO0;
        UU+= htmp[i0].w*SI0*SI0;
        IQ+= htmp[i0].w*CO0;
        IU+= htmp[i0].w*SI0;
        QU+= htmp[i0].w*SI0*CO0;

      }
      if (Param->OUT_NOPOL[0]%2==1) {
        flgpix[k]=1;
        for (i0=0;i0<ndata;i0++) {
          flg_rg[htmp[i0].ib][htmp[i0].rg-globalBeginRing]=1;
          l_nmatpix++;
        }
      }
      else {
        double cond=cond_3_3_thres(II,IQ,IU,IQ,QQ,QU,IU,QU,UU);
        
        if (cond < Param->seuilcond) {
          flgpix[k]=1;
          for (i0=0;i0<ndata;i0++) {
            flg_rg[htmp[i0].ib][htmp[i0].rg-globalBeginRing]=1;
            //flg_rg2[htmp[i0].ib][htmp[i0].hrg-globalBeginRing*CUTRG]=1;
            l_nmatpix++;
          }
        }
        else {
          flgpix[k]=0;

        // Print detail for pix that does not meet cond requirement
//        fprintf(stderr,"[DBG COND] Pix#%ld is flagged out! II=%g IQ=%g IU=%g QQ=%g QU=%g UU=%g \n",
//                k+begpix[rank], II, IQ, IU, QQ, QU, UU);
        }
      }
    }
    else {
      //      fprintf(stderr,"[DBG COND] Pix#%ld is flagged out! ndata<=2\n", k+begpix[rank]);
      flgpix[k]=0;
    }
  }


  //inv_stat_pix = (int *) malloc(l_nmatpix*sizeof(int));
  //the_stat_pix = (int *) malloc(sizeof(int)*nnbpix);

  //for (j=0;j<nnbpix;j++) the_stat_pix[j]=-1;

  l_nmatpix=0;
  for (j=0;j<nnbpix;j++) {
    if (flgpix[j]==1) {
      //inv_stat_pix[l_nmatpix]=j+begpix[rank];
      //the_stat_pix[j]=l_nmatpix;
      l_nmatpix++;
    }
  }
  if (rank==0) fprintf(stderr,"l_nmatpix[%d] %ld\n",rank,(long) l_nmatpix);

  if (rank==0) {
    long lb;
    nmatpix=l_nmatpix;
    for (rrk=1;rrk<mpi_size;rrk++) {
      MPI_Recv(&lb,sizeof(long), MPI_BYTE, rrk,450, MPI_COMM_WORLD,&statu);
      nmatpix+=lb;
    }
  }
  else MPI_Send(&l_nmatpix, sizeof(long), MPI_BYTE, 0, 450, MPI_COMM_WORLD);
  MPI_Bcast(&nmatpix, sizeof(long), MPI_BYTE, 0, MPI_COMM_WORLD);

  for (ib=0;ib<nbolo;ib++) {
    PIOBYTE *l_flg_rg = (PIOBYTE *) malloc(sizeof(PIOBYTE)*
                                           (globalRangeRing));

    for (rrk=0;rrk<mpi_size;rrk++) {
      if (rrk==rank) memcpy(l_flg_rg ,flg_rg[ib] ,sizeof(PIOBYTE)*
                            (globalRangeRing));
      MPI_Bcast(l_flg_rg,sizeof(PIOBYTE)*
                (globalRangeRing), MPI_BYTE, rrk, MPI_COMM_WORLD);
      for (i=0;i<(globalRangeRing);i++) if (l_flg_rg[i]==1) flg_rg[ib][i]=1;

    }
    free(l_flg_rg);
  }

  
  //for (ib=0;ib<nbolo;ib++) {
  //  PIOBYTE *l_flg_rg = (PIOBYTE *) malloc(sizeof(PIOBYTE)*
  //                                       (globalRangeRing)*CUTRG);
  //
  //  for (rrk=0;rrk<mpi_size;rrk++) {
  //    if (rrk==rank) memcpy(l_flg_rg ,flg_rg2[ib] ,sizeof(PIOBYTE)*
  //                        (globalRangeRing)*CUTRG);
  //    MPI_Bcast(l_flg_rg,sizeof(PIOBYTE)*
  //            (globalRangeRing)*CUTRG, MPI_BYTE, rrk, MPI_COMM_WORLD);
  //    for (i=0;i<(globalRangeRing)*CUTRG;i++) if (l_flg_rg[i]==1) flg_rg2[ib][i]=1;
  //
  //  }
  //  free(l_flg_rg);
  //}



  GetProcMem(&vmem,&phymem);
  if (rank%64==0) fprintf(stderr,"Rank: %ld[%d] MEM %.1lf[%.1lf]MB line=%d\n",
              (long) rank, getpid(),
              (double) vmem/1024./1024.,
              (double) phymem/1024./1024.,__LINE__);

  rgord = (PIOINT **) malloc(nbolo*sizeof(PIOINT *));
  PIOINT **rgordinv = (PIOINT **) malloc(nbolo*sizeof(PIOINT *));
  newnr  = (PIOLONG *) malloc((nbolo+1)*sizeof(PIOLONG));
  rgord2 = (PIOINT **) malloc(nbolo*sizeof(PIOINT *));
  PIOINT **rgordinv2 = (PIOINT **) malloc(nbolo*sizeof(PIOINT *));
  newnr2  = (PIOLONG *) malloc((nbolo+1)*sizeof(PIOLONG));

  newnr[0]=0;
  for (ib=0;ib<nbolo;ib++) {
    rgord[ib] = (PIOINT *) malloc(sizeof(PIOINT)*
                                  (globalRangeRing));
    rgordinv[ib] = (PIOINT *) malloc(sizeof(PIOINT)*
                                     (globalRangeRing));

    newnr[ib+1]=0;
    for (i=0;i<(globalRangeRing);i++) if (flg_rg[ib][i]>0) {
      rgord[ib][i]=newnr[ib+1];
      rgordinv[ib][newnr[ib+1]]=i;
      newnr[ib+1]++;
    }
  }
  for (ib=1;ib<nbolo+1;ib++) newnr[ib]+=newnr[ib-1];

  newnr2[0]=0;
  for (ib=0;ib<nbolo;ib++) {
    rgord2[ib] = (PIOINT *) malloc(sizeof(PIOINT)*
                                  (globalRangeRing)*CUTRG);
    rgordinv2[ib] = (PIOINT *) malloc(sizeof(PIOINT)*
                                     (globalRangeRing)*CUTRG);

    newnr2[ib+1]=0;
    for (i=0;i<(globalRangeRing)*CUTRG;i++) if (flg_rg[ib][i/CUTRG]>0) {
      rgord2[ib][i]=newnr2[ib+1];
      rgordinv2[ib][newnr2[ib+1]]=i;
      newnr2[ib+1]++;
    }
  }
  for (ib=1;ib<nbolo+1;ib++) newnr2[ib]+=newnr2[ib-1];

  nmatco=nbolo;
  if (Param->flag_Theo_CO!=_PAR_TRUE) nmatco=0;
  nfreefree=nbolo;
  if (Param->flag_Theo_FREEFREE!=_PAR_TRUE) nfreefree=0;
  nmatdust=nbolo;
  if (Param->flag_Theo_Dust_I!=_PAR_TRUE) nmatdust=0;

  if (rank==0) fprintf(stderr,"NMATCO %ld NMATDUST %ld\n FREEFREE %ld\n",(long)  nmatco,(long)  nmatdust,(long)  nfreefree);

#ifdef DONSIDE 
  PIOLONG nmat= nmatpix*2+newnr[nbolo]+nbolo*npixbeam+nbolo*(DEGREE+nadu3)+nmatco+nmatdust+nfreefree;
#else
  PIOLONG nmat= nmatpix*2+newnr[nbolo]+nbolo*(DEGREE+nadu3)+nmatco+nmatdust+nfreefree;
#endif


  if (rank==0) {
    fprintf(stderr,"RK%d SHOULD DETERMINE %ld VALUES ",rank,(long) nmat);
    for (i=0;i<nbolo;i++) fprintf(stderr,"%ld ",(long) newnr[i]);
    fprintf(stderr,"\n");
  }

  if (Param->TESTPOL==7) {
    double dpsico2[30];
    double dpsisi2[30];
    for (k=0;k<nbolo;k++) dpsico2[k]=cos((k-(nbolo-1)/2.)*0.0174);
    for (k=0;k<nbolo;k++) dpsisi2[k]=sin((k-(nbolo-1)/2.)*0.0174);
    for (k=0;k<nnbpix;k++)  {
      long ndata = loc_nhpix[k];
      hpix *htmp = loc_hpix[k];
      int l1;
      double l_th,l_ph;
      pix2ang_ring(Nside,k+begpix[rank],&l_th,&l_ph);
      for (l1=0;l1<ndata;l1++) {
            
        double CO1=eta[htmp[l1].ib]*(dpsico2[htmp[l1].ib]*htmp[l1].co
                                        -dpsisi2[htmp[l1].ib]*htmp[l1].si);
        double SI1=eta[htmp[l1].ib]*(dpsico2[htmp[l1].ib]*htmp[l1].si
                                    +dpsisi2[htmp[l1].ib]*htmp[l1].co);
        
          htmp[l1].listp[0] = htmp[l1].dip+htmp[l1].fsl+1E-4*(CO1*cos(l_th*4)+SI1*sin(l_th*3));
      }
    }
  }

  if (Param->REMHDIP==1) {
    PIODOUBLE *templatemap;
    templatemap = (PIODOUBLE *) malloc(sizeof(PIODOUBLE)*(edpix[rank]-begpix[rank]+1));
    PIOLONG nmask = noDMC_readObject_PIODOUBLE(Param->TEMPLATEMAP,begpix[rank],edpix[rank]-begpix[rank]+1,templatemap);
    if (nmask!=0) {
      fprintf(stderr,"Impossible to read TEMPLATEMAP %ld %ld\n",(long) (edpix[rank]-begpix[rank]+1), (long) nmask);
      exit(0);
    }
    
    for (k=0;k<nnbpix;k++)  {
      long ndata = loc_nhpix[k];
      hpix *htmp = loc_hpix[k];
      int l1;
      double tmpmap=templatemap[k];
      for (l1=0;l1<ndata;l1++) {
        
          htmp[l1].freefree  = htmp[l1].dip;
          htmp[l1].dip  = tmpmap+htmp[l1].dip;
      }
    }
    free(templatemap);
  }

#if 0 //NOT YET USED
  else {
    PIODOUBLE *templatemap;
    templatemap = (PIODOUBLE *) malloc(sizeof(PIODOUBLE)*(edpix[rank]-begpix[rank]+1));
    PIOLONG nmask = noDMC_readObject_PIODOUBLE(Param->TEMPLATEMAP,begpix[rank],edpix[rank]-begpix[rank]+1,templatemap);
    if (nmask!=0) {
      fprintf(stderr,"Impossible to read templatemap\n");
      exit(0);
    }
    for (k=0;k<nnbpix;k++)  {
      long ndata = loc_nhpix[k];
      hpix *htmp = loc_hpix[k];
      int l1;
      double tmpmap=templatemap[k];
      for (l1=0;l1<ndata;l1++) {
        
          //htmp[l1].dipmod  = tmpmap+htmp[l1].dip;
          htmp[l1].dipmod  = 1;//htmp[l1].dip;
      }
    }
    free(templatemap);
  }
#endif

 //========================================================================
  // if ADU option is open do the histogram and transfrom the adu field in the data  
  
  double *histo_gi = (double *) malloc(32000*sizeof(double)*nbolo);
  double *histo2_gi = (double *) malloc(32000*sizeof(double)*nbolo);
  double *histon_gi = (double *) malloc(32000*sizeof(double)*nbolo);
  int *invgi       = (int *)    malloc(32000*sizeof(int)*nbolo);
  GAINSTEP = Param->GAINSTEP;
  //GAINSTEP = nadu3;
  double *xgi      = (double *) malloc(GAINSTEP*sizeof(double)*nbolo);
  double *nxgi     = (double *) malloc(GAINSTEP*sizeof(double)*nbolo);

  double mat_dip[4*15];
  double vec_dip[2*15];
    
  memset(mat_dip,0,4*sizeof(double)*nbolo);
  memset(vec_dip,0,2*sizeof(double)*nbolo);

  memset(histo_adu,0,nbolo*NADUSTEP*sizeof(double));

  memset(histo_gi ,0,nbolo*32000*sizeof(double));
  memset(histo2_gi ,0,nbolo*32000*sizeof(double));
  memset(histon_gi ,0,nbolo*32000*sizeof(double));

  if (Param->OUT_NOPOL[0]%2==0) {
    for (k=0;k<nnbpix;k++)  {
      if (flgpix[k]>0) {
        long ndata = loc_nhpix[k];
        hpix *htmp = loc_hpix[k];
        int l1;
        double DI=0,DQ=0,DU=0;
        double II2=0,IQ2=0,IU2=0,QQ2=0,QU2=0,UU2=0;
        for (l1=0;l1<ndata;l1++) {
          long ri1=htmp[l1].rg-globalBeginRing;
          if (flg_rg[htmp[l1].ib][ri1]!=0) {
            
            double CO1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].co
                                              -dpsisi[htmp[l1].ib]*htmp[l1].si);
            double SI1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].si
                                              +dpsisi[htmp[l1].ib]*htmp[l1].co);
            
            DI+=htmp[l1].w*htmp[l1].dip;
            DQ+=htmp[l1].w*htmp[l1].dip*CO1;
            DU+=htmp[l1].w*htmp[l1].dip*SI1;
            
            II2+=htmp[l1].w;
            IQ2+=htmp[l1].w*CO1;
            IU2+=htmp[l1].w*SI1;
            QQ2+=htmp[l1].w*CO1*CO1;
            QU2+=htmp[l1].w*SI1*CO1;
            UU2+=htmp[l1].w*SI1*SI1;
          }
        }
        solvemap(&DI,&DQ,&DU,II2,IQ2,IU2,QQ2,QU2,UU2);
        
        for (l1=0;l1<ndata;l1++) {
          long ri1=htmp[l1].rg-globalBeginRing;
          if (flg_rg[htmp[l1].ib][ri1]!=0) {
            double CO1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].co
                                              -dpsisi[htmp[l1].ib]*htmp[l1].si);
            double SI1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].si
                                              +dpsisi[htmp[l1].ib]*htmp[l1].co);
            
            double li=htmp[l1].w,lco=CO1*htmp[l1].w,lsi=SI1*htmp[l1].w;
            solvemap(&li,&lco,&lsi,II2,IQ2,IU2,QQ2,QU2,UU2);
            double divi=NEP_tab[htmp[l1].ib]*htmp[l1].hit*(1-(li+CO1*lco+SI1*lsi));
            divi=divi*divi;
            
            long l3;
            
            for (l3=0;l3<ndata;l3++) if (l3!=l1) {
                long ri3=htmp[l3].rg-globalBeginRing;
                if (flg_rg[htmp[l3].ib][ri3]!=0) {
               li=htmp[l3].w;lco=CO1*htmp[l3].w;lsi=SI1*htmp[l3].w;
               solvemap(&li,&lco,&lsi,II2,IQ2,IU2,QQ2,QU2,UU2);
               double sigma=NEP_tab[htmp[l3].ib]*htmp[l3].hit*(li+CO1*lco+SI1*lsi);
               divi+=sigma*sigma;
                }
              }
            double ww=1/divi;
            double tmp=(htmp[l1].dip); //-(DI+CO1*DQ+SI1*DU));
            
            histo_gi[htmp[l1].rg+htmp[l1].ib*32000]+=ww*tmp;
            histo2_gi[htmp[l1].rg+htmp[l1].ib*32000]+=ww*tmp*tmp;
            histon_gi[htmp[l1].rg+htmp[l1].ib*32000]+=ww;
            
            histo_adu[htmp[l1].adu+htmp[l1].ib*NADUSTEP]+=ww;

            mat_dip[0+4*htmp[l1].ib]+=ww*htmp[l1].dip*htmp[l1].dip;
            mat_dip[1+4*htmp[l1].ib]+=ww*htmp[l1].dip;
            mat_dip[2+4*htmp[l1].ib]+=ww*htmp[l1].dip;
            mat_dip[3+4*htmp[l1].ib]+=ww;
            
            vec_dip[0+2*htmp[l1].ib]+=ww*htmp[l1].listofpix[0]*htmp[l1].dip;
            vec_dip[1+2*htmp[l1].ib]+=ww*htmp[l1].listofpix[0];
            
          }
        }
      }
    }
  }// NO POL METHOD
  else {
    for (k=0;k<nnbpix;k++)  {
      if (flgpix[k]>0) {
        long ndata = loc_nhpix[k];
        hpix *htmp = loc_hpix[k];
        int l1;
        double DI=0;
        double II2=0;
        for (l1=0;l1<ndata;l1++) {
          long ri1=htmp[l1].rg-globalBeginRing;
          if (flg_rg[htmp[l1].ib][ri1]!=0) {
            
            
            DI+=htmp[l1].w*htmp[l1].dip;
            
            II2+=htmp[l1].w;
          }
        }
        DI=DI/II2;
        
        for (l1=0;l1<ndata;l1++) {
          long ri1=htmp[l1].rg-globalBeginRing;
          if (flg_rg[htmp[l1].ib][ri1]!=0) {

            double divi=NEP_tab[htmp[l1].ib]*htmp[l1].hit;
            
            divi=divi*divi;
            double ww=1/divi;
            
            double tmp=(htmp[l1].dip); 
            
            histo_gi[ htmp[l1].rg+htmp[l1].ib*32000]+=ww*tmp;
            histo2_gi[htmp[l1].rg+htmp[l1].ib*32000]+=ww*tmp*tmp;
            histon_gi[htmp[l1].rg+htmp[l1].ib*32000]+=ww;
            
            histo_adu[htmp[l1].adu+htmp[l1].ib*NADUSTEP]+=ww;

            if (htmp[l1].rg==300) {
              mat_dip[0+4*htmp[l1].ib]+=ww*htmp[l1].dip*htmp[l1].dip;
              mat_dip[1+4*htmp[l1].ib]+=ww*htmp[l1].dip;
              mat_dip[2+4*htmp[l1].ib]+=ww*htmp[l1].dip;
              mat_dip[3+4*htmp[l1].ib]+=ww;
              
              vec_dip[0+2*htmp[l1].ib]+=ww*(htmp[l1].adu)*htmp[l1].dip;
              vec_dip[1+2*htmp[l1].ib]+=ww*(htmp[l1].adu);
            }
          }
        }
      }
    }
    
  }


  /*===============================================================================================
    COMPUTE CUTTING RING FOR OPTIMAL MAP-MAKING
    ===============================================================================================*/

  if (CUTRG>1) {
    int l_rank;
    for (l_rank=0;l_rank<mpi_size;l_rank++) {
      long nring=globalRankInfo.EndRing[l_rank]-globalRankInfo.BeginRing[l_rank]+1;
      PIOLONG * histo_phase = (PIOLONG *) malloc(sizeof(PIOLONG)*NBPH*nbolo*nring); 
      PIOLONG * l_histo_phase = (PIOLONG *) malloc(sizeof(PIOLONG)*NBPH*nbolo*nring);
      int l1;
      memset(histo_phase,0,sizeof(PIOLONG)*NBPH*nbolo*nring); 
      for (k=0;k<nnbpix;k++)  {
        if (flgpix[k]>0) {
          long ndata = loc_nhpix[k];
          hpix *htmp = loc_hpix[k];
          for (l1=0;l1<ndata;l1++) {
            long ri1=htmp[l1].rg-globalBeginRing;
            if (flg_rg[htmp[l1].ib][ri1]!=0) {
              if (htmp[l1].rg>=globalRankInfo.BeginRing[l_rank]&&htmp[l1].rg<=globalRankInfo.EndRing[l_rank]) {
                histo_phase[htmp[l1].hrg+(htmp[l1].rg-globalRankInfo.BeginRing[l_rank]+htmp[l1].ib*nring)*NBPH]+=ndata-1;
              }
            }
          }
        }
      }
      if (rank==0) {
        for (rrk=1;rrk<mpi_size;rrk++) {
          MPI_Recv(l_histo_phase,sizeof(double)*NBPH*nbolo*nring, MPI_BYTE, rrk,351, MPI_COMM_WORLD,&statu);
          for (l1=0;l1<NBPH*nbolo*nring;l1++) histo_phase[l1]+=l_histo_phase[l1];
        }
        for (ib=0;ib<nbolo;ib++) {
          for (j=0;j<nring;j++) {
            for (l1=1;l1<NBPH;l1++) histo_phase[l1+(j+ib*nring)*NBPH]+=histo_phase[l1-1+(j+ib*nring)*NBPH];
            if (histo_phase[NBPH-1+(j+ib*nring)*NBPH]==0) {
              for (l1=0;l1<NBPH;l1++) histo_phase[l1+(j+ib*nring)*NBPH]=(l1*CUTRG)/NBPH;
            }
            else {
              for (l1=0;l1<NBPH;l1++) histo_phase[l1+(j+ib*nring)*NBPH]=(histo_phase[l1+(j+ib*nring)*NBPH]*CUTRG)/histo_phase[NBPH-1+(j+ib*nring)*NBPH];
            }
            for (l1=0;l1<NBPH;l1++) if (histo_phase[l1+(j+ib*nring)*NBPH]>CUTRG-1) histo_phase[l1+(j+ib*nring)*NBPH]=CUTRG-1;
          }
        }
        
        fprintf(stderr,"Ring %d %d Done\n",(int) globalRankInfo.BeginRing[l_rank], (int) globalRankInfo.EndRing[l_rank]);
        
      }
      else {
        MPI_Send(histo_phase, sizeof(double)*NBPH*nbolo*nring, MPI_BYTE, 0, 351, MPI_COMM_WORLD);
      }
        
      MPI_Bcast(histo_phase, sizeof(double)*NBPH*nbolo*nring, MPI_BYTE, 0, MPI_COMM_WORLD);
        
      for (k=0;k<nnbpix;k++)  {
        long ndata = loc_nhpix[k];
        hpix *htmp = loc_hpix[k];
        int l1;
        for (l1=0;l1<ndata;l1++) {
          long ri1=htmp[l1].rg-globalBeginRing;
          if (flg_rg[htmp[l1].ib][ri1]!=0) {
            if (htmp[l1].rg>=globalRankInfo.BeginRing[l_rank]&&htmp[l1].rg<=globalRankInfo.EndRing[l_rank]) {
              htmp[l1].hrg=htmp[l1].rg*CUTRG+histo_phase[htmp[l1].hrg+(htmp[l1].rg-globalRankInfo.BeginRing[l_rank]+htmp[l1].ib*nring)*NBPH];
            }
          }
        }
      }
      free(histo_phase);
      free(l_histo_phase);
    }
  }

  //============================================================
  // remove dipole from the first harmonic
  //
  
  if (rank==0) {
    MPI_Status statu;
    double l_mat[4*14];
    double l_vec[2*14];
    
    for (rrk=1;rrk<mpi_size;rrk++) {
      MPI_Recv(l_mat,sizeof(double)*4*nbolo, MPI_BYTE, rrk,71, MPI_COMM_WORLD,&statu);
      MPI_Recv(l_vec,sizeof(double)*2*nbolo, MPI_BYTE, rrk,71, MPI_COMM_WORLD,&statu);
      for (j=0;j<4*nbolo;j++) mat_dip[j]+=l_mat[j];
      for (j=0;j<2*nbolo;j++) vec_dip[j]+=l_vec[j];
    }

    for (j=0;j<nbolo;j++) {
      lusol(mat_dip+4*j,vec_dip+2*j,2);

      fprintf(stderr,"GAIN DIP BOL %ld : %lf\n",j,vec_dip[2*j]);
    }
  }
  else  {
    MPI_Send(mat_dip, sizeof(double)*4*nbolo, MPI_BYTE, 0, 71, MPI_COMM_WORLD);
    MPI_Send(vec_dip, sizeof(double)*2*nbolo, MPI_BYTE, 0, 71, MPI_COMM_WORLD);
  }


  if (Param->flag_ADU==_PAR_TRUE) {
    for (ib=0;ib<nbolo;ib++) {
      if (rank==0) {
        MPI_Status statu;
        int rrk;
        double *l_histo = (double *) malloc((NADUSTEP)*sizeof(double));
        for (rrk=1;rrk<mpi_size;rrk++) {
          MPI_Recv(l_histo,sizeof(double)*NADUSTEP, MPI_BYTE, rrk,51, MPI_COMM_WORLD,&statu);
          for (j=0;j<NADUSTEP;j++) histo_adu[j+ib*NADUSTEP]+=l_histo[j];
          MPI_Recv(l_histo,sizeof(double)*32000, MPI_BYTE, rrk,54, MPI_COMM_WORLD,&statu);
          for (j=0;j<32000;j++) histo_gi[j+ib*32000]+=l_histo[j];
          MPI_Recv(l_histo,sizeof(double)*32000, MPI_BYTE, rrk,55, MPI_COMM_WORLD,&statu);
          for (j=0;j<32000;j++) histo2_gi[j+ib*32000]+=l_histo[j];
          MPI_Recv(l_histo,sizeof(double)*32000, MPI_BYTE, rrk,56, MPI_COMM_WORLD,&statu);
          for (j=0;j<32000;j++) histon_gi[j+ib*32000]+=l_histo[j];
        }
        free(l_histo);

        PIOSTRING saveg;
        int tmpi=0;

        memset(xadu+ib*GAINSTEPADU,0,sizeof(double)*GAINSTEPADU);
        memset(yadu+ib*GAINSTEPADU,0,sizeof(double)*GAINSTEPADU);
        memset(nxadu+ib*GAINSTEPADU,0,sizeof(double)*GAINSTEPADU);
        double cumul=0;

        for (j=0;j<4096*4096;j++) {
          int xj=getx12(j);
          int yj=gety12(j);
          if (xj<ADURGSTEP&&yj<ADUSTEP) {
            histo_adu[xj+yj*ADURGSTEP+ib*NADUSTEP]+=cumul;
            cumul=histo_adu[xj+yj*ADURGSTEP+ib*NADUSTEP];
          }
        }
        double step=(cumul)/(GAINSTEPADU);
        cumul=0;
        for (j=0;j<4096*4096;j++) {
          int xj=getx12(j);
          int yj=gety12(j);
          if (xj<ADURGSTEP&&yj<ADUSTEP) {
            invadu[xj+yj*ADURGSTEP+ib*NADUSTEP]=tmpi;
            double delta=histo_adu[xj+yj*ADURGSTEP+ib*NADUSTEP]-cumul;
            xadu[tmpi+ib*GAINSTEPADU]+=(xj)*(delta);
            yadu[tmpi+ib*GAINSTEPADU]+=(yj)*(delta);
            nxadu[tmpi+ib*GAINSTEPADU]+=(delta);
            cumul=histo_adu[xj+yj*ADURGSTEP+ib*NADUSTEP];
            if (histo_adu[xj+yj*ADURGSTEP+ib*NADUSTEP]>step*(1+tmpi)) {
              tmpi++;
            }
          }
        }
        for (i=0;i<GAINSTEPADU;i++) {
          xadu[i+ib*GAINSTEPADU]/=nxadu[i+ib*GAINSTEPADU];
          yadu[i+ib*GAINSTEPADU]/=nxadu[i+ib*GAINSTEPADU];
        }
        
        sprintf(saveg,"%s_HISTO_ADU",Param->Out_Offset[ib]);
        fprintf(stderr,"Write HISTO_ADU  %lld\n",(long long) (PIOWriteVECT(saveg,histo_adu+ib*NADUSTEP,0,sizeof(PIODOUBLE)*NADUSTEP)/sizeof(PIODOUBLE)));

        sprintf(saveg,"%s_INV_ADU",Param->Out_Offset[ib]);
        fprintf(stderr,"Write INV_ADU  %lld\n",(long long) (PIOWriteVECT(saveg,invadu+ib*NADUSTEP,0,sizeof(int)*NADUSTEP)/sizeof(int)));

        sprintf(saveg,"%s_XADU",Param->Out_Offset[ib]);
        fprintf(stderr,"Write XADU  %lld\n",(long long) (PIOWriteVECT(saveg,xadu+ib*GAINSTEPADU,0,sizeof(PIODOUBLE)*GAINSTEPADU)/sizeof(PIODOUBLE)));
        
        sprintf(saveg,"%s_YADU",Param->Out_Offset[ib]);
        fprintf(stderr,"Write YADU  %lld\n",(long long) (PIOWriteVECT(saveg,yadu+ib*GAINSTEPADU,0,sizeof(PIODOUBLE)*GAINSTEPADU)/sizeof(PIODOUBLE)));
        
        sprintf(saveg,"%s_HISTO_GAIN",Param->Out_Offset[ib]);
        for (i=0;i<32000;i++) if (histon_gi[i+ib*32000]>0) histo_gi[i+ib*32000]=
                                                             histo2_gi[i+ib*32000]
                                                             -histo_gi[i+ib*32000]*histo_gi[i+ib*32000]/histon_gi[i+ib*32000];
        fprintf(stderr,"Write HISTO_GAIN  %lld\n",(long long) (PIOWriteVECT(saveg,histo_gi+ib*32000,0,32000*sizeof(double))/sizeof(double)));
        for (i=1;i<32000;i++) histo_gi[i+ib*32000]+=histo_gi[i-1+ib*32000];

        step=(histo_gi[31999+ib*32000]+1)/GAINSTEP;

        memset(xgi+ib*GAINSTEP,0,sizeof(double)*GAINSTEP);
        memset(nxgi+ib*GAINSTEP,0,sizeof(double)*GAINSTEP);
        invgi[0+ib*32000]=0;
        for (i=1;i<32000;i++) {
          invgi[i+ib*32000]=tmpi;
          xgi[tmpi+ib*GAINSTEP]+=(i)*(histo_gi[i+ib*32000]-histo_gi[i-1+ib*32000]);
          nxgi[tmpi+ib*GAINSTEP]+=(histo_gi[i+ib*32000]-histo_gi[i-1+ib*32000]);
          if (histo_gi[i+ib*32000]>step*(1+tmpi)) {
            if (tmpi<GAINSTEP-1) tmpi++;
          }
        }
        for (i=0;i<32000;i++) {
          if (invgi[i+ib*32000]>=GAINSTEP) invgi[i+ib*32000]=GAINSTEP-1;
        }
          
        for (i=0;i<GAINSTEP;i++) {
          xgi[i+ib*GAINSTEP]/=nxgi[i+ib*GAINSTEP];
        }
         
        sprintf(saveg,"%s_INV_GI",Param->Out_Offset[ib]);
        fprintf(stderr,"Write INV_GI  %lld\n",(long long) PIOWriteVECT(saveg,invgi+ib*32000,0,32000*sizeof(int))/sizeof(int));
        
        sprintf(saveg,"%s_XGI",Param->Out_Offset[ib]);
        fprintf(stderr,"Write XGI  %lld\n",(long long) PIOWriteVECT(saveg,xgi+ib*GAINSTEP,0,sizeof(PIODOUBLE)*GAINSTEP)/sizeof(PIODOUBLE));
      }
      
      else {
        MPI_Send(histo_adu+ib*NADUSTEP, sizeof(double)*NADUSTEP, MPI_BYTE, 0, 51, MPI_COMM_WORLD);
        MPI_Send(histo_gi+ib*32000, sizeof(double)*32000, MPI_BYTE, 0, 54, MPI_COMM_WORLD);
        MPI_Send(histo2_gi+ib*32000, sizeof(double)*32000, MPI_BYTE, 0, 55, MPI_COMM_WORLD);
        MPI_Send(histon_gi+ib*32000, sizeof(double)*32000, MPI_BYTE, 0, 56, MPI_COMM_WORLD);
      }
    }
    
    MPI_Bcast(invadu,sizeof(int)*NADUSTEP*nbolo, MPI_BYTE, 0, MPI_COMM_WORLD);
    MPI_Bcast(xadu,sizeof(double)*GAINSTEPADU*nbolo, MPI_BYTE, 0, MPI_COMM_WORLD);
    MPI_Bcast(yadu,sizeof(double)*GAINSTEPADU*nbolo, MPI_BYTE, 0, MPI_COMM_WORLD);
    
    MPI_Bcast(invgi,sizeof(int)*32000*nbolo, MPI_BYTE, 0, MPI_COMM_WORLD);
    MPI_Bcast(xgi,sizeof(double)*GAINSTEP*nbolo, MPI_BYTE, 0, MPI_COMM_WORLD);
    
    for (k=0;k<nnbpix;k++)  {
      long ndata = loc_nhpix[k];
      hpix *htmp = loc_hpix[k];
      int l1;
      for (l1=0;l1<ndata;l1++) {
        htmp[l1].gi=invgi[htmp[l1].gi+32000*htmp[l1].ib];
        htmp[l1].adu=invadu[htmp[l1].adu+NADUSTEP*htmp[l1].ib];
      }
    }
  }

  free(invadu);
  free(histo_adu);
  free(histo_gi);
  free(histo2_gi);
  free(histon_gi);

  SSI = (double *) malloc(sizeof(double)*nnbpix);
  SSQ = (double *) malloc(sizeof(double)*nnbpix);
  SSU = (double *) malloc(sizeof(double)*nnbpix);
  SSI2 = (double *) malloc(sizeof(double)*nnbpix);
  SSQ2 = (double *) malloc(sizeof(double)*nnbpix);
  SSU2 = (double *) malloc(sizeof(double)*nnbpix);
  II = (double *) malloc(sizeof(double)*nnbpix);
  IQ = (double *) malloc(sizeof(double)*nnbpix);
  IU = (double *) malloc(sizeof(double)*nnbpix);
  QQ = (double *) malloc(sizeof(double)*nnbpix);
  UU = (double *) malloc(sizeof(double)*nnbpix);
  QU = (double *) malloc(sizeof(double)*nnbpix);
  
  dthetai = (double *) malloc(sizeof(double)*nbolo*nnbpix);
  dthetaq = (double *) malloc(sizeof(double)*nbolo*nnbpix);
  dthetau = (double *) malloc(sizeof(double)*nbolo*nnbpix);
  
  dii = (double *) malloc(sizeof(double)*nbolo*GAINSTEP); //*nnbpix);
  dqq = (double *) malloc(sizeof(double)*nbolo*GAINSTEP); //*nnbpix);
  duu = (double *) malloc(sizeof(double)*nbolo*GAINSTEP); //*nnbpix);
  
  //ddegi = (double *) malloc(sizeof(double)*nbolo*DEGREE*nnbpix);
  //ddegq = (double *) malloc(sizeof(double)*nbolo*DEGREE*nnbpix);
  //ddegu = (double *) malloc(sizeof(double)*nbolo*DEGREE*nnbpix);
  
  dcoi = (double *) malloc(sizeof(double)*nbolo*nnbpix);
  dcoq = (double *) malloc(sizeof(double)*nbolo*nnbpix);
  dcou = (double *) malloc(sizeof(double)*nbolo*nnbpix);
  
  dfri = (double *) malloc(sizeof(double)*nbolo*nnbpix);
  dfrq = (double *) malloc(sizeof(double)*nbolo*nnbpix);
  dfru = (double *) malloc(sizeof(double)*nbolo*nnbpix);
  
  ddusti = (double *) malloc(sizeof(double)*nbolo*nnbpix);
  ddustq = (double *) malloc(sizeof(double)*nbolo*nnbpix);
  ddustu = (double *) malloc(sizeof(double)*nbolo*nnbpix);
  
  dpixi = (double *) malloc(sizeof(double)*(npixbeam+1)*nbolo*nnbpix);
  dpixq = (double *) malloc(sizeof(double)*(npixbeam+1)*nbolo*nnbpix);
  dpixu = (double *) malloc(sizeof(double)*(npixbeam+1)*nbolo*nnbpix);
  
  cdip= (double *) malloc(sizeof(double)*nbolo*GAINSTEP);
  cco= (double *) malloc(sizeof(double)*nbolo);
  ccfree= (double *) malloc(sizeof(double)*nbolo);
  ctheta= (double *) malloc(sizeof(double)*nbolo);
  cdust= (double *) malloc(sizeof(double)*nbolo);
  cpix= (double *) malloc(sizeof(double)*nbolo*(npixbeam+1));
  ctmp= (double *) malloc(sizeof(double)*(newnr[nbolo]+nbolo*GAINSTEP));



  long l_off_nmatpix=0;
  {
    long rrk=0;
    for (rrk=0;rrk<mpi_size;rrk++) {
      long l_trans_nmatpix;
      if (rrk==rank) l_trans_nmatpix=l_nmatpix;
      MPI_Bcast(&l_trans_nmatpix,sizeof(long), MPI_BYTE, rrk, MPI_COMM_WORLD);
      if (rank>rrk) l_off_nmatpix+=l_trans_nmatpix;
    }
  }
  
  if (rank==0) fprintf(stderr,"BEG-ED %d %ld %ld - %ld %ld\n",rank,(long) nmatpix,(long) begpix[rank],(long) edpix[rank],
          (long) l_nmatpix);
  MPI_Barrier(MPI_COMM_WORLD);

  //double *tvec;
  //double *vvec;
  //double *totmat;

  long l1;
  
  // decoupe vecteur matrice par proc
#ifdef DORGG
  double *gtab2 = (double *) calloc (newnr[nbolo], sizeof (double));
#else
  double *gtab2 = (double *) calloc (GAINSTEP*nbolo, sizeof (double));
#endif


  long maxsizemat=newnr[nbolo]+nbolo*(GAINSTEPADU+npixbeam+1)+nmatco+nmatdust+nfreefree;
  if (GAINSTEP>GAINSTEPADU) maxsizemat=newnr[nbolo]+nbolo*(GAINSTEP+npixbeam+1)+nmatco+nmatdust+nfreefree;
  if (newnr2[nbolo]>maxsizemat&&CUTRG>(1)) {
    maxsizemat=newnr2[nbolo];
  }
  if (rank==0) fprintf(stderr,"MAXSIZE %ld %ld\n",(long) maxsizemat,(long) newnr[nbolo]+nbolo*(GAINSTEP+npixbeam)+nmatco+nmatdust+nfreefree);

  
  //vvec=(double *) malloc((maxsizemat)*sizeof (double));
  //tvec=(double *) malloc((maxsizemat)*sizeof (double));

  if (rank==0) fprintf(stderr,"NbGain %ld\n",newnr2[nbolo]);

  
    //if (Param->CALCODUST==1) {
    //  totmat=(double *) malloc(8*sizeof(double)*nnbpix);
    //}
    qmat=(double *) malloc(sizeof(double)*nnbpix);
    umat=(double *) malloc(sizeof(double)*nnbpix);
    
#ifdef COMATOUT
    double *rcomat=(double *) malloc(sizeof(double)*nnbpix);
    double *co_mat=(double *) malloc(sizeof(double)*nnbpix);
    double *co_qmat=(double *) malloc(sizeof(double)*nnbpix);
    double *co_umat=(double *) malloc(sizeof(double)*nnbpix);
#endif
    
    
#ifdef DUSTMATOUT
    double *dust_mat=(double *) malloc(sizeof(double)*nnbpix);
    double *dust_qmat=(double *) malloc(sizeof(double)*nnbpix);
    double *dust_umat=(double *) malloc(sizeof(double)*nnbpix);
#endif
    x2 =     (double *) malloc (maxsizemat*sizeof (double));
    x2old =  (double *) malloc (maxsizemat*sizeof (double));
    x2init = (double *) malloc (maxsizemat*sizeof (double));
    b2 =     (double *) malloc (maxsizemat*sizeof (double));
    d2 =     (double *) malloc (maxsizemat*sizeof (double));
    q2 =     (double *) malloc (maxsizemat*sizeof (double));
    r2 =     (double *) malloc (maxsizemat*sizeof (double));
    s2 =     (double *) malloc (maxsizemat*sizeof (double));
    hit2 =   (double *) malloc (maxsizemat*sizeof (double));


    //if (Param->TESTPOL==-1) {
    //  for (i=0;i<Param->n_MAP;i++) {
    //  sprintf(mapout[i],"%s_I%d",Param->MAP[i],(int) nseed);
    //  if (rank==0) fprintf(stderr,"MAP %s\n",mapout[i]);
    //  }
    //}
    //else 
    for (i=0;i<Param->n_MAP;i++) strcpy(mapout[i],Param->MAP[i]);
    
  
    if (rank==0) fprintf(stderr,"BEG %d %ld %ld - %ld %ld\n",rank,(long) nmatpix,(long) begpix[rank],(long) edpix[rank],
          (long) l_nmatpix);

    MPI_Barrier(MPI_COMM_WORLD);


  /*======================================================
    =
    =    LOOP ON SEEDs
    =
    =*/

  for (int iter = 0; iter < number_of_iterations; iter++) {
    for (k=0;k<nnbpix;k++)  {
      long ndata = loc_nhpix[k];
      hpix *htmp = loc_hpix[k];
      int l1;
      double l_th,l_ph;
      pix2ang_ring(Nside,k+begpix[rank],&l_th,&l_ph);
      for (l1=0;l1<ndata;l1++) {
        htmp[l1].sig = htmp[l1].listp[iter];
      }
    }
  
    itbogo=0;
  
    /*======================================================
      =
      =    solve matrix
      =
      =*/
  
    // ca plante pas
    memset(x2old  ,0,(maxsizemat)*sizeof (double));
  
    /*======================================================
      =
      =    PUT to 0
      =
      =*/
  
    // ca plante
    //if (Param->TESTPOL==-1) {
    //  for (k=0;k<nnbpix;k++)  {
    //    long ndata = loc_nhpix[k];
    //    hpix *htmp = loc_hpix[k];
    //    int l1;
    //    for (l1=0;l1<ndata;l1++) {
    //    htmp[l1].sig  = htmp[l1].dip+
    //      Param->NEP[htmp[l1].ib]/Param->Calibration[htmp[l1].ib]
    //      *sqrt(-2*log( drand48()))*cos(2*M_PI*drand48());
    //    }
    //  }
    //}
    // ca plante
    docutrg=nitbogo-1;
  
    ittt=1;
    g=gtab2;
    memset(x2     ,0,maxsizemat*sizeof (double));
    memset(x2old  ,0,maxsizemat*sizeof (double));
    memset(x2init ,0,maxsizemat*sizeof (double));
    memset(b2     ,0,maxsizemat*sizeof (double));
    memset(d2     ,0,maxsizemat*sizeof (double));
    memset(q2     ,0,maxsizemat*sizeof (double));
    memset(r2     ,0,maxsizemat*sizeof (double));
    memset(s2     ,0,maxsizemat*sizeof (double));
    memset(hit2   ,0,maxsizemat*sizeof (double));
  #if DORGG
    for (i=0;i<newnr[nbolo];i++) g[i]=1.;
  #else
    for (i=0;i<GAINSTEP*nbolo;i++) g[i]=1.;
  #endif
  
    nadu=GAINSTEPADU;
    double *reshisto=(double *) malloc(sizeof(double)*GAINSTEPADU*nbolo);
    memset(reshisto,0,sizeof(double)*GAINSTEPADU*nbolo);
  
    if (rank==0) {
      fprintf(stderr,"NADU %ld\n", (long) nadu);
    }
    PrintFreeMemOnNodes( rank, mpi_size, "before solving matrix");
  
    for (k=0;k<nnbpix;k++)  {
      long ndata = loc_nhpix[k];
      hpix *htmp = loc_hpix[k];
      for (l1=0;l1<ndata;l1++) {
        reshisto[htmp[l1].adu+htmp[l1].ib*GAINSTEPADU]+=htmp[l1].w;
        if (htmp[l1].adu<0||htmp[l1].adu>GAINSTEPADU-1) {
          fprintf(stderr,"ERRER .ADU\n");
          exit(0);
        }
      }
    }
  
    if (rank==0) {
      double *l_reshisto=(double *) malloc(sizeof(double)*nadu*nbolo);
      MPI_Status statu;
      int rrk;
      for (rrk=1;rrk<mpi_size;rrk++) {
        MPI_Recv(l_reshisto,sizeof(double)*nadu*nbolo, MPI_BYTE, rrk,751, MPI_COMM_WORLD,&statu);
        for (j=0;j<nadu*nbolo;j++) reshisto[j]+=l_reshisto[j];
      }
      for (j=0;j<nadu;j++) {
        fprintf(stderr,"RESHIST ");
        for (i=0;i<nbolo;i++) {
          fprintf(stderr,"%lg\t",reshisto[j+nadu*i]);
        }
        fprintf(stderr,"\n");
      }
      free(l_reshisto);
    }
    else {
      MPI_Send(reshisto, sizeof(double)*nadu*nbolo, MPI_BYTE, 0, 751, MPI_COMM_WORLD);
    }
    
    free(reshisto);
  
    MPI_Barrier(MPI_COMM_WORLD);
    
    // Exchange dip against sig in the XI2
    
  
    for (k=0;k<nnbpix;k++)  {
      long ndata = loc_nhpix[k];
      hpix *htmp = loc_hpix[k];
      II[k]=0;
      IQ[k]=0;
      IU[k]=0;
      QQ[k]=0;
      UU[k]=0;
      QU[k]=0;
      
      for (l1=0;l1<ndata;l1++) {
        long ri1=htmp[l1].rg-globalBeginRing;
        
        double CO1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].co
                                          -dpsisi[htmp[l1].ib]*htmp[l1].si);
        double SI1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].si
                                          +dpsisi[htmp[l1].ib]*htmp[l1].co);
        if (flg_rg[htmp[l1].ib][ri1]!=0) {
          II[k]+=htmp[l1].w;
          IQ[k]+=htmp[l1].w*CO1;
          IU[k]+=htmp[l1].w*SI1;
          QQ[k]+=htmp[l1].w*CO1*CO1;
          QU[k]+=htmp[l1].w*SI1*CO1;
          UU[k]+=htmp[l1].w*SI1*SI1;
        }
      }
      if (ndata>0&&flgpix[k]>0) {
        
        double SI=0;
        double SQ=0;
        double SU=0;
        for (l1=0;l1<ndata;l1++) {
          long ri1=htmp[l1].rg-globalBeginRing;
          
          double CO1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].co
                                            -dpsisi[htmp[l1].ib]*htmp[l1].si);
          double SI1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].si
                                            +dpsisi[htmp[l1].ib]*htmp[l1].co);
          if (flg_rg[htmp[l1].ib][ri1]!=0) {
            
            SI+=htmp[l1].w*(htmp[l1].dip+htmp[l1].fsl); 
            SQ+=htmp[l1].w*(htmp[l1].dip+htmp[l1].fsl)*CO1; 
            SU+=htmp[l1].w*(htmp[l1].dip+htmp[l1].fsl)*SI1; 
  
  
  #ifdef UPDATE_DIP
            double tmp=htmp[l1].sig-htmp[l1].fsl;
            matdip[iri1*4]  +=htmp[l1].w*htmp[l1].dip*htmp[l1].dip;
            matdip[iri1*4+1]+=htmp[l1].w*htmp[l1].dip;
            matdip[iri1*4+2]+=htmp[l1].w*htmp[l1].dip;
            matdip[iri1*4+3]+=htmp[l1].w;
            vecdip[iri1*2]  +=htmp[l1].w*tmp*htmp[l1].dip;
            vecdip[iri1*2+1]+=htmp[l1].w*tmp;
  #endif
  
          }
        }
        
        solvemap(&SI,&SQ,&SU,II[k],IQ[k],IU[k],QQ[k],QU[k],UU[k]);
        
        SSI2[k]=SI;
        SSQ2[k]=SQ;
        SSU2[k]=SU;
      }
    }
  
    
  
  
    double *gain = (double *) malloc(sizeof(double)*nbolo*GAINSTEP);
    double *gainadu = (double *) malloc(sizeof(double)*nbolo*GAINSTEPADU);
  
    for (i=0;i<nbolo;i++) {
      for (j=0;j<GAINSTEP;j++) gain[i*GAINSTEP+j]=1; //+3E-3*i+1E-3*cos(j/4.);
      for (j=0;j<GAINSTEPADU;j++) gainadu[i*GAINSTEPADU+j]=0;
    }
  
    memset(x2,0,(newnr[nbolo]+nbolo*npixbeam+nmatco+nmatdust+nfreefree)*sizeof(double));
    
    double *xoff= (double *) malloc(sizeof(double)*newnr2[nbolo]);
    
    gainoff=0;
    //frprmn(gain,nbolo*GAINSTEP,ftol,&iter,&fret,funcxigain,dfuncxigain);
    //dfpmin(gain,nbolo*GAINSTEP,ftol,&iter,&fret,funcxigain,dfuncxigain);
    
    //minimize(x2,gain);
  
    //if (CUTRG>1) minimize_optimize(x2,xoff,gain);
  
    
    if (GAINSTEP<GAINSTEPADU) nmatres=newnr[nbolo]+nbolo*(GAINSTEPADU+npixbeam)+nmatdust+nmatco+nfreefree;
    else nmatres=newnr[nbolo]+nbolo*(GAINSTEP+npixbeam)+nmatdust+nmatco+nfreefree;
    
    double *x3= (double *) malloc(sizeof(double)*(nmatres+nbolo)); //+nbolo in case of fittetha
  
    GetProcMem(&vmem,&phymem);
    if (rank==0) fprintf(stderr,"Rank: %ld[%d] MEM %.1lf[%.1lf]MB line=%d\n",
                (long) rank, getpid(),
                (double) vmem/1024./1024.,
                (double) phymem/1024./1024.,__LINE__);
    int itt;
    memset(x3,0,sizeof(double)*nmatres);
  
    if (DODISTOR!=0) { 
      //=========================================================
      GAINSTEP=1;
      for (k=0;k<nnbpix;k++) {
        long ndata = loc_nhpix[k];
        hpix *htmp = loc_hpix[k];
        for (l1=0;l1<ndata;l1++) {
          htmp[l1].hrg=htmp[l1].gi;
          htmp[l1].gi=0;
          //htmp[l1].sig=htmp[l1].dip+htmp[l1].fsl+1E-7*(rgord[htmp[l1].ib][htmp[l1].rg-globalBeginRing])
          //+exp(-(yadu[htmp[l1].adu]-1600)*(yadu[htmp[l1].adu]-1600)*0.0001)*1E-5;
        }
      }
    }
    itt=0;
    double resxi=1;
    while ((itt<DODISTOR&&sqrt(resxi)>1E-10)||(DODISTOR==0&&itt<7)) {
      //minimize(x2,gain);
      memset(newnr[nbolo]+x3,0,sizeof(double)*nbolo*GAINSTEP);
  
      if (Param->FITTHETA==1&&itt>0) {  
        for (k=0;k<nnbpix;k++) {
          long ndata = loc_nhpix[k];
          hpix *htmp = loc_hpix[k];
          
          for (l1=0;l1<ndata;l1++) {
            htmp[l1].listofpix[npixbeam-1]=htmp[l1].model/eta_dest[htmp[l1].ib];
          }
        }
      }
  
      if (Param->OUT_NOPOL[0]%2==1) minimize_gain_nopol(x3,gainadu,gain);
      else minimize_gain_gi(x3,gainadu,gain);
  
      resxi=0;
      for (i=0;i<nbolo;i++){
        for (j=0;j<GAINSTEP;j++) {
          resxi+=x3[newnr[nbolo]+i*GAINSTEP+j]*x3[newnr[nbolo]+i*GAINSTEP+j];
          gain[i*GAINSTEP+j]-=x3[newnr[nbolo]+i*GAINSTEP+j];
        }
      }
      resxi/=(nbolo*GAINSTEP);
  
      if (rank==0) {
        fprintf(stderr,"GI XIGAIN %.10lg\n",sqrt(resxi));
        if (Param->FITTHETA==1&&itt>0) { 
          for (i=GAINSTEP*nbolo+newnr[nbolo];i<(GAINSTEP+npixbeam-1)*nbolo+newnr[nbolo];i++) {
            if ((i-newnr[nbolo])%nbolo==0) fprintf(stderr,"X3 ");
            fprintf(stderr,"%lg ",x3[i]);
            if ((i-newnr[nbolo])%nbolo==nbolo-1) fprintf(stderr,"\n");
          } 
          for (i=(GAINSTEP+npixbeam-1)*nbolo+newnr[nbolo];i<(GAINSTEP+npixbeam)*nbolo+newnr[nbolo];i++) {
            if ((i-newnr[nbolo])%nbolo==0) fprintf(stderr,"ANG ");
            fprintf(stderr,"%lg(%.2lg) ",x3[i],asin(x3[i])/M_PI*180.);
            if ((i-newnr[nbolo])%nbolo==nbolo-1) fprintf(stderr,"\n");
          }
          for (i=(GAINSTEP+npixbeam)*nbolo+newnr[nbolo];i<(GAINSTEP+npixbeam)*nbolo+nmatco+nmatdust+nfreefree+newnr[nbolo];i++) {
            if ((i-newnr[nbolo])%nbolo==0) fprintf(stderr,"X3 ");
            fprintf(stderr,"%lg ",x3[i]);
            if ((i-newnr[nbolo])%nbolo==nbolo-1) fprintf(stderr,"\n");
          }
        }
        else {
          for (i=GAINSTEP*nbolo+newnr[nbolo];i<(GAINSTEP+npixbeam)*nbolo+nmatco+nmatdust+nfreefree+newnr[nbolo];i++) {
            if ((i-newnr[nbolo])%nbolo==0) fprintf(stderr,"X3 ");
            fprintf(stderr,"%lg ",x3[i]);
            if ((i-newnr[nbolo])%nbolo==nbolo-1) fprintf(stderr,"\n");
          }
        }
      }
      if (Param->FITTHETA==1||DODISTOR!=0) {
        double *totadu=NULL;
        double *ntotadu=NULL;
        
        if (DODISTOR!=0) {
          totadu = (double *) malloc(sizeof(double)*GAINSTEPADU*nbolo);
          ntotadu = (double *) malloc(sizeof(double)*GAINSTEPADU*nbolo);
          
          memset(totadu,0,sizeof(double)*GAINSTEPADU*nbolo);
          memset(ntotadu,0,sizeof(double)*GAINSTEPADU*nbolo);
        }
  
        for (k=0;k<nnbpix;k++) {
          long ndata = loc_nhpix[k];
          hpix *htmp = loc_hpix[k];
          
          double SII=0;
          double SIQ=0;
          double SIU=0;
          double SQQ=0;
          double SUU=0;
          double SQU=0;
          
          double SI=0;
          double SQ=0;
          double SU=0;
          
          int l1,m;
          
          for (l1=0;l1<ndata;l1++) {
            long ri1=htmp[l1].rg-globalBeginRing;
            long iri1=rgord[htmp[l1].ib][ri1]+newnr[htmp[l1].ib];
            if (flg_rg[htmp[l1].ib][ri1]!=0) {
              
              double CO1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].co
                                                -dpsisi[htmp[l1].ib]*htmp[l1].si);
              double SI1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].si
                                                +dpsisi[htmp[l1].ib]*htmp[l1].co);
              
              double gg2=gain[htmp[l1].ib];
              
              double sig_corr = htmp[l1].sig*gg2-htmp[l1].fsl;
              if (REMHDIP==0) sig_corr-=htmp[l1].dip;
              else sig_corr-=htmp[l1].freefree;
              
              sig_corr-=x3[iri1];
              
              // ATTENTION GAINSTEP EST FORCE A 1
              sig_corr-=x3[newnr[nbolo]+htmp[l1].ib]*htmp[l1].dip+gainadu[htmp[l1].adu+htmp[l1].ib*GAINSTEPADU];
              
              if (nmatco!=0) {
                sig_corr -=htmp[l1].comap*x3[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP)+htmp[l1].ib];
              }
              if (nmatdust!=0) {
                sig_corr -=htmp[l1].dustmap*x3[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP)+nmatco+htmp[l1].ib];
              }
              if (nfreefree!=0) {
                sig_corr -=htmp[l1].freefree*x3[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP)+nmatco+nmatdust+htmp[l1].ib];
              }
              for (m=0;m<npixbeam;m++)  {
                sig_corr -= htmp[l1].listofpix[m]*x3[newnr[nbolo]+nbolo*(GAINSTEP)+m*nbolo+htmp[l1].ib];
              }
              SI +=htmp[l1].w*sig_corr;
              SQ +=htmp[l1].w*CO1*sig_corr;
              SU +=htmp[l1].w*SI1*sig_corr;
              
              SII +=htmp[l1].w;
              SIQ +=htmp[l1].w*CO1;
              SIU +=htmp[l1].w*SI1;
              SQQ +=htmp[l1].w*CO1*CO1;
              SQU +=htmp[l1].w*CO1*SI1;
              SUU +=htmp[l1].w*SI1*SI1;
            }
          }
          
          if (Param->OUT_NOPOL[0]%2==0) {
            double cond=cond_3_3_thres(SII,SIQ,SIU,
                                       SIQ,SQQ,SQU,
                                       SIU,SQU,SUU);
            
            if (cond<Param->seuilcond) {
              solvemap(&SI,&SQ,&SU,SII,SIQ,SIU,SQQ,SQU,SUU);
              
              for (l1=0;l1<ndata;l1++) {
                long ri1=htmp[l1].rg-globalBeginRing;
                long iri1=rgord[htmp[l1].ib][ri1]+newnr[htmp[l1].ib];
                if (flg_rg[htmp[l1].ib][ri1]!=0) {
                  double CO1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].co
                                                    -dpsisi[htmp[l1].ib]*htmp[l1].si);
                  double SI1=eta_dest[htmp[l1].ib]*(dpsico[htmp[l1].ib]*htmp[l1].si
                                                    +dpsisi[htmp[l1].ib]*htmp[l1].co);
                  
                  
                  if (flgpix[k]>0) {
                    if (DODISTOR!=0) {
                      double gg2=gain[htmp[l1].ib];
                      
                      double sig_corr = (htmp[l1].sig*gg2-htmp[l1].fsl);
                      if (REMHDIP==0) sig_corr-=htmp[l1].dip;
                      else sig_corr-=htmp[l1].freefree;
                      
                      sig_corr-=x3[iri1];
                      
                      // ATTENTION GAINSTEP EST FORCE A 1
                      sig_corr-=x3[newnr[nbolo]+htmp[l1].ib]*htmp[l1].dip+gainadu[htmp[l1].adu+htmp[l1].ib*GAINSTEPADU];
                      
                      if (nmatco!=0) {
                        sig_corr -=htmp[l1].comap*x3[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP)+htmp[l1].ib];
                      }
                      if (nmatdust!=0) {
                        sig_corr -=htmp[l1].dustmap*x3[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP)+nmatco+htmp[l1].ib];
                      }
                      if (nfreefree!=0) {
                        sig_corr -=htmp[l1].freefree*x3[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP)+nmatco+nmatdust+htmp[l1].ib];
                      }
                      for (m=0;m<npixbeam;m++)  {
                        sig_corr -= htmp[l1].listofpix[m]*x3[newnr[nbolo]+nbolo*(GAINSTEP)+m*nbolo+htmp[l1].ib];
                      }
                      
                      totadu[htmp[l1].adu+htmp[l1].ib*GAINSTEPADU]  += htmp[l1].w*(sig_corr-SI-CO1*SQ-SI1*SU);
                      ntotadu[htmp[l1].adu+htmp[l1].ib*GAINSTEPADU] += htmp[l1].w;
                    }
                  }
                  if (Param->FITTHETA==1) htmp[l1].model=-SI1*SQ+CO1*SU;
                }
              }
            }
          } // END POL
          else {
            if (flgpix[k]>0) {
              if (DODISTOR!=0) {
                if (SII>0) {
                  SI/=SII;
                  for (l1=0;l1<ndata;l1++) {
                    long ri1=htmp[l1].rg-globalBeginRing;
                    long iri1=rgord[htmp[l1].ib][ri1]+newnr[htmp[l1].ib];
                    if (flg_rg[htmp[l1].ib][ri1]!=0) {
                      
                      double gg2=gain[htmp[l1].ib];
                      
                      double sig_corr = (htmp[l1].sig*gg2-htmp[l1].fsl);
                      if (REMHDIP==0) sig_corr-=htmp[l1].dip;
                      else sig_corr-=htmp[l1].freefree;
                      
                      sig_corr-=x3[iri1];
                      
                      // ATTENTION GAINSTEP EST FORCE A 1
                      sig_corr-=x3[newnr[nbolo]+htmp[l1].ib]*htmp[l1].dip+gainadu[htmp[l1].adu+htmp[l1].ib*GAINSTEPADU];
                      
                      if (nmatco!=0) {
                        sig_corr -=htmp[l1].comap*x3[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP)+htmp[l1].ib];
                      }
                      if (nmatdust!=0) {
                        sig_corr -=htmp[l1].dustmap*x3[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP)+nmatco+htmp[l1].ib];
                      }
                      if (nfreefree!=0) {
                        sig_corr -=htmp[l1].freefree*x3[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP)+nmatco+nmatdust+htmp[l1].ib];
                      }
                      for (m=0;m<npixbeam;m++)  {
                        sig_corr -= htmp[l1].listofpix[m]*x3[newnr[nbolo]+nbolo*(GAINSTEP)+m*nbolo+htmp[l1].ib];
                      }
                      totadu[htmp[l1].adu+htmp[l1].ib*GAINSTEPADU]  += htmp[l1].w*(sig_corr-SI);
                      ntotadu[htmp[l1].adu+htmp[l1].ib*GAINSTEPADU] += htmp[l1].w;
                    }
                  }
                }
              }
            }
          }
        }
  
        if (DODISTOR!=0) {
          resxi=0;
          
          if (rank==0) {
            double *l_offtab=(double *) malloc(sizeof(double)*GAINSTEPADU*nbolo);
            MPI_Status statu;
            int rrk;
            for (rrk=1;rrk<mpi_size;rrk++) {
              MPI_Recv(l_offtab,sizeof(double)*GAINSTEPADU*nbolo, MPI_BYTE, rrk,2951, MPI_COMM_WORLD,&statu);
              for (j=0;j<GAINSTEPADU*nbolo;j++) totadu[j]+=l_offtab[j];
              MPI_Recv(l_offtab,sizeof(double)*GAINSTEPADU*nbolo, MPI_BYTE, rrk,2952, MPI_COMM_WORLD,&statu);
              for (j=0;j<GAINSTEPADU*nbolo;j++) ntotadu[j]+=l_offtab[j];
            }
            free(l_offtab);
            for (j=0;j<GAINSTEPADU*nbolo;j++) {
              if (ntotadu[j]>0) {
                resxi+=(totadu[j]/ntotadu[j])*(totadu[j]/ntotadu[j]);
                gainadu[j]+=totadu[j]/ntotadu[j];
                totadu[j]/=ntotadu[j];
              }
            }
          }
          else {
            MPI_Send(totadu, sizeof(double)*GAINSTEPADU*nbolo, MPI_BYTE, 0, 2951, MPI_COMM_WORLD);
            MPI_Send(ntotadu, sizeof(double)*GAINSTEPADU*nbolo, MPI_BYTE, 0, 2952, MPI_COMM_WORLD);
          }
          MPI_Bcast(gainadu,sizeof(double)*GAINSTEPADU*nbolo, MPI_BYTE, 0, MPI_COMM_WORLD);
        
          if (rank==0) {
            fprintf(stderr,"GI XIADU %.10lg\n",sqrt(resxi));
          }
        }
        if (itt==0&&Param->FITTHETA==1) {
          npixbeam++;
          DOFITTHETA=1; //do not average angle
          if (GAINSTEP<GAINSTEPADU) nmatres=newnr[nbolo]+nbolo*(GAINSTEPADU+npixbeam)+nmatdust+nmatco+nfreefree;
          else nmatres=newnr[nbolo]+nbolo*(GAINSTEP+npixbeam)+nmatdust+nmatco+nfreefree;
          memset(x3,0,sizeof(double)*nmatres);
          //for (i=0;i<nbolo;i++) {
          //  for (j=0;j<GAINSTEP;j++) gain[i*GAINSTEP+j]=1; //+3E-3*i+1E-3*cos(j/4.);
          //  for (j=0;j<GAINSTEPADU;j++) gainadu[i*GAINSTEPADU+j]=0;
          //}
        }
      }
  
      itt++;
      MPI_Barrier(MPI_COMM_WORLD);
    }
  
    if (CUTRG>1) {
      minimize_optimize(x3,xoff,gain);
    }
  
    MPI_Barrier(MPI_COMM_WORLD);
  
    nmatres=newnr[nbolo]+nbolo*(GAINSTEP+npixbeam)+nmatdust+nmatco+nfreefree;
    if (rank==0) {
  
      PIOSTRING commm;
      PIOSTRING saveg;
      // Write offset
      PIODOUBLE *tmpoff = (PIODOUBLE *) malloc(sizeof(PIODOUBLE)*(globalRangeRing));
      
      for (i=0;i<nbolo;i++) {
        sprintf(commm,"begin=%lld;end=%lld",(long long) globalBeginRing,(long long) globalEndRing);
        for (j=0;j<globalRangeRing;j++) {
          tmpoff[j]=-10000;
        }
        for (j=newnr[i];j<newnr[i+1];j++) tmpoff[rgordinv[i][j-newnr[i]]]=x3[j];
        if (stim_first_seed+iter==0) sprintf(saveg,"%s_OFF",Param->Out_Offset[i]);
        else sprintf(saveg,"%s_%d_OFF",Param->Out_Offset[i],stim_first_seed+iter);
  
        fprintf(stderr,"Write OFF  %lld\n",(long long) PIOWriteVECT(saveg,tmpoff,sizeof(PIODOUBLE)*globalBeginRing,sizeof(PIODOUBLE)*(globalRangeRing)));
        
        for (j=0;j<globalRangeRing;j++) tmpoff[j]=gain[invgi[j+globalBeginRing+i*32000]+GAINSTEP*i];
        if (stim_first_seed+iter==0) sprintf(saveg,"%s_GAIN",Param->Out_Offset[i]);
        else sprintf(saveg,"%s_%d_GAIN",Param->Out_Offset[i],stim_first_seed+iter);
  
        fprintf(stderr,"Write GAIN  %lld\n",(long long) PIOWriteVECT(saveg,tmpoff,sizeof(PIODOUBLE)*globalBeginRing,sizeof(PIODOUBLE)*(globalRangeRing)));
  
        if (stim_first_seed+iter==0) sprintf(saveg,"%s_GADU",Param->Out_Offset[i]);
        else sprintf(saveg,"%s_%d_GADU",Param->Out_Offset[i],stim_first_seed+iter);
        fprintf(stderr,"Write GADU  %lld\n",(long long) PIOWriteVECT(saveg,gainadu+i*GAINSTEPADU,0,sizeof(PIODOUBLE)*(GAINSTEPADU)));
  
        if (CUTRG>1) {
          PIODOUBLE *tmpcutoff = (PIODOUBLE *) malloc(sizeof(PIODOUBLE)*(globalRangeRing)*CUTRG);
          for (j=0;j<CUTRG*(globalRangeRing);j++) {
            tmpcutoff[j]=-10000;
          }
          if (stim_first_seed+iter==0) sprintf(saveg,"%s_X3",Param->Out_Offset[i]);
          else sprintf(saveg,"%s_%d_X3",Param->Out_Offset[i],stim_first_seed+iter);
  
          for (j=newnr2[i];j<newnr2[i+1];j++) tmpcutoff[rgordinv2[i][j-newnr2[i]]]=xoff[j];
          fprintf(stderr,"Save X3 CUTRG\n");
          sprintf(commm,"begin=%lld;end=%lld",(long long) globalBeginRing*CUTRG,(long long) globalEndRing*CUTRG);
          fprintf(stderr,"Write OPTOFF  %lld\n",(long long) PIOWriteVECT(saveg,tmpcutoff,globalBeginRing*CUTRG*sizeof(PIODOUBLE),
                                                                         (globalRangeRing)*CUTRG*sizeof(PIODOUBLE)));
          free(tmpcutoff);
        }
      }
      free(tmpoff);
  
      if (stim_first_seed+iter==0) sprintf(saveg,"%s_X2",Param->Out_Offset[0]);
      else sprintf(saveg,"%s_%d_X2",Param->Out_Offset[0],stim_first_seed+iter);
  
      //for (i=newnr[nbolo]+GAINSTEP*nbolo;i<nmatres;i++) fprintf(stderr,"X2 %d %lg\n",(int) i,x3[i]);
      fprintf(stderr,"Save X2 NO CUTRG\n");
      fprintf(stderr,"Write MAT  %lld\n",(long long) PIOWriteVECT(saveg,x3+newnr[nbolo],0,(nbolo*(GAINSTEP+npixbeam)+nmatdust+nmatco+nfreefree)*sizeof(PIODOUBLE)));
  
    }
  
    double avvgain=0;
    for (i=0;i<nbolo*GAINSTEP;i++) avvgain+=x3[newnr[nbolo]+i];
    avvgain/=((double)(nbolo*GAINSTEP));
  
    
    if (rank==0)  {
      fprintf(stderr,"AVVGAIN: %lf\n",avvgain);
    }
  
    if (number_of_iterations==iter+1) { 
      free(x2);
      free(x2old);
      free(x2init);
      free(b2);
      free(d2);
      free(q2);
      free(r2);
      free(s2);
      free(hit2); 
      
      free(dthetai );
      free(dthetaq );
      free(dthetau );
      free(dii     );
      free(dqq     );
      free(duu     );
      free(dcoi );
      free(dcoq );
      free(dcou );
      free(dfri );
      free(dfrq );
      free(dfru );
      free(ddusti);
      free(ddustq);
      free(ddustu);
      free(dpixi );
      free(dpixq );
      free(dpixu );
      free(cdip );
      free(cco);
      free(ccfree);
      free(ctheta);
      free(cdust);
      free(cpix );
      free(ctmp);
    }
  
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank%64==0) {
      GetProcMem(&vmem,&phymem);
      fprintf(stderr,"Rank: %ld[%d] MEM %.1lf[%.1lf]MB line=%d\n",
                (long) rank, getpid(),
                (double) vmem/1024./1024.,
                (double) phymem/1024./1024.,__LINE__);
    }
  
  
    /*======================================================
      =
      =    build maps
      =
      =*/
  
  
    double *MAPS_CORR[3]; // corrected I, Q, U maps
    double *MAPS_COV[6];  // corresponding II, IQ, IU, QQ, QU, UU covariance matrix
    
    for (i=0; i<6; i++) {
      MAPS_COV[i] = (double *) malloc(sizeof(double)*nnbpix);
      assert( MAPS_COV[i] != NULL);
    }
    for (i=0; i<3; i++) {
      MAPS_CORR[i] = (double *) malloc(sizeof(double)*nnbpix);
      assert( MAPS_CORR[i] != NULL);
    }
  
    for (int detset=0; detset<Param->n_MAP; detset++) {
      // convert MAPRINGS codes to list of survey maps to produce
      long surveys = 0; // one bit in <surveys> = one map to produce
      PIOSTRING mapname; // map name suffix describing the ring range / survey
      int isurv = 0; // to know if we're in the first loop iteration
      long current_survey = 0; // currently processed survey bit
  
      if (Param->MAPRINGS[detset] & FULL) {
        surveys += MAPFULL;
      }
      if (Param->MAPRINGS[detset] & HM12) {
        surveys += HM1 + HM2;
      }
      if (Param->MAPRINGS[detset] & S12345) {
        surveys += S1 + S2 + S3 + S4 + S5;
      }
      if (Param->MAPRINGS[detset] & YEAR12) {
        surveys += YEAR1 + YEAR2;
      }
      if (Param->MAPRINGS[detset] & FULLODDEVEN) {
        surveys += FULLODD + FULLEVEN;
      }
      if (Param->MAPRINGS[detset] & HM12ODDEVEN) {
        surveys += HM1ODD + HM1EVEN + HM2ODD + HM2EVEN;
      }
  
      while (surveys != 0) {
        int nside_out = Param->Nside;
        int reduced_nside = nside_out;
        if (nside_out == 2048) {
          // reduced_nside is used to reduce nside for maps with bad conditioning
          reduced_nside = 1024;
        }
        for (i=0; i<3; i++) {
          memset(MAPS_CORR[i],0,sizeof(double)*nnbpix);
        }
        for (i=0; i<6; i++) {
          memset(MAPS_COV[i],0,sizeof(double)*nnbpix);
        }
        GetProcMem(&vmem,&phymem);
        if (rank%64==0) fprintf(stderr,"Rank: %ld[%d] MEM %.1lf[%.1lf]MB line=%d\n",
                                (long) rank, getpid(),
                                (double) vmem/1024./1024.,
                                (double) phymem/1024./1024.,__LINE__);
        
        MPI_Barrier(MPI_COMM_WORLD);
  
        if (surveys & MAPFULL) {
          current_survey = MAPFULL;
          strcpy( mapname, "full");
          
        } else if (surveys & HM1) {
          current_survey = HM1;
          strcpy( mapname, "hm1");
          if (freq == 100) nside_out = reduced_nside;
          
        } else if (surveys & HM2) {
          current_survey = HM2;
          strcpy( mapname, "hm2");
          if (freq == 100) nside_out = reduced_nside;
          
        } else if (surveys & S1) {
          current_survey = S1;
          strcpy( mapname, "s1");
          
        } else if (surveys & S2) {
          current_survey = S2;
          strcpy( mapname, "s2");
          
        } else if (surveys & S3) {
          current_survey = S3;
          strcpy( mapname, "s3");
          
        } else if (surveys & S4) {
          current_survey = S4;
          strcpy( mapname, "s4");
          
        } else if (surveys & S5) {
          current_survey = S5;
          strcpy( mapname, "s5");
          
        } else if (surveys & YEAR1) {
          current_survey = YEAR1;
          strcpy( mapname, "year1");
          
        } else if (surveys & YEAR2) {
          current_survey = YEAR2;
          strcpy( mapname, "year2");
          
        } else if (surveys & FULLODD) {
          current_survey = FULLODD;
          strcpy( mapname, "fullodd");
          
        } else if (surveys & FULLEVEN) {
          current_survey = FULLEVEN;
          strcpy( mapname, "fulleven");
  
        } else if (surveys & HM1ODD) {
          current_survey = HM1ODD;
          strcpy( mapname, "hm1odd");
          nside_out = reduced_nside;
          
        } else if (surveys & HM1EVEN) {
          current_survey = HM1EVEN;
          strcpy( mapname, "hm1even");
          nside_out = reduced_nside;
  
        } else if (surveys & HM2ODD) {
          current_survey = HM2ODD;
          strcpy( mapname, "hm2odd");
          nside_out = reduced_nside;
          
        } else if (surveys & HM2EVEN) {
          current_survey = HM2EVEN;
          strcpy( mapname, "hm2even");
          nside_out = reduced_nside;
        }
        surveys &= ~current_survey; // clear current_survey bit
  
        for (k=0; k<nnbpix; k++) { // for each pixel in the rank
          long ndata = loc_nhpix[k];
          long m;
          hpix *htmp = loc_hpix[k];
          long ri,ri2;
          
          for (i=0;i<ndata;i++) {
            ri=htmp[i].rg-globalBeginRing;
            ri2=htmp[i].hrg-globalBeginRing*CUTRG;
            ib=htmp[i].ib;
            
            // select only bolometers in detset
            int use_bolo = Param->bolomask[detset*nbolo+ib];
  
            if (use_bolo) {
              // and select only needed rings
              switch (current_survey) {
                case HM1:
                  if (htmp[i].surv > 10) use_bolo = 0;
                  break;
                case HM2:
                  if (htmp[i].surv < 10) use_bolo = 0;
                  break;
                case S1:
                  if (htmp[i].surv%10 != 1) use_bolo = 0;
                  break;
                case S2:
                  if (htmp[i].surv%10 != 2) use_bolo = 0;
                  break;
                case S3:
                  if (htmp[i].surv%10 != 3) use_bolo = 0;
                  break;
                case S4:
                  if (htmp[i].surv%10 != 4) use_bolo = 0;
                  break;
                case S5:
                  if (htmp[i].surv%10 != 5) use_bolo = 0;
                  break;
                case YEAR1:
                  if (htmp[i].surv%10!=1 && htmp[i].surv%10!=2) use_bolo=0;
                  break;
                case YEAR2:
                  if (htmp[i].surv%10!=3 && htmp[i].surv%10!=4) use_bolo=0;
                  break;
                case FULLODD:
                  if (htmp[i].rg%2 == 0) use_bolo = 0;
                  break;
                case FULLEVEN:
                  if (htmp[i].rg%2 == 1) use_bolo = 0;
                  break;
                case HM1ODD:
                  if (htmp[i].surv > 10) use_bolo = 0;
                  if (htmp[i].rg%2 == 0) use_bolo = 0;
                  break;
                case HM1EVEN:
                  if (htmp[i].surv > 10) use_bolo = 0;
                  if (htmp[i].rg%2 == 1) use_bolo = 0;
                  break;
                case HM2ODD:
                  if (htmp[i].surv < 10) use_bolo = 0;
                  if (htmp[i].rg%2 == 0) use_bolo = 0;
                  break;
                case HM2EVEN:
                  if (htmp[i].surv < 10) use_bolo = 0;
                  if (htmp[i].rg%2 == 1) use_bolo = 0;
              }
            }
  
            if (flg_rg[ib][ri]!=0 && use_bolo==1) {
              ri=rgord[ib][ri]+newnr[htmp[i].ib];
              if (CUTRG>1) {
                ri2=rgord2[ib][ri2]+newnr2[htmp[i].ib];
              }
              double gg2 = gain[htmp[i].gi+htmp[i].ib*GAINSTEP];
              double sig_corr = htmp[i].sig*gg2-htmp[i].fsl-gainadu[htmp[i].adu+htmp[i].ib*GAINSTEPADU];
              if (REMHDIP==0) {
                sig_corr -= htmp[i].dip;
              }
              else {
                sig_corr -= htmp[i].freefree;
              }
              
              if (CUTRG>1) {
                sig_corr-=xoff[ri2];
              }
              else {
                sig_corr-=x3[ri];
              }
              
              sig_corr-=x3[newnr[nbolo]+htmp[i].gi+htmp[i].ib*GAINSTEP]*htmp[i].dip;        
              
              if (nmatco!=0) {
                sig_corr -=htmp[i].comap*x3[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP)+htmp[i].ib];
              }
              if (nmatdust!=0) {
                sig_corr -=htmp[i].dustmap*x3[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP)+nmatco+htmp[i].ib];
              }
              if (nfreefree!=0) {
                sig_corr -=htmp[i].freefree*x3[newnr[nbolo]+nbolo*(npixbeam+GAINSTEP)+nmatco+nmatdust+htmp[i].ib];
              }
              
              for (m=0;m<npixbeam;m++)  {
                sig_corr -= htmp[i].listofpix[m]*x3[newnr[nbolo]+nbolo*(GAINSTEP)+m*nbolo+htmp[i].ib];
              }
              
              double CO1 = eta[htmp[i].ib]*(dpsico[htmp[i].ib]*htmp[i].co
                           - dpsisi[htmp[i].ib]*htmp[i].si);
              double SI1 = eta[htmp[i].ib]*(dpsico[htmp[i].ib]*htmp[i].si
                           + dpsisi[htmp[i].ib]*htmp[i].co);
                          
              MAPS_CORR[0][k] += htmp[i].w*sig_corr;
              MAPS_COV [0][k] += htmp[i].w;
              if (Param->OUT_NOPOL[detset]%2==1 && htmp[i].vi!=UNSEENPIX) {
                MAPS_CORR[1][k] += htmp[i].w*(sig_corr-htmp[i].vi);
                MAPS_COV [1][k] += htmp[i].w;
              }
              else {
                MAPS_CORR[1][k] += htmp[i].w*sig_corr*CO1;
                MAPS_COV [1][k] += htmp[i].w*CO1;
              }
              MAPS_CORR[2][k] += htmp[i].w*sig_corr*SI1;
              MAPS_COV [2][k] += htmp[i].w*SI1;
              MAPS_COV [3][k] += htmp[i].w*CO1*CO1;
              MAPS_COV [4][k] += htmp[i].w*CO1*SI1;
              MAPS_COV [5][k] += htmp[i].w*SI1*SI1;
            }
          }  
        }
        
        MPI_Barrier(MPI_COMM_WORLD);
        double l_th, l_ph; // local_theta, local_phi
        int p; // some map counter
  
  #if 0
        /*==========================================================
        SAUVEGARDE MAP_CORR AVANT INVERSION POUR REPROJECTION APRES CALCUL
        ===========================================================*/
        for (p=0;p<3;p++) {
          char suffix[10];
          if (p==0) strcpy(suffix,"IS");
          if (Param->OUT_NOPOL[detset]%2!=1 || p==0) {
            if (p==1) strcpy(suffix,"QS");
            if (p==2) strcpy(suffix,"US");
            
            PIOSTRING OUTMAP;
                  if (Param->flag_stim_paramfiles == 1) {
                    sprintf(OUTMAP,"%s_%03d_%s_%s",mapout[detset], (int) iter + stim_first_seed, mapname, suffix);
                  } else {
                    sprintf(OUTMAP,"%s_%s_%s", mapout[detset], mapname, suffix);
                  }     
            PIOWriteMAP(OUTMAP,MAPS_CORR[p],begpix[rank],begpix[rank]+nnbpix-1);
            MPI_Barrier(MPI_COMM_WORLD);
          }
        }
  #endif
  
        for (k=0;k<nnbpix;k++) {
          long ndata = loc_nhpix[k];
          hpix *htmp = loc_hpix[k];
          if (detset==0 && isurv==0) {
            for (i=0;i<ndata;i++) htmp[i].vi=UNSEENPIX;
          }
          
          pix2ang_ring(Nside,k+begpix[rank],&l_th,&l_ph);
  
          if (MAPS_COV[0][k]>0) {
            if (Param->OUT_NOPOL[detset]%2==1) {
              MAPS_CORR[0][k]/=MAPS_COV[0][k];
              if (MAPS_COV[1][k]>0) {
                MAPS_CORR[1][k]/=MAPS_COV[1][k]; // get polar cleanned version
              }
              else {
                MAPS_CORR[1][k]=UNSEENPIX;
              }
              MAPS_CORR[2][k]=UNSEENPIX;
            }
            else {
              solvemap(MAPS_CORR[0]+k,MAPS_CORR[1]+k,MAPS_CORR[2]+k,
                       MAPS_COV[0][k],MAPS_COV[1][k],MAPS_COV[2][k],
                       MAPS_COV[3][k],MAPS_COV[4][k],MAPS_COV[5][k]);
              double cond=cond_3_3_thres(MAPS_COV[0][k],MAPS_COV[1][k],MAPS_COV[2][k],
                                         MAPS_COV[1][k],MAPS_COV[3][k],MAPS_COV[4][k],
                                         MAPS_COV[2][k],MAPS_COV[4][k],MAPS_COV[5][k]);
              if (cond >= Param->seuilcond) {
                for (p=0;p<6;p++) {
                  MAPS_COV[p][k]=0;
                }
                for (p=0;p<3;p++) {
                  MAPS_CORR[p][k]=UNSEENPIX;
                }
              }
              else {
                if (detset==0 && isurv==0) {
                  for (i=0;i<ndata;i++) {
                    double CO1 = eta[htmp[i].ib]*(dpsico[htmp[i].ib]*htmp[i].co
                                 - dpsisi[htmp[i].ib]*htmp[i].si);
                    double SI1 = eta[htmp[i].ib]*(dpsico[htmp[i].ib]*htmp[i].si
                                 + dpsisi[htmp[i].ib]*htmp[i].co);
              
                    htmp[i].vi=CO1*MAPS_CORR[1][k]+SI1*MAPS_CORR[2][k];
                  }
                }
              }
            }
          }
          else {
            for (p=0;p<6;p++) {
              MAPS_COV[p][k]=0;
            }
            for (p=0;p<3;p++) {
              MAPS_CORR[p][k]=UNSEENPIX;
            }
          }
        }
        if (rank==0) fprintf(stderr,"%d %d : %ld %ld %ld\n",__LINE__,rank,(long) begpix[rank],(long) edpix[rank],(long) nnbpix);
        
        // only save II,QQ,UU,etc. when not running stim or in first realisation of stim
        if (stim_first_seed + iter == 0) {
          for (p=0;p<6;p++) {
            char suffix[10];
            if (p==0) strcpy(suffix,"II");
            if (Param->OUT_NOPOL[detset]%2!=1 || p==0) {
              if (p==1) strcpy(suffix,"IQ");
              if (p==2) strcpy(suffix,"IU");
              if (p==3) strcpy(suffix,"QQ");
              if (p==4) strcpy(suffix,"QU");
              if (p==5) strcpy(suffix,"UU");
              
              PIOSTRING OUTMAP;
              sprintf(OUTMAP,"%s_%s_%s", mapout[detset], mapname, suffix);
              PIOWriteMAP(OUTMAP,MAPS_COV[p],begpix[rank],begpix[rank]+nnbpix-1);
              MPI_Barrier(MPI_COMM_WORLD);
            }
          }
        }
        
        for (p=0;p<3;p++) {
          char suffix[10];
          if (Param->OUT_NOPOL[detset]%2==0 || p==0 || (p==1 && Param->OUT_NOPOL[detset]%2==1 && Param->OUT_NOPOL[0]%2==0)) {
            if (p==0) strcpy(suffix,"I");
            if (p==1) strcpy(suffix,"Q");
            if (p==1 && Param->OUT_NOPOL[detset]%2==1 && Param->OUT_NOPOL[0]%2==0) strcpy(suffix,"C");
            if (p==2) strcpy(suffix,"U");
            
            PIOSTRING OUTMAP;
            if (Param->flag_stim_paramfiles == 1) {
              sprintf(OUTMAP,"%s_%03d_%s_%s",mapout[detset], (int) iter + stim_first_seed, mapname, suffix);
            } else {
              sprintf(OUTMAP,"%s_%s_%s", mapout[detset], mapname, suffix);
            }     
            PIOWriteMAP(OUTMAP,MAPS_CORR[p],begpix[rank],begpix[rank]+nnbpix-1);
            MPI_Barrier(MPI_COMM_WORLD);
          }
        }
        
        if (rank==0) fprintf(stderr,"line: %d rank: %d file: %s\n",__LINE__,rank,__FILE__);
        MPI_Barrier(MPI_COMM_WORLD);
        isurv++;
      }
    }
  
    for (i=0;i<6;i++) {
      free(MAPS_COV[i]);
    }
    for (i=0;i<3;i++) {
      free(MAPS_CORR[i]);
    }
  } // while (iter++ < number_of_iterations)

  MPI_Barrier(MPI_COMM_WORLD);

  if (rank==0) {
    now = time( NULL);
    fprintf(stderr, "%s: --------------------------\n", __FILE__ );
    fprintf(stderr, "%s: Finished sucessfully at %s",   __FILE__, ctime( &now));
    fprintf(stderr, "%s: --------------------------\n", __FILE__ );
  }
  
  MPI_Finalize();        /* free parameters info */

  exit (0);
}

