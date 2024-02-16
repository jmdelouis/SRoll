// this is for srand48(), drand48(), realpath()
#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE 700
#endif

// MAP NAME DEFINITION
#define MAX_OUT_NAME_LENGTH (2048)

#define CNN_NSIDE (32)


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
#include <assert.h>
#include <errno.h>
#include <omp.h>
#include <libgen.h>
#include "spline.h"

#include "no_dmc_piolib_type_def.h"
#include "no_dmc_metadata.h"
#include "no_dmc_data_access.h"
#include "no_dmc_util.h"
#include "no_dmc_debug.h"
#include "no_dmc_version.h"
#include "chealpix.h"

#define PY_SSIZE_T_CLEAN
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#if 0
#include <numpy/arrayobject.h>
#endif

#include "troll_param.h"

#include "stim_parLoader.h"
#include "stim.h"
#include "stim_tools.h"
#if 1
#define OPTIMPI
#endif

#ifndef DOOFFSET
#define DOOFFSET 1
#endif

PIOINT **rgordinv;

float *tpparam=NULL; // parameters for neural network initialise to NULL at the beginning.
int CNN_NB_PARAM=0;

int verbose=0;
PIOINT Nside;
long NORM_GAIN=0;
long REMOVE_CAL=0;

int NUMBEROFITER=500;
int do_offset=DOOFFSET;
int NORMFITPOL=0;

int PIOWriteVECT(const char *path,void *value,int off,int size);
double det(double *mat, int n);

int *realpix;
int *irealpix;

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

typedef struct {
  int ipx;
  int hit;
} hpint;

typedef struct {
  PyObject * sys;
  PyObject * py_path;
  PyObject * ModuleString;
  PyObject * Module;
  PyObject * Dict;
  PyObject * init;
  PyObject * initFromFile;
  PyObject * initFrom_Files;
  PyObject * init_net;
  PyObject * init_net_data;
  PyObject * allocf32;
  PyObject * alloci32;
  PyObject * Clean;
  PyObject * grad;
  PyObject * agrad;
  PyObject * gloss;
  PyObject * close;
  PyObject * pred;
  PyObject * corr;
  PyObject * initpar;
  PyObject * getparam;
  PyObject * initconvw;
  PyObject * getconvw;
  PyObject * initconvb;
  PyObject * getconvb;
  PyObject * initflw;
  PyObject * getflw;
  PyObject * initflb;
  PyObject * getflb;
  PyObject * pltvec;
  PyObject * plthisto;
  PyObject * run_foscat;
} WrapPython;

int python_rank;
int mpi_python_size;
int tensorflow_rank;
int mpi_tensorflow_size;
MPI_Comm python_comm;
MPI_Comm tensorflow_comm;
int NB_EXTERNAL=0;

int compar_int(const void *a, const void *b)
{
  hpint *pa = (hpint *) a;
  hpint *pb = (hpint *) b;
  return(pb->hit-pa->hit);
}

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
    +((idx/4194304)%2)*Nside;
  return(res);
}

int gety12(int idx)
{
  int res=((idx/2)%2)
    +((idx/8)%2)*2
    +((idx/32)%2)*4
    +((idx/128)%2)*8
    +((idx/512)%2)*16
    +((idx/Nside)%2)*32
    +((idx/8192)%2)*64
    +((idx/32768)%2)*128
    +((idx/131072)%2)*256
    +((idx/524288)%2)*512
    +((idx/2097152)%2)*1024
    +((idx/8388608)%2)*Nside;
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
    +((x/Nside)%2)*4194304
    +((y)%2)*2
    +((y/2)%2)*8
    +((y/4)%2)*32
    +((y/8)%2)*128
    +((y/16)%2)*512
    +((y/32)%2)*Nside
    +((y/64)%2)*8192
    +((y/128)%2)*32768
    +((y/256)%2)*131072
    +((y/512)%2)*524288
    +((y/1024)%2)*2097152
    +((y/Nside)%2)*8388608;
  return(res);
}

#define MAXRAND 1000000
// ------------------------------------------------------------------------
PyObject *EXECPYTHON(PyObject *TheObject)
{
  if (TheObject==NULL) {
    fprintf(stderr,"EXECPYTHON : error object null\n");
    PyErr_Print();
    MPI_Finalize(); 
    exit(0);
  }
  return(TheObject);
}
PyObject *CALLPYTHON(PyObject *Dict,const char *name)
{
  PyObject *res=PyDict_GetItemString(Dict,name);
  if (PyCallable_Check(res)) {
    return(res); 
  }
  else {
    fprintf(stderr,"Problem while loading %s function\n",name);
    PyErr_Print();
    MPI_Finalize(); 
    exit(0);
  }
  return(NULL);
}

// ------------------------------------------------------------------------
typedef struct {
  int degree;
  int nodes;
  double *norm;
} spl1d;

// ------------------------------------------------------------------------
typedef struct {
   int degree;
   int nodes;
   double *norm;
   spl1d *spline1;
   spl1d *spline2;
} scirc;


scirc *myspline;  //define splines

// ------------------------------------------------------------------------
int Fact_spline(long int x)
{
   if (x <= 1)
     return 1;
   return (x * Fact_spline(x-1));
}
// ------------------------------------------------------------------------
//degree 3 par default
scirc *circ_spline(const int nodes,const int degree)
{
   int i;
   scirc *out = (scirc *) malloc(sizeof(scirc));
   out->degree=degree;
   out->nodes=nodes;
   out->norm=(double *) malloc(sizeof(double)*(out->degree+1));

   for (i=0;i<out->degree+1;i++) {
     out->norm[i]=pow(-1,i)*(out->degree+1)/
       (Fact_spline(out->degree+1-i)*Fact_spline(i));
   }
   return(out);
}
// ------------------------------------------------------------------------
spl1d *spline1d(int nodes,int degree)
{
  int i;
  spl1d *out = (spl1d *) malloc(sizeof(spl1d));
  out->degree=degree;
  out->nodes=nodes;
  out->norm=(double *) malloc(sizeof(double)*(out->degree+1));
  
  for (i=0;i<out->degree+1;i++) {
    out->norm[i]=pow(-1,i)*(out->degree+1)/
      (Fact_spline(out->degree+1-i)*Fact_spline(i));
  }
  return(out);
}
// ------------------------------------------------------------------------
void free_spline1d(spl1d *spline)
{
  free(spline);
}

// ------------------------------------------------------------------------
//degree 3 par default
scirc *circ_spline2D(const int nodes1,const int nodes2,const int degree)
{
   scirc *out = (scirc *) malloc(sizeof(scirc));
   out->spline1=(spl1d *) circ_spline(nodes1,degree);
   out->spline2=spline1d(nodes2,degree);
   return(out);
}
// ------------------------------------------------------------------------
void free_circ_spline(scirc *spline)
{
   free(spline);
}
// ------------------------------------------------------------------------
double yplus_spline1d(spl1d *spline,double x)
{
  if (x<0.0) return(0.0);
  if (spline->degree==0) {
    if (x==0.0) return(0.5);
    else return(1.0);
  }
  return(pow(x,spline->degree));
}
// ------------------------------------------------------------------------
void calc_spline1d(spl1d *spline,double x,float *y)
{
  int i,j;

  memset(y,0,spline->nodes*sizeof(double));
  for (i=0;i<spline->nodes;i++) {
    double tmp=0;
    double tx=(spline->nodes-1)*x-i;
    if (x<0) tx=-i;
    if (x>1.0) tx=(spline->nodes-1)-i;
    for (j=0;j<spline->degree+1;j++) {
      tmp+=spline->norm[j]*yplus_spline1d(spline,tx-j+(spline->degree+1)/2);
    }
    if (tmp<0) tmp=0.0;
    y[i]+=tmp;
  }
  double tmp=0;
  for (i=0;i<spline->nodes;i++) {
    tmp+=y[i];
  }
  for (i=0;i<spline->nodes;i++) {
    y[i]/=tmp;
  }
}

// ------------------------------------------------------------------------
double yplus_circ_spline(scirc *spline,double x)
{
   if (x<0.0) return(0.0);
   if (spline->degree==0) {
     if (x==0.0) return(0.5);
     else return(1.0);
   }
   return(pow(x,spline->degree));
}
// ------------------------------------------------------------------------
void calc_circ_spline(scirc *spline,double x,float *y)
{
   int i,j;
 
   memset(y,0,spline->nodes*sizeof(double));
   for (i=0;i<spline->nodes+spline->degree/2+1;i++) {
     double tmp=0;
     double tx=spline->nodes*fmod(x+2*M_PI,2*M_PI)/(M_PI*2)-i;
     for (j=0;j<spline->degree+1;j++) {
      tmp+=spline->norm[j]*yplus_circ_spline(spline,tx-j+(spline->degree+1)/2);
     }
     if (tmp<0) tmp=0.0;
     y[i%spline->nodes]+=tmp;
   }
   for (i=0;i<spline->degree/2;i++) {
     double tmp=0;
     double tx=spline->nodes*fmod(x+2*M_PI,2*M_PI)/(M_PI*2)+1+i;
     for (j=0;j<spline->degree+1;j++) {
      tmp+=spline->norm[j]*yplus_circ_spline(spline,tx-j+(spline->degree+1)/2);
     }
     if (tmp<0) tmp=0.0;
     y[spline->nodes-1-i]+=tmp;
   }

}


// ------------------------------------------------------------------------
void calc_spline(scirc *spline,double x,float *y)
{
   int i,j;
 
   memset(y,0,spline->nodes*sizeof(double));
   for (i=0;i<spline->nodes;i++) {
     double tmp=0;
     double tx=spline->nodes*x-i;
     for (j=0;j<spline->degree;j++) {
      tmp+=spline->norm[j]*yplus_circ_spline(spline,tx-j+(spline->degree+1)/2);
     }
     if (tmp<0) tmp=0.0;
     y[i]+=tmp;
   }
}

// ------------------------------------------------------------------------
void calc_circ_spline2D(scirc *spline,double x,double y,float *v)
{
  float *v1 = (float *) malloc(spline->spline1->nodes*sizeof(float));
  float *v2 = (float *) malloc(spline->spline2->nodes*sizeof(float));
  
  calc_circ_spline((scirc *)spline->spline1,x,v1);
  calc_spline1d(spline->spline2,y,v2);

  for (int i=0;i<spline->spline1->nodes;i++) {
    for (int j=0;j<spline->spline2->nodes;j++) {
      v[i+spline->spline1->nodes*j]=v1[i]*v2[j];
    }
  }
  
  free(v1);
  free(v2);
}

// --------------------------------------------------------------------------------
void InitPython(WrapPython *mywrap, char *myfunct, int rank)
{
  char thepath[1024];
  int i,lastslash=-1;
  strcpy(thepath,myfunct);
  for (i=0;i<(int)strlen(myfunct);i++) 
    if (myfunct[i]=='/') lastslash=i;
  if (lastslash==-1) {
    sprintf(thepath,".");
    lastslash=0;
  }
  else {
    thepath[lastslash]='\0';
    lastslash++;
  }
  if (rank==0) {
    fprintf(stderr,"Call python module in path\n%s\n",thepath);
  }
  Py_Initialize();

  PyRun_SimpleString("sys.path.append(os.getcwd())");

  mywrap->sys = EXECPYTHON(PyImport_ImportModule("sys"));
  
  mywrap->py_path = EXECPYTHON(PyObject_GetAttrString(mywrap->sys, "path"));
  
#ifdef PYTHON3
  PyList_Append(mywrap->py_path, PyUnicode_FromString(thepath));
#else
  PyList_Append(mywrap->py_path, PyString_FromString(thepath));
#endif

  strcpy(thepath,myfunct+lastslash);
  lastslash=-1;
  for (i=0;i<(int)strlen(thepath);i++) 
    if (thepath[i]=='.') lastslash=i;
  if (lastslash!=-1) thepath[lastslash]='\0';

  if (rank==0) {
    fprintf(stderr,"Call module : %s\n",thepath);
  }

  
#ifdef PYTHON3
  mywrap->ModuleString = PyUnicode_FromString(thepath);
#else
  mywrap->ModuleString = PyString_FromString(thepath);
#endif

    
  mywrap->Module = EXECPYTHON(PyImport_Import(mywrap->ModuleString));
  
  mywrap->Dict = EXECPYTHON(PyModule_GetDict(mywrap->Module));
  mywrap->run_foscat = CALLPYTHON(mywrap->Dict, "run");

#if 0
#ifdef PYTHON3
  mywrap->ModuleString = PyUnicode_FromString(thepath);
#else
  mywrap->ModuleString = PyString_FromString(thepath);
#endif
  mywrap->Module = EXECPYTHON(PyImport_Import(mywrap->ModuleString));
  mywrap->Dict = EXECPYTHON(PyModule_GetDict(mywrap->Module));

  
  
  mywrap->init     = CALLPYTHON(mywrap->Dict, "init_shape");
  mywrap->initFrom_Files = CALLPYTHON(mywrap->Dict, "init_shape_from_files");
  mywrap->init_net = CALLPYTHON(mywrap->Dict, "init_network");
  mywrap->init_net_data = CALLPYTHON(mywrap->Dict, "init_network_data");
  mywrap->allocf32 = CALLPYTHON(mywrap->Dict, "alloc_table_float32");
  mywrap->alloci32 = CALLPYTHON(mywrap->Dict, "alloc_table_int32");
  mywrap->Clean    = CALLPYTHON(mywrap->Dict, "free_table");
  mywrap->grad     = CALLPYTHON(mywrap->Dict, "calc_grad");
  mywrap->agrad    = CALLPYTHON(mywrap->Dict, "apply_grad");
  mywrap->gloss    = CALLPYTHON(mywrap->Dict, "get_loss");
  mywrap->close    = CALLPYTHON(mywrap->Dict, "close_session");
  mywrap->pred     = CALLPYTHON(mywrap->Dict, "get_prediction");
  mywrap->corr     = CALLPYTHON(mywrap->Dict, "get_correction");
  mywrap->initpar  = CALLPYTHON(mywrap->Dict, "alloc_param");
  mywrap->getparam = CALLPYTHON(mywrap->Dict, "get_param");
  mywrap->initconvw = CALLPYTHON(mywrap->Dict, "alloc_convw");
  mywrap->getconvw  = CALLPYTHON(mywrap->Dict, "get_convw");
  mywrap->initconvb = CALLPYTHON(mywrap->Dict, "alloc_convb");
  mywrap->getconvb  = CALLPYTHON(mywrap->Dict, "get_convb");
  mywrap->initflw = CALLPYTHON(mywrap->Dict, "alloc_flw");
  mywrap->getflw  = CALLPYTHON(mywrap->Dict, "get_flw");
  mywrap->initflb = CALLPYTHON(mywrap->Dict, "alloc_flb");
  mywrap->getflb  = CALLPYTHON(mywrap->Dict, "get_flb");
  mywrap->pltvec    = CALLPYTHON(mywrap->Dict, "plt_vec");
  mywrap->plthisto  = CALLPYTHON(mywrap->Dict, "plt_histo");
#endif
}
// --------------------------------------------------------------------------------
void CleanPython(WrapPython *mywrap)
{
  // Clean up
  Py_DECREF(mywrap->Module);
  Py_DECREF(mywrap->ModuleString);
  
  
  Py_Finalize();
}
// --------------------------------------------------------------------------------
ssize_t pwrite(int fildes, const void *buf, size_t nbyte,
              off_t offset);
// --------------------------------------------------------------------------------
int isPowerOfTwo (unsigned int x)
{
  return ((x != 0) && ((x & (~x + 1)) == x));
}
// --------------------------------------------------------------------------------

#if !defined(RINGSIZE)
#define RINGSIZE (27664l)
#endif

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
    int err=fscanf(fp,"%d %s %c %d %d %d %d %d %lu %lu %lu %lu %lu %lu %lu %ld %ld %ld %ld %ld %ld %lu %lu %ld",
		   &pid,comm,&state,&ppid,&pgrp,&session,&tty_nr,&tpgid,&flags,&minflt,&cminflt,&majflt,
		   &cmajflt,&utime,&stime,&cutime,&cstime,&priority,&nice,&xxx,&itrealvalue,&starttime,&vsize,
		   &rss);
    fclose(fp);
    if (err==EOF) {
      fprintf(stderr,"Problem while trying to get proc info\n");
    }
    *vmem=vsize;
    *phymem=rss*4096;
}
// --------------------------------------------------------------------------------
float GetLoadAvg()
{
    char Line[256];
    float a,b,c;

    FILE *fp=fopen("/proc/loadavg","r");
    char *err=fgets(Line,256,fp);
    fclose(fp);
    if (err==NULL) {
      fprintf(stderr,"Problem while trying to get load average info\n");
    }
    sscanf(Line,"%f %f %f",&a,&b,&c);
    return(a);
}
// --------------------------------------------------------------------------------

#include "troll_parLoader.h"

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

#define UNSEENPIX ((double) (-1.6375000E+30))

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

void invertMatrix(double * mat,double *vec,int n,int rank);

// --------------------------------------------------------------------------------
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
      //fprintf(stderr,"Singular matrix in routine LUDCMP\n");
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
// --------------------------------------------------------------------------------
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
// --------------------------------------------------------------------------------
int lusol(double *a,double *b,int n)
{
  double d;
  int *indx= (int *) malloc(n*sizeof(int));
  if (ludcmp(a,&d,indx,n)!=0) return(-1);
  lubksb(a,b,indx,n);
  free(indx);
  return(0);
}
// --------------------------------------------------------------------------------
void invert(double *a,double *y,int n)
{
  if (n==1) {
    *y=1/(*a);
  }
  double d;
  double *col= (double *) malloc(n*sizeof(double));
  double *l_a= (double *) malloc(n*n*sizeof(double));
  memcpy(l_a,a,n*n*sizeof(double));
  int i,j;
  int *indx= (int *) malloc(n*sizeof(int));
  memset(indx,0,n*sizeof(int));
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
// --------------------------------------------------------------------------------
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
// --------------------------------------------------------------------------------
double Dcond(double *a,int n)
{
  double *b=(double *) malloc(n*n*sizeof(double));
  memset(b,0,n*n*sizeof(double));
  invert(a,b,n);
  double cond=norm(a,n)*norm(b,n);
  free(b);
  return(cond);
}

// --------------------------------------------------------------------------------
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

// --------------------------------------------------------------------------------
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

// ----------------------------------------------------------------------------------------------------- 
double cond_thres(double * mat,double *res,int maxchannels){
 
  double determinant =0.0;
  memset(res,0,sizeof(maxchannels*maxchannels*sizeof(double)));
  determinant = det(mat,maxchannels);    

  if (abs(determinant)<1E-7) return(1E30);
  
  double cond=0.0,cond2=0.0;
  
  if (maxchannels>1) {
    invert(mat,res,maxchannels);
    if (isnan(res[0])) return(1E30);
    cond=norm(mat,maxchannels)*norm(res,maxchannels);
  }else {
    res[0]=1/mat[0];
    cond=abs(res[0]);
    return(cond);
  }
  
  //return(cond);
    
  if (cond<1E10) {
    cond2=norm(res,maxchannels); 
  }
  if (cond<=0.0) return(1E30);
  return(cond2);
}

// ----------------------------------------------------------------------------------------------------- 
void lmatrice(double *mat,double *lmat, int n, int l){
 int ld=0;
 int k=n-1;
 for(int i=0;i<n;i++){
  if(i!=l){
   for(int j=1;j<n;j++){
    lmat[ld+(j-1)*k]=mat[i+j*n];
   }
    ld++;
  }
 }
}
// ----------------------------------------------------------------------------------------------------- 
double det(double *mat, int n){
  
 double resultat=0.0;
 int k=n-1;
 double signe=1.0;
 double lmat[k*k];
 memset(lmat,0,k*k*sizeof(double));

 if(n==1){
  return mat[0];
 }
 for(int i=0;i<n;i++) {
  lmatrice(mat,lmat,n,i);
  resultat=resultat+signe*mat[i]*det(lmat,k);
  signe=-signe;
 }
 return resultat;

}
//-------------------------------------------------------------------------------------------------------
double cond_3_3_thres(double a,double b,double c,
                      double d,double e,double f,
                      double g,double h,double i){
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

//-------------------------------------------------------------------------------------------------------
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

// maximum number of stim iterations in one troll run
#ifndef MAXSIMU
#define MAXSIMU (1)
#endif

// * maximum number of External HPR templates per bolometer
// * from main():
//     npixShpr = Param->n_External / nbolo;
//     assert( npixShpr <= MAXEXTERNALHPR);
// * can be set at compilation time with '-DMAXEXTERNALHPR=xxx'
// * values used in RD12 are 4 at 100ghz-217ghz, 7 at 353ghz and 9 at 545ghz-857ghz
// * need to remember why it was set to 11 for RD12...

// TROLL update:
// MAXEXTERNALHPR = 5 TF (H0, H1, H2, H3, H4 + TT1@353) + ANGLE + POLEFF + TDUST + CO13 + SYNCROTRON = 10
// see "extra external tempaltes" comment around line 6616
// updated assert( npixShpr+5 <= MAXEXTERNALHPR);
// 545/857 : 
// MAXEXTERNALHPR = 8 TF + 1FSL + ANGLE + POLEFF + CO13  = 12

#ifndef MAXEXTERNALHPR
#define MAXEXTERNALHPR (10)
#endif
#ifndef MAXEXTERNALSHPR
#define MAXEXTERNALSHPR (4)
#endif
#ifndef MAXEXTERNALMAP
#define MAXEXTERNALMAP (3)
#endif
#ifndef MAXCHAN
#define MAXCHAN (16)
#endif

typedef struct {
  PIOFLOAT sig;
  PIOFLOAT listp[MAXSIMU];
  PIOINT   nShpr;
  PIOFLOAT listofShpr[MAXEXTERNALSHPR];
  PIOINT   listofShpr_idx[MAXEXTERNALSHPR];
  PIOFLOAT listofmap[MAXEXTERNALMAP];
  PIOINT ipix;
  PIOINT rg;
  PIOFLOAT corr_cnn;
  PIOINT gi;
  PIOFLOAT hpr_cal;
  PIOFLOAT Sub_HPR;
  PIOFLOAT w;
  PIOBYTE  surv;
  PIOBYTE  ib;
  PIOFLOAT hit;
  PIOFLOAT model;
  PIOFLOAT channels[MAXCHAN];
  PIOINT ipix_foscat;
  //add for new version
  int diag_idx;
  double m;
  double R_ij;
  double alpha[MAXCHAN];
} hpix;

int MAXCHANNELS  = 0;
long DODISTOR=0;
#define _PIOMALLOC malloc
#define _PIOFREE free

// ---------------------------------------------------------------------------------------------------
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
          if (rank==0&&iter%10==0) fprintf (stderr,"gcmat_mpi() iter = %d - delta0 = %lg - delta_new = %lg\n", iter, delta0, delta_new);
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

//-------------------------------------------------------------------------------------------------------
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

void matvec(double *mat,double *vec,double *out,int n)
{
  int i,j;
  
  for (i=0;i<n;i++) {
    double tmp = 0;
    for (j=0;j<n;j++) {
        tmp+=mat[j+n*i]*vec[j];
    }
    out[i]=tmp;
  }
}
//-------------------------------------------------------------------------------------------------------
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
//-------------------------------------------------------------------------------------------------------
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
//-------------------------------------------------------------------------------------------------------
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
//-------------------------------------------------------------------------------------------------------
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
//-------------------------------------------------------------------------------------------------------
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

// ------------------------------------------------------------------------
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
// ------------------------------------------------------------------------


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

double **rawcst;

// TODO TODO
long nnbpix;
hpix **loc_hpix;
PIOINT *loc_nhpix;

#if 0
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
#endif
PIOBYTE *flgpix;
PIOLONG globalBeginRing;
PIOLONG globalEndRing;

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

PIOBYTE **flg_rg;
long npixShpr;
long npixmap;
long ittt,itbogo;
long nmatres;
int testzero=-1;

double normaoff=0;
double *NEP_tab = NULL;

long DOCO13=0;
long NDOCNN=0;
int DOCNN[MAXEXTERNALHPR];
PIOBYTE COMP_CNN=0;
WrapPython MyPythonBackend;

PIOINT DOSYNCHRO=0;
double *x2;
double *x2old;
double *x2init;
double *b2;
double *d2;
double *q2;
double *r2;
double *s2;
double *hit2;
double *hit_WW;
double *hit_L1;
double *hit_L2;

long nmatpix;
#define FITGAIN
#ifdef FITGAIN
  int nitbogo=4;
#else
  int nitbogo=1;
#endif
int DOTDUST;
int MAXFREQ=0;
PIODOUBLE *imatrice;
PIOINT **rgord;
troll_parContent *Param;
//int *the_stat_pix;
double delta0;
PIOLONG  globalRangeRing;
double *eta;
double *eta_dest;
PIODOUBLE *cond;
PIODOUBLE *g;

typedef struct {
  PIOLONG *BeginRing;
  PIOLONG *EndRing;
} rankInfo;

rankInfo globalRankInfo;

double **cache_xi2;
double **cache_gain_xi2;
long n_cache_xi2=0;

int NOMOREFITTED=0;

#ifdef UPDATE_DIP
#define UPDATE_DIP
#endif

// -------------------------------------------------------------------------------------------------------------
void proj_data(double *b2,int nnbpix,int rank,double nmatres,int GAINSTEP2){

  long ir,irt;  
  
  for(int pix =0;pix<nnbpix;pix++){
    double sum_Rij = 0.0,tmp = 0.0,sum_channels[MAXCHAN];
    
    long ndata = loc_nhpix[pix];                // nombre de donnees dans le pixel pix
    hpix *htmp = loc_hpix[pix]; 

    memset(sum_channels,0,MAXCHANNELS*sizeof(double)); 

    if (ndata>1&&flgpix[pix]>0) {

      //calcul sum_Rij && calcul sum_channels
      for(int l1=0;l1<ndata;l1++){
        long ri1=htmp[l1].rg-globalBeginRing;
        if (flg_rg[htmp[l1].ib][ri1]!=0) {
          sum_Rij+= htmp[l1].R_ij; 
        }               
        for(int k =0;k<MAXCHANNELS;k++){
          sum_channels[k] += htmp[l1].channels[k]* htmp[l1].R_ij;
        }
      }
 
      // Calcul de b2
      for(int l1=0;l1<ndata;l1++){
        long ri1=htmp[l1].rg-globalBeginRing;

        if (flg_rg[htmp[l1].ib][ri1]!=0) {
          ir=rgord[htmp[l1].ib][ri1]+newnr[htmp[l1].ib];
                 
          tmp = htmp[l1].R_ij;                            
          for(int k =0;k<MAXCHANNELS;k++){
            tmp-= sum_channels[k]*htmp[l1].alpha[k] ;
          }

	  if (do_offset==1) {
	    b2[ir]+=tmp;
	  }
	  
          for(int m=0;m<npixmap;m++){
            b2[newnr[nbolo]+htmp[l1].ib+(m+GAINSTEP2)*nbolo]+=htmp[l1].listofmap[m]*tmp;
          }
          for(int m=0;m<htmp[l1].nShpr;m++){
            b2[newnr[nbolo]+htmp[l1].ib+(htmp[l1].listofShpr_idx[m]+npixmap+GAINSTEP2)*nbolo]+=htmp[l1].listofShpr[m]*tmp;
          }
	  

          if(GAINSTEP2 != 0){
            b2[newnr[nbolo]+htmp[l1].gi+htmp[l1].ib*GAINSTEP2]+= tmp * htmp[l1].model;
          }

        }
      }
      
  
      //Calcul hit2
      // le but est de calculer la variance de chacun des paramètre dans chaque pixel.
      // hit 2 est la somme des variances sur tous les pixels.
      
      for(int l1=0;l1<ndata;l1++){
        long ri1=htmp[l1].rg-globalBeginRing;  
        if (flg_rg[htmp[l1].ib][ri1]!=0) {
	  if (do_offset==1) {
	    ir=rgord[htmp[l1].ib][ri1]+newnr[htmp[l1].ib];      
	    hit_WW[ir]+= htmp[l1].w; // poids stat de chaque valeurs 
	    hit_L2[ir]+= htmp[l1].w; // poids stat de chaque valeurs 
	    //hit_L1[ir]+= htmp[l1].w; // poids stat de chaque valeurs 
	  }
          for(int m = 0;m<npixmap;m++){
            irt = newnr[nbolo]+htmp[l1].ib+nbolo*(m+GAINSTEP2);
            hit_L2[irt]+= htmp[l1].w*htmp[l1].listofmap[m]*htmp[l1].listofmap[m]; // poids stat de chaque valeurs 
            hit_L1[irt]+= htmp[l1].w*htmp[l1].listofmap[m]; // poids stat de chaque valeurs 
            hit_WW[irt]+= htmp[l1].w; // poids stat de chaque valeurs 
          }
          for(int m = 0;m<htmp[l1].nShpr;m++){
            irt = newnr[nbolo]+htmp[l1].ib+nbolo*(htmp[l1].listofShpr_idx[m]+GAINSTEP2+npixmap);
            hit_L2[irt]+= htmp[l1].w*htmp[l1].listofShpr[m]*htmp[l1].listofShpr[m]; // poids stat de chaque valeurs 
            hit_L1[irt]+= htmp[l1].w*htmp[l1].listofShpr[m]; // poids stat de chaque valeurs 
            hit_WW[irt]+= htmp[l1].w; // poids stat de chaque valeurs 
          }
	  
          if(GAINSTEP2 != 0) {
	    irt=newnr[nbolo]+htmp[l1].gi+htmp[l1].ib*GAINSTEP2;
	    hit_L2[irt]+=htmp[l1].w*htmp[l1].model*htmp[l1].model;
	    hit_L1[irt]+=htmp[l1].w*htmp[l1].model;
	    hit_WW[irt]+=htmp[l1].w;
	  }
        }
      }
      //Remise à 0 et calcul de hit2 si hit_WW est different de 0
      //On considere que c'est plus rapide de remetre à 0 comme cela.
      
      for(int l1=0;l1<ndata;l1++){
        long ri1=htmp[l1].rg-globalBeginRing;  
        if (flg_rg[htmp[l1].ib][ri1]!=0) {
	  if (do_offset==1) {
	    ir=rgord[htmp[l1].ib][ri1]+newnr[htmp[l1].ib];
	    if (hit_WW[ir]>0&&hit_L2[ir]!=hit_L1[ir]) {
	      hit2[ir] += hit_L2[ir] - hit_L1[ir]*hit_L1[ir]/hit_WW[ir];
	    }
	    hit_WW[ir] =0; // poids stat de chaque valeurs 
	    hit_L1[ir] =0; // poids stat de chaque valeurs 
	    hit_L2[ir] =0; // poids stat de chaque valeurs
	  }
	  
          for(int m = 0;m<npixmap;m++){
            irt = newnr[nbolo]+htmp[l1].ib+nbolo*(m+GAINSTEP2);
	    if (hit_WW[irt]>0&&hit_L2[irt]!=hit_L1[irt]) {
	      hit2[irt] += (hit_L2[irt] - hit_L1[irt]*hit_L1[irt]/hit_WW[irt]);
	    }
            hit_L2[irt] = 0; // poids stat de chaque valeurs 
            hit_L1[irt] = 0; // poids stat de chaque valeurs 
            hit_WW[irt] = 0; // poids stat de chaque valeurs 
          }
          for(int m = 0;m<htmp[l1].nShpr;m++){
            irt = newnr[nbolo]+htmp[l1].ib+nbolo*(htmp[l1].listofShpr_idx[m]+GAINSTEP2+npixmap);
	    if (hit_WW[irt]>0&&hit_L2[irt]!=hit_L1[irt]) {
	      hit2[irt] += (hit_L2[irt] - hit_L1[irt]*hit_L1[irt]/hit_WW[irt]);
	    }
            hit_L2[irt] = 0; // poids stat de chaque valeurs 
            hit_L1[irt] = 0; // poids stat de chaque valeurs 
            hit_WW[irt] = 0; // poids stat de chaque valeurs 
          }
	  
          if(GAINSTEP2 != 0) {
	    irt=newnr[nbolo]+htmp[l1].gi+htmp[l1].ib*GAINSTEP2;
	    if (hit_WW[irt]>0&&hit_L2[irt]!=hit_L1[irt]) {
	      hit2[irt] += (hit_L2[irt] - hit_L1[irt]*hit_L1[irt]/hit_WW[irt]);
	    }
	    hit_L2[irt] = 0;
	    hit_L1[irt] = 0;
	    hit_WW[irt] = 0;
	  }
        }
      }
    }
  }
}
// -------------------------------------------------------------------------------------------------------------
void proj_grad(double * q2,double nmatres,double * x,int nnbpix,int rank,int GAINSTEP2){
  //plap
  double val =0.0,sum =0.0;
  long ir,irt;
  double tmp= 0.0;
  double *csum = (double *) malloc(sizeof(double)*MAXCHANNELS);

  memset(csum,0,sizeof(double)*MAXCHANNELS);

  for(int pix =0;pix<nnbpix;pix++){

    long ndata = loc_nhpix[pix];                // nombre de donnees dans le pixel k
    hpix *htmp = loc_hpix[pix]; 
    double s_X[MAXCHAN]; // somme des X
    double sum_channels[MAXCHAN];
    double Rij_bis=0.0;
#ifdef CALCMATRIX      
    double MAT[MAXCHAN*MAXCHAN];
    
    //init matrice and s_X
    memset(MAT,0,MAXCHANNELS*MAXCHANNELS*sizeof(double));
#endif
    memset(s_X,0,MAXCHANNELS*sizeof(double)); 
    memset(sum_channels,0,MAXCHANNELS*sizeof(double)); 

    if (ndata>1&&flgpix[pix]>0) {       

#ifdef CALCMATRIX      
      //Calcul matrice 
      for(int l1=0;l1<ndata;l1++){ 
        for(int m=0;m<MAXCHANNELS;m++){
          for(int l =0;l<MAXCHANNELS;l++){
            MAT[m+MAXCHANNELS*l]+= htmp[l1].w * htmp[l1].channels[l]*htmp[l1].channels[m];
          }
        }
       } 
#endif  
 
      //Calcul s_X  && calcul sum_channels
      for(int l1=0;l1<ndata;l1++){ 
   
        long ri1=htmp[l1].rg-globalBeginRing;
        if (flg_rg[htmp[l1].ib][ri1]!=0) {
	  
          ir=rgord[htmp[l1].ib][ri1]+newnr[htmp[l1].ib];
	  if (do_offset==1) {
	    val = x[ir];
          }
	  
          for(int m = 0;m<npixmap;m++){
            irt = newnr[nbolo]+htmp[l1].ib+nbolo*(m+GAINSTEP2);
            val+= x[irt]*htmp[l1].listofmap[m];            
          }
	  
          for(int m = 0;m<htmp[l1].nShpr;m++){
            irt = newnr[nbolo]+htmp[l1].ib+nbolo*(htmp[l1].listofShpr_idx[m]+npixmap+GAINSTEP2);
            val+= x[irt]*htmp[l1].listofShpr[m];            
          }
          
          if(GAINSTEP2 != 0)  val+= x[newnr[nbolo]+htmp[l1].gi+htmp[l1].ib*GAINSTEP2]*htmp[l1].model;

          for(int l =0;l<MAXCHANNELS;l++){ 
	    if (do_offset==1) {
	      csum[l]+=x[ir]*htmp[l1].channels[l];
	    }
            s_X[l]+= htmp[l1].w *htmp[l1].channels[l]*val;
          }  
        }
      }

      invertMatrix(imatrice+pix*MAXCHANNELS*MAXCHANNELS,s_X,MAXCHANNELS,rank);         
  
      //calcul sum_Rij_bis
      for(int l1=0;l1<ndata;l1++){  

        long ri1=htmp[l1].rg-globalBeginRing;
        if (flg_rg[htmp[l1].ib][ri1]!=0) {
          ir=rgord[htmp[l1].ib][ri1]+newnr[htmp[l1].ib];

	  if (do_offset==1) {
	    val = x[ir]; // calcul val
	  }
	  
          for(int m = 0;m<npixmap;m++){
            irt = newnr[nbolo]+htmp[l1].ib+nbolo*(m+GAINSTEP2);
            val+= x[irt]*htmp[l1].listofmap[m];  
          }   
          for(int m = 0;m<htmp[l1].nShpr;m++){
            irt = newnr[nbolo]+htmp[l1].ib+nbolo*(htmp[l1].listofShpr_idx[m]+npixmap+GAINSTEP2);
            val+= x[irt]*htmp[l1].listofShpr[m];  
          }
          if(GAINSTEP2 != 0)  val+= x[newnr[nbolo]+htmp[l1].gi+htmp[l1].ib*GAINSTEP2]*htmp[l1].model;

          //R_ij_bis = val-s_X;
          Rij_bis = val;
          for(int k=0;k<MAXCHANNELS;k++){
            Rij_bis -= s_X[k]*htmp[l1].channels[k];
          }

          for(int k=0;k<MAXCHANNELS;k++){
            sum_channels[k]+= htmp[l1].channels[k]*Rij_bis;
          }
        }
      }


      for(int l1=0;l1<ndata;l1++){
        long ri1=htmp[l1].rg-globalBeginRing;  
       
         if (flg_rg[htmp[l1].ib][ri1]!=0) {
          ir=rgord[htmp[l1].ib][ri1]+newnr[htmp[l1].ib]; 

	  if (do_offset==1) {
	    val = x[ir]; // calcul val            
	  }
	  
          for(int m = 0;m<npixmap;m++){
            irt = newnr[nbolo]+htmp[l1].ib+nbolo*(m+GAINSTEP2);
            val+= x[irt]*htmp[l1].listofmap[m];
          }
          for(int m = 0;m<htmp[l1].nShpr;m++){
            irt = newnr[nbolo]+htmp[l1].ib+nbolo*(htmp[l1].listofShpr_idx[m]+npixmap+GAINSTEP2);
            val+= x[irt]*htmp[l1].listofShpr[m];
          }
          if(GAINSTEP2 != 0)  val+= x[newnr[nbolo]+htmp[l1].gi+htmp[l1].ib*GAINSTEP2 ]*htmp[l1].model;
          
          //Calcul Rij_bis
          Rij_bis = val;
          for(int k=0;k<MAXCHANNELS;k++){
            Rij_bis -= htmp[l1].channels[k]*s_X[k];
          }
   
          tmp = Rij_bis;
          for(int k=0;k<MAXCHANNELS;k++){                
            tmp -=sum_channels[k]*htmp[l1].alpha[k];
          }

	  if (do_offset==1) {
	    q2[ir]+= tmp;        //ie  ((val-s_X[k])-htmp[l1].alpha[k]*htmp[l1].channels[k]*R_ij_bis);
	  }
       
          for(int m=0;m<npixmap;m++){           
            q2[newnr[nbolo]+htmp[l1].ib+(m+GAINSTEP2)*nbolo]+= htmp[l1].listofmap[m]*tmp;          
          }
          for(int m=0;m<htmp[l1].nShpr;m++){           
            q2[newnr[nbolo]+htmp[l1].ib+(htmp[l1].listofShpr_idx[m]+npixmap+GAINSTEP2)*nbolo]+= htmp[l1].listofShpr[m]*tmp;          
          }

          if(GAINSTEP2 != 0){ 
            q2[newnr[nbolo]+htmp[l1].gi+htmp[l1].ib*GAINSTEP2] +=htmp[l1].model*tmp;            
          }

        }
      }
    }

  }


  if(rank == 0){
    double msum=1E4;
    
    for(int n = 0;n<Param->n_val_mean;n++){
      double sum2 = 0.0;
      for(int b = 0;b<nbolo*(npixShpr);b++){
	sum2+=(x[newnr[nbolo]+b+(GAINSTEP2)*nbolo]-Param->val_mean[n])*Param->do_mean[b+n*(nbolo*(npixShpr))];
      }
      for(int b = 0;b<nbolo*(npixShpr);b++){
	q2[newnr[nbolo]+b+(GAINSTEP2)*nbolo]+= sum2*Param->w_mean[n]*Param->do_mean[b+n*(nbolo*(npixShpr))];
      }
    }
    // NORMLISATION IS MANDATORY FOR MAP TEMPLATE
    for(int m =0;m<npixmap;m++){
        double sum2 = 0.0;
	for(int b = 0;b<nbolo;b++){
          sum2+=x[newnr[nbolo]+b+(m+GAINSTEP2)*nbolo];
        }
        
	for(int b = 0;b<nbolo;b++){
          q2[newnr[nbolo]+b+(m+GAINSTEP2)*nbolo]+= sum2*msum;
        }
    }
    
    // NORMLISATION GAIN MOYEN
    if(GAINSTEP2 != 0 && NORM_GAIN==1) {
      double sum2 = 0.0;
      for(int b = 0;b<nbolo*GAINSTEP2;b++){
        sum2+=x[newnr[nbolo]+b];
      }
        
      for(int b = 0;b<nbolo*GAINSTEP2;b++){
        q2[newnr[nbolo]+b]+= sum2*msum;
      }
    }

    if (do_offset==1) {
      // NORMLISATION OFFSET
      //offset
      for(int l1 = 0;l1 <newnr[nbolo];l1++){
	sum += x[l1]*msum;
      }
      for(int l1 = 0;l1 <newnr[nbolo];l1++){
	q2[l1] += sum*msum;
      }
#if 1
      // NORMALIZE OFFSET AGAINT CHANNELS
      for (int i=0;i<MAXCHANNELS;i++) {
	for (int j=0;j<newnr[nbolo];j++) q2[j]+=csum[i];
      }   
#endif
    }
  }

  free(csum);
}
// --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void invertMatrix(double * mat,double *vec,int n,int rank){
 

#ifdef CALCMATRIX      
      double *tmp = malloc(n*n*sizeof(double)); 
      memset(tmp,0,n*n*sizeof(double));    
      invert(mat,tmp,n);  
      memcpy(mat,tmp,n*n*sizeof(double));

#endif
      
      double *tmp_vec = malloc(n*sizeof(double));
      memset(tmp_vec,0,n*sizeof(double));  
      matvec(mat,vec,tmp_vec,n);
      memcpy(vec,tmp_vec,n*sizeof(double)); 
      
      free(tmp_vec);
#ifdef CALCMATRIX   
      free(tmp);
#endif

 
}
// #########################################################################################################################################################################################
void minimize_gain_tf(double *ix2,double *gaingi){
  //gaingi -> gain de chaque detecteur entré
  // xi2 -> vecteur des offsets + amplitude de templates

  // Init 
  long i,k,l1;
  int itermax = NUMBEROFITER;
  double tol = 1E-20;
  PIOLONG GAINSTEP2;
  int rank;
  int size;

  int rank_size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank_size);
  rank=rank_size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);

#ifdef CALCMATRIX   
  double MAT[MAXCHAN*MAXCHAN];
#endif

  struct timeval tp1,tp2;
  gettimeofday(&tp1,NULL);

  GAINSTEP2=GAINSTEP;

  nmatres=newnr[nbolo]+nbolo*(npixmap+npixShpr+GAINSTEP2);

  MPI_Bcast(&nmatres, sizeof(long), MPI_BYTE, 0, MPI_COMM_WORLD);

  double *x_tab = (double *) malloc(sizeof(double)*(nmatres));
  double *projX =(double *) malloc(sizeof(double)*(nmatres));
  double *new_p =(double *) malloc(sizeof(double)*(nmatres)); 
  double *new_x =(double *) malloc(sizeof(double)*(nmatres));
  double *new_r = (double *) malloc(sizeof(double)*(nmatres));
  double *res = (double *) malloc(sizeof(double)*(nmatres));


  if (rank==0) {
    fprintf(stderr,"==============================\n\nminimize_gain_nopol(): ITERATION NOPOL %ld/%d %ld \n\n==============================\n",itbogo,Param->NITT,(long) GAINSTEP);
    fprintf(stderr,"GAIN ");
    for (i=0;i<nbolo;i++) fprintf(stderr,"%lg ",gaingi[i*GAINSTEP]);
    fprintf(stderr,"\n");
    if (GAINSTEP>1) {
      fprintf(stderr,"...\nGAIN ");
      for (i=0;i<nbolo;i++) fprintf(stderr,"%lg ",gaingi[i*GAINSTEP+(GAINSTEP-1)]);
      fprintf(stderr,"\n");
    }
  }

  if (itbogo==0) delta0=0;      // init delta0
  MPI_Barrier(MPI_COMM_WORLD); //Synchronization MPI process
   
  //init x    
  memcpy(x_tab,ix2,sizeof(double)*nmatres);
  //memset(x_tab,0,sizeof(double)*nmatres);
 
  double SI[MAXCHAN]; // signal 
  double val_tmp = 0.0;


  // Init
  for (k=0;k<nnbpix;k++)  {                   // Pour chaque pixel in proc
    long ndata = loc_nhpix[k];               // nombre de donnees dans le pixel k
    hpix *htmp = loc_hpix[k]; 

    //init matrice pour pixel k 
#ifdef CALCMATRIX   
    memset(MAT,0,MAXCHANNELS*MAXCHANNELS*sizeof(double)); 
#endif
    memset(SI,0,MAXCHANNELS*sizeof(double));

    if (ndata>0) {
      for (l1=0;l1<ndata;l1++) {      
        long ri1=htmp[l1].rg-globalBeginRing; // check if ring valide
        if (flg_rg[htmp[l1].ib][ri1]!=0) {
        
#ifdef CALCMATRIX   
          for(int m=0;m<MAXCHANNELS;m++){
            for(int l =0;l<MAXCHANNELS;l++){
              MAT[m+MAXCHANNELS*l]+= htmp[l1].w * htmp[l1].channels[l]*htmp[l1].channels[m];   
            }
          }
#endif

          //if(rank == 0) for(int i = 0;i<MAXCHANNELS*MAXCHANNELS;i++) fprintf(stderr,"MAT[%d] = %lg\n",i,MAT[i]);
    
          //Calcul de carte        
          //long ri1=htmp[l1].rg-globalBeginRing; // check valide ring 
          double g1=gaingi[htmp[l1].gi+htmp[l1].ib*GAINSTEP]; //get gain 
          //long ir = rgord[htmp[l1].ib][ri1]+newnr[htmp[l1].ib];


          val_tmp =htmp[l1].sig*g1- htmp[l1].Sub_HPR-htmp[l1].corr_cnn;
	  if (REMOVE_CAL==1) {
	    val_tmp -=  htmp[l1].hpr_cal;
	  }
	  
          for(int l =0;l<MAXCHANNELS;l++){           
            SI[l]+= htmp[l1].w * htmp[l1].channels[l]*val_tmp;
          }

          htmp[l1].m = val_tmp; // calcul de m 
	}
      }
    }
   
     
#ifdef CALCMATRIX   
    //Inversion de matrice
    invertMatrix(MAT,SI,MAXCHANNELS,rank);
#else
    invertMatrix(imatrice+k*MAXCHANNELS*MAXCHANNELS,SI,MAXCHANNELS,rank);
#endif
    /*
      if(rank == 0){
      for(int i=0;i<MAXCHANNELS;i++){
      for(int j=0;j<MAXCHANNELS;j++){   
      fprintf(stderr,"MAT[%d] = %lf\n",i,MAT[j+MAXCHANNELS*i]);
      }
      }
      }
      exit(0);
    */
    // JMD: COMPUTE M INSIDE THE MASK BUT NOT Q2
    if (flgpix[k]>0) {
	    
      //Calcul ----- R_ij and alpha 
      for (l1=0;l1<ndata;l1++) { // parcours des donnees du pixel
	long ri1=htmp[l1].rg-globalBeginRing; 
	//long ir1 = rgord[htmp[l1].ib][ri1]+newnr[htmp[l1].ib];
	
	if (flg_rg[htmp[l1].ib][ri1]!=0) {                  //if pixel ok
	  htmp[l1].R_ij = htmp[l1].m ; //- SI;            
	  for(int i =0;i<MAXCHANNELS;i++){
	    htmp[l1].alpha[i] = 0.0;                  
	    htmp[l1].R_ij -= SI[i]*htmp[l1].channels[i];  
	    for(int j=0;j<MAXCHANNELS;j++){
#ifdef CALCMATRIX   
	      htmp[l1].alpha[i] += htmp[l1].w*MAT[j+MAXCHANNELS*i]*htmp[l1].channels[j];                //calcul alpha
#else
	      htmp[l1].alpha[i] += htmp[l1].w*imatrice[j+MAXCHANNELS*i+k*MAXCHANNELS*MAXCHANNELS]*htmp[l1].channels[j];                //calcul alpha
#endif
	    }  
	  }        
	}
      }
    }
  }

  
  
   
  //Init b2 and q2 to 0
  memset(b2,0,sizeof(double)*(nmatres));
  memset(hit2,0,sizeof(double)*(nmatres));
  memset(hit_WW,0,sizeof(double)*(nmatres));
  memset(hit_L1,0,sizeof(double)*(nmatres));
  memset(hit_L2,0,sizeof(double)*(nmatres));
  memset(q2,0,sizeof(double)*(nmatres));


  proj_data(b2,nnbpix,rank,nmatres,GAINSTEP2); //calcul b2
  proj_grad(q2,nmatres,x_tab,nnbpix,rank,GAINSTEP2); // Calcul q2 


  // Recuperation de b2
  #ifdef OPTIMPI
  {
    double *lb = (double *) malloc(sizeof(double)*(nmatres));
    MPI_Reduce(b2,lb,nmatres,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD); //parallisation MPI -- recuperation b
    memcpy(b2,lb,sizeof(double)*(nmatres)); // copy b2 dans lb
    free(lb);
  }
  #else
  if (rank==0) {
    double *lb = (double *) malloc(sizeof(double)*(nmatres));
    for (rrk=1;rrk<mpi_size;rrk++) {
      MPI_Recv(lb,sizeof(double)*(nmatres), MPI_BYTE, rrk,1031, MPI_COMM_WORLD,&statu);
      for (l=0;l<nmatres;l++) b2[l]+=lb[l];
    }
    free(lb);
  }
  else MPI_Send(b2, sizeof(double)*(nmatres), MPI_BYTE, 0, 1031, MPI_COMM_WORLD);
  #endif


  #ifdef OPTIMPI //recup and send hit2
  {
    double *lb = (double *) malloc(sizeof(double)*(nmatres));
    MPI_Reduce(hit2,lb,nmatres,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    memcpy(hit2,lb,sizeof(double)*(nmatres));
    free(lb);
  }
  #else
  if (rank==0) {
    double *lb = (double *) malloc(sizeof(double)*(nmatres));
    for (rrk=1;rrk<mpi_size;rrk++) {
      MPI_Recv(lb,sizeof(double)*(nmatres), MPI_BYTE, rrk,1033, MPI_COMM_WORLD,&statu);
      for (l=0;l<nmatres;l++) hit2[l]+=lb[l];
    }
    free(lb);
  }
  else {
    MPI_Send(hit2, sizeof(double)*(nmatres), MPI_BYTE, 0, 1033, MPI_COMM_WORLD); // nmatres = nb total d'inconnus = ring + n templates
  }
  #endif
  // END Recuperation de b2 and hit2


  //Send q2
  #ifdef OPTIMPI
  {
    double *lb = (double *) malloc(sizeof(double)*(nmatres));
    MPI_Reduce(q2,lb,nmatres,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    memcpy(q2,lb,sizeof(double)*(nmatres));
    free(lb);
  }
  #endif
  // End send q2


  if (rank==0)  fprintf(stderr,"B2 %lg\n",b2[0]);
  if (rank==0)  fprintf(stderr,"H2 %lg\n",hit2[0]);
  if (rank==0)  fprintf(stderr,"Q2 %lg\n",q2[0]);

  //init r2 and d2 to 0
  memset(r2,0,sizeof(double)*(nmatres)); 
  memset(d2,0,sizeof(double)*(nmatres));


  //Calcul r and p
  if (rank==0) {
    for (i=0; i < nmatres; i++){
      r2[i] = b2[i] - q2[i]; //r = b - Ax0 = Ax - Ax0  
      if (hit2[i]>1.0) d2[i] = r2[i]/hit2[i]; //d2 => p
      //if (rank==0)  fprintf(stderr,"[DEBUG] i = %d q2[] %lf b2[] =%lf r2[] = %lf ,d2[] = %lf , hit2[] = %lf \n",i,q2[i],b2[i],r2[i],d2[i],hit2[i]);

    }
  }


  MPI_Bcast(d2, sizeof(double)*nmatres, MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Bcast(r2, sizeof(double)*nmatres, MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Bcast(hit2, sizeof(double)*nmatres, MPI_BYTE, 0, MPI_COMM_WORLD);  

    
  memset(new_x,0,sizeof(double)*nmatres);
  memset(new_r,0,sizeof(double)*nmatres);
  memset(new_p,0,sizeof(double)*nmatres);
  memset(res,0,sizeof(double)*nmatres);

  
  //init
  double sum = 0.0 ,tmp = 0.0,alpha_tmp = 0.0 ,tmp_delta=0.0,beta= 0.0,delta = 0.0,alpha=0.0,tot_time = 0.0,time_exc =0.0;
    
  
  if (itbogo==0) delta0 = 0;

  if(rank==0) fprintf(stderr,"\n---------\n");
  for(int n=0;n<itermax;n++){
      
      gettimeofday(&tp1,NULL);
      
      //calcul projX
      memset(projX,0,sizeof(double)*nmatres);
      proj_grad(projX,nmatres,d2,nnbpix,rank,GAINSTEP2);

      //send projX
      double *lb = (double *) malloc(sizeof(double)*(nmatres));
      memset(lb,0,sizeof(double)*nmatres);
      MPI_Reduce(projX,lb,nmatres,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      memcpy(projX,lb,sizeof(double)*(nmatres));
      free(lb);

      MPI_Bcast(projX, sizeof(double)*nmatres, MPI_BYTE, 0, MPI_COMM_WORLD);


      //Calcul delta
      tmp = 0.0;
      tmp_delta = 0.0;
      for(int k = 0;k < nmatres;k++){
        tmp_delta += r2[k]*r2[k];
      }
      
      delta = tmp_delta;
      
      if (n==0) delta0 = delta;

      //calcul alpha
      alpha_tmp = 0.0;
      tmp = 0.0;
      
      for(int k = 0;k < nmatres;k++){
        alpha_tmp+= d2[k]*projX[k];
        if (hit2[k]>1.0) tmp += (r2[k]/hit2[k])*r2[k];
      }
      //alpha = delta /alpha_tmp;
      alpha = tmp /alpha_tmp;
      // manage case where alpha_tmp=0
      if (isnormal(alpha)) {
	      MPI_Bcast(&alpha, sizeof(double), MPI_BYTE, 0, MPI_COMM_WORLD);
      }
      else {
      	alpha=0;
      }
      //calcul new_x
      for(int k = 0;k < nmatres;k++){
        new_x[k] = x_tab[k]+alpha*d2[k];
      }
     
      //calcul new_r
      for(int k = 0;k <nmatres;k++){
        new_r[k] = r2[k]-alpha*projX[k];
      }

      //test delta
      if(delta<tol){
        gettimeofday(&tp2,NULL);
        time_exc = (double)(tp2.tv_sec-tp1.tv_sec)+(1E-6)*(tp2.tv_usec-tp1.tv_usec);
        tot_time += time_exc;

        if(rank ==0) fprintf(stderr,"\n==> End delta = %lg time = %3lfs\n",delta,tot_time);
        //for(int i =0;i<nmatres;i++) res[i] = new_x[i];

        memcpy(ix2,new_x,sizeof(double)*(nmatres));
        MPI_Bcast(ix2, sizeof(double)*(nmatres), MPI_BYTE, 0, MPI_COMM_WORLD);
        
        itbogo ++;       
        
        return;     
      }

      //calcul beta
      tmp = 0.0;
      sum =0.0;
      for(int i =0;i<nmatres;i++){
        if (hit2[i]>1.0) {
	  tmp += (new_r[i]/hit2[i])*new_r[i];
	  sum += (r2[i]/hit2[i])*r2[i];
	}
      }

      beta =tmp/sum;
      
    
      MPI_Bcast(&beta, sizeof(double), MPI_BYTE, 0, MPI_COMM_WORLD);

      //calcul new_p
      for(int k =0;k<nmatres;k++){
        if (hit2[k]>1.0) new_p[k] = (new_r[k]/hit2[k])+beta*d2[k];
      }

      //mise a jour r,p,x
      memcpy(r2,new_r,sizeof(double)*(nmatres));
      memcpy(d2,new_p,sizeof(double)*(nmatres));
      memcpy(x_tab,new_x,sizeof(double)*(nmatres));

      gettimeofday(&tp2,NULL);
      time_exc = (double)(tp2.tv_sec-tp1.tv_sec)+(1E-6)*(tp2.tv_usec-tp1.tv_usec);
      tot_time += time_exc;
      if (rank==0&&n%10==0) fprintf(stderr,"iter: %d/%d beta = %12lg  alpha = %12lg  delta = %12lg  %12lfs\n",n,itermax,beta,alpha,delta,time_exc);

      
  }
  if(rank==0) fprintf(stderr,"tot_time = %lg\n",tot_time);

  itbogo ++;
  memcpy(ix2,new_x,sizeof(double)*(nmatres));
  MPI_Bcast(ix2, sizeof(double)*(nmatres), MPI_BYTE, 0, MPI_COMM_WORLD);
      

}
#if 0
// --------------------------------------------------------------------------------------------------
void buildmap(double * map,double *signal,int begpix,int endpix){
 
  int rank;
  int size;
  int mpi_size;
  int rank_size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank_size);
  rank=rank_size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  mpi_size=size;

  int getbeginfo=0;
  int *allbeg;
  int *allend;
  int *all_realpix;
  float *all_map;

  // Convert array from double to float for smaller file to be written
  float *value = (float *)malloc(sizeof(float)*(endpix-begpix+1));
  for (int i = 0; i < (endpix-begpix+1); ++i) {
    value[i] = (float)signal[i];
  }


  if (getbeginfo==0) {
    getbeginfo=1;
    int i,rrk;
    if (rank==0) {
      allbeg = (int *) malloc(sizeof(int)*mpi_size);
      allend = (int *) malloc(sizeof(int)*mpi_size);
    }
    MPI_Gather(&begpix,sizeof(int),MPI_BYTE,allbeg,sizeof(int),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Gather(&endpix,sizeof(int),MPI_BYTE,allend,sizeof(int),MPI_BYTE,0,MPI_COMM_WORLD);

    if (rank==0) {
      for (rrk=0;rrk<mpi_size;rrk++) {
	      if (maxsize<allend[rrk]-allbeg[rrk]+1) maxsize=allend[rrk]-allbeg[rrk]+1;
      }
      all_realpix = (int *) malloc(sizeof(int)*mpi_size*maxsize);
      all_map = (float *) malloc(sizeof(float)*mpi_size*maxsize);
    }

    MPI_Bcast(&maxsize,sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD);

    int *l_idx =(int *)malloc(sizeof(int)*maxsize);
    for (i=begpix;i<=endpix;i++) l_idx[i-begpix]=realpix[i-begpix];

    MPI_Gather(l_idx,sizeof(int)*maxsize,MPI_BYTE,all_realpix,sizeof(int)*maxsize,MPI_BYTE,0,MPI_COMM_WORLD);

    free(l_idx);
  }
  float *l_map =(float *)malloc(sizeof(float)*maxsize);
  for (int k=begpix;k<=endpix;k++) l_map[k-begpix]=value[k-begpix];
  MPI_Gather(l_map,sizeof(float)*maxsize,MPI_BYTE,all_map,sizeof(float)*maxsize,MPI_BYTE,0,MPI_COMM_WORLD);
  free(l_map);


  if (rank==0) {
    int i,rrk;
    for (rrk=0;rrk<mpi_size;rrk++) {
      int l_beg,l_end;
      l_beg=allbeg[rrk];
      l_end=allend[rrk];
      for (i=l_beg;i<=l_end;i++) map[all_realpix[i-l_beg+rrk*maxsize]]=all_map[i-l_beg+rrk*maxsize];
    }
  }

}
// --------------------------------------------------------------------------------------------------
void foscat(double *x3,double *gain,int nside,int begpix,int endpix,int * do_templates){

  MPI_Status statu;
  int rank;
  int size;
  int mpi_size;
  int rank_size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank_size);
  rank=rank_size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  mpi_size=size;
  int all_ndata = 0;

  int getbeginfo=0;
  int *allbeg;
  int *allend;
  int *all_realpix;
  float *all_map;
  int maxsize=0;
  
  PIOLONG GAINSTEP2 = GAINSTEP;

  MPI_Barrier(MPI_COMM_WORLD); //Synchronization MPI process


  double ** signal = (double **) malloc(sizeof(double*)*MAXCHANNELS);
  for (int i = 0;i<MAXCHANNELS;i++){
    signal[i]=(double *) malloc(sizeof(double)*nnbpix);
  }
  

  // Init matrice and vecteur
  double *matrix = malloc(MAXCHANNELS*MAXCHANNELS*sizeof(double)); 
  double *Imatrix = malloc(MAXCHANNELS*MAXCHANNELS*sizeof(double));
  double *vector = malloc(MAXCHANNELS*sizeof(double));
  
  for (int k=0;k<nnbpix;k++) {
	
    for(int i =0;i<MAXCHANNELS;i++){
      signal[i][k]=UNSEENPIX;  // ie hp.UNSEEN
    }
    
    long ndata = loc_nhpix[k];
    hpix *htmp = loc_hpix[k];

	    
    //calcul map     
    memset(matrix,0,MAXCHANNELS*MAXCHANNELS*sizeof(double));
    memset(Imatrix,0,MAXCHANNELS*MAXCHANNELS*sizeof(double));
    memset(vector,0,MAXCHANNELS*sizeof(double));
	    
    for (int l1=0;l1<ndata;l1++) {
      // select only bolometers in detset
      
      long ri1=htmp[l1].rg-globalBeginRing;
	     
	      if (flg_rg[htmp[l1].ib][ri1]!=0) {
	        long iri1=rgord[htmp[l1].ib][ri1]+newnr[htmp[l1].ib];
	        //calcul signal corriger
	        double g1=gain[htmp[l1].gi+htmp[l1].ib*GAINSTEP];
	        double sig_corr = htmp[l1].sig*g1 - htmp[l1].Sub_HPR-htmp[l1].corr_cnn;
	        double sig_corr2 = sig_corr;
		if (do_offset==1) {
		  sig_corr-=x3[iri1];
		}
	        
	        if (REMOVE_CAL==1) {
	          sig_corr-=htmp[l1].hpr_cal;
	          sig_corr2-=htmp[l1].hpr_cal;
	        }
	        
	        for(int m =0;m<npixmap;m++){
	          sig_corr-=x3[newnr[nbolo]+htmp[l1].ib+(m+GAINSTEP2)*nbolo]*htmp[l1].listofmap[m];
	        }
	        for(int m =0;m<htmp[l1].nShpr;m++){
	          sig_corr-=x3[newnr[nbolo]+htmp[l1].ib+(htmp[l1].listofShpr_idx[m]+npixmap+GAINSTEP2)*nbolo]*htmp[l1].listofShpr[m];
	        }
		
	        sig_corr-=x3[newnr[nbolo]+htmp[l1].ib*GAINSTEP2]*htmp[l1].model;  
	        //calcul matrix & vector
	        for(int i = 0;i<MAXCHANNELS;i++){  
	          vector[i]+= htmp[l1].w *htmp[l1].channels[i]*sig_corr;
	        }
	        for(int i = 0;i<MAXCHANNELS;i++){        
	          for(int j = 0;j<MAXCHANNELS;j++){
	            matrix[i+j*MAXCHANNELS] += htmp[l1].w *htmp[l1].channels[j]*htmp[l1].channels[i];
	          }
	        } 
	      }
      }
	             
	     cond[k]=cond_thres(matrix,Imatrix,MAXCHANNELS);  

	     if (cond[k] < Param->seuilcond) {	               
	       invertMatrix(Imatrix,vector,MAXCHANNELS,rank);	       
	       for(int i = 0;i<MAXCHANNELS;i++){
	         signal[i][k]= vector[i];
	       }
	     }
      
	  }

  
  // #####  synch mpi rank and add map ##### // 

  int map_size = Nside*Nside*12; 

  //fprintf(stderr,"[DEBUG] rank = %d begpix %d endpix %d // %d %d\n",rank,begpix,endpix,begpix+nnbpix-1,map_size);
  
  double ** maps = (double **) malloc(sizeof(double*)*MAXCHANNELS);
  for (int i = 0;i<MAXCHANNELS;i++){
    maps[i]=(double *) malloc(sizeof(double)*map_size);
    for (int ii = 0;ii<map_size;ii++){
      maps[i][ii]=UNSEENPIX;
    }
  }

  for(int i =0;i<MAXCHANNELS;i++) buildmap(maps[i],signal[i],begpix,endpix);
  
  MPI_Barrier(MPI_COMM_WORLD);


  double ** new_templates = (double **) malloc(sizeof(double*)*2);
  for (int i =0;i<2;i++){
    new_templates[i]=(double *) malloc(sizeof(double)*map_size);
    for (int ii = 0;ii<map_size;ii++){
      new_templates[i][ii]=UNSEENPIX;
    }
  }

  if(rank == 0){

    PyObject* py_list = PyList_New(map_size);
    PyObject* py_nside,*py_rank,*py_ite,*func,*result;

    // For debug only 1 iteration sould be MAXCHANNELS send Q maps[1]
                                                                                                                                                                                                          
    for(int i =0;i<2;i++){
      for (int k = 0; k < map_size; k++) {
        PyObject* py_double = PyFloat_FromDouble(maps[i+1][k]);
        PyList_SET_ITEM(py_list,k, py_double);
      }

      //PyObject* model_path = PyUnicode_FromString("/export/home1/tfoulquier/out_maps/test_rstep2_353psb_nside64_I.fits");

      py_nside = Py_BuildValue("i",nside);
      py_ite = Py_BuildValue("i",i);
      py_rank = Py_BuildValue("i",rank);

      // Call the function in the Python script that will receive the array    
      func = MyPythonBackend.run_foscat;
      result = PyObject_CallObject(func, Py_BuildValue("(OOOO)", py_list,py_nside,py_rank,py_ite));

      
      //get result map from foscat      
      for(int m = 0;m<map_size;m++){
        new_templates[i][m] = PyFloat_AsDouble(PyList_GetItem(result,m));        
      }
    }

  
    
    // Clean up
    Py_DECREF(py_list);
    Py_DECREF(func);
    Py_DECREF(result);
    Py_DECREF(py_nside);
    Py_DECREF(py_ite);
    Py_DECREF(py_rank);
        
    
  }


  //Send templates // for debug only 1 iteration (only send Q) 
  for(int i =0;i<2;i++) MPI_Bcast(new_templates[i], sizeof(double)*(map_size), MPI_BYTE, 0, MPI_COMM_WORLD);

 
  // Update templates     
  for (int k=0;k<nnbpix;k++) {
    int index_pix = realpix[k];
    long ndata = loc_nhpix[k];
    hpix *htmp = loc_hpix[k];
    for (int l1=0;l1<ndata;l1++) {
      long ri1=htmp[l1].rg-globalBeginRing;
  
      
      if (flg_rg[htmp[l1].ib][ri1]!=0) {
          for (int m=0;m<npixShpr;m++)  {
              if (do_templates[m]==1) { // update PDUST
                htmp[l1].listofhpr[m]=htmp[l1].channels[1]*new_templates[0][index_pix]+htmp[l1].channels[2]*new_templates[1][index_pix];                
              }
              if (do_templates[m]==2) { // update ODUST
                htmp[l1].listofhpr[m]= -htmp[l1].channels[2]*new_templates[0][index_pix]+htmp[l1].channels[1]*new_templates[1][index_pix]; 
              }
          }      
      }
    }
  }

  free(new_templates);
  free(maps);
}
#endif
// --------------------------------------------------------------------------------------------------


double avv_poleff=-1;
PyObject *cnn_coef=NULL;
PyObject *py_signal;
PyObject *py_weights;
PyObject *py_TCO1;
PyObject *py_TSI1;
PyObject *py_hidx;
PyObject *py_idx;
PyObject *mynetwork;
PyObject *py_realpix;
PyObject *myrun;
PyObject *py_coef;
PyObject *py_smap;

int *recvcount;
int *recvcounts;
long nnbpixmax;
long all_ndatamax;
long offpix=0;

void localreduce(double *in,double *tmp,int nval)
{
  memset(tmp,0,sizeof(double)*nval);
  MPI_Allreduce(in,tmp,nval,MPI_DOUBLE, MPI_SUM,MPI_COMM_WORLD);
  memcpy(in,tmp,nval*sizeof(double));
}


void fit_cnn(double *x3,double *gain,int out_itt)
{
  
#if 0
  int rank,i,k,j,l,rrk,itt;
  int size;
  MPI_Status statu;
  int mpi_size;
  int rank_size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank_size);
  rank=rank_size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  mpi_size=size;
  int ib;
  int **rgcnn=(int **) malloc(sizeof(int *)*nbolo);

  for (i=0;i<nbolo;i++) {
    rgcnn[i]=(int *) malloc(sizeof(int)*(globalEndRing+1));
    FILE *l_fp=fopen(Param->rgcnn[i],"r");
	if(l_fp != NULL ){
	    fread(rgcnn[i],sizeof(int)*(globalEndRing+1),1,l_fp);
    }else{
		fprintf(stderr,"ERROR NULL File => %s\n",Param->rgcnn[i]);
	}
	fclose(l_fp);
  }

  double *A0 = (double *) malloc(sizeof(double)*12*CNN_NSIDE*CNN_NSIDE);
  double *A1 = (double *) malloc(sizeof(double)*12*CNN_NSIDE*CNN_NSIDE);
  double *A2 = (double *) malloc(sizeof(double)*12*CNN_NSIDE*CNN_NSIDE);
  double *B0 = (double *) malloc(sizeof(double)*12*CNN_NSIDE*CNN_NSIDE);
  double *B1 = (double *) malloc(sizeof(double)*12*CNN_NSIDE*CNN_NSIDE);
  double *B2 = (double *) malloc(sizeof(double)*12*CNN_NSIDE*CNN_NSIDE);
  double *C0 = (double *) malloc(sizeof(double)*12*CNN_NSIDE*CNN_NSIDE);
  double *C1 = (double *) malloc(sizeof(double)*12*CNN_NSIDE*CNN_NSIDE);
  double *C2 = (double *) malloc(sizeof(double)*12*CNN_NSIDE*CNN_NSIDE);
  double *D0 = (double *) malloc(sizeof(double)*12*CNN_NSIDE*CNN_NSIDE);
  double *D1 = (double *) malloc(sizeof(double)*12*CNN_NSIDE*CNN_NSIDE);
  double *D2 = (double *) malloc(sizeof(double)*12*CNN_NSIDE*CNN_NSIDE);
  double *V00 = (double *) malloc(sizeof(double)*12*CNN_NSIDE*CNN_NSIDE);
  double *V01 = (double *) malloc(sizeof(double)*12*CNN_NSIDE*CNN_NSIDE);
  double *V02 = (double *) malloc(sizeof(double)*12*CNN_NSIDE*CNN_NSIDE);
  double *V10 = (double *) malloc(sizeof(double)*12*CNN_NSIDE*CNN_NSIDE);
  double *V11 = (double *) malloc(sizeof(double)*12*CNN_NSIDE*CNN_NSIDE);
  double *V12 = (double *) malloc(sizeof(double)*12*CNN_NSIDE*CNN_NSIDE);

  memset(A0,0,sizeof(double)*12*CNN_NSIDE*CNN_NSIDE);
  memset(A1,0,sizeof(double)*12*CNN_NSIDE*CNN_NSIDE);
  memset(A2,0,sizeof(double)*12*CNN_NSIDE*CNN_NSIDE);
  memset(B0,0,sizeof(double)*12*CNN_NSIDE*CNN_NSIDE);
  memset(B1,0,sizeof(double)*12*CNN_NSIDE*CNN_NSIDE);
  memset(B2,0,sizeof(double)*12*CNN_NSIDE*CNN_NSIDE);
  memset(C0,0,sizeof(double)*12*CNN_NSIDE*CNN_NSIDE);
  memset(C1,0,sizeof(double)*12*CNN_NSIDE*CNN_NSIDE);
  memset(C2,0,sizeof(double)*12*CNN_NSIDE*CNN_NSIDE);
  memset(D0,0,sizeof(double)*12*CNN_NSIDE*CNN_NSIDE);
  memset(D1,0,sizeof(double)*12*CNN_NSIDE*CNN_NSIDE);
  memset(D2,0,sizeof(double)*12*CNN_NSIDE*CNN_NSIDE);
  memset(V00,0,sizeof(double)*12*CNN_NSIDE*CNN_NSIDE);
  memset(V01,0,sizeof(double)*12*CNN_NSIDE*CNN_NSIDE);
  memset(V02,0,sizeof(double)*12*CNN_NSIDE*CNN_NSIDE);
  memset(V10,0,sizeof(double)*12*CNN_NSIDE*CNN_NSIDE);
  memset(V11,0,sizeof(double)*12*CNN_NSIDE*CNN_NSIDE);
  memset(V12,0,sizeof(double)*12*CNN_NSIDE*CNN_NSIDE);

  int *reduce_realpix = (int *) malloc(sizeof(int)*nnbpix);
  int cnn_down = (Nside/CNN_NSIDE)*(Nside/CNN_NSIDE);
  for (i=0;i<nnbpix;i++){
    long ipix;
    ring2nest(Nside,realpix[i],&ipix);
    reduce_realpix[i]=ipix/cnn_down;
  }
  
  int CNN_XSIZE=Param->CNN_XSIZE;
  int CNN_YSIZE=Param->CNN_YSIZE;
  float coefresidu=Param->CNN_RESIDU;

  float *incorr = (float *) malloc(sizeof(float)*nbolo*CNN_XSIZE*CNN_YSIZE);
  float *wincorr = (float *) malloc(sizeof(float)*nbolo*CNN_XSIZE*CNN_YSIZE);

  memset(incorr,0,sizeof(float)*nbolo*CNN_XSIZE*CNN_YSIZE);
  memset(wincorr,0,sizeof(float)*nbolo*CNN_XSIZE*CNN_YSIZE);

  long all_ndata=0;
  int l1,m;

  for (k=0;k<nnbpix;k++) {
   // if (flgpix[k]>0) {
      long ndata = loc_nhpix[k];
      hpix *htmp = loc_hpix[k];
      for (l1=0;l1<ndata;l1++) {
	long ri1=htmp[l1].rg-globalBeginRing;
	if (flg_rg[htmp[l1].ib][ri1]!=0) {
	  all_ndata+=1;
	}
      }
  //  }
  }

  PyObject *arglist;

  PyObject *mynetworkSUB_HPR;
  
  if (out_itt==0) {
    if (python_rank==0) {
      arglist = Py_BuildValue("(l)", (long) 4*12*32*32);

      py_smap  = EXECPYTHON(PyObject_CallObject(MyPythonBackend.allocf32, arglist));
      
      arglist = Py_BuildValue("(s[s]llsss)", Param->CNN_WEIGHTS,pixnames[0],(long)rank,(long) Param->CNN_LEARN_PARAM,Param->CNN_TMPID,Param->MAP_CNN,Param->INST_CNN);
      if (nbolo==4)
	arglist = Py_BuildValue("(s[s,s,s,s]llsss)", Param->CNN_WEIGHTS,pixnames[0],pixnames[1],
				pixnames[2],pixnames[3],(long)rank,(long) Param->CNN_LEARN_PARAM,Param->CNN_TMPID,Param->MAP_CNN,Param->INST_CNN);
      if (nbolo==8)
	arglist = Py_BuildValue("(s[s,s,s,s,s,s,s,s]llsss)", Param->CNN_WEIGHTS,pixnames[0],pixnames[1],
				pixnames[2],pixnames[3],pixnames[4],pixnames[5],pixnames[6],pixnames[7],
				rank,(long) Param->CNN_LEARN_PARAM,Param->CNN_TMPID,Param->MAP_CNN,Param->INST_CNN);
      if (nbolo==11)
	arglist = Py_BuildValue("(s[s,s,s,s,s,s,s,s,s,s,s]llsss)", Param->CNN_WEIGHTS,pixnames[0],pixnames[1],
				pixnames[2],pixnames[3],pixnames[4],pixnames[5],pixnames[6],pixnames[7],
				pixnames[8],pixnames[9],pixnames[10],
				rank,(long) Param->CNN_LEARN_PARAM,Param->CNN_TMPID,Param->MAP_CNN,Param->INST_CNN);
      if (nbolo==12)
	arglist = Py_BuildValue("(s[s,s,s,s,s,s,s,s,s,s,s,s]llsss)", Param->CNN_WEIGHTS,pixnames[0],pixnames[1],
				pixnames[2],pixnames[3],pixnames[4],pixnames[5],pixnames[6],pixnames[7],
				pixnames[8],pixnames[9],pixnames[10],pixnames[11],
				rank,(long) Param->CNN_LEARN_PARAM,Param->CNN_TMPID,Param->MAP_CNN,Param->INST_CNN);



      mynetwork=EXECPYTHON(PyObject_CallObject(MyPythonBackend.initFrom_Files, arglist));
      //mynetwork=EXECPYTHON(PyObject_CallObject(MyPythonBackend.initFromFile, arglist));
      //Py_DECREF(arglist);
    }

    //fprintf(stderr,"nnbpix %d %d\n",(int) rank,(int) nnbpix);
    //fprintf(stderr,"all_ndata %d %d\n",(int) rank,(int) all_ndata);

    
    MPI_Reduce(&nnbpix,&nnbpixmax,1,MPI_LONG, MPI_SUM,0,python_comm);
    MPI_Reduce(&all_ndata,&all_ndatamax,1,MPI_LONG, MPI_SUM,0,python_comm);

    if (python_rank==0) {
      fprintf(stderr,"MAX %d %d NNBPIX %ld\n",python_rank,rank,(long) nnbpixmax);
      fprintf(stderr,"MAX %d %d ALL_NDATA %ld\n",python_rank,rank,(long) all_ndatamax);
    }
      
    recvcount = (int *) malloc(sizeof(int)*mpi_python_size);
    recvcounts = (int *) malloc(sizeof(int)*mpi_python_size);
    
    if (python_rank==0) {
      arglist = Py_BuildValue("(l)", (long) nnbpixmax);
      py_realpix = EXECPYTHON(PyObject_CallObject(MyPythonBackend.alloci32, arglist));
    }
      

    MPI_Allgather(&nnbpix,1,MPI_INT,recvcount,1,MPI_INT,python_comm);
    
     
    recvcounts[0]=0;
    for (i=1;i<mpi_python_size;i++) recvcounts[i]=recvcounts[i-1]+recvcount[i-1];

    offpix=recvcounts[python_rank];
    
    int *ires;
    if (python_rank==0) {
      ires=PyArray_DATA((PyArrayObject *)py_realpix);
    }
    else ires=NULL;
    
    MPI_Gatherv(realpix,nnbpix,MPI_INT,ires,
    		recvcount,recvcounts,MPI_INT,0,python_comm);

    //convert realpix in nested
    if (python_rank==0) {
      for (i=0;i<nnbpixmax;i++){
	long ipix;
	ring2nest(Nside,ires[i],&ipix);
	ires[i]=ipix;
      }
    }
    
    MPI_Allgather(&all_ndata,1,MPI_INT,recvcount,1,MPI_INT,python_comm);

    recvcounts[0]=0;
    for (i=1;i<mpi_python_size;i++) {
      recvcounts[i]=recvcounts[i-1]+recvcount[i-1];
      //if (rank==0) fprintf(stderr,"%d %d %d\n",i,recvcount[i],recvcounts[i]);
    }

 
    if (python_rank==0) {
      arglist = Py_BuildValue("(l)", (long) all_ndatamax);

      py_signal  = EXECPYTHON(PyObject_CallObject(MyPythonBackend.allocf32, arglist));
      py_weights = EXECPYTHON(PyObject_CallObject(MyPythonBackend.allocf32, arglist));
      py_TCO1    = EXECPYTHON(PyObject_CallObject(MyPythonBackend.allocf32, arglist));
      py_TSI1    = EXECPYTHON(PyObject_CallObject(MyPythonBackend.allocf32, arglist));
      py_hidx    = EXECPYTHON(PyObject_CallObject(MyPythonBackend.alloci32, arglist));
      py_idx     = EXECPYTHON(PyObject_CallObject(MyPythonBackend.alloci32, arglist));
    }    
  }

  float *signal   ;
  float *weights  ;
  float *TCO1     ;
  float *TSI1     ;
  int   *hidx     ;
  int   *idx      ;
  signal  = (float * ) malloc(sizeof(float)*all_ndata);
  weights = (float * ) malloc(sizeof(float)*all_ndata);
  TCO1    = (float * ) malloc(sizeof(float)*all_ndata);
  TSI1    = (float * ) malloc(sizeof(float)*all_ndata);
  hidx    = (int * )   malloc(sizeof(int)*all_ndata);
  idx     = (int * )   malloc(sizeof(int)*all_ndata);
  
  memset(signal,0,sizeof(float)*all_ndata);
  memset(weights,0,sizeof(float)*all_ndata);
  memset(TCO1,0,sizeof(float)*all_ndata);
  memset(TSI1,0,sizeof(float)*all_ndata);
  memset(hidx,0,sizeof(int)*all_ndata);
  memset(idx,0,sizeof(int)*all_ndata);

  all_ndata=0;

  double sig2[4];
  sig2[0]=0.0;
  sig2[1]=0.0;
  sig2[2]=0.0;
  sig2[3]=0.0;

  for (k=0;k<nnbpix;k++) {
    if (flgpix[k]>0&&cond[k]<Param->seuilcond) {
      double *vector = (double *) malloc(sizeof(double)*MAXCHANNELS*MAXCHANNELS);
      long ndata = loc_nhpix[k];
      hpix *htmp = loc_hpix[k];
     
      for (l1=0;l1<ndata;l1++) {
	long ri1=htmp[l1].rg-globalBeginRing;
	long iri1=rgord[htmp[l1].ib][ri1]+newnr[htmp[l1].ib];
	if (flg_rg[htmp[l1].ib][ri1]!=0) {
	  long iri1=rgord[htmp[l1].ib][ri1]+newnr[htmp[l1].ib];
	  //calcul signal corriger 
	  double sig_corr = htmp[l1].m;
	  
	  if (do_offset==1) {
	    sig_corr -= x3[iri1];
	  }
	 
	  for(int m =0;m<npixmap;m++){
	    sig_corr-=x3[newnr[nbolo]+htmp[l1].ib+(m+GAINSTEP2)*nbolo]*htmp[l1].listofmap[m]
	  }
	  for(int m =0;m<htmp[l1].nShpr;m++){
	    sig_corr-=x3[newnr[nbolo]+htmp[l1].ib+(htmp[l1].listofShpr_idx[m]+npixmap+GAINSTEP2)*nbolo]*htmp[l1].listofShpr[m]
	  }
	  sig_corr-=x3[newnr[nbolo]+htmp[l1].ib*GAINSTEP2]*htmp[l1].model;
	  
	  //calcul matrix & vector
	  for(int i = 0;i<MAXCHANNELS;i++){  
	    vector[i]+= htmp[l1].w *htmp[l1].channels[i]*sig_corr;
	  }  
	}
      }
      
      invertMatrix(imatrice+k*MAXCHANNELS*MAXCHANNELS,vector,MAXCHANNELS,rank);
      
      for (l1=0;l1<ndata;l1++) {
	long ri1=htmp[l1].rg-globalBeginRing;
	long iri1=rgord[htmp[l1].ib][ri1]+newnr[htmp[l1].ib];
	if (flg_rg[htmp[l1].ib][ri1]!=0) {
	  long iri1=rgord[htmp[l1].ib][ri1]+newnr[htmp[l1].ib];
	  //calcul signal corriger 
	  double sig_corr = htmp[l1].m+htmp[l1].corr_cnn; 
	  
	  if (do_offset==1) {
	    sig_corr -= x3[iri1];
	  }
	     
	  double temp=0;
	  for(int m =0;m<htmp[l1].nShpr;m++){
	    sig_corr-=x3[newnr[nbolo]+htmp[l1].ib+(htmp[l1].listofShpr_idx[m]+npixmap+GAINSTEP2)*nbolo]*htmp[l1].listofShpr[m];
	  }
	  
	  for(int m =0;m<npixShpr;m++){
	    if (DOCNN[m]==0) {
	      sig_corr-=x3[newnr[nbolo]+htmp[l1].ib+(m+GAINSTEP2)*nbolo]*htmp[l1].listofhpr[m]
	     }
	    else {
	      if (DOCNN[m]!=2) {
		temp+=htmp[l1].listofhpr[m]*x3[newnr[nbolo]+nbolo*(GAINSTEP)+m*nbolo+htmp[l1].ib];
	      }
	    }
	  }
	  sig_corr-=x3[newnr[nbolo]+htmp[l1].ib*GAINSTEP2]*htmp[l1].model;
	    
	  int iring=rgcnn[htmp[l1].ib][htmp[l1].rg];
	  int ipha=0.0 ; //(int)(CNN_XSIZE*(htmp[l1].External)/(2*M_PI));
	  if (ipha==CNN_XSIZE) ipha=CNN_XSIZE-1;
	  if (ipha<0||ipha>CNN_XSIZE-1) {
	    fprintf(stderr,"CNN_XSIZE PBS %d\n",ipha);
	    exit(0);
	  }
	    
	  // if nopol mapq and mapu have been set to 0 
	  signal[all_ndata]=sig_corr; //-coefresidu*(mapi+CO1*mapq+SI1*mapu);
	  if (flgpix[k]) {
	      sig2[0]+=sig_corr;
	      sig2[1]+=htmp[l1].w;
	      sig2[3]+=htmp[l1].w*(sig_corr-coefresidu*(mapi+CO1*mapq+SI1*mapu));
	      for (int k=0;k<MAXCHANNELS;k++) {
		sig2[0]-=vector[k]*htmp[l1].channels[k];
		sig2[3]-=coefresidu*htmp[l1].w*vector[k]*htmp[l1].channels[k];
	      }
	      sig2[0]=htmp[l1].w*sig2[0]*sig2[0];
	   
	      weights[all_ndata]=htmp[l1].w;

	      // compute polarised map for cross-ST input
	      double a=sig_corr-temp-htmp[l1].corr_cnn-SI/SII;
	      double b=CO1-SIQ/SII;
	      double c=SI1-SIU/SII;
	    
	      V00[reduce_realpix[k]] += htmp[l1].w*a*b;
	      V10[reduce_realpix[k]] += htmp[l1].w*a*c;
	      A0[reduce_realpix[k]]  += htmp[l1].w*b*b;
	      B0[reduce_realpix[k]]  += htmp[l1].w*b*c;
	      C0[reduce_realpix[k]]  += htmp[l1].w*c*b;
	      D0[reduce_realpix[k]]  += htmp[l1].w*c*c;

	      //if (htmp[l1].surv%10==1||htmp[l1].surv%10==3) {
	      if (htmp[l1].surv<10) {
		V01[reduce_realpix[k]] += htmp[l1].w*a*b;
		V11[reduce_realpix[k]] += htmp[l1].w*a*c;
		A1[reduce_realpix[k]]  += htmp[l1].w*b*b;
		B1[reduce_realpix[k]]  += htmp[l1].w*b*c;
		C1[reduce_realpix[k]]  += htmp[l1].w*c*b;
		D1[reduce_realpix[k]]  += htmp[l1].w*c*c;
	      }

	      //if (htmp[l1].surv%10==2||htmp[l1].surv%10==4) {
	      if (htmp[l1].surv>=10) {
		V02[reduce_realpix[k]] += htmp[l1].w*a*b;
		V12[reduce_realpix[k]] += htmp[l1].w*a*c;
		A2[reduce_realpix[k]]  += htmp[l1].w*b*b;
		B2[reduce_realpix[k]]  += htmp[l1].w*b*c;
		C2[reduce_realpix[k]]  += htmp[l1].w*c*b;
		D2[reduce_realpix[k]]  += htmp[l1].w*c*c;
	      }
	      
	      
	    }
	    
	    TCO1[all_ndata]=CO1;
	    TSI1[all_ndata]=SI1;
	    hidx[all_ndata]=k+offpix;
	    int l_idx=htmp[l1].ib*CNN_YSIZE*CNN_XSIZE+iring*CNN_XSIZE+ipha;

      //fprintf(stderr,"l_idx %d nbolo %d CNN_YSIZE %d CNN_XSIZE %d \n",(int) l_idx,(int)nbolo,(int) CNN_YSIZE,(int) CNN_XSIZE);

	    if (l_idx<0||l_idx>=nbolo*CNN_YSIZE*CNN_XSIZE) {      
	      fprintf(stderr,"IDX ERROR %d %d %d\n",(int) htmp[l1].ib,(int) iring,(int) ipha);
	      exit(0);
	    }
	    idx[all_ndata]=l_idx;
	    incorr[l_idx]+=htmp[l1].w*(sig_corr-(mapi+CO1*mapq+SI1*mapu));
	    wincorr[l_idx]+=htmp[l1].w;
	    all_ndata=all_ndata+1;
	 // }
	}
      }
   // }
  }

  double *ltmp=(double *) malloc(sizeof(double)*12*CNN_NSIDE*CNN_NSIDE);
  localreduce(A0,ltmp,12*CNN_NSIDE*CNN_NSIDE);
  localreduce(A1,ltmp,12*CNN_NSIDE*CNN_NSIDE);
  localreduce(A2,ltmp,12*CNN_NSIDE*CNN_NSIDE);
  localreduce(B0,ltmp,12*CNN_NSIDE*CNN_NSIDE);
  localreduce(B1,ltmp,12*CNN_NSIDE*CNN_NSIDE);
  localreduce(B2,ltmp,12*CNN_NSIDE*CNN_NSIDE);
  localreduce(C0,ltmp,12*CNN_NSIDE*CNN_NSIDE);
  localreduce(C1,ltmp,12*CNN_NSIDE*CNN_NSIDE);
  localreduce(C2,ltmp,12*CNN_NSIDE*CNN_NSIDE);
  localreduce(D0,ltmp,12*CNN_NSIDE*CNN_NSIDE);
  localreduce(D1,ltmp,12*CNN_NSIDE*CNN_NSIDE);
  localreduce(D2,ltmp,12*CNN_NSIDE*CNN_NSIDE);
  localreduce(V00,ltmp,12*CNN_NSIDE*CNN_NSIDE);
  localreduce(V01,ltmp,12*CNN_NSIDE*CNN_NSIDE);
  localreduce(V02,ltmp,12*CNN_NSIDE*CNN_NSIDE);
  localreduce(V10,ltmp,12*CNN_NSIDE*CNN_NSIDE);
  localreduce(V11,ltmp,12*CNN_NSIDE*CNN_NSIDE);
  localreduce(V12,ltmp,12*CNN_NSIDE*CNN_NSIDE);
  free(ltmp);

  for (i=0;i<12*CNN_NSIDE*CNN_NSIDE;i++) {
    double det0=A0[i]*D0[i]-B0[i]*C0[i];
    double det1=A1[i]*D1[i]-B1[i]*C1[i];
    double det2=A2[i]*D2[i]-B2[i]*C2[i];

    if (det0>0) {
      double a=( D0[i]*V00[i] - B0[i]*V10[i])/det0;
      double b=(-C0[i]*V00[i] + A0[i]*V10[i])/det0;
      V00[i]=a;
      V10[i]=b;
    }
    else {
      V00[i]=UNSEENPIX;
      V10[i]=UNSEENPIX;
    }
    
    if (det1>0) {
      double a=( D1[i]*V01[i] - B1[i]*V11[i])/det1;
      double b=(-C1[i]*V01[i] + A1[i]*V11[i])/det1;
      V01[i]=a;
      V11[i]=b;
    }
    else {
      V01[i]=UNSEENPIX;
      V11[i]=UNSEENPIX;
    }

    if (det2>0) {
      double a=( D2[i]*V02[i] - B2[i]*V12[i])/det2;
      double b=(-C2[i]*V02[i] + A2[i]*V12[i])/det2;
      V02[i]=a;
      V12[i]=b;
    }
    else {
      V02[i]=UNSEENPIX;
      V12[i]=UNSEENPIX;
    }
  }
#if 1
  if (python_rank==0) {
      float *tp_val = PyArray_DATA((PyArrayObject *) py_smap);
      
      for (i=0;i<12*32*32;i++) tp_val[i*2]              = V01[i];
      for (i=0;i<12*32*32;i++) tp_val[i*2+1]            = V11[i];
      for (i=0;i<12*32*32;i++) tp_val[i*2+12*32*32*2]   = V02[i];
      for (i=0;i<12*32*32;i++) tp_val[i*2+1+12*32*32*2] = V12[i];
      if (rank==0) {
	char l_path[1024];
	sprintf(l_path,"MAPSDUST2_%d.dat",out_itt);
	FILE *fp=fopen(l_path,"w");
	fwrite(tp_val,sizeof(float)*4*12*CNN_NSIDE*CNN_NSIDE,1,fp);
	fclose(fp);
      }
  }
  
#else
  if (rank==0) {
    FILE*fp;

    fp=fopen("MAPQ0.dat","w");
    fwrite(V00,sizeof(double)*12*CNN_NSIDE*CNN_NSIDE,1,fp);
    fclose(fp);
    fp=fopen("MAPU0.dat","w");
    fwrite(V10,sizeof(double)*12*CNN_NSIDE*CNN_NSIDE,1,fp);
    fclose(fp);
    fp=fopen("MAPQ1.dat","w");
    fwrite(V01,sizeof(double)*12*CNN_NSIDE*CNN_NSIDE,1,fp);
    fclose(fp);
    fp=fopen("MAPU1.dat","w");
    fwrite(V11,sizeof(double)*12*CNN_NSIDE*CNN_NSIDE,1,fp);
    fclose(fp);
    fp=fopen("MAPQ2.dat","w");
    fwrite(V02,sizeof(double)*12*CNN_NSIDE*CNN_NSIDE,1,fp);
    fclose(fp);
    fp=fopen("MAPU2.dat","w");
    fwrite(V12,sizeof(double)*12*CNN_NSIDE*CNN_NSIDE,1,fp);
    fclose(fp);
    exit(0);
  }
#endif

  
  //Keep for debug
  float *l_incorr = (float *) malloc(sizeof(float)*nbolo*CNN_XSIZE*CNN_YSIZE);
  float *l_wincorr = (float *) malloc(sizeof(float)*nbolo*CNN_XSIZE*CNN_YSIZE);

  memset(l_incorr ,0, sizeof(float)*nbolo*CNN_XSIZE*CNN_YSIZE);
  memset(l_wincorr ,0, sizeof(float)*nbolo*CNN_XSIZE*CNN_YSIZE);

  MPI_Reduce(incorr,l_incorr,nbolo*CNN_XSIZE*CNN_YSIZE,MPI_FLOAT, MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(wincorr,l_wincorr,nbolo*CNN_XSIZE*CNN_YSIZE,MPI_FLOAT, MPI_SUM,0,MPI_COMM_WORLD);

  double ssig2[4];
  sig2[2]=all_ndata;
  sig2[1]=0.0;
  sig2[3]=0.0;
  for (i=0;i<nbolo*CNN_XSIZE*CNN_YSIZE;i++) {
    if (l_wincorr[i]>0) {
      sig2[1]+=l_wincorr[i];
      sig2[3]+=l_incorr[i];
      l_incorr[i]/=l_wincorr[i];
    }
    else  {
      l_incorr[i]=0.0;
    }
  }

  MPI_Bcast(l_incorr,nbolo*CNN_XSIZE*CNN_YSIZE,MPI_FLOAT,0,MPI_COMM_WORLD);
	    
  MPI_Allreduce(sig2,ssig2,4,MPI_DOUBLE, MPI_SUM,MPI_COMM_WORLD);

  double offregul=ssig2[3]/ssig2[1];

#if 0
  if (rank==0) {
    fprintf(stderr,"SAVE CNN CORRECTION in %s_CNN_IN_CORR\n",Param->Out_VEC[0]);
  
    PIOSTRING saveval;
    sprintf(saveval,"%s_CNN_IN_CORR",Param->Out_VEC[0]);
    FILE *fp=fopen(saveval,"w");
    if(fp != NULL){
      fwrite(l_incorr,1,sizeof(float)*nbolo*CNN_XSIZE*CNN_YSIZE,fp);
    }else{
      fprintf(stderr,"ERROR NULL file => %s \n",Param->Out_VEC[0]);
    }
    fclose(fp);
    sprintf(saveval,"%s_CNN_IN_NCORR",Param->Out_VEC[0]);
    fp=fopen(saveval,"w");
    fwrite(l_wincorr,1,sizeof(float)*nbolo*CNN_XSIZE*CNN_YSIZE,fp);
    fclose(fp);
  }
#endif
#if 0
  if (rank==0) {
    fprintf(stderr,"SAVE FOR CNN %ld DATA INSIDE %s_CNNINFO_*\n",(long) all_ndata,Param->Out_VEC[0]);
  }

  PIOSTRING saveval;
  sprintf(saveval,"%s_CNNINFO_RPIX_%d",Param->Out_VEC[0],rank);
  FILE *fp=fopen(saveval,"w");
  fwrite(realpix,1,sizeof(int)*nnbpix,fp);
  fclose(fp);

  sprintf(saveval,"%s_CNNINFO_SIG_%d",Param->Out_VEC[0],rank);
  fp=fopen(saveval,"w");
  fwrite(signal,1,sizeof(float)*all_ndata,fp);
  fclose(fp);

  sprintf(saveval,"%s_CNNINFO_W_%d",Param->Out_VEC[0],rank);
  fp=fopen(saveval,"w");
  fwrite(weights,1,sizeof(float)*all_ndata,fp);
  fclose(fp);
  
  sprintf(saveval,"%s_CNNINFO_CO_%d",Param->Out_VEC[0],rank);
  fp=fopen(saveval,"w");
  fwrite(TCO1,1,sizeof(float)*all_ndata,fp);
  fclose(fp);
  
  sprintf(saveval,"%s_CNNINFO_SI_%d",Param->Out_VEC[0],rank);
  fp=fopen(saveval,"w");
  fwrite(TSI1,1,sizeof(float)*all_ndata,fp);
  fclose(fp);
  
  sprintf(saveval,"%s_CNNINFO_HIDX_%d",Param->Out_VEC[0],rank);
  fp=fopen(saveval,"w");
  fwrite(hidx,1,sizeof(int)*all_ndata,fp);
  fclose(fp);
  
  sprintf(saveval,"%s_CNNINFO_IDX_%d",Param->Out_VEC[0],rank);
  fp=fopen(saveval,"w");
  fwrite(idx,1,sizeof(int)*all_ndata,fp);
  fclose(fp);

  MPI_Barrier(MPI_COMM_WORLD);
  exit(0);
#endif
  
  for (k=0;k<all_ndata;k++) {
    if (weights[k]!=0) {
      signal[k]-=offregul;
    }
    else signal[k]=0.0;
  }
  
  free(incorr);
  free(l_wincorr);
  free(wincorr);

  double rap=sqrt(ssig2[0]/ssig2[1])*sqrt((nbolo*CNN_YSIZE*CNN_XSIZE)/ssig2[2]);

  if (rank==0) {
    fprintf(stderr,"SIGMA: %lf %lf\n",(double) rap,(double) sqrt((nbolo*CNN_YSIZE*CNN_XSIZE)/ssig2[2]));
  }

  rap=1.0;
  
  int ncorr=0,testcorr=0;
  if (avv_poleff==-1) {
    testcorr=-1;
    avv_poleff=0;
  }
  for (m=0;m<npixShpr;m++)  {
    if (DOCNN[m]!=0) ncorr++;
  }

  if (out_itt==0) {
    if (python_rank==0) {
      arglist = Py_BuildValue("(l)", (long) (ncorr*nbolo));
      if (NOMOREFITTED==0) {
	      cnn_coef  = EXECPYTHON(PyObject_CallObject(MyPythonBackend.allocf32, arglist));
	      float *l_coef=PyArray_DATA((PyArrayObject *) cnn_coef);
	      for (m=0;m<npixShpr;m++)  {
	        if (DOCNN[m]!=0) {
	          for (j=0;j<nbolo;j++) {
	            l_coef[j+(DOCNN[m]-1)*nbolo]=x3[newnr[nbolo]+nbolo*(GAINSTEP)+m*nbolo+j];
	          }
	        }
	      }
      }

      py_coef = cnn_coef;
      float *tp_coef = PyArray_DATA((PyArrayObject *) py_coef);
  
      for (m=0;m<npixShpr;m++)  {
	      if (DOCNN[m]!=0) {
	        if (rank==0) {
	          fprintf(stderr,"iPOLEFF %3d ",DOCNN[m]);
	        }
	        for (j=0;j<nbolo;j++) {
	          if (rank==0) {
	            fprintf(stderr,"%.3f ",tp_coef[j+(DOCNN[m]-1)*nbolo]);
	          }
	        }
	        if (rank==0) {
	          fprintf(stderr,"\n");
	        }
	      }
      }
    }
    
    float *ires=NULL;

    if (python_rank==0) ires=PyArray_DATA((PyArrayObject *) py_signal);
    MPI_Gatherv(signal,all_ndata,MPI_FLOAT,ires,
    		recvcount,recvcounts,MPI_FLOAT,0,python_comm);
    if (python_rank==0)  ires=PyArray_DATA((PyArrayObject *) py_weights);
    MPI_Gatherv(weights,all_ndata,MPI_FLOAT,ires,
    		recvcount,recvcounts,MPI_FLOAT,0,python_comm);
    if (python_rank==0)  ires=PyArray_DATA((PyArrayObject *) py_TCO1);
    MPI_Gatherv(TCO1,all_ndata,MPI_FLOAT,ires,
    		recvcount,recvcounts,MPI_FLOAT,0,python_comm);
    if (python_rank==0)  ires=PyArray_DATA((PyArrayObject *) py_TSI1);
    MPI_Gatherv(TSI1,all_ndata,MPI_FLOAT,ires,
    		recvcount,recvcounts,MPI_FLOAT,0,python_comm);
    if (python_rank==0)  ires=PyArray_DATA((PyArrayObject *) py_hidx);
    MPI_Gatherv(hidx,all_ndata,MPI_INT,ires,
    		recvcount,recvcounts,MPI_INT,0,python_comm);
    if (python_rank==0)  ires=PyArray_DATA((PyArrayObject *) py_idx);
    MPI_Gatherv(idx,all_ndata,MPI_INT,ires,
    		recvcount,recvcounts,MPI_INT,0,python_comm);

    if (python_rank==0) { 
      arglist = PyTuple_New(10);
      PyTuple_SetItem(arglist,0,mynetwork);
      PyTuple_SetItem(arglist,1,py_signal);
      PyTuple_SetItem(arglist,2,py_weights);
      PyTuple_SetItem(arglist,3,py_TCO1);
      PyTuple_SetItem(arglist,4,py_TSI1);
      PyTuple_SetItem(arglist,5,py_hidx);
      PyTuple_SetItem(arglist,6,py_idx);
      PyTuple_SetItem(arglist,7,py_realpix);
      PyTuple_SetItem(arglist,8,py_coef);
      PyTuple_SetItem(arglist,9,py_smap);
      myrun = EXECPYTHON(PyObject_CallObject(MyPythonBackend.init_net, arglist));
    }

    
  }
  else {
    float *ires=NULL;
    
    if (python_rank==0) ires=PyArray_DATA((PyArrayObject *) py_signal);
    MPI_Gatherv(signal,all_ndata,MPI_FLOAT,ires,
    		recvcount,recvcounts,MPI_FLOAT,0,python_comm);

    if (python_rank==0) {
      arglist = PyTuple_New(6);
      PyTuple_SetItem(arglist,0,myrun);
      PyTuple_SetItem(arglist,1,py_signal);
      PyTuple_SetItem(arglist,2,py_weights);
      PyTuple_SetItem(arglist,3,py_hidx);
      PyTuple_SetItem(arglist,4,py_realpix);
      PyTuple_SetItem(arglist,5,py_smap);
      myrun = EXECPYTHON(PyObject_CallObject(MyPythonBackend.init_net_data, arglist));
    }
  }

  free(signal  );
  free(weights );
  free(TCO1    );
  free(TSI1    );
  free(hidx    );
  free(idx     );

  if (python_rank==0) {
    struct timeval tp1,tp2;
    gettimeofday(&tp1,NULL);
    double deltaloss=1.0;
    itt=0;
    while (itt<Param->CNN_ITT) {
      int k;
      Py_ssize_t nval;
      PyArrayObject *value;
      
      arglist = PyTuple_New(1);
      PyTuple_SetItem(arglist,0,myrun);
      PyObject *mygradient = EXECPYTHON(PyObject_CallObject(MyPythonBackend.grad, arglist));
      PyTypeObject *keys = (PyTypeObject *) PyDict_Keys(mygradient);
      nval = PyList_GET_SIZE(keys);
      for (k=0;k<nval;k++) {
	char keystr[128];
	sprintf(keystr,"%03d",k);
	value = (PyArrayObject *)EXECPYTHON(PyDict_GetItemString(mygradient,keystr));
	long nnn=1;
	int l;
	for (l=0;l<PyArray_NDIM(value);l++) nnn*=PyArray_DIM(value,l);
	float *res = (float *) malloc(sizeof(float)*nnn);
	MPI_Allreduce(PyArray_DATA(value),res,nnn,MPI_FLOAT,MPI_SUM,tensorflow_comm);
	memcpy(PyArray_DATA(value),res,sizeof(float)*nnn);
      }
      
      PyObject *l_arglist = PyTuple_New(2);
      PyTuple_SetItem(l_arglist,0,myrun);
      PyTuple_SetItem(l_arglist,1,mygradient);
      EXECPYTHON(PyObject_CallObject(MyPythonBackend.agrad, l_arglist));
      
      if (itt%10==0) {
	PyObject *res=EXECPYTHON(PyObject_CallObject(MyPythonBackend.gloss, arglist)); 
	if (res==NULL) {
	  PyErr_Print();
	  MPI_Finalize();
	  exit(0)
	}
	float *loss=PyArray_DATA((PyArrayObject *)res);
	float l_loss[2];
	
	MPI_Allreduce(loss,l_loss,2,MPI_FLOAT,MPI_SUM,tensorflow_comm);
	gettimeofday(&tp2,NULL);
	double dt=(double)(tp2.tv_sec-tp1.tv_sec)+(1E-6)*(tp2.tv_usec-tp1.tv_usec);
	if (itt==0) deltaloss=sqrt(l_loss[0]/l_loss[1]);
	else deltaloss=fabs(deltaloss-sqrt(l_loss[0]/l_loss[1]));
	if (rank==0) {
	  fprintf(stderr,"Itt %d loss=%10.4g Dloss=%10.4g Lr=%.4f Lpwst=%.4g Dt=%.4f\n",(int) itt,sqrt(l_loss[0]/l_loss[1]),deltaloss,loss[2],loss[3],dt);
	}
	deltaloss=sqrt(l_loss[0]/l_loss[1]);
	gettimeofday(&tp1,NULL);
	//Py_DECREF(res);
      }
      itt++;
    }
  }
  // la prediction doit avoir la taille de l'echantillon : all_ndata - change value -- have to be the size of timeline
  // tous les processeurs applique prdiction
  //nbolo*CNN_SIZE*CNN_YSIZE

  //Py_DECREF(arglist);
  PyObject *py_num;
  if (python_rank==0) {
    arglist = Py_BuildValue("(l)", (long) 2);
    py_num  = EXECPYTHON(PyObject_CallObject(MyPythonBackend.alloci32, arglist));
    //Py_DECREF(arglist);
  }
  float *icorrection=NULL;
  float *correction=NULL;
  float *iprediction[2];
  float *prediction[2];
  
  if (python_rank==0) {
    int *iter_docnn = (int *) (PyArray_DATA((PyArrayObject *)py_num));
    *iter_docnn = 0;
      
    arglist = PyTuple_New(2);
    PyTuple_SetItem(arglist,0,myrun);
    PyTuple_SetItem(arglist,1,py_num);
    
    PyObject *thepred=EXECPYTHON(PyObject_CallObject(MyPythonBackend.corr, arglist));
    icorrection = (float *) malloc(sizeof(float)*all_ndatamax);
    memcpy(icorrection,PyArray_DATA((PyArrayObject *)thepred),sizeof(float)*all_ndatamax);
    arglist = PyTuple_New(1);
    PyTuple_SetItem(arglist,0,thepred);
    PyObject *nulpred=EXECPYTHON(PyObject_CallObject(MyPythonBackend.Clean, arglist));
    
    *iter_docnn = 1;
    arglist = PyTuple_New(2);
    PyTuple_SetItem(arglist,0,myrun);
    PyTuple_SetItem(arglist,1,py_num);
    thepred=EXECPYTHON(PyObject_CallObject(MyPythonBackend.corr, arglist));
    iprediction[0] = (float *) malloc(sizeof(float)*all_ndatamax);
    memcpy(iprediction[0],PyArray_DATA((PyArrayObject *)thepred),sizeof(float)*all_ndatamax);
    arglist = PyTuple_New(1);
    PyTuple_SetItem(arglist,0,thepred);
    nulpred=EXECPYTHON(PyObject_CallObject(MyPythonBackend.Clean, arglist));
    
    *iter_docnn = 2;
    arglist = PyTuple_New(2);
    PyTuple_SetItem(arglist,0,myrun);
    PyTuple_SetItem(arglist,1,py_num);
    thepred=EXECPYTHON(PyObject_CallObject(MyPythonBackend.corr, arglist));
    iprediction[1] = (float *) malloc(sizeof(float)*all_ndatamax);
    memcpy(iprediction[1],PyArray_DATA((PyArrayObject *)thepred),sizeof(float)*all_ndatamax);
    arglist = PyTuple_New(1);
    PyTuple_SetItem(arglist,0,thepred);
    nulpred=EXECPYTHON(PyObject_CallObject(MyPythonBackend.Clean, arglist));
  }
  
  correction = (float *) malloc(sizeof(float)*all_ndata);
    
  // and now distribute data
  MPI_Scatterv(icorrection, recvcount,recvcounts,MPI_FLOAT,
	       correction, all_ndata,MPI_FLOAT,0,python_comm);

  if (python_rank==0) {
    free(icorrection);
  }

  for (i=0;i<2;i++) {
    prediction[i] = (float *) malloc(sizeof(float)*all_ndata);
    
    // and now distribute data
    MPI_Scatterv(iprediction[i], recvcount,recvcounts,MPI_FLOAT,
		 prediction[i], all_ndata,MPI_FLOAT,0,python_comm);

    if (python_rank==0) {
      free(iprediction[i]);
    }
  }
  //Py_DECREF(py_num);

  
  // and now replace the template by the CNN version

  long iii=0;
  for (k=0;k<nnbpix;k++) {
    long ndata = loc_nhpix[k];
    hpix *htmp = loc_hpix[k];
    for (l1=0;l1<ndata;l1++) {
      
      long ri1=htmp[l1].rg-globalBeginRing;
      
      if (flg_rg[htmp[l1].ib][ri1]!=0) {
	
	int iring = rgcnn[htmp[l1].ib][htmp[l1].rg];
	
	int ipha=(int)(0); //CNN_XSIZE*(htmp[l1].External)/(2*M_PI));
	
	if (ipha==CNN_XSIZE) ipha=CNN_XSIZE-1;

	for (m=0;m<npixShpr;m++)  {
	    if (DOCNN[m]==1) {
	      htmp[l1].listofhpr[m]=sqrt(-2*log( drand48()))*cos(2*M_PI*drand48()); // should be set to something otherwise it generates Nan during the minimization
	    }
	    if (DOCNN[m]==2) {
	      htmp[l1].listofhpr[m]=prediction[0][iii]; // should be set to something otherwise it generates Nan during the minimization
	    }
	    if (DOCNN[m]==3) {
	      htmp[l1].listofhpr[m]=prediction[1][iii]; // should be set to something otherwise it generates Nan during the minimization
	    }
	}
	//	htmp[l1].corr_cnn=l_incorr[htmp[l1].ib*CNN_YSIZE*CNN_XSIZE+iring*CNN_XSIZE+ipha];
	  //htmp[l1].corr_cnn=l_incorr[htmp[l1].ib*CNN_YSIZE*CNN_XSIZE+iring*CNN_XSIZE+ipha];
	//}
	htmp[l1].corr_cnn=correction[iii]; //htmp[l1].ib*CNN_YSIZE*CNN_XSIZE+iring*CNN_XSIZE+ipha];
	iii++;
      }
    }
  }
#if 1
  
  free(correction);
  for (i=0;i<2;i++) {
    free(prediction[i]);
  }
  
#endif
  for (i=0;i<nbolo;i++) {
    free(rgcnn[i]);
  }
  free(rgcnn);
  
  free(l_incorr);

  if (NOMOREFITTED==0) {
    int mm;
    NOMOREFITTED=1;
    for (mm=0;mm<npixShpr;mm++)  {
      if (DOCNN[mm]==1) { //NOMOREFITTED) {
	for (i=0;i<nbolo;i++) {
	  x3[newnr[nbolo]+nbolo*(GAINSTEP)+mm*nbolo+i]=0;
	}
      }
    }
  }
 
  
  free(A0);
  free(A1);
  free(A2);
  free(B0);
  free(B1);
  free(B2);
  free(C0);
  free(C1);
  free(C2);
  free(D0);
  free(D1);
  free(D2);
  free(V00);
  free(V01);
  free(V02);
  free(V10);
  free(V11);
  free(V12);
#endif
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
  double *map = (double *) malloc(sizeof(double)*12*Nside*Nside);

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
  err=write(fp,map,12*Nside*Nside*sizeof(double));
  close(fp);
  free(map);

  return(err);
}

int getbeginfo=0;
int *allbeg;
int *allend;
int *all_realpix;
float *all_map;
int maxsize=0;

////////////////////////////////////////////////////////////////////////////////
// get hostname without forking a process to not mess with MPI
// HOSTNAME environment variable always contains the master node HOSTNAME

void GetHostname( PIOSTRING hostname) {

  FILE *f = fopen("/proc/sys/kernel/hostname", "r");
  int err=fscanf( f, "%s", hostname);
  fclose( f);
  if (err==EOF) {
    fprintf(stderr,"Problem while trying to get proc info\n");
  }
}
 
////////////////////////////////////////////////////////////////////////////////
// print free memory on each node

void PrintFreeMemOnNodes( int mpi_rank, int mpi_size, char* msg) {

  int mpi_nodes=0;
  char *temps;
  PIOSTRING hostname;
  time_t now = time( NULL);

  GetHostname( hostname);
  temps = getenv( "SLURM_NNODES"); // NERSC
  if (temps) {
    mpi_nodes = atoi( temps);
  }
  temps = getenv( "OMPI_MCA_orte_num_nodes"); // magique4
  if (temps) {
    mpi_nodes = atoi( temps);
  }
  if (mpi_nodes == 0) { // magique3?
    mpi_nodes = 64;
  }

  MPI_Barrier( MPI_COMM_WORLD);
  if (mpi_rank == 0) {
    fprintf( stderr, "\nFree memory per node %s, at %s", msg, ctime( &now));
  }

  for (int i=0; i<mpi_nodes; i++) {
    if (i == mpi_rank) {
      fprintf( stderr, "rank %02d/%d, node %02d/%02d, %s: %.1fGB free\n", mpi_rank, mpi_size-1, i, mpi_nodes-1, hostname, GetFreeMemGB());
    }
  }
  MPI_Barrier( MPI_COMM_WORLD);
}

////////////////////////////////////////////////////////////////////////////////
// get free memory in GB without forking a process to not mess with MPI

float GetFreeMemGB() {

  PIOSTRING fline;
  long memfree=0, buffers=0, cached=0;
  FILE *f = fopen("/proc/meminfo", "r");
  while (fgets( fline, PIOSTRINGMAXLEN, f) != NULL) {
    if (strstr( fline, "MemFree:")) {
      sscanf( fline, "MemFree: %ld kB", &memfree);
    }
    if (strstr( fline, "Buffers:")) {
      sscanf( fline, "Buffers: %ld kB", &buffers);
    }
    if (strstr( fline, "Cached:")) {
      sscanf( fline, "Cached: %ld kB", &cached);
    }
  }
  fclose( f);
  return( (memfree+buffers+cached) / 1024.0 / 1024.0);
}

int PIOWriteMAP(const char *path, double *value_in_double,int beg,int end)
{
  int rank,mpi_size;
  int map_size = Nside*Nside*12; //12*Nside*Nside;

  int rank_size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank_size);
  rank=rank_size;
  MPI_Comm_size(MPI_COMM_WORLD,&rank_size);
  mpi_size=rank_size;
  int k;

  ssize_t err=0;
  // Convert array from double to float for smaller file to be written
  float *value = (float *)malloc(sizeof(float)*(end-beg+1));
  for (int i = 0; i < (end-beg+1); ++i) {
    value[i] = (float)value_in_double[i];
  }
  float *map=NULL;
  if (rank==0)  {
    map=(float *)malloc(sizeof(float)*map_size);
  }

  if (getbeginfo==0) {
    getbeginfo=1;
    int i,rrk;
    if (rank==0) {
      allbeg = (int *) malloc(sizeof(int)*mpi_size);
      allend = (int *) malloc(sizeof(int)*mpi_size);
    }
    MPI_Gather(&beg,sizeof(int),MPI_BYTE,allbeg,sizeof(int),MPI_BYTE,0,MPI_COMM_WORLD);
    MPI_Gather(&end,sizeof(int),MPI_BYTE,allend,sizeof(int),MPI_BYTE,0,MPI_COMM_WORLD);

    if (rank==0) {
      for (rrk=0;rrk<mpi_size;rrk++) {
	if (maxsize<allend[rrk]-allbeg[rrk]+1) maxsize=allend[rrk]-allbeg[rrk]+1;
      }
      all_realpix = (int *) malloc(sizeof(int)*mpi_size*maxsize);
      all_map = (float *) malloc(sizeof(float)*mpi_size*maxsize);
    }

    MPI_Bcast(&maxsize,sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD);

    int *l_idx =(int *)malloc(sizeof(int)*maxsize);
    for (i=beg;i<=end;i++) l_idx[i-beg]=realpix[i-beg];

    MPI_Gather(l_idx,sizeof(int)*maxsize,MPI_BYTE,all_realpix,sizeof(int)*maxsize,MPI_BYTE,0,MPI_COMM_WORLD);

    free(l_idx);
  }
  float *l_map =(float *)malloc(sizeof(float)*maxsize);
  for (k=beg;k<=end;k++) l_map[k-beg]=value[k-beg];
  MPI_Gather(l_map,sizeof(float)*maxsize,MPI_BYTE,all_map,sizeof(float)*maxsize,MPI_BYTE,0,MPI_COMM_WORLD);
  free(l_map);


  if (rank==0) {
    int i,rrk;
    for (rrk=0;rrk<mpi_size;rrk++) {
      int l_beg,l_end;
      l_beg=allbeg[rrk];
      l_end=allend[rrk];
      for (i=l_beg;i<=l_end;i++) map[all_realpix[i-l_beg+rrk*maxsize]]=all_map[i-l_beg+rrk*maxsize];
    }
  }

  if (rank==0) {
    PIOSTRING fitspath;
    sprintf( fitspath, "%s.fits", path);
    if (remove( fitspath) == 0) {
      fprintf( stderr, "removed existing MAP: %s\n", fitspath);
    }
    fprintf( stderr, "WRITE MAP: %s\n", fitspath);
    write_healpix_map( map, Nside, fitspath, 0, "G");
    if (map!=NULL)
      free(map);
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

/* ---------------------------------------------------------------------------------*/
PyObject * init_PyFunction(char* path,char *funcname){
    /* Load python parameters form file gived in parameters path */    

    // init variables
    PyObject *pName, *pModule, *pDict, *pFunc;

    // Initialize the Python Interpreter
    Py_Initialize();

    strip_ext(path);
    
    //add current path to python
    PyRun_SimpleString("import sys");
    PyRun_SimpleString("import os");
    PyRun_SimpleString("sys.path.append(os.getcwd())");

    // Build the name object
    pName = PyUnicode_FromString(path);
    
    // Load the module object
    pModule = PyImport_Import(pName);

    if ((pModule = PyImport_Import(pName)) == NULL) {
      fprintf(stderr,"Error: PyImport_Import\n");
      exit(0);
    }
 
    // pDict is a borrowed reference 
    pDict = PyModule_GetDict(pModule);

    // pFunc is also a borrowed reference 
    pFunc = PyDict_GetItemString(pDict, funcname);

    if (!PyCallable_Check(pFunc)){
        PyErr_Print();
	MPI_Finalize();
        exit(0);
    }
    
    return pFunc;
}
      
/* ---------------------------------------------------------------------------------*/
long exec_PyFunction(PyObject *pFunc){
  PyObject * result;
  PyObject *pArgs;
    if (PyCallable_Check(pFunc)){
      pArgs = PyTuple_Pack(1, PyUnicode_FromString("nom_du_fichier"));
      result=PyObject_CallObject(pFunc,pArgs);
        
    } else {
        PyErr_Print();
	MPI_Finalize();
        exit(0);
    }
  
    Py_DECREF(pArgs);

    return PyLong_AsLong(result);
}
// Fonction pour copier les données d'un tableau Python d'entiers à un tableau C
int copy_int_array(PyObject *pArray, PIOINT *array) {
  
  if (!PyList_Check(pArray) && !PyTuple_Check(pArray)) {
    fprintf(stderr, "L'objet n'est ni une liste ni un tuple\n");
    exit(0);
    return -1;
  }

  Py_ssize_t n = PyList_Check(pArray) ? PyList_Size(pArray) : PyTuple_Size(pArray);

  for (Py_ssize_t i = 0; i < n; i++) {
    PyObject *item = PyList_Check(pArray) ? PyList_GetItem(pArray, i) : PyTuple_GetItem(pArray, i);
    if (!PyLong_Check(item)) {
      fprintf(stderr, "L'élément n'est pas un entier\n");
      return -1;
    }
    array[i] = (PIOINT) PyLong_AsLong(item);
  }

  return (int) n;
}
// Fonction pour copier les données d'un tableau Python d'entiers à un tableau C
int copy_float_array(PyObject *pArray, PIOFLOAT *array) {
  
  if (!PyList_Check(pArray) && !PyTuple_Check(pArray)) {
    fprintf(stderr, "The provided object is not a list or a tuple\n");
    exit(0);
    return -1;
  }

  Py_ssize_t n = PyList_Check(pArray) ? PyList_Size(pArray) : PyTuple_Size(pArray);

  for (Py_ssize_t i = 0; i < n; i++) {
    PyObject *item = PyList_Check(pArray) ? PyList_GetItem(pArray, i) : PyTuple_GetItem(pArray, i);
    if (!PyFloat_Check(item)) {
      fprintf(stderr, "L'élément n'est pas un double\n");
      free(array);
      return -1;
    }
    array[i] = (PIOFLOAT) PyFloat_AsDouble(item);
  }

  return (int) n;
}

int Get_NumberOfChannels(PyObject *projFunc)
{
  PyObject *pValue=NULL;
  int n=0;
  
  if (projFunc != NULL){
    pValue = PyObject_CallMethod(projFunc,"getnumber_of_channels",NULL);

    // Traitement de la valeur de retour
    if (pValue != NULL) {
      // Assurer que pValue est un tuple
      n=(int) PyLong_AsLong(pValue);
      Py_DECREF(pValue);
    }
    else {
      PyErr_Print();
      MPI_Finalize();
      exit(0);
    }
    
  }
  else {
    PyErr_Print();
    MPI_Finalize();
    exit(0);
  }
  
  return n;
}


int Get_NumberOfDiag(PyObject *diagFunc)
{
  PyObject *pValue=NULL;
  int n=0;
  
  if (diagFunc != NULL){
    pValue = PyObject_CallMethod(diagFunc,"getnumber_of_index",NULL);

    // Traitement de la valeur de retour
    if (pValue != NULL) {
      // Assurer que pValue est un tuple
      n=(int) PyLong_AsLong(pValue);
      Py_DECREF(pValue);
    }
    else {
      PyErr_Print();
      MPI_Finalize();
      exit(0);
    }
    
  }
  else {
    PyErr_Print();
    MPI_Finalize();
    exit(0);
  }
  
  return n;
}

double init_channels(hpix * h,PyObject *projFunc,double psi,PIOFLOAT *External,double rgnorm,int ipix,int idx_bolo,double hit,int idx_in_ring,PIOFLOAT *sig)
{
  
  PyObject *pValue=NULL;
  double ohit=0;
  
  if (projFunc != NULL){
    PyObject *pArg1 = PyFloat_FromDouble((double)psi); // val est un flottant
    PyObject *pList = PyList_New(NB_EXTERNAL);
    
    for (int i = 0; i < NB_EXTERNAL; i++) {
      PyList_SetItem(pList, i, PyFloat_FromDouble((double) External[i*RINGSIZE]));
    }

    PyObject *pArg3 = PyFloat_FromDouble((double)rgnorm); // val est un flottant
    PyObject *pArg4 = PyLong_FromLong((long) ipix); // val est un flottant
    PyObject *pArg5 = PyLong_FromLong((long)idx_bolo); // val est un flottant
    PyObject *pArg6 = PyFloat_FromDouble((double) hit); // val est un flottant
    PyObject *pArg7 = PyLong_FromLong((long)idx_in_ring); // val est un flottant
    PyObject *pArg8 = PyFloat_FromDouble((double) *sig); // val est un flottant
    
    pValue = PyObject_CallMethod(projFunc, "eval", "(OOOOOOOO)",
				 pArg1, pList, pArg3,
				 pArg4, pArg5, pArg6,
				 pArg7, pArg8);

    // Traitement de la valeur de retour
    if (pValue != NULL) {
      if (PyTuple_Check(pValue) && PyTuple_Size(pValue) == 3) {
	// Extraire les deux tableaux du tuple
	*sig = (double) PyFloat_AsDouble(PyTuple_GetItem(pValue, 0));
	ohit = (double) PyFloat_AsDouble(PyTuple_GetItem(pValue, 1));
	int err=copy_float_array(PyTuple_GetItem(pValue, 2),h->channels);
	if (err!=MAXCHANNELS) {
	  fprintf(stderr, "Projection function does not provide the good number of channels: expected %d  get %d\n",MAXCHANNELS,err);
	}
      }
      else {
	  fprintf(stderr, "Projection function does not provide the good number of argument (3): expected sig,hit,chan_value\n");
	  fprintf(stderr, " Where sig and hit are a double (respectively signal and hit count) and\n");
	  fprintf(stderr, " chan_value is a list of double values with len(chan_value)=%d",MAXCHANNELS);
	  exit(0);
      }
      Py_DECREF(pValue);
    }
    else {
      fprintf(stderr, "Problem while trying to compute the projection\n");
      PyErr_Print();
      MPI_Finalize();  
      exit(0);
    }

    // Nettoyage des arguments
    Py_DECREF(pArg1);
    Py_DECREF(pList);
    Py_DECREF(pArg3);
    Py_DECREF(pArg4);
    Py_DECREF(pArg5);
    Py_DECREF(pArg6);
    Py_DECREF(pArg7);
    Py_DECREF(pArg8);
    
  } else {
    PyErr_Print();
    MPI_Finalize(); 
    exit(0);
  }

    return ohit;
  
}

int calc_sparse_hpr(PyObject *sparseFunc,
		    long rg,
		    long ib,
		    long hpix,
		    double psi,
		    PIOFLOAT *External,
		    PIOFLOAT *oval,
		    PIOINT *oval_idx)
{
  PyObject *pValue=NULL;
  int n=-1;

  if (sparseFunc != NULL){
    // Créer des arguments pour la fonction Python
    // Création des arguments
    PyObject *pArg1 = PyLong_FromLong((long)rg); // rg est un entier
    PyObject *pArg2 = PyLong_FromLong((long)ib); // ib est un entier
    PyObject *pArg3 = PyLong_FromLong((long)hpix); // idx est un entier
    PyObject *pArg4 = PyFloat_FromDouble((double)psi); // val est un flottant
    PyObject *pList = PyList_New(NB_EXTERNAL);
    
    for (int i = 0; i < NB_EXTERNAL; i++) {
      PyList_SetItem(pList, i, PyFloat_FromDouble((double) External[i*RINGSIZE]));
    }

    // Appel de la méthode 'eval's
    pValue = PyObject_CallMethod(sparseFunc, "eval", "(OOOOO)", pArg1, pArg2, pArg3, pArg4, pList);

    // Traitement de la valeur de retour
    if (pValue != NULL) {
      // Assurer que pValue est un tuple
      if (PyTuple_Check(pValue) && PyTuple_Size(pValue) == 2) {
	// Extraire les deux tableaux du tuple
	int n1=copy_int_array(PyTuple_GetItem(pValue, 0), oval_idx);
	int n2=copy_float_array(PyTuple_GetItem(pValue, 1), oval);
	if (n1!=n2) {
	  fprintf(stderr, "Sparse function should provide an equal number of invex and value, here Sroll received %d %d\n",n1,n2);
	}
	n=n1;
	Py_DECREF(pValue);
      }
      else {
	fprintf(stderr,"Problem while executing the method get value inside SparseFunc class\n");
	PyErr_Print();
	MPI_Finalize(); 
	exit(0);
      }
    }

    // Nettoyage des arguments
    Py_DECREF(pArg1);
    Py_DECREF(pArg2);
    Py_DECREF(pArg3);
    Py_DECREF(pArg4);
    Py_DECREF(pList);
    
  } else {
    PyErr_Print();
    MPI_Finalize(); 
    exit(0);
  }
  
  return n;
}


int calc_diag_hpr(PyObject *diagFunc,
		  long rg,
		  long ib,
		  long hpix,
		  double psi,
		  PIOFLOAT *External)
{
  PyObject *pValue=NULL;
  int n=-1;

  if (diagFunc != NULL){
    // Créer des arguments pour la fonction Python
    // Création des arguments
    PyObject *pArg1 = PyLong_FromLong((long)rg); // rg est un entier
    PyObject *pArg2 = PyLong_FromLong((long)ib); // ib est un entier
    PyObject *pArg3 = PyLong_FromLong((long)hpix); // idx est un entier
    PyObject *pArg4 = PyFloat_FromDouble((double)psi); // val est un flottant
    PyObject *pList = PyList_New(NB_EXTERNAL);
    
    for (int i = 0; i < NB_EXTERNAL; i++) {
      PyList_SetItem(pList, i, PyFloat_FromDouble((double) External[i*RINGSIZE]));
    }

    // Appel de la méthode 'eval's
    pValue = PyObject_CallMethod(diagFunc, "get_diag_idx", "(OOOOO)", pArg1, pArg2, pArg3, pArg4, pList);

    
    // Traitement de la valeur de retour
    if (pValue != NULL) {
      n=(int) PyLong_AsLong(pValue);
      Py_DECREF(pValue);
    }
    else {
      fprintf(stderr,"%d\n",(int) rg);
      fprintf(stderr,"%d\n",(int) ib);
      fprintf(stderr,"%d\n",(int) hpix);
      fprintf(stderr,"Problem while executing the method get_diag_idx inside DiagFunc class\n");
      PyErr_Print();
      MPI_Finalize(); 
      exit(0);
    }

    // Nettoyage des arguments
    Py_DECREF(pArg1);
    Py_DECREF(pArg2);
    Py_DECREF(pArg3);
    Py_DECREF(pArg4);
    Py_DECREF(pList);
    
  } else {
    fprintf(stderr,"no DiagFunc class defined\n");
    PyErr_Print();
    MPI_Finalize(); 
    exit(0);
  }
  
  return n;
}

/* ---------------------------------------------------------------------------------*/
void ClosepFunc(PyObject *pFunc){
  Py_XDECREF(pFunc);
}

///////////////////////////////////////////////////////////////////////////////////////////////
/////////                                 MAIN                                         ////////
///////////////////////////////////////////////////////////////////////////////////////////////

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
    char * err=realpath( argv[0], rpath);
    if (err) 
      fprintf( stderr, "%s: --------------------------\n", __FILE__ );
    now = time( NULL);
    fprintf( stderr, "%s: --------------------------\n", __FILE__ );
    fprintf( stderr, "%s: Starting %s\n",                __FILE__, rpath);
    fprintf( stderr, "%s: at %s",                        __FILE__, ctime( &now));
    fprintf( stderr, "%s: with %d MPI ranks\n", __FILE__, mpi_size); //, omp_get_max_threads());
    fprintf( stderr, "%s: --------------------------\n", __FILE__ );
#if 0
    /* Ensure mpi_rank is a power of two (required by some algo...) */
    if (!isPowerOfTwo(mpi_size)) {
      fprintf(stderr,"INVALID mpi_size(%d), shall be a power of two!\n", mpi_size);
      return 1;
    }
#endif
    fprintf( stderr, "\nstarting troll with MAXSIMU=%d, MAXEXTERNALHPR=%d\n", MAXSIMU, MAXEXTERNALHPR);
  }

  //------ read params ----
  troll_parContent par;
  GetHostname( hostname);

  PyObject *pyParam = troll_readParam(&par, argv[1] );
  if (pyParam==NULL) {
    fprintf(stderr, "Unable to parse the parameter file.\n");
    exit(-1);
  }

  Param = &par;

  if (Param->flag_verbose==_PAR_TRUE) verbose=Param->verbose;
  if (verbose==1) fprintf(stderr,"%s %d %d\n",__FILE__,__LINE__,rank);


  if (Param->flag_do_offset==_PAR_TRUE) do_offset=Param->do_offset;
  
  PyObject *projFunc=NULL;
  if (Param->flag_projection==_PAR_TRUE) {
    PyObject *pArgs = NULL;
    PyObject *pClass = NULL;
    PyObject *pName = PyUnicode_FromString(argv[1]);
    PyObject *pModule = PyImport_Import(pName);
    Py_DECREF(pName);

    if (pModule == NULL) {
      PyErr_Print();
      fprintf(stderr, "Échec du chargement du module\n");
      return 1;
    }

    pClass = PyObject_GetAttrString(pModule, Param->projection);
    if (pClass == NULL || !PyCallable_Check(pClass)) {
      PyErr_Print();
      fprintf(stderr, "Échec de la récupération de la classe\n");
      Py_XDECREF(pClass);
      Py_DECREF(pModule);
      return 1;
    }
    // Créer un tuple pour les arguments du constructeur
    pArgs = PyTuple_New(2);
    PyObject *pyRank = PyLong_FromLong((long) rank); // val est un flottant
    
    PyTuple_SetItem(pArgs, 0, pyParam);  // Le tuple prend la propriété de 'param'
    PyTuple_SetItem(pArgs, 1, pyRank);  // Provide the mpi rank 
    projFunc = PyObject_CallObject(pClass,pArgs); 

    if (projFunc == NULL) {
      PyErr_Print();
      fprintf(stderr, "Échec de la création de l'instance de la classe projection\n");
      Py_DECREF(pModule);
      return 1;
    }
  }
  
  if (verbose==1) fprintf(stderr,"%s %d %d\n",__FILE__,__LINE__,rank);
  
  if (projFunc==NULL) {
    MAXCHANNELS = 1;
  }
  else {
    MAXCHANNELS = Get_NumberOfChannels(projFunc);
    
    if (MAXCHANNELS>MAXCHAN) {
      fprintf(stderr,"The number of channels from python class %s is %d. Troll has been compiled with MAXCHAN=%d\n",Param->projection,
	      (int) MAXCHANNELS,(int) MAXCHAN);
      exit(0);
    }
  }
  
  if (verbose==1) fprintf(stderr,"%s %d %d\n",__FILE__,__LINE__,rank);

  if (rank==0)  fprintf(stderr,"Projection uses %d channels\n",MAXCHANNELS);

  int nside128=128; // IF DONside=1 then all 128 maps are used in Nside
  if (Param->flag_TEMPLATE_NSIDE) {
    nside128=Param->TEMPLATE_NSIDE;
  }

  GAINSTEP = Param->GAINSTEP;

  if (rank==0) fprintf(stderr,"RINGSIZE = %d \n",(int) (RINGSIZE));
  
  /*-------------------------------------------------------------------------*/
  /*   SAVE PARAMETER FILE                                                   */
  /*-------------------------------------------------------------------------*/
  if (rank==0) {
    char commandtest[PIOSTRINGMAXLEN*64];
    sprintf(commandtest,"cp %s.py %s.py",argv[1],Param->Out_VEC[0]);
    int err=system(commandtest);
    if (err) {
      fprintf(stderr,"Error while copy the parameters\n");
    }
  }
  
  if (verbose==1) fprintf(stderr,"%s %d %d\n",__FILE__,__LINE__,rank);
  
  NUMBEROFITER = Param->N_IN_ITT;

  /*-------------------------------------------------------------------------*/
  /* parameters consistency checks                                           */
  /* and default values / legacy behavior for optional parameters            */
  /*-------------------------------------------------------------------------*/

  // global number of bolometers
  nbolo = Param->n_Ptg;
  
  assert( Param->n_OUT_NOPOL == Param->n_Out_MAP);
  assert( Param->n_Sub_HPR == Param->n_SUB_HPRCOEF);

  if (verbose==1) fprintf(stderr,"%s %d %d\n",__FILE__,__LINE__,rank);
  
  // get CNN parameters if needed
  memset(DOCNN,0,MAXEXTERNALHPR*sizeof(int));
  
  if (Param->flag_DOCNN==_PAR_TRUE) {
    if (Param->n_DOCNN>MAXEXTERNALHPR) {
      fprintf(stderr,"To big DOCNN table [%d], should be smaller than [%d] as the max number of TF\n",(int)(Param->n_DOCNN),(int)(MAXEXTERNALHPR));
      exit(0);
    }
    COMP_CNN=1;
    for (i=0;i<Param->n_DOCNN;i++) {
      DOCNN[i]=Param->DOCNN[i];
      if (DOCNN[i]==1) COMP_CNN+=2;
    }
    
    int num_tensorflow=mpi_size;
    if (Param->flag_CNN_CORE==_PAR_TRUE) num_tensorflow=Param->CNN_CORE;
    
    /* Determine my color for tensorflow serveur*/ 
    int color = rank % num_tensorflow;
    
    if (verbose==1) fprintf(stderr,"%s %d %d\n",__FILE__,__LINE__,rank);
         
    /* Split the intercommunicator */ 
    MPI_Comm_split (MPI_COMM_WORLD, color, rank, &python_comm );

    MPI_Comm_rank(python_comm, &python_rank );
    MPI_Comm_size(python_comm, &mpi_python_size );

    int ntensorflow=0;
    for (i=0;i<mpi_size;i++) if ((i/num_tensorflow)==0) ntensorflow++;
    int *tabtens = (int *) malloc(sizeof(int)*ntensorflow);
    ntensorflow=0;
    for (i=0;i<mpi_size;i++) {
      if ((i/num_tensorflow)==0) {
	      tabtens[ntensorflow]=i;
	      //if (rank==0) fprintf(stderr,"%d %d\n",ntensorflow,tabtens[ntensorflow]);
	      ntensorflow++;
      }
    }
    if (rank==0) fprintf(stderr,"USE %d tensoflow rank\n",ntensorflow);
	
    MPI_Group world_group;
    MPI_Comm_group(MPI_COMM_WORLD, &world_group);
    MPI_Group prime_group;
    MPI_Group_incl(world_group,ntensorflow,tabtens,&prime_group );

    // Create a new communicator based on the group
    MPI_Comm_create_group(MPI_COMM_WORLD, prime_group, 0, &tensorflow_comm);
    
    free(tabtens);
    
    if (python_rank==0) {
      MPI_Comm_rank(tensorflow_comm, &tensorflow_rank ); 
      MPI_Comm_size(tensorflow_comm, &mpi_tensorflow_size );
    }

    /*
      TO DEBUG
    if (python_rank==0) {
      InitPython(&MyPythonBackend,Param->CALLFOSCAT,rank);
    }
    */
    
  }


  // default don't BUILDTF
  if (Param->flag_BUILDTF == 0) {
    Param->BUILDTF = 0;
    Param->flag_BUILDTF = 1;
  }

  if (Param->n_bolomask == 0) {
    // if bolomask is empty, produce a map with all bolometers
    assert( Param->n_Out_MAP == 1);
    Param->n_bolomask = nbolo;
    Param->bolomask = malloc( nbolo * sizeof( PIOINT));
    assert( Param->bolomask != NULL);
    for (i=0; i<nbolo; i++) {
      Param->bolomask[i] = 1;
    }
  }

  if (verbose==1) fprintf(stderr,"%s %d %d\n",__FILE__,__LINE__,rank);
  
  if (Param->flag_MAPRINGS == 0) {
    Param->flag_MAPRINGS = 1;
    Param->n_MAPRINGS = Param->n_Out_MAP;
    Param->MAPRINGS = malloc( Param->n_MAPRINGS * sizeof( PIOLONG));
    assert( Param->MAPRINGS != NULL);
    for (i=0; i<Param->n_MAPRINGS; i++) {
      Param->MAPRINGS[i] = 1;
    }
  }
  assert( Param->n_MAPRINGS == Param->n_Out_MAP);

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

  if (verbose==1) fprintf(stderr,"%s %d %d\n",__FILE__,__LINE__,rank);
  
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

  if (verbose==1) fprintf(stderr,"%s %d %d\n",__FILE__,__LINE__,rank);
  
  /* Display ring dispatching between ranks */
 //if (Param->verbose > 0) {
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
 // }

  /*--- End : MPI: Ring dispatching --- */


  //PIOINT Nside=Param->Nside;
  Nside = Param->Nside;
  if (Param->flag_stim_paramfiles == 1) {
    assert( Param->n_stim_paramfiles == nbolo);
  }

#ifndef DOMAP
  PIOLONG LMAX=Param->LMAX*(Param->LMAX+1)+Param->LMAX;
#endif
  MPI_Barrier(MPI_COMM_WORLD);

  PIOSTRING *mapout = (PIOSTRING *) malloc(sizeof(PIOSTRING)*Param->n_Out_MAP);

  if (verbose==1) fprintf(stderr,"%s %d %d\n",__FILE__,__LINE__,rank);
  
  NORM_GAIN = Param->NORM_GAIN;
  REMOVE_CAL = Param->REMOVE_CAL;

  assert( Param->n_NEP == nbolo);
  assert( Param->n_Calibration == nbolo);
  NEP_tab = malloc( nbolo * sizeof( double));

  double avvnep=0;
  for (i=0;i<nbolo;i++) avvnep+=Param->NEP[i]/Param->Calibration[i];
  for (i=0;i<nbolo;i++) NEP_tab[i]=nbolo*Param->NEP[i]/Param->Calibration[i]/avvnep;
  if (rank==0) {
    for (i=0;i<nbolo;i++) fprintf( stderr,"NEP_tab[%ld]=%lf\n", i, NEP_tab[i]);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  
  PIOINT  *badring;
  PIOINT  *ibadring;
  int rg_max;
#ifdef GAIN_RATIO
  PIODOUBLE  **gain_ratio = (PIODOUBLE **) malloc(nbolo*sizeof(PIODOUBLE *));
  PIOLONG     **gain_ratio_off = (PIOLONG **) malloc(nbolo*sizeof(PIOLONG *));
#endif

  srand48((long)(rank+Param->SEED[0]));
  
  PIOINT *l_nhpix  = (PIOINT *) malloc(sizeof(PIOINT)*12*Nside*Nside);
  hpix **l_hpix = (hpix **) malloc(sizeof(hpix *)*12*Nside*Nside);
  memset(l_nhpix,0,sizeof(PIOINT)*12*Nside*Nside);
#if 0
  assert( Param->n_External % nbolo == 0);
  npixShpr = Param->n_External / nbolo;
  assert( npixShpr <= MAXEXTERNALHPR);
#endif
  
  assert( Param->n_External_MAP % MAXCHANNELS == 0);
  npixmap = Param->n_External_MAP/MAXCHANNELS;
  assert( npixmap <= MAXEXTERNALMAP); 

  /*======================================================
    =
    =      read data
    =
    =*/

  long vmem,phymem;
  
  PIOLONG ib;
  eta=(double *) malloc(sizeof(double)*nbolo); // (1-crosspol)/(1+crosspol)
  eta_dest=(double *) malloc(sizeof(double)*nbolo);

  if (rank==0) fprintf(stderr,"Avv GAIN is equal to 0 if ==1 : NORM_GAIN : %d\n",(int) Param->NORM_GAIN);

  PIOFLOAT **skymodel = (PIOFLOAT **) malloc(sizeof(PIOFLOAT *)*MAXCHANNELS);
  for (i=0;i<MAXCHANNELS;i++) skymodel[i]=(PIOFLOAT *) malloc(sizeof(PIOFLOAT)*12*Nside*Nside);

  PIOFLOAT **mapmodel = NULL;

  if (verbose==1) fprintf(stderr,"%s %d %d\n",__FILE__,__LINE__,rank);
  
  if (Param->flag_External_MAP==_PAR_TRUE) {

    mapmodel = (PIOFLOAT **) malloc(sizeof(PIOFLOAT)*Param->n_External_MAP);
    
    for (int l=0;l<Param->n_External_MAP;l++) {
      mapmodel[l]= (PIOFLOAT *) malloc(sizeof(PIOFLOAT)*12*nside128*nside128);

      PIOLONG nsa  = noDMC_readObject_PIOFLOAT(Param->External_MAP[l],0,12*nside128*nside128,mapmodel[l]);
      if (nsa<0) {
	      fprintf(stderr, "Impossible to the model map: %s %d\n",Param->External_MAP[l],(int) nsa);

	      exit ( -1);
      }
    }
  }
  
  PIOFLOAT *addpol=NULL;
  if (Param->flag_ADDPOL==_PAR_TRUE) {
    addpol= (PIOFLOAT *) malloc(sizeof(PIOFLOAT)*12*128*128*2);
    PIOLONG nsa  = noDMC_readObject_PIOFLOAT(Param->ADDPOL,0,12*128*128*2,addpol);
    if (nsa<0) {
      fprintf(stderr, "Impossible to read ADDPOL: %s %d\n",Param->ADDPOL,(int) nsa);

      exit ( -1);
    }
  }
 
  PyObject *sparseFunc=NULL;
  PyObject *diagFunc=NULL;
  PyObject *pArgs = NULL;
  PyObject *pClass = NULL;
  
  if (verbose==1) fprintf(stderr,"%s %d %d\n",__FILE__,__LINE__,rank);
  
  npixShpr=0;
  if (Param->flag_SparseFunc==_PAR_TRUE) {
    PyObject *pName = PyUnicode_FromString(argv[1]);
    PyObject *pModule = PyImport_Import(pName);
    Py_DECREF(pName);

    if (pModule == NULL) {
      PyErr_Print();
      fprintf(stderr, "Problem while creating the sparseFunc from class %s\n",argv[1]);
      return 1;
    }

    pClass = PyObject_GetAttrString(pModule, Param->SparseFunc);
    if (pClass == NULL || !PyCallable_Check(pClass)) {
      PyErr_Print();
      fprintf(stderr, "Problem while loading the classe SparseFunc\n");
      Py_XDECREF(pClass);
      Py_DECREF(pModule);
      return 1;
    }
    // Créer un tuple pour les arguments du constructeur
    pArgs = PyTuple_New(2);
    PyObject *pyRank = PyLong_FromLong((long) rank); // val est un flottant
    PyTuple_SetItem(pArgs, 0, pyParam);  // Le tuple prend la propriété de 'param'
    PyTuple_SetItem(pArgs, 1, pyRank);  // Provide the mpi rank 
    sparseFunc = PyObject_CallObject(pClass,pArgs); 

    if (sparseFunc == NULL) {
      PyErr_Print();
      fprintf(stderr, "Problem while creating the class instance\n");
      Py_DECREF(pModule);
      return 1;
    }
  }

  if (verbose==1) fprintf(stderr,"%s %d %d\n",__FILE__,__LINE__,rank);
  
  int nb_diag=0;
  if (Param->flag_DiagFunc==_PAR_TRUE) {
    PyObject *pName = PyUnicode_FromString(argv[1]);
    PyObject *pModule = PyImport_Import(pName);
    Py_DECREF(pName);

    if (pModule == NULL) {
      PyErr_Print();
      fprintf(stderr, "Problem while creating the DiagFunc from class %s\n",argv[1]);
      return 1;
    }

    pClass = PyObject_GetAttrString(pModule, Param->DiagFunc);
    if (pClass == NULL || !PyCallable_Check(pClass)) {
      PyErr_Print();
      fprintf(stderr, "Problem while loading the DiagFunc class\n");
      Py_XDECREF(pClass);
      Py_DECREF(pModule);
      return 1;
    }
    // Créer un tuple pour les arguments du constructeur
    pArgs = PyTuple_New(2);
    PyObject *pyRank = PyLong_FromLong((long) rank); // val est un flottant
    PyTuple_SetItem(pArgs, 0, pyParam);  // Le tuple prend la propriété de 'param'
    PyTuple_SetItem(pArgs, 1, pyRank);  // Provide the mpi rank 
    diagFunc = PyObject_CallObject(pClass,pArgs); 

    if (sparseFunc == NULL) {
      PyErr_Print();
      fprintf(stderr, "Problem while creating the class instance\n");
      Py_DECREF(pModule);
      return 1;
    }
    nb_diag=Get_NumberOfDiag(diagFunc);
  }
  
  int number_of_iterations = 1;

  if (verbose==1) fprintf(stderr,"%s %d %d\n",__FILE__,__LINE__,rank);
  
  NB_EXTERNAL=Param->n_External/nbolo;
  if (NB_EXTERNAL==0) NB_EXTERNAL=1;
  
  for (ib=0;ib<nbolo;ib++) {
    if (Param->n_in_template_map>0) {
      for (int nchan=0;nchan<MAXCHANNELS;nchan++) {
	      PIOLONG nsa  = noDMC_readObject_PIOFLOAT(Param->in_template_map[ib*MAXCHANNELS+nchan],0,12*Nside*Nside,
	      skymodel[nchan]);  
	      if (nsa<0) {
	        fprintf(stderr, "Impossible to read in_template_map: %s %d\n",Param->in_template_map[ib],(int) nsa);
	        exit ( -1);
	      }
      }
    }
    else {
      for (int nchan=0;nchan<MAXCHANNELS;nchan++) {
	      memset(skymodel[nchan],0,12*Nside*Nside*sizeof(PIOFLOAT));
      }
    }

    eta[ib]=(1.-Param->CrossPol[ib])/(1.+Param->CrossPol[ib]);
    if (Param->D_NOPOL) eta_dest[ib]=0;
    else eta_dest[ib]=eta[ib];

    double sxi= (Param->Calibration[ib]/Param->NEP[ib])*(Param->Calibration[ib]/Param->NEP[ib]);
    if (rank==0) fprintf(stderr,"SXI %d %lg\n",(int) ib,sxi);
    sprintf(Command,"begin=%lld;end=%lld",
            (long long) (globalRankInfo.BeginRing[rank]),
            (long long) (globalRankInfo.EndRing[rank]));

    int ring_count = globalEndRing-globalBeginRing+1;
    badring = (PIOINT *) malloc(sizeof(PIOINT)*ring_count);
    if (badring == NULL) {
      fprintf(stderr, "**** ERROR **** Memory allocation error !!!!\n");
      fprintf(stderr, "** rank = %d  globalRankInfo.BeginRing[rank]="PIOLONG_FMT"  globalRankInfo.EndRing[rank]="PIOLONG_FMT"\n", rank, globalRankInfo.BeginRing[rank], globalRankInfo.EndRing[rank]);
    }
    memset(badring,0,sizeof(PIOINT)*ring_count);
    int tperr;
    int bad_rings;
    if (Param->n_Badring>0) {
      tperr = noDMC_readObject_PIOINT(Param->Badring[ib],globalBeginRing,ring_count,badring);
      if (tperr<0) {
	fprintf(stderr, "Impossible to read Badring[%ld]: %s %d\n",ib,Param->Badring[ib],tperr);
	exit ( -1);
      }
    }

    ibadring=(PIOINT *) malloc(sizeof(PIOINT)*ring_count);
    // compute ring taking into account the badring
    {
      rg_max=0;
      bad_rings = 0;
      for (i=0;i<ring_count;i++) {
	      ibadring[i]=rg_max;
	      if (badring[i]==0) {
	        rg_max++;
	      }else{
          bad_rings++;
        }
      }
    }
    
    if (rank==0) fprintf(stderr,"ring_count = %d\n",ring_count);
    if (rank==0) fprintf(stderr,"bad_rings = %d\n",bad_rings);
    
    if (rank==0) fprintf(stderr,"RG_MAX %ld\n",(long) rg_max);
    if (rank==0) fprintf(stderr,"NB_DIAG %ld\n",(long) nb_diag);
    /*=========================================================================================
      Compute spline in Time:
      =========================================================================================*/

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

      if (rank==0) {
        GetProcMem(&vmem,&phymem);
        fprintf(stderr,"\nafter stim bolo loop: used VMEM %.1lf[PHYS %.1lf]MB\n",
            (double) vmem/1024./1024., (double) phymem/1024./1024.);
      }

    }

    if (verbose==1) fprintf(stderr,"%s %d %d\n",__FILE__,__LINE__,rank);
    
    for (rg=globalRankInfo.BeginRing[rank];rg<=globalRankInfo.EndRing[rank];rg+=Param->RSTEP) {

      if (badring[rg-globalBeginRing]==0) {

        PIOFLOAT *h,*Sub_HPR;
        PIODOUBLE *ph,*th,*psi;
        PIOFLOAT *hprtab_cal;
        PIOBYTE surv=-1;
	
        if (rg>=0&&rg<=globalEndRing/2)    surv=1;
        if (rg>=globalEndRing/2)  surv=12;
	

        PIOFLOAT *y[MAXSIMU];
        // get rid of gcc "maybe-uninitialized" warning depending on optimsation level
        for (iter = 0; iter < number_of_iterations; iter++) {
          y[iter] = NULL;
        }

        sprintf(Command,"ring=%lld",(long long)(rg));
        h = (PIOFLOAT *) malloc(sizeof(PIOFLOAT)*RINGSIZE);
        assert( (h != NULL) && "Not enough memory to allocate hitcount HPR array");
        tperr = noDMC_readObject_PIOFLOAT(Param->Hit[ib],rg*RINGSIZE,RINGSIZE,h);
        if (tperr<0) {
          fprintf(stderr, "Impossible to read Hit[%ld]: %s %d\n",ib,Param->Hit[ib],tperr);          exit ( -1);
        }

        if (Param->flag_Sub_HPR == _PAR_TRUE) {
          Sub_HPR = (PIOFLOAT *) malloc(sizeof(PIOFLOAT)*RINGSIZE);
          assert( (Sub_HPR != NULL) && "Not enough memory to allocate SUB_HPR HPR array");
          tperr = noDMC_readObject_PIOFLOAT(Param->Sub_HPR[ib],rg*RINGSIZE,RINGSIZE,Sub_HPR);
          if (tperr<0) {
            fprintf(stderr, "Impossible to read Sub_HPR[%ld]: %s %d\n",ib,Param->Sub_HPR[ib],tperr);
            exit ( -1);
          }
        } else {
          Sub_HPR = NULL;
        }

        if (Param->flag_stim_paramfiles == 0) {
          y[0] = (PIOFLOAT *) malloc(sizeof(PIOFLOAT)*RINGSIZE);
          assert( (y[0] != NULL) && "Not enough memory to allocate signal HPR array");
          tperr = noDMC_readObject_PIOFLOAT(Param->Signal[ib],rg*RINGSIZE,RINGSIZE,y[0]);
          if (tperr<0) {
            fprintf(stderr, "Impossible to read Signal[%ld]: %s %d\n",ib,Param->Signal[ib],tperr);
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

        hprtab_cal = (PIOFLOAT *) malloc(sizeof(PIOFLOAT)*RINGSIZE);
        assert( (hprtab_cal != NULL) && "Not enough memory to allocate total dipole HPR array");
	if (GAINSTEP>0) {
	  tperr = noDMC_readObject_PIOFLOAT(Param->HPR_Calib[ib],rg*RINGSIZE,RINGSIZE,hprtab_cal);
	  if (tperr<0) {
	    fprintf(stderr, "Impossible to read HPR_Calib[%ld]: %s %d\n",ib,Param->HPR_Calib[ib],tperr);
	    exit ( -1);
	  }
	}
	else {
	  memset(hprtab_cal,0,sizeof(PIOFLOAT)*RINGSIZE);
	}

	for (i=0; i<RINGSIZE; i++) {
	  y[0][i] *= Param->Calibration[ib];
	}

        if (Param->ADDDIP == 1) {
          // for projection only, signal must contain the total dipole, add it if it's not here
          for (i=0; i<RINGSIZE; i++) {
            y[0][i] = y[0][i] + hprtab_cal[i] * Param->Calibration[ib];
          }
        }

        if (Param->flag_addHPR_name == 1) {
          PIOFLOAT *addhpr = (PIOFLOAT *) malloc( sizeof( PIOFLOAT) * RINGSIZE);
          for (int i = 0; i < Param->n_addHPR_name / nbolo; i++) {
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

        PIOFLOAT *External;
        if (Param->flag_External==_PAR_TRUE) {
          External = (PIOFLOAT *) _PIOMALLOC(sizeof(PIOFLOAT)*RINGSIZE*NB_EXTERNAL);
	  for (i=0;i<NB_EXTERNAL;i++) {
	    tperr = noDMC_readObject_PIOFLOAT(Param->External[ib+i*nbolo],rg*RINGSIZE,RINGSIZE,External+i*RINGSIZE);
	    if (tperr<0) {
	      fprintf(stderr, "Impossible to read External[%ld]: %s %d %ld\n",ib,Param->External[ib+i*nbolo],tperr,(long) rg);
	      exit ( -1);
	    }
	  }
        }
        else {
          External = (PIOFLOAT *) _PIOMALLOC(sizeof(PIOFLOAT)*RINGSIZE);
          for (i=0;i<RINGSIZE;i++) External[i]=0;
        }

        //
        //====================================================================================

        ph = (PIODOUBLE *) malloc(sizeof(PIODOUBLE)*RINGSIZE);
        tperr = noDMC_readObject_PIODOUBLE(Param->Ptg[ib],rg*RINGSIZE,RINGSIZE,ph);
        if (tperr<0) {
          fprintf( stderr, "Impossible to read Ptg[%ld]: %s %d\n", ib, Param->Ptg[ib], tperr);
          exit ( -1);
        }
	char tpname[MAX_OUT_NAME_LENGTH];
        sprintf(tpname,"%s_TUPLE_1",Param->Ptg[ib]);
        th = (PIODOUBLE *) malloc(sizeof(PIODOUBLE)*RINGSIZE);
        tperr = noDMC_readObject_PIODOUBLE(tpname,rg*RINGSIZE,RINGSIZE,th);
        if (tperr<0) {
          fprintf(stderr, "Impossible to read Ptg[%ld]: %s %d\n", ib, tpname, tperr);
          exit ( -1);
        }

        sprintf(tpname,"%s_TUPLE_2",Param->Ptg[ib]);
        psi = (PIODOUBLE *) malloc(sizeof(PIODOUBLE)*RINGSIZE);
        tperr = noDMC_readObject_PIODOUBLE(tpname,rg*RINGSIZE,RINGSIZE,psi);
        if (tperr<0) {
          fprintf(stderr, "Impossible to read Ptg[%ld]: %s %d\n", ib, tpname, tperr);
          exit ( -1);
        }
        if (Param->flag_delta_psi == _PAR_TRUE) {
          //fprintf(stderr,"[DEBUG] Param->delta_psi[ib] = %f \n",Param->delta_psi[ib]);
          float delta_psi_radians = Param->delta_psi[ib] * M_PI / 180.0;
          for (i=0; i<RINGSIZE; i++) {
            psi[i] += delta_psi_radians;
            // psi must be in [-pi, +pi]
            if (psi[i] >  M_PI) psi[i] -= 2*M_PI;
            if (psi[i] < -M_PI) psi[i] += 2*M_PI;
          }
        }
	/*========================= ~ NORTH ECLIPTIC POLE COORDINATE (max hit count smoothed at 30 degree)
	  th,ph=(1.0899100645823487, 1.6705050780074633)
	  iarg=13523074
	  vec=(-0.08825391070996505, 0.8821818245983525, 0.46256510416666663)
	  lat,lon=(27.552753250600432, 95.712890625)
	  
	  External=0 is defined as the nearest point of the ring to this direction
	*/
	
#define EPOLEX (-0.08825391070996505)
#define EPOLEY (0.8821818245983525)
#define EPOLEZ (0.46256510416666663)
	
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
	    fprintf(stderr,"LUSOL POINTING %s %d %d\n",__FILE__,__LINE__,rank);
	  }
	  v0[0]=rres[0];v0[1]=rres[1];v0[2]=-1;
	  tmp=0;
	  for (k=0;k<3;k++) tmp+=v0[k]*v0[k];
	  for (k=0;k<3;k++) v0[k]/=sqrt(tmp);
	  double navr=0,avr=0,avr2=0;
	  //fprintf(stderr,"EXTERNAL v0 %ld %lf %lf %lf\n",(long) rg,v0[0],v0[1],v0[2]);
	  //fprintf(stderr,"EXTERNAL v1 %ld %lf %lf %lf\n",(long) rg,v1[0],v1[1],v1[2]);
	  //fprintf(stderr,"EXTERNAL v2 %ld %lf %lf %lf\n",(long) rg,v2[0],v2[1],v2[2]);
	  
	  PIOFLOAT vExternal0=1E30;
	  
	  for (i=0;i<RINGSIZE;i++) if (h[i]>0) {
	      double vec[3];
	      ang2vec(th[i],ph[i],vec);
	      double xx=vec[0]*v0[0]+vec[1]*v0[1]+vec[2]*v0[2];
	      double yy=vec[0]*v1[0]+vec[1]*v1[1]+vec[2]*v1[2];
	      
	      tmp=sqrt(xx*xx+yy*yy);
	      navr+=1;
	      avr+=tmp;
	      avr2+=tmp*tmp;
	      tmp=(vec[0]-(EPOLEX))*(vec[0]-(EPOLEX))+
		(vec[1]-(EPOLEY))*(vec[1]-(EPOLEY))+
		(vec[2]-(EPOLEZ))*(vec[2]-(EPOLEZ));
	      if (tmp<vExternal0) {
		vExternal0=tmp;
	      } 
	    }

	  long nrg_htmp = 0;

          for (i=0;i<RINGSIZE;i++) if (h[i]>0) {
            long ipix;
            long ipix128;
            long l_ipix128;
            nrg_htmp++;
            double vecpix[3];
            ang2pix_ring( Nside, th[i], ph[i], &ipix);
            ang2pix_ring( nside128, th[i], ph[i], &ipix128); // Beware taht nside128 could be different from 128
            ang2pix_ring( 128, th[i], ph[i], &l_ipix128);

	    ang2vec(th[i], ph[i],vecpix);
	    hpix *tp_hpix;
            
	    if (l_nhpix[ipix]==0) {
	      l_hpix[ipix] = (hpix *) malloc(sizeof(hpix));
	    }
	    else {
	      l_hpix[ipix] = (hpix *) realloc(l_hpix[ipix],(l_nhpix[ipix]+1)*sizeof(hpix));
	    }
	    tp_hpix=l_hpix[ipix]+l_nhpix[ipix];

	    memset(tp_hpix,0,sizeof(hpix));
	    
	    for (int iter = 0; iter < number_of_iterations; iter++) {
	      tp_hpix->listp[iter]=(y[iter][i]-Param->Monop[ib])/Param->Calibration[ib];
	    }
	    
	    tp_hpix->sig=tp_hpix->listp[0];
   
	    tp_hpix->hpr_cal   = hprtab_cal[i];
	    
	    tp_hpix->hit    = 1/sqrt(h[i]);
	    if (Sub_HPR != NULL) {
	      tp_hpix->Sub_HPR = Sub_HPR[i]*Param->SUB_HPRCOEF[ib];
	    } else {
	      tp_hpix->Sub_HPR = 0.0;
	    }
	    
	    tp_hpix->corr_cnn = 0;
	    
	    tp_hpix->hit = init_channels(tp_hpix,
					 projFunc,
					 psi[i],
					 External+i,
					 (rg-globalBeginRing)/((double) (globalEndRing-globalBeginRing)),
					 ipix,
					 ib,
					 tp_hpix->hit,
					 i,
					 &(tp_hpix->sig));

	    if (sparseFunc!=NULL) {
	      
	      tp_hpix->nShpr = calc_sparse_hpr(sparseFunc,
					       rg,
					       ib,
					       ipix,
					       psi[i],
					       External+i,
					       tp_hpix->listofShpr,
					       tp_hpix->listofShpr_idx);
	      
	      
	      for (j=0;j<tp_hpix->nShpr;j++) {
		if (npixShpr<tp_hpix->listofShpr_idx[j]+1) {
		  npixShpr=tp_hpix->listofShpr_idx[j]+1;
		}
	      }
	    }
	    if (diagFunc!=NULL) {
	      tp_hpix->diag_idx = calc_diag_hpr(diagFunc,
						rg,
						ib,
						ipix,
						psi[i],
						External+i);
	      
	    }
	    
	    
	    for (j=0;j<npixmap;j++) {
	      tp_hpix->listofmap[j]=0;
	      for (int k=0;k<MAXCHANNELS;k++) {
		      tp_hpix->listofmap[j]+=tp_hpix->channels[k]*mapmodel[j][ipix128];
	      }
	    }
	    
	    tp_hpix->model  = hprtab_cal[i];
	    for (int k=0;k<MAXCHANNELS;k++) {
	      tp_hpix->model  += tp_hpix->channels[k]*skymodel[k][ipix]; 
	    }

	    tp_hpix->w=h[i]*sxi;
	         
	    tp_hpix->surv = surv;
	    tp_hpix->ib = ib;
	    tp_hpix->rg = rg;
#define NBPH (4096)
	    tp_hpix->gi = 0;
	    tp_hpix->ipix = ipix;
	    l_nhpix[ipix]+=1;
	    if (isnan(tp_hpix->sig)||isnan(tp_hpix->w)||isnan(tp_hpix->hpr_cal)) {
	      fprintf(stderr,"NAN NAN NAN %lf %lf %lf %ld %ld\n",tp_hpix->hpr_cal,tp_hpix->sig,tp_hpix->w,(long)rg,(long)i);
	    }
	    
            }
	  



          if ((rank==0) && (rg==globalRankInfo.BeginRing[rank]) ) {
            GetProcMem(&vmem,&phymem);
	    fprintf(stderr,"Com %s Rank: %ld[%d] MEM %.1lf[%.1lf]MB Nd=%ld  %lg %lg \n",
		    Command, (long) rank, getpid(),
		    (double) vmem/1024./1024.,
		    (double) phymem/1024./1024.,
		    (long) nrg_htmp,
		    avr/navr,sqrt(avr2/navr-(avr/navr)*(avr/navr)));
          }
        }

        free(h);
        if (Param->flag_stim_paramfiles == 0) {
          free(y[0]);
        }
        free(Sub_HPR);
        free(hprtab_cal);
        free(External);
        free(th);
        free(ph);
        free(psi);
      } // if not badring
    } // end ring loop


    if (Param->flag_stim_paramfiles == 1) {
      for (int iter = 0; iter < number_of_iterations; iter++) {
        free( stim_hpr[iter]);
      }
    }

    if (rank==0) {
      GetProcMem(&vmem,&phymem);
      fprintf(stderr,"\nafter troll detector (%ld/%ld): used VMEM %.1lf[PHYS %.1lf]MB\n",
          ib, nbolo-1, (double) vmem/1024./1024., (double) phymem/1024./1024.);
    }

  } // end bolometer loop

  if (addpol!=NULL) {
    free(addpol);
  }

  if (skymodel!=NULL) {
    for (int l=0;l<MAXCHANNELS;l++) {
      free(skymodel[l]);
    }
  }

  if (mapmodel!=NULL) {
    for (int l=0;l<Param->n_External_MAP;l++) {
      free(mapmodel[l]);
    }
    free(mapmodel);
  }
  
  PrintFreeMemOnNodes( rank, mpi_size, "before pixel balancing");

  if (sparseFunc!=NULL) {
    long lnpix=npixShpr;
    long mnpix;
    MPI_Allreduce(&lnpix,&mnpix,1,MPI_LONG, MPI_MAX,MPI_COMM_WORLD);
    npixShpr=mnpix;
    ClosepFunc(sparseFunc);
    
    Py_DECREF(pClass);
    Py_DECREF(pArgs);
  }

  if (rank==0) fprintf(stderr,"NSHPR %d\n",(int) npixShpr);
  if (rank==0) fprintf(stderr,"NMAP %d\n",(int) npixmap);
 
  /*======================================================
    =
    =      compute load balancing between proc
    =
    =*/

  long *begpix = (long *) malloc(sizeof(long)*mpi_size);
  long *edpix = (long *) malloc(sizeof(long)*mpi_size);
  long rrk;

  PIOINT *mask;

  if (Param->flag_Mask==_PAR_TRUE) {
    PIOSTRING commask;
    sprintf(commask,"begin=0;end=%lld",(long long) 12*Nside*Nside);
    mask = (PIOINT *) malloc(sizeof(PIOINT)*(12*Nside*Nside));
    long resmask=(long) noDMC_readObject_PIOINT(Param->Mask,0,12*Nside*Nside,mask);
    if (rank==0) fprintf(stderr,"Mask %ld\n",(long) resmask);
  }
  else {
    mask = (PIOINT *) _PIOMALLOC(sizeof(PIOINT)*12*Nside*Nside);
    memset(mask,1,sizeof(PIOINT)*nnbpix);
  }

  if (rank==0) fprintf(stderr,"Sort pixel to make it run faster\n");
  nnbpix=12*Nside*Nside/mpi_size;
 
  for (i=0;i<mpi_size;i++) {
    begpix[i]=i*nnbpix;
    edpix[i]=(i+1)*nnbpix-1;
  }
  edpix[mpi_size-1]=12*Nside*Nside-1;
  
  int *stat_pix = (int *) malloc( 2*12*Nside*Nside*sizeof(int));
  memset(stat_pix,0,2*12*Nside*Nside*sizeof(int));
  realpix = (int *) malloc(sizeof(int)*nnbpix);
  irealpix = (int *) malloc(sizeof(int)*nnbpix);
  int *newmask = (int *) malloc(sizeof(int)*nnbpix);
  //int *inv_stat_pix;
  {
    long l_nr=12*Nside*Nside;
    for (i=0;i<l_nr;i++) {
      stat_pix[i]=l_nhpix[i];
    }

#ifdef OPTIMPI
    {
      int *l_stat_pix = (int *) malloc(12*Nside*Nside*sizeof(int));
      MPI_Reduce(stat_pix,l_stat_pix,12*Nside*Nside,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);

      memcpy(stat_pix,l_stat_pix,12*Nside*Nside*sizeof(int));
      free(l_stat_pix);
    }
#endif

    if (rank==0) {
#ifndef OPTIMPI
      MPI_Status statu;
      int *l_stat_pix = (int *) malloc(12*Nside*Nside*sizeof(int));
      for (rrk=1;rrk<mpi_size;rrk++) {
        MPI_Recv(l_stat_pix,sizeof(int)*12*Nside*Nside, MPI_BYTE, rrk,31, MPI_COMM_WORLD,&statu);
        for (j=0;j<12*Nside*Nside;j++) stat_pix[j]+=l_stat_pix[j];
      }
      free(l_stat_pix);
#endif
      fprintf(stderr,"regrid\n");
      hpint *statp = (hpint *) malloc(12*Nside*Nside*sizeof(hpint));

      for (j=0;j<12*Nside*Nside;j++) {
	      statp[j].ipx=j;
	      statp[j].hit=stat_pix[j]*mask[j];
      }

      fprintf(stderr,"sort\n");
      qsort(statp,12*Nside*Nside,sizeof(hpint),compar_int);

      for (rrk=0;rrk<mpi_size;rrk++) {
	for (i=0;i<nnbpix/4;i++) {

	  stat_pix[i           +(mpi_size-1-rrk)*nnbpix] = statp[i*mpi_size*4+rrk].ipx;
	  stat_pix[i+nnbpix/4  +(mpi_size-1-rrk)*nnbpix] = statp[i*mpi_size*4+2*mpi_size-1-rrk].ipx;
	  stat_pix[i+nnbpix/2  +rrk*nnbpix]              = statp[i*mpi_size*4+2*mpi_size+rrk].ipx;
	  stat_pix[i+3*nnbpix/4+rrk*nnbpix]              = statp[i*mpi_size*4+2*mpi_size+2*mpi_size-1-rrk].ipx;

#if 0
	  if (i==0) fprintf(stderr,"%d %d %d %d %d %d\n",(int) (mpi_size-1-rrk),
			    statp[i*mpi_size*2+rrk].hit,statp[i*mpi_size*2+2*mpi_size-1-rrk].hit,
			    (int) rrk,
			    statp[i*mpi_size*4+2*mpi_size+rrk].ipx,statp[i*mpi_size*4+2*mpi_size+2*mpi_size-1-rrk].ipx);
#endif
	  stat_pix[12*Nside*Nside+statp[i*mpi_size*4+rrk].ipx]=(mpi_size-1-rrk);
	  stat_pix[12*Nside*Nside+statp[i*mpi_size*4+2*mpi_size-1-rrk].ipx]=(mpi_size-1-rrk);
	  stat_pix[12*Nside*Nside+statp[i*mpi_size*4+2*mpi_size+rrk].ipx]=rrk;
	  stat_pix[12*Nside*Nside+statp[i*mpi_size*4+2*mpi_size+2*mpi_size-1-rrk].ipx]=rrk;
	}
      }
      free(statp);
#if 0
      FILE *fp=fopen("statpix.dat","w");
      fwrite(stat_pix,sizeof(int)*12*Nside*Nside,1,fp);
      fclose(fp);
#endif
    }
#ifndef OPTIMPI
    else  {
      MPI_Send(stat_pix, sizeof(int)*12*Nside*Nside, MPI_BYTE, 0, 31, MPI_COMM_WORLD);
    }
#endif
    // NOT A BCAST TO BE CHANGED TO MPI_SCATTER FOR OPTIMALITY
    MPI_Bcast(stat_pix,sizeof(int)*12*Nside*Nside*2, MPI_BYTE, 0, MPI_COMM_WORLD);

    for (i=0;i<nnbpix;i++) realpix[i]=stat_pix[i+rank*nnbpix];
    for (i=0;i<nnbpix;i++) irealpix[i]=stat_pix[12*Nside*Nside+i+rank*nnbpix];
    for (i=0;i<nnbpix;i++) newmask[i]=mask[realpix[i]];
  }
  free(mask);
  mask=newmask;
  // AND NOW IN FRONT OF YOUR EYES MIX PIXELS TO GET THEM EFFICIENTLY COMPUTED!!!

  if (rank==0) fprintf(stderr,"rebin pixel\n");
  hpix **new_hpix = (hpix **) malloc(sizeof(hpix *)*12*Nside*Nside);
  long l_nr=12*Nside*Nside;
  for (i=0;i<l_nr;i++) {
    new_hpix[i]=l_hpix[stat_pix[i]];
    if (l_nhpix[stat_pix[i]]>0) {
      for (j=0;j<l_nhpix[stat_pix[i]];j++) {
	      new_hpix[i][j].ipix=i;
      }
    }
    //new_hpix[i]=l_hpix[i];
  }
  free(l_hpix);
  l_hpix=new_hpix;

  PIOINT *new_nhpix = (PIOINT *) malloc(sizeof(PIOINT)*12*Nside*Nside);

  for (i=0;i<l_nr;i++) {
    //new_nhpix[i]=l_nhpix[i];
    new_nhpix[i]=l_nhpix[stat_pix[i]];
  }
  free(l_nhpix);
  l_nhpix=new_nhpix;
  free(stat_pix);

  GetProcMem(&vmem,&phymem);
  if (rank==0) fprintf(stderr,"Rank: %ld[%d] MEM %.1lf[%.1lf]MB line=%d\n",
              (long) rank, getpid(),
              (double) vmem/1024./1024.,
              (double) phymem/1024./1024.,__LINE__);
  MPI_Barrier(MPI_COMM_WORLD);

  /*======================================================
    =
    =      order data per pixel and per proc
    =
    =*/
  
#if 1
  int *begbuf = (int *) malloc(mpi_size*sizeof(int));
  hpix * tbs = NULL;
  long ntbs=0;
  
  for (int k=0;k<mpi_size;k++) {
    begbuf[k]=ntbs*sizeof(hpix);
    for (i=begpix[k];i<edpix[k]+1;i++) {
      if (l_nhpix[i]>0) {
	if (tbs==NULL) {
	  tbs = (hpix *) malloc(l_nhpix[i]*sizeof(hpix));
	}
	else {
	  tbs = (hpix *) realloc(tbs,(ntbs+l_nhpix[i])*sizeof(hpix));
	}
	memcpy(tbs+ntbs,l_hpix[i],l_nhpix[i]*sizeof(hpix));
	ntbs+=l_nhpix[i];
	free(l_hpix[i]);
      }
    }
  }
  free(l_hpix);
  free(l_nhpix);
  
  begbuf[mpi_size]=ntbs*sizeof(hpix);
       	    
  /* exchange all buffer information */
  int *allbuf = (int *) malloc(mpi_size*(mpi_size+1)*sizeof(int));
  
  MPI_Allgather(begbuf,mpi_size+1,MPI_INT,allbuf,(mpi_size+1),MPI_INT,MPI_COMM_WORLD);

  int *sdispls = (int *) malloc(sizeof(int)*mpi_size);
  int *sendcounts = (int *) malloc(sizeof(int)*mpi_size);
  int *recvcounts = (int *) malloc(sizeof(int)*mpi_size);
  int *rdispls = (int *) calloc(mpi_size,sizeof(int));
  
  for (i=0;i<size;i++) {
    sendcounts[i]=begbuf[i+1]-begbuf[i];
    sdispls[i]=begbuf[i];
    recvcounts[i]=allbuf[i*(size+1)+rank+1]-allbuf[i*(size+1)+rank];
  }
  
  free(begbuf);
  free(allbuf);
  
  rdispls[0]=0;
  for (i=1;i<size;i++) {
    rdispls[i]=rdispls[i-1]+recvcounts[i-1];
  }

  long ldata=(rdispls[size-1]+recvcounts[size-1])/sizeof(hpix);

  hpix *otbs = (hpix *) malloc(ldata*sizeof(hpix));

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Alltoallv(tbs,sendcounts,sdispls,MPI_BYTE,
		otbs,recvcounts,rdispls,MPI_BYTE,
		MPI_COMM_WORLD);

  free(sdispls);
  free(sendcounts);
  free(recvcounts);
  free(rdispls);
  free(tbs);
  tbs=otbs;
  
  nnbpix = edpix[rank]-begpix[rank]+1;

  loc_hpix = (hpix **) malloc(sizeof(hpix *)*nnbpix);
  loc_nhpix = (PIOINT *) malloc(sizeof(PIOINT)*nnbpix);
  memset(loc_nhpix,0,sizeof(PIOINT)*nnbpix);

  for (int kk=0;kk<ldata;kk++) {
    j=tbs[kk].ipix-begpix[rank];
    if (loc_nhpix[j]==0) {
      loc_hpix[j] = (hpix *) malloc(sizeof(hpix));
    }
    else {
      loc_hpix[j] = (hpix *) realloc(loc_hpix[j],(1+loc_nhpix[j])*sizeof(hpix));
    }
    memcpy(loc_hpix[j]+loc_nhpix[j],tbs+kk,sizeof(hpix));
    loc_nhpix[j] +=1;
  }
  
  free(tbs);
#else
  nnbpix = edpix[rank]-begpix[rank]+1;

  loc_hpix = (hpix **) malloc(sizeof(hpix *)*nnbpix);
  loc_nhpix = (PIOINT *) malloc(sizeof(PIOINT)*nnbpix);
  memset(loc_nhpix,0,sizeof(PIOINT)*nnbpix);

  MPI_Barrier(MPI_COMM_WORLD);
  
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
          //fprintf (stderr,"%d <- %d\n",(int)rank,(int)(rk0+(k+exch/2)%(exch)));
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
          for (j=begpix[js+k];j<=edpix[js+k];j++){
	    if (l_nhpix[j]>0) {
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


  memcpy(loc_nhpix,l_nhpix+begpix[rank],sizeof(PIOINT)*(nnbpix));
  for (j=0;j<nnbpix;j++) loc_hpix[j]=l_hpix[j+begpix[rank]];

  free(l_hpix);
  free(l_nhpix);
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
  fprintf(stderr,"Rank: %ld[%d] MEM %.1lf[%.1lf]MB line=%d\n",
              (long) rank, getpid(),
              (double) vmem/1024./1024.,
              (double) phymem/1024./1024.,__LINE__);
  PIOLONG k;
  long i0;
  flgpix = (PIOBYTE *) malloc(sizeof(PIOBYTE)*nnbpix);
  memset(flgpix,0,sizeof(PIOBYTE)*nnbpix);

  imatrice = (PIODOUBLE *) malloc(nnbpix*MAXCHANNELS*MAXCHANNELS*sizeof(PIODOUBLE));
  memset(imatrice,0,nnbpix*MAXCHANNELS*MAXCHANNELS*sizeof(PIODOUBLE));
  
  cond = (PIODOUBLE *) malloc(nnbpix*sizeof(PIODOUBLE));
  memset(imatrice,0,nnbpix*sizeof(PIODOUBLE));
  
  long l_nmatpix=0;
  fprintf(stderr,"Rank: %ld[%d] MEM %.1lf[%.1lf]MB line=%d\n",
              (long) rank, getpid(),
              (double) vmem/1024./1024.,
	(double) phymem/1024./1024.,__LINE__);
  for (k=0;k<nnbpix;k++) {
      
    cond[k] = 0.0;
    double matrice[MAXCHAN*MAXCHAN];
    memset(matrice,0,MAXCHANNELS*MAXCHANNELS*sizeof(double)); 
    
    long ndata = loc_nhpix[k];
    if (ndata>2) {
      hpix *htmp = loc_hpix[k];
      for (i0=0;i0<ndata;i0++) { 
	      for(int i = 0;i<MAXCHANNELS;i++){        
	        for(int j = 0;j<MAXCHANNELS;j++){
	          matrice[i+j*MAXCHANNELS] += htmp[i0].w *htmp[i0].channels[j]*htmp[i0].channels[i];
	        }
	      } 
      }     
      
      //test if matrice inversible 
      cond[k]=cond_thres(matrice,imatrice+MAXCHANNELS*MAXCHANNELS*k,MAXCHANNELS);
      
      if (cond[k] < Param->seuilcond) {
	      flgpix[k]=1;
	      for (i0=0;i0<ndata;i0++) {
	        flg_rg[htmp[i0].ib][htmp[i0].rg-globalBeginRing]=1;
	      }
      }else {
	      flgpix[k]=0;
	      // Print detail for pix that does not meet cond requirement
	      //fprintf(stderr,"[DBG COND] Pix#%ld is flagged out! II=%g IQ=%g IU=%g QQ=%g QU=%g UU=%g \n",k+begpix[rank], II, IQ, IU, QQ, QU, UU);
	      cond[k]=1E30;
      }
    }else {
      //fprintf(stderr,"[DBG COND] Pix#%ld is flagged out! ndata<=2\n", k+begpix[rank]);
      flgpix[k]=0;
      cond[k]=1E30;
    }
  }
  fprintf(stderr,"Rank: %ld[%d] MEM %.1lf[%.1lf]MB line=%d\n",
              (long) rank, getpid(),
              (double) vmem/1024./1024.,
              (double) phymem/1024./1024.,__LINE__);
  
  l_nmatpix=0;
  for (j=0;j<nnbpix;j++) {
    if (flgpix[j]==1) {
      l_nmatpix++;
    }
  }


  if (rank==0) fprintf(stderr,"l_nmatpix[%d] %ld / %ld \n",rank,(long) l_nmatpix ,nnbpix );

#ifdef OPTIMPI
  long lb;
  MPI_Reduce(&l_nmatpix,&lb,1,MPI_LONG,MPI_SUM,0,MPI_COMM_WORLD);
  nmatpix=lb;
#else
  if (rank==0) {
    long lb;
    nmatpix=l_nmatpix;
    for (rrk=1;rrk<mpi_size;rrk++) {
      MPI_Recv(&lb,sizeof(long), MPI_BYTE, rrk,450, MPI_COMM_WORLD,&statu);
      nmatpix+=lb;
    }
  }
  else MPI_Send(&l_nmatpix, sizeof(long), MPI_BYTE, 0, 450, MPI_COMM_WORLD);
#endif
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

  GetProcMem(&vmem,&phymem);
  if (rank%64==0) fprintf(stderr,"Rank: %ld[%d] MEM %.1lf[%.1lf]MB line=%d\n",
              (long) rank, getpid(),
              (double) vmem/1024./1024.,
              (double) phymem/1024./1024.,__LINE__);

  rgord = (PIOINT **) malloc(nbolo*sizeof(PIOINT *));
  rgordinv = (PIOINT **) malloc(nbolo*sizeof(PIOINT *));
  
  newnr  = (PIOLONG *) malloc((nbolo+1)*sizeof(PIOLONG));
  
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
  
  if (do_offset==0) {
    for (ib=1;ib<nbolo+1;ib++) newnr[ib]=0;
  }
  
  if (rank==0) {
    fprintf(stderr,"RK%d SHOULD DETERMINE %ld VALUES ",rank,(long) nnbpix);
    for (i=0;i<nbolo+1;i++) fprintf(stderr,"%ld ",(long) newnr[i]);
    fprintf(stderr,"\n");
  }
  
 //========================================================================


  double *histo_gi = (double *) malloc(32000*sizeof(double)*nbolo);
  double *histo2_gi = (double *) malloc(32000*sizeof(double)*nbolo);
  double *histon_gi = (double *) malloc(32000*sizeof(double)*nbolo);

  double mat_dip[4*nbolo];
  double vec_dip[2*nbolo];

  memset(mat_dip,0,4*sizeof(double)*nbolo);
  memset(vec_dip,0,2*sizeof(double)*nbolo);

  memset(histo_gi ,0,nbolo*32000*sizeof(double));
  memset(histo2_gi ,0,nbolo*32000*sizeof(double));
  memset(histon_gi ,0,nbolo*32000*sizeof(double));

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
	  
	  
	  DI+=htmp[l1].w*htmp[l1].model;
	  
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
	  
	  double tmp=(htmp[l1].model);
	  
	  histo_gi[ htmp[l1].rg+htmp[l1].ib*32000]+=ww*tmp;
	  histo2_gi[htmp[l1].rg+htmp[l1].ib*32000]+=ww*tmp*tmp;
	  histon_gi[htmp[l1].rg+htmp[l1].ib*32000]+=ww;
	  
	  mat_dip[0+4*htmp[l1].ib]+=ww*htmp[l1].model*htmp[l1].model;
	  mat_dip[1+4*htmp[l1].ib]+=ww*htmp[l1].model;
	  mat_dip[2+4*htmp[l1].ib]+=ww*htmp[l1].model;
	  mat_dip[3+4*htmp[l1].ib]+=ww;
	  
	  vec_dip[0+2*htmp[l1].ib]+=ww*(htmp[l1].sig)*htmp[l1].model;
	  vec_dip[1+2*htmp[l1].ib]+=ww*(htmp[l1].sig);
	}
      }
    }
  }


  //============================================================
  // remove dipole from the first harmonic
  //

  double l_mat[4*100];
  double l_vec[2*100];
  MPI_Reduce(mat_dip,l_mat,4*nbolo,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(vec_dip,l_vec,2*nbolo,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

  if (rank==0) {

    for (j=0;j<nbolo;j++) {

      lusol(mat_dip+4*j,vec_dip+2*j,2);

      fprintf(stderr,"GAIN DIP BOL (%ld/%ld): %lf\n", j, nbolo-1, vec_dip[2*j]);
    }
  }

  free(histo_gi);
  free(histo2_gi);
  free(histon_gi);

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

  // decoupe vecteur matrice par proc

  long maxsizemat=newnr[nbolo]+nbolo*(GAINSTEP+npixmap+npixShpr);
  
  if (rank==0) fprintf(stderr,"MAXSIZE %ld %ld\n",(long) maxsizemat,(long) newnr[nbolo]+nbolo*(GAINSTEP+npixmap+npixShpr));

  x2 =     (double *) malloc (maxsizemat*sizeof (double));
  x2old =  (double *) malloc (maxsizemat*sizeof (double));
  x2init = (double *) malloc (maxsizemat*sizeof (double));
  b2 =     (double *) malloc (maxsizemat*sizeof (double));
  d2 =     (double *) malloc (maxsizemat*sizeof (double));
  q2 =     (double *) malloc (maxsizemat*sizeof (double));
  r2 =     (double *) malloc (maxsizemat*sizeof (double));
  s2 =     (double *) malloc (maxsizemat*sizeof (double));
  hit2 =   (double *) malloc (maxsizemat*sizeof (double));
  hit_WW =   (double *) malloc (maxsizemat*sizeof (double));
  hit_L2 =   (double *) malloc (maxsizemat*sizeof (double));
  hit_L1 =   (double *) malloc (maxsizemat*sizeof (double));


  for (i=0;i<Param->n_Out_MAP;i++) strcpy(mapout[i],Param->Out_MAP[i]);

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

    ittt=1;
    memset(x2     ,0,maxsizemat*sizeof (double));
    memset(x2old  ,0,maxsizemat*sizeof (double));
    memset(x2init ,0,maxsizemat*sizeof (double));
    memset(b2     ,0,maxsizemat*sizeof (double));
    memset(d2     ,0,maxsizemat*sizeof (double));
    memset(q2     ,0,maxsizemat*sizeof (double));
    memset(r2     ,0,maxsizemat*sizeof (double));
    memset(s2     ,0,maxsizemat*sizeof (double));
    memset(hit2   ,0,maxsizemat*sizeof (double));
    memset(hit_WW   ,0,maxsizemat*sizeof (double));
    memset(hit_L1   ,0,maxsizemat*sizeof (double));
    memset(hit_L2   ,0,maxsizemat*sizeof (double));


    MPI_Barrier(MPI_COMM_WORLD);
    // Exchange dip against sig in the XI2

    double *gain;
    if (GAINSTEP==0) {
      gain = (double *) malloc(sizeof(double)*nbolo);

      for (i=0;i<nbolo;i++) {
	gain[i]=1; //+3E-3*i+1E-3*cos(j/4.);
      }
    }
    else {
      gain = (double *) malloc(sizeof(double)*nbolo*GAINSTEP);

      for (i=0;i<nbolo;i++) {
	for (j=0;j<GAINSTEP;j++) gain[i*GAINSTEP+j]=1; //+3E-3*i+1E-3*cos(j/4.);
      }
    }

    memset(x2,0,(newnr[nbolo]+nbolo*(GAINSTEP+npixmap+npixShpr))*sizeof(double));

    gainoff=0;

    nmatres=newnr[nbolo]+nbolo*(GAINSTEP+npixmap+npixShpr);

    double *x3= (double *) malloc(sizeof(double)*(nmatres)); 

    GetProcMem(&vmem,&phymem);
    if (rank==0) fprintf(stderr,"Rank: %ld[%d] MEM %.1lf[%.1lf]MB line=%d\n",
                (long) rank, getpid(),
                (double) vmem/1024./1024.,
                (double) phymem/1024./1024.,__LINE__);
    int itt;

    memset(x3,0,sizeof(double)*nmatres);

    itt=0;
    double resxi=1;

    while (itt<Param->NITT) {
      memset(newnr[nbolo]+x3,0,sizeof(double)*nbolo*GAINSTEP);
     
      minimize_gain_tf(x3,gain);
     
      resxi=0;
      if (GAINSTEP>0) {
	for (i=0;i<nbolo;i++){
	  for (j=0;j<GAINSTEP;j++) {
	    resxi+=x3[newnr[nbolo]+i*GAINSTEP+j]*x3[newnr[nbolo]+i*GAINSTEP+j];
	  }
	}
	resxi/=(nbolo*GAINSTEP);
      }


      if (GAINSTEP>0) {
	for (i=0;i<nbolo;i++){
	  for (j=0;j<GAINSTEP;j++) {
	    if (itt<Param->NITT-1) gain[i*GAINSTEP+j]-=x3[newnr[nbolo]+i*GAINSTEP+j];
	  }
	}
      }

      if (rank==0) {
        fprintf(stderr,"GI XIGAIN %.10lg\n",sqrt(resxi));
	if (do_offset==1) {
	  fprintf(stderr,"MEAN_OFF=[");
	  for (i=0;i<nbolo;i++) {
	    double  osum=0;
	    for (int ll=newnr[i];ll<newnr[i+1];ll++) osum+=x3[ll];
	    fprintf(stderr,"%lg,",osum/(newnr[i+1]-newnr[i]));
	  }
	  fprintf(stderr,"]\n");
	}
	if (do_offset==1) {
	  fprintf(stderr,"MEAN2_OFF=[");
	  for (i=0;i<nbolo;i++) {
	    double  osum=0;
	    for (int ll=newnr[i];ll<newnr[i+1];ll++) osum+=x3[ll]*x3[ll];
	    fprintf(stderr,"%lg,",sqrt(osum)/(newnr[i+1]-newnr[i]));
	  }
	}
	fprintf(stderr,"]\n");
	if (GAINSTEP>0) {
	  for (i=newnr[nbolo];i<nbolo+newnr[nbolo];i++) {
	    if ((i-newnr[nbolo])%nbolo==0) fprintf(stderr,"GAIN=[");
	    fprintf(stderr,"%lg,",x3[i]);
	    if ((i-newnr[nbolo])%nbolo==nbolo-1) fprintf(stderr,"]\n");
	  }
	}
	for (i=GAINSTEP*nbolo+newnr[nbolo];i<(GAINSTEP)*nbolo+newnr[nbolo];i++) {
	  if ((i-newnr[nbolo])%nbolo==0) fprintf(stderr,"TF=[");
	  fprintf(stderr,"%lg,",x3[i]);
	  if ((i-newnr[nbolo])%nbolo==nbolo-1) fprintf(stderr,"]\n");
	}	
	for (i=(GAINSTEP)*nbolo+newnr[nbolo];i<(GAINSTEP+npixmap)*nbolo+newnr[nbolo];i++) {
	  if ((i-newnr[nbolo])%nbolo==0) fprintf(stderr,"TMAP=[");
	  fprintf(stderr,"%lg,",x3[i]);
	  if ((i-newnr[nbolo])%nbolo==nbolo-1) fprintf(stderr,"]\n");
	}
	int nbSTF=npixShpr;
	if (nbSTF>200) nbSTF=200;
	for (i=(GAINSTEP+npixmap)*nbolo+newnr[nbolo];i<(GAINSTEP+npixmap+nbSTF)*nbolo+newnr[nbolo];i++) {
	  if ((i-newnr[nbolo])%nbolo==0) fprintf(stderr,"STF=[");
	  fprintf(stderr,"%lg,",x3[i]);
	  if ((i-newnr[nbolo])%nbolo==nbolo-1) fprintf(stderr,"]\n");
	}	
      }
      itt++;
      MPI_Barrier(MPI_COMM_WORLD);
    }


    MPI_Barrier(MPI_COMM_WORLD);

    nmatres=newnr[nbolo]+nbolo*(GAINSTEP+npixmap+npixShpr);
    if (rank==0) {

      PIOSTRING commm;
      PIOSTRING saveg;
      if (do_offset==1) {
	// Write offset
	PIODOUBLE *tmpoff = (PIODOUBLE *) malloc(sizeof(PIODOUBLE)*(globalRangeRing));
	
	for (i=0;i<nbolo;i++) {
	  sprintf(commm,"begin=%lld;end=%lld",(long long) globalBeginRing,(long long) globalEndRing);
	  for (j=0;j<globalRangeRing;j++) {
	    tmpoff[j]=-10000;
	  }
	  for (j=newnr[i];j<newnr[i+1];j++) tmpoff[rgordinv[i][j-newnr[i]]]=x3[j];
	  if (stim_first_seed+iter==0) sprintf(saveg,"%s_OFF",Param->Out_VEC[i]);
	  else sprintf(saveg,"%s_%d_OFF",Param->Out_VEC[i],stim_first_seed+iter);
	  
	  fprintf(stderr,"Write OFF  %lld\n",(long long) PIOWriteVECT(saveg,tmpoff,sizeof(PIODOUBLE)*globalBeginRing,sizeof(PIODOUBLE)*(globalRangeRing)));
	}
	free(tmpoff);
      }

      if (stim_first_seed+iter==0) sprintf(saveg,"%s_X2",Param->Out_VEC[0]);
      else sprintf(saveg,"%s_%d_X2",Param->Out_VEC[0],stim_first_seed+iter);

      //for (i=newnr[nbolo]+GAINSTEP*nbolo;i<nmatres;i++) fprintf(stderr,"X2 %d %lg\n",(int) i,x3[i]);
      fprintf(stderr,"Write MAT  %lld\n",(long long) PIOWriteVECT(saveg,x3+newnr[nbolo],0,(nbolo*(GAINSTEP+npixmap+npixShpr))*sizeof(PIODOUBLE)));

    }

    double avvgain=0;
    if (GAINSTEP>0) {
      for (i=0;i<nbolo*GAINSTEP;i++) avvgain+=x3[newnr[nbolo]+i];
      avvgain/=((double)(nbolo*GAINSTEP));

      if (rank==0)  {
	fprintf(stderr,"AVVGAIN: %lf\n",avvgain);
      }
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
      free(hit_WW);
      free(hit_L1);
      free(hit_L2);
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
  
  PIOLONG GAINSTEP2 = GAINSTEP;

  double ** map = (double **) malloc(sizeof(double*)*MAXCHANNELS);
  for (i = 0;i<MAXCHANNELS;i++){
    map[i]=(double *) malloc(sizeof(double)*nnbpix);
  }
  double ** rmap = (double **) malloc(sizeof(double*)*MAXCHANNELS);
  for (i = 0;i<MAXCHANNELS;i++){
    rmap[i]=(double *) malloc(sizeof(double)*nnbpix);
  }
  double * cmap = (double *) malloc(sizeof(double)*nnbpix);
  
  double *diag_avv=NULL;
  double *diag_avv2=NULL;
  double *diag_n=NULL;
  
  if (diagFunc!=NULL) {
    diag_avv  = (double *) malloc(sizeof(double)*nb_diag);
    diag_avv2  = (double *) malloc(sizeof(double)*nb_diag);
    diag_n = (double *) malloc(sizeof(double)*nb_diag);
  }
  
  double ** imap = (double **) malloc(sizeof(double*)*MAXCHANNELS);
  for (i = 0;i<MAXCHANNELS;i++){
    imap[i]=(double *) malloc(sizeof(double)*nnbpix);
  }
  double ** matmap = (double **) malloc(sizeof(double*)*MAXCHANNELS*MAXCHANNELS);
  for (i = 0;i<MAXCHANNELS*MAXCHANNELS;i++){
    matmap[i]=(double *) malloc(sizeof(double)*nnbpix);
  }
  
  // Init matrice and vecteur
  double *matrix = malloc(MAXCHANNELS*MAXCHANNELS*sizeof(double)); 
  double *Imatrix = malloc(MAXCHANNELS*MAXCHANNELS*sizeof(double));
  double *vector = malloc(MAXCHANNELS*sizeof(double));
  double *rvector = malloc(MAXCHANNELS*sizeof(double));
  
  if (rank%64==0) {
    GetProcMem(&vmem,&phymem);
    fprintf(stderr,"Rank: %ld[%d] MEM %.1lf[%.1lf]MB line=%d\n",
	    (long) rank, getpid(),
	    (double) vmem/1024./1024.,
	    (double) phymem/1024./1024.,__LINE__);
  }
  
  for (int detset=0; detset<Param->n_Out_MAP; detset++) {
    for (int isurv=0; isurv<Param->MAPRINGS[detset]; isurv++) {
      PIOSTRING mapname;
      strcpy(mapname, Param->name_surv[isurv]);
      
      if (diagFunc!=NULL) {
	memset(diag_avv,0,sizeof(double)*nb_diag);
	memset(diag_avv2,0,sizeof(double)*nb_diag);
	memset(diag_n,0,sizeof(double)*nb_diag);
      }
      int l1;
      for (k=0;k<nnbpix;k++) {
	
      for(int i =0;i<MAXCHANNELS;i++){
	map[i][k]=UNSEENPIX;  // ie hp.UNSEEN
	imap[i][k]=UNSEENPIX;  // ie hp.UNSEEN
	rmap[i][k]=UNSEENPIX;  // ie hp.UNSEEN
      }
      cmap[k]=UNSEENPIX;  // ie hp.UNSEEN
	
      for(int i =0;i<MAXCHANNELS*MAXCHANNELS;i++){
	matmap[i][k]=UNSEENPIX;  // ie hp.UNSEEN
      }
      
	    
      long ndata = loc_nhpix[k];
      hpix *htmp = loc_hpix[k];
      //fprintf(stderr,"%s %d\n",__FILE__,__LINE__);
      
      //buildmap
      
      memset(matrix,0,MAXCHANNELS*MAXCHANNELS*sizeof(double));
      memset(Imatrix,0,MAXCHANNELS*MAXCHANNELS*sizeof(double));
      memset(vector,0,MAXCHANNELS*sizeof(double));
      memset(rvector,0,MAXCHANNELS*sizeof(double));
      
      for (l1=0;l1<ndata;l1++) {
	// select only bolometers in detset
	int use_bolo = 0;
	
	if (Param->bolomask[detset*nbolo+htmp[l1].ib]==1 && htmp[l1].rg>Param->beg_surv[isurv]&&htmp[l1].rg<=Param->end_surv[isurv]) {
	  use_bolo=1;
	}
	
	long ri1=htmp[l1].rg-globalBeginRing;
	
	if (flg_rg[htmp[l1].ib][ri1]!=0 && use_bolo==1) {
	  long iri1=rgord[htmp[l1].ib][ri1]+newnr[htmp[l1].ib];
	  //calcul signal corriger
	  double g1=gain[htmp[l1].gi+htmp[l1].ib*GAINSTEP];
	  double rsig = htmp[l1].sig*g1;
	  
	  double sig_corr = htmp[l1].sig*g1 - htmp[l1].Sub_HPR-htmp[l1].corr_cnn;
	  double sig_corr2 = sig_corr;

	  if (do_offset==1) {
	    sig_corr-=x3[iri1];
	  }
	  
	  if (REMOVE_CAL==1) {
	    sig_corr-=htmp[l1].hpr_cal;
	    sig_corr2-=htmp[l1].hpr_cal;
	  }
	  
	  
	  for(int m =0;m<npixmap;m++){
	    sig_corr-=x3[newnr[nbolo]+htmp[l1].ib+(m+GAINSTEP2)*nbolo]*htmp[l1].listofmap[m];
	  }
	  
	  for(int m =0;m<htmp[l1].nShpr;m++){
	    sig_corr-=x3[newnr[nbolo]+htmp[l1].ib+(htmp[l1].listofShpr_idx[m]+GAINSTEP2+npixmap)*nbolo]*htmp[l1].listofShpr[m];
	  }
	  
	  sig_corr-=x3[newnr[nbolo]+htmp[l1].ib*GAINSTEP2]*htmp[l1].model;  
	  //calcul matrix & vector
	  for(int i = 0;i<MAXCHANNELS;i++){  
	    vector[i]+= htmp[l1].w *htmp[l1].channels[i]*sig_corr;
	  }
	  for(int i = 0;i<MAXCHANNELS;i++){  
	    rvector[i]+= htmp[l1].w *htmp[l1].channels[i]*rsig;
	  }
	  for(int i = 0;i<MAXCHANNELS;i++){        
	    for(int j = 0;j<MAXCHANNELS;j++){
	      matrix[i+j*MAXCHANNELS] += htmp[l1].w *htmp[l1].channels[j]*htmp[l1].channels[i];
	    }
	  } 
	}
      }
	             
      cond[k]=cond_thres(matrix,Imatrix,MAXCHANNELS);  
      
      if (cond[k] < Param->seuilcond) {	               
	invertMatrix(Imatrix,vector,MAXCHANNELS,rank);	       
	for(int i = 0;i<MAXCHANNELS;i++){
	  map[i][k]= vector[i];
	}      
	invertMatrix(Imatrix,rvector,MAXCHANNELS,rank);	       
	for(int i = 0;i<MAXCHANNELS;i++){
	  rmap[i][k]= rvector[i];
	}
	for(int i = 0;i<MAXCHANNELS*MAXCHANNELS;i++){
	  matmap[i][k]=Imatrix[i];
	}

	/*================================================================================================
	  COMPUTE CHI2 MAP 
	  ================================================================================================*/
	
	double avv=0,navv=0,avv2=0;
	
	for (l1=0;l1<ndata;l1++) {
	  // select only bolometers in detset
	  int use_bolo = 0;
	  
	  if (Param->bolomask[detset*nbolo+htmp[l1].ib]==1 && htmp[l1].rg>Param->beg_surv[isurv]&&htmp[l1].rg<=Param->end_surv[isurv]) {
	    use_bolo=1;
	  }
	  
	  long ri1=htmp[l1].rg-globalBeginRing;
	  
	  if (flg_rg[htmp[l1].ib][ri1]!=0 && use_bolo==1) {
	    long iri1=rgord[htmp[l1].ib][ri1]+newnr[htmp[l1].ib];
	    //calcul signal corriger
	    double g1=gain[htmp[l1].gi+htmp[l1].ib*GAINSTEP];
	    double sig_corr = htmp[l1].sig*g1 - htmp[l1].Sub_HPR-htmp[l1].corr_cnn;
	    double sig_corr2 = sig_corr;

	    if (do_offset==1) {
	      sig_corr-=x3[iri1];
	    }
	    
	    if (REMOVE_CAL==1) {
	      sig_corr-=htmp[l1].hpr_cal;
	      sig_corr2-=htmp[l1].hpr_cal;
	    }
	    
	    
	    for(int m =0;m<npixmap;m++){
	      sig_corr-=x3[newnr[nbolo]+htmp[l1].ib+(m+GAINSTEP2)*nbolo]*htmp[l1].listofmap[m];
	    }
	    
	    for(int m =0;m<htmp[l1].nShpr;m++){
	      sig_corr-=x3[newnr[nbolo]+htmp[l1].ib+(htmp[l1].listofShpr_idx[m]+GAINSTEP2+npixmap)*nbolo]*htmp[l1].listofShpr[m];
	    }
	    
	    sig_corr-=x3[newnr[nbolo]+htmp[l1].ib*GAINSTEP2]*htmp[l1].model;  
	    //calcul matrix & vector
	    for(int i = 0;i<MAXCHANNELS;i++){  
	      sig_corr-= htmp[l1].channels[i]*map[i][k];
	    }
	    avv=avv+htmp[l1].w *sig_corr;
	    avv2=avv2+htmp[l1].w *sig_corr*sig_corr;
	    navv=navv+htmp[l1].w;
	    if (diagFunc!=NULL) {
	      diag_avv[htmp[l1].diag_idx]+=htmp[l1].w *sig_corr;
	      diag_avv2[htmp[l1].diag_idx]+=htmp[l1].w *sig_corr*sig_corr;
	      diag_n[htmp[l1].diag_idx]+=htmp[l1].w;
	    }
	  }
	}

	cmap[k]=sqrt(avv2/navv-(avv/navv)*(avv/navv));
      }
    }

    if (diagFunc!=NULL) {
      double *l_avv = (double *) malloc(sizeof(double)*(nb_diag));
      double *l_avv2 = (double *) malloc(sizeof(double)*(nb_diag));
      double *l_n = (double *) malloc(sizeof(double)*(nb_diag));
      MPI_Reduce(diag_avv,l_avv,nb_diag,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Reduce(diag_avv2,l_avv2,nb_diag,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Reduce(diag_n,l_n,nb_diag,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

      for (k=0;k<nb_diag;k++) diag_avv[k]=sqrt(l_avv2[k]/l_n[k]-(l_avv[k]/l_n[k])*(l_avv[k]/l_n[k]));

      if  (rank==0) {
	char TEST_OUTMAP[MAX_OUT_NAME_LENGTH];
	sprintf(TEST_OUTMAP,"%s_%s_DIAG", mapout[detset],mapname);
	PIOWriteVECT(TEST_OUTMAP,diag_avv,0,sizeof(PIODOUBLE)*nb_diag);
      }

      free(l_avv);
      free(l_avv2);
      free(l_n);
    }
    
    
    for(int i = 0;i<MAXCHANNELS;i++){
      char TEST_OUTMAP[MAX_OUT_NAME_LENGTH];
      sprintf(TEST_OUTMAP,"%s_%s_%d", mapout[detset],mapname,i);
      fprintf(stderr,"%s %d %d\n",__FILE__,__LINE__,rank);
      PIOWriteMAP(TEST_OUTMAP,map[i],begpix[rank],begpix[rank]+nnbpix-1);
      MPI_Barrier(MPI_COMM_WORLD);
    }
   
    for(int i = 0;i<MAXCHANNELS;i++){
      char TEST_OUTMAP[MAX_OUT_NAME_LENGTH];
      fprintf(stderr,"%s %d %d\n",__FILE__,__LINE__,rank);
      sprintf(TEST_OUTMAP,"%s_%s_%d_RAW", mapout[detset],mapname,i);
      PIOWriteMAP(TEST_OUTMAP,rmap[i],begpix[rank],begpix[rank]+nnbpix-1);
      MPI_Barrier(MPI_COMM_WORLD);
    }
   
    {
      char TEST_OUTMAP[MAX_OUT_NAME_LENGTH];
      fprintf(stderr,"%s %d %d\n",__FILE__,__LINE__,rank);
      sprintf(TEST_OUTMAP,"%s_%s_STD", mapout[detset],mapname);
      PIOWriteMAP(TEST_OUTMAP,cmap,begpix[rank],begpix[rank]+nnbpix-1);
      MPI_Barrier(MPI_COMM_WORLD);
    }
#if 0
    for(int i = 0;i<MAXCHANNELS*MAXCHANNELS;i++){
      char TEST_OUTMAP[MAX_OUT_NAME_LENGTH];
      sprintf(TEST_OUTMAP,"%s_%s_MAT%d", mapout[detset],mapname,i);
      PIOWriteMAP(TEST_OUTMAP,matmap[i],begpix[rank],begpix[rank]+nnbpix-1);
      MPI_Barrier(MPI_COMM_WORLD);
    }
#endif
    {
      char TEST_OUTMAP[MAX_OUT_NAME_LENGTH];
      fprintf(stderr,"%s %d %d\n",__FILE__,__LINE__,rank);
      sprintf(TEST_OUTMAP,"%s_%s_COND", mapout[detset],mapname);
      PIOWriteMAP(TEST_OUTMAP,cond,begpix[rank],begpix[rank]+nnbpix-1);
      MPI_Barrier(MPI_COMM_WORLD);
    }
    }
  }

  for(int i = 0;i<MAXCHANNELS;i++){
    free(map[i]);
    free(rmap[i]);
  }
  free(cmap);
  free(map);
  free(rmap);

  if (diagFunc!=NULL) {
    free(diag_avv);
    free(diag_avv2);
    free(diag_n);
  }
  
  if (rank==0) {
    now = time( NULL);
    fprintf(stderr, "\n%s: --------------------------\n", __FILE__ );
    fprintf(stderr, "%s: Finished successfully at %s",   __FILE__, ctime( &now));
    fprintf(stderr, "%s: --------------------------\n", __FILE__ );
  }
  }
  MPI_Finalize();        /* free parameters info */

  exit (0);
 }


