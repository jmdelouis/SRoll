// this is for srand48(), drand48(), realpath()
#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE 700
#endif

#define MAXPIXDEF (9)

// MAP NAME DEFINITION
#define MAX_OUT_NAME_LENGTH (2048)

#define CNN_NSIDE (32)

//should not compute calcmatrix
#if 0
#define CALCMATRIX
#endif

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
int *rank_map;

float *tpparam=NULL; // parameters for neural network initialise to NULL at the beginning.
int CNN_NB_PARAM=0;

int verbose=0;
PIOLONG Nside;
long NORM_GAIN=0;
long REMOVE_CAL=0;

int NUMBEROFITER=500;
double LIMIT_ITER=1E-30;
int do_offset=DOOFFSET;
int NORMFITPOL=0;
long RINGSIZE=27664; // default value used by Planck HFI

PyObject *sparseFunc=NULL;
PyObject *CorrTODFunc=NULL;
PIOINT TestNormalizeSparse=0;

int PIOWriteVECT(const char *path,void *value,int off,int size);
double det(double *mat, int n);

long *realpix;
int rank_zero=0;

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
  int ipx;
  float hit;
} hpfloat;

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

void sum_double_array(PyObject *pArray, PIODOUBLE *array,int maxsize);

int check_dir(const char *dir) {
    struct stat st;

    // Check if the directory exists
    if (stat(dir, &st) == 0) {
        // Check if it's actually a directory
        if (S_ISDIR(st.st_mode)) {
            return 0;  // Directory exists
        } else {
            fprintf(stderr, "%s exists but is not a directory.\n", dir);
            return -1;  // Path exists but is not a directory
        }
    }

    // If the directory does not exist, create it
    if (mkdir(dir, 0755) == -1) {  // 0755 is the standard permission for directories
        fprintf(stderr, "Error creating directory %s: %s\n", dir, strerror(errno));
        return -1;  // Return -1 in case of an error
    }

    return 0;  // Directory successfully created
}

int compar_int(const void *a, const void *b)
{
  hpint *pa = (hpint *) a;
  hpint *pb = (hpint *) b;
  return(pb->hit-pa->hit);
}

int compar_float(const void *a, const void *b)
{
  hpfloat *pa = (hpfloat *) a;
  hpfloat *pb = (hpfloat *) b;
  return((int) (pb->hit-pa->hit));
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
  if (rank==rank_zero) {
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

  if (rank==rank_zero) {
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

#define BADCOND (1E30)
// ----------------------------------------------------------------------------------------------------- 
double cond_thres(double * mat,double *res,int maxchannels){
 
  double determinant =0.0;
  memset(res,0,sizeof(maxchannels*maxchannels*sizeof(double)));
  determinant = det(mat,maxchannels);    

  if (fabs(determinant)<1E-30) return(BADCOND);
  
  double cond=0.0,cond2=0.0;
  
  if (maxchannels>1) {
    invert(mat,res,maxchannels);
    if (isnan(res[0])) return(BADCOND);
    cond=norm(mat,maxchannels)*norm(res,maxchannels);
  }else {
    res[0]=1/mat[0];
    cond=fabs(res[0]);
    return(cond);
  }
  
  //return(cond);
    
  if (cond<1E10) {
    cond2=norm(res,maxchannels); 
  }
  if (cond<=0.0) return(BADCOND);
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

  if (invert_3_3(mat,res)==-1) return(BADCOND);

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
#ifndef MAXCHAN
#define MAXCHAN (16)
#endif
#ifndef MAXEXTERNAL 
#define MAXEXTERNAL (4)
#endif

typedef struct {
  PIOFLOAT sig;
  PIOFLOAT listp[MAXSIMU];
  PIOINT   nShpr;
  PIOFLOAT listofShpr[MAXEXTERNALSHPR];
  PIOINT   listofShpr_idx[MAXEXTERNALSHPR];
  PIOLONG  ipix;
  PIOINT   rg;
  PIOINT   mpi_rank;
  PIOFLOAT corr_cnn;
  PIOINT   gi;
  PIOFLOAT hpr_cal;
  PIOFLOAT Sub_HPR;
  PIOFLOAT w;
  PIOBYTE  surv;
  PIOBYTE  ib;
  PIOFLOAT hit;
  PIOFLOAT model;
  PIOFLOAT inc;
  PIOFLOAT External[MAXEXTERNAL];
  PIOFLOAT channels[MAXCHAN];
  PIOINT irank;
  //add for new version
  int diag_idx;
  double m;
  double R_ij;
  double alpha[MAXCHAN];
} hpix;

hpix ** ptr_l_hpix;

int compar_hpix_rank(const void *a, const void *b)
{
  hpix *pa = (hpix *) a;
  hpix *pb = (hpix *) b;
  return(pa->mpi_rank-pb->mpi_rank);
}

int compar_hpix_ipix(const void *a, const void *b)
{
  hpix *pa = (hpix *) a;
  hpix *pb = (hpix *) b;
  return(pa->ipix-pb->ipix);
}

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


  if (rank==rank_zero) fprintf (stderr, "iter = %d - delta0 = %lg - delta_new = %lg\n", iter, delta0, delta_new);

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
          if (rank==rank_zero&&iter%10==0) fprintf (stderr,"gcmat_mpi() iter = %d - delta0 = %lg - delta_new = %lg\n", iter, delta0, delta_new);
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
  if (rank==rank_zero) fprintf (stderr,"gcmat_mpi2() iter = %d - delta0 = %lg - delta_new = %lg\n",
          iter, delta0, delta_new);
  if (rank==rank_zero) fprintf (stderr,"CG in iter = %d (max=%d)\n", iter, itermax);
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

hpix *com_hpix;
int memdisk=0;
hpix *call_hpix_buffer(long rank_buffer,int rank);
hpix *free_hpix_buffer(long rank_buffer,int rank);
void save_hpix_buffer(long rank_buffer,hpix *ptr,long nval,int rank);

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
      if (rank==rank_zero) fprintf(stderr,"END ABS fret[%ld] = %.10lg %.10lg %lg %lg %lg %lg : %lg %lg %lg\n",(long) its,fp,*fret,
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

    if (rank==rank_zero) fprintf(stderr,"fret[%ld] %lg =\t%.10lg\t%.10lg\t%lg\t[ %lg %lg %lg %lg ]\n",(long) its,diffp,fp,*fret,fp-*fret,
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
long npixShpr=0;
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
long l_rank_ptr_hpix=0;

#ifdef UPDATE_DIP
#define UPDATE_DIP
#endif

// -------------------------------------------------------------------------------------------------------------
void proj_data(double *b2,double *hit2,int nnbpix,int rank,double nmatres,int GAINSTEP2){

  long ir,irt;  

  double *sum_channels =(double *) malloc(sizeof(double)*nnbpix*MAXCHANNELS);
  memset(sum_channels,0,sizeof(double)*nnbpix*MAXCHANNELS);

  for(int ibuffer=0;ibuffer<l_rank_ptr_hpix;ibuffer++) {
    hpix *l_htmp=call_hpix_buffer(ibuffer,rank);
    for(int pix =0;pix<loc_nhpix[ibuffer];pix++){
      hpix *htmp = l_htmp+pix; 

      if (flgpix[htmp->ipix]>0) {
	// calcul sum_channels
        long ri1=htmp->rg-globalBeginRing;
        if (flg_rg[htmp->ib][ri1]!=0) {
	  for(int k =0;k<MAXCHANNELS;k++){
	    sum_channels[htmp->ipix*MAXCHANNELS+k] += htmp->channels[k]* htmp->R_ij;
	  }
	}
      }
    }
  }

  for(int ibuffer=0;ibuffer<l_rank_ptr_hpix;ibuffer++) {
    hpix *l_htmp=call_hpix_buffer(ibuffer,rank);
    for(int pix =0;pix<loc_nhpix[ibuffer];pix++){
      hpix *htmp = l_htmp+pix; 
      
      if (flgpix[htmp->ipix]>0) {
	double tmp=0;
	//calcul sum_Rij && calcul sum_channels
	
        long ri1=htmp->rg-globalBeginRing;
        if (flg_rg[htmp->ib][ri1]!=0) {
	  ir=rgord[htmp->ib][ri1]+newnr[htmp->ib];
	  
	  tmp = htmp->R_ij;                            
	  for(int k =0;k<MAXCHANNELS;k++){
	    tmp-= sum_channels[k+htmp->ipix*MAXCHANNELS]*htmp->alpha[k] ;
	  }
	  
	  if (do_offset==1) {
	    b2[ir]+=tmp;
	  }
	  
	  for(int m=0;m<htmp->nShpr;m++){
	    b2[newnr[nbolo]+htmp->listofShpr_idx[m]+(GAINSTEP2)*nbolo]+=htmp->listofShpr[m]*tmp;
	  }
	  
	  
	  if(GAINSTEP2 != 0){
	    b2[newnr[nbolo]+htmp->gi+htmp->ib*GAINSTEP2]+= tmp * htmp->model;
	  }
	  
	  if (do_offset==1) {
	    ir=rgord[htmp->ib][ri1]+newnr[htmp->ib];      
	    hit2[ir]+= htmp->w; // poids stat de chaque valeurs 
	  }
	  for(int m = 0;m<htmp->nShpr;m++){
	    irt = newnr[nbolo]+htmp->listofShpr_idx[m]+nbolo*(GAINSTEP2);
	    hit2[irt]+= htmp->listofShpr[m]*htmp->listofShpr[m]*htmp->w; // poids stat de chaque valeurs 
	  }
	  
	  if(GAINSTEP2 != 0) {
	    irt=newnr[nbolo]+htmp->gi+htmp->ib*GAINSTEP2;
	    hit2[irt]+=htmp->w*htmp->model*htmp->model;
	  }
	}
      }
    }
  }
  
  free(sum_channels);
}

// -------------------------------------------------------------------------------------------------------------
void proj_grad(double * q2,double nmatres,double * x,int nnbpix,int rank,int GAINSTEP2){
  //plap
  double sum =0.0;
  long ir,irt;
  double *csum = (double *) malloc(sizeof(double)*MAXCHANNELS);
  double *s_X = (double *) malloc(sizeof(double)*nnbpix*MAXCHANNELS);
  double *sum_channels = (double *) malloc(sizeof(double)*nnbpix*MAXCHANNELS);

  memset(csum,0,sizeof(double)*MAXCHANNELS);
  memset(s_X,0,sizeof(double)*MAXCHANNELS*nnbpix);
  memset(sum_channels,0,sizeof(double)*MAXCHANNELS*nnbpix);

  for(int ibuffer=0;ibuffer<l_rank_ptr_hpix;ibuffer++) {
    hpix *l_htmp=call_hpix_buffer(ibuffer,rank);
    for(int pix =0;pix<loc_nhpix[ibuffer];pix++){
      hpix *htmp = l_htmp+pix; 

      if (flgpix[htmp->ipix]>0) {
	// calcul sum_channels
        long ri1=htmp->rg-globalBeginRing;
        if (flg_rg[htmp->ib][ri1]!=0) {
          ir=rgord[htmp->ib][ri1]+newnr[htmp->ib];
	  double val=0.0;
	  if (do_offset==1) {
	    val = x[ir];
          }
	  
          for(int m = 0;m<htmp->nShpr;m++){
            irt = newnr[nbolo]+htmp->listofShpr_idx[m]+nbolo*(GAINSTEP2);
            val+= x[irt]*htmp->listofShpr[m];            
          }
          
          if(GAINSTEP2 != 0)  val+= x[newnr[nbolo]+htmp->gi+htmp->ib*GAINSTEP2]*htmp->model;

          for(int l =0;l<MAXCHANNELS;l++){ 
	    if (do_offset==1) {
	      csum[l]+=x[ir]*htmp->channels[l];
	    }
            s_X[l+htmp->ipix*MAXCHANNELS]+= htmp->w *htmp->channels[l]*val;
          }  
        }
      }
    }
  }

  for(int pix =0;pix<nnbpix;pix++){
    if (flgpix[pix]>0) {
      invertMatrix(imatrice+pix*MAXCHANNELS*MAXCHANNELS,s_X+pix*MAXCHANNELS,MAXCHANNELS,rank); 
    }
  }

  for(int ibuffer=0;ibuffer<l_rank_ptr_hpix;ibuffer++) {
    hpix *l_htmp=call_hpix_buffer(ibuffer,rank);
    for(int pix =0;pix<loc_nhpix[ibuffer];pix++){
      hpix *htmp = l_htmp+pix; 

      if (flgpix[htmp->ipix]>0) {
	// calcul sum_channels
        long ri1=htmp->rg-globalBeginRing;
        if (flg_rg[htmp->ib][ri1]!=0) {
          ir=rgord[htmp->ib][ri1]+newnr[htmp->ib];
	  double val=0.0;
	  if (do_offset==1) {
	    val = x[ir];
          }
	  
          for(int m = 0;m<htmp->nShpr;m++){
            irt = newnr[nbolo]+htmp->listofShpr_idx[m]+nbolo*(GAINSTEP2);
            val+= x[irt]*htmp->listofShpr[m];            
          }
          
          if(GAINSTEP2 != 0)  val+= x[newnr[nbolo]+htmp->gi+htmp->ib*GAINSTEP2]*htmp->model;

          //Calcul Rij_bis
          double Rij_bis = val;
          for(int k=0;k<MAXCHANNELS;k++){
            Rij_bis -= htmp->channels[k]*s_X[k+htmp->ipix*MAXCHANNELS];
          }

          for(int k=0;k<MAXCHANNELS;k++){
            sum_channels[k+htmp->ipix*MAXCHANNELS]+= htmp->channels[k]*Rij_bis;
          }
	}
      }
    }
  }

  for(int ibuffer=0;ibuffer<l_rank_ptr_hpix;ibuffer++) {
    hpix *l_htmp=call_hpix_buffer(ibuffer,rank);
    for(int pix =0;pix<loc_nhpix[ibuffer];pix++){
      hpix *htmp = l_htmp+pix; 

      if (flgpix[htmp->ipix]>0) {
	// calcul sum_channels
        long ri1=htmp->rg-globalBeginRing;
        if (flg_rg[htmp->ib][ri1]!=0) {
          ir=rgord[htmp->ib][ri1]+newnr[htmp->ib];
	  double val=0.0;
	  if (do_offset==1) {
	    val = x[ir];
          }
	  
          for(int m = 0;m<htmp->nShpr;m++){
            irt = newnr[nbolo]+htmp->listofShpr_idx[m]+nbolo*(GAINSTEP2);
            val+= x[irt]*htmp->listofShpr[m];            
          }
          
          if(GAINSTEP2 != 0)  val+= x[newnr[nbolo]+htmp->gi+htmp->ib*GAINSTEP2]*htmp->model;

          //Calcul Rij_bis
          double Rij_bis = val;
          for(int k=0;k<MAXCHANNELS;k++){
            Rij_bis -= htmp->channels[k]*s_X[k+htmp->ipix*MAXCHANNELS];
          }
   
          double tmp = Rij_bis;
	  
          for(int k=0;k<MAXCHANNELS;k++){                
            tmp -=sum_channels[k+htmp->ipix*MAXCHANNELS]*htmp->alpha[k];
          }

	  if (do_offset==1) {
	    q2[ir]+= tmp;        //ie  ((val-s_X[k])-htmp[l1].alpha[k]*htmp[l1].channels[k]*R_ij_bis);
	  }
       
          for(int m=0;m<htmp->nShpr;m++){           
            q2[newnr[nbolo]+htmp->listofShpr_idx[m]+(GAINSTEP2)*nbolo]+= htmp->listofShpr[m]*tmp;          
          }

          if(GAINSTEP2 != 0){ 
            q2[newnr[nbolo]+htmp->gi+htmp->ib*GAINSTEP2] += htmp->model*tmp;            
          }

	}
      }
    }
  }

  if(rank == rank_zero){ // to compute normalisation on precessor with less rings
    double msum=1E4;

    if (TestNormalizeSparse==1) {
      PyObject *pX = PyList_New(npixShpr);
      
      for(int b = 0;b<npixShpr;b++){
	PyList_SetItem(pX, b, PyFloat_FromDouble(x[newnr[nbolo]+b+(GAINSTEP2)*nbolo]));
      }
      // Appel de la méthode sparse
      PyObject *pValue = PyObject_CallMethod(sparseFunc, "normalize", "(O)", pX);
      
      // Traitement de la valeur de retour
      if (pValue != NULL) {
	// Assurer que pValue est un tuple
	sum_double_array(pValue,q2+newnr[nbolo]+(GAINSTEP2)*nbolo,npixShpr);
	Py_DECREF(pValue);
      }
      else {
	fprintf(stderr,"Problem while loading the 'normalize' method of the class sparse\n");
	exit(0);
      }
      // Nettoyage des arguments
      Py_DECREF(pX);
    }
    

    for(int n = 0;n<Param->n_val_mean;n++){
      double sum2 = 0.0;
      for(int b = 0;b<npixShpr;b++){
	sum2+=(x[newnr[nbolo]+b+(GAINSTEP2)*nbolo]-Param->val_mean[n])*Param->do_mean[b+n*(npixShpr)];
      }
      for(int b = 0;b<npixShpr;b++){
	q2[newnr[nbolo]+b+(GAINSTEP2)*nbolo]+= sum2*Param->w_mean[n]*Param->do_mean[b+n*(npixShpr)];
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
      // NORMALIZE OFFSET AGAINST CHANNELS
      for (int i=0;i<MAXCHANNELS;i++) {
	for (int j=0;j<newnr[nbolo];j++) q2[j]+=csum[i];
      }   
#endif
    }
  }

  free(sum_channels);
  free(csum);
  free(s_X);
}

// --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void apply_invertMatrix(double * mat,double *vec,int n,int rank){
 

      double *tmp_vec = malloc(n*sizeof(double));
      memset(tmp_vec,0,n*sizeof(double));  
      matvec(mat,vec,tmp_vec,n);
      memcpy(vec,tmp_vec,n*sizeof(double)); 
      
      free(tmp_vec);
 
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

double tol=-1.0;

// #########################################################################################################################################################################################
void minimize_gain_tf(double *ix2,double *gaingi){
  //gaingi -> gain de chaque detecteur entré
  // xi2 -> vecteur des offsets + amplitude de templates

  // Init 
  long i,k;
  int itermax = NUMBEROFITER;
  PIOLONG GAINSTEP2;
  int rank;
  int size;

  int rank_size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank_size);
  rank=rank_size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  struct timeval tp1,tp2;
  gettimeofday(&tp1,NULL);

  GAINSTEP2=GAINSTEP;

  nmatres=newnr[nbolo]+nbolo*(GAINSTEP2)+npixShpr;


  MPI_Bcast(&nmatres, sizeof(long), MPI_BYTE, rank_zero, MPI_COMM_WORLD);

  double *x_tab = (double *) malloc(sizeof(double)*(nmatres));
  double *projX =(double *) malloc(sizeof(double)*(nmatres));
  double *new_p =(double *) malloc(sizeof(double)*(nmatres)); 
  double *new_x =(double *) malloc(sizeof(double)*(nmatres));
  double *new_r = (double *) malloc(sizeof(double)*(nmatres));
  double *res = (double *) malloc(sizeof(double)*(nmatres));


  if (rank==rank_zero) {
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
 
  double *SI = (double *) malloc(sizeof(double)*(nnbpix*MAXCHANNELS));
  memset(SI,0,sizeof(double)*(nnbpix*MAXCHANNELS));

  for(int ibuffer=0;ibuffer<l_rank_ptr_hpix;ibuffer++) {
    hpix *l_htmp=call_hpix_buffer(ibuffer,rank);
    for(int pix =0;pix<loc_nhpix[ibuffer];pix++){
      hpix *htmp = l_htmp+pix; 

      if (flgpix[htmp->ipix]>0) {
        long ri1=htmp->rg-globalBeginRing;
        if (flg_rg[htmp->ib][ri1]!=0) {
	  double g1=gaingi[htmp->gi+htmp->ib*GAINSTEP];
          double val_tmp =htmp->sig*g1- htmp->Sub_HPR-htmp->corr_cnn;
	  if (REMOVE_CAL==1) {
	    val_tmp -=  htmp->hpr_cal;
	  }
	  
          for(int l =0;l<MAXCHANNELS;l++){           
            SI[l+htmp->ipix*MAXCHANNELS]+= htmp->w * htmp->channels[l]*val_tmp;
          }

          htmp->m = val_tmp; // calcul de m 
	}
      }
    }
    save_hpix_buffer(ibuffer,l_htmp,loc_nhpix[ibuffer],rank);
  }

  // Init
  for (k=0;k<nnbpix;k++)  { 
    if (flgpix[k]>0) {
      invertMatrix(imatrice+k*MAXCHANNELS*MAXCHANNELS,SI+k*MAXCHANNELS,MAXCHANNELS,rank);
    }
  }

  for(int ibuffer=0;ibuffer<l_rank_ptr_hpix;ibuffer++) {
    hpix *l_htmp=call_hpix_buffer(ibuffer,rank);
    for(int pix =0;pix<loc_nhpix[ibuffer];pix++){
      hpix *htmp = l_htmp+pix; 

      if (flgpix[htmp->ipix]>0) {
        long ri1=htmp->rg-globalBeginRing;
        if (flg_rg[htmp->ib][ri1]!=0) {               //if pixel ok
	  htmp->R_ij = htmp->m ; //- SI;            
	  for(int i =0;i<MAXCHANNELS;i++){
	    htmp->alpha[i] = 0.0;                  
	    htmp->R_ij -= SI[i+htmp->ipix*MAXCHANNELS]*htmp->channels[i];  
	    for(int j=0;j<MAXCHANNELS;j++){ 
	      //calcul alpha
	      htmp->alpha[i] += htmp->w*imatrice[j+MAXCHANNELS*i+htmp->ipix*MAXCHANNELS*MAXCHANNELS]*htmp->channels[j];  
	    }        
	  }
	}
      }
    }
    save_hpix_buffer(ibuffer,l_htmp,loc_nhpix[ibuffer],rank);
  }
  
  free(SI);

  //Init b2 and q2 to 0
  memset(b2,0,sizeof(double)*(nmatres));
  memset(hit2,0,sizeof(double)*(nmatres));
  memset(q2,0,sizeof(double)*(nmatres));

  proj_data(b2,hit2,nnbpix,rank,nmatres,GAINSTEP2); //calcul b2
  proj_grad(q2,nmatres,x_tab,nnbpix,rank,GAINSTEP2); // Calcul q2 

  // Recuperation de b2
  {
    double *lb = (double *) malloc(sizeof(double)*(nmatres));
    MPI_Reduce(b2,lb,nmatres,MPI_DOUBLE,MPI_SUM,rank_zero,MPI_COMM_WORLD); //parallisation MPI -- recuperation b
    memcpy(b2,lb,sizeof(double)*(nmatres)); // copy b2 dans lb
    free(lb);
  }

  //recup and send hit2
  {
    double *lb = (double *) malloc(sizeof(double)*(nmatres));
    MPI_Reduce(hit2,lb,nmatres,MPI_DOUBLE,MPI_SUM,rank_zero,MPI_COMM_WORLD);
    memcpy(hit2,lb,sizeof(double)*(nmatres));
    free(lb);
  }
  // END Recuperation de b2 and hit2


  //Send q2
  {
    double *lb = (double *) malloc(sizeof(double)*(nmatres));
    MPI_Reduce(q2,lb,nmatres,MPI_DOUBLE,MPI_SUM,rank_zero,MPI_COMM_WORLD);
    memcpy(q2,lb,sizeof(double)*(nmatres));
    free(lb);
  }
  // End send q2


  if (rank==rank_zero) { 
    fprintf(stderr,"B2 %lg\n",b2[0]);
    fprintf(stderr,"H2 %lg\n",hit2[0]);
    fprintf(stderr,"Q2 %lg\n",q2[0]);
    for (int i=0;i<nmatres;i++) if (hit2[i]==0) fprintf(stderr,"H2 %d %lg\n",i,hit2[i]);
  }

  //init r2 and d2 to 0
  memset(r2,0,sizeof(double)*(nmatres)); 
  memset(d2,0,sizeof(double)*(nmatres));


  //Calcul r and p
  if (rank==rank_zero) {
    for (i=0; i < nmatres; i++){
      r2[i] = b2[i] - q2[i]; //r = b - Ax0 = Ax - Ax0  
      if (hit2[i]>0.0) d2[i] = r2[i]/hit2[i]; //d2 => p
      //if (rank==rank_zero)  fprintf(stderr,"[DEBUG] i = %d q2[] %lf b2[] =%lf r2[] = %lf ,d2[] = %lf , hit2[] = %lf \n",i,q2[i],b2[i],r2[i],d2[i],hit2[i]);

    }
  }


  MPI_Bcast(d2, sizeof(double)*nmatres, MPI_BYTE, rank_zero, MPI_COMM_WORLD);
  MPI_Bcast(r2, sizeof(double)*nmatres, MPI_BYTE, rank_zero, MPI_COMM_WORLD);
  MPI_Bcast(hit2, sizeof(double)*nmatres, MPI_BYTE, rank_zero, MPI_COMM_WORLD);  

    
  memset(new_x,0,sizeof(double)*nmatres);
  memset(new_r,0,sizeof(double)*nmatres);
  memset(new_p,0,sizeof(double)*nmatres);
  memset(res,0,sizeof(double)*nmatres);

  
  //init
  double sum = 0.0 ,tmp = 0.0,alpha_tmp = 0.0 ,tmp_delta=0.0,beta= 0.0,delta = 0.0,alpha=0.0,tot_time = 0.0,time_exc =0.0;
    
  
  if (itbogo==0) delta0 = 0;
#if 1
  int n=0;
  while(n<itermax) {
      
      gettimeofday(&tp1,NULL);
      
      //calcul projX
      memset(projX,0,sizeof(double)*nmatres);
      proj_grad(projX,nmatres,d2,nnbpix,rank,GAINSTEP2);

      //send projX
      double *lb = (double *) malloc(sizeof(double)*(nmatres));
      memset(lb,0,sizeof(double)*nmatres);
      MPI_Reduce(projX,lb,nmatres,MPI_DOUBLE,MPI_SUM,rank_zero,MPI_COMM_WORLD);
      memcpy(projX,lb,sizeof(double)*(nmatres));
      free(lb);

      MPI_Bcast(projX, sizeof(double)*nmatres, MPI_BYTE, rank_zero, MPI_COMM_WORLD);

      //Calcul delta
      tmp = 0.0;
      tmp_delta = 0.0;
      for(int k = 0;k < nmatres;k++){
        tmp_delta += r2[k]*r2[k];
      }
      
      delta = tmp_delta;
      
      if (n==0) {
	delta0 = delta;
	if (tol==-1) tol=LIMIT_ITER*delta0;

	if(rank==rank_zero) fprintf(stderr,"\n---------\n");
	if(rank==rank_zero) fprintf(stderr,"\n D0=%lg LIMIT=%lg\n",delta0,tol);
	if(rank==rank_zero) fprintf(stderr,"\n---------\n");
      }
      //calcul alpha
      alpha_tmp = 0.0;
      tmp = 0.0;
      
      for(int k = 0;k < nmatres;k++){
        alpha_tmp+= d2[k]*projX[k];
        if (hit2[k]>0.0) tmp += (r2[k]/hit2[k])*r2[k];
      }
      //alpha = delta /alpha_tmp;
      alpha = tmp /alpha_tmp;
      // manage case where alpha_tmp=0
      if (isnormal(alpha)) {
	      MPI_Bcast(&alpha, sizeof(double), MPI_BYTE, rank_zero, MPI_COMM_WORLD);
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

        if(rank == rank_zero) fprintf(stderr,"\n==> End delta = %lg time = %3lfs\n",delta,tot_time);
        //for(int i =0;i<nmatres;i++) res[i] = new_x[i];

        memcpy(ix2,new_x,sizeof(double)*(nmatres));
        MPI_Bcast(ix2, sizeof(double)*(nmatres), MPI_BYTE, rank_zero, MPI_COMM_WORLD);
        
        itbogo ++;       
        
        return;     
      }

      //calcul beta
      tmp = 0.0;
      sum =0.0;
      for(int i =0;i<nmatres;i++){
        if (hit2[i]>0.0) { //>1 ?????
	  tmp += (new_r[i]/hit2[i])*new_r[i];
	  sum += (r2[i]/hit2[i])*r2[i];
	}
      }

      beta =tmp/sum;
      
    
      MPI_Bcast(&beta, sizeof(double), MPI_BYTE, rank_zero, MPI_COMM_WORLD);

      //calcul new_p
      for(int k =0;k<nmatres;k++){
        if (hit2[k]>0.0) new_p[k] = (new_r[k]/hit2[k])+beta*d2[k];
      }

      //mise a jour r,p,x
      memcpy(r2,new_r,sizeof(double)*(nmatres));
      memcpy(d2,new_p,sizeof(double)*(nmatres));
      memcpy(x_tab,new_x,sizeof(double)*(nmatres));

      gettimeofday(&tp2,NULL);
      time_exc = (double)(tp2.tv_sec-tp1.tv_sec)+(1E-6)*(tp2.tv_usec-tp1.tv_usec);
      tot_time += time_exc;
      if (rank==rank_zero&&n%10==0) fprintf(stderr,"iter: %d/%d beta = %12lg  alpha = %12lg  delta = %12lg  %12lfs\n",n,itermax,beta,alpha,delta,time_exc);
      n=n+1;
  }
#endif
  if(rank==rank_zero) fprintf(stderr,"tot_time = %lg\n",tot_time);

  itbogo ++;
  memcpy(ix2,new_x,sizeof(double)*(nmatres));
  MPI_Bcast(ix2, sizeof(double)*(nmatres), MPI_BYTE, rank_zero, MPI_COMM_WORLD);
      

}

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

PIOLONG last_call_ptr=-1;

hpix *call_hpix_buffer(long rank_buffer,int rank)
{
  
  if (memdisk) {
    if (ptr_l_hpix[rank_buffer]==NULL) {
      if (last_call_ptr!=-1&&ptr_l_hpix[last_call_ptr]!=NULL) {
	free(ptr_l_hpix[last_call_ptr]);
	ptr_l_hpix[last_call_ptr]=NULL;
      }
      char fname[2048];
      sprintf(fname,"%s/TEMP_%d_%ld.dat",Param->TEMP_DISK,rank,rank_buffer);
      
      // Ouvrir le fichier en mode binaire
      FILE* file = fopen(fname, "rb");
      if (!file) {
	fprintf(stderr,"File : %s\n",fname);
        perror("Erreur lors de l'ouverture du fichier");
        return NULL;
      }
      
      // Aller à la fin du fichier pour connaître sa taille
      fseek(file, 0, SEEK_END);
      long file_size = ftell(file);
      rewind(file);
      
      // Calculer le nombre d'éléments de type double
      long nval = file_size / sizeof(hpix);
      
      // Allouer de la mémoire pour le buffer
      hpix* buffer = (hpix*) malloc(sizeof(hpix)*nval);
      if (!buffer) {
        perror("Erreur lors de l'allocation mémoire");
        fclose(file);
        return NULL;
      }
      
      // Lire le contenu du fichier dans le buffer
      size_t result = fread(buffer, sizeof(hpix), nval, file);
      if (result != nval) {
        perror("Erreur lors de la lecture du fichier");
        free(buffer);
        fclose(file);
        return NULL;
      }

      // Fermer le fichier
      fclose(file);
      loc_nhpix[rank_buffer]=nval;
      ptr_l_hpix[rank_buffer]=buffer;
      last_call_ptr=rank_buffer;
      return buffer;
    }
    else {
      last_call_ptr=rank_buffer;
      return ptr_l_hpix[rank_buffer];
    }
  }
  else {
    return ptr_l_hpix[rank_buffer];
  }
  return NULL;
}

hpix *free_hpix_buffer(long rank_buffer,int rank)
{
  if (memdisk) {
    if (ptr_l_hpix[rank_buffer]==NULL) {
      // read from disk HPTR
    }
    else {
      free(ptr_l_hpix[rank_buffer]);
    }
  }
  else {
    if (ptr_l_hpix[rank_buffer]==NULL) {
      // read from disk HPTR
    }
    else {
      free(ptr_l_hpix[rank_buffer]);
    }
  }
  return NULL;
}

void save_hpix_buffer(long rank_buffer,hpix *ptr,long nval,int rank)
{
  if (memdisk) {
    char fname[2048];
    sprintf(fname,"%s/TEMP_%d_%ld.dat",Param->TEMP_DISK,rank,rank_buffer);
    FILE *fp=fopen(fname,"wb");
    long nout=fwrite(ptr,sizeof(hpix),nval,fp);
    if (nout!=nval) {
      fprintf(stderr,"Problem while writing the memory cache in %s %ld %ld\n",fname,nout,nval);
      exit(0);
    }
    fclose(fp);
    free(ptr);
    
    loc_nhpix[rank_buffer]=nval;
    ptr_l_hpix[rank_buffer]=NULL;
    // read from disk HPTR
  }
  else {
    if (ptr!=NULL) {
      ptr_l_hpix[rank_buffer]=ptr;
      loc_nhpix[rank_buffer]=nval;
    }
    else {
      ptr_l_hpix[rank_buffer]= (hpix *) _PIOMALLOC(sizeof(hpix)*nval);
      loc_nhpix[rank_buffer]=0;
    }
  }
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
  if (mpi_rank == rank_zero) {
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
  if (rank==rank_zero)  {
    map=(float *)malloc(sizeof(float)*map_size);
  }

  if (getbeginfo==0) {
    getbeginfo=1;
    int i,rrk;
    if (rank==rank_zero) {
      allbeg = (int *) malloc(sizeof(int)*mpi_size);
      allend = (int *) malloc(sizeof(int)*mpi_size);
    }
    MPI_Gather(&beg,sizeof(int),MPI_BYTE,allbeg,sizeof(int),MPI_BYTE,rank_zero,MPI_COMM_WORLD);
    MPI_Gather(&end,sizeof(int),MPI_BYTE,allend,sizeof(int),MPI_BYTE,rank_zero,MPI_COMM_WORLD);

    if (rank==rank_zero) {
      for (rrk=0;rrk<mpi_size;rrk++) {
	if (maxsize<allend[rrk]-allbeg[rrk]+1) maxsize=allend[rrk]-allbeg[rrk]+1;
      }
      all_realpix = (int *) malloc(sizeof(int)*mpi_size*maxsize);
      all_map = (float *) malloc(sizeof(float)*mpi_size*maxsize);
    }

    MPI_Bcast(&maxsize,sizeof(int), MPI_BYTE, rank_zero, MPI_COMM_WORLD);

    int *l_idx =(int *)malloc(sizeof(int)*maxsize);
    for (i=beg;i<=end;i++) l_idx[i-beg]=realpix[i-beg];

    MPI_Gather(l_idx,sizeof(int)*maxsize,MPI_BYTE,all_realpix,sizeof(int)*maxsize,MPI_BYTE,rank_zero,MPI_COMM_WORLD);

    free(l_idx);
  }
  float *l_map =(float *)malloc(sizeof(float)*maxsize);
  for (k=beg;k<=end;k++) l_map[k-beg]=value[k-beg];
  MPI_Gather(l_map,sizeof(float)*maxsize,MPI_BYTE,all_map,sizeof(float)*maxsize,MPI_BYTE,rank_zero,MPI_COMM_WORLD);
  free(l_map);


  if (rank==rank_zero) {
    int i,rrk;
    for (rrk=0;rrk<mpi_size;rrk++) {
      int l_beg,l_end;
      l_beg=allbeg[rrk];
      l_end=allend[rrk];
      for (i=l_beg;i<=l_end;i++) map[all_realpix[i-l_beg+rrk*maxsize]]=all_map[i-l_beg+rrk*maxsize];
    }
  }

  if (rank==rank_zero) {
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
  if (remove(path)==0) {
    fprintf(stderr,"%s has been remove before writing data\n",path);
  }
  int fp=open(path,O_WRONLY|O_CREAT,0664);
  int err=pwrite(fp,value,size,off);
  close(fp);
  return(err);
}

/* ---------------------------------------------------------------------------------*/
int CheckMethods(PyObject *pClass,const char *name) {
  int test_method=0;
  // Récupérer la liste des attributs de la classe
  PyObject *pDir = PyObject_Dir(pClass);
  
  // Vérifier si la liste est valide
  if (pDir == NULL) {
    printf("Erreur: Impossible de récupérer la liste des attributs de la classe\n");
    return 0;
  }
  
  // Parcourir la liste des attributs
  Py_ssize_t size = PyList_Size(pDir);
  for (Py_ssize_t i = 0; i < size; ++i) {
    // Récupérer le nom de l'attribut
    PyObject *pAttrName = PyList_GetItem(pDir, i);
    const char *attrName = PyUnicode_AsUTF8(pAttrName);
    
    // Vérifier si l'attribut est une méthode
    PyObject *pAttr = PyObject_GetAttrString(pClass, attrName);
    if (pAttr != NULL && PyCallable_Check(pAttr)) {
      if (strcmp(name,attrName)==0) {
	//fprintf(stderr,"Méthode : %s\n", attrName);
	test_method=1;
      }
    }
    Py_XDECREF(pAttr);
  }
  
  // Libérer les ressources
  Py_DECREF(pDir);

  return(test_method);
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
        exit(0);
    }
  
    Py_DECREF(pArgs);

    return PyLong_AsLong(result);
}
// Fonction pour copier les données d'un tableau Python d'entiers à un tableau C
int copy_int_array(PyObject *pArray, PIOINT *array,int maxsize) {
  
  if (!PyList_Check(pArray) && !PyTuple_Check(pArray)) {
    fprintf(stderr, "L'objet n'est ni une liste ni un tuple\n");
    exit(0);
    return -1;
  }

  Py_ssize_t n = PyList_Check(pArray) ? PyList_Size(pArray) : PyTuple_Size(pArray);
  
  if (maxsize!=-1) {
    if (n>maxsize) {
      fprintf(stderr,"Table read from python (%d) is bigger than allocated memory (%d)\n",(int)n,maxsize);
      exit(0);
    }
  }
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
int copy_long_array(PyObject *pArray, PIOLONG *array,int maxsize) {
  
  if (!PyList_Check(pArray) && !PyTuple_Check(pArray)) {
    fprintf(stderr, "L'objet n'est ni une liste ni un tuple\n");
    exit(0);
    return -1;
  }

  Py_ssize_t n = PyList_Check(pArray) ? PyList_Size(pArray) : PyTuple_Size(pArray);
  
  if (maxsize!=-1) {
    if (n>maxsize) {
      fprintf(stderr,"Table read from python (%d) is bigger than allocated memory (%d)\n",(int)n,maxsize);
      exit(0);
    }
  }
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
void sum_double_array(PyObject *pArray, PIODOUBLE *array,int maxsize) {
  
  if (!PyList_Check(pArray) && !PyTuple_Check(pArray)) {
    fprintf(stderr, "The provided object is not a list or a tuple\n");
    exit(0);
  }

  Py_ssize_t n = PyList_Check(pArray) ? PyList_Size(pArray) : PyTuple_Size(pArray);

  if (maxsize!=-1) {
    if (n>maxsize) {
      fprintf(stderr,"Table read from python (%d) is bigger than allocated memory (%d)\n",(int)n,maxsize);
      exit(0);
    }
  }
  for (Py_ssize_t i = 0; i < n; i++) {
    PyObject *item = PyList_Check(pArray) ? PyList_GetItem(pArray, i) : PyTuple_GetItem(pArray, i);
    if (!PyFloat_Check(item)) {
      fprintf(stderr, "L'élément n'est pas un double %d\n",__LINE__);
      exit(0);
    }
    array[i] += PyFloat_AsDouble(item);
  }

}
// Fonction pour copier les données d'un tableau Python d'entiers à un tableau C
int copy_float_array(PyObject *pArray, PIOFLOAT *array,int maxsize) {
  
  if (!PyList_Check(pArray) && !PyTuple_Check(pArray)) {
    fprintf(stderr, "The provided object is not a list or a tuple\n");
    return -1;
  }

  Py_ssize_t n = PyList_Check(pArray) ? PyList_Size(pArray) : PyTuple_Size(pArray);

  if (maxsize!=-1) {
    if (n>maxsize) {
      fprintf(stderr,"Table read from python (%d) is bigger than allocated memory (%d)\n",(int)n,maxsize);
      return -1;
    }
  }
  for (Py_ssize_t i = 0; i < n; i++) {
    PyObject *item = PyList_Check(pArray) ? PyList_GetItem(pArray, i) : PyTuple_GetItem(pArray, i);
    if (!PyFloat_Check(item)) {
      fprintf(stderr, "L'élément n'est pas un double %d\n",__LINE__);
      return -1;
    }
    array[i] = (PIOFLOAT) PyFloat_AsDouble(item);
  }

  return (int) n;
}

// Fonction pour copier les données d'un tableau Python d'entiers à un tableau C
int copy_double_array(PyObject *pArray, PIODOUBLE *array,int maxsize) {
  
  if (!PyList_Check(pArray) && !PyTuple_Check(pArray)) {
    fprintf(stderr, "The provided object is not a list or a tuple\n");
    return -1;
  }

  Py_ssize_t n = PyList_Check(pArray) ? PyList_Size(pArray) : PyTuple_Size(pArray);

  if (maxsize!=-1) {
    if (n>maxsize) {
      fprintf(stderr,"Table read from python (%d) is bigger than allocated memory (%d)\n",(int)n,maxsize);
      return -1;
    }
  }
  for (Py_ssize_t i = 0; i < n; i++) {
    PyObject *item = PyList_Check(pArray) ? PyList_GetItem(pArray, i) : PyTuple_GetItem(pArray, i);
    if (!PyFloat_Check(item)) {
      fprintf(stderr, "L'élément n'est pas un double %d\n",__LINE__);
      return -1;
    }
    array[i] = (PIODOUBLE) PyFloat_AsDouble(item);
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
      exit(0);
    }
    
  }
  else {
    PyErr_Print();
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
      exit(0);
    }
    
  }
  else {
    PyErr_Print();
    exit(0);
  }
  
  return n;
}

int Get_hidx(PyObject *projFunc,double ph,double th,double psi,int idx_bolo,
	     int idx_in_ring,double rg_norm,PIOFLOAT *External,
	     PIOFLOAT *o_widx,PIOLONG *o_hidx,int rank)
{
  
  PyObject *pValue=NULL;
  int n=0;

  if (projFunc != NULL){
    PyObject *pArg1 = PyFloat_FromDouble((double)ph); // val est un flottant
    PyObject *pArg2 = PyFloat_FromDouble((double)th); // val est un flottant
    PyObject *pArg3 = PyFloat_FromDouble((double)psi); // val est un flottant
    PyObject *pArg4 = PyLong_FromLong((long)idx_bolo); // val est un flottant
    PyObject *pArg5 = PyFloat_FromDouble((double)rg_norm); // val est un flottant
    PyObject *pArg6 = PyLong_FromLong((long)idx_in_ring); // val est un flottant

    PyObject *pList = PyList_New(NB_EXTERNAL);
    for (int i = 0; i < NB_EXTERNAL; i++) {
      PyList_SetItem(pList, i, PyFloat_FromDouble((double) External[i*RINGSIZE]));
    }

    pValue = PyObject_CallMethod(projFunc, "get_healpix_idx", "(OOOOOOO)",
				 pArg1, pArg2, pArg3,
				 pArg4, pArg5, pArg6,
				 pList);

    // Traitement de la valeur de retour
    if (pValue != NULL) {
      if (PyTuple_Check(pValue) && PyTuple_Size(pValue) == 2) {
	int n1=copy_long_array(PyTuple_GetItem(pValue, 0), o_hidx,MAXPIXDEF);
	int n2=copy_float_array(PyTuple_GetItem(pValue, 1), o_widx,MAXPIXDEF);
	if (n1!=n2) {
	  fprintf(stderr, "get_heapix_idx function does not return the same number of values %d %d\n",n1,n2);
	}
	n=n1;
	Py_DECREF(pValue);
      }
      else {
	fprintf(stderr,"%d %d\n",(int) PyTuple_Check(pValue),(int) PyTuple_Size(pValue));
	fprintf(stderr,"Problem while executing the get_heapix_idx method get value inside projection class\n");
	PyErr_Print();
	exit(0);
      }
    }
    else {
      fprintf(stderr,"%d %lf %lf %lf %d %lf %d %lf\n",rank,
	      ph,th,psi,(int) idx_bolo,rg_norm,(int)idx_in_ring,External[RINGSIZE]);
      fprintf(stderr, "Problem while trying to compute the healpix coordinate in the projection class %d\n",rank);
      PyErr_Print();
      exit(0);
    }

    // Nettoyage des arguments
    Py_DECREF(pArg1);
    Py_DECREF(pArg2);
    Py_DECREF(pArg3);
    Py_DECREF(pArg4);
    Py_DECREF(pArg5);
    Py_DECREF(pArg6);
    Py_DECREF(pList);
    
  } else {
    PyErr_Print();
    exit(0);
  }
  
  return n;
}

int init_channels(hpix * h,PyObject *projFunc,double psi,PIOFLOAT *External,double rgnorm,long ipix,int idx_bolo,int idx_in_ring,PIOFLOAT *sig,PIOFLOAT *calib,PIOFLOAT *hit,int rank)
{
  
  PyObject *pValue=NULL;
  int is_valid=0;
  
  if (projFunc != NULL){
    PyObject *pArg1 = PyFloat_FromDouble((double)psi); // val est un flottant
    PyObject *pList = PyList_New(NB_EXTERNAL);
    
    for (int i = 0; i < NB_EXTERNAL; i++) {
      PyList_SetItem(pList, i, PyFloat_FromDouble((double) External[i*RINGSIZE]));
    }

    PyObject *pArg3 = PyFloat_FromDouble((double)rgnorm); // val est un flottant
    PyObject *pArg4 = PyLong_FromLong((long) ipix); // val est un flottant
    PyObject *pArg5 = PyLong_FromLong((long)idx_bolo); // val est un flottant
    PyObject *pArg6 = PyFloat_FromDouble((double) *hit); // val est un flottant
    PyObject *pArg7 = PyLong_FromLong((long)idx_in_ring); // val est un flottant
    PyObject *pArg8 = PyFloat_FromDouble((double) *sig); // val est un flottant
    PyObject *pArg9 = PyFloat_FromDouble((double) *calib); // val est un flottant
    
    pValue = PyObject_CallMethod(projFunc, "eval", "(OOOOOOOOO)",
				 pArg1, pList, pArg3,
				 pArg4, pArg5, pArg6,
				 pArg7, pArg8, pArg9);

    // Traitement de la valeur de retour
    if (pValue != NULL) {
      if (PyTuple_Check(pValue) && PyTuple_Size(pValue) == 5) {
	// Extraire les deux tableaux du tuple
	is_valid = (int) PyLong_AsLong(PyTuple_GetItem(pValue, 0));
	*sig = (double) PyFloat_AsDouble(PyTuple_GetItem(pValue, 1));
	*hit = (double) PyFloat_AsDouble(PyTuple_GetItem(pValue, 2));
	*calib = (double) PyFloat_AsDouble(PyTuple_GetItem(pValue, 3));
	int err=copy_float_array(PyTuple_GetItem(pValue, 4),h->channels,MAXCHAN);
	if (err!=MAXCHANNELS) {
	  fprintf(stderr, "Projection function does not provide the good number of channels: expected %d  get %d\n",MAXCHANNELS,err);
	  exit(0);
	}
      }
      else {
	fprintf(stderr,"double psi %lf\n",psi);
	for (int k=0;k<NB_EXTERNAL;k++)
	  fprintf(stderr,"double External[%d] %f\n",k,External[k*RINGSIZE]);
	fprintf(stderr,"double rgnorm %lf\n",rgnorm);
	fprintf(stderr,"long ipix %ld\n",ipix);
	fprintf(stderr,"int idx_bolo %d\n",idx_bolo);
	fprintf(stderr,"int idx_in_ring %d\n",idx_in_ring);
	fprintf(stderr, "Projection function does not provide the good number of argument (4): expected is_valid,sig,hit,chan_value\n");
	fprintf(stderr, " Where sig and hit are a double (respectively signal and hit count) and\n");
	fprintf(stderr, " chan_value is a list of double values with len(chan_value)=%d",MAXCHANNELS);
	exit(0);
      }
      Py_DECREF(pValue);
    }
    else {
      fprintf(stderr,"double psi %lf\n",psi);
      for (int k=0;k<NB_EXTERNAL;k++)
	fprintf(stderr,"double External[%d] %f\n",k,External[k*RINGSIZE]);
      fprintf(stderr,"double rgnorm %lf\n",rgnorm);
      fprintf(stderr,"long ipix %ld\n",ipix);
      fprintf(stderr,"int idx_bolo %d\n",idx_bolo);
      fprintf(stderr,"int idx_in_ring %d\n",idx_in_ring);
      fprintf(stderr, "Problem while trying to compute the projection %d %d\n",rank,__LINE__);
      PyErr_Print();
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
    Py_DECREF(pArg9);
    
  } else {
    PyErr_Print();
    exit(0);
  }
  
  return is_valid;
  
}

int npar_corrtod(PyObject *corrtodFunc,
		 PIOINT HealIdx)
{
  PyObject *pValue=NULL;
  int n=0;
  if (corrtodFunc != NULL){
    // Créer des arguments pour la fonction Python
    // Création des arguments

    PyObject *pyHidx   = PyLong_FromLong((long) HealIdx);
    
    // Appel de la méthode 'eval's
    pValue = PyObject_CallMethod(CorrTODFunc, "npar_correction", "(O)",pyHidx);

    // Traitement de la valeur de retour
    if (pValue != NULL) {
      n=(int) PyLong_AsLong(pValue);
      Py_DECREF(pValue);
    }
    else {
      fprintf(stderr,"Problem while evaluating the correction class\n");
      PyErr_Print(); 
      exit(0);
    }

    // Nettoyage des arguments
    Py_DECREF(pyHidx);
  } else {
    PyErr_Print();
    exit(0);
  }

  return n;
}

void eval_corrtod(PyObject *corrtodFunc,
		 PIODOUBLE inc_ref,
		 PIODOUBLE rg_ref,
		 PIOINT HealIdx,
		 double *delta)
{
  PyObject *pValue=NULL;
  
  if (corrtodFunc != NULL){
    // Créer des arguments pour la fonction Python
    // Création des arguments

    PyObject *pyInc_ref= PyFloat_FromDouble((double) inc_ref);
    PyObject *pyRg_ref = PyFloat_FromDouble((double) rg_ref);
    PyObject *pyHidx   = PyLong_FromLong((long) HealIdx);
    
    // Appel de la méthode 'eval's
    pValue = PyObject_CallMethod(CorrTODFunc, "eval_correction", "(OOO)",
				 pyInc_ref,pyRg_ref,pyHidx);

    // Traitement de la valeur de retour
    if (pValue != NULL) {
      // Assurer que pValue est un tuple
      int err=copy_double_array(pValue,delta,MAXCHANNELS);
      if (err!=MAXCHANNELS) {
	fprintf(stderr,"Problem while evaluating the correction class (did not get the proper number of data) %d instead of %d\n",
		(int) err,(int) MAXCHANNELS);
	PyErr_Print(); 
	exit(0);
      }
      Py_DECREF(pValue);
    }
    else {
      fprintf(stderr,"Problem while evaluating the correction class\n");
      PyErr_Print(); 
      exit(0);
    }

    // Nettoyage des arguments
    Py_DECREF(pyInc_ref);
    Py_DECREF(pyRg_ref);
    Py_DECREF(pyHidx);
  } else {
    PyErr_Print();
    exit(0);
  }
}

int calc_corrtod_hpr(PyObject *corrtodFunc,
		     PIODOUBLE *Signal,
		     PIODOUBLE *Hit,
		     PIODOUBLE *Inc,
		     PIOINT *rg,
		     PIOINT *ib,
		     PIOINT hidx,
		     PIODOUBLE *External,
		     PIODOUBLE *o_signal,
		     PIODOUBLE *o_hit,
		     PIOLONG n_values)
{
  PyObject *pValue=NULL;

  if (corrtodFunc != NULL){
    // Créer des arguments pour la fonction Python
    // Création des arguments

    PyObject *pSignal   = PyList_New(n_values);
    PyObject *pHit      = PyList_New(n_values);
    PyObject *pInc      = PyList_New(n_values);
    PyObject *prg       = PyList_New(n_values);
    PyObject *pib       = PyList_New(n_values);
    PyObject *phidx     = PyList_New(n_values);
    PyObject *pExternal = PyList_New(n_values*NB_EXTERNAL);

    for (int i = 0; i < n_values; i++) {
      PyList_SetItem(pSignal, i, PyFloat_FromDouble((double) Signal[i]));
      PyList_SetItem(pHit, i, PyFloat_FromDouble((double) Hit[i]));
      PyList_SetItem(pInc, i, PyFloat_FromDouble((double) Inc[i]));
      PyList_SetItem(prg, i, PyLong_FromLong((long) rg[i]));
      PyList_SetItem(pib, i, PyLong_FromLong((long) ib[i]));
      PyList_SetItem(phidx, i, PyLong_FromLong((long) hidx));
      for (int j = 0; j < NB_EXTERNAL; j++) {
	PyList_SetItem(pExternal, i*NB_EXTERNAL+j, PyFloat_FromDouble((double) External[i*NB_EXTERNAL+j]));
      }
    }

    // Appel de la méthode 'eval's
    pValue = PyObject_CallMethod(CorrTODFunc, "eval", "(OOOOOOO)",pSignal,pHit,pInc,prg,pib,phidx,pExternal);

    // Traitement de la valeur de retour
    if (pValue != NULL) {
      // Assurer que pValue est un tuple
      if (PyTuple_Check(pValue) && PyTuple_Size(pValue) == 2) {
	// Extraire les deux tableaux du tuple
	int n1=copy_double_array(PyTuple_GetItem(pValue, 0), o_signal,n_values);
	int n2=copy_double_array(PyTuple_GetItem(pValue, 1), o_hit,n_values);
	if (n1!=n_values||n2!=n_values) {
	  fprintf(stderr, "CorrTOD function should provide an equal number of invex and value, here Sroll received %d %d expected %ld\n",n1,n2,n_values);
	}
	Py_DECREF(pValue);
      }
      else {
	fprintf(stderr,"Problem 1 while executing the method get value inside CorrTODFunc class\n");
	PyErr_Print(); 
	exit(0);
      }
    }
    else {
      fprintf(stderr,"Problem 2 while executing the method get value inside CorrTODFunc class\n");
      PyErr_Print(); 
      exit(0);
    }

    // Nettoyage des arguments
    Py_DECREF(pSignal);
    Py_DECREF(pHit);
    Py_DECREF(pInc);
    Py_DECREF(prg);
    Py_DECREF(pib);
    Py_DECREF(phidx);
    Py_DECREF(pExternal);
  } else {
    PyErr_Print();
    exit(0);
  }
  
  return n_values;
}



int calc_sparse_hpr(PyObject *sparseFunc,
		    long rg,
		    long ib,
		    long hpix,
		    long idx,
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
    PyObject *pArg4 = PyLong_FromLong((long)idx); // idx est un entier
    PyObject *pArg5 = PyFloat_FromDouble((double)psi); // val est un flottant
    PyObject *pList = PyList_New(NB_EXTERNAL);
    
    for (int i = 0; i < NB_EXTERNAL; i++) {
      PyList_SetItem(pList, i, PyFloat_FromDouble((double) External[i*RINGSIZE]));
    }

    // Appel de la méthode 'eval's
    pValue = PyObject_CallMethod(sparseFunc, "eval", "(OOOOOO)", pArg1, pArg2, pArg3, pArg4,  pArg5, pList);

    // Traitement de la valeur de retour
    if (pValue != NULL) {
      // Assurer que pValue est un tuple
      if (PyTuple_Check(pValue) && PyTuple_Size(pValue) == 2) {
	// Extraire les deux tableaux du tuple
	int n1=copy_int_array(PyTuple_GetItem(pValue, 0), oval_idx,MAXEXTERNALSHPR);
	int n2=copy_float_array(PyTuple_GetItem(pValue, 1), oval,MAXEXTERNALSHPR);
	if (n1!=n2) {
	  fprintf(stderr, "Sparse function should provide an equal number of invex and value, here Sroll received %d %d\n",n1,n2);
	}
	n=n1;
	Py_DECREF(pValue);
      }
      else {
	fprintf(stderr,"Problem while executing the method get value inside SparseFunc class\n");
	PyErr_Print(); 
	exit(0);
      }
    }

    // Nettoyage des arguments
    Py_DECREF(pArg1);
    Py_DECREF(pArg2);
    Py_DECREF(pArg3);
    Py_DECREF(pArg4);
    Py_DECREF(pArg5);
    Py_DECREF(pList);
    
  } else {
    PyErr_Print();
    exit(0);
  }
  
  return n;
}


int calc_diag_hpr(PyObject *diagFunc,
		  long rg,
		  long ib,
		  long hpix,
		  long idx,
		  double psi,
		  double sig,
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
    PyObject *pArg4 = PyLong_FromLong((long)idx); // idx est un entier
    PyObject *pArg5 = PyFloat_FromDouble((double)psi); // val est un flottant
    PyObject *pArg6 = PyFloat_FromDouble((double)sig); // val est un flottant
    PyObject *pList = PyList_New(NB_EXTERNAL);
    
    for (int j= 0; j < NB_EXTERNAL; j++) {
      PyList_SetItem(pList, j, PyFloat_FromDouble((double) External[j*RINGSIZE]));
    }

    // Appel de la méthode 'eval's
    pValue = PyObject_CallMethod(diagFunc, "get_diag_idx", "(OOOOOOO)", pArg1, pArg2, pArg3, pArg4, pArg5, pArg6, pList);

    
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
      exit(0);
    }

    // Nettoyage des arguments
    Py_DECREF(pArg1);
    Py_DECREF(pArg2);
    Py_DECREF(pArg3);
    Py_DECREF(pArg4);
    Py_DECREF(pArg5);
    Py_DECREF(pArg6);
    Py_DECREF(pList);
    
  } else {
    fprintf(stderr,"no DiagFunc class defined\n");
    PyErr_Print();
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
  
  if (rank==rank_zero) {
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
    MPI_Finalize();        /* free parameters info */
    exit(-1);
  }

  Param = &par;

  long MAXMPIBUFFER=Param->MAXMPIBUFFER;
  if (Param->EndRing-Param->BeginRing+1<mpi_size) {
    if (rank==rank_zero)
      fprintf(stderr, "Number of MPI rank %d should be smaller than number pointing period %d.\n",
	    mpi_size,(int) (Param->EndRing-Param->BeginRing+1));
    MPI_Finalize();        /* free parameters info */
    exit(-1);
  }
  if (Param->flag_verbose==_PAR_TRUE) verbose=Param->verbose;
  if (verbose==1) fprintf(stderr,"%s %d %d\n",__FILE__,__LINE__,rank);
  if (Param->flag_do_offset==_PAR_TRUE) do_offset=Param->do_offset;
  
  RINGSIZE = Param->RINGSIZE;

  if (rank==rank_zero) fprintf(stderr,"RINGSIZE = %d \n",(int) (RINGSIZE));
  if (rank==rank_zero) fprintf(stderr,"Nside = %d\n",(int) Param->Nside);

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
  Param->regrid=0;

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
      if (rank==rank_zero) {
        fprintf(stderr,"ERROR: too few data to be processed ("PIOLONG_FMT" rings) regarding the available ranks (%d)\n",
            globalRangeRing, mpi_size);
      }
      return 1;
    } else { // Take all available procs
      balancing_correction = globalRangeRing - (rings_per_rank * mpi_size);
    }

    if (rank==rank_zero) {
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

  if (Param->flag_TEMP_DISK==_PAR_TRUE) {
    memdisk=1;
    if (check_dir(Param->TEMP_DISK)!=0) {
      fprintf(stderr,"Impossible to create the cache directory TEMP_DISK=%s\n",
	      Param->TEMP_DISK);
      exit(0);
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

  PyObject *projFunc=NULL;
  int calc_hidx_proj=0;
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
    pArgs = PyTuple_New(4);

    PyObject *pyRank = PyLong_FromLong((long) rank); // val est un flottant
    PyObject *pyBeg = PyLong_FromLong((long) globalRankInfo.BeginRing[rank]);
    PyObject *pyEnd = PyLong_FromLong((long) globalRankInfo.EndRing[rank]);

    PyTuple_SetItem(pArgs, 0, pyParam);  // Le tuple prend la propriété de 'param'
    PyTuple_SetItem(pArgs, 1, pyRank);  // Provide the mpi rank 
    PyTuple_SetItem(pArgs, 2, pyBeg);  // Provide the Beginring of the given rank 
    PyTuple_SetItem(pArgs, 3, pyEnd);  // Provide the Endring of the given rank 

    projFunc = PyObject_CallObject(pClass,pArgs); 

    if (projFunc == NULL) {
      PyErr_Print();
      fprintf(stderr, "Échec de la création de l'instance de la classe projection\n");
      Py_DECREF(pModule);
      return 1;
    }

    if (CheckMethods(pClass,"get_healpix_idx")==1) {
      if (rank==0) fprintf(stderr,"Use specific pixel map computation from Proj class\n");
      calc_hidx_proj=1;
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

  if (rank==rank_zero)  fprintf(stderr,"Projection uses %d channels\n",MAXCHANNELS);

  GAINSTEP = Param->GAINSTEP;
  
  /*-------------------------------------------------------------------------*/
  /*   SAVE PARAMETER FILE                                                   */
  /*-------------------------------------------------------------------------*/
  if (rank==rank_zero) {
    char commandtest[PIOSTRINGMAXLEN*64];
    sprintf(commandtest,"cp %s.py %s.py",argv[1],Param->Out_VEC[0]);
    int err=system(commandtest);
    if (err) {
      fprintf(stderr,"Error while copy the parameters\n");
    }
  }
  
  if (verbose==1) fprintf(stderr,"%s %d %d\n",__FILE__,__LINE__,rank);
  
  NUMBEROFITER = Param->N_IN_ITT;
  LIMIT_ITER = Param->S_IN_ITT;

  /*-------------------------------------------------------------------------*/
  /* parameters consistency checks                                           */
  /* and default values / legacy behavior for optional parameters            */
  /*-------------------------------------------------------------------------*/

  // global number of bolometers
  nbolo = Param->n_Ptg;
  
  assert( Param->n_Sub_HPR == Param->n_SUB_HPRCOEF);

  if (verbose==1) fprintf(stderr,"%s %d %d\n",__FILE__,__LINE__,rank);
  
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

  //PIOINT Nside=Param->Nside;
  Nside = Param->Nside;
  if (Param->flag_stim_paramfiles == 1) {
    assert( Param->n_stim_paramfiles == nbolo);
  }

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
  if (rank==rank_zero) {
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
  if (rank==rank_zero) {
    fprintf(stderr,"MAXMPIBUFFER %ld\n",MAXMPIBUFFER);
    if (MAXMPIBUFFER<1024*128) {
      fprintf(stderr,"MAXMPIBUFFER is usually around %ld be sure of your choice\n",(long) (1024*1024));
    }
  }
  if (MAXMPIBUFFER==0) {
    fprintf(stderr,"MAXMPIBUFFER==0 does not work!! set a value to MAXMPIBUFFER\n");
    exit(0);
  }

  ptr_l_hpix = (hpix **) malloc(sizeof(hpix **)*MAXMPIBUFFER);
  memset(ptr_l_hpix,0,sizeof(hpix **)*MAXMPIBUFFER);
  loc_nhpix = (PIOINT *) malloc(sizeof(PIOINT)*MAXMPIBUFFER);
  memset(loc_nhpix,0,sizeof(PIOINT)*MAXMPIBUFFER);

  long alloc_ptr_hpix=MAXMPIBUFFER;
  long rank_ptr_hpix=0;
  ptr_l_hpix[rank_ptr_hpix] = (hpix *) malloc(sizeof(hpix)*MAXMPIBUFFER);
  com_hpix = ptr_l_hpix[rank_ptr_hpix];
  long n_l_hpix=0;
  

  /*======================================================
    =
    =      read data
    =
    =*/

  long vmem,phymem;
  
  PIOLONG ib;

  if (rank==rank_zero) fprintf(stderr,"Avv GAIN is equal to 0 if ==1 : NORM_GAIN : %d\n",(int) Param->NORM_GAIN);
  
  if (verbose==1) fprintf(stderr,"%s %d %d\n",__FILE__,__LINE__,rank);
  
  PIOFLOAT *addpol=NULL;
  if (Param->flag_ADDPOL==_PAR_TRUE) {
    addpol= (PIOFLOAT *) malloc(sizeof(PIOFLOAT)*12*128*128*2);
    PIOLONG nsa  = noDMC_readObject_PIOFLOAT(Param->ADDPOL,0,12*128*128*2,addpol);
    if (nsa<0) {
      fprintf(stderr, "Impossible to read ADDPOL: %s %d\n",Param->ADDPOL,(int) nsa);

      exit ( -1);
    }
  }
 
  PyObject *diagFunc=NULL;
  PyObject *pArgs = NULL;
  PyObject *pClass = NULL;

  PIOINT TestUpdateSparse=0;
  PIOINT TestCorrTODEval_corr=0;
  
  if (verbose==1) fprintf(stderr,"%s %d %d\n",__FILE__,__LINE__,rank);
  
  int TestCorrTOD=0;
  if (Param->flag_CorrTOD==_PAR_TRUE) {
    PyObject *pName = PyUnicode_FromString(argv[1]);
    PyObject *pModule = PyImport_Import(pName);
    Py_DECREF(pName);

    if (pModule == NULL) {
      PyErr_Print();
      fprintf(stderr, "Problem while creating the CorrTOD from class %s\n",argv[1]);
      return 1;
    }

    pClass = PyObject_GetAttrString(pModule, Param->CorrTOD);
    if (pClass == NULL || !PyCallable_Check(pClass)) {
      PyErr_Print();
      fprintf(stderr, "Problem while loading the classe CorrTOD\n");
      Py_XDECREF(pClass);
      Py_DECREF(pModule);
      return 1;
    }

    if (CheckMethods(pClass,"eval")==0) {
      PyErr_Print();
      fprintf(stderr, "CorrTOD class does not have an eval method\n");
      Py_DECREF(pModule);
      Py_DECREF(CorrTODFunc);
      return 1;
    }

    // check if sparse funtion need to be upgraded
    TestCorrTODEval_corr=CheckMethods(pClass,"eval_correction");
    
    // Créer un tuple pour les arguments du constructeur
    pArgs = PyTuple_New(2);

    PyObject *pyRank = PyLong_FromLong((long) rank); // val est un flottant

    PyTuple_SetItem(pArgs, 0, pyParam);  // Le tuple prend la propriété de 'param'
    PyTuple_SetItem(pArgs, 1, pyRank);  // Provide the mpi rank 

    CorrTODFunc = PyObject_CallObject(pClass,pArgs); 

    if (CorrTODFunc == NULL) {
      PyErr_Print();
      fprintf(stderr, "Problem while creating the sparse class instance\n");
      Py_DECREF(pModule);
      return 1;
    }
    TestCorrTOD=1;
  }

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

    if (CheckMethods(pClass,"eval")==0) {
      PyErr_Print();
      fprintf(stderr, "Sparse class does not have an eval method\n");
      Py_DECREF(pModule);
      Py_DECREF(sparseFunc);
      return 1;
    }

    // check if sparse funtion need to be upgraded
    TestUpdateSparse=CheckMethods(pClass,"update_eval");
    TestNormalizeSparse=CheckMethods(pClass,"normalize");
    
    // Créer un tuple pour les arguments du constructeur
    pArgs = PyTuple_New(4);

    PyObject *pyRank = PyLong_FromLong((long) rank); // val est un flottant
    PyObject *pyBeg = PyLong_FromLong((long) globalRankInfo.BeginRing[rank]);
    PyObject *pyEnd = PyLong_FromLong((long) globalRankInfo.EndRing[rank]);

    PyTuple_SetItem(pArgs, 0, pyParam);  // Le tuple prend la propriété de 'param'
    PyTuple_SetItem(pArgs, 1, pyRank);  // Provide the mpi rank 
    PyTuple_SetItem(pArgs, 2, pyBeg);  // Provide the Beginring of the given rank 
    PyTuple_SetItem(pArgs, 3, pyEnd);  // Provide the Endring of the given rank 

    sparseFunc = PyObject_CallObject(pClass,pArgs); 

    if (sparseFunc == NULL) {
      PyErr_Print();
      fprintf(stderr, "Problem while creating the sparse class instance\n");
      Py_DECREF(pModule);
      return 1;
    }

    // check update eval method to avoid problem much latter
    if (TestUpdateSparse==1) {
   
      PyObject *pIdx = PyList_New(MAXCHANNELS);
      PyObject *pWw  = PyList_New(MAXCHANNELS);
      PyObject *pChan= PyList_New(MAXCHANNELS);
      PyObject *pVal = PyList_New(MAXCHANNELS);
      
      for (int i = 0; i < MAXCHANNELS; i++) {
	PyList_SetItem(pIdx, i, PyLong_FromLong(i));
	PyList_SetItem(pWw, i, PyFloat_FromDouble(i));
	PyList_SetItem(pChan, i, PyFloat_FromDouble(i));
	PyList_SetItem(pVal, i, PyFloat_FromDouble(i));
      }
      
      // Appel de la méthode 'eval's
      PyObject *pValue = PyObject_CallMethod(sparseFunc, "update_eval", "(OOOO)", pIdx,pWw,pChan,pVal);
      
      // Traitement de la valeur de retour
      if (pValue != NULL) {
	Py_DECREF(pValue);
      }
      
      // Nettoyage des arguments
      Py_DECREF(pIdx);
      Py_DECREF(pWw);
      Py_DECREF(pChan);
      Py_DECREF(pVal);
    }
    
    if (CheckMethods(pClass,"getnumber_of_sparse")==1) {
      
      PyObject *pValue = PyObject_CallMethod(sparseFunc,"getnumber_of_sparse",NULL);

      // Traitement de la valeur de retour
      if (pValue != NULL) {
	// Assurer que pValue est un tuple
	npixShpr=(int) PyLong_AsLong(pValue);
	if (rank==rank_zero)
	  fprintf(stderr,"Number of sparse value %ld\n",(long) npixShpr);
	Py_DECREF(pValue);
      }
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
    pArgs = PyTuple_New(4);

    PyObject *pyRank = PyLong_FromLong((long) rank); // val est un flottant
    PyObject *pyBeg = PyLong_FromLong((long) globalRankInfo.BeginRing[rank]);
    PyObject *pyEnd = PyLong_FromLong((long) globalRankInfo.EndRing[rank]);

    PyTuple_SetItem(pArgs, 0, pyParam);  // Le tuple prend la propriété de 'param'
    PyTuple_SetItem(pArgs, 1, pyRank);  // Provide the mpi rank 
    PyTuple_SetItem(pArgs, 2, pyBeg);  // Provide the Beginring of the given rank 
    PyTuple_SetItem(pArgs, 3, pyEnd);  // Provide the Endring of the given rank 

    diagFunc = PyObject_CallObject(pClass,pArgs); 

    if (sparseFunc == NULL) {
      PyErr_Print();
      fprintf(stderr, "Problem while creating the class instance\n");
      Py_DECREF(pModule);
      return 1;
    }
    nb_diag=Get_NumberOfDiag(diagFunc);

    if (CheckMethods(diagFunc,"get_diag_idx")==0) {
      PyErr_Print();
      fprintf(stderr, "Diag class does not have a 'get_diag_idx' method\n");
      Py_DECREF(diagFunc);
      return 1;
    }
  }

  int *stat_pix=NULL;

  if (Param->regrid==1) {
    stat_pix = (int *) malloc(2*12*Nside*Nside*sizeof(int));
    if (rank==rank_zero) fprintf(stderr,"Sort pixel to make it run faster\n");
    memset(stat_pix,0,2*12*Nside*Nside*sizeof(int));
  }
  int number_of_iterations = 1;

  if (verbose==1) fprintf(stderr,"%s %d %d\n",__FILE__,__LINE__,rank);
  
  NB_EXTERNAL=Param->n_External/nbolo;
  if (NB_EXTERNAL==0) NB_EXTERNAL=1;

  int prtnan=0;

  for (ib=0;ib<nbolo;ib++) {

    double sxi= (Param->Calibration[ib]/Param->NEP[ib])*(Param->Calibration[ib]/Param->NEP[ib]);
    if (rank==rank_zero) fprintf(stderr,"SXI %d %lg\n",(int) ib,sxi);
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
    
    if (rank==rank_zero) fprintf(stderr,"ring_count = %d\n",ring_count);
    if (rank==rank_zero) fprintf(stderr,"bad_rings = %d\n",bad_rings);
    
    if (rank==rank_zero) fprintf(stderr,"RG_MAX %ld\n",(long) rg_max);
    if (rank==rank_zero) fprintf(stderr,"NB_DIAG %ld\n",(long) nb_diag);
    /*=========================================================================================
      Compute spline in Time:
      =========================================================================================*/

    int iter;
    PIOFLOAT *stim_hpr[MAXSIMU];
    for (iter = 0; iter < number_of_iterations; iter++) {
      stim_hpr[iter] = NULL;
    }

    if (Param->flag_stim_paramfiles == 1) {

      if (rank==rank_zero) {
        GetProcMem(&vmem,&phymem);
        fprintf(stderr,"\nbefore stim bolo loop: used VMEM %.1lf[PHYS %.1lf]MB\n",
            (double) vmem/1024./1024., (double) phymem/1024./1024.);
      }

      if (rank==rank_zero) {
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
        PIOBYTE surv=-1;
	
	if (verbose==1) fprintf(stderr,"%s %d %d\n",__FILE__,__LINE__,rank);
	
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
          fprintf(stderr, "Impossible to read Hit[%ld]: %s %d\n",ib,Param->Hit[ib],tperr);
	  exit ( -1);
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

        PIOFLOAT *External;
        if (Param->flag_External==_PAR_TRUE) {
          External = (PIOFLOAT *) _PIOMALLOC(sizeof(PIOFLOAT)*RINGSIZE*NB_EXTERNAL);
	  for (i=0;i<NB_EXTERNAL;i++) {
	    tperr = noDMC_readObject_PIOFLOAT(Param->External[i+ib*NB_EXTERNAL],rg*RINGSIZE,RINGSIZE,External+i*RINGSIZE);
	    if (tperr<0) {
	      fprintf(stderr, "Impossible to read External[%ld]: %s %d %ld\n",ib,Param->External[i],tperr,(long) rg);
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
	long nrg_htmp=0;
	
	if (verbose==1) fprintf(stderr,"%s %d %d\n",__FILE__,__LINE__,rank);
	
	hpix *l_tp_hpix = (hpix *) malloc(sizeof(hpix));
	
	for (i=0;i<RINGSIZE;i++) {
	  if (h[i]!=0.0) { 
	    long ipix;
	    int lll;
	    PIOFLOAT o_widx[MAXPIXDEF];
	    PIOLONG o_hidx[MAXPIXDEF];
	    
	    int npt=1;
	    if (calc_hidx_proj==0) {
	      ang2pix_ring( Nside, th[i], ph[i], &ipix);
	      o_widx[0]=1.0;
	      o_hidx[0]=ipix;
	    }
	    else {
	      npt=Get_hidx(projFunc,ph[i],th[i],psi[i],ib,i,
			   (rg-globalBeginRing)/((double) (globalEndRing-globalBeginRing)),
			   External+i,
			   o_widx,
			   o_hidx,
			   rank);
	    }
	    
	    for (lll=0;lll<npt;lll++) {
	      hpix *tp_hpix = l_tp_hpix;
	      
	      memset(tp_hpix,0,sizeof(hpix));
	      
	      ipix=o_hidx[lll];
	      
	      if (ipix<0||ipix>=12*Nside*Nside) {
		fprintf(stderr,"Problem with pointing, healpix index %ld not in domain [%ld,%ld]\n",
			(long) ipix,
			(long)0, (long) 12*Nside*Nside);
		exit(0);
	      }
	      
	      for (int l_iter = 0; l_iter < number_of_iterations; l_iter++) {
		tp_hpix->listp[l_iter]=(y[l_iter][i]-Param->Monop[ib])/Param->Calibration[ib];
	      }
	      
	      tp_hpix->sig=tp_hpix->listp[0];
	      
	      if (Sub_HPR != NULL) {
		tp_hpix->Sub_HPR = Sub_HPR[i]*Param->SUB_HPRCOEF[ib];
	      } else {
		tp_hpix->Sub_HPR = 0.0;
	      }
	      
	      tp_hpix->corr_cnn = 0.0;

	      tp_hpix->hit =  h[i]*o_widx[lll];
	      

	      int is_valid = init_channels(tp_hpix,
					   projFunc,
					   psi[i],
					   External+i,
					   (rg-globalBeginRing)/((double) (globalEndRing-globalBeginRing)),
					   ipix,
					   ib,
					   i,
					   &(tp_hpix->sig),
					   &(tp_hpix->hpr_cal),
					   &(tp_hpix->hit),rank);

	      if (is_valid) {
		
		nrg_htmp++;
		if (n_l_hpix+1-MAXMPIBUFFER*rank_ptr_hpix>MAXMPIBUFFER) {
		  if (rank_ptr_hpix+1>alloc_ptr_hpix) {
		    fprintf(stderr,"Alloc to new buffer\n");
		    ptr_l_hpix = (hpix **) realloc(ptr_l_hpix,sizeof(hpix **)*(2*alloc_ptr_hpix));
		    alloc_ptr_hpix*=2;
		  }
		  if (memdisk) {
		    save_hpix_buffer(rank_ptr_hpix,ptr_l_hpix[rank_ptr_hpix],MAXMPIBUFFER,rank);
		  }
		  rank_ptr_hpix++;
		  // If not saved realloc memory
		  ptr_l_hpix[rank_ptr_hpix] = (hpix *) malloc(sizeof(hpix)*MAXMPIBUFFER);
		  com_hpix = ptr_l_hpix[rank_ptr_hpix];
		}

		hpix *l_hpix = com_hpix+n_l_hpix-MAXMPIBUFFER*rank_ptr_hpix;

		memcpy(l_hpix,tp_hpix,sizeof(hpix));

		tp_hpix=l_hpix;

		tp_hpix->inc=psi[i];

		for (j=0;j<NB_EXTERNAL;j++) tp_hpix->External[j]=External[i+j*RINGSIZE];
		
		tp_hpix->model=tp_hpix->hpr_cal;
		
		if (sparseFunc!=NULL) {
		  
		  tp_hpix->nShpr = calc_sparse_hpr(sparseFunc,
						   rg,
						   ib,
						   ipix,
						   i,
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
						    i,
						    psi[i],
						    tp_hpix->sig,
						    External+i);
		  
		}
		
		tp_hpix->w = tp_hpix->hit*sxi;
		tp_hpix->surv = surv;
		tp_hpix->ib = ib;
		tp_hpix->rg = rg;
		tp_hpix->gi = 0;
		tp_hpix->ipix = ipix;

		if (Param->regrid==1) {
		  stat_pix[ipix]+=1;
		}
		
		if (isnan(tp_hpix->sig)||isnan(tp_hpix->w)||isnan(tp_hpix->hpr_cal)) {
		  if (prtnan<10) {
		    fprintf(stderr,"NAN NAN NAN %lf %lf %lf %ld %ld\n",tp_hpix->hpr_cal,
			    tp_hpix->sig,tp_hpix->w,
			    (long)rg,(long)i);
		    prtnan++;
		  }
		}  
		else n_l_hpix+=1;
	      }
	    }
	  }
	}
	free(l_tp_hpix);
	  
	if (verbose==1) fprintf(stderr,"%s %d %d\n",__FILE__,__LINE__,rank);

	//if ((rank==rank_zero) && (rg==globalRankInfo.BeginRing[rank]) ) {
	{
	  GetProcMem(&vmem,&phymem);
	  fprintf(stderr,"Com %d %s Rank: %ld[%d] MEM %.1lf[%.1lf]MB Nd=%ld \n",
		  (int) ib, Command, (long) rank, getpid(),
		  (double) vmem/1024./1024.,
		  (double) phymem/1024./1024.,
		  (long) nrg_htmp);
	}

        free(h);
        if (Param->flag_stim_paramfiles == 0) {
          free(y[0]);
        }
	if (Sub_HPR!=NULL)
	  free(Sub_HPR);
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

    {
      GetProcMem(&vmem,&phymem);
      fprintf(stderr,"\nafter troll detector (%ld/%ld): used VMEM %.1lf[PHYS %.1lf]MB\n",
          ib, nbolo-1, (double) vmem/1024./1024., (double) phymem/1024./1024.);
    }

  } // end bolometer loop

  if (addpol!=NULL) {
    free(addpol);
  }

  if (memdisk) {
    save_hpix_buffer(rank_ptr_hpix,ptr_l_hpix[rank_ptr_hpix],n_l_hpix-rank_ptr_hpix*MAXMPIBUFFER,rank);
  }
  
  //PrintFreeMemOnNodes( rank, mpi_size, "before pixel balancing");

  if (sparseFunc!=NULL) {
    long lnpix=npixShpr;
    long mnpix;
    MPI_Allreduce(&lnpix,&mnpix,1,MPI_LONG, MPI_MAX,MPI_COMM_WORLD);
    npixShpr=mnpix;
    if (TestUpdateSparse==0&&TestNormalizeSparse==0)
      ClosepFunc(sparseFunc);
  }

  if (rank==rank_zero) fprintf(stderr,"Number Of Values %d\n",(int) n_l_hpix);
  if (rank==rank_zero) fprintf(stderr,"Number Of Sparse Value %d\n",(int) npixShpr);
  if (rank==rank_zero) fprintf(stderr,"Number of Detector %d\n",(int) nbolo);
  if (rank==rank_zero) fprintf(stderr,"Number of val means %d\n",(int) Param->n_val_mean);
#if 0
  if (Param->n_do_mean!=npixShpr*Param->n_val_mean) {
    fprintf(stderr, "do_mean param should have the size [%ld (Number Of Sparse Value x Number of Detector x Number of val weigts)] but found %ld\n",
	    (long) (npixShpr*Param->n_val_mean),
	    (long) (Param->n_do_mean));
    exit(0);
  }
#endif
 
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
    if (rank==rank_zero) fprintf(stderr,"Mask %ld\n",(long) resmask);
  }
  else {
    mask=NULL;
  }

  nnbpix=12*Nside*Nside/mpi_size;
 
  for (i=0;i<mpi_size;i++) {
    begpix[i]=i*nnbpix;
    edpix[i]=(i+1)*nnbpix-1;
  }
  edpix[mpi_size-1]=12*Nside*Nside-1;
  
  nnbpix=edpix[rank]-begpix[rank]+1;

  realpix = (long *) malloc(sizeof(long)*nnbpix);
  int *newmask = (int *) malloc(sizeof(int)*nnbpix);
  if (Param->regrid==1) {
      
    {
      int *l_stat_pix = (int *) malloc(12*Nside*Nside*sizeof(int));
      MPI_Reduce(stat_pix,l_stat_pix,12*Nside*Nside,MPI_INT,MPI_SUM,rank_zero,MPI_COMM_WORLD);
      
      memcpy(stat_pix,l_stat_pix,12*Nside*Nside*sizeof(int));
      free(l_stat_pix);
    }
    
    if (rank==rank_zero) {
      fprintf(stderr,"regrid\n");
      hpint *statp = (hpint *) malloc(12*Nside*Nside*sizeof(hpint));
      
      if (mask!=NULL) {
	for (j=0;j<12*Nside*Nside;j++) {
	  statp[j].ipx=j;
	  statp[j].hit=stat_pix[j]*mask[j];
	}
      }
      else {
	for (j=0;j<12*Nside*Nside;j++) {
	  statp[j].ipx=j;
	  statp[j].hit=stat_pix[j];
	}
      }
      
      fprintf(stderr,"sort\n");
      qsort(statp,12*Nside*Nside,sizeof(hpint),compar_int);
      
      for (j=0;j<12*Nside*Nside;j++) {
	statp[j].hit=j%(2*mpi_size);
      }
      for (j=0;j<12*Nside*Nside;j++) {
	if (statp[j].hit>mpi_size-1) statp[j].hit=mpi_size-1-statp[j].hit;
      }
      
      qsort(statp,12*Nside*Nside,sizeof(hpint),compar_int);
      
      for (rrk=0;rrk<mpi_size;rrk++) {
	for (i=0;i<edpix[rrk]-begpix[rrk]+1;i++) {
	  
	  stat_pix[i+begpix[rrk]] = statp[i+begpix[rrk]].ipx;
	  
	  stat_pix[12*Nside*Nside+statp[i+begpix[rrk]].ipx]=i+begpix[rrk];
	}
      }
      free(statp);
#if 1
      FILE *fp=fopen("statpix.dat","w");
      fwrite(stat_pix,sizeof(int)*12*Nside*Nside*2,1,fp);
      fclose(fp);
#endif
    }
      
    // NOT A BCAST TO BE CHANGED TO MPI_SCATTER FOR OPTIMALITY
    MPI_Bcast(stat_pix,sizeof(int)*12*Nside*Nside*2, MPI_BYTE, rank_zero, MPI_COMM_WORLD);
    
    for (i=0;i<nnbpix;i++) realpix[i]=stat_pix[i+begpix[rank]];
    if (mask!=NULL) {
      for (i=0;i<nnbpix;i++) newmask[i]=mask[realpix[i]];
    }
    else {
      for (i=0;i<nnbpix;i++) newmask[i]=1.0;
    }
  
    if (mask!=NULL) free(mask);
    mask=newmask;
    // AND NOW IN FRONT OF YOUR EYES MIX PIXELS TO GET THEM EFFICIENTLY COMPUTED!!!
    
    if (rank==rank_zero) fprintf(stderr,"rebin pixel\n");
    for (long lll=0;lll<rank_ptr_hpix+1;lll++) {
      long ed_buffer=MAXMPIBUFFER;
      
      if (lll==rank_ptr_hpix) ed_buffer=n_l_hpix-lll*MAXMPIBUFFER;
      if (lll>rank_ptr_hpix) ed_buffer=0;
      hpix *local_ptr=call_hpix_buffer(lll,rank);
      
      for (long k=0;k<ed_buffer;k++) {
	hpix *l_hpix = local_ptr+k;
	l_hpix->ipix = stat_pix[12*Nside*Nside+l_hpix->ipix-begpix[rank]];
      }
    }
    
    free(stat_pix);
  }
  else {
    for (i=0;i<nnbpix;i++) realpix[i]=rank+i*mpi_size;
    if (mask!=NULL) {
      for (i=0;i<nnbpix;i++) newmask[i]=mask[realpix[i]];
    
      free(mask);
    }
    else {
      for (i=0;i<nnbpix;i++) newmask[i]=1.0;
    }
    mask=newmask;

    if (rank==rank_zero) fprintf(stderr,"rebin pixel\n");
    for (long lll=0;lll<rank_ptr_hpix+1;lll++) {
      long ed_buffer=MAXMPIBUFFER;
      
      if (lll==rank_ptr_hpix) ed_buffer=n_l_hpix-lll*MAXMPIBUFFER;
      if (lll>rank_ptr_hpix) ed_buffer=0;
      if (ed_buffer>0) {
	hpix *local_ptr=call_hpix_buffer(lll,rank);
      
	for (long k=0;k<ed_buffer;k++) {
	  hpix *l_hpix = local_ptr+k;
	  l_hpix->ipix = begpix[l_hpix->ipix%mpi_size]+ (int) (l_hpix->ipix/mpi_size);
	}
      }
    }
    
  }
  
  GetProcMem(&vmem,&phymem);
  if (rank==rank_zero) fprintf(stderr,"Rank: %ld[%d] MEM %.1lf[%.1lf]MB line=%d\n",
			       (long) rank, getpid(),
			       (double) vmem/1024./1024.,
			       (double) phymem/1024./1024.,__LINE__);
  MPI_Barrier(MPI_COMM_WORLD);

  /*======================================================
    =
    =      order data per pixel and per proc
    =
    =*/
  
  
  nnbpix = edpix[rank]-begpix[rank]+1;
    

  long ntot_pts=0;
#if 0
  if (Param->regrid==1) { // do not allocate rank_map
    rank_map = (int *) malloc(sizeof(int)*12*Nside*Nside);
    memset(rank_map,0,sizeof(int)*12*Nside*Nside);
    for (long k=0;k<mpi_size;k++) {
      for (int l=begpix[k];l<edpix[k]+1;l++) rank_map[l]=realpix[k];
    }
  }
#endif
  
  MPI_Allreduce(&rank_ptr_hpix,&l_rank_ptr_hpix,1,MPI_LONG, MPI_MAX,MPI_COMM_WORLD);

  for (long lll=0;lll<l_rank_ptr_hpix+1;lll++) {
    long ed_buffer=MAXMPIBUFFER;
    hpix *ltbs=NULL;

    if (rank==rank_zero) fprintf(stderr,"Start Exchange[%ld/%ld]\n",(long) lll,
				 (long) l_rank_ptr_hpix+1);

    if (lll==rank_ptr_hpix) {
      ed_buffer=n_l_hpix-lll*MAXMPIBUFFER;
    }
    if (lll>rank_ptr_hpix) ed_buffer=0;
    else ltbs=call_hpix_buffer(lll,rank);

    MPI_Barrier(MPI_COMM_WORLD);
    if (rank==rank_zero) fprintf(stderr,"Exchange[%ld/%d]\n",(long) lll,__LINE__);

    long ntbs=ed_buffer;

    if (mpi_size>1) { 
#if 0
      if (Param->regrid==1) {
	for (long k=0;k<ntbs;k++) {
	  hpix *l_hpix = ltbs+k;
	  l_hpix->mpi_rank=rank_map[l_hpix->ipix];
	}
      }
      else {
#endif
	for (long k=0;k<ntbs;k++) {
	  hpix *l_hpix = ltbs+k;
	  int l_pix=0;
	  // could be optimize with dichotomie but should be taken into account odd/even mpi_size
	  while (edpix[l_pix]<l_hpix->ipix&&l_pix<mpi_size-1) l_pix++;
	  l_hpix->mpi_rank=l_pix;
	}
#if 0
      }
#endif
      if (ed_buffer>0)
	qsort(ltbs,ed_buffer,sizeof(hpix),compar_hpix_rank);


      MPI_Barrier(MPI_COMM_WORLD);
      if (rank==rank_zero) fprintf(stderr,"Exchange[%ld/%d]\n",(long) lll,__LINE__);

      hpix *otbs=NULL;

      long *begbuf = (long *) malloc((mpi_size+1)*sizeof(long));
      memset(begbuf,0,(mpi_size+1)*sizeof(long));
      
      for (long k=0;k<ntbs;k++) {
	hpix *l_hpix = ltbs+k;
	begbuf[l_hpix->mpi_rank+1]=(k+1)*sizeof(hpix);
      }

      // case some process have no data
      for (long k=1;k<mpi_size+1;k++) {
	if (begbuf[k]==0) {
	  begbuf[k]=begbuf[k-1];
	}
      }
      
      MPI_Barrier(MPI_COMM_WORLD);
      if (rank==rank_zero) fprintf(stderr,"Exchange[%ld/%d]\n",(long) lll,__LINE__);

      /* exchange all buffer information */
      long *allbuf = (long *) malloc(mpi_size*(mpi_size+1)*sizeof(long));
      MPI_Allgather(begbuf,mpi_size+1,MPI_LONG,allbuf,(mpi_size+1),MPI_LONG,MPI_COMM_WORLD);
      int *sdispls = (int *) malloc(sizeof(int)*mpi_size);
      int *sendcounts = (int *) malloc(sizeof(int)*mpi_size);
      int *recvcounts = (int *) malloc(sizeof(int)*mpi_size);
      int *rdispls = (int *) calloc(mpi_size,sizeof(int));
      
      MPI_Barrier(MPI_COMM_WORLD);
      if (rank==rank_zero) fprintf(stderr,"Exchange[%ld/%d]\n",(long) lll,__LINE__);

      for (i=0;i<mpi_size;i++) {
	sendcounts[i]=begbuf[i+1]-begbuf[i];
	sdispls[i]=begbuf[i];
	recvcounts[i]=allbuf[i*(mpi_size+1)+rank+1]-allbuf[i*(mpi_size+1)+rank];
      }
      free(begbuf);
      free(allbuf);
      rdispls[0]=0;
      for (i=1;i<mpi_size;i++) {
	rdispls[i]=rdispls[i-1]+recvcounts[i-1];
      }
      ntbs=(rdispls[mpi_size-1]+recvcounts[mpi_size-1])/sizeof(hpix);

      MPI_Barrier(MPI_COMM_WORLD);
      if (rank==rank_zero) fprintf(stderr,"Exchange[%ld/%d]\n",(long) lll,__LINE__);

      if (ntbs==0) {
	otbs = (hpix *) malloc(sizeof(hpix));
	memset(otbs,0,sizeof(hpix));
      }
      else{
	otbs = (hpix *) malloc(ntbs*sizeof(hpix));
	memset(otbs,0,ntbs*sizeof(hpix));
      }

      MPI_Alltoallv(ltbs,sendcounts,sdispls,MPI_BYTE,
		    otbs,recvcounts,rdispls,MPI_BYTE,
		    MPI_COMM_WORLD);
      
      free(sdispls);
      free(sendcounts);
      free(recvcounts);
      free(rdispls);
      ltbs=otbs;
      MPI_Barrier(MPI_COMM_WORLD);
      if (rank==rank_zero) fprintf(stderr,"Exchange[%ld/%d]\n",(long) lll,__LINE__);

      if (lll<=rank_ptr_hpix) free_hpix_buffer(lll,rank);
    }
    
    save_hpix_buffer(lll,ltbs,ntbs,rank);
    //ptr_l_hpix[lll]=ltbs;
    //loc_nhpix[lll]=ntbs;
    ntot_pts+=ntbs;
    
    if (rank==rank_zero) fprintf(stderr,"Exchange[%ld/%ld] %ld values\n",(long) lll,
				 (long) l_rank_ptr_hpix+1,(long) ntbs);
    
  }
#if 0
  if (Param->regrid==1) {
    free(rank_map);
  }
#endif

  l_rank_ptr_hpix=l_rank_ptr_hpix+1;

  if (rank==rank_zero)
    fprintf(stderr,"\n===============================================\n\n");
  MPI_Barrier(MPI_COMM_WORLD);
  
  fprintf(stderr,"   NTOT[%d] = %ld\n",(int) rank,(long) ntot_pts);
  
  MPI_Barrier(MPI_COMM_WORLD);
  if (rank==rank_zero) {
    fprintf(stderr,"\n\n===============================================\n");
    fprintf(stderr,"\n NBuffer = %ld\n",(long) l_rank_ptr_hpix);
    fprintf(stderr,"\n\n===============================================\n");
  }
  
  
  GetProcMem(&vmem,&phymem);
  if (rank%64==0) fprintf(stderr,"Rank: %ld[%d] MEM %.1lf[%.1lf]MB line=%d\n",
			  (long) rank, getpid(),
			  (double) vmem/1024./1024.,
			  (double) phymem/1024./1024.,__LINE__);
  
  PrintFreeMemOnNodes( rank, mpi_size, "before matrix computation");
  
  /*======================================================
    =
    =     Compute pixel index inside the process
    =
    =*/
  

  for(int ibuffer=0;ibuffer<l_rank_ptr_hpix;ibuffer++) {
    hpix *l_htmp=call_hpix_buffer(ibuffer,rank);
    for(int pix =0;pix<loc_nhpix[ibuffer];pix++){
      hpix *htmp = l_htmp+pix; 
      htmp->ipix-=begpix[rank];
    }
    save_hpix_buffer(ibuffer,l_htmp,loc_nhpix[ibuffer],rank);
  }

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


  flgpix = (PIOBYTE *) malloc(sizeof(PIOBYTE)*nnbpix);
  memset(flgpix,0,sizeof(PIOBYTE)*nnbpix);
  
  imatrice = (PIODOUBLE *) malloc(nnbpix*MAXCHANNELS*MAXCHANNELS*sizeof(PIODOUBLE));
  memset(imatrice,0,nnbpix*MAXCHANNELS*MAXCHANNELS*sizeof(PIODOUBLE));
  cond = (PIODOUBLE *) malloc(nnbpix*sizeof(PIODOUBLE));
  
  long l_nmatpix=0;
 
  fprintf(stderr,"Rank: %ld[%d] MEM %.1lf[%.1lf]MB line=%d\n",
    (long) rank, getpid(),
    (double) vmem/1024./1024.,
    (double) phymem/1024./1024.,__LINE__);
  
  for(int ibuffer=0;ibuffer<l_rank_ptr_hpix;ibuffer++) {
    hpix *l_htmp=call_hpix_buffer(ibuffer,rank);
    for(int pix =0;pix<loc_nhpix[ibuffer];pix++){
      hpix *htmp = l_htmp+pix;
      for(int i = 0;i<MAXCHANNELS;i++){
	double x = htmp->w*htmp->channels[i];
	for(int j = 0;j<MAXCHANNELS;j++){
	  imatrice[i+j*MAXCHANNELS+htmp->ipix*MAXCHANNELS*MAXCHANNELS] += x*htmp->channels[j];
	}
      }
    }
  }
  
  for (k=0;k<nnbpix;k++) {    
    cond[k] = 0.0;
    double matrice[MAXCHAN*MAXCHAN];
    memcpy(matrice,imatrice+k*MAXCHANNELS*MAXCHANNELS,MAXCHANNELS*MAXCHANNELS*sizeof(double)); 
    
    //test if matrice inversible 
    cond[k]=cond_thres(matrice,imatrice+MAXCHANNELS*MAXCHANNELS*k,MAXCHANNELS);
    if (cond[k] < Param->seuilcond) {
      flgpix[k]=1;
    }else {
      flgpix[k]=0;
      // Print detail for pix that does not meet cond requirement
      //fprintf(stderr,"[DBG COND] Pix#%ld is flagged out! II=%g IQ=%g IU=%g QQ=%g QU=%g UU=%g \n",k+begpix[rank], II, IQ, IU, QQ, QU, UU);
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  
  for(int ibuffer=0;ibuffer<l_rank_ptr_hpix;ibuffer++) {
    hpix *l_htmp=call_hpix_buffer(ibuffer,rank);
    for(int pix =0;pix<loc_nhpix[ibuffer];pix++){
      hpix *htmp = l_htmp+pix; 
      if (flgpix[htmp->ipix]==1) {
	flg_rg[htmp->ib][htmp->rg-globalBeginRing]=1;
      }
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
  
  if (rank==rank_zero) fprintf(stderr,"l_nmatpix[%d] %ld / %ld \n",rank,(long) l_nmatpix ,nnbpix );

  long lb;
  MPI_Reduce(&l_nmatpix,&lb,1,MPI_LONG,MPI_SUM,rank_zero,MPI_COMM_WORLD);
  nmatpix=lb;
  MPI_Bcast(&nmatpix, sizeof(long), MPI_BYTE, rank_zero, MPI_COMM_WORLD);

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
    for (ib=0;ib<nbolo+1;ib++) newnr[ib]=0;
  }
  
  if (rank==rank_zero) {
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

  for(int ibuffer=0;ibuffer<l_rank_ptr_hpix;ibuffer++) {
    hpix *l_htmp=call_hpix_buffer(ibuffer,rank);
    for(int pix =0;pix<loc_nhpix[ibuffer];pix++){
      hpix *htmp = l_htmp+pix; 
      if (flgpix[htmp->ipix]>0) {
	long ri1=htmp->rg-globalBeginRing;
	if (flg_rg[htmp->ib][ri1]!=0) {
	  
	  double divi=NEP_tab[htmp->ib]*htmp->hit;
	  
	  divi=divi*divi;
	  double ww=1/divi;
	  
	  double tmp=(htmp->model);
	  
	  histo_gi[ htmp->rg+htmp->ib*32000]+=ww*tmp;
	  histo2_gi[htmp->rg+htmp->ib*32000]+=ww*tmp*tmp;
	  histon_gi[htmp->rg+htmp->ib*32000]+=ww;
	  
	  mat_dip[0+4*htmp->ib]+=ww*htmp->model*htmp->model;
	  mat_dip[1+4*htmp->ib]+=ww*htmp->model;
	  mat_dip[2+4*htmp->ib]+=ww*htmp->model;
	  mat_dip[3+4*htmp->ib]+=ww;
	  
	  vec_dip[0+2*htmp->ib]+=ww*(htmp->sig)*htmp->model;
	  vec_dip[1+2*htmp->ib]+=ww*(htmp->sig);

	  if (isnan(ww)||isnan(htmp->sig)||isnan(htmp->model)) {
	    fprintf(stderr,"NAN AAA %lf %lf %lf\n",(double) ww,(double) htmp->sig,(double) htmp->model);
	  }
	}
      }
    }
  }


  //============================================================
  // remove dipole from the first harmonic
  //

  double l_mat[4*100];
  double l_vec[2*100];
  MPI_Reduce(mat_dip,l_mat,4*nbolo,MPI_DOUBLE,MPI_SUM,rank_zero,MPI_COMM_WORLD);
  MPI_Reduce(vec_dip,l_vec,2*nbolo,MPI_DOUBLE,MPI_SUM,rank_zero,MPI_COMM_WORLD);

  if (rank==rank_zero) {

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

  if (rank==rank_zero) fprintf(stderr,"BEG-ED %d %ld %ld - %ld %ld\n",rank,(long) nmatpix,(long) begpix[rank],(long) edpix[rank],
          (long) l_nmatpix);
  MPI_Barrier(MPI_COMM_WORLD);

  //double *tvec;
  //double *vvec;
  //double *totmat;

  // decoupe vecteur matrice par proc

  long maxsizemat=newnr[nbolo]+nbolo*(GAINSTEP)+npixShpr;
  
  if (rank==rank_zero) fprintf(stderr,"MAXSIZE %ld %ld\n",(long) maxsizemat,(long) newnr[nbolo]+nbolo*(GAINSTEP)+npixShpr);

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

  if (rank==rank_zero) fprintf(stderr,"BEG %d %ld %ld - %ld %ld\n",rank,(long) nmatpix,(long) begpix[rank],(long) edpix[rank],
		       (long) l_nmatpix);

  MPI_Barrier(MPI_COMM_WORLD);

  
  /*======================================================
    =
    =    LOOP ON SEEDs
    = 
    =*/
  
  for (int iter = 0; iter < number_of_iterations; iter++) {

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

    memset(x2,0,(newnr[nbolo]+nbolo*(GAINSTEP)+npixShpr)*sizeof(double));

    gainoff=0;

    nmatres=newnr[nbolo]+nbolo*(GAINSTEP)+npixShpr;

    double *x3= (double *) malloc(sizeof(double)*(nmatres)); 

    GetProcMem(&vmem,&phymem);
    if (rank==rank_zero) fprintf(stderr,"Rank: %ld[%d] MEM %.1lf[%.1lf]MB line=%d\n",
                (long) rank, getpid(),
                (double) vmem/1024./1024.,
                (double) phymem/1024./1024.,__LINE__);
    int itt;

    memset(x3,0,sizeof(double)*nmatres);

    itt=0;
    double resxi=1;

    while (itt<Param->NITT) {

      //=======================================================================================================
      // Start update if needed
      //=======================================================================================================
      if (TestUpdateSparse==1) {
	if (rank==rank_zero)
	  fprintf(stderr,"update_eval method available Update model to be fitted\n");

	// Init matrice and vecteur
	double *vector = malloc(MAXCHANNELS*sizeof(double)*nnbpix);
	memset(vector,0,MAXCHANNELS*sizeof(double)*nnbpix);

	//=======================================================================================================
	// Compute the map
	//=======================================================================================================
	for(int ibuffer=0;ibuffer<l_rank_ptr_hpix;ibuffer++) {
	  hpix *l_htmp=call_hpix_buffer(ibuffer,rank);
	  for(int pix =0;pix<loc_nhpix[ibuffer];pix++){
	    hpix *htmp = l_htmp+pix; 
	    if (flgpix[htmp->ipix]>0) {
	      long ri1=htmp->rg-globalBeginRing;
	      if (flg_rg[htmp->ib][ri1]!=0) {
		long iri1=rgord[htmp->ib][ri1]+newnr[htmp->ib];
		//calcul signal corriger
		double g1=gain[htmp->gi+htmp->ib*GAINSTEP];
	      
		double sig_corr = htmp->sig*g1 - htmp->Sub_HPR-htmp->corr_cnn;
	      
		if (do_offset==1) {
		  sig_corr-=x3[iri1];
		}
		
		if (REMOVE_CAL==1) {
		  sig_corr-=htmp->hpr_cal;
		}
		
		if(GAINSTEP!=0)  
		  sig_corr-=x3[newnr[nbolo]+htmp->gi+htmp->ib*GAINSTEP]*htmp->model;
		
		
		for(int m =0;m<htmp->nShpr;m++){
		  sig_corr-=x3[newnr[nbolo]+htmp->listofShpr_idx[m]+(GAINSTEP)*nbolo]*htmp->listofShpr[m];
		}
		
		//calcul matrix & vector
		for(int i = 0;i<MAXCHANNELS;i++){  
		  vector[i+htmp->ipix*MAXCHANNELS]+= htmp->w *htmp->channels[i]*sig_corr;
		}
	      } 
	    }
	  }
	}
	for (k=0;k<nnbpix;k++) {
	  if (flgpix[k]>0) {
	    apply_invertMatrix(imatrice+MAXCHANNELS*MAXCHANNELS*k,vector+MAXCHANNELS*k,MAXCHANNELS,rank);
	  }
	}
	
	for(int ibuffer=0;ibuffer<l_rank_ptr_hpix;ibuffer++) {
	  hpix *l_htmp=call_hpix_buffer(ibuffer,rank);
	  for(int pix =0;pix<loc_nhpix[ibuffer];pix++){
	    hpix *htmp = l_htmp+pix; 
	    if (flgpix[htmp->ipix]>0) {
	      long ri1=htmp->rg-globalBeginRing;
	      if (flg_rg[htmp->ib][ri1]!=0) {
		PyObject *pChan= PyList_New(MAXCHANNELS);
		PyObject *pVal = PyList_New(MAXCHANNELS);
		for (int i = 0; i < MAXCHANNELS; i++) {
		  PyList_SetItem(pVal, i, PyFloat_FromDouble((double) vector[i+htmp->ipix*MAXCHANNELS]));
		}
		PyObject *pIdx = PyList_New(htmp->nShpr);
		PyObject *pWw  = PyList_New(htmp->nShpr);
		for(int m =0;m<htmp->nShpr;m++){
		  PyList_SetItem(pIdx, m, PyLong_FromLong((long) htmp->listofShpr_idx[m]));
		  PyList_SetItem(pWw, m, PyFloat_FromDouble((double) htmp->listofShpr[m]));
		}
		for (int i = 0; i < MAXCHANNELS; i++) {
		  PyList_SetItem(pChan, i, PyFloat_FromDouble((double) htmp->channels[i]));
		}
		
		PyObject *pValue = PyObject_CallMethod(sparseFunc, "update_eval", "(OOOO)", pIdx,pWw,pChan,pVal);
		
		// Traitement de la valeur de retour
		if (pValue != NULL) {
		  // Assurer que pValue est un tuple
		  if (PyTuple_Check(pValue) && PyTuple_Size(pValue) == 2) {
		    // Extraire les deux tableaux du tuple
		    int n1=copy_int_array(PyTuple_GetItem(pValue, 0), htmp->listofShpr_idx,MAXEXTERNALSHPR);
		    int n2=copy_float_array(PyTuple_GetItem(pValue, 1), htmp->listofShpr,MAXEXTERNALSHPR);
		    if (n1!=n2) {
		      fprintf(stderr, "Update Sparse function should provide an equal number of invex and value, here Sroll received %d %d\n",n1,n2);
		      exit(0);
		    }
		    htmp->nShpr=n1;
		    Py_DECREF(pValue);
		  }
		}
		Py_DECREF(pIdx);
		Py_DECREF(pWw);
		
		// Nettoyage des arguments
		Py_DECREF(pChan);
		Py_DECREF(pVal);
	      }
	    }
	  }
	}
	free(vector);
      }
      //=======================================================================================================
      // Start TOD correction if needed
      //=======================================================================================================
      if (TestCorrTOD==1) {
	if (rank==rank_zero)
	  fprintf(stderr,"corr_tod_eval method available clean tod\n");
	//extract corrected data for timeline interpretation
	double **l_signal=(double **) malloc(sizeof(double*)*nnbpix);
	double **l_weights=(double **) malloc(sizeof(double*)*nnbpix);
	double **l_inc=(double **) malloc(sizeof(double*)*nnbpix);
	int **l_time=(int **) malloc(sizeof(int*)*nnbpix);
	int **l_did=(int **) malloc(sizeof(int*)*nnbpix);
	double **l_external=(double **) malloc(sizeof(double*)*nnbpix);
	hpix ***l_pointer=(hpix ***) malloc(sizeof(hpix **)*nnbpix);
	long *l_ndata=(long *) malloc(sizeof(long)*nnbpix);

	memset(l_ndata,0,sizeof(long)*nnbpix);

	for(int ibuffer=0;ibuffer<l_rank_ptr_hpix;ibuffer++) {
	  hpix *l_htmp=call_hpix_buffer(ibuffer,rank);
	  for(int pix =0;pix<loc_nhpix[ibuffer];pix++){
	    hpix *htmp = l_htmp+pix; 
	    if (flgpix[htmp->ipix]>0) {
	      long ri1=htmp->rg-globalBeginRing;
	      if (flg_rg[htmp->ib][ri1]!=0) {
		long iri1=rgord[htmp->ib][ri1]+newnr[htmp->ib];
		//calcul signal corriger
		double g1=gain[htmp->gi+htmp->ib*GAINSTEP];
	      
		double sig_corr = htmp->sig*g1 - htmp->Sub_HPR-htmp->corr_cnn;
	      
		if (do_offset==1) {
		  sig_corr-=x3[iri1];
		}
		
		if (REMOVE_CAL==1) {
		  sig_corr-=htmp->hpr_cal;
		}
		
		if(GAINSTEP!=0)  
		  sig_corr-=x3[newnr[nbolo]+htmp->gi+htmp->ib*GAINSTEP]*htmp->model;
		
		
		for(int m =0;m<htmp->nShpr;m++){
		  sig_corr-=x3[newnr[nbolo]+htmp->listofShpr_idx[m]+(GAINSTEP)*nbolo]*htmp->listofShpr[m];
		}
		
		if (l_ndata[htmp->ipix]==0) {
		  l_signal[htmp->ipix]  = (double *) malloc(sizeof(double));
		  l_weights[htmp->ipix] = (double *) malloc(sizeof(double));
		  l_inc[htmp->ipix]     = (double *) malloc(sizeof(double));
		  l_time[htmp->ipix]    = (int *) malloc(sizeof(int));
		  l_did[htmp->ipix]     = (int *) malloc(sizeof(int));
		  l_pointer[htmp->ipix]     = (hpix **) malloc(sizeof(hpix *));
		  l_external[htmp->ipix]    = (double *) malloc(sizeof(double)*NB_EXTERNAL);
		  
		}
		else {
		  l_signal[htmp->ipix]  = (double *) realloc(l_signal[htmp->ipix] ,sizeof(double)*(l_ndata[htmp->ipix]+1));
		  l_weights[htmp->ipix] = (double *) realloc(l_weights[htmp->ipix],sizeof(double)*(l_ndata[htmp->ipix]+1));
		  l_inc[htmp->ipix]     = (double *) realloc(l_inc[htmp->ipix],sizeof(double)*(l_ndata[htmp->ipix]+1));
		  l_time[htmp->ipix]    = (int *) realloc(l_time[htmp->ipix]      ,sizeof(int)*(l_ndata[htmp->ipix]+1));
		  l_did[htmp->ipix]     = (int *)   realloc(l_did[htmp->ipix]     ,sizeof(int)*(l_ndata[htmp->ipix]+1));
		  l_pointer[htmp->ipix] = (hpix **) realloc(l_pointer[htmp->ipix] ,sizeof(hpix *)*(l_ndata[htmp->ipix]+1));
		  l_external[htmp->ipix]    = (double *) realloc(l_external[htmp->ipix],
								    sizeof(double)*(l_ndata[htmp->ipix]+1)*NB_EXTERNAL);
		}
		l_signal[htmp->ipix][l_ndata[htmp->ipix]]=sig_corr;
		l_weights[htmp->ipix][l_ndata[htmp->ipix]]=htmp->w;
		l_inc[htmp->ipix][l_ndata[htmp->ipix]]=htmp->inc;
		l_time[htmp->ipix][l_ndata[htmp->ipix]]=ri1;
		l_did[htmp->ipix][l_ndata[htmp->ipix]]=htmp->ib;
		l_pointer[htmp->ipix][l_ndata[htmp->ipix]]=htmp;
		for (int o=0;o<NB_EXTERNAL;o++) {
		  l_external[htmp->ipix][l_ndata[htmp->ipix]*NB_EXTERNAL+o]=htmp->External[o];
		}
		l_ndata[htmp->ipix]+=1;
	      } 
	    }
	  }
	}
	int already_computed=(int)(nnbpix/100);
	for (k=0;k<nnbpix;k++) 
	  {
	    if (rank==0&&k==already_computed) {
	      fprintf(stderr,"Compute TOD correction %.2f%%\n",((double)(100*k))/nnbpix);
	      already_computed+=(int)(nnbpix/100);
	    }
	    if (l_ndata[k]>0) {
	      double *o_signal = (double *) malloc(l_ndata[k]*sizeof(double));
	      double *o_hit = (double *) malloc(l_ndata[k]*sizeof(double));
	      
	      int ncorrtod= calc_corrtod_hpr(CorrTODFunc,
					     l_signal[k],
					     l_weights[k],
					     l_inc[k],
					     l_time[k],
					     l_did[k],
					     realpix[k],
					     l_external[k],
					     o_signal,
					     o_hit,
					     l_ndata[k]);
	      for (int o=0;o<ncorrtod;o++) {
		hpix *htmp=l_pointer[k][o];
		htmp->corr_cnn=o_signal[o];
		htmp->w=o_hit[o];
	      }
	      
	      free(o_signal);
	      
	    }
	  }


#if 0
	char FILENAME[1024];
	sprintf(FILENAME,"/home1/scratch/jmdeloui/SIGNAL_%d.dat",(int) rank);
	FILE *fp=fopen(FILENAME,"wb");
	for (k=0;k<nnbpix;k++) if (l_ndata[k]>0) {
	    fwrite(l_signal[k],l_ndata[k]*sizeof(double),1,fp);
	  }
	fclose(fp);
	sprintf(FILENAME,"/home1/scratch/jmdeloui/DID_%d.dat",(int) rank);
	fp=fopen(FILENAME,"wb");
	for (k=0;k<nnbpix;k++) if (l_ndata[k]>0) {
	    fwrite(l_did[k],l_ndata[k]*sizeof(double),1,fp);
	  }
	fclose(fp);
	sprintf(FILENAME,"/home1/scratch/jmdeloui/WEIGHTS_%d.dat",(int) rank);
	fp=fopen(FILENAME,"wb");
	for (k=0;k<nnbpix;k++) if (l_ndata[k]>0) {
	    fwrite(l_weights[k],l_ndata[k]*sizeof(double),1,fp);
	  }
	fclose(fp);
	sprintf(FILENAME,"/home1/scratch/jmdeloui/TIME_%d.dat",(int) rank);
	fp=fopen(FILENAME,"wb");
	for (k=0;k<nnbpix;k++) if (l_ndata[k]>0) {
	    fwrite(l_time[k],l_ndata[k]*sizeof(double),1,fp);
	  }
	fclose(fp);
	sprintf(FILENAME,"/home1/scratch/jmdeloui/INC%d_%d.dat",(int) i,(int) rank);
	fp=fopen(FILENAME,"wb");
	for (k=0;k<nnbpix;k++) if (l_ndata[k]>0) {
	    fwrite(l_inc[k],l_ndata[k]*sizeof(double),1,fp);
	  }
	fclose(fp);
	sprintf(FILENAME,"/home1/scratch/jmdeloui/IDX_%d.dat",(int) rank);
	fp=fopen(FILENAME,"wb");
	for (k=0;k<nnbpix;k++) if (l_ndata[k]>0) {
	    int *pixindex = (int *) malloc(l_ndata[k]*sizeof(int));
	    for (int l=0;l<l_ndata[k];l++) pixindex[l]=realpix[k];
	    fwrite(pixindex,l_ndata[k]*sizeof(int),1,fp);
	    free(pixindex);
	  }
	fclose(fp);
#endif	    

	free(l_signal);
	free(l_weights);
	free(l_time);
	free(l_did);
	free(l_pointer);
	free(l_external);
	free(l_inc);
      } // END testCorrTod

      memset(newnr[nbolo]+x3,0,sizeof(double)*nbolo*GAINSTEP);
     
      if (nmatres>0)
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

      if (rank==rank_zero) {
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
	for (i=(GAINSTEP)*nbolo+newnr[nbolo];i<(GAINSTEP)*nbolo+newnr[nbolo];i++) {
	  if ((i-newnr[nbolo])%nbolo==0) fprintf(stderr,"TMAP=[");
	  fprintf(stderr,"%lg,",x3[i]);
	  if ((i-newnr[nbolo])%nbolo==nbolo-1) fprintf(stderr,"]\n");
	}
	int nbSTF=npixShpr;
	if (nbSTF>200) nbSTF=200;
	for (i=(GAINSTEP)*nbolo+newnr[nbolo];i<(GAINSTEP)*nbolo+newnr[nbolo]+nbSTF;i++) {
	  if ((i-newnr[nbolo])%nbolo==0) fprintf(stderr,"STF=[");
	  fprintf(stderr,"%lg,",x3[i]);
	  if ((i-newnr[nbolo])%nbolo==nbolo-1) fprintf(stderr,"]\n");
	}	
      }

      itt++;
      MPI_Barrier(MPI_COMM_WORLD);
    }


    MPI_Barrier(MPI_COMM_WORLD);

    nmatres=newnr[nbolo]+nbolo*(GAINSTEP)+npixShpr;
    if (rank==rank_zero) {

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
      fprintf(stderr,"Write MAT  %lld\n",(long long) PIOWriteVECT(saveg,x3+newnr[nbolo],0,(nbolo*(GAINSTEP)+npixShpr)*sizeof(PIODOUBLE)));

    }

    double avvgain=0;
    if (GAINSTEP>0) {
      for (i=0;i<nbolo*GAINSTEP;i++) avvgain+=x3[newnr[nbolo]+i];
      avvgain/=((double)(nbolo*GAINSTEP));

      if (rank==rank_zero)  {
	fprintf(stderr,"AVVGAIN: %lf\n",avvgain);
      }
    }

    if (iter==number_of_iterations-1) {
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
  
  double *diag_avv=NULL;
  double *diag_avv2=NULL;
  double *diag_n=NULL;
  
  if (diagFunc!=NULL) {
    diag_avv  = (double *) malloc(sizeof(double)*nb_diag);
    diag_avv2  = (double *) malloc(sizeof(double)*nb_diag);
    diag_n = (double *) malloc(sizeof(double)*nb_diag);
  }
  
  // Init matrice and vecteur
  double *vector = malloc(MAXCHANNELS*sizeof(double)*nnbpix);
  double *rvector = malloc(MAXCHANNELS*sizeof(double)*nnbpix);
  double *avv = malloc(MAXCHANNELS*sizeof(double)*nnbpix);
  double *navv = malloc(MAXCHANNELS*sizeof(double)*nnbpix);
  double *avv2 = malloc(MAXCHANNELS*sizeof(double)*nnbpix);
  
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
      
      for (k=0;k<nnbpix;k++) cond[k]=UNSEENPIX;  // ie hp.UNSEEN
	 
      memset(vector,0,MAXCHANNELS*sizeof(double)*nnbpix);
      memset(rvector,0,MAXCHANNELS*sizeof(double)*nnbpix);
      memset(avv,0,MAXCHANNELS*sizeof(double)*nnbpix);
      memset(navv,0,MAXCHANNELS*sizeof(double)*nnbpix);
      memset(avv2,0,MAXCHANNELS*sizeof(double)*nnbpix);
      memset(imatrice,0,MAXCHANNELS*MAXCHANNELS*sizeof(double)*nnbpix);
      
      for(int ibuffer=0;ibuffer<l_rank_ptr_hpix;ibuffer++) {
	hpix *l_htmp=call_hpix_buffer(ibuffer,rank);
	for(int pix =0;pix<loc_nhpix[ibuffer];pix++){
	  hpix *htmp = l_htmp+pix; 
	  if (flgpix[htmp->ipix]>0) {
	    long ri1=htmp->rg-globalBeginRing;
	    if (flg_rg[htmp->ib][ri1]!=0&&Param->bolomask[detset*nbolo+htmp->ib]==1 && 
		htmp->rg>Param->beg_surv[isurv]&&htmp->rg<=Param->end_surv[isurv]) {

	      long iri1=rgord[htmp->ib][ri1]+newnr[htmp->ib];
	      //calcul signal corriger
	      double g1=gain[htmp->gi+htmp->ib*GAINSTEP];
	      double rsig = htmp->sig;
	      
	      double sig_corr = htmp->sig*g1-htmp->Sub_HPR;
	      double sig_corr2 = sig_corr;

	      sig_corr -= htmp->corr_cnn;
	      
	      if (do_offset==1) {
		sig_corr-=x3[iri1];
	      }
	      
	      if (REMOVE_CAL==1) {
		sig_corr-=htmp->hpr_cal;
		sig_corr2-=htmp->hpr_cal;
		rsig-=htmp->hpr_cal;
	      }
	      
	      
	      if(GAINSTEP2!=0)  
		sig_corr-= x3[newnr[nbolo]+htmp->gi+htmp->ib*GAINSTEP2]*htmp->model;
	      
	      
	      for(int m =0;m<htmp->nShpr;m++){
		sig_corr-=x3[newnr[nbolo]+htmp->listofShpr_idx[m]+(GAINSTEP2)*nbolo]*htmp->listofShpr[m];
	      }
	      
	      //calcul matrix & vector
	      for(int i = 0;i<MAXCHANNELS;i++){  
		vector[htmp->ipix+i*nnbpix]+= htmp->w *htmp->channels[i]*sig_corr;
	      }
	      for(int i = 0;i<MAXCHANNELS;i++){  
		rvector[htmp->ipix+i*nnbpix]+= htmp->w *htmp->channels[i]*rsig;
	      }
	      for(int i = 0;i<MAXCHANNELS;i++){
		for(int j = 0;j<MAXCHANNELS;j++){
		  imatrice[i+j*MAXCHANNELS+htmp->ipix*MAXCHANNELS*MAXCHANNELS] += htmp->w *htmp->channels[j]*htmp->channels[i];
		}
	      }
	    }
	  }
	}
      }

      for (k=0;k<nnbpix;k++) {   
	cond[k] = 0.0;
	double matrice[MAXCHAN*MAXCHAN];
	memcpy(matrice,imatrice+k*MAXCHANNELS*MAXCHANNELS,MAXCHANNELS*MAXCHANNELS*sizeof(double)); 
    
	//test if matrice inversible 
	cond[k]=cond_thres(matrice,imatrice+MAXCHANNELS*MAXCHANNELS*k,MAXCHANNELS);
    
	if (cond[k] < Param->seuilcond) {
	  double tmp[MAXCHAN];
	  double rtmp[MAXCHAN];
	  double *mat=imatrice+MAXCHANNELS*MAXCHANNELS*k;
	  for (i=0;i<MAXCHANNELS;i++) {
	    tmp[i] = 0.0;
	    rtmp[i] = 0.0;
	    for (j=0;j<MAXCHANNELS;j++) {
	      tmp[i]+=mat[j+MAXCHANNELS*i]*vector[j*nnbpix+k];
	      rtmp[i]+=mat[j+MAXCHANNELS*i]*rvector[j*nnbpix+k];
	    }
	  }
	  for (i=0;i<MAXCHANNELS;i++) {
	    vector[i*nnbpix+k]=tmp[i];
	    rvector[i*nnbpix+k]=rtmp[i];
	  }
	}
	else {
	  cond[k]=UNSEENPIX;
	  for (i=0;i<MAXCHANNELS;i++) {
	    vector[i*nnbpix+k]=UNSEENPIX;
	    rvector[i*nnbpix+k]=UNSEENPIX;
	  }
	}
      }

      /*================================================================================================
	COMPUTE CHI2 MAP 
	================================================================================================*/
      for(int ibuffer=0;ibuffer<l_rank_ptr_hpix;ibuffer++) {
	hpix *l_htmp=call_hpix_buffer(ibuffer,rank);
	for(int pix =0;pix<loc_nhpix[ibuffer];pix++){
	  hpix *htmp = l_htmp+pix; 
	  if (cond[htmp->ipix]!=UNSEENPIX) {
	    long ri1=htmp->rg-globalBeginRing;
	    if (flg_rg[htmp->ib][ri1]!=0&&Param->bolomask[detset*nbolo+htmp->ib]==1 && 
		htmp->rg>Param->beg_surv[isurv]&&htmp->rg<=Param->end_surv[isurv]) {

	      long iri1=rgord[htmp->ib][ri1]+newnr[htmp->ib];
	      //calcul signal corriger
	      double g1=gain[htmp->gi+htmp->ib*GAINSTEP];
	      
	      double sig_corr = htmp->sig*g1-htmp->Sub_HPR;
	      double sig_corr2 = sig_corr;

	      sig_corr -= htmp->corr_cnn;
	      
	      if (do_offset==1) {
		sig_corr-=x3[iri1];
	      }
	      
	      if (REMOVE_CAL==1) {
		sig_corr-=htmp->hpr_cal;
		sig_corr2-=htmp->hpr_cal;
	      }
	      
	      
	      if(GAINSTEP2!=0)  
		sig_corr-= x3[newnr[nbolo]+htmp->gi+htmp->ib*GAINSTEP2]*htmp->model;
	      
	      
	      for(int m =0;m<htmp->nShpr;m++){
		sig_corr-=x3[newnr[nbolo]+htmp->listofShpr_idx[m]+(GAINSTEP2)*nbolo]*htmp->listofShpr[m];
	      }
	      
	      //calcul matrix & vector
	      for(int i = 0;i<MAXCHANNELS;i++){  
		sig_corr-=htmp->channels[i]*vector[htmp->ipix+i*nnbpix];
	      }

	      avv[htmp->ipix]=avv[htmp->ipix]+htmp->w *sig_corr;
	      avv2[htmp->ipix]=avv2[htmp->ipix]+htmp->w *sig_corr*sig_corr;
	      navv[htmp->ipix]=navv[htmp->ipix]+htmp->w;
	      if (diagFunc!=NULL) {
		diag_avv[htmp->diag_idx]  += htmp->w *sig_corr;
		diag_avv2[htmp->diag_idx] += htmp->w *sig_corr*sig_corr;
		diag_n[htmp->diag_idx]    += htmp->w;
	      }
	    }
	  }
	}
      }
    

      for (k=0;k<nnbpix;k++) {  
	if (cond[k]!=UNSEENPIX) {
	  avv[k]=sqrt(avv2[k]/navv[k]-(avv[k]/navv[k])*(avv[k]/navv[k]));
	}
	else avv[k]=UNSEENPIX;
      }
      
      if (diagFunc!=NULL) {
	
	double *l_avv = (double *) malloc(sizeof(double)*(nb_diag));
	double *l_avv2 = (double *) malloc(sizeof(double)*(nb_diag));
	double *l_n = (double *) malloc(sizeof(double)*(nb_diag));
	
	MPI_Reduce(diag_avv,l_avv,nb_diag,MPI_DOUBLE,MPI_SUM,rank_zero,MPI_COMM_WORLD);
	MPI_Reduce(diag_avv2,l_avv2,nb_diag,MPI_DOUBLE,MPI_SUM,rank_zero,MPI_COMM_WORLD);
	MPI_Reduce(diag_n,l_n,nb_diag,MPI_DOUBLE,MPI_SUM,rank_zero,MPI_COMM_WORLD);
	
	for (k=0;k<nb_diag;k++) diag_avv[k]=sqrt(l_avv2[k]/l_n[k]-(l_avv[k]/l_n[k])*(l_avv[k]/l_n[k]));
	
	if  (rank==rank_zero) {
	  char TEST_OUTMAP[MAX_OUT_NAME_LENGTH];
	  sprintf(TEST_OUTMAP,"%s_%s_DIAG", mapout[detset],mapname);
	  PIOWriteVECT(TEST_OUTMAP,diag_avv,0,sizeof(PIODOUBLE)*nb_diag);
	  sprintf(TEST_OUTMAP,"%s_%s_DIAG_W", mapout[detset],mapname);
	  PIOWriteVECT(TEST_OUTMAP,l_n,0,sizeof(PIODOUBLE)*nb_diag);
	  sprintf(TEST_OUTMAP,"%s_%s_DIAG_L1", mapout[detset],mapname);
	  PIOWriteVECT(TEST_OUTMAP,l_avv,0,sizeof(PIODOUBLE)*nb_diag);
	  sprintf(TEST_OUTMAP,"%s_%s_DIAG_L2", mapout[detset],mapname);
	  PIOWriteVECT(TEST_OUTMAP,l_avv,0,sizeof(PIODOUBLE)*nb_diag);
	}
	
	free(l_avv);
	free(l_avv2);
	free(l_n);
      }
      
	      
      if (TestCorrTODEval_corr) {
	double *l_vector = (double *) malloc(MAXCHANNELS*sizeof(double)*nnbpix);
	memset(l_vector,0,MAXCHANNELS*sizeof(double)*nnbpix);
	for (int pix=0;pix<nnbpix;pix++) {
	  l_vector[pix]=npar_corrtod(CorrTODFunc,realpix[pix]);
	}
	if (verbose==1) fprintf(stderr,"%s %d %d\n",__FILE__,__LINE__,rank);
	{
	  char TEST_OUTMAP[MAX_OUT_NAME_LENGTH];
	  sprintf(TEST_OUTMAP,"%s_%s_PAR", mapout[detset],mapname);
	  PIOWriteMAP(TEST_OUTMAP,l_vector,begpix[rank],begpix[rank]+nnbpix-1);
	  MPI_Barrier(MPI_COMM_WORLD);
	}
	
	double *delta = (double *) malloc(MAXCHANNELS*sizeof(double));
	for (int istep=0;istep<Param->n_inc_sub_ref;istep++) {
	  memset(l_vector,0,MAXCHANNELS*sizeof(double)*nnbpix);
	  for (int pix=0;pix<nnbpix;pix++) {
	    int err=npar_corrtod(CorrTODFunc,realpix[pix]);
	    if (err>0) {
	      eval_corrtod(CorrTODFunc,
			   Param->inc_sub_ref[istep],
			   Param->rg_sub_ref[istep],
			   realpix[pix],
			   delta);
	      for(int i = 0;i<MAXCHANNELS;i++){
		l_vector[pix+i*nnbpix]=vector[pix+i*nnbpix]-delta[i];
	      }
	    }
	    else {
	      for(int i = 0;i<MAXCHANNELS;i++){
		l_vector[pix+i*nnbpix]=UNSEENPIX;
	      }
	    }
	  }
	  if (verbose==1) fprintf(stderr,"%s %d %d\n",__FILE__,__LINE__,rank);
	  for(int i = 0;i<MAXCHANNELS;i++){
	    char TEST_OUTMAP[MAX_OUT_NAME_LENGTH];
	    sprintf(TEST_OUTMAP,"%s_%s_%d_%s", mapout[detset],mapname,i,Param->name_sub[istep]);
	    PIOWriteMAP(TEST_OUTMAP,l_vector+i*nnbpix,begpix[rank],begpix[rank]+nnbpix-1);
	    MPI_Barrier(MPI_COMM_WORLD);
	  }
	}
	free(l_vector);
	free(delta);
      }
      
      if (verbose==1) fprintf(stderr,"%s %d %d\n",__FILE__,__LINE__,rank);
      for(int i = 0;i<MAXCHANNELS;i++){
	char TEST_OUTMAP[MAX_OUT_NAME_LENGTH];
	sprintf(TEST_OUTMAP,"%s_%s_%d", mapout[detset],mapname,i);
	PIOWriteMAP(TEST_OUTMAP,vector+i*nnbpix,begpix[rank],begpix[rank]+nnbpix-1);
	MPI_Barrier(MPI_COMM_WORLD);
      }
	
      if (verbose==1) fprintf(stderr,"%s %d %d\n",__FILE__,__LINE__,rank);
      for(int i = 0;i<MAXCHANNELS;i++){
	char TEST_OUTMAP[MAX_OUT_NAME_LENGTH];
	sprintf(TEST_OUTMAP,"%s_%s_%d_RAW", mapout[detset],mapname,i);
	PIOWriteMAP(TEST_OUTMAP,rvector+i*nnbpix,begpix[rank],begpix[rank]+nnbpix-1);
	MPI_Barrier(MPI_COMM_WORLD);
      }
      
      if (verbose==1) fprintf(stderr,"%s %d %d\n",__FILE__,__LINE__,rank);
      {
	char TEST_OUTMAP[MAX_OUT_NAME_LENGTH];
	sprintf(TEST_OUTMAP,"%s_%s_STD", mapout[detset],mapname);
	PIOWriteMAP(TEST_OUTMAP,avv,begpix[rank],begpix[rank]+nnbpix-1);
	MPI_Barrier(MPI_COMM_WORLD);
      }
      
      // Save the cond matrix
      if (verbose==1) fprintf(stderr,"%s %d %d\n",__FILE__,__LINE__,rank);
      {
	char TEST_OUTMAP[MAX_OUT_NAME_LENGTH];
	sprintf(TEST_OUTMAP,"%s_%s_COND", mapout[detset],mapname);
	PIOWriteMAP(TEST_OUTMAP,cond,begpix[rank],begpix[rank]+nnbpix-1);
	MPI_Barrier(MPI_COMM_WORLD);
      }
    }
  }
  free(avv);
  free(avv2);
  free(navv);
  free(vector);
  free(rvector);
  if (diagFunc!=NULL) {
    free(diag_avv);
    free(diag_avv2);
    free(diag_n);
  }

  
  
  }

  if (rank==rank_zero) {
    now = time( NULL);
    fprintf(stderr, "\n%s: --------------------------\n", __FILE__ );
    fprintf(stderr, "%s: Finished successfully at %s",   __FILE__, ctime( &now));
    fprintf(stderr, "%s: --------------------------\n", __FILE__ );
  }
  MPI_Finalize();        /* free parameters info */

  exit (0);
 }
	      
    

	      
