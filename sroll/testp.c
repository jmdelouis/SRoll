#define PY_SSIZE_T_CLEAN
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <stdio.h>
#include <numpy/arrayobject.h>
#include <mpi.h>
#include <stdio.h>

#define MAXRAND 1000000

PyObject *EXECPYTHON(PyObject *TheObject)
{
  if (TheObject==NULL) {
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

int main(int argc, char** argv) {
  int itt,k;
  // Initialize the MPI environment
  MPI_Init(&argc, &argv);
  
  // Get the number of processes
  int mpi_size;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  
  // Get the rank of the process
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  char *path=argv[1];
  int nbolo=atoi(argv[2]);
  int xcnn=atoi(argv[3]);
  int ycnn=atoi(argv[4]);
  
  char saveval[1024];
  struct stat buf;
  sprintf(saveval,"%s_CNNINFO_RPIX_%d",path,rank);
  stat(saveval,&buf);
  long nnbpix=buf.st_size/4;
  
  sprintf(saveval,"%s_CNNINFO_SIG_%d",path,rank);
  stat(saveval,&buf);
  long all_ndata=buf.st_size/4;
  
  Py_Initialize();
  
  PyObject * sys = EXECPYTHON(PyImport_ImportModule("sys"));
  PyObject * py_path = EXECPYTHON(PyObject_GetAttrString(sys, "path"));
  PyList_Append(py_path, PyString_FromString("/scratch/cnt0028/ias1717/SHARED/bware/sroll22_trick/SrollEx/sroll"));
  
  PyObject * ModuleString = PyString_FromString((char*) "py_function");
  PyObject * Module = EXECPYTHON(PyImport_Import(ModuleString));
  PyObject * Dict = EXECPYTHON(PyModule_GetDict(Module));
  
  PyObject * init     = CALLPYTHON(Dict, "init_shape");
  PyObject * init_net = CALLPYTHON(Dict, "init_network");
  PyObject * allocf32 = CALLPYTHON(Dict, "alloc_table_float32");
  PyObject * alloci32 = CALLPYTHON(Dict, "alloc_table_int32");
  PyObject * Clean    = CALLPYTHON(Dict, "free_table");
  PyObject * grad     = CALLPYTHON(Dict, "calc_grad");
  PyObject * agrad    = CALLPYTHON(Dict, "apply_grad");
  PyObject * gloss    = CALLPYTHON(Dict, "get_loss");
  PyObject * pred     = CALLPYTHON(Dict, "get_prediction");
  
  PyObject *arglist = Py_BuildValue("(llll)", (long) 8,(long) 256,(long) 256,(long) rank);
  PyObject *mynetwork;
  mynetwork=EXECPYTHON(PyObject_CallObject(init, arglist));
  Py_DECREF(arglist); 
  
  arglist = Py_BuildValue("(l)", (long) nnbpix);
  PyObject *realpix = EXECPYTHON(PyObject_CallObject(alloci32, arglist));
  Py_DECREF(arglist); 

  float *l_random = (float *) malloc(sizeof(float)*MAXRAND);
  for (k=0;k<MAXRAND;k++) {
    l_random[k]=drand48()-0.5;
  }
  MPI_Bcast(l_random, sizeof(float)*MAXRAND, MPI_BYTE, 0, MPI_COMM_WORLD);

  arglist = Py_BuildValue("(l)", (long) MAXRAND);
  PyObject *randombase = EXECPYTHON(PyObject_CallObject(allocf32, arglist));
  Py_DECREF(arglist); 
  float *ptr_random = PyArray_DATA((PyArrayObject *) randombase);
  for (k=0;k<MAXRAND;k++) {
    ptr_random[k]=l_random[k];
  }
  free(l_random);
  
  arglist = Py_BuildValue("(l)", (long) all_ndata);
  PyObject *signal  = EXECPYTHON(PyObject_CallObject(allocf32, arglist));
  PyObject *weights = EXECPYTHON(PyObject_CallObject(allocf32, arglist));
  PyObject *TCO1    = EXECPYTHON(PyObject_CallObject(allocf32, arglist));
  PyObject *TSI1    = EXECPYTHON(PyObject_CallObject(allocf32, arglist));
  PyObject *MAT0    = EXECPYTHON(PyObject_CallObject(allocf32, arglist));
  PyObject *MAT1    = EXECPYTHON(PyObject_CallObject(allocf32, arglist));
  PyObject *MAT2    = EXECPYTHON(PyObject_CallObject(allocf32, arglist));
  PyObject *hidx    = EXECPYTHON(PyObject_CallObject(alloci32, arglist));
  PyObject *idx     = EXECPYTHON(PyObject_CallObject(alloci32, arglist));
  
  Py_DECREF(arglist); 

  FILE *fp=fopen(saveval,"r");
  fread(PyArray_DATA((PyArrayObject *)signal),1,sizeof(float)*all_ndata,fp);
  fclose(fp);
  
  sprintf(saveval,"%s_CNNINFO_W_%d",path,rank);
  fp=fopen(saveval,"r");
  fread(PyArray_DATA((PyArrayObject *)weights),1,sizeof(float)*all_ndata,fp);
  fclose(fp);
  
  sprintf(saveval,"%s_CNNINFO_CO_%d",path,rank);
  fp=fopen(saveval,"r");
  fread(PyArray_DATA((PyArrayObject *)TCO1),1,sizeof(float)*all_ndata,fp);
  fclose(fp);
  
  sprintf(saveval,"%s_CNNINFO_SI_%d",path,rank);
  fp=fopen(saveval,"r");
  fread(PyArray_DATA((PyArrayObject *)TSI1),1,sizeof(float)*all_ndata,fp);
  fclose(fp);
  
  sprintf(saveval,"%s_CNNINFO_MAT0_%d",path,rank);
  fp=fopen(saveval,"r");
  fread(PyArray_DATA((PyArrayObject *)MAT0),1,sizeof(float)*all_ndata,fp);
  fclose(fp);
  
  sprintf(saveval,"%s_CNNINFO_MAT1_%d",path,rank);
  fp=fopen(saveval,"r");
  fread(PyArray_DATA((PyArrayObject *)MAT1),1,sizeof(float)*all_ndata,fp);
  fclose(fp);
  
  sprintf(saveval,"%s_CNNINFO_MAT2_%d",path,rank);
  fp=fopen(saveval,"r");
  fread(PyArray_DATA((PyArrayObject *)MAT2),1,sizeof(float)*all_ndata,fp);
  fclose(fp);
  
  sprintf(saveval,"%s_CNNINFO_HIDX_%d",path,rank);
  fp=fopen(saveval,"r");
  fread(PyArray_DATA((PyArrayObject *)hidx),1,sizeof(int)*all_ndata,fp);
  fclose(fp);
  
  sprintf(saveval,"%s_CNNINFO_IDX_%d",path,rank);
  fp=fopen(saveval,"r");
  fread(PyArray_DATA((PyArrayObject *)idx),1,sizeof(int)*all_ndata,fp);
  fclose(fp);
  
  fprintf(stderr,"NNBPIX %ld\n",(long) nnbpix);
  sprintf(saveval,"%s_CNNINFO_RPIX_%d",path,rank);
  fp=fopen(saveval,"r");
  fread(PyArray_DATA((PyArrayObject *)realpix),1,sizeof(int)*nnbpix,fp);
  fclose(fp);
  
  arglist = PyTuple_New(12);
  PyTuple_SetItem(arglist,0,mynetwork);
  PyTuple_SetItem(arglist,1,signal);
  PyTuple_SetItem(arglist,2,weights);
  PyTuple_SetItem(arglist,3,TCO1);
  PyTuple_SetItem(arglist,4,TSI1);
  PyTuple_SetItem(arglist,5,MAT0);
  PyTuple_SetItem(arglist,6,MAT1);
  PyTuple_SetItem(arglist,7,MAT2);
  PyTuple_SetItem(arglist,8,hidx);
  PyTuple_SetItem(arglist,9,idx);
  PyTuple_SetItem(arglist,10,realpix);
  PyTuple_SetItem(arglist,11,randombase);
  
  PyObject *myrun = EXECPYTHON(PyObject_CallObject(init_net, arglist));
  Py_DECREF(arglist); 

  for (itt=0;itt<300;itt++) {
    int k;
    Py_ssize_t nval;
    PyArrayObject *value;

    arglist = PyTuple_New(1);
    PyTuple_SetItem(arglist,0,myrun);
    PyObject *mygradient = EXECPYTHON(PyObject_CallObject(grad, arglist)); 
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
      MPI_Allreduce(PyArray_DATA(value),res,nnn,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
      //fprintf(stderr,"%s %ld %f\n",keystr,(long) nnn,res[0]);
      for (l=0;l<nnn;l++){
	res[l]/=mpi_size;
      }
      memcpy(PyArray_DATA(value),res,sizeof(float)*nnn);
    }


    PyObject *l_arglist = PyTuple_New(2);
    PyTuple_SetItem(l_arglist,0,myrun);
    PyTuple_SetItem(l_arglist,1,mygradient);
    EXECPYTHON(PyObject_CallObject(agrad, l_arglist));

    if (itt%10==0) {
      PyObject *res=EXECPYTHON(PyObject_CallObject(gloss, arglist)); 
      if (res==NULL) {
	PyErr_Print();
      }
      float *loss=PyArray_DATA((PyArrayObject *)res);
      float l_loss[2];
      MPI_Reduce(loss,l_loss,2,MPI_FLOAT,MPI_SUM,0,MPI_COMM_WORLD);
      if (rank==0) {
	fprintf(stderr,"Itt %d loss=%.4f Lr=%.4f\n",(int) itt,sqrt(l_loss[0]/mpi_size),loss[1]);
      }
      Py_DECREF(res);
    }
  }
  
  if (rank==0) {
    PyObject *thepred=EXECPYTHON(PyObject_CallObject(pred, arglist)); 
    float *correction = (float *) malloc(sizeof(float)*8*256*256);
    memcpy(correction,PyArray_DATA((PyArrayObject *)thepred),sizeof(float)*8*256*256);
    Py_DECREF(thepred);
  }

  Py_DECREF(arglist);
  
  // Clean up
  Py_DECREF(Module);
  Py_DECREF(ModuleString);
  
  
  Py_Finalize();
}
