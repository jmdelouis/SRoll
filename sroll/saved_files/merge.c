#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include "sroll_param.h"
#include <math.h>
#include <sys/time.h>
#include <unistd.h>   
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

int PIOMergeMAP(const char *path,int mpi_size)
{
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
    double *value = (double *) malloc(sizefile);
    int fp=open(thepath,O_RDONLY,0664);
    err=read(fp,value,sizefile);
    close(fp);
     int ii;
    for (ii=0;ii<sizefile/sizeof(double);ii++) map[ii+nnn]=value[ii];
    nnn+=sizefile/sizeof(double);
  }
  
  int fp=open(path,O_WRONLY|O_CREAT,0664);
  err=write(fp,map,12*2048*2048*sizeof(double));
  close(fp);
  free(map);

  sprintf(thepath,"rm -rf %s_dir_%d",path,mpi_size);
  system(thepath);

  return(err);
}



#include "sroll_parLoader.h"

int main(int argc,char *argv[])  
{
  
  sroll_parContent* Param;
  
  fprintf(stderr,"--------------------------\n" );
  fprintf(stderr,"Starting\n" );
  fprintf(stderr, "--------------------------\n" );
  
  /* read params */
  sroll_parContent par;
  int res = sroll_readParam(&par, argv[1] );
  if (res != 0) {
    fprintf(stderr, "Unable to parse the parameter file.\n");
    exit(res);
  }
  Param = &par;

  int nbproc=atoi(argv[2]);
  long i;

  PIOSTRING *mapout = (PIOSTRING *) malloc(sizeof(PIOSTRING)*Param->n_MAP);
  for (i=0;i<Param->n_MAP;i++) strcpy(mapout[i],Param->MAP[i]);

  long maptype;
  long p;
  for (maptype=0;maptype<Param->n_MAP;maptype++) {

    for (p=0;p<6;p++) {
      char suffix[10];
      if (p==0) strcpy(suffix,"II");
      if (Param->OUT_NOPOL[maptype]!=1||p==0) {
        if (p==1) strcpy(suffix,"IQ");
        if (p==2) strcpy(suffix,"IU");
        if (p==3) strcpy(suffix,"QQ");
        if (p==4) strcpy(suffix,"QU");
        if (p==5) strcpy(suffix,"UU");

        PIOSTRING OUTMAP;
        sprintf(OUTMAP,"%s_survey_0_%s",mapout[maptype],suffix);
        PIOMergeMAP(OUTMAP,nbproc);
        for (i=0;i<2;i++) {
          sprintf(OUTMAP,"%s_year_%d_%s",mapout[maptype],(int) i,suffix);
          PIOMergeMAP(OUTMAP,nbproc);
        }
      }
    }
    for (p=0;p<3;p++) {
      char suffix[10];
      if (Param->OUT_NOPOL[maptype]!=1||p==0) {
        if (p==0) strcpy(suffix,"I");
        if (p==1) strcpy(suffix,"Q");
        if (p==2) strcpy(suffix,"U");

        PIOSTRING OUTMAP;
        sprintf(OUTMAP,"%s_corr_survey_0_%s",mapout[maptype],suffix);
        PIOMergeMAP(OUTMAP,nbproc);
        for (i=0;i<2;i++) {
          sprintf(OUTMAP,"%s_corr_year_%d_%s",mapout[maptype],(int) i,suffix);
          PIOMergeMAP(OUTMAP,nbproc);
        }
      }
    }


  }

  fprintf(stderr, "fit_raw: --------------------------\n" );
  fprintf(stderr, "fit_raw: Finished sucessfully\n" );
  fprintf(stderr, "fit_raw: --------------------------\n" );

  exit ( 0 );
}
