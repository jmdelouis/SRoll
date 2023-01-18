
#define _XOPEN_SOURCE 500

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>
#include <errno.h>
#include <fcntl.h>
#include <math.h>
#include <assert.h>
#include <mpi.h>

#include "stim_tools.h"

#include "no_dmc_data_access.h"
#include "no_dmc_metadata.h"
#include "ring_index.h"


static PIOSTRING ELAPSEDTIME;


////////////////////////////////////////////////////////////////////////////////
// RngStreams helper functions
// create_rngstreams() + delete_rngstreams() takes less than 10ms

// create one RngStream per ring, takes 4.5MB
void create_rngstreams( stimParameters *stimPar) {
  for (int ring=0; ring<MAXRSTREAMS; ring++) {
    stimPar->rstream[ring] = RngStream_CreateStream( NULL);
    assert(stimPar->rstream[ring] != NULL);
  }
}

// delete all RngStreams
void delete_rngstreams( stimParameters *stimPar) {
  for (int ring=0; ring<MAXRSTREAMS; ring++) {
    RngStream_DeleteStream( stimPar->rstream[ring]);
  }
}

// draw a random number with a white noise distribution
float whitenoise( stimParameters *stimPar, int ring) {
/*
  if (stimPar->rstream[ring] == NULL) {
    fprintf( stderr, "rstream=NULL in whitenoise(), rank=%d, ring=%d", stimPar->mpi_rank, ring);
    exit(-1);
  }
*/
  return( sqrt( -2 * log( RANDOM( stimPar, ring))) * cos( 2 * M_PI * RANDOM( stimPar, ring)));
}

////////////////////////////////////////////////////////////////////////////////
// read 4K lines object used in sim_adu.c and correct_adc.c
// 4K lines object example: $DMCDATA/RAW_4K/02_HARMS_LM4K_GP21

int read_4k_vect( char      *objname,   // in: DMC VECT object name
                  int       begin_ring, // in: first ring to read
                  int       ring_count, // in: number of rings to read
                  stim_Data *data)      // out: structure containing 4K lines data
{

  if (strcmp( objname, "0") == 0) {
    return( 0);
  }

  if (data->n_harmonics != 0) {
    return( 0);
  }

  keyword   fourk_keyw;
  PIOSTRING idx_harmonics;

  // read <n_harmonics> and <idx_harmonics> from metadata keywords
  assert( getMetadataKeywordFor( objname, "n_harmonics", &fourk_keyw) == 0);
  data->n_harmonics = atol( fourk_keyw.Val);
  assert( (data->n_harmonics > 0) && (data->n_harmonics < 20));
  assert( getMetadataKeywordFor( objname, "idx_harmonics", &fourk_keyw) == 0);
  strncpy( idx_harmonics, fourk_keyw.Val, PIOSTRINGMAXLEN);

  // decode 4K harmonics values
  data->fourk_harmonics = malloc( data->n_harmonics * sizeof( float));
  assert ( data->fourk_harmonics != NULL);
  char * pchar;
  pchar = strtok ( idx_harmonics, ":");
  int i = 0;
  while (pchar != NULL) {
    assert( i < data->n_harmonics);
    data->fourk_harmonics[i++] = atof( pchar); // * 180.0 / 9.0;
    pchar = strtok( NULL, ":");
  }

  // read 4K lines: one complex number (float, float) per harmonic per ring
  // for interval BeginRing-EndRing
  float *fourk_piofloat = malloc( data->n_harmonics * ring_count * 2 * sizeof( float));
  assert( fourk_piofloat != NULL);
  assert( noDMC_readObject_PIOFLOAT( objname, data->n_harmonics * begin_ring * 2,
                                              data->n_harmonics * ring_count * 2,
                                              fourk_piofloat) >= 0);

  data->fourk_amplitudes = malloc( data->n_harmonics * ring_count * sizeof( float complex));
  assert( data->fourk_amplitudes != NULL);
  for (i = 0; i < data->n_harmonics * ring_count; i++) {
    data->fourk_amplitudes[i] = fourk_piofloat[2*i] + I * fourk_piofloat[2*i+1];
  }
  free( fourk_piofloat);

  return( 0);
}


////////////////////////////////////////////////////////////////////////////////
// compute 4K lines contribution for given ring

void compute_4K_lines( int       ring,             // in: ring index in fourk_amplitudes array
                       long      begin_ring_index, // in: BEGINRINGINDEX of absolute ring number
                       stim_Data *data,            // in: structure containing 4K lines data
                       double    *ring_4k)         // out: 4K lines value per fast sample
{
  for (int isamp = 0; isamp < 9; isamp++) {
    long pp = begin_ring_index + 2*isamp + PARITY( begin_ring_index);
    for (int j = 0; j < 80; j++) {
      // IDL: for i=0,n_elements(fourk_harms)-1 do raw4K+=double(amp[i,ring-ccc])*cos(((dindgen(80))/40d0 + (pp mod 18288))*fourk_harms(i)/9d0*2d0*!dpi) + imaginary(amp[i,ring-ccc])*sin(((dindgen(80))/40d0 + (pp mod 18288))*fourk_harms(i)/9d0*2d0*!dpi)
      ring_4k[isamp*80 + j] = 0.0;
      for (int harm = 0; harm < data->n_harmonics; harm++) {
        ring_4k[isamp*80 + j] += creal( data->fourk_amplitudes[ring * data->n_harmonics + harm]) * cos( ((double)j/40.0 + (pp % 18288)) * data->fourk_harmonics[harm]/9.0*2.0*M_PI)
                               + cimag( data->fourk_amplitudes[ring * data->n_harmonics + harm]) * sin( ((double)j/40.0 + (pp % 18288)) * data->fourk_harmonics[harm]/9.0*2.0*M_PI);
      }
      assert(!isnan(ring_4k[isamp*80 + j]));
    }
  }
}


////////////////////////////////////////////////////////////////////////////////
// find_interval_*: finds the interval from an array where a value lies.
// <xarray> must be monotonic
// use *_increasing or *_decreasing according to <xarray> direction

long find_interval_increasing( double x, double *xarray, long array_size) {

  if (x < xarray[0]) return( 0);
  if (x > xarray[array_size-2]) return( array_size-2);

  long jump_size = array_size / 2;
  long pos = jump_size;
  while (1) {
    if (jump_size > 1) jump_size /= 2;
    if (x < xarray[pos])          pos -= jump_size;
    else if (x > xarray[pos + 1]) pos += jump_size;
    else return( pos);
  }
}

long find_interval_decreasing( double x, double *xarray, long array_size) {

  if (x > xarray[0]) return( 0);
  if (x < xarray[array_size-2]) return( array_size-2);

  long jump_size = array_size / 2;
  long pos = jump_size;
  while (1) {
    if (jump_size > 1) jump_size /= 2;
    if (x > xarray[pos])          pos -= jump_size;
    else if (x < xarray[pos + 1]) pos += jump_size;
    else return( pos);
  }
}


////////////////////////////////////////////////////////////////////////////////
// linearly interpolate (xin, yin) at <xout> abscissae into <yout>
// <xin> must be monotonic
// <xout> may not be monotonic
// <yout> must be allocated by the caller
// if xout[n] is outside of xin[], it's extrapolated from the first or last xin[] interval

double linint( double x, double x1, double y1, double x2, double y2) {
  return( (y2-y1)/(x2-x1)*(x-x1)+y1);
}

void linear_interpolate( double *xin, double *yin, double *xout, double *yout, long n_in, long n_out, int trace) {

  long in_idx=0;

  long (*find_interval)(double, double *, long);

  if (trace) printf("linear interpolate: xin=%p, xout=%p, yin=%p, yout=%p\n", xin, xout, yin, yout);

  if (xin[n_in-1] > xin[0]) {
    find_interval = &find_interval_increasing;
  } else {
    find_interval = &find_interval_decreasing;
  }

  for (long out_idx = 0; out_idx < n_out; out_idx++) {
    in_idx = find_interval( xout[out_idx], xin, n_in);
    assert( in_idx >= 0 && in_idx < n_in-1);
    yout[out_idx] = linint( xout[out_idx], xin[in_idx], yin[in_idx], xin[in_idx+1], yin[in_idx+1]);
  }  
}


////////////////////////////////////////////////////////////////////////////////
// get pysical and virtual memory used

void stim_GetProcMem( long *vmem,long *phymem)
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


////////////////////////////////////////////////////////////////////////////////
// get total, used and free memory in MB as return by the free command, including swap

void stim_FreeMem( long *mem_total, long *mem_used, long *mem_free) {

  *mem_total = 0;
  *mem_used  = 0;
  *mem_free  = 0;
  
  char tmp[256]={0x0};
  FILE *shellcommand = popen( "free -tb | tail -n 1", "r");
  while( fgets( tmp, sizeof(tmp),shellcommand)!=NULL) {
    sscanf( tmp, "Total: %ld %ld %ld", mem_total, mem_used, mem_free);
  }
  pclose( shellcommand);
  // convert from Bytes to MegaBytes
  *mem_total /= 1024*1024;
  *mem_used  /= 1024*1024;
  *mem_free  /= 1024*1024;
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


////////////////////////////////////////////////////////////////////////////////
// get hostname without forking a process to not mess with MPI
// HOSTNAME environment variable always contains the master node HOSTNAME

void GetHostname( PIOSTRING hostname) {

  FILE *f = fopen("/proc/sys/kernel/hostname", "r");
  fscanf( f, "%s", hostname);
  fclose( f);
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
// bin a TOI ring into PBRsize bins

int TOI2PBR( float  *signal,
             PIOINT *flag,
             double *phase,
             long   signal_length,
             int    PBRsize,
             double *pbr_signal,
             int    *pbr_hitcount) {

  long i;
  int bin;

  // initialise output arrays
  for (bin = 0; bin < PBRsize; bin++) {
    pbr_signal[bin] = 0.0;
    pbr_hitcount[bin] = 0;
  }

  // sum signal in bins
  for (i = 0; i < signal_length; i++) {
    if ((flag == NULL) || (flag[i] != -1)) {
      bin = phase[i] * PBRsize / (2 * M_PI);
      pbr_signal[bin]   += signal[i];
      pbr_hitcount[bin] += 1;
    }
  }

  // average signal in bins
  for (bin = 0; bin < PBRsize; bin++) {
    if (pbr_hitcount[bin] > 0) {
      pbr_signal[bin] /= pbr_hitcount[bin];
    } else {
      pbr_signal[bin] = NAN;
    }
  }

  return( 0);
}


////////////////////////////////////////////////////////////////////////////////
// standard deviation and optional mean, ignoring nans

void stddev( float *data, long ndata, double *mean, double *sigma) {

/*
In [1]: import numpy
In [2]: x = [65, 63, 67, 64, 68, 62, 70, 66, 68, 67, 69, 71, 66, 65, 70]
In [3]: numpy.std(x, ddof=1)
Out[3]: 2.6583202716502514
In [4]: numpy.std(x)
Out[4]: 2.5681813712344295
*/

  long i;
  long   count = 0;
  double sum = 0.0;
  double sumsqr = 0.0;
  double  fmean, deviation;
  
  for (i = 0; i < ndata; i++) {
    if (!isnan(data[i])) {
      sum += data[i];
      count += 1;
    }
  }
  fmean = sum / count;

  for (i = 0 ; i < ndata; i++) {
    if (!isnan(data[i])) {
      deviation = data[i] - fmean;
      sumsqr += deviation * deviation;
    }
  }

  if (mean != NULL) *mean = fmean;
  *sigma = sqrt( sumsqr / (count-1));
}


////////////////////////////////////////////////////////////////////////////////
// linear combination of vectors

void lincombvec( float *vector_out, float **vectors_in, float *factors, long vec_length, int vec_count) {

  long  i;
  int   j;
  float tmp;

  assert( vec_count > 0);
  
  for (i = 0; i < vec_length; i++) {
    tmp = 0.0;
    for (j = 0; j < vec_count; j++) {
      tmp += factors[j] * vectors_in[j][i];
    }
    vector_out[i] = tmp;
  }
}

// from http://coding.debuntu.org/c-implementing-str_replace-replace-all-occurrences-substring
/* !!! UNTESTED !!! */

char *str_replace ( const char *string, const char *substr, const char *replacement ){
  char *tok = NULL;
  char *newstr = NULL;
  char *oldstr = NULL;
  char *head = NULL;
 
  newstr = malloc( strlen( string) + 1);
  if ( newstr == NULL ) return NULL;
  strcpy( newstr, string);
  
  /* if either substr or replacement is NULL, duplicate string a let caller handle it */
  if ( substr == NULL || replacement == NULL ) return newstr;
  head = newstr;
  while ( (tok = strstr ( head, substr ))){
    oldstr = newstr;
    newstr = malloc ( strlen ( oldstr ) - strlen ( substr ) + strlen ( replacement ) + 1 );
    /*failed to alloc mem, free old string and return NULL */
    if ( newstr == NULL ){
      free (oldstr);
      return NULL;
    }
    memcpy ( newstr, oldstr, tok - oldstr );
    memcpy ( newstr + (tok - oldstr), replacement, strlen ( replacement ) );
    memcpy ( newstr + (tok - oldstr) + strlen( replacement ), tok + strlen ( substr ), strlen ( oldstr ) - strlen ( substr ) - ( tok - oldstr ) );
    memset ( newstr + strlen ( oldstr ) - strlen ( substr ) + strlen ( replacement ) , 0, 1 );
    /* move back head right after the last replacement */
    head = newstr + (tok - oldstr) + strlen( replacement );
    free (oldstr);
  }
  return newstr;
}


////////////////////////////////////////////////////////////////////////////////

char *get_elapsed_time( struct timeval *t0) {

  struct timeval t1;
  gettimeofday( &t1, NULL);
  time_t ts = (t1.tv_sec - t0->tv_sec) % 60;
  time_t tm = (t1.tv_sec - t0->tv_sec - ts) / 60 % 60;
  time_t th = ((t1.tv_sec - t0->tv_sec - ts) / 60 - tm) / 60;
  sprintf( ELAPSEDTIME, "%02ldh:%02ldm:%02lds ( %.1f sec.)", th, tm, ts, (float)(t1.tv_sec-t0->tv_sec) + (float)(t1.tv_usec-t0->tv_usec)/1.0e6);
  return ELAPSEDTIME;
}

void print_elapsed_time( struct timeval *t0, char *msg_prefix, int trace_level) {

  if (trace_level > 0) {
    fprintf( stderr, "%s %s\n", msg_prefix, get_elapsed_time( t0));
  }
}


////////////////////////////////////////////////////////////////////////////////
// search for a pixname or boloid in a filename and returns the pixname and BC if found,
// return -1 if not found

char *get_pixname( char *objname) {
  for (int i=0; i<NUMBEROFBOLO; i++) {
    if (strstr( objname, PIXNAMES(i)) || strstr( objname, BOLOIDS(i))) {
      return (char *) PIXNAMES(i);
    }
  }
  // not found
  fprintf( stderr, "%s::%s() ERROR: pixname not found in <%s>\n", __FILE__, __FUNCTION__, objname);
  exit( -1);
}

////////////////////////////////////////////////////////////////////////////////
// return 1 if <str> starts with <prefix> else 0

int startswith( char *str, char *prefix) {
  return (strstr( str, prefix) == str);
}

  
////////////////////////////////////////////////////////////////////////////////
// write data to binary file

void writebin( char *filename, void *data, long count, long offset) {

  int fp = open( filename, O_WRONLY|O_CREAT, 0664);
  if (fp == -1) {
    fprintf( stderr, "ERROR %d in open( %s, O_WRONLY|O_CREAT)\n%s\n", errno, filename, strerror( errno));
    exit( -1);
  }
  else {
    size_t ioret = pwrite( fp, data, count, offset);
    if (ioret == -1) {
      fprintf( stderr, "ERROR %d in write( %s)\n%s\n", errno, filename, strerror( errno));
      exit( -1);
    }
    else if (ioret != count) {
      fprintf( stderr, "ERROR: incomplete write( %s), %ld bytes written out of %ld\n%s\n", filename, ioret, count, strerror( errno));
      exit( -1);
    }
    fprintf( stderr, "written %ld bytes to %s\n", ioret, filename);
    close( fp);
  }
}

// change in place a string to lower case
char *lowercase( char *str) {
  for (int i = 0; i<strlen( str); i++){
    str[i] = tolower( str[i]);
  }
  return( str);
}

// return 1 if the given pixname or corresponding boloid is present in the object name
int check_pixname( char *objname, char *pixname) {
  if ( strstr( objname, pixname) || strstr( objname, toBoloID( pixname))) {
    return 1;
  }
  else {
    fprintf( stderr, "%s::%s(): <%s> not found in <%s>\n", __FILE__, __FUNCTION__, pixname, objname);
    return 0;
  }
}


