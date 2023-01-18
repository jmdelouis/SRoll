#ifndef _STIMTOOLS_H_
#define _STIMTOOLS_H_


#include <stdlib.h>
#include <sys/time.h>
#include <complex.h>

#include "no_dmc_data_access.h"
#include "stim_param.h"

#define STIM_TRACE( tl, msg, ...) {if (stimPar->trace_level >= (tl)) fprintf( stderr, "%s %s (line %d) - "msg"\n", stimPar->msg_prefix, __FILE__, __LINE__, ##__VA_ARGS__);}

#define RANDOM( stimPar, ring) (RngStream_RandU01( stimPar->rstream[ring]))

// create one RngStream per ring, takes 4.5MB
void create_rngstreams( stimParameters *stimPar);

// delete all RngStreams
void delete_rngstreams( stimParameters *stimPar);

// draw a random number with a white noise distribution
float whitenoise( stimParameters *stimPar, int ring);


// read 4K lines object used in sim_adu.c and correct_adc.c
// fourk_harmonics and fourk_amplitudes are allocated by the function and must be freed by the caller
int read_4k_vect( char      *objname,   // in: DMC VECT object name
                  int       begin_ring, // in: first ring to read
                  int       ring_count, // in: number of rings to read
                  stim_Data *data);     // out: structure containing 4K lines data


// compute 4K lines contribution for given ring
void compute_4K_lines( int       ring,             // in: ring index in fourk_amplitudes array
                       long      begin_ring_index, // in: BEGINRINGINDEX of absolute ring number
                       stim_Data *data,            // in: structure containing 4K lines data
                       double    *ring_4k);        // out: 4K lines value per fast sample

// linearly interpolate (xin, yin) at xout abscissae into yout
void linear_interpolate( double *xin,  // in: abscissae of function to interpolate
                         double *yin,  // in: values of function to interpolate
                         double *xout, // in: abscissae where to interpolate function
                         double *yout, // out: interpolated function values
                         long n_in,    // in: number of values in xin and yin
                         long n_out,   // in: number of values in xout and yout
                         int  trace);  // in: set to 1 to add debug traces to stdout


// get pysical and virtual memory used
void stim_GetProcMem( long *vmem, long *phymem);


// get total, used and free memory in MB as return by the free command, including swap
void stim_FreeMem( long *mem_total, long *mem_used, long *mem_free);

// get free memory in GB not forking a new process
float GetFreeMemGB();

// get host name not forking a new process
void GetHostname( PIOSTRING hostname);

// print free memory on each MPI node
void PrintFreeMemOnNodes( int mpi_rank, int mpi_size, char* msg);

// bin a TOI ring into PBRsize bins
int TOI2PBR( float  *signal,        // in: signal
             PIOINT *flag,          // in: hpr idx, =-1 when flagged. If NULL, all signal samples are used
             double *phase,         // in: TOI containing the signal phase between 0 and 2PI
             long   signal_length,  // in: number of samples in signal
             int    PBRsize,        // in: number of bins in output
             double *pbr_signal,    // out: allocated phase binned ring
             int    *pbr_hitcount); // out: allocated hit count per bin

// standard deviation
void stddev( float *data, long ndata, double *mean, double *sigma);

// linear combination of vectors
void lincombvec( float *vector_out, float **vectors_in, float *factors, long vec_length, int vec_count);

// string replace
char *str_replace( const char *string, const char *substr, const char *replacement);

// print elpased time between two timestamps
char *get_elapsed_time( struct timeval *t0);
void print_elapsed_time( struct timeval *t0, char *msg_prefix, int trace_level);

// return a pixname or boloid found in an objectname, exit on failure
char *get_pixname( char *objname);

// return 1 if <str> starts with <prefix> else 0
int startswith( char *str, char *prefix);

// write data to binary file
void writebin( char *filename, void *data, long count, long offset);

// change in place a string to lower case
char *lowercase( char *str);

// return 1 if the given pixname or corresponding boloid is present in the object name
int check_pixname( char *objname, char *pixname);

#endif
