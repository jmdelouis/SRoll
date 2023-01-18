#define _GNU_SOURCE

#include <unistd.h>

#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <fcntl.h>
#include <assert.h>
#include <errno.h>

#include "binary_file_toi.h"


/*******************************************************************************
* binary_file_toi.c specifications

- object names are a full path to a filename prefix
  eg: /pscratch1/delouis/SIM_data/dmc/MISS03/DATA/e2e_common_TOI/143-1a_skynops_Oct15

- the object is stored using two files, one with the '.written' suffix,
  one with the '.<datatype>.bin' suffix

- the '.written' file is a binary file containing an ordered list of pairs of PIOLONG (long int) giving
  the first and last sample numbers of gapless written data. When data is written to the TOI,
  the written list is updated and the new written interval is merged with the existing ones,
  maintaining the list ordered.

- the other file contains flat binary data, each sample is written at the offset
  correspondong to its sample number starting from 0. <datatype> can be float32 (PIOFLOAT),
  float64 (PIODOUBLE), int32 (PIOINT) or int64 (PIOLONG).

- when data is read from the TOI, the '.written' file is checked to see if all the requested samples
  are written. If not, an error (assert) occurs.

- when data is written to a TOI, the files are created if they don't exist. There is no
  explicit creation procedure.

- read and write functions have the same parameters as noDMC_readObject_*() for ease of use.

*******************************************************************************/

#define BFT_DEBUG 0


////////////////////////////////////////////////////////////////////////////////

// print chunk list to stdout for debugging
void print_chunk_list( char *print_me_first,
                       written_chunk_struct chunk_list[],
                       int  chunk_count)
{
  fprintf( stdout, print_me_first);
  int chunk_idx;
  for (chunk_idx = 0; chunk_idx < chunk_count; chunk_idx++) {
    fprintf( stdout, "(%d) %ld-%ld\n", chunk_idx, chunk_list[chunk_idx].first_sample, chunk_list[chunk_idx].last_sample);
  }
  if (chunk_count == 0) {
    fprintf( stdout, "empty chunk list\n");
  }
}


////////////////////////////////////////////////////////////////////////////////

long BFT_type_size( const char *data_type) {

  if        (!strcmp( data_type, "float32")) {
    return( 4l);
  } else if (!strcmp( data_type, "float64")) {
    return( 8l);
  } else if (!strcmp( data_type, "byte")) {
    return( 1l);
  } else if (!strcmp( data_type, "int32")) {
    return( 4l);
  } else if (!strcmp( data_type, "int64")) {
    return( 8l);
  }
  printf( "%s: unknow data type: %s\n", __FILE__, data_type);
  exit( -1);
}


////////////////////////////////////////////////////////////////////////////////

// add a chunk to a list of chunks, merging overlapping chunks and keeping the list sorted
int add_chunk_to_list( long first_sample,
                       long last_sample,
                       written_chunk_struct chunk_list[],
                       int  chunk_count) {

  if (BFT_DEBUG) {
    fprintf( stdout, "add_chunk_to_list(): adding %ld-%ld to\n", first_sample, last_sample);
    print_chunk_list( NULL, chunk_list, chunk_count);
  }

  if (chunk_count == 0) {
    // chunk list is empty, adding a new one
    chunk_list[0].first_sample = first_sample;
    chunk_list[0].last_sample  = last_sample;
    chunk_count = 1;
    return chunk_count;
  }

  // first_sample or last_sample either fall inside an existing chunk or between two chunks
  // that's 4 cases to manage
  int first_inside_chunk = -1;          // fisrt_sample is inside this chunk number
  int first_before_chunk = chunk_count; // fisrt_sample is before this chunk number
  int last_inside_chunk  = -1;          // last_sample is inside this chunk number
  int last_before_chunk  = chunk_count; // last_sample is before this chunk number
  int chunk_idx;

  for (chunk_idx = 0; chunk_idx < chunk_count; chunk_idx++) {
    if ((first_sample < chunk_list[chunk_idx].first_sample) &&
        (chunk_idx < first_before_chunk)) {
      first_before_chunk = chunk_idx;
    }
    if ((first_sample >= chunk_list[chunk_idx].first_sample) &&
        (first_sample <= chunk_list[chunk_idx].last_sample + 1)) {
      first_inside_chunk = chunk_idx;
    }
    if ((last_sample < chunk_list[chunk_idx].first_sample) &&
        (chunk_idx < last_before_chunk)) {
      last_before_chunk = chunk_idx;
    }
    if ((last_sample >= chunk_list[chunk_idx].first_sample - 1) &&
        (last_sample <= chunk_list[chunk_idx].last_sample)) {
      last_inside_chunk = chunk_idx;
    }
  }

  if (BFT_DEBUG) {
    fprintf( stdout, "first_inside_chunk = %d\n", first_inside_chunk);
    fprintf( stdout, "first_before_chunk = %d\n", first_before_chunk);
    fprintf( stdout, "last_inside_chunk  = %d\n",  last_inside_chunk);
    fprintf( stdout, "last_before_chunk  = %d\n",  last_before_chunk);
  }

  if (first_inside_chunk == -1) {
    if ((last_inside_chunk == -1) && (first_before_chunk == last_before_chunk)) {
      // this is the only case where we just add a chunk to the list:
      // chunk from first_sample to last_sample has no intersection with chunk list
      for (chunk_idx = chunk_count; chunk_idx > first_before_chunk; chunk_idx--) {
        assert( (chunk_idx >= 1) && (chunk_idx <= chunk_count));
        chunk_list[chunk_idx].first_sample = chunk_list[chunk_idx - 1].first_sample;
        chunk_list[chunk_idx].last_sample  = chunk_list[chunk_idx - 1].last_sample;
      }
      chunk_list[first_before_chunk].first_sample = first_sample;
      chunk_list[first_before_chunk].last_sample  = last_sample;
      chunk_count += 1;
      return chunk_count;
    }
    // extend first_before_chunk to include first_sample
    chunk_list[first_before_chunk].first_sample = first_sample;
    first_inside_chunk = first_before_chunk;
  }

  // first sample is now inside first_inside_chunk
  int first_chunk_todel, last_chunk_todel;
  if (last_inside_chunk != -1) {
    // merge chunks first_inside_chunk and last_inside_chunk,
    chunk_list[first_inside_chunk].last_sample = chunk_list[last_inside_chunk].last_sample;
    first_chunk_todel = first_inside_chunk + 1;
    last_chunk_todel  = last_inside_chunk;
  } else {
    chunk_list[first_inside_chunk].last_sample = last_sample;
    first_chunk_todel = first_inside_chunk + 1;
    last_chunk_todel  = last_before_chunk - 1;
  }

  int chunks_to_del = last_chunk_todel - first_chunk_todel + 1;
  if (chunks_to_del > 0) {
    for (chunk_idx = last_chunk_todel + 1; chunk_idx < chunk_count; chunk_idx++) {
      chunk_list[chunk_idx - chunks_to_del].first_sample = chunk_list[chunk_idx].first_sample;
      chunk_list[chunk_idx - chunks_to_del].last_sample  = chunk_list[chunk_idx].last_sample;
    }
    chunk_count -= chunks_to_del;
  }
  return chunk_count;
}


////////////////////////////////////////////////////////////////////////////////

// read BinaryFileToi data, returns 0 if ok, failed assert if nok
int BFT_readObject( const char *object_name,
                    const char *data_type,
                    long first_sample,
                    long sample_count,
                    void *data)
{
#if 0
  long last_sample = first_sample + sample_count - 1;
  int chunk_idx;

  // get a file descriptor on .written file
  PIOSTRING written_filename;
  snprintf( written_filename, PIOSTRINGMAXLEN, "%s.written", object_name);
  int written_filedesc = open(written_filename, O_RDONLY);
  assert( written_filedesc >= 0);

  // get and check .written file size
  struct stat written_stat;
  memset( &written_stat, 0, sizeof( written_stat));
  assert( fstat( written_filedesc, &written_stat) != -1);
  assert( written_stat.st_size % sizeof( written_chunk_struct) == 0);

  // read chunk list from .written file
  int chunk_count = written_stat.st_size / sizeof( written_chunk_struct);
  written_chunk_struct *chunk_list = NULL;
  chunk_list = (written_chunk_struct *) malloc(sizeof(written_chunk_struct) * chunk_count);
  assert( chunk_list != NULL);
  assert( pread( written_filedesc, chunk_list, sizeof(written_chunk_struct) * chunk_count, 0) == sizeof(written_chunk_struct) * chunk_count);

  // check that first_sample and last_sample are in the same chunk of gapless data
  for (chunk_idx = 0; chunk_idx < chunk_count; chunk_idx++) {
    if (first_sample >= chunk_list[chunk_idx].first_sample) {
      assert( last_sample <= chunk_list[chunk_idx].last_sample);
      break;
    }
  }
  // if no gapless data chunk found, some requested samples are not written, fail.
  assert( chunk_idx < chunk_count);
#endif
  // read binary file
  PIOSTRING bin_filename;
  if (strstr( object_name, ".bin") != NULL) {
    snprintf( bin_filename, PIOSTRINGMAXLEN, "%s", object_name);
  } else {
    snprintf( bin_filename, PIOSTRINGMAXLEN, "%s.%s.bin", object_name, data_type);
  }
  int bin_filedesc = open( bin_filename, O_RDONLY);
  assert( bin_filedesc >= 0);
  long bytespersample = BFT_type_size( data_type);
  assert( pread64( bin_filedesc, data, sample_count * bytespersample, first_sample * bytespersample) == sample_count * bytespersample); 
  close( bin_filedesc);

  return 0;
}


////////////////////////////////////////////////////////////////////////////////

// write BinaryFileToi data, returns sample_count (>=0) if ok, -1 if nok
int BFT_writeObject( const char *object_name,
                     const char *data_type,
                     long first_sample,
                     long sample_count,
                     void *data)
{
  // get a file descriptor on .written file, create the file if it doesn't exist
  PIOSTRING written_filename;
  snprintf( written_filename, PIOSTRINGMAXLEN, "%s.written", object_name);
  int written_filedesc = open( written_filename, O_RDWR | O_CREAT, S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH);
  assert( written_filedesc >= 0);

  // lock .written file
  struct flock lock;
  memset( &lock, 0, sizeof( lock));
  lock.l_type = F_WRLCK;
  assert( fcntl( written_filedesc, F_SETLKW, &lock) != -1);

  // get and check .written file size
  struct stat written_stat;
  memset( &written_stat, 0, sizeof( written_stat));
  assert( fstat( written_filedesc, &written_stat) != -1);
  assert( written_stat.st_size % sizeof( written_chunk_struct) == 0);

  // allocate one more chunk than are present in the file to have space to add one
  int chunk_count = written_stat.st_size / sizeof( written_chunk_struct);
  written_chunk_struct *chunk_list = NULL;
  chunk_list = (written_chunk_struct *) malloc(sizeof(written_chunk_struct) * (chunk_count + 1));
  assert( chunk_list != NULL);

  // read chunk list from .written file
  assert( pread( written_filedesc, chunk_list, sizeof(written_chunk_struct) * chunk_count, 0) == sizeof(written_chunk_struct) * chunk_count);

  // add written chunk to chunk list
  chunk_count = add_chunk_to_list( first_sample, first_sample + sample_count - 1, chunk_list, chunk_count);

  // overwrite new chunk list to .written file
  assert( ftruncate( written_filedesc, 0) == 0);
  assert( pwrite( written_filedesc, chunk_list, sizeof(written_chunk_struct) * chunk_count, 0) == sizeof(written_chunk_struct) * chunk_count);

  // write data to binary file (create it if it doesn't exist)
  PIOSTRING bin_filename;
  snprintf( bin_filename, PIOSTRINGMAXLEN, "%s.%s.bin", object_name, data_type);
  int bin_filedesc = open( bin_filename, O_RDWR | O_CREAT, S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH);
  assert( bin_filedesc >= 0);
  long bytespersample = BFT_type_size( data_type);
  size_t pwrite_result = pwrite64( bin_filedesc, data, sample_count * bytespersample, first_sample * bytespersample);
  if (pwrite_result == -1) {
    fprintf( stderr, "ERROR %d in pwrite64() in BFT_writeObject(filename=%s, count=%ld, offset=%ld)\n", errno, object_name, sample_count * bytespersample, first_sample * bytespersample);
  } else if (pwrite_result != sample_count * bytespersample) {
    fprintf( stderr, "ERROR: only %ld bytes written out of %ld in pwrite64() in BFT_writeObject(filename=%s, count=%ld, offset=%ld)\n", pwrite_result, sample_count * bytespersample, object_name, sample_count * bytespersample, first_sample * bytespersample);
  }

  // clean and leave
  close( bin_filedesc);
  free( chunk_list);
  close( written_filedesc);

  if (pwrite_result != sample_count * bytespersample) {
    return -1;
  } else {
    return pwrite_result;
  }
}

