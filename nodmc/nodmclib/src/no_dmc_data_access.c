/******************************************************************************
 * "no_dmc_data_access.c"
 *
 * author:  Christan Madsen
 * date:    2014-12-22 (initial)
 * version: STABLE
 *****************************************************************************/

// Used for strndup() definition
// Note that the first define is for glibc < 2.10 and second define is for glibc >= 2.10
// RMK: on magique 3 we have: glibc 2.4
#define _GNU_SOURCE 1
#define POSIX_C_SOURCE 200809L

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#include "no_dmc_data_access.h"
#include "no_dmc_metadata.h"
#include "no_dmc_piolib_type_def.h"
#include "no_dmc_util.h"
#include "no_dmc_debug.h"
#include "fitsio2.h"


// Chunk definition constants (see <object>.0.pio)
#define PIOCHUNKINFOSIZE (32)

#define NON_EXISTING_CHUNK_INDEX (-1)

// In case of non existing chunk we use the following default value to fill data (depending on type):
// Note for PIOCOMPLEX type: since this is 2 floats we use the default float value
#define DEFAULT_DATA_VALUE_BYTE (0)
#define DEFAULT_DATA_VALUE_FLAG (0) // TODO: to see if relevant!???
#define DEFAULT_DATA_VALUE_SHORT (0)
#define DEFAULT_DATA_VALUE_INT (0)
#define DEFAULT_DATA_VALUE_LONG (0)
#define DEFAULT_DATA_VALUE_FLOAT (0.)
#define DEFAULT_DATA_VALUE_DOUBLE (0.)
#define DEFAULT_DATA_VALUE_STRING ""

// Flag context
#define DEFAULT_FLAG_VALUE (0)


// Used to check an error. If the condition in first argument is true, the program is exited
// and shows the message in second argument.
#define ERRORCHECK( test, msg, ...) {if (test) {fprintf( stderr, "** ERROR ** in %s calling '%s' (line %d): "msg"\n", __FILE__, __FUNCTION__, __LINE__, ##__VA_ARGS__);return -1;}}

// Enum used by the generic read function to know what to read: data or flag content.
typedef enum {READ_DATA, READ_FLAG_WRITTEN} whatToRead;


// Structures used to store info read from index file (<object>.0.pio)
// The file is composed of a header and a variable number of chunks.

// Note: Each field of this structure is stored as a 8 bytes value (see doc)
typedef struct {
  PIOLONG chunk_num;
  PIOLONG version_id;
  PIOLONG filenum;
  PIOLONG chunk_ord;
} chunk_index;

typedef struct {
  char          *header;
  /* The header seems to be composed as follow (unit is bytes, 1 indexed) length is inside perenthesis:
  1-20   (20) : obj identifier
  21     (1)  : separator ('-')
  22-29  (8)  : file identifier
  30     (1)  : separator ('-')
  31-54  (24) : obj type
  55     (1)  : separator ('-')
  56-65  (10) : data type + length
  66     (1)  : separator ('-')
  67-86  (20) : database name
  87     (1)  : separator ('-')
  88-128 (41) : spare data...
  */
  char          *data_type_str; // The label associated to the data type. ex: "DOUBLE", "FLOAT", "INT"...
  PIOSHORT       data_type_size; // size in Bytes for the data type considered. ex: for PIODOUBLE this is 8, for PIOFLOAT this is 4...
  chunk_index   *chunk_index_list; // List of chunk_index
  PIOLONG        chunk_index_list_size;
} pio_chunk_indexes;


typedef struct {
  PIOLONG *indexes_in_chunk_index_list; // indexes in the 'chunk_index_list' of the 'pio_chunk_indexes'. A negative value indicate a non existing chunk def.
  PIOLONG indexes_in_chunk_index_list_size;
  PIOLONG first_chunk_offset; // Offset from which we must read in the first chunk (= what we must NOT read)
  PIOLONG last_chunk_limit_offset; // The limit offset UP TO which we must read in the last chunk (= what we must read, included)
} ordered_chunk_list_to_be_read;


/*
 * Init structure of type 'pio_chunk_indexes'.
 */
int initPioChunkIndexes(pio_chunk_indexes *pio_chunk_indexes_p) {
  debug_print("* %s()\n", __FUNCTION__);
  if (pio_chunk_indexes_p == NULL) {
    fprintf(stderr, "ERROR: unable to init pio_chunk_indexes_p: NULL pointer!");
    return -1;
  }

  pio_chunk_indexes_p->header = NULL;
  pio_chunk_indexes_p->data_type_str = NULL;
  pio_chunk_indexes_p->data_type_size = 0;
  pio_chunk_indexes_p->chunk_index_list = NULL;
  pio_chunk_indexes_p->chunk_index_list_size = 0;

  return 0;
}


/*
 * Free memory for structure of type 'pio_chunk_indexes'.
 */
int freePioChunkIndexes(pio_chunk_indexes *pio_chunk_indexes_p) {
  debug_print("* %s()\n", __FUNCTION__);
  if (pio_chunk_indexes_p == NULL) {
    fprintf(stderr, "ERROR: unable to free pio_chunk_indexes_p: NULL pointer!");
    return -1;
  }

  if (pio_chunk_indexes_p->header != NULL) {
    free(pio_chunk_indexes_p->header);
    pio_chunk_indexes_p->header = NULL;
  }
  if (pio_chunk_indexes_p->data_type_str != NULL) {
    free(pio_chunk_indexes_p->data_type_str);
    pio_chunk_indexes_p->data_type_str = NULL;
  }
  if (pio_chunk_indexes_p->chunk_index_list != NULL) {
    free(pio_chunk_indexes_p->chunk_index_list);
    pio_chunk_indexes_p->chunk_index_list = NULL;
    pio_chunk_indexes_p->chunk_index_list_size = 0;
  }

  return 0;
}


/*
 * Init for structure of type 'ordered_chunk_list_to_be_read'.
 */
int initOrderedChunkListToBeRead(ordered_chunk_list_to_be_read *ordered_chunk_list_to_be_read_p) {
  if (ordered_chunk_list_to_be_read_p == NULL) {
    fprintf(stderr, "ERROR: unable to init ordered_chunk_list_to_be_read_p: NULL pointer!");
    return -1;
  }

  ordered_chunk_list_to_be_read_p->indexes_in_chunk_index_list = NULL;
  ordered_chunk_list_to_be_read_p->indexes_in_chunk_index_list_size = 0;
  ordered_chunk_list_to_be_read_p->first_chunk_offset = 0;
  ordered_chunk_list_to_be_read_p->last_chunk_limit_offset = 0;

  return 0;
}


/*
 * Free memory for structure of type 'ordered_chunk_list_to_be_read'.
 */
int freeOrderedChunkListToBeRead(ordered_chunk_list_to_be_read *ordered_chunk_list_to_be_read_p) {
  if (ordered_chunk_list_to_be_read_p == NULL) {
    fprintf(stderr, "ERROR: unable to free ordered_chunk_list_to_be_read_p: NULL pointer!");
    return -1;
  }

  if (ordered_chunk_list_to_be_read_p->indexes_in_chunk_index_list != NULL) {
    free(ordered_chunk_list_to_be_read_p->indexes_in_chunk_index_list);
    ordered_chunk_list_to_be_read_p->indexes_in_chunk_index_list = NULL;
    ordered_chunk_list_to_be_read_p->indexes_in_chunk_index_list_size = 0;
  }

  return 0;
}


/*
 * Read the <object>.0.pio (chunk index file) corresponding to requested DMC object.
 * Note: chunk index file is describe at page 48 of "HFI_DMC_ADD.doc".
 * This function rely on the information present in the previously generated
 * 'no dmc metadata' to find the chunk definition file (<object>.0.pio).
 * Since this function only read indexes it is not critical for performance
 * issue, so that it can use "fread()" function.
 */
int getChunkIndexInfo(const metadata *metadata_p, pio_chunk_indexes *current_pio_chunk_indexes_p) {
  debug_print("* %s()\n", __FUNCTION__);

  char *backendname     = metadata_p->Backendname;
  PIOLONG headerSize    = metadata_p->Iooffset;
  PIOLONG flagchunkSize = metadata_p->Flagchunksize;

  debug_print("  backendname   = %s\n", backendname);
  debug_print("  headerSize    = "PIOLONG_FMT"\n", headerSize);
  debug_print("  flagchunkSize = "PIOLONG_FMT"\n", flagchunkSize);

  current_pio_chunk_indexes_p->header = calloc(MAX_STR_LEN, sizeof(char));
  if (current_pio_chunk_indexes_p->header == NULL) {
    perror("Error");
    return -1;
  }

  // Retrieve the full path to object file
  debug_print("Retrieve the fullpath for the INDEX file ('%s')\n", INDEX_FILE_SUFFIX);
  char index_fullpath[MAX_STR_LEN];
  if (getFullpathToFileFromBackendname(backendname, INDEX_FILE_SUFFIX, 0, index_fullpath) != 0) {
    fprintf(stderr, "ERROR: unable to retrieve fullpath to object file '%s'\n", INDEX_FILE_SUFFIX);
    freePioChunkIndexes(current_pio_chunk_indexes_p);
    return -1;
  }

  // Open the file and extract info from it
  FILE *file = fopen(index_fullpath, "rb" );
  if (file == NULL) {
    fprintf(stderr, "ERROR while trying to open CHUNK INDEX file: '%s'\n", index_fullpath);
    perror("Error");
    freePioChunkIndexes(current_pio_chunk_indexes_p);
    return -1;
  }

  // a) the header (occupy the first 'headerSize' bytes)
  fgets(current_pio_chunk_indexes_p->header, headerSize, file);
  debug_print("    header         = '%s'\n", current_pio_chunk_indexes_p->header);
  // Extract from the header some info
  // (See description of .pio 'header' for index)
  char *str_data_type_str = strndup(current_pio_chunk_indexes_p->header + 55, 6);
  if (str_data_type_str == NULL) {
    perror("Error");
    fclose(file);
    return -1;
  }
  debug_print("    str_data_type_str = '%s'\n", str_data_type_str);
  current_pio_chunk_indexes_p->data_type_str = trim(str_data_type_str);
  debug_print("    data_type_str     = '%s'\n", current_pio_chunk_indexes_p->data_type_str);

  char *str_data_type_size = strndup(current_pio_chunk_indexes_p->header + 62, 3);
  if (str_data_type_size == NULL) {
    perror("Error");
    fclose(file);
    return -1;
  }
  debug_print("    str_data_type_size = '%s'\n", str_data_type_size);
  current_pio_chunk_indexes_p->data_type_size = strtol(str_data_type_size, NULL, 10);
  debug_print("    data_type_size = %d\n", current_pio_chunk_indexes_p->data_type_size);
  free(str_data_type_size); // sub string is no more usefull

  // b) the chunk_index info
  // Guess the number of chunk def from the file size. Note that each chunk info are stored on 32 bytes.
  // So the total file size correspond to [ headerSize + n * 32 ] where 'n' is the number of chunk def.
  off_t tot = getSizeOfFile(file);
  current_pio_chunk_indexes_p->chunk_index_list_size = (tot - headerSize) / PIOCHUNKINFOSIZE;
  debug_print("    size of 0.pio  = %zd\n", tot);
  debug_print("    nb chunk info  = "PIOLONG_FMT"\n", current_pio_chunk_indexes_p->chunk_index_list_size);

  // Allocate memory accordingly
  current_pio_chunk_indexes_p->chunk_index_list = malloc(current_pio_chunk_indexes_p->chunk_index_list_size * sizeof(chunk_index));
  if (current_pio_chunk_indexes_p->chunk_index_list == NULL) {
    perror("Error");
    fclose(file);
    return -1;
  }

  // Start to read just after the header
  if (fseek(file, headerSize, SEEK_SET) != 0) {
    perror("Error");
    fclose(file);
    return -1;
  }

  // Read all chunk indexes in a uniq call
  if (fread(current_pio_chunk_indexes_p->chunk_index_list, sizeof(chunk_index),
        current_pio_chunk_indexes_p->chunk_index_list_size, file) != current_pio_chunk_indexes_p->chunk_index_list_size) {
    perror("Error");
    fclose(file);
    return -1;
  }

  // DEBUG print only (can be deleted or commented!)
  /*
  if (DEBUG) {
    int i;
    // Loop on chunk index def until we reach the end
    for (i=0; i < current_pio_chunk_indexes_p->chunk_index_list_size; i++) {
      debug_print("    chunk_index def (iter %d)\n", i);

      debug_print("      Chunk.chunk_num  = "PIOLONG_FMT"\n", current_pio_chunk_indexes_p->chunk_index_list[i].chunk_num);
      debug_print("      Chunk.version_id = "PIOLONG_FMT"\n", current_pio_chunk_indexes_p->chunk_index_list[i].version_id);
      debug_print("      Chunk.filenum    = "PIOLONG_FMT"\n", current_pio_chunk_indexes_p->chunk_index_list[i].filenum);
      debug_print("      Chunk.chunk_ord  = "PIOLONG_FMT"\n", current_pio_chunk_indexes_p->chunk_index_list[i].chunk_ord);
    }
  }
  */

  // Closing file
  fclose(file);

  return 0;
}


/*
 * This function allow to compute the ordered list of chunk that will need to be read (regarding user request).
 * Note that for practicity the result is gathered in a structure of type: 'ordered_chunk_list_to_be_read'.
 */
int getOrderedChunkListToBeRead(const metadata *metadata_p, PIOLONG offset, PIOLONG nbsample, const pio_chunk_indexes current_pio_chunk_indexes, ordered_chunk_list_to_be_read *output_ordered_chunk_list_to_be_read) {
  debug_print("* %s()\n", __FUNCTION__);

  debug_print("  offset      = "PIOLONG_FMT"\n", offset);
  debug_print("  nbsample    = "PIOLONG_FMT"\n", nbsample);
  PIOLONG flagchunksize = metadata_p->Flagchunksize;

  // 0) Compute/Retrieve required info
  // Number of sample in a complete chunk
  PIOLONG nb_sample_per_chunk = 8 * flagchunksize;
  debug_print("  nb_sample_per_chunk     = "PIOLONG_FMT"\n", nb_sample_per_chunk);

  // a) Determine first chunk and the corresponding offset
  PIOLONG first_chunk = offset / nb_sample_per_chunk; // This value is 0 indexed!
  output_ordered_chunk_list_to_be_read->first_chunk_offset = offset % nb_sample_per_chunk;
  debug_print("  first_chunk             = "PIOLONG_FMT"\n", first_chunk);
  debug_print("  first_chunk_offset      = "PIOLONG_FMT"\n", output_ordered_chunk_list_to_be_read->first_chunk_offset);

  // b) Determine last chunk (with its limit offset)
  PIOLONG last_chunk = (offset+nbsample-1) / nb_sample_per_chunk; // This value is 0 indexed!
  output_ordered_chunk_list_to_be_read->last_chunk_limit_offset = (offset+nbsample-1) % nb_sample_per_chunk; // This value is 0 indexed!
  debug_print("  last_chunk              = "PIOLONG_FMT"\n", last_chunk);
  debug_print("  last_chunk_limit_offset = "PIOLONG_FMT"\n", output_ordered_chunk_list_to_be_read->last_chunk_limit_offset);

  // c) Allocate appropriate memory to store indexes
  PIOLONG total_nb_chunk_to_read = last_chunk - first_chunk + 1;
  debug_print("  total_nb_chunk_to_read  = "PIOLONG_FMT"\n", total_nb_chunk_to_read);
  output_ordered_chunk_list_to_be_read->indexes_in_chunk_index_list = malloc(total_nb_chunk_to_read * sizeof(PIOLONG));
  if (output_ordered_chunk_list_to_be_read->indexes_in_chunk_index_list == NULL) {
    perror("Error");
    return -1;
  }
  output_ordered_chunk_list_to_be_read->indexes_in_chunk_index_list_size = 0; // Still there is no data in the list

  // d) Loop between first and last chunk_num to retrieve the corresponding indexes in 'pio_chunk_indexes'
  PIOLONG chunk_num;
  for (chunk_num = first_chunk; chunk_num <= last_chunk; chunk_num++) {
    debug_print("  seeking for chunk_num "PIOLONG_FMT" in 0.pio indexes...\n", chunk_num);
    // Loop in pio_chunk_indexes in order to find the corresponding index
    PIOLONG i;
    int found = 0;
    for (i = 0; i < current_pio_chunk_indexes.chunk_index_list_size; i++) {
      if (chunk_num == current_pio_chunk_indexes.chunk_index_list[i].chunk_num) {
        found = 1;
        debug_print("    found at index "PIOLONG_FMT"\n", i);
        output_ordered_chunk_list_to_be_read->indexes_in_chunk_index_list[output_ordered_chunk_list_to_be_read->indexes_in_chunk_index_list_size] = i;
        // Update size
        output_ordered_chunk_list_to_be_read->indexes_in_chunk_index_list_size++;
        break;
      }
    }
    if (found != 1) { // Case of non existing chunk num (-> we will fill returned data with default value but for now we just identify inexistance with a negativ index)
      debug_print("    Unable to find chunk_num "PIOLONG_FMT" in 0.pio file -> non existing but will be handle as 0.\n", chunk_num);
      output_ordered_chunk_list_to_be_read->indexes_in_chunk_index_list[output_ordered_chunk_list_to_be_read->indexes_in_chunk_index_list_size] = NON_EXISTING_CHUNK_INDEX;
      // Update size
      output_ordered_chunk_list_to_be_read->indexes_in_chunk_index_list_size++;
    }
  }

  debug_print("    Total chunk_num found = "PIOLONG_FMT"\n", output_ordered_chunk_list_to_be_read->indexes_in_chunk_index_list_size);

  return 0;
}


/*
 * Read the required chunks according to user request and store result in 'data' parameter.
 * Note that due to offset and request sample the first/last chunks might NOT be read entirely!
 * When a chunk does not exist then we fill data with a default value (see DEFAULT_DATA_VALUE_* where "*" is the value type).
 */
int readChunks(const metadata *metadata_p, const pio_chunk_indexes *chunk_ind, const ordered_chunk_list_to_be_read *chunk_list, void *data) {
  debug_print("* %s()\n", __FUNCTION__);
  PIOLONG headerSize  = metadata_p->Iooffset;

  PIOLONG chunkLength_s = metadata_p->Flagchunksize * 8; // Length of a chunk (in number of samples)
  debug_print("  chunkLength_s = "PIOLONG_FMT"\n", chunkLength_s);

  debug_print("  Start reading chunk(s)...\n");
  PIOLONG cur_pos = 0; // Current position in output data (index = sample index)
  PIOLONG i;
  for (i = 0; i < chunk_list->indexes_in_chunk_index_list_size; i++) {
    PIOLONG chunkIndexInIndexList = chunk_list->indexes_in_chunk_index_list[i];
    debug_print("   >chunk index pos: "PIOLONG_FMT"\n", i);

    // 0- Determine 'start_s' and 'length_s' of data to be read. Values are in number of samples.
    PIOLONG start_s  = 0;
    PIOLONG length_s = chunkLength_s; // By default the whole chunk
    // Handle case of FIRST chunk
    if (i == 0) {
      start_s = chunk_list->first_chunk_offset;
    }
    // Handle case of LAST chunk
    if (i == (chunk_list->indexes_in_chunk_index_list_size - 1)) {
      length_s = (chunk_list->last_chunk_limit_offset+1); // +1 since index is included
    }

    PIOLONG size_to_read_s = length_s - start_s; // Size to read in number of samples

    debug_print("    start_s       : "PIOLONG_FMT"\n", start_s);
    debug_print("    length_s      : "PIOLONG_FMT"\n", length_s);
    debug_print("    size_to_read_s: "PIOLONG_FMT"\n", size_to_read_s);

    // We must handle the case of non existing data (when index is NON_EXISTING_CHUNK_INDEX)
    if (chunkIndexInIndexList == NON_EXISTING_CHUNK_INDEX) {
      PIOLONG j;
      // Fill with default value for required length
      debug_print("    Non existing chunk -> fill '"PIOLONG_FMT"' samples with default value.\n", size_to_read_s);

      // IMPORTANT: here we are specific to the required data type (ex: PIODOUBLE, PIOFLOAT...)
      if (strcmp(chunk_ind->data_type_str, "DOUBLE") == 0) {
        for (j = 0; j < size_to_read_s; j++) {
          ((PIODOUBLE *)data)[cur_pos+j] = DEFAULT_DATA_VALUE_DOUBLE;
        }
      } else if (strcmp(chunk_ind->data_type_str, "FLOAT") == 0) {
        for (j = 0; j < size_to_read_s; j++) {
          ((PIOFLOAT *)data)[cur_pos+j] = DEFAULT_DATA_VALUE_FLOAT;
        }
      } else if (strcmp(chunk_ind->data_type_str, "INT") == 0) {
        for (j = 0; j < size_to_read_s; j++) {
          ((PIOINT *)data)[cur_pos+j] = DEFAULT_DATA_VALUE_INT;
        }
      } else if (strcmp(chunk_ind->data_type_str, "LONG") == 0) {
        for (j = 0; j < size_to_read_s; j++) {
          ((PIOLONG *)data)[cur_pos+j] = DEFAULT_DATA_VALUE_LONG;
        }
      } else if (strcmp(chunk_ind->data_type_str, "SHORT") == 0) {
        for (j = 0; j < size_to_read_s; j++) {
          ((PIOSHORT *)data)[cur_pos+j] = DEFAULT_DATA_VALUE_SHORT;
        }
      } else if (strcmp(chunk_ind->data_type_str, "FLAG") == 0) {
        for (j = 0; j < size_to_read_s; j++) {
          ((PIOFLAG *)data)[cur_pos+j] = DEFAULT_DATA_VALUE_FLAG;
        }
      } else if (strcmp(chunk_ind->data_type_str, "BYTE") == 0) {
        for (j = 0; j < size_to_read_s; j++) {
          ((PIOBYTE *)data)[cur_pos+j] = DEFAULT_DATA_VALUE_BYTE;
        }
      } else if (strcmp(chunk_ind->data_type_str, "STRING") == 0) {
        for (j = 0; j < size_to_read_s; j++) {
          strncpy(((PIOSTRING *)data)[cur_pos+j], DEFAULT_DATA_VALUE_STRING, 255);
        }
      } else if (strcmp(chunk_ind->data_type_str, "COMPLE") == 0) {
        PIOCOMPLEX default_complex = {.real = DEFAULT_DATA_VALUE_FLOAT, .imaginary = DEFAULT_DATA_VALUE_FLOAT};
        for (j = 0; j < size_to_read_s; j++) {
          // Here we used default float value in a PIOCOMPLEX structure
          ((PIOCOMPLEX *)data)[cur_pos+j] = default_complex;
        }
      } else {
        fprintf(stderr, "ERROR: data format '%s' NOT SUPPORTED in %s()!\n", chunk_ind->data_type_str, __FUNCTION__);
        return -1;
      }
    } else { // Case of existing chunk
      chunk_index current_chunk_ind = chunk_ind->chunk_index_list[chunkIndexInIndexList];
      debug_print("    Corresponding chunk num: "PIOLONG_FMT"\n", current_chunk_ind.chunk_num);

      // 1- Retrieve the full path to object file
      debug_print("Retrieve the fullpath for the DATA file ('%s')\n", DATA_FILE_SUFFIX);
      char data_fullpath[MAX_STR_LEN];
      if (getFullpathToFileFromBackendname(metadata_p->Backendname, DATA_FILE_SUFFIX, current_chunk_ind.filenum+1, data_fullpath) != 0) {
        fprintf(stderr, "ERROR: unable to retrieve fullpath to object file '%s'\n", DATA_FILE_SUFFIX);
        return -1;
      }

      debug_print("    Copy "PIOLONG_FMT" samples of size %d\n", size_to_read_s, chunk_ind->data_type_size);

      // 2- Open file and read appropriate data

      // Compute the offset to be used for reading file
      off_t in_file_offset = headerSize + (current_chunk_ind.chunk_ord * chunkLength_s * chunk_ind->data_type_size) + (start_s*chunk_ind->data_type_size);

      // Compute the number of byte to be read in file
      size_t nb_byte_to_read_from_file = chunk_ind->data_type_size * size_to_read_s; // Size in Bytes

      // Check max number of bytes to be read
      if (nb_byte_to_read_from_file > SSIZE_MAX) {
        fprintf(stderr, "ERROR: The requested size (in byte) to be read is over capacity! (max=%ld  request=%ld)\n", SSIZE_MAX, nb_byte_to_read_from_file);
        return -1;
      }

      // Open the file
      int fd = open(data_fullpath, O_RDONLY);
      if (fd < 0) {
        fprintf(stderr, "ERROR while trying to open '%s'\n", data_fullpath);
        perror("Error");
        return -1;
      }

      // Read file and store result to 'data'
      ssize_t res_read = pread64(fd, data+(cur_pos*chunk_ind->data_type_size), nb_byte_to_read_from_file, in_file_offset);

      // Check return of read depending on scenario
      if (res_read < 0) {
        perror("Error");
        close(fd);
        return -1;
      }
      if (res_read != nb_byte_to_read_from_file) {
        fprintf(stderr, "ERROR while trying to read '%s' (not all requested bytes have been read!)\n", data_fullpath);
        close(fd);
        return -1;
      }

      // Close file
      close(fd);
    }

    cur_pos += size_to_read_s; // Update position index for output data
  }
  debug_print("  Read END\n");

  return 0;
}


/*
 * Allow to extract the selected bit (at position 'i') from the Byte 'data'.
 */
char extract(char data, int i) {
  return (data >> i) & 1;
}


/*
 * Read the required chunks of flags according to user request and store result in 'flags' parameter.
 * Note that due to offset and request sample the first/last chunks might NOT be read entirely!
 * when a chunk does not exist then we fill data with a default value.
 */
int readChunksFlags(const metadata *metadata_p, const pio_chunk_indexes *chunk_ind, const ordered_chunk_list_to_be_read *chunk_list, PIOFLAG *flags, whatToRead what_to_read) {
  debug_print("* %s()\n", __FUNCTION__);
  PIOLONG headerSize  = metadata_p->Iooffset;

  PIOLONG chunkLength_s = metadata_p->Flagchunksize * 8; // Length of a chunk (in number of samples). IMPORTANT: each byte code for 8 flags.
  debug_print("  chunkLength_s = "PIOLONG_FMT"\n", chunkLength_s);

  debug_print("  Start reading chunk(s)...\n");
  PIOLONG cur_pos = 0; // Current position in output data (index = sample index)
  PIOLONG i;
  for (i = 0; i < chunk_list->indexes_in_chunk_index_list_size; i++) {
    PIOLONG chunkIndexInIndexList = chunk_list->indexes_in_chunk_index_list[i];
    debug_print("   >chunk index pos: "PIOLONG_FMT"\n", i);

    // 0- Determine 'start_s' and 'length_s' of data to be read. Values are in number of samples.
    PIOLONG start_s  = 0;
    PIOLONG length_s = chunkLength_s; // By default the whole chunk
    // Handle case of FIRST chunk
    if (i == 0) {
      start_s = chunk_list->first_chunk_offset;
    }
    // Handle case of LAST chunk
    if (i == (chunk_list->indexes_in_chunk_index_list_size - 1)) {
      length_s = (chunk_list->last_chunk_limit_offset+1); // +1 since index is included
    }

    PIOLONG size_to_read_s = length_s - start_s; // Size to read in number of samples
    debug_print("    start_s       : "PIOLONG_FMT"\n", start_s);
    debug_print("    length_s      : "PIOLONG_FMT"\n", length_s);
    debug_print("    size_to_read_s: "PIOLONG_FMT"\n", size_to_read_s);

    // We must handle the case of non existing data (when index is NON_EXISTING_CHUNK_INDEX)
    if (chunkIndexInIndexList == NON_EXISTING_CHUNK_INDEX) {
      PIOLONG j;
      // Fill with default value for required length
      debug_print("    Non existing chunk -> fill '"PIOLONG_FMT"' samples with default value.\n", size_to_read_s);

      for (j = 0; j < size_to_read_s; j++) {
        flags[cur_pos+j] = DEFAULT_FLAG_VALUE;
      }
      cur_pos += size_to_read_s;

    } else { // Case of existing chunk
      chunk_index current_chunk_ind = chunk_ind->chunk_index_list[chunkIndexInIndexList];
      debug_print("    Corresponding chunk num: "PIOLONG_FMT"\n", current_chunk_ind.chunk_num);

      // 1- Retrieve the full path to object file
      const char *suffix = (what_to_read == READ_FLAG_WRITTEN) ? FLAG_FILE_SUFFIX:DATA_FILE_SUFFIX;
      debug_print("Retrieve the fullpath for the DATA or FLAG file ('%s')\n", suffix);
      char fullpath[MAX_STR_LEN];
      if (getFullpathToFileFromBackendname(metadata_p->Backendname, suffix, current_chunk_ind.filenum+1, fullpath) != 0) {
        fprintf(stderr, "ERROR: unable to retrieve fullpath to object file '%s'\n", suffix);
        return -1;
      }

      // 2- Open file and read appropriate data
      // Open the file
      int fd = open(fullpath, O_RDONLY);
      if (fd < 0) {
        fprintf(stderr, "ERROR while trying to open '%s'\n", fullpath);
        perror("Error");
        return -1;
      }

      // Remember that one byte contains 8 flags (one per bit)!!!
      // So we read the whole chunk of bytes then we extract bits info from them

      // Start to read after the header at the required offset!
      off_t seek_pos = headerSize + (current_chunk_ind.chunk_ord * (chunkLength_s/8)) + (start_s/8);
      debug_print("    need to seek to = %ld\n", seek_pos);

      PIOBYTE mod_start_s = start_s % 8; // the offset in the first byte
      debug_print("    mod_start_s = %d\n", mod_start_s);

      PIOLONG nb_byte_to_read = ((size_to_read_s-1)/8 + 1);
      if ((mod_start_s+(size_to_read_s-1)%8) >= 8) { // Detect if required to read an extra byte
        nb_byte_to_read += 1;
      }
      debug_print("    Need to retrieve "PIOLONG_FMT" samples => corresponding to "PIOLONG_FMT" BYTE\n", size_to_read_s, nb_byte_to_read);

      // Check max number of bytes to be read
      if (nb_byte_to_read > SSIZE_MAX) {
        fprintf(stderr, "ERROR: The requested size (in byte) to be read is over capacity! (max=%ld  request=%ld)\n", SSIZE_MAX, nb_byte_to_read);
        close(fd);
        return -1;
      }

      // Allocate tmp memory
      char *tmpFlagBytes = malloc(nb_byte_to_read * sizeof(char));
      if (tmpFlagBytes == NULL) {
        perror("Error");
        close(fd);
        return -1;
      }

      // Read the whole bunch of required bytes
      ssize_t res_read = pread64(fd, tmpFlagBytes, nb_byte_to_read, seek_pos);

      // Check return of read depending on scenario
      if (res_read < 0) {
        perror("Error");
        close(fd);
        return -1;
      }
      if (res_read != nb_byte_to_read) {
        fprintf(stderr, "ERROR while trying to read '%s' (not all requested bytes have been read!)\n", fullpath);
        close(fd);
        return -1;
      }

      // Extract flag values from bytes
      // Note that some bits may be not relevant in the first/last BYTE!
      PIOLONG processStartByte = 0;
      PIOBYTE mod_stop_s = (mod_start_s + size_to_read_s - 1) % 8; // This index is INCLUDED!!!
      debug_print("    mod_stop_s = %d\n", mod_stop_s);

      PIOBYTE bpos;
      if (nb_byte_to_read == 1) { // Case of only one BYTE to read: so we must care of 'start bit' and 'end bit'
        // Case of first BYTE
        debug_print("    *** first byte = '%x'\n", tmpFlagBytes[0]);
        for (bpos = mod_start_s; bpos <= mod_stop_s; bpos++) {
          flags[cur_pos] = extract(tmpFlagBytes[0], bpos);
          debug_print("    *** flags[]     = %d\n", flags[cur_pos]);
          cur_pos += 1;
        }
      } else { // Case of at least 2 BYTES to proceed
        // Case of first BYTE
        for (bpos = mod_start_s; bpos < 8; bpos++) {
          flags[cur_pos] = extract(tmpFlagBytes[0], bpos);
          debug_print("    *** flags["PIOLONG_FMT"] (in first byte)  = %d\n", cur_pos, flags[cur_pos]);
          cur_pos += 1;
        }

        // Case of inner bytes
        for (processStartByte = 1; processStartByte < (nb_byte_to_read-1); processStartByte++) {
          for (bpos = 0; bpos < 8; bpos++) {
            flags[cur_pos] = extract(tmpFlagBytes[processStartByte], bpos);
            cur_pos += 1;
          }
        }

        // Case of last byte
        for (bpos = 0; bpos <= mod_stop_s; bpos++) {
          flags[cur_pos] = extract(tmpFlagBytes[nb_byte_to_read-1], bpos);
          cur_pos += 1;
        }
      }

      // Free tmp memory
      free(tmpFlagBytes);

      // Closing file
      close(fd);
    }

  }
  debug_print("  Read END\n");

  return 0;
}

/*===========================================================================*/
/*===========================================================================*/
/*===========================================================================*/


/*
 * Reads data from the column named <column_name> in a binary table written in a FITS file.
 * The data is read starting from <first_sample> and with length <count>, then it is written
 * in <data>.
 * If <column_name>=NULL, the data is read from the first column of the first HDU containing data.
 * <expected_data_type> is the data type that the user expects to find at the requested location,
 * and and error is returned if it doesn't corresponfd to the data type encountered at the
 * requested location.
 * If <expected_data_type>=NULL, there is no verification.
 *
 * Used in the function "noDMC_readObject_GEN" below.
 */
int read_data_in_fits(const char *filename, const char *column_name, PIOLONG first_sample, PIOLONG count, void *data, const char *expected_data_type)
{
  int status = 0;

  debug_print("* %s()\n", __FUNCTION__);

  // Open the FITS file for reading
  fitsfile *fits_ptr;
  fits_open_file(&fits_ptr, filename, READONLY, &status);
  ERRORCHECK( status != 0, "occured when calling 'fits_open_file'. Error status: %d", status);

  // Get the number of HDUs present in the file
  int nhdu;
  fits_get_num_hdus(fits_ptr, &nhdu, &status);
  ERRORCHECK( status != 0, "occured when calling 'fits_get_num_hdus'. Error status: %d", status);


  // Browse all the HDUs to find the requested data
  int hdu_type;
  LONGLONG offset_in_file;
  int i;
  for(i=2; i<nhdu+1; i++) {

    // Point to the HDU number <i> (with the Primary HDU numeroted 1)
    fits_movabs_hdu(fits_ptr, i, &hdu_type, &status);
    ERRORCHECK( status != 0, "occured when calling 'fits_movabs_hdu'. Error status: %d", status);
    // Skips non binary table HDUs
    if(hdu_type != 2) {
      fprintf(stderr, "WARNING [%s()]: HDU %d has been skipped because it is not a binary table\n", __FUNCTION__, i);
      break;
    }

    // Reads the number of columns and rows in the HDU <i>
    int ncols;
    long nrows;
    fits_get_num_cols(fits_ptr, &ncols, &status);
    ERRORCHECK( status != 0, "occured when calling 'fits_get_num_cols'. Error status: %d", status);
    fits_get_num_rows(fits_ptr, &nrows, &status);
    ERRORCHECK( status != 0, "occured when calling 'fits_get_num_rowsll'. Error status: %d", status);



    // Check the case when no column has been requested, if there is some data in the HDU
    if(column_name == NULL && nrows != 0 && ncols !=0) {

      // Read the parameters of the first column
      char ttype[81], tunit[81], dtype[81];
      LONGLONG dimension;
      fits_get_bcolparmsll(fits_ptr, 1, ttype, tunit, dtype, &dimension, 0, 0, 0, 0, &status);
      ERRORCHECK( status != 0, "occured when calling 'fits_get_bcolparmsll'. Error status: %d", status);
      // Determine the cfitsio typecode corresponding to the data type of the first column
      int typecode;
      fits_binary_tform(dtype, &typecode, 0, 0, &status);
      ERRORCHECK( status != 0, "occured when calling 'fits_binary_tform'. Error status: %d", status);

      // Search for a keyword 'first_sample' and read its value if it is found
      fits_read_key( fits_ptr, TLONGLONG, "first_sample", &offset_in_file, 0, &status);
      if(status == KEY_NO_EXIST) {
        status         = 0;
        offset_in_file = 0;
      }
      ERRORCHECK( status != 0, "occured when calling 'fits_read_key'. Error status: %d", status);

      // Check if the requested <first_sample> is coherent with the 'first_sample' keyword's value
      LONGLONG offset = first_sample - offset_in_file;
      ERRORCHECK( offset < 0, "'first_sample' shall be greater than the first_sample written in the .fits file");

      // Verify that the found data type corresponds to the expected one
      if(typecode == TFLOAT) {
        ERRORCHECK( (expected_data_type != NULL) && (strcmp(expected_data_type, "FLOAT") != 0), "The data present in the first readable column has not expected type '%s', but has type '%s'", expected_data_type, "FLOAT");
      }
      else if(typecode == TDOUBLE) {
        ERRORCHECK( (expected_data_type != NULL) && (strcmp(expected_data_type, "DOUBLE") != 0), "The data present in the first readable column has not expected type '%s', but has type '%s'", expected_data_type, "DOUBLE");
      }
      else if(typecode == TBYTE) {
        ERRORCHECK( (expected_data_type != NULL) && (strcmp(expected_data_type, "BYTE") != 0), "The data present in the first readable column has not expected type '%s', but has type '%s'", expected_data_type, "BYTE");
      }
      else if(typecode == TLONG) {
        ERRORCHECK( (expected_data_type != NULL) && (strcmp(expected_data_type, "INT") != 0), "The data present in the first readable column has not expected type '%s', but has type '%s'", expected_data_type, "INT");

// TODO: understand why we need to change 'typecode' from TLONG to TINT to correctly read the file
        typecode = TINT;

      }
      else if(typecode == TLONGLONG) {
        ERRORCHECK( (expected_data_type != NULL) && (strcmp(expected_data_type, "LONG") != 0), "The data present in the first readable column has not expected type '%s', but has type '%s'", expected_data_type, "LONG");
      }
      else {
        ERRORCHECK( 1, "The data present in the first readable column has fitsio typecode %d which is not yet supported", typecode);
      }

      // Read the requested data and write it in <data>
      LONGLONG frow  = offset/dimension+1;
      LONGLONG felem = offset%dimension+1;

      fits_read_col(fits_ptr, typecode, 1, frow, felem, count, 0, data, 0, &status);
      ERRORCHECK( status != 0, "occured when calling 'fits_read_col'. Error status: %d", status);

      // Close the FITS file
      fits_close_file(fits_ptr, &status);
        ERRORCHECK( status != 0, "occured when calling 'fits_close_file'. Error status: %d", status);

      return 0;
    }


    // When there is one, search for the requested column
    else if(column_name != NULL && ncols !=0) {

      // Read the parameters of each column, and check if the name (TTYPE) corresponds to
      // the requested column
      char ttype[81], tunit[81], dtype[81];
      LONGLONG dimension;
      int j;
      for(j=1; j < ncols+1; j++) {
        fits_get_bcolparmsll(fits_ptr, j, ttype, tunit, dtype, &dimension, 0, 0, 0, 0, &status);
        ERRORCHECK( status != 0, "occured when calling 'fits_get_bcolparmsll'. Error status: %d", status);
        if(strcmp(column_name, ttype) == 0) {
          break;
          }
      }
      // Read the data if the requested column was found
      if(j != ncols+1) {
        // Check if the requested column holds data
        ERRORCHECK( nrows == 0, "requested column was found in HDU %d but is empty (nrows=0)", i);

        // Determine the cfitsio typecode corresponding to the data type of the first column
        int typecode;
        fits_binary_tform(dtype, &typecode, 0, 0, &status);
        ERRORCHECK( status != 0, "occured when calling 'fits_binary_tform'. Error status: %d", status);

        // Search for a keyword 'first_sample' and read its value if it is found
        fits_read_key( fits_ptr, TLONGLONG, "first_sample", &offset_in_file, 0, &status);
        if(status == KEY_NO_EXIST) {
            status         = 0;
            offset_in_file = 0;
        }
        ERRORCHECK( status != 0, "occured when calling 'fits_read_key'. Error status: %d", status);

        // Check if the requested <first_sample> is coherent with the 'first_sample' keyword's value
        LONGLONG offset = first_sample - offset_in_file;
        ERRORCHECK( offset < 0, "'first_sample' shall be greater than the first_sample written in the .fits file");

        // Verify that the found data type corresponds to the expected one
        if(typecode == TFLOAT) {
            ERRORCHECK( (expected_data_type != NULL) && (strcmp(expected_data_type, "FLOAT") != 0), "The data present in the first readable column has not expected type '%s', but has type '%s'", expected_data_type, "FLOAT");
        }
        else if(typecode == TDOUBLE) {
            ERRORCHECK( (expected_data_type != NULL) && (strcmp(expected_data_type, "DOUBLE") != 0), "The data present in the first readable column has not expected type '%s', but has type '%s'", expected_data_type, "DOUBLE");
        }
        else if(typecode == TBYTE) {
            ERRORCHECK( (expected_data_type != NULL) && (strcmp(expected_data_type, "BYTE") != 0), "The data present in the first readable column has not expected type '%s', but has type '%s'", expected_data_type, "BYTE");
        }
        else if(typecode == TLONG) {
            ERRORCHECK( (expected_data_type != NULL) && (strcmp(expected_data_type, "INT") != 0), "The data present in the first readable column has not expected type '%s', but has type '%s'", expected_data_type, "INT");

// TODO: understand why we need to change 'typecode' from TLONG to TINT to correctly read the file
        typecode = TINT;

        }
        else if(typecode == TLONGLONG) {
            ERRORCHECK( (expected_data_type != NULL) && (strcmp(expected_data_type, "LONG") != 0), "The data present in the first readable column has not expected type '%s', but has type '%s'", expected_data_type, "LONG");
        }
        else {
            ERRORCHECK( 1, "The data present in the first readable column has fitsio typecode %d which is not yet supported", typecode);
        }

        // Read the requested data and write it in <data>
        LONGLONG frow  = offset/dimension+1;
        LONGLONG felem = offset%dimension+1;
        fits_read_col(fits_ptr, typecode, j, frow, felem, count, 0, data, 0, &status);
        ERRORCHECK( status != 0, "occured when calling 'fits_read_col'. Error status: %d", status);

          // Close the FITS file
          fits_close_file(fits_ptr, &status);
            ERRORCHECK( status != 0, "occured when calling 'fits_close_file'. Error status: %d", status);

        return 0;
      }
    }
  }

  // Returns an error if either no column was requested but no column holds any valid data,
  // or if a column was requested but wasn't found
  if(i == nhdu+1 && column_name == NULL) {
    fprintf(stderr, "ERROR [%s()]: there is no readable data in the provided FITS file\n", __FUNCTION__);
  }
  else if(i == nhdu+1 && column_name != NULL) {
    fprintf(stderr, "ERROR [%s()]: the requested column was not found in the provided FITS file\n", __FUNCTION__);
  }

  // Close the FITS file
  fits_close_file(fits_ptr, &status);
  ERRORCHECK( status != 0, "occured when calling 'fits_close_file'. Error status: %d", status);

  return -1;
}




/*
 * This function is able to dispatche reading accordingly depending of the datatype (data on x byte OR data on bit).
 */
int noDMC_readObject_GEN(const char *object_name, PIOLONG offset, PIOLONG nbsample, void *data, const char *current_function_data_type_label, PIOSHORT output_data_type_size, whatToRead what_to_read) {

/* SM - 21 jan 2020
This function should be completely rewritten as follows:
- if object_name doesn't exist, add ".fits" at the end,
- if object_name is a file:
  - if it ends with ".fits" or contains ".fits#", call "read_data_in_fits()" which
    must be updated to handle the possible #column_name at the end of the file name
  - else call "read_data_from_binary()" which must be written with the block starting with
    " else { // read a flat binary file"
- if object_name is a directory, call a function "read_data_from_nodmc()", to be written,
  which contains everything from: "// 1- Try to retrieve the corresponding no dmc metadata"
- else return a "file not found" error
*/

  int ret;
  debug_print("* %s()\n", __FUNCTION__);
  debug_print("  object_name = %s\n", object_name);
  debug_print("  offset      = "PIOLONG_FMT"\n", offset);
  debug_print("  nbsample    = "PIOLONG_FMT"\n", nbsample);
  debug_print("  output data type str = %s\n", current_function_data_type_label);
  debug_print("  output data type     = %d\n", output_data_type_size);
  debug_print("  what_to_read         = %s\n", what_to_read == READ_DATA ? "READ_DATA":"READ_FLAG");

  // 0 - Check parameters
  if (object_name == NULL || data == NULL) {
    fprintf(stderr, "ERROR [%s()]: 'object_name' and 'data' parameters shall not be null!\n", __FUNCTION__);
    return -1;
  }
  if (offset < 0) {
    fprintf(stderr, "ERROR [%s()]: 'offset' parameter shall be >= 0\n", __FUNCTION__);
    return -1;
  }
  if (nbsample < 1) {
    fprintf(stderr, "ERROR [%s()]: 'nbsample' parameter shall be >= 1\n", __FUNCTION__);
    return -1;
  }


  // 0 - check if object_name is a file or a directory
  // if it's a file, read it as a flat binary file or as a FITS file

  struct stat s;
  uint len_filename = strlen(object_name);
  char object_name_fits[len_filename+5+1];
  object_name_fits[len_filename+5] = '\0';
  strncpy(object_name_fits, object_name, len_filename);
  strncpy(&object_name_fits[len_filename], ".fits", 5);
  struct stat s_fits;
  if ((stat( object_name, &s) == -1) && (strstr(object_name, ".fits") == NULL)) {
    s.st_mode = -1; // used to avoid bugs
    if (stat( object_name_fits, &s_fits) == -1) {
      fprintf(stderr, "ERROR [%s()]: '%s' not found\n", __FUNCTION__, object_name);
      return -1;
    }
  }

// commented out by SM on 20 jan 2020 because it produces segfaults on fits files.

/*
  else if (strstr(object_name, ".fits") != NULL) {// if a column is asked using '#', check that the given FITS file exists
    char *column_name;
    uint len_column_name;
    column_name = strstr(object_name, ".fits#");
    if(column_name != NULL) {
      column_name += 6;
      len_column_name = 0;
    }
    else{
      len_column_name = strlen(column_name);
    }
    char filename[len_filename - len_column_name +1];
    strncpy(filename, object_name, len_filename-len_column_name);
    filename[len_filename] = '\0';
    if (stat( filename, &s) == -1) {//changes info in s
      fprintf(stderr, "ERROR [%s()]: '%s' not found\n", __FUNCTION__, object_name);
      return -1;
    }
  }
*/

  if (S_ISREG(s.st_mode)) {
    if(strstr(object_name, ".fits") != NULL) { // read a FITS file
      char *column_name;
      column_name = strstr(object_name, ".fits#");
      if(column_name != NULL) {
        column_name += 6;
        len_filename -= strlen(column_name) + 1;
      }
      char filename[len_filename +1];
      strncpy(filename, object_name, len_filename);
      filename[len_filename] = '\0';

      return read_data_in_fits(filename, column_name, offset, nbsample, data, current_function_data_type_label);
    }
    else { // read a flat binary file
      int bytespersample;
      if (strcmp(current_function_data_type_label, "FLOAT") == 0) {
        bytespersample = 4;
      } else if (strcmp(current_function_data_type_label, "DOUBLE") == 0) {
        bytespersample = 8;
      } else if (strcmp(current_function_data_type_label, "BYTE") == 0) {
        bytespersample = 1;
      } else if (strcmp(current_function_data_type_label, "INT") == 0) {
        bytespersample = 4;
      } else if (strcmp(current_function_data_type_label, "LONG") == 0) {
        bytespersample = 8;
      } else {
        fprintf(stderr, "ERROR [%s()]: %s not implemented yet for flat binary files, please do it.\n", __FUNCTION__, current_function_data_type_label);
          return -1;
      }

      int bin_filedesc = open( object_name, O_RDONLY);
      if (bin_filedesc < 0) {
        fprintf( stderr, "ERROR [%s()]: error while opening <%s>.\n", __FUNCTION__, object_name);
        return -1;
      }
      ssize_t read_ret = pread64( bin_filedesc, data, nbsample * bytespersample, offset * bytespersample);
      close( bin_filedesc);
      if (read_ret != nbsample * bytespersample) {
        fprintf( stderr, "ERROR [%s()]: read %ld byte instead of %ld from <%s>.\n", __FUNCTION__, read_ret, nbsample * bytespersample, object_name);
        return -1;
      }
      return 0;
    }
  }
  else if (S_ISREG(s_fits.st_mode)) {
    return read_data_in_fits(object_name_fits, NULL, offset, nbsample, data, current_function_data_type_label);
  }

  // 1- Try to retrieve the corresponding no dmc metadata
  metadata meta_obj; // the metadata used to store info

  // Retrieve the "no dmc metadata"
  ret = getMetadataFor(object_name, &meta_obj);
  if (ret != 0) {
    return ret;
  }

  // 2- Retrieving info from chunk index file (<object>.0.pio)
  // The *.pio and *.pio.flg files are located in meta_obj.Backendname directory
  debug_print("  meta_obj.Backendname = %s\n", meta_obj.Backendname);
  pio_chunk_indexes current_pio_chunk_indexes;
  // Init the struct
  initPioChunkIndexes(&current_pio_chunk_indexes);
  // Retrieve chunk info
  ret = getChunkIndexInfo(&meta_obj, &current_pio_chunk_indexes);
  if (ret != 0) {
    freeMetadata(&meta_obj);
    return ret;
  }

  // Make some additional assertion about type of data to be read and current read function (must be exactly the same).
  // These test only apply when we try to read data not flags!
  if (what_to_read == READ_DATA) {
    if ((output_data_type_size != current_pio_chunk_indexes.data_type_size) ||
        // This second test on data type label is even more strict
        // Note that we restrict comparison up to 6 char since 'data_type_str' obetain from pio file is 6 char length max! For exemple a COMPLEX is labeled "COMPLE"
        (strncmp(current_function_data_type_label, current_pio_chunk_indexes.data_type_str, 6) != 0)) {
      fprintf(stderr, "\n*************\n**  ERROR  ** You are using noDMC_readObject_PIO%s() function to read another type of data (%s)!\n*************\n\n", current_function_data_type_label, current_pio_chunk_indexes.data_type_str);
      freePioChunkIndexes(&current_pio_chunk_indexes);
      freeMetadata(&meta_obj);
      return -1;
    }
  }

  // 3- Using the info from chunk index file to read appropriate data and copy it to user var 'data'

  // 3a- Determine the appropriate chunk to be read
  ordered_chunk_list_to_be_read current_ordered_chunk_list_to_be_read; // Structure to store what chunk indexes to read
initOrderedChunkListToBeRead(&current_ordered_chunk_list_to_be_read);
  ret = getOrderedChunkListToBeRead(&meta_obj, offset, nbsample, current_pio_chunk_indexes, &current_ordered_chunk_list_to_be_read);
  if (ret != 0) {
    freePioChunkIndexes(&current_pio_chunk_indexes);
    freeMetadata(&meta_obj);
    return ret;
  }

  //DEBUG!
  if (DEBUG) {
    debug_print("  ***************************\n");
    debug_print("  Which chunks are to be read (in order) [total: "PIOLONG_FMT" chunks]:\n", current_ordered_chunk_list_to_be_read.indexes_in_chunk_index_list_size);
    int i;
    for (i = 0; i < current_ordered_chunk_list_to_be_read.indexes_in_chunk_index_list_size; i++) {
      debug_print("  "PIOLONG_FMT"", current_ordered_chunk_list_to_be_read.indexes_in_chunk_index_list[i]);
      if (i == 0) {
        debug_printNP(" (offset: "PIOLONG_FMT")", current_ordered_chunk_list_to_be_read.first_chunk_offset);
      }
      if (i == (current_ordered_chunk_list_to_be_read.indexes_in_chunk_index_list_size - 1)) {
        debug_printNP(" (limit offset: "PIOLONG_FMT")", current_ordered_chunk_list_to_be_read.last_chunk_limit_offset);
      }
      debug_printNP("\n");
    }
    debug_print("  ***************************\n");
  }

  // 3b- Read chunks and copy data.
  if (strcmp(current_function_data_type_label, "FLAG") == 0) { // Case of data on bit (8 per byte)!!!!
    ret = readChunksFlags(&meta_obj, &current_pio_chunk_indexes, &current_ordered_chunk_list_to_be_read, data, what_to_read);
  } else { // All other cases (here we are sure that we need to read data)
    ret = readChunks(&meta_obj, &current_pio_chunk_indexes, &current_ordered_chunk_list_to_be_read, data);
  }

  // Free all allocated memory (even in case of error)
  freePioChunkIndexes(&current_pio_chunk_indexes);
freeOrderedChunkListToBeRead(&current_ordered_chunk_list_to_be_read);
  freeMetadata(&meta_obj);

  return ret;
}


/*
 *
 */
int noDMC_readObject_PIOBYTE(const char *object_name, PIOLONG offset, PIOLONG nbsample, PIOBYTE *data) {
  debug_print("* %s()\n", __FUNCTION__);

  // Custom def
  char *current_function_data_type_label = "BYTE";
  PIOSHORT output_data_type_size = sizeof(PIOBYTE);

  return noDMC_readObject_GEN(object_name, offset, nbsample, data, current_function_data_type_label, output_data_type_size, READ_DATA);
}


/*
 *
 */
int noDMC_readObject_PIOFLAG(const char *object_name, PIOLONG offset, PIOLONG nbsample, PIOFLAG *data) {
  debug_print("* %s()\n", __FUNCTION__);

  // Custom def
  char *current_function_data_type_label = "FLAG";
  PIOSHORT output_data_type_size = sizeof(PIOFLAG);

  return noDMC_readObject_GEN(object_name, offset, nbsample, data, current_function_data_type_label, output_data_type_size, READ_DATA);
}


/*
 *
 */
int noDMC_readObject_PIOSHORT(const char *object_name, PIOLONG offset, PIOLONG nbsample, PIOSHORT *data) {
  debug_print("* %s()\n", __FUNCTION__);

  // Custom def
  char *current_function_data_type_label = "SHORT";
  PIOSHORT output_data_type_size = sizeof(PIOSHORT);

  return noDMC_readObject_GEN(object_name, offset, nbsample, data, current_function_data_type_label, output_data_type_size, READ_DATA);
}


/*
 *
 */
int noDMC_readObject_PIOINT(const char *object_name, PIOLONG offset, PIOLONG nbsample, PIOINT *data) {
  debug_print("* %s()\n", __FUNCTION__);

  // Custom def
  char *current_function_data_type_label = "INT";
  PIOSHORT output_data_type_size = sizeof(PIOINT);

  return noDMC_readObject_GEN(object_name, offset, nbsample, data, current_function_data_type_label, output_data_type_size, READ_DATA);
}


/*
 *
 */
int noDMC_readObject_PIOLONG(const char *object_name, PIOLONG offset, PIOLONG nbsample, PIOLONG *data) {
  debug_print("* %s()\n", __FUNCTION__);

  // Custom def
  char *current_function_data_type_label = "LONG";
  PIOSHORT output_data_type_size = sizeof(PIOLONG);

  return noDMC_readObject_GEN(object_name, offset, nbsample, data, current_function_data_type_label, output_data_type_size, READ_DATA);
}


/*
 *
 */
int noDMC_readObject_PIOFLOAT(const char *object_name, PIOLONG offset, PIOLONG nbsample, PIOFLOAT *data) {
  debug_print("* %s()\n", __FUNCTION__);

  // Custom def
  char *current_function_data_type_label = "FLOAT";
  PIOSHORT output_data_type_size = sizeof(PIOFLOAT);

  return noDMC_readObject_GEN(object_name, offset, nbsample, data, current_function_data_type_label, output_data_type_size, READ_DATA);
}


/*
 *
 */
int noDMC_readObject_PIODOUBLE(const char *object_name, PIOLONG offset, PIOLONG nbsample, PIODOUBLE *data) {
  debug_print("* %s()\n", __FUNCTION__);

  // Custom def
  char *current_function_data_type_label = "DOUBLE";
  PIOSHORT output_data_type_size = sizeof(PIODOUBLE);

  return noDMC_readObject_GEN(object_name, offset, nbsample, data, current_function_data_type_label, output_data_type_size, READ_DATA);
}


/*
 *
 */
int noDMC_readObject_PIOSTRING(const char *object_name, PIOLONG offset, PIOLONG nbsample, PIOSTRING *data) {
  debug_print("* %s()\n", __FUNCTION__);

  // Custom def
  char *current_function_data_type_label = "STRING";
  PIOSHORT output_data_type_size = sizeof(PIOSTRING);

  return noDMC_readObject_GEN(object_name, offset, nbsample, data, current_function_data_type_label, output_data_type_size, READ_DATA);
}


/*
 *
 */
int noDMC_readObject_PIOCOMPLEX(const char *object_name, PIOLONG offset, PIOLONG nbsample, PIOCOMPLEX *data) {
  debug_print("* %s()\n", __FUNCTION__);

  // Custom def
  char *current_function_data_type_label = "COMPLEX";
  PIOSHORT output_data_type_size = sizeof(PIOCOMPLEX);

  return noDMC_readObject_GEN(object_name, offset, nbsample, data, current_function_data_type_label, output_data_type_size, READ_DATA);
}


/*
 *
 */
int noDMC_readFlagWritten(const char *object_name, PIOLONG offset, PIOLONG nbsample, PIOFLAG *flags) {
  debug_print("* %s()\n", __FUNCTION__);

  // Custom def
  char *current_function_data_type_label = "FLAG";
  PIOSHORT output_data_type_size = sizeof(PIOFLAG);

  return noDMC_readObject_GEN(object_name, offset, nbsample, flags, current_function_data_type_label, output_data_type_size, READ_FLAG_WRITTEN);
}


/*
 * Helper function to read stim input signal as PIOFLOAT of PIOINT in one call
 */

int noDMC_readObject_PIOINTorPIOFLOAT(const char *object_name, PIOLONG offset, PIOLONG nbsample, PIOFLOAT *data) {

  debug_print("* %s()\n", __FUNCTION__);

  // check if object_name is a dir or a file
  struct stat s;
  if (stat( object_name, &s) == -1) {
    fprintf(stderr, "ERROR [%s()]: '%s' not found\n", __FUNCTION__, object_name);
    return -1;
  }

  // default data type when something fails...
  PIOSTRING datatype;
  strcpy( datatype, "PIOFLOAT");

  if (S_ISREG(s.st_mode)) {
    // it's a flat binary file
    if (strstr( object_name, "int32.bin") != NULL) {
      strcpy( datatype, "PIOINT");
    }
  }
  else {
    // it's a nodmclib object
    metadata metdat;
    if (getMetadataFor(object_name, &metdat) == 0) {
      strncpy( datatype, metdat.Datatype, PIOSTRINGMAXLEN);
    }
  freeMetadata( &metdat);
  }

  if (strcmp( datatype, "PIOFLOAT") == 0) {
    return noDMC_readObject_GEN(object_name, offset, nbsample, data, "FLOAT", sizeof( PIOFLOAT), READ_DATA);
  }

  if (strcmp( datatype, "PIOINT") == 0) {
    PIOINT *tmp = (PIOINT *) malloc( nbsample * sizeof( PIOINT));
    int ret = noDMC_readObject_GEN(object_name, offset, nbsample, tmp, "INT", sizeof( PIOINT), READ_DATA);
    assert( ret == 0);
    for (long i = 0; i < nbsample; i++) {
      data[i] = tmp[i];
    }
    free( tmp);
    return ret;
  }

  assert( "noDMC_readObject_PIOINTorPIOFLOAT() can only read PIOINT or PIOFLOAT.");
  return -1;

}


