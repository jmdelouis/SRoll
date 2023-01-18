/******************************************************************************
 * "no_dmc_metadata.c"
 * This file aims to group all functions used to read the "no dmc metadata".
 * 
 * author:  Christan Madsen
 * date:    2014-12-22 (initial)
 * version: STABLE
 *****************************************************************************/

// Used for strtok_r() definition
// Note that the first define is for glibc < 2.10 and second define is for glibc >= 2.10
// RMK: on magique 3 we have: glibc 2.4
#define _GNU_SOURCE 1
#define POSIX_C_SOURCE 200809L

#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <libgen.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>

#include "no_dmc_metadata.h"
#include "no_dmc_util.h"
#include "no_dmc_debug.h"


// CONSTANTS for metadata
static const char METADATA_FILENAME[] = "no_dmc_metadata.txt"; // The filename in which metadata are stored
//static const char METADATA_SUPPORTED_VERSION[] = "0.0.2";

#define MAX_SIZE_KEYWORD_LIST 10000


//------------------------------------------------------------------------------

void extractIndexAndRingInfo(char *str, PIOLONG *index_p, PIOINT *ring_p) {
  // For memo, the line can be as: "BeginIndex : 1371347725 (ring 240)" or "BeginIndex : 0"
  char *delimiter_ring = "("; // The char used to know if there is a ring info inside [Begin/End]Index line

  debug_print("full line = '%s'\n", str);

  char *index = strtok(str, delimiter_ring);
  char *ring  = strtok(NULL, "");

  // a) Extract index info
  debug_print("index = '%s'\n", index);
  *index_p = strtol(index, NULL, 10);

  // b) Detect if there is ring info
  if (ring != NULL) { 
    // Extract ring info
    debug_print("We find a ring info = '%s'\n", ring);
    ring += strlen("ring"); // remove "ring" word
    ring = replace(ring, ')', ' '); // remove last parenthesis
    debug_print("clean ring info = '%s'\n", ring);
    *ring_p = strtol(ring, NULL, 10);
  }
  else {
    // Set default value
    debug_print("No ring info!\n");
    *ring_p = -1;
  }
}

//------------------------------------------------------------------------------

int getFullpathToFileFromBackendname(const char *backendname, const char *file_suffix, PIOLONG file_num, char *fullpath) {
  // Retrieve the object filename (=basename)
  char *tmpBackendname = strdupnodmc(backendname); // required for use in basename()
  if (tmpBackendname == NULL) {
    perror("Error");
    return -1;
  }
  char *obj_name = basename(tmpBackendname);
  debug_print("  obj_name = '%s'\n", obj_name);
  
  // Construct the index filename
  char filename[MAX_STR_LEN];
  snprintf(filename, MAX_STR_LEN-1, "%s."PIOLONG_FMT"%s", obj_name, file_num, file_suffix);
  debug_print("  filename = '%s'\n", filename);

  // Prepend the fullpath
  snprintf(fullpath, MAX_STR_LEN-1, "%s/%s", backendname, filename);
  debug_print("  fullpath = '%s'\n", fullpath);

  // Free mem
  if (tmpBackendname) {
    free(tmpBackendname);
    tmpBackendname = NULL;
  }

  return 0;
}

//------------------------------------------------------------------------------

int getMetadataFor(const char *dmc_obj, metadata *metadata_p) {
  char *key;
  char *value;
  char *delimiter = ":";

  char fullpath[MAX_STR_LEN];
  snprintf(fullpath, sizeof fullpath, "%s/%s", dmc_obj, METADATA_FILENAME); 
  //debug_print("fullpath = %s\n", fullpath);

  FILE *file = fopen(fullpath, "r" );
  if (file != NULL) {
    char line[MAX_STR_LEN];
    while (fgets(line, sizeof line, file) != NULL) { /* read a line */
      // Case of the 'Keyword' declaration line (no ':' delimiter)
      if (startsWith(line, "Kw=") == 0) {
        //debug_print("Find the keyword declaration!\n");
        // Concatenate the keyword on a newline
        strcat(metadata_p->Keyword_list, "\n");
        strcat(metadata_p->Keyword_list, trim(line));
        continue;
      }

      // 'key' will point to the part before the delimiter
      key = strtok(line, delimiter);
      key = trim(key);
      //debug_print("read key   = '%s'\n", key);

      // 'value' will point to the part after the delimiter
      value = strtok(NULL, "");
      value = trim(value);
      //debug_print("read value = '%s'\n", value);

      char *newStr = strdupnodmc(value);
      if (newStr == NULL) {
        perror("Error");
        fclose(file);
        return -1;
      }

      // Update metadata content
      if (strcmp(key, "noDMCmetadata_date") == 0) {
        metadata_p->noDMCmetadata_date = newStr;
      } else if (strcmp(key, "noDMCmetadata_version") == 0) {
        metadata_p->noDMCmetadata_version = newStr;
/*
        if (strcmp(metadata_p->noDMCmetadata_version, METADATA_SUPPORTED_VERSION) != 0) {
          fprintf(stderr, "ERROR: Unsupported metadata version ('%s') in file '%s'. Required version: '%s'\n",
              metadata_p->noDMCmetadata_version, fullpath, METADATA_SUPPORTED_VERSION);
          fclose(file);
          return -1;
        }
*/
      } else if (strcmp(key, "TOItype") == 0) {
        metadata_p->PIOtype = newStr;
      } else if (strcmp(key, "Datatype") == 0) {
        metadata_p->Datatype = newStr;
      } else if (strcmp(key, "BeginIndex") == 0) {
        extractIndexAndRingInfo(newStr, &metadata_p->BeginIndex, &metadata_p->BeginRing);
        free(newStr);
      } else if (strcmp(key, "EndIndex") == 0) {
        extractIndexAndRingInfo(newStr, &metadata_p->EndIndex, &metadata_p->EndRing);
        free(newStr);
      } else if (strcmp(key, "Author") == 0) {
        metadata_p->Author = newStr;
      } else if (strcmp(key, "Date") == 0) {
        metadata_p->Date = newStr;
      } else if (strcmp(key, "Keyword list") == 0) {
        metadata_p->Keyword_list = malloc(MAX_SIZE_KEYWORD_LIST * sizeof(char));
        if (metadata_p->Keyword_list == NULL) {
          perror("Error");
          free(newStr);
          fclose(file);
          return -1;
        }
        strncpy(metadata_p->Keyword_list, newStr, MAX_SIZE_KEYWORD_LIST);
        if (MAX_SIZE_KEYWORD_LIST > 0) { // Ensure the dest string contains as last char a null byte
          metadata_p->Keyword_list[MAX_SIZE_KEYWORD_LIST - 1]= '\0';
        }
        free(newStr);
      } else if (strcmp(key, "Location") == 0) {
        // Do nothing: this info is not relevant!
        free(newStr);
      } else if (strcmp(key, "Backendname") == 0) {
        debug_print("Backendname = '%s'\n", newStr);
        metadata_p->Backendname = newStr;
        //fprintf(stderr,"*** TEST: metadata_p->Backendname = '%s'\n", metadata_p->Backendname);
      } else if (strcmp(key, "Iooffset") == 0) {
        metadata_p->Iooffset = strtol(newStr, NULL, 10);
        free(newStr);
      } else if (strcmp(key, "Flagchunksize") == 0) {
        metadata_p->Flagchunksize = strtol(newStr, NULL, 10);
        free(newStr);
/*
      } else {
        fprintf(stderr,"WARNING: unknow key in metadata : '%s'\n", key);
        free(newStr);
*/
      }
    }
    fclose(file);
  } else {
    fprintf(stderr, "ERROR while trying to open METADATA file: '%s'\n", fullpath);
    perror("Error");
    return -1;
  }

  // Make some additional test to determine if the metadata may be outdated:
  // ie. in the case where the object has changed since the metadata generation
  // To check this, we compare the metadata generation date and the 0.pio file
  // date. The metadata date should be after 0.pio date, otherwise we print a
  // warning message to the user!
  char index_fullpath[MAX_STR_LEN];
  getFullpathToFileFromBackendname(metadata_p->Backendname, INDEX_FILE_SUFFIX, 0, index_fullpath);

  struct stat s_pio;
  struct stat s_met;
  if (stat(index_fullpath, &s_pio) != 0 || stat(fullpath, &s_met) != 0) {
    fprintf(stderr, "ERROR while trying to retrieve file modification date\n");
    perror("Error");
/*
  } else if (s_pio.st_mtime > s_met.st_mtime) {
    debug_print("Last file modification OBJECT  : %s", ctime(&s_pio.st_mtime));
    debug_print("Last file modification METADATA: %s", ctime(&s_met.st_mtime));
    fprintf(stderr, "!!! WARNING !!!: metadata file '%s' might contains outdated information since the object has changed since then.\n--> Please consider to re-generate metadata for this object!\n", fullpath);
*/
  }

  return 0;
}

//------------------------------------------------------------------------------

int freeMetadata(metadata *metadata_p) {
  if (metadata_p == NULL) {
    fprintf(stderr, "ERROR: unable to free metadata : NULL pointer!");
    return -1;
  }

  if (metadata_p->noDMCmetadata_date != NULL) {
    free(metadata_p->noDMCmetadata_date);
    metadata_p->noDMCmetadata_date = NULL;
  }
  if (metadata_p->noDMCmetadata_version != NULL) {
    free(metadata_p->noDMCmetadata_version);
    metadata_p->noDMCmetadata_version = NULL;
  }
  if (metadata_p->PIOtype != NULL) {
    free(metadata_p->PIOtype);
    metadata_p->PIOtype = NULL;
  }
  if (metadata_p->Datatype != NULL) {
    free(metadata_p->Datatype);
    metadata_p->Datatype = NULL;
  }
  if (metadata_p->Author != NULL) {
    free(metadata_p->Author);
    metadata_p->Author = NULL;
  }
  if (metadata_p->Date != NULL) {
    free(metadata_p->Date);
    metadata_p->Date = NULL;
  }
  if (metadata_p->Keyword_list != NULL) {
    free(metadata_p->Keyword_list);
    metadata_p->Keyword_list = NULL;
  }
  if (metadata_p->Backendname != NULL) {
    free(metadata_p->Backendname);
    metadata_p->Backendname = NULL;
  }

  return 0;
}

//------------------------------------------------------------------------------

int getMetadataKeywordFor(const char *dmc_obj, const char *keywordName, keyword *keyword_p) {
  metadata tmpmet;

  // 0- Init output keyword struct
  keyword_p->Kw = NULL;
  keyword_p->Val = NULL;
  keyword_p->Type = NULL;
  keyword_p->Com = NULL;

  // 1- First try to retrieve the full metadata
  int ret = getMetadataFor(dmc_obj, &tmpmet);
  if (ret != 0) {
    return ret;
  }

  // 2- Loop on each lines and try to extract the one that match: Kw='keywordName' (if any)
  // Need to copy string since strtok is destructive!
  char *tmpFullKeywords = strdupnodmc(tmpmet.Keyword_list);
  if (tmpFullKeywords == NULL) {
    perror("Error");
    freeMetadata(&tmpmet);
    return -1;
  }
  
  char *saveptr;
  char *line = strtok_r(tmpFullKeywords, "\n", &saveptr);
  while(line) {
    debug_print("line from keyword = \"%s\"\n", line);

    // Check that this is the line we want
    char tmpkey[MAX_STR_LEN];
    sprintf(tmpkey, "Kw='%s'", keywordName);
    if (startsWith(line, tmpkey) == 0) {
      int i;
      // Loop on all the 4 info composing a keyword
      char *tmpStr = line;
      char *saveptrInfo;
      for(i = 0; i < 4; i++) {
        tmpStr = strtok_r(tmpStr, "='", &saveptrInfo);
        debug_print("** tmpStr = \"%s\"\n", tmpStr);
        char *tmpContent = strtok_r(NULL, "'", &saveptrInfo); 
        debug_print("** tmpContent = \"%s\"\n", tmpContent);
        tmpStr = NULL;

        // Duplicate the string for storing in struct
        char *newStr = strdupnodmc(tmpContent);
        if (newStr == NULL) {
          perror("Error");
          return -1;
        }
        // Trim it (remove any space)
        newStr = trim(newStr);

        // Store it accordingly in keyword struct
        switch(i) {
          case 0:
            keyword_p->Kw = newStr;
            break;
          case 1:
            keyword_p->Val = newStr;
            break;
          case 2:
            keyword_p->Type = newStr;
            break;
          case 3:
            keyword_p->Com = newStr;
            break;
        }
      }

      // Since all done just stop the loop
      break;
    }

    // Otherwise, go for next line
    line = strtok_r(NULL, "\n", &saveptr);
  }

  // Check that we found the keyword, otherwise it is an error!
  if(keyword_p->Kw == NULL) {
    fprintf(stderr, "ERROR: There is no keyword named '%s' associated to object '%s'\n", keywordName, dmc_obj);
    ret = -1;
  }

  // Finally free memory
  freeMetadata(&tmpmet);

  return ret;
}

//------------------------------------------------------------------------------

int freeKeyword(keyword *keyword_p) {
  if (keyword_p == NULL) {
    fprintf(stderr, "ERROR: unable to free keyword structure: NULL pointer!");
    return -1;
  }

  if (keyword_p->Kw != NULL) {
    free(keyword_p->Kw);
    keyword_p->Kw = NULL;
  }
  if (keyword_p->Val != NULL) {
    free(keyword_p->Val);
    keyword_p->Val = NULL;
  }
  if (keyword_p->Type != NULL) {
    free(keyword_p->Type);
    keyword_p->Type = NULL;
  }
  if (keyword_p->Com != NULL) {
    free(keyword_p->Com);
    keyword_p->Com = NULL;
  }

  return 0;
}
