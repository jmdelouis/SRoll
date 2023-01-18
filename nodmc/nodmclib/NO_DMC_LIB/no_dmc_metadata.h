/******************************************************************************
 * "no_dmc_metadata.h"
 * 
 * author:  Christan Madsen
 * date:    2014-12-22
 * version: STABLE
 *****************************************************************************/

#ifndef _NO_DMC_METADATA_H_
#define _NO_DMC_METADATA_H_

#include "no_dmc_piolib_type_def.h"



// Other constants
static const char INDEX_FILE_SUFFIX[] = ".pio"; // The suffix of index file
static const char DATA_FILE_SUFFIX[]  = ".pio"; // The suffix of data file
static const char FLAG_FILE_SUFFIX[]  = ".pio.flg"; // The suffix of flag file



// Define the structure used for containing all the available metadata information
// that will be required in order be able to read directly the .pio files.
typedef struct {
  char *noDMCmetadata_date;
  char *noDMCmetadata_version;
  char *PIOtype; // Object type
  char *Datatype;
  PIOLONG BeginIndex; // Note: see 'BeginRing' for the optional (ring xx) info!
  PIOLONG EndIndex; // Note: see 'EndRing' for the optional (ring xx) info!
  PIOINT BeginRing; // This an optional value (ie. may not be set! -> in this case the value is -1)
  PIOINT EndRing; // This an optional value (ie. may not be set! -> in this case the value is -1)
  char *Author;
  char *Date;
  char *Keyword_list; // The string containing the keywords
  char *Backendname; // Where object data are stored
  PIOLONG Iooffset; // header size in byte of the .pio files
  PIOLONG Flagchunksize; // flag chunk size (in number of sample value) of .pio.flg files. Note that data chunk size is 8 times this one (in .pio file).
} metadata;

// Define the structure used for a keyword which can be present in metadata
// See function getMetadataKeywordFor()
typedef struct {
  char *Kw; // name (identifier)
  char *Val; // value (as string)
  char *Type; // type
  char *Com; // optional commentary associated to this keyword
} keyword;


/*
 * Allow to retrieve metadata corresponding to 'dmc_obj' into struct
 * 'metadata_p' pass in arg.
 * Note the caller is responsible to free memory after use by calling function
 * 'freeMetadata()'.
 * IMPORTANT: note that if the function detect an old metadata version, it
 * will stop immediatly and return an error code.
 * IMPORTANT2: note that if the function detect that the object has changed
 * since the no_dmc_metadata file generation (which may lead to wrong
 * information but do not prevent to read object content) we print a warning
 * message to user and process as usual!
 *
 * @param dmc_obj the fullpath to the dmc object name.
 * @param metadata_p the metadata structure in which info will be stored.
 *
 * @return 0 in case of success, non zero value otherwise.
 */
int getMetadataFor(const char *dmc_obj, metadata *metadata_p);


/*
 * Act exactly as getMetadataFor() function but instead of retrieving the full
 * metadata structure, it only retrieve (if available) the keyword value
 * associated to the keyword name pass in arg. Therefore this is only a subset
 * of the field named "Keyword_list".
 * Note the caller is responsible to free memory (in case of successfull
 * return) by calling freeKeyword().
 *
 * See getMetadataFor() for documentation.
 *
 * @param dmc_obj the fullpath to the dmc object name.
 * @param keywordName the keyword name to be retrieve.
 * @param keyword_p   a valid pointer to a keyword. If any value is
 *                    retrieve it will be copied in this structure.
 *                    Note: the current function is responsible for the
 *                    memory allocation of each element inside the struct.
 *
 * @return 0 in case of success, non zero value otherwise.
 */
int getMetadataKeywordFor(const char *dmc_obj, const char *keywordName, keyword *keyword_p);


/*
 * Free memory associated to the metadata structure content.
 *
 * @param metadata_p the metadata for which we want to free memory.
 *
 * @return 0 in case of success, non zero value otherwise.
 */
int freeMetadata(metadata *metadata_p);

/*
 * Free memory associated to the keyword structure content.
 *
 * @param keyword_p the keyword for which we want to free memory.
 *
 * @return 0 in case of success, non zero value otherwise.
 */
int freeKeyword(keyword *keyword_p);


int getFullpathToFileFromBackendname(const char *backendname, const char *file_suffix, PIOLONG file_num, char *fullpath);



#endif
