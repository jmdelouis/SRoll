/******************************************************************************
 * "test_no_dmc_lib.c"
 *
 * Unit/Functionnal testing of the "no dmc library".
 *
 * Note that this prog is not responsible to compare the read results regarding
 * the ones get with legacy DMC, see python script "compare_dmc.py" for this
 * purpose.
 *
 * Author: C. Madsen
 * Date: 2015-01-20
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include "no_dmc_metadata.h"
#include "no_dmc_data_access.h"
#include "no_dmc_piolib_type_def.h"
#include "no_dmc_util.h"
#include "no_dmc_debug.h"
#include "no_dmc_version.h"


// Size in byte of each reference types: In this test file, we used it only for
// control! (see check_piotype_length())
#define SIZEOF_PIOBYTE (1)
#define SIZEOF_PIOFLAG (1)
#define SIZEOF_PIOSHORT (2)
#define SIZEOF_PIOINT (4)
#define SIZEOF_PIOLONG (8)
#define SIZEOF_PIOFLOAT (4)
#define SIZEOF_PIODOUBLE (8)
#define SIZEOF_PIOSTRING (256)
#define SIZEOF_PIOCOMPLEX (8)

#define ID_PIOBYTE    "ID_PIOBYTE" 
#define ID_PIOFLAG    "ID_PIOFLAG"
#define ID_PIOSHORT   "ID_PIOSHORT"
#define ID_PIOINT     "ID_PIOINT"
#define ID_PIOLONG    "ID_PIOLONG"
#define ID_PIOFLOAT   "ID_PIOFLOAT"
#define ID_PIODOUBLE  "ID_PIODOUBLE"
#define ID_PIOSTRING  "ID_PIOSTRING"
#define ID_PIOCOMPLEX "ID_PIOCOMPLEX"


#define TEST_OK 0
#define TEST_ERROR -1

// Define Env config, so that input object can be accessed on different host target
// For now only M3 and M4 are supported
#define DATA_PATH_M3 "/data/dmc/MISS03/DATA/"
// This M4 path def rely on product on M3!!!
#define DATA_PATH_M4 "/m3gpfs3/datadmc/dmc/MISS03/DATA/"
//#define DATA_PATH_M4 "/pscratch1/delouis/SIM_data/dmc/MISS03/DATA/"


// This structure is used to group together all info required to process a test.
typedef struct {
  int      expectedResult; // Allow to handle error test: if we expecte an error result and we get it then the test is CORRECT! @see TEST_OK, TEST_ERROR
  char    *dmc_obj;    // full path to dmc object
  PIOLONG  offset;     // offset from which to read inside the object
  PIOLONG  nbSample;   // Number of sample to be read
  char    *datatype;   // data type id
  char    *keywords;   // keywords to be extracted (separated by a ',')
} test_scenario;


//-----------------------------------------------------------------------------

// Check methods from the no_dmc_util file.
int check_util() {
  int ret = 0;

  // BEGIN test for "replace()"

  char mystr[]      = "##this is##a examp#le";
  char mystr_REF1[] = "  this is  a examp le";
  char mystr_REF2[] = "__this_is__a_examp_le";

  debug_print("[TEST] mystr (input)  = '%s'\n", mystr);
  replace(mystr, '#', ' ');
  debug_print("[TEST] mystr (output) = '%s'\n", mystr);

  if (strcmp(mystr, mystr_REF1) != 0) {
    fprintf(stderr, "  Check replace(): FAILED!\n");
    return 1;
  }

  debug_print("[TEST] mystr (input)  = '%s'\n", mystr);
  replace(mystr, ' ', '_');
  debug_print("[TEST] mystr (output) = '%s'\n", mystr);

  if (strcmp(mystr, mystr_REF2) != 0) {
    fprintf(stderr, "  Check replace(): FAILED!\n");
    return 1;
  }

  fprintf(stderr, "  Check replace(): OK\n");

  // END test for "replace()"

  return ret;
}


//-----------------------------------------------------------------------------

// Check that all pio types has the expected size (which shall be constant on all platform!).
int check_piotype_length() {
  int ret = 0;

  // Check PIOBYTE
  if (sizeof(PIOBYTE) == SIZEOF_PIOBYTE) {
    fprintf(stderr, "  Check sizeof PIOBYTE (%d): OK\n", SIZEOF_PIOBYTE);
  } else {
    fprintf(stderr, "  Check sizeof PIOBYTE: ERROR!\n");
    ret = 1;
  }
  // Check PIOFLAG
  if (sizeof(PIOFLAG) == SIZEOF_PIOFLAG) {
    fprintf(stderr, "  Check sizeof PIOFLAG (%d): OK\n", SIZEOF_PIOFLAG);
  } else {
    fprintf(stderr, "  Check sizeof PIOFLAG: ERROR!\n");
    ret = 1;
  }
  // Check PIOSHORT
  if (sizeof(PIOSHORT) == SIZEOF_PIOSHORT) {
    fprintf(stderr, "  Check sizeof PIOSHORT (%d): OK\n", SIZEOF_PIOSHORT);
  } else {
    fprintf(stderr, "  Check sizeof PIOSHORT: ERROR!\n");
    ret = 1;
  }
  // Check PIOINT
  if (sizeof(PIOINT) == SIZEOF_PIOINT) {
    fprintf(stderr, "  Check sizeof PIOINT (%d): OK\n", SIZEOF_PIOINT);
  } else {
    fprintf(stderr, "  Check sizeof PIOINT: ERROR!\n");
    ret = 1;
  }
  // Check PIOLONG
  if (sizeof(PIOLONG) == SIZEOF_PIOLONG) {
    fprintf(stderr, "  Check sizeof PIOLONG (%d): OK\n", SIZEOF_PIOLONG);
  } else {
    fprintf(stderr, "  Check sizeof PIOLONG: ERROR!\n");
    ret = 1;
  }
  // Check PIOFLOAT
  if (sizeof(PIOFLOAT) == SIZEOF_PIOFLOAT) {
    fprintf(stderr, "  Check sizeof PIOFLOAT (%d): OK\n", SIZEOF_PIOFLOAT);
  } else {
    fprintf(stderr, "  Check sizeof PIOFLOAT: ERROR!\n");
    ret = 1;
  }
  // Check PIODOUBLE
  if (sizeof(PIODOUBLE) == SIZEOF_PIODOUBLE) {
    fprintf(stderr, "  Check sizeof PIODOUBLE (%d): OK\n", SIZEOF_PIODOUBLE);
  } else {
    fprintf(stderr, "  Check sizeof PIODOUBLE: ERROR!\n");
    ret = 1;
  }
  // Check PIOSTRING
  if (sizeof(PIOSTRING) == SIZEOF_PIOSTRING) {
    fprintf(stderr, "  Check sizeof PIOSTRING (%d): OK\n", SIZEOF_PIOSTRING);
  } else {
    fprintf(stderr, "  Check sizeof PIOSTRING: ERROR!\n");
    ret = 1;
  }
  // Check PIOCOMPLEX
  if (sizeof(PIOCOMPLEX) == SIZEOF_PIOCOMPLEX) {
    fprintf(stderr, "  Check sizeof PIOCOMPLEX (%d): OK\n", SIZEOF_PIOCOMPLEX);
  } else {
    fprintf(stderr, "  Check sizeof PIOCOMPLEX: ERROR!\n");
    ret = 1;
  }

  return ret;
}


//-----------------------------------------------------------------------------

// Check the "no dmc" version.
int check_nodmc_version() {
  int ret = 0;

  char* ver = noDMC_getVersion();
  if (ver != NULL) {
    fprintf(stderr, "  Check \"no dmc lib\" version (%s): OK\n", ver);
    free(ver);
  } else {
    fprintf(stderr, "  ERROR unable to retrieve the \"no dmc lib\" version!\n");
    ret = 1;
  }

  return ret;
}


//-----------------------------------------------------------------------------

// Allow to get size for the type specified in arg.
// This is used to automate the test using as input only a 'type' string.
size_t getSizeOfType(char *type) {
  int res = 0;
  size_t size = 0;

  if (strcmp(type, ID_PIOBYTE) == 0) {
    size = sizeof(PIOBYTE);
  } else if (strcmp(type, ID_PIOFLAG) == 0) {
    size = sizeof(PIOFLAG);
  } else if (strcmp(type, ID_PIOSHORT) == 0) {
    size = sizeof(PIOSHORT);
  } else if (strcmp(type, ID_PIOINT) == 0) {
    size = sizeof(PIOINT);
  } else if (strcmp(type, ID_PIOLONG) == 0) {
    size = sizeof(PIOLONG);
  } else if (strcmp(type, ID_PIOFLOAT) == 0) {
    size = sizeof(PIOFLOAT);
  } else if (strcmp(type, ID_PIODOUBLE) == 0) {
    size = sizeof(PIODOUBLE);
  } else if (strcmp(type, ID_PIOSTRING) == 0) {
    size = sizeof(PIOSTRING);
  } else if (strcmp(type, ID_PIOCOMPLEX) == 0) {
    size = sizeof(PIOCOMPLEX);
  } else {
    res = -1;
  }

  if (res != 0) {
    fprintf(stderr, "Unable to get size for type %s!\n", type);
    exit(res);
  }

  return size;
}



//-----------------------------------------------------------------------------

// Check 'res' value and if non zero print error message and exit, unless it was expected!
void check_res(int res, int expected_res) {
  bool success = (res == expected_res);

  fprintf(stderr, ">> Result = %d   (expected %d)   %s\n\n", res, expected_res, success?"SUCCESS":"FAILED");

  if (!success) {
    fprintf(stderr, "TEST FAILED detected in last call!\n");
    exit(res);
  }
}


//-----------------------------------------------------------------------------


// Dump read data (for debug purpose).
void dumpData(const void *data, const char *type, PIOLONG size) {

  fprintf(stderr, "-----------------------------------\n");
  fprintf(stderr, "DATA DUMP:\n");
  fprintf(stderr, "  * data type to be used for dump = %s\n", type);
  fprintf(stderr, "  * data size = "PIOLONG_FMT"\n", size);

  const PIOLONG MAX_SIZE_FOR_DUMP = 100;
  PIOLONG maxLoop = size;
  if (size > MAX_SIZE_FOR_DUMP) {
    maxLoop = MAX_SIZE_FOR_DUMP;
    fprintf(stderr, "  * [!] Too many data to dump: skip dump after "PIOLONG_FMT" elements\n", maxLoop);
  }

  PIOLONG i;
  if (strcmp(type, ID_PIOBYTE) == 0) {
    for (i = 0; i < maxLoop; i++) {
      // TODO: to be tested!
      fprintf(stderr, "  val #"PIOLONG_FMT" = \\x%02x\n", i, ((PIOBYTE *)data)[i]);
    }
  } else if (strcmp(type, ID_PIOFLAG) == 0) {
    for (i = 0; i < maxLoop; i++) {
      fprintf(stderr, "  flag #"PIOLONG_FMT" = '%s'\n", i, ((PIOFLAG *)data)[i] == 0 ? "False":"True");
    }
  } else if (strcmp(type, ID_PIOSHORT) == 0) {
    for (i = 0; i < maxLoop; i++) {
      // TODO: to be tested!
      fprintf(stderr, "  val #"PIOLONG_FMT" = %hi\n", i, ((PIOSHORT *)data)[i]);
    }
  } else if (strcmp(type, ID_PIOINT) == 0) {
    for (i = 0; i < maxLoop; i++) {
      fprintf(stderr, "  val #"PIOLONG_FMT" = "PIOINT_FMT"\n", i, ((PIOINT *)data)[i]);
    }
  } else if (strcmp(type, ID_PIOLONG) == 0) {
    for (i = 0; i < maxLoop; i++) {
      fprintf(stderr, "  val #"PIOLONG_FMT" = "PIOLONG_FMT"\n", i, ((PIOLONG *)data)[i]);
    }
  } else if (strcmp(type, ID_PIOFLOAT) == 0) {
    for (i = 0; i < maxLoop; i++) {
      // TODO: to be tested!
      fprintf(stderr, "  val #"PIOLONG_FMT" = %f\n", i, ((PIOFLOAT *)data)[i]);
    }
  } else if (strcmp(type, ID_PIODOUBLE) == 0) {
    for (i = 0; i < maxLoop; i++) {
      fprintf(stderr, "  val #"PIOLONG_FMT" = %g\n", i, ((PIODOUBLE *)data)[i]);
    }
  } else if (strcmp(type, ID_PIOSTRING) == 0) {
    for (i = 0; i < maxLoop; i++) {
      fprintf(stderr, "  val #"PIOLONG_FMT" = '%s'\n", i, ((PIOSTRING *)data)[i]);
    }
  } else if (strcmp(type, ID_PIOCOMPLEX) == 0) {
    for (i = 0; i < maxLoop; i++) {
      // TODO: to be tested!
      fprintf(stderr, "  val #"PIOLONG_FMT" = (%f,%f)\n", i,
          ((PIOCOMPLEX *)data)[i].real, ((PIOCOMPLEX *)data)[i].imaginary);
    }
  } else {
    fprintf(stderr, "INTERNAL TEST ERROR: Unsupported data type for dump! (%s)\n", type);
    exit(-1);
  }

  fprintf(stderr, "-----------------------------------\n");
}


//-----------------------------------------------------------------------------
// TEST "no_dmc_metadata.c"
// Test the metadata reading
//-----------------------------------------------------------------------------
int test_no_dmc_metadata(test_scenario ts) {
  metadata m_test; // the metadata used to store info

  // Call metadata read function for test obj
  int ret = getMetadataFor(ts.dmc_obj, &m_test);
  if (ret != 0) {
    return ret;
  }

  // Print retrieved info
  fprintf(stderr, "METADATA read for DMC object '%s'\n", ts.dmc_obj);
  fprintf(stderr, "  noDMCmetadata_date    = '%s'\n", m_test.noDMCmetadata_date);
  fprintf(stderr, "  noDMCmetadata_version = '%s'\n", m_test.noDMCmetadata_version);
  fprintf(stderr, "  PIOtype               = '%s'\n", m_test.PIOtype);
  fprintf(stderr, "  Datatype              = '%s'\n", m_test.Datatype);
  fprintf(stderr, "  BeginIndex            = "PIOLONG_FMT"\n", m_test.BeginIndex);
  fprintf(stderr, "  EndIndex              = "PIOLONG_FMT"\n", m_test.EndIndex);
  fprintf(stderr, "  BeginRing             = %d\n", m_test.BeginRing);
  fprintf(stderr, "  EndRing               = %d\n", m_test.EndRing);
  fprintf(stderr, "  Author                = '%s'\n", m_test.Author);
  fprintf(stderr, "  Date                  = '%s'\n", m_test.Date);
  fprintf(stderr, "  Keyword_list          = '%s'\n", m_test.Keyword_list);
  fprintf(stderr, "  Backendname           = '%s'\n", m_test.Backendname);
  fprintf(stderr, "  Iooffset              = "PIOLONG_FMT"\n", m_test.Iooffset);
  fprintf(stderr, "  Flagchunksize         = "PIOLONG_FMT"\n", m_test.Flagchunksize);

  // Free mem
  freeMetadata(&m_test);

  return 0;
}


//-----------------------------------------------------------------------------
// TEST "no_dmc_metadata.c"
// Test the metadata keyword reading
// This method will try to retrieve all required keywords values.
// Return 0 in case of all success, non zero if at least one failed.
//-----------------------------------------------------------------------------
int test_no_dmc_metadata_keyword(test_scenario ts) {
  int final_ret = 0; // The final result
  
  fprintf(stderr, "Will test for keywords {%s}:\n", ts.keywords);

  char *tmp; // Copy string: required for strtok!
  char *key; // Used to store the keyword name
  char *saveptr; // Required since we want to continue the strtok after having
                 // call some internal method that also use it!

  // Loop on all available keywords to be checked
  tmp = strdup(ts.keywords);
  char *memo = tmp;
  for(tmp = strdup(ts.keywords); ; tmp = NULL) {
    key = strtok_r(tmp, ",", &saveptr);
    if (key == NULL) { // When we reach the end (no more keyword)
      break;
    }

    // Remove spaces around retrieved string
    key = trim(key);
    //debug_print("[TEST] key  = '%s'\n", key);

    fprintf(stderr, "  Checking metadata keyword '%s':", key);

    // Try to retrieve the keyword info from metadata
    keyword k_test; // the keyword used to store info
    int ret = getMetadataKeywordFor(ts.dmc_obj, key, &k_test);
    if (ret != 0) {
      fprintf(stderr, "  --> FAILED!\n");
      final_ret = ret;
    } else {
      fprintf(stderr, "  --> OK\n");
      debug_print("[TEST] Kw   = '%s'\n", k_test.Kw);
      debug_print("[TEST] Val  = '%s'\n", k_test.Val);
      debug_print("[TEST] Type = '%s'\n", k_test.Type);
      debug_print("[TEST] Com  = '%s'\n", k_test.Com);

      // Free keyword struct
      freeKeyword(&k_test);
    }

  }
  free(memo);

  return final_ret;
}

//-----------------------------------------------------------------------------
// TEST "no_dmc_lib_access.c"
// Test capacity to read data using previously retrieve info.
//-----------------------------------------------------------------------------
int test_no_dmc_lib_access(test_scenario ts) {
  int ret;
  void *data;

  // Print debug info about current type to be used
  //fprintf(stderr, "[TEST] sizeof(PIODOUBLE) = %ld\n", sizeof(PIODOUBLE));


  // Allocate memory for data to be read
  data = calloc(ts.nbSample, getSizeOfType(ts.datatype));
  if (data == NULL) {
    fprintf(stderr, "Error: unable to allocate memory for data\n");
    return -1;
  }

  // Retrieve data using the new library
  // Select the appropriate function according to data type to be read
  if (strcmp(ts.datatype, ID_PIOBYTE) == 0) {
    ret = noDMC_readObject_PIOBYTE(ts.dmc_obj, ts.offset, ts.nbSample, data);
  } else if (strcmp(ts.datatype, ID_PIOFLAG) == 0) {
    ret = noDMC_readObject_PIOFLAG(ts.dmc_obj, ts.offset, ts.nbSample, data);
  } else if (strcmp(ts.datatype, ID_PIOSHORT) == 0) {
    ret = noDMC_readObject_PIOSHORT(ts.dmc_obj, ts.offset, ts.nbSample, data);
  } else if (strcmp(ts.datatype, ID_PIOINT) == 0) {
    ret = noDMC_readObject_PIOINT(ts.dmc_obj, ts.offset, ts.nbSample, data);
  } else if (strcmp(ts.datatype, ID_PIOLONG) == 0) {
    ret = noDMC_readObject_PIOLONG(ts.dmc_obj, ts.offset, ts.nbSample, data);
  } else if (strcmp(ts.datatype, ID_PIOFLOAT) == 0) {
    ret = noDMC_readObject_PIOFLOAT(ts.dmc_obj, ts.offset, ts.nbSample, data);
  } else if (strcmp(ts.datatype, ID_PIODOUBLE) == 0) {
    ret = noDMC_readObject_PIODOUBLE(ts.dmc_obj, ts.offset, ts.nbSample, data);
  } else if (strcmp(ts.datatype, ID_PIOSTRING) == 0) {
    ret = noDMC_readObject_PIOSTRING(ts.dmc_obj, ts.offset, ts.nbSample, data);
  } else if (strcmp(ts.datatype, ID_PIOCOMPLEX) == 0) {
    ret = noDMC_readObject_PIOCOMPLEX(ts.dmc_obj, ts.offset, ts.nbSample, data);
  } else {
    fprintf(stderr, "INTERNAL TEST ERROR: data type '%s' NOT supported!\n", ts.datatype);
    return -2;
  }

  // Check return
  if (ret != 0) {
    free(data);
    return ret;
  }
  
  // For debug purpose: dump retrieve data with the no_dmc_lib
  dumpData(data, ts.datatype, ts.nbSample);

  // free mem
  free(data);

  return 0;
}


//-----------------------------------------------------------------------------
// TEST "no_dmc_lib_readFlag.c"
// Test capacity to read flag using previously retrieve info.
//-----------------------------------------------------------------------------
int test_no_dmc_lib_readFlag(test_scenario ts) {
  int ret;
  PIOFLAG *flags;

  // Allocate memory for data to be read
  flags = calloc(ts.nbSample, sizeof(PIOFLAG));
  if (flags == NULL) {
    fprintf(stderr, "Error: unable to allocate memory for flags\n");
    return -1;
  }

  // Retrieve flags using the new library
  ret = noDMC_readFlagWritten(ts.dmc_obj, ts.offset, ts.nbSample, flags);

  // Check return
  if (ret != 0) {
    free(flags);
    return ret;
  }

  // For debug purpose: dump retrieve flags with the no_dmc_lib
  dumpData(flags, ID_PIOFLAG, ts.nbSample);

  // free mem
  free(flags);

  return 0;
}


//-----------------------------------------------------------------------------

// Allow to set appropriatly the define DATA_PATH
char *retrieveEnvDataPath() {
  char *val = getenv("HOSTNAME"); 
  char *ret = malloc(PIOSTRINGMAXLEN);
  ret[0] = '\0';
  
  //fprintf(stderr, "retrieveEnv() val = %s\n", val);

  // M3 config
  if (startsWith(val, "ln") == 0) {
    return strncat(ret, DATA_PATH_M3, PIOSTRINGMAXLEN);
  }
  // M4 config
  else if (startsWith(val, "log") == 0) {
    return strncat(ret, DATA_PATH_M4, PIOSTRINGMAXLEN);
  } else {
    fprintf(stderr, "ERROR: unknown env, can't proceed test. Only M3 and M4 are supported\n");
    exit(1);
  }
}


//*****************************************************************************
// MAIN
//*****************************************************************************
int main() {
  int res;

  char *dataPathBase = retrieveEnvDataPath();

  // Define the test scenario
  test_scenario test_scenario_list[] = {
    // Nominal cases: (= cases where we expect SUCCESS)
    {TEST_OK, "PSM190_FFP8_MAP/100-1a_fg_map_I", 0, 50331648, ID_PIODOUBLE, ""},
    {TEST_OK, "CL_MC/F8-WN-100-D1-JKYRxWN-100-D2-JKYR", 131071, 10, ID_PIODOUBLE, ""},
    {TEST_OK, "MAP_JMD_2048_PROD_MC/FULL_NOISE_100GHz_RD11_F0_v64_corr_survey_0_I", 0, 1, ID_PIODOUBLE, ""},
    {TEST_OK, "PSM190_FFP8_MAP/100-1a_fg_map_I", 2048, 2, ID_PIODOUBLE, ""},
    {TEST_OK, "DR2_MAP_GALACTIC_2048/mask_galaxy_0.80_apodised", 25600000, 20, ID_PIOFLOAT, ""},
// DOES NOT EXIST ANYMORE...    {TEST_OK, "e2eJan15_TOI/143-1a_00", 1371347725, 5, ID_PIOFLOAT, ""}, // First 5 values
// DOES NOT EXIST ANYMORE...    {TEST_OK, "e2eJan15_TOI/143-1a_00", 1371347720, 100, ID_PIOFLOAT, ""}, // 1 out then 4 inside
// DOES NOT EXIST ANYMORE...    {TEST_OK, "e2e_common_TOI/100-1a_skynops", 1348657332, 5000, ID_PIOFLOAT, ""}, // FLOAT - First real record (non zero)
// DOES NOT EXIST ANYMORE...    {TEST_OK, "e2e_common_TOI/100-1a_skynops", 1348599807, 5000, ID_PIOFLOAT, ""}, // FLOAT - before first chunk
    {TEST_OK, "ECC_v44/NAME", 910, 120, ID_PIOSTRING, ""}, // PIOSTRING (the only .pio files to have flagchunksize set to 16 instead of 16384)
//TOSEE for now do not support the chunk format!    {TEST_OK, "Sa_HFI_C_Bolo_ADCDec14/HFI_00_C", 1, 10, ID_PIOINT}, // This object is empty (begin=end=-1)! Therefore we shall not be able to read even a sample from it, and therefore the lib shall return only zeros... CAUTION: This object shall be of raw chunk type (the one with flag just after data)
//TOSEE for now do not support the chunk format!    {TEST_OK, "Sa_HFI_C_Bolo/HFI_00_C/"), 1, 10, ID_PIOINT}, // Same as previous but non empty!
    {TEST_OK, "RAW_4K/02_HARMS_LM4K_GP21", 700, 80327, ID_PIOCOMPLEX, ""},
    {TEST_OK, "RAW_4K/02_HARMS_LM4K_GP21", 81026, 10, ID_PIOCOMPLEX, "idx_harmonics, method, n_harmonics, Unit"}, // Read 1 valid value and 9 outside (ie. default)!

    // vvv BELOW vvv Error cases: (= cases where we expect an error) Note: the expected error is specify in comment
    //TOSEE ???? {TEST_ERROR, "PCCS_100_RC2/NAME", 3, 253, ID_PIOSTRING},
  };
  int test_scenario_list_size = 8;

  // Prepend path to data for each scenario
  int l;
  for (l = 0; l < test_scenario_list_size; l++) {
    char *str = malloc(PIOSTRINGMAXLEN);
    sprintf(str, "%s%s", dataPathBase, test_scenario_list[l].dmc_obj);
    test_scenario_list[l].dmc_obj = str;
  }


  fprintf(stderr, "***********************************************\n");
  fprintf(stderr, "** This is the test program for 'no dmc lib' **\n");
  fprintf(stderr, "***********************************************\n");
  fprintf(stderr, "\n\n");

  fprintf(stderr, "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
  fprintf(stderr, "$$ PRELIMINARY TESTS \n");
  fprintf(stderr, "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
  fprintf(stderr, "\n");

  // check_util()
  fprintf(stderr, "===============================================\n");
  fprintf(stderr, "Testing: check_util()\n");
  fprintf(stderr, "-----------------------------------------------\n");
  res = check_util();
  fprintf(stderr, "-----------------------------------------------\n");
  check_res(res, TEST_OK); // Always expect to succeed!

  // check_piotype_length()
  fprintf(stderr, "===============================================\n");
  fprintf(stderr, "Testing: check_piotype_length()\n");
  fprintf(stderr, "-----------------------------------------------\n");
  res = check_piotype_length();
  fprintf(stderr, "-----------------------------------------------\n");
  check_res(res, TEST_OK); // Always expect to succeed!

  // check_nodmc_version()
  fprintf(stderr, "===============================================\n");
  fprintf(stderr, "Testing: check_nodmc_version()\n");
  fprintf(stderr, "-----------------------------------------------\n");
  res = check_nodmc_version();
  fprintf(stderr, "-----------------------------------------------\n");
  check_res(res, TEST_OK); // Always expect to succeed!


  // Loop on all test scenario (unless force_test is set to a specific index)
  int force_test = -1; // -1 mean : loop on all test. Otherwise only test the specified test sample. IMPORTANT: test number is zero indexed! See 'test_scenario_list[]' for full list of test sample.
  int i;
  int start = 0;
  int stop  = test_scenario_list_size;
  if (force_test != -1) {
    start = force_test;
    stop  = force_test+1;
  }
  int nb_test_to_process = stop - start;
  
  for (i = start; i < stop; i++) {
    fprintf(stderr, "\n\n");
    fprintf(stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    fprintf(stderr, "!! SCENARIO #%d (%d/%d)\n", i, (i-start+1), nb_test_to_process);
    fprintf(stderr, "!!   dmc_obj  = %s\n", test_scenario_list[i].dmc_obj);
    fprintf(stderr, "!!   offset   = "PIOLONG_FMT"\n", test_scenario_list[i].offset);
    fprintf(stderr, "!!   nbSample = "PIOLONG_FMT"\n", test_scenario_list[i].nbSample);
    fprintf(stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n");

    // test_no_dmc_metadata()
    fprintf(stderr, "===============================================\n");
    fprintf(stderr, "Testing: test_no_dmc_metadata()\n");
    fprintf(stderr, "-----------------------------------------------\n");
    res = test_no_dmc_metadata(test_scenario_list[i]);
    fprintf(stderr, "-----------------------------------------------\n");
    check_res(res, TEST_OK); // Always expect to succeed!

    // Only test if something to do
    if (strcmp(test_scenario_list[i].keywords, "") != 0) {
      // test_no_dmc_metadata_keyword()
      fprintf(stderr, "===============================================\n");
      fprintf(stderr, "Testing: test_no_dmc_metadata_keyword()\n");
      fprintf(stderr, "-----------------------------------------------\n");
      res = test_no_dmc_metadata_keyword(test_scenario_list[i]);
      fprintf(stderr, "-----------------------------------------------\n");
      check_res(res, TEST_OK); // Always expect to succeed!
    }

    // test_no_dmc_lib_access()
    fprintf(stderr, "===============================================\n");
    fprintf(stderr, "Testing: test_no_dmc_lib_access()\n");
    fprintf(stderr, "-----------------------------------------------\n");
    res = test_no_dmc_lib_access(test_scenario_list[i]);
    fprintf(stderr, "-----------------------------------------------\n");
    check_res(res, test_scenario_list[i].expectedResult);

    // test_no_dmc_lib_readFlag()
    fprintf(stderr, "===============================================\n");
    fprintf(stderr, "Testing: test_no_dmc_lib_readFlag()\n");
    fprintf(stderr, "-----------------------------------------------\n");
    res = test_no_dmc_lib_readFlag(test_scenario_list[i]);
    fprintf(stderr, "-----------------------------------------------\n");
    check_res(res, test_scenario_list[i].expectedResult);
  }


  fprintf(stderr, "\n\n");
  fprintf(stderr, "vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv\n");
  fprintf(stderr, "Total: proceed %d test --> ALL SUCCESS !\n", nb_test_to_process);
  fprintf(stderr, "vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv\n");

  // Free memory
  for (l = 0; l < test_scenario_list_size; l++) {
    free(test_scenario_list[l].dmc_obj);
    test_scenario_list[l].dmc_obj = NULL;
  }
  free(dataPathBase);
  dataPathBase = NULL;

  return 0;
}
