/******************************************************************************
 * "testParLoader.c"
 *
 * This is the test program for the 'parLoader' facility.
 *
 * Author: C. Madsen
 * Date: 2015-04-03
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include "test_parLoader.h"

//-----------------------------------------------------------------------------

#define ID_PIOINT    "ID_PIOINT"
#define ID_PIOLONG   "ID_PIOLONG"
#define ID_PIOFLOAT  "ID_PIOFLOAT"
#define ID_PIODOUBLE "ID_PIODOUBLE"
#define ID_PIOSTRING "ID_PIOSTRING"

#define TEST_OK (0)
#define TEST_ERROR (1)

//-----------------------------------------------------------------------------

// Declare external function from the parLoader so that we can reproduce the
// same read behaviour.

extern PIODOUBLE myRead_PIODOUBLE(PIOSTRING value);
extern PIOLONG   myRead_PIOLONG(PIOSTRING value);
extern PIOFLOAT  myRead_PIOFLOAT(PIOSTRING value);
extern PIOINT    myRead_PIOINT(PIOSTRING value);

// Same as in parLoader.c
#define MAX_LINE_LENGTH 256

//-----------------------------------------------------------------------------

// This structure is used to group together all info required to process a UNIT test.
typedef struct {
  int  expectedResult;  // Allow to handle error test: if we expecte an error
                        // result and we get it then the test is CORRECT!
                        // @see TEST_OK, TEST_ERROR
  char *refStringValue; // the ref value string as it would be read in the parameter files...
  char *datatype;       // data type id
} test_unit;

//-----------------------------------------------------------------------------

// This structure is used to group together all info required to process a SCENARIO test.
typedef struct {
  int  expectedResult; // Allow to handle error test: if we expecte an error result and we get it then the test is CORRECT! @see TEST_OK, TEST_ERROR
  char *samplefile;    // Path to sample file (containing the parameters definition.
  char *description;   // Description of the sample file and/or expected result.
} test_scenario;

//-----------------------------------------------------------------------------

// Check 'res' value and if non zero print error message and exit, unless it was expected!
void check_res(int res, int expected_res, bool verbose) {
  bool success = (res == expected_res);

  if (verbose) {
    fprintf(stderr, ">> Result = %d   (expected %d)   %s\n\n", res, expected_res, success?"SUCCESS":"FAILED");
  } else {
    fprintf(stderr, " --> %s\n", success?"SUCCESS":"FAILED");
  }
    
  if (!success) {
    fprintf(stderr, "TEST FAILED detected in last call!\n");
    exit(res);
  }
}

//-----------------------------------------------------------------------------

/**
 * This function is responsible to put numbers representation (a string) in the
 * same shape so that they can be compared with each other!
 * In particular it process the string by:
 * - trim (remove space before/after)
 * - removing leading and trailing zeros ("00123.456000" -> "123.456")
 * - remove unnecessary dot ("12.0" -> "12"; "0." -> "0")
 */
char *removeLeadingTrailingZero(char *str) {
  char *res = (char *) malloc(sizeof(char) * MAX_LINE_LENGTH);
  int start = 0;
  int stop  = strlen(str)-1;

  //fprintf(stderr, "  str = '%s' (length=%zd)\n", str, strlen(str));

  // 0) Trim leading/trailing space
  while (*(str+start) == ' ' && strlen(str+start) > 0) {
    start++;
  }
  while (*(str+stop) == ' ' && stop >= start) {
    stop--;
  }
  
  //fprintf(stderr, "  after trim: str size = %d (%d - %d + 1)\n", stop-start+1, stop, start);
  if ((stop-start+1) == 0) {
    fprintf(stderr, "INTERNAL ERROR: Empty string!\n");
    exit(1);
  }

  // 1) Remove leading/trailing zeros

  // In all case we can remove leading zeros... (BEWARE of the value '0')
  while (*(str+start) == '0' && strlen(str+start) > 1) {
    start++;
  }

  //Find decimal (if any)
  char *p = strchr (str,'.');

  if (p != NULL) {
    //fprintf(stderr, "  decimal detected!\n");

    // Remove trailing zeros... (even last one after the dot)
    while (*(str+stop) == '0') {
      stop--;
    }
  }
  //fprintf(stderr, "  start = %d   stop = %d\n", start, stop);

  // 2) Remove optional dot '.'

  if (*(str+stop) == '.') {
    stop--;
  }

  // ** Make a copy in res of the result string
  if ((stop-start+1) == 0) {
    strncpy(res, "0", 1);
    res[1] = '\0';
  } else {
    // Return a copy
    strncpy(res, str+start, stop-start+1);
    // Set end of char
    res[stop-start+1] = '\0';

    // Avoid case of ".0" (transforme it to '0')
    if (strcmp(res, ".0") == 0) {
      res[0] = '0';
      res[1] = '\0';
    }
  }

  return res;
}

//-----------------------------------------------------------------------------

/**
 * This function compare the two string representation of number.
 * Note: it reshape number representations to be able to check.
 *
 * @return 0 if equal, non zero value otherwise.
 */
int compare(char *val, char *ref) {
  int res;
  //fprintf(stderr, "\n");

  char *val_new = removeLeadingTrailingZero(val);
  char *ref_new = removeLeadingTrailingZero(ref);

  //fprintf(stderr, "val = '%s' \t--> '%s'\n", val, val_new);
  //fprintf(stderr, "ref = '%s' \t--> '%s'\n", ref, ref_new);
  if (strncmp(val_new, ref_new, MAX_LINE_LENGTH) != 0) {
    //fprintf(stderr, "--> NOT THE SAME!!!\n");
    res = 1;
  } else {
    //fprintf(stderr, "--> OK\n");
    res = 0;
  }
  
  free(val_new);
  free(ref_new);

  return res;
}

//-----------------------------------------------------------------------------

/*
 * FOR DEBUG PURPOSE ONLY!!!
 * Allow to dump some content of Param.
 */
/*
void dumpParam(test_parContent *par) {
  fprintf(stderr, "---------------------------------\n");
  fprintf(stderr, "-- DUMP of Param\n");
  fprintf(stderr, "---------------------------------\n");

  char *param_name = NULL;
  int i;
  
  // Param of type PIOSTRING (list), optional
  param_name = "fsl";
  if (par->flag_fsl) {
    fprintf(stderr, "* \"%s\" (#"PIOLONG_FMT" elements):\n", param_name, par->n_fsl);
    for (i = 0; i < par->n_fsl; i++) {
      fprintf(stderr, "  %d:  %s\n", i+1, par->fsl[i]);
    }
    //fprintf(stderr, "\n");
  } else {
    fprintf(stderr, "* \"%s\" has not been defined!\n", param_name);
  }

  // Param of type PIOSTRING (list), optional
  param_name = "Signal_noPS";
  if (par->flag_Signal_noPS) {
    fprintf(stderr, "* \"%s\" (#"PIOLONG_FMT" elements):\n", param_name, par->n_Signal_noPS);
    for (i = 0; i < par->n_Signal_noPS; i++) {
      fprintf(stderr, "  %d:  %s\n", i+1, par->Signal_noPS[i]);
    }
    //fprintf(stderr, "\n");
  } else {
    fprintf(stderr, "* \"%s\" has not been defined!\n", param_name);
  }

  // Param of type PIOINT, mandatory
  param_name = "NADU";
  fprintf(stderr, "* \"%s\": "PIOINT_FMT"\n", param_name, par->NADU);

  // Param of type PIODOUBLE (list), mandatory
  param_name = "Calibration";
  fprintf(stderr, "* \"%s\" (#"PIOLONG_FMT" elements):\n", param_name, par->n_Calibration);
  for (i = 0; i < par->n_Calibration; i++) {
    fprintf(stderr, "  %d:  %.12g\n", i+1, par->Calibration[i]);
  }
  //fprintf(stderr, "\n");

  // Param of type PIOFLOAT (list), mandatory
  param_name = "FSLCOEF";
  fprintf(stderr, "* \"%s\" (#"PIOLONG_FMT" elements):\n", param_name, par->n_FSLCOEF);
  for (i = 0; i < par->n_FSLCOEF; i++) {
    fprintf(stderr, "  %d:  %f\n", i+1, par->FSLCOEF[i]);
  }
  //fprintf(stderr, "\n");

  // Param of type PIOLONG (list), mandatory
  param_name = "BeginRing";
  fprintf(stderr, "* \"%s\" (#"PIOLONG_FMT" elements):\n", param_name, par->n_BeginRing);
  for (i = 0; i < par->n_BeginRing; i++) {
    fprintf(stderr, "  %d:  "PIOLONG_FMT"\n", i+1, par->BeginRing[i]);
  }
  //fprintf(stderr, "\n");




  // Other....


  fprintf(stderr, "---------------------------------\n");
}
*/

//-----------------------------------------------------------------------------


//*****************************************************************************
// MAIN
//*****************************************************************************
int main(int argc, char *argv[]) {
  int res;

  fprintf(stderr, "********************************************************\n");
  fprintf(stderr, "** This is the test program for 'parLoader' facility. **\n");
  fprintf(stderr, "********************************************************\n");


  //DEBUG:
  //fprintf(stderr, "res = %d\n", compare("", ""));
  //fprintf(stderr, "res = %d\n", compare("123", "123 "));
  //fprintf(stderr, "res = %d\n", compare("123", "0123 "));
  //return 0;


  //***************************************************************************
  // *** PHASE 1 *** (unit testing)
  //***************************************************************************

  fprintf(stderr, "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
  fprintf(stderr, "$$ UNIT TESTS \n");
  fprintf(stderr, "$$ Purpose of theses tests are to check that the conversion\n");
  fprintf(stderr, "$$ from string to value are valid.\n");
  fprintf(stderr, "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
  fprintf(stderr, "\n");

  // In this phase we loop over a set of unit test.
  // For each we use the loader conversion tool (string -> value) then we
  // compare the value (using a reverse process value -> string) with the
  // reference string.

  // Define the test unit list
  test_unit test_unit_list[] = {
    // Nominal cases: (= cases where we expect SUCCESS)

    // -- PIODOUBLE --
    {TEST_OK, "0", ID_PIODOUBLE},
    {TEST_OK, "0.", ID_PIODOUBLE},
    {TEST_OK, ".0", ID_PIODOUBLE},
    {TEST_OK, "1.", ID_PIODOUBLE},
    {TEST_OK, "1.0", ID_PIODOUBLE},
    {TEST_OK, "10.00", ID_PIODOUBLE},
    {TEST_OK, "10.20", ID_PIODOUBLE},
    {TEST_OK, ".0123", ID_PIODOUBLE}, // Value starting by a dot
    {TEST_OK, "0.0915", ID_PIODOUBLE},
    {TEST_OK, "0.09150", ID_PIODOUBLE}, // Check that leading/trailing zero is OK
    {TEST_OK, "0.0915000000", ID_PIODOUBLE}, // Check that leading/trailing zero is OK
    {TEST_OK, "1.78260000", ID_PIODOUBLE}, // Check trailing zeros
    {TEST_OK, "9.99999999999e-13", ID_PIODOUBLE}, // Test with many 9
    {TEST_OK, "1.88008633662e-13", ID_PIODOUBLE}, // Number on 12 digits (MAX digits!)
    {TEST_OK, "1.38569754004e-16", ID_PIODOUBLE}, // Number on 12 digits (MAX digits!)
    {TEST_OK, "6.44604950677e-15", ID_PIODOUBLE}, // Number on 12 digits (MAX digits!)
    {TEST_OK, "1.78260067e-13", ID_PIODOUBLE}, // Number on 9 digits

    // -- PIOLONG --
    {TEST_OK, "0", ID_PIOLONG},
    {TEST_OK, "-9223372036854775807", ID_PIOLONG}, // Min value!
    {TEST_OK, "9223372036854775807", ID_PIOLONG}, // Max value!
    {TEST_OK, "240", ID_PIOLONG},

    // -- PIOFLOAT --
    {TEST_OK, "0", ID_PIOFLOAT},
    {TEST_OK, "0.", ID_PIOFLOAT},
    {TEST_OK, ".0", ID_PIOFLOAT},
    {TEST_OK, "1.", ID_PIOFLOAT},
    {TEST_OK, "1.0", ID_PIOFLOAT},
    {TEST_OK, "10.00", ID_PIOFLOAT},
    {TEST_OK, "10.20", ID_PIOFLOAT},
    {TEST_OK, ".0123", ID_PIOFLOAT}, // Value starting by a dot
    {TEST_OK, "0.0915", ID_PIOFLOAT},
    {TEST_OK, "0.09150", ID_PIOFLOAT}, // Check that leading/trailing zero is OK

    // -- PIOINT --
    {TEST_OK, "0", ID_PIOINT},
    {TEST_OK, "1234", ID_PIOINT},
    {TEST_OK, "-1", ID_PIOINT},
    {TEST_OK, "-10", ID_PIOINT},
    {TEST_OK, "-2147483647", ID_PIOINT}, // Min value!
    {TEST_OK, "2147483647", ID_PIOINT}, // Max value!

  };

  int test_unit_list_size = 37;
 
  // Loop on each test... valid or invalid...
  for (int i = 0; i < test_unit_list_size; i++) {
    test_unit tu = test_unit_list[i];
    char str[MAX_LINE_LENGTH];

    char tmp[MAX_LINE_LENGTH+2];
    strcpy(tmp, "\"");
    strcat(tmp, tu.refStringValue);
    strcat(tmp, "\"");
    fprintf(stderr, "Testing TU#%02d: %22s as %13s", i+1, tmp, tu.datatype);

    if (strcmp(tu.datatype, ID_PIODOUBLE) == 0) {
      // 1) Read value from ref string using parloader function
      PIODOUBLE tmp_piodouble = myRead_PIODOUBLE(tu.refStringValue);
      // 2) Convert it back to string to compare it
      snprintf(str, MAX_LINE_LENGTH, "%.12g", tmp_piodouble);
    } else if (strcmp(tu.datatype, ID_PIOLONG) == 0) {
      PIOLONG tmp_piolong = myRead_PIOLONG(tu.refStringValue);
      snprintf(str, MAX_LINE_LENGTH, PIOLONG_FMT, tmp_piolong);
    } else if (strcmp(tu.datatype, ID_PIOFLOAT) == 0) {
      PIOFLOAT tmp_piofloat = myRead_PIOFLOAT(tu.refStringValue);
      snprintf(str, MAX_LINE_LENGTH, "%f", tmp_piofloat);
    } else if (strcmp(tu.datatype, ID_PIOINT) == 0) {
      PIOINT tmp_pioint = myRead_PIOINT(tu.refStringValue);
      snprintf(str, MAX_LINE_LENGTH, PIOINT_FMT, tmp_pioint);
    } else {
      fprintf(stderr, "INTERNAL ERROR! (unknown datatype)\n");
      return 1;
    }

    // 3) Finally compare the two string representation
    res = compare(str, tu.refStringValue);

    // DEBUG
    /*
    if (res != 0) {
      fprintf(stderr, "\n[DBG] ref='%s' str='%s'\n", tu.refStringValue, str);
    }
    */

    check_res(res, tu.expectedResult, false); // Check according test unit expected result!
  }



  //***************************************************************************
  // *** PHASE 2 *** (real case testing)
  //***************************************************************************

  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
  fprintf(stderr, "$$ REAL CASES TESTS \n");
  fprintf(stderr, "$$ Purpose of theses tests are to check that the loader\n");
  fprintf(stderr, "$$ is operationnal in real case scenario.\n");
  fprintf(stderr, "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
  fprintf(stderr, "\n");

  // In this phase we loop over a real set of parameter files in order to check
  // that we are able to read them or able to detect their errors.

  // Define the test scenario
  test_scenario test_scenario_list[] = {
    // Nominal cases: (= cases where we expect SUCCESS)
    { TEST_OK, "test/samples/Param_OK_100_RD12.txt", "" },
    { TEST_OK, "test/samples/Param_OK_143_RD12.txt", "" },
    { TEST_OK, "test/samples/Param_OK_217_RD12.txt", "" },
    { TEST_OK, "test/samples/Param_OK_353_RD12.txt", "" },
    { TEST_OK, "test/samples/Param_OK_100_RD12_SIM_64.txt", "" },
    { TEST_OK, "test/samples/Param_OK_but_poorly_written.txt", "Well this file is poorly written BUT shall succeed!" },
    { TEST_OK, "test/samples/Param_OK_line_very_long.txt", "Ensure that 253 chars long line is ok" },
    { TEST_OK, "test/samples/Param_OK_numberOf_set_to_zero_for_optional_param.txt", "'number_of' can be set to 0 for optional parameter." },

    // Error cases: (= cases where we expect ERROR)
    { TEST_ERROR, "test/samples/Param_INVALID_missing_elt_01.txt", "Missing some elements in a poorly written file." },
    { TEST_ERROR, "test/samples/Param_INVALID_missing_value_01.txt", "Missing value for an element in a list" },
    { TEST_ERROR, "test/samples/Param_INVALID_missing_value_02.txt", "Missing value for 'number_of' parameter (empty spaces)" },
    { TEST_ERROR, "test/samples/Param_INVALID_missing_value_03.txt", "Missing value for 'number_of' parameter" },
    { TEST_ERROR, "test/samples/Param_INVALID_missing_value_04.txt", "Invalid value." },
    { TEST_ERROR, "test/samples/Param_INVALID_line_too_long_01.txt", "Line of comment is too long." },
    { TEST_ERROR, "test/samples/Param_INVALID_line_too_long_02.txt", "Line of param is too long." },
    { TEST_ERROR, "test/samples/Param_INVALID_numberOf_zero_but_elements.txt", "If 'number_of' set to 0 for optional parameter we expect no parameter value!" },
    { TEST_ERROR, "test/samples/Param_INVALID_number_of.txt", "Element 'number_of' is not allowed for element which are not defined as a list!" },
  };

  int test_scenario_list_size = 17;

  //char *filename = "samples/TMP_DEBUG_JMD/Param_143_RD12_SIM_512.txt";

  // Loop on each sample... valid or invalid...
  for (int i = 0; i < test_scenario_list_size; i++) {
    test_scenario ts = test_scenario_list[i];
    char *filename = ts.samplefile;
    test_parContent par;
    
    fprintf(stderr, "===============================================\n");
    fprintf(stderr, "Test Case TC#%02d\n", i+1);
    fprintf(stderr, "Testing param loading from: \"%s\"\n", filename);
    fprintf(stderr, "Description: %s\n", ts.description);
    fprintf(stderr, "-----------------------------------------------\n");
    res = test_readParam(&par, filename);

    // For debug
    /*
    if (i == 8) {
      fprintf(stderr, "\n");
      dumpParam(&par);
      fprintf(stderr, "\n");
    }
    */

    fprintf(stderr, "-----------------------------------------------\n");
    check_res(res, ts.expectedResult, true); // Check according test scenario expected result!
  }


  //***************************************************************************
  // *** PHASE 3 *** (specific case testing)
  //***************************************************************************
  char *filename = "test/samples/Param_TEST_FLAG.txt";
  test_parContent par;

  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
  fprintf(stderr, "$$ SPECIFIC CASES TESTS \n");
  fprintf(stderr, "$$ Purpose of theses tests are to check some specific situations.\n");
  fprintf(stderr, "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
  fprintf(stderr, "\n");

  fprintf(stderr, "===============================================\n");
  fprintf(stderr, "Test Specific Case TSC#01\n");
  fprintf(stderr, "Testing flag setting from: \"%s\"\n", filename);
  fprintf(stderr, "Description: This test allow to check that flag of optional parameters are correctly set in any situation (presence or absence of the optional parameter).\n");
  fprintf(stderr, "-----------------------------------------------\n");
  res = test_readParam(&par, filename);

  fprintf(stderr, "-----------------------------------------------\n");
  fprintf(stderr, "Loading param file");
  check_res(res, TEST_OK, false); // Check that file can be read

  // Now make the flag test:
  bool flag_test_success = (par.flag_AVGR0 == _PAR_FALSE) &&
                           (par.flag_ALMMAP == _PAR_FALSE) &&
                           (par.flag_VarGain == _PAR_FALSE) &&
                           (par.flag_Signal_noPS == _PAR_TRUE) &&
                           (par.flag_fsl == _PAR_TRUE) &&
                           (par.flag_verbose == _PAR_TRUE) &&
                           (par.flag_dmc_output_path == _PAR_TRUE);
  
  if (flag_test_success) {
    res = TEST_OK;
  } else {
    res = TEST_ERROR;
    fprintf(stderr, "DEBUG flags:\n");
    fprintf(stderr, "  > par.flag_AVGR0 = %d (expecting %d)\n", par.flag_AVGR0, _PAR_FALSE);
    fprintf(stderr, "  > par.flag_ALMMAP = %d (expecting %d)\n", par.flag_ALMMAP, _PAR_FALSE);
    fprintf(stderr, "  > par.flag_VarGain = %d (expecting %d)\n", par.flag_VarGain, _PAR_FALSE);
    fprintf(stderr, "  > par.flag_Signal_noPS = %d (expecting %d)\n", par.flag_Signal_noPS, _PAR_TRUE);
    fprintf(stderr, "  > par.flag_fsl = %d (expecting %d)\n", par.flag_fsl, _PAR_TRUE);
    fprintf(stderr, "  > par.flag_verbose = %d (expecting %d)\n", par.flag_verbose, _PAR_TRUE);
    fprintf(stderr, "  > par.flag_dmc_output_path = %d (expecting %d)\n", par.flag_dmc_output_path, _PAR_TRUE);
  }

  check_res(res, TEST_OK, true); // Check the whole flag test



  fprintf(stderr, "\n\nALL SUCCESS !\n");

  return 0;
}
