##############################################################################
# fillParLoader.py
#
# This python script get generic files of the parLoader and fill them with
# appropriate data then save them to target files.
#
# Date    : 2015-04-06
# version : BETA
##############################################################################

class ParamInfo:
  """Class representing a parameter"""
  name = None
  islist = False
  mandatory = True
  type = None

# ---------------------------------------------------------------------------
def getTypeParam(param):

  if isinstance(param,str):
    paramType= "PIOSTRING"
  elif isinstance(param,int):
    paramType= "PIOINT"
  elif isinstance(param,float) :
    paramType = "PIOFLOAT"
  elif isinstance(param,double):
    paramType = "PIODOUBLE"
  elif isinstance(param,long):
    paramType = "PIOLONG"


  return paramType
#-----------------------------------------------------------------------------
def fillParamInfo_old(headerParContent):
  paramInfoList = [] # The structure to store parameters info
  # We will only consider values between these symbols
  BEGIN_STATE_EXPR = "typedef struct"
  END_STATE_EXPR = "parContent;"
  validState = False
  PARAM_PREFIX_FLAG = "flag_"
  PARAM_PREFIX_SIZE = "n_"

  fileH = open(headerParContent)

  # Loop on each line of the header file
  for line in fileH:
    # Ignore empty or comment lines
    sline = line.strip()
    if sline == '' or sline.startswith('*') or sline.startswith('/*') or sline.startswith('//'):
      continue
  
    if (validState == False) and (BEGIN_STATE_EXPR in line):
      validState = True
      #print("validState becomes True!")
      continue
  
    if (validState == True) and (END_STATE_EXPR in line):
      validState = False
      #print("validState becomes False!")
      break # we know there is no more to read!

    # Ignore line outside valid range
    if validState == False:
      continue

    #print("[DEBUG] We found a VALID LINE: %s" % line) 
    tmpParamInfo = ParamInfo()
    # Set default values
    #tmpParamInfo.islist = False
    #tmpParamInfo.mandatory = True

    # Remove everything after ';'
    sline = sline.split(';',1)[0]
    #print("[DEBUG] split = '%s'" % sline)

    # Remove '*' but keep that it indicate a list!
    if '*' in sline:
      tmpParamInfo.islist = True
      sline = sline.replace('*','')
      #print("[DEBUG] this is a list")

    # Retrieve pio type
    # Retrieve param name
    p_type, p_name = sline.split()
    tmpParamInfo.name = p_name
    tmpParamInfo.type = p_type
    
    #print("[DEBUG] p_type = '%s', p_name = '%s'" % (p_type, p_name))

    # Check for special prefix on param name

    # - PARAM_PREFIX_FLAG
    if p_name.startswith(PARAM_PREFIX_FLAG) and (p_type == "PIOBYTE"):
      del tmpParamInfo # the tmp object is of no use
      # In this case we must find the corresponding parameter in the list and update its info
      tmpp = [x for x in paramInfoList if x.name == p_name[len(PARAM_PREFIX_FLAG):]][0]
      tmpp.mandatory = False
      continue
    
    # - PARAM_PREFIX_SIZE
    if p_name.startswith(PARAM_PREFIX_SIZE):
      del tmpParamInfo # the tmp object is of no use
      # In this case we must find the corresponding parameter in the list and update its info
      tmpp = [x for x in paramInfoList if x.name == p_name[len(PARAM_PREFIX_SIZE):]][0]
      
      if (tmpp.islist == False):
        print("ERROR: param '%s' should be a pointer to store list!" % p_name[len(PARAM_PREFIX_FLAG):])
        exit(1)
      #tmpp.list = True # not usefull since we already know that...
      continue
    
    # Add newly updated param
    paramInfoList.append(tmpParamInfo)

  return paramInfoList
  
#-----------------------------------------------------------------------------
def fillParamInfo(headerParContent):
  paramInfoList = [] # The structure to store parameters info
  # We will only consider values between these symbols
  BEGIN_STATE_EXPR = "typedef struct"
  END_STATE_EXPR = "parContent;"
  validState = False
  PARAM_PREFIX_FLAG = "flag_"
  PARAM_PREFIX_SIZE = "n_"

  #plip

  #debug path  
  path_py_param ="/home3/homedir7/perso/tfoulqui/workspace/srollexx_work/sroll4/"
  sys.path.append(path_py_param)
  #sys.path.append(headerParContent)

  #import python params
  from params_unit_test import main
  pyValues = main()

  #parcours dict params python
  for k,v in pyValues.items():
    #init ParamInfo
    tmpParamInfo = ParamInfo()
    #set default value 
    tmpParamInfo.islist = False
    #if list
    if isinstance(v,list):
      tmpParamInfo.islist = True

    tmpParamInfo.name = k
    #tmpParamInfo.type = type(v) 
    if isinstance(v,list):
      tmpParamInfo.type=getTypeParam(v[0])
    else:
      tmpParamInfo.type=getTypeParam(v)

    # Add newly updated param
    paramInfoList.append(tmpParamInfo)

  return paramInfoList

#-----------------------------------------------------------------------------

def codeC_paramDef_list(paramInfoList):
  """Convert the already retrieve info about params to the corresponding C code."""

# Example:

#  paramDef paramDef_list[] = {
#    {"fsl", false, false}, // PIOSTRING (list), optional
#    {"NADU", true, false}, // PIOINT, mandatory
#    {"Calibration", true, false}, // PIODOUBLE (list), mandatory
#  };
#  int paramDef_list_size = 3;

  codeC = "paramDef __TAG__paramDef_list[] = {\n"

  for p in paramInfoList:
    codeC += "  {\"%s\", %s, %s, %s},\n" % (p.name, repr(p.islist).lower(), repr(p.mandatory).lower(), "false")

  # Add final code
  codeC += "};\nint __TAG__paramDef_list_size = %d;\n" % len(paramInfoList)
  
  return codeC

#-----------------------------------------------------------------------------

def getCodeC_affectation_forType(name, type, isList):

  result = ""

  additionalSpecifier = ""
  prefixValue = "*"
  indent = ""
  # if isList then this is a list! So we need of additional specifier
  if isList:
    additionalSpecifier = "[i]"
    prefixValue = ""
    indent = "  "

  if type != "PIOSTRING":
    result += indent+"    errno = 0;\n"

  result +=  {
    #'PIODOUBLE': indent+"    param->%s%s = strtod(%svalue%s, NULL);\n" % (name, additionalSpecifier, prefixValue, additionalSpecifier)
    'PIODOUBLE': indent+"    param->%s%s = myRead_PIODOUBLE(%svalue%s);\n" % (name, additionalSpecifier, prefixValue, additionalSpecifier)
                ,
    #'PIOLONG':   indent+"    param->%s%s = strtol(%svalue%s, NULL, 10);\n" % (name, additionalSpecifier, prefixValue, additionalSpecifier)
    'PIOLONG':   indent+"    param->%s%s = myRead_PIOLONG(%svalue%s);\n" % (name, additionalSpecifier, prefixValue, additionalSpecifier)
                ,
    #'PIOFLOAT':  indent+"    param->%s%s = strtof(%svalue%s, NULL);\n" % (name, additionalSpecifier, prefixValue, additionalSpecifier)
    'PIOFLOAT':  indent+"    param->%s%s = myRead_PIOFLOAT(%svalue%s);\n" % (name, additionalSpecifier, prefixValue, additionalSpecifier)
                ,
    'PIOSTRING': indent+"    strcpy(param->%s%s, %svalue%s);\n" % (name, additionalSpecifier, prefixValue, additionalSpecifier)
                ,
    #'PIOINT':    indent+"    param->%s%s = strtol(%svalue%s, NULL, 10);\n" % (name, additionalSpecifier, prefixValue, additionalSpecifier)
    'PIOINT':    indent+"    param->%s%s = myRead_PIOINT(%svalue%s);\n" % (name, additionalSpecifier, prefixValue, additionalSpecifier)
  }[type]

  if type != "PIOSTRING":
    result +=  indent+"    if (errno != 0) {\n"\
              +indent+"      fprintf(stderr, \"ERROR: '%s': Unable to convert value '%s' to target type %s\\n\", %s, \"%s\");\n" % (name, "%s", "%s", prefixValue+"value"+additionalSpecifier, type)\
              +indent+"      return 1;\n"\
              +indent+"    }\n"

  return result

#-----------------------------------------------------------------------------

def codeC_param_affectation(paramInfoList):
  """Convert the already retrieve info about params to the corresponding C code."""

# Example:

#  // Copy into correct entry in parameters struct
#  if (strcmp(name, "fsl") == 0) {
#    //strncpy(param->fsl, (char*)value, 100); //TODO to check for an array os PIOSTRING
#    param->fsl = (PIOSTRING *)value; //TODO to check for an array os PIOSTRING
#    param->flag_fsl = _PAR_TRUE;
#    param->n_fsl = list_size;
#  } else if (strcmp(name, "NADU") == 0) {
#    param->NADU = strtol(*value, NULL, 10);
#  } else if (strcmp(name, "Calibration") == 0) {
#    param->n_Calibration = list_size;
#    param->Calibration = malloc(list_size * sizeof(PIODOUBLE));
#    if (param->Calibration == NULL) {
#      perror("Error");
#      return 1;
#    }
#    int i;
#    for (i = 0; i < list_size; i++) {
#      errno = 0;
#      param->Calibration[i] = strtod(value[i], NULL);
#      if (errno != 0) {
#        fprintf(stderr, "ERROR: \"%s%d\": Unable to convert value '%s' to target type %s\n", "Calibration", i, value[i], "PIODOUBLE");
#        return 1;
#      }
#    }
#  } else {
#    fprintf(stderr, "ERROR: %s: Unknown param name!\n", name);
#    return 1;
#  }


  codeC = "// Copy into correct entry in parameters struct\n"

  for i, p in enumerate(paramInfoList):
    conditionInstruction = "else if"
    if i == 0: # For the very first test
      conditionInstruction = "if"
      
    # Open the if statement
    codeC += "  %s (strcmp(name, \"%s\") == 0) {\n" % (conditionInstruction, p.name)

    # Handle case of not mandatory (ie. flag_ param)
    if not p.mandatory:
      codeC += "    param->flag_%s = _PAR_TRUE;\n" % (p.name)
    # Handle case of list param
    if p.islist:
      codeC += "    param->n_%s = list_size;\n" % (p.name)
      codeC += "    param->%s = malloc(list_size * sizeof(%s));\n" % (p.name, p.type)\
              +"    if (param->%s == NULL) {\n" % (p.name)\
              +"      perror(\"Error\");\n"\
              +"      return 1;\n"\
              +"    }\n"
      codeC += "    int i;\n"\
              +"    for (i = 0; i < list_size; i++) {\n"\
              + getCodeC_affectation_forType(p.name, p.type, True)\
              +"    }\n"
    # Handle case of simple param
    else:
      codeC += getCodeC_affectation_forType(p.name, p.type, False) # 'isList' set to False if not list!
    
    # Close the if statement
    codeC += "  }\n"

  # Add final code
  codeC += ('  else {\n'
            '    fprintf(stderr, \"ERROR: %s: Unknown param name!\\n\", name);\n'
            '    return 1;\n'
            '  }\n'
            )
  
  return codeC


#-----------------------------------------------------------------------------

def replace_auto_gen_param_def_full_list(s, auto_content):
  return s.replace("__AUTO_GEN_PARAM_DEF_FULL_LIST__", auto_content)

#-----------------------------------------------------------------------------

def replace_auto_gen_param_affectation(s, auto_content):
  return s.replace("__AUTO_GEN_PARAM_AFFECTATION__", auto_content)

#-----------------------------------------------------------------------------

def replaceHeaderTag(s):
  import os.path
  return s.replace("__INCLUDE_PARCONTENT_HEADER__", "#include \"%s\"" % os.path.basename(headerParContent))

#-----------------------------------------------------------------------------

def replaceTagPrefix(s, upperCase=False):
  if upperCase:
    return s.replace("__TAG_UP__", tagPrefix.upper())
  else:
    return s.replace("__TAG__", tagPrefix)

#-----------------------------------------------------------------------------

def savefileWithContent(fileoutput, content):
  print("Saving file \"%s\"..." % fileoutput)

  f = open(fileoutput, 'w')
  f.write(content)
  f.flush()
  f.close()

#-----------------------------------------------------------------------------

def printUsage(cmd):
  print("Usage: %s </path_to/header_file_parContent.h> <tagPrefix>" % cmd)

#-----------------------------------------------------------------------------


##############################################################################
# MAIN
# This code allow to execute the module as a standalone script.
##############################################################################
if __name__ == "__main__":
  import sys
  #import shutil
  import os.path

  print("")
  print("************************************************")
  print("* %s" % sys.argv[0])
  print("* Goal: generate a fully functionnal parLoader")
  print("************************************************")
  print("")

  if len(sys.argv) != 3:
    print("ERROR: wrong number of argument!")
    printUsage(sys.argv[0])
    exit(1)

  SUFFIX_GEN = "_GEN"

  headerParContent = sys.argv[1]
  headerParLoaderGEN = "parLoader.h_GEN"
  sourceParLoaderGEN = "parLoader.c_GEN"
  tagPrefix = sys.argv[2]
  outputHeaderParLoader = os.path.dirname(headerParContent) + "/" + tagPrefix + headerParLoaderGEN[:-len(SUFFIX_GEN)]
  outputSourceParLoader = os.path.dirname(headerParContent) + "/" + tagPrefix + sourceParLoaderGEN[:-len(SUFFIX_GEN)]

  print(">> Using as input header for parContent: \"%s\"" % headerParContent)
  print(">> Using as generic header parLoader: \"%s\"" % headerParLoaderGEN)
  print(">> Using as generic source parLoader: \"%s\"" % sourceParLoaderGEN)
  print(">> Using as tag prefix for file and function: \"%s\"" % tagPrefix)
  print(">> Will generate header: \"%s\"" % outputHeaderParLoader)
  print(">> Will generate source: \"%s\"" % outputSourceParLoader)
  print("")

  # 1- Retrieve info
  print("Retrieving data from header file...")
  paramInfoList = fillParamInfo(headerParContent)
  """  
  print("[DEBUG] Finally we found %d parameters:" % len(paramInfoList))
  print("   %-20s %-10s %-6s %-6s" % ("NAME", "TYPE", "ISLIST", "MANDATORY"))
  for x in paramInfoList:
    print("p: %-20s %-10s %-6s %-6s " % (x.name, x.type, x.islist, x.mandatory))
  """
  auto_gen_param_def_full_list = codeC_paramDef_list(paramInfoList)
  auto_gen_param_affectation = codeC_param_affectation(paramInfoList)

  # 2- Load GEN file then replace its content appropriatly for saving to ouptut
  print("Loading input GEN files...")
  h = open(headerParLoaderGEN).read()
  s = open(sourceParLoaderGEN).read()

  print("Replacing content...")
  s = replace_auto_gen_param_def_full_list(s, auto_gen_param_def_full_list)
  s = replace_auto_gen_param_affectation(s, auto_gen_param_affectation)

  h = replaceHeaderTag(h)
  h = replaceTagPrefix(h, upperCase=True) # For the #ifndef...
  h = replaceTagPrefix(h) # For others...
  s = replaceHeaderTag(s)
  s = replaceTagPrefix(s) # For tag_prefix
  
  print("Writting result to output...")
  savefileWithContent(outputHeaderParLoader, h)
  savefileWithContent(outputSourceParLoader, s)
  
  print("")
  print("Done!")
  
