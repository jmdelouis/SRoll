SHELL = /bin/sh

#CFLAGS := -g -O0 -Wall $(CFG) -std=c99 -fgnu89-inline -fopenmp -ffunction-sections -fdata-sections -fprofile-arcs -ftest-coverage
#CFLAGS :=  -O3 -Wall $(CFG) -std=c99 -fgnu89-inline -fopenmp -ffunction-sections -fdata-sections
CFLAGS := -O3 -Wall $(CFG) -std=c99 -fgnu89-inline -fopenmp -w



# Ensure 'mpicc' (ex on magique4) is available otherwise use 'cc' (ex on edison at NERSC)
ifeq "$(shell which mpicc 2> /dev/null)" ""
CC := gcc
OPTIONAL_FLAGS := -dynamic -Wl,--print-gc-sections
else
CC := mpicc -Wunused-result
OPTIONAL_FLAGS := -Wl,--print-gc-sections -Wunused-result
endif

TAR := tar

ARCHIVE_FILENAME := sroll_`date +"%Y-%m-%d"`.tgz

BUILDDIR := .


# --- FFTW package ---
FFTW_BASE_DIR := $(SROLLDIR)/fftw
FFTW_NAME := fftw-3.3.4
FFTWDIR := $(FFTW_NAME)
FFTWDIR_INCLUDE := $(FFTW_BASE_DIR)/$(FFTWDIR)/api
FFTW_TAR_FILE := $(FFTW_NAME).tgz
FFTW_LIB := $(FFTW_BASE_DIR)/$(FFTWDIR)/.libs/libfftw3.a
FFTW_OMPLIB := $(FFTW_BASE_DIR)/$(FFTWDIR)/threads/.libs/libfftw3_omp.a

# --- CFITSIO package ---
CFITSIO_ROOT := $(SROLLDIR)/cfitsio
CFITSIO_TAR  := $(CFITSIO_ROOT)/cfitsio3410.tar.gz
CFITSIO_DIR  := $(CFITSIO_ROOT)/cfitsio
CFITSIO_LIB  := $(CFITSIO_DIR)/libcfitsio.a

# --- Python wrapper ---
ifeq ("$(SROLLHOST)","OCCIGEN")
DIRPYTHONPATH := /opt/software/occigen/tools/python/2.7.12/intel/17.0/nompi
DIRPYTHONINC := $(DIRPYTHONPATH)/include/python2.7
DIRNUMPYINC  := $(DIRPYTHONPATH)/lib/python2.7/site-packages/numpy/core/include
DIRPYTHONLIB := $(DIRPYTHONPATH)/lib64/python2.7/config
LIBPYTHONLIB := python2.7
OPTIONPYTHON := -DPYTHON2
endif 

# -- Test debug path python garoupe --
ifeq ("$(SROLLHOST)","GAROUPE")
DIRPYTHONPATH := /export/home/tfoulquier/env_sroll #/home3/homedir7/perso/tfoulqui/workspace/my_env #/export/home/jmdeloui/anaconda3/envs/healpy
DIRPYTHONINC := /export/home/tfoulquier/env_sroll/include/python3.6m #/home3/homedir7/perso/tfoulqui/workspace/my_env/include/python2.7 #/export/home/jmdeloui/anaconda3/envs/healpy/include/python3.8
DIRNUMPYINC  := /export/home/tfoulquier/env_sroll/lib/python3.6/site-packages/numpy/core/include #/home3/homedir7/perso/tfoulqui/workspace/my_env/lib/python2.7/site-packages/numpy/core/include #/python3.8
DIRPYTHONLIB := /export/home/tfoulquier/env_sroll/lib/python3.6/config-3.6m-x86_64-linux-gnu #/home3/homedir7/perso/tfoulqui/workspace/my_env/lib/python2.7/config #/python3.8/config
LIBPYTHONLIB := python3.6 #python2.7 #/python3.8
OPTIONPYTHON := -DPYTHON3 #-DPYTHON2 #-DPYTHON3
endif


# -- Test debug path python datarmor --
ifeq ("$(SROLLHOST)","DATARMOR")
#DEBUG 
#DIRPYTHONPATH :=/home1/datahome/jmdeloui/conda-env/healpy
#DIRPYTHONINC :=/home1/datahome/jmdeloui/conda-env/healpy/include/python3.8
#DIRNUMPYINC  :=/home1/datahome/jmdeloui/conda-env/healpy/lib/python3.8/site-packages/numpy/core/include
#DIRPYTHONLIB := $(DIRPYTHONPATH)/lib
#LIBPYTHONLIB := python3.8
#OPTIONPYTHON := -DPYTHON3

PYTHONPATH := /home1/datahome/tfoulqui/.conda/envs/sroll_env
DIRPYTHONINC := $(PYTHONPATH)/include/python2.7
DIRNUMPYINC  := $(PYTHONPATH)/lib/python2.7/site-packages/numpy/core/include
DIRPYTHONLIB := $(PYTHONPATH)/lib/python2.7
LIBPYTHONLIB := python2.7
endif

# --- NO_DMC_LIB package ---
NO_DMC_LIB_NAME := no_dmc
NO_DMC_LIB_BASE_DIR := $(SROLLDIR)/nodmc/nodmclib
NO_DMC_LIB_DIR_INCLUDE := $(NO_DMC_LIB_BASE_DIR)/NO_DMC_LIB
NO_DMC_LIB_DIR := $(NO_DMC_LIB_BASE_DIR)/Linux-x86_64
NO_DMC_LIB_FILE := $(NO_DMC_LIB_DIR)/libno_dmc.so

# --- PAR_LOADER package ---
PAR_LOADER_BASE_DIR := ../nodmc/parloader
# Config for sroll!
PAR_LOADER_INPUT_HEADER_DEF := $(shell pwd)/sroll_param.h
PAR_LOADER_TAG_PREFIX := sroll_
# Config for troll!
PAR_LOADER_INPUT_HEADER_TROLL_DEF := $(shell pwd)/troll_param.h
PAR_LOADER_TAG_PREFIX_TROLL := troll_
# Config for hprclean!
PAR_LOADER_INPUT_HEADER_HPRCLEAN_DEF := $(shell pwd)/hprclean_param.h
PAR_LOADER_TAG_PREFIX_HPRCLEAN := hprclean_

# --- STIM package ---
STIM_LIB_NAME := stim
STIM_BASE_DIR := ../stim
STIM_DIR_INCLUDE := $(STIM_BASE_DIR)
STIM_LIB_DIR := $(STIM_BASE_DIR)
STIM_LIB_FILE := $(STIM_LIB_DIR)/stim.so


ODIR := $(BUILDDIR)
_OBJ := chealpix.o $(PAR_LOADER_TAG_PREFIX)parLoader.o spline.o
OBJ := $(patsubst %,$(ODIR)/%,$(_OBJ))
_OBJTROLL := chealpix.o $(PAR_LOADER_TAG_PREFIX_TROLL)parLoader.o spline.o
OBJTROLL := $(patsubst %,$(ODIR)/%,$(_OBJTROLL))
_OBJHPRCLEAN := $(PAR_LOADER_TAG_PREFIX_HPRCLEAN)parLoader.o spline.o HPR_Cleaner.o
OBJHPRCLEAN := $(patsubst %,$(ODIR)/%,$(_OBJHPRCLEAN))

CFGSROLL := -I$(FFTWDIR_INCLUDE) -I$(CFITSIO_DIR) -I$(NO_DMC_LIB_DIR_INCLUDE) -I$(DIRPYTHONINC) -I$(DIRNUMPYINC) -L$(DIRPYTHONLIB) -L$(NO_DMC_LIB_DIR) -I$(STIM_DIR_INCLUDE) -L$(STIM_LIB_DIR) -I../RngStreams -I. 


.PHONY: all clean dist fftw cfitsio levels cleanAll alwaysCheckModules


# build hprclean before sroll, because sroll uses stimlib which needs HPR_Clean.o
all: sroll srollit srollit545 merge troll troll_14tf


# DEBUG rule, to print a variable computed in this Makefile, for example OBJ, do: make print-OBJ
print-%  : ; @echo $* = $($*)


# This rule allow to ensure the proper (re)compilation of modules that sroll depend on.
# For now the relevant modules are: no_dmc_lib and stim.
alwaysCheckModules: fftw cfitsio sroll_parLoader.c
	$(MAKE) -C $(NO_DMC_LIB_BASE_DIR)
	$(MAKE) -C $(STIM_BASE_DIR)

fftw: $(FFTW_OMPLIB)

$(FFTW_OMPLIB) :
	cd $(FFTW_BASE_DIR) && $(TAR) zxf $(FFTW_TAR_FILE)
	cd $(FFTW_BASE_DIR)/$(FFTWDIR) && ./configure --enable-threads --enable-openmp
	$(MAKE) -C $(FFTW_BASE_DIR)/$(FFTWDIR)

cfitsio: $(CFITSIO_LIB)

$(CFITSIO_LIB) :
	cd $(CFITSIO_ROOT) && $(TAR) zxf $(CFITSIO_TAR)
	cd $(CFITSIO_DIR) && ./configure
	$(MAKE) -C $(CFITSIO_DIR) all shared


####
#### LevelS
####
LSDIR := $(SROLLDIR)/LevelS/LevelS-2015-01-09
levels: $(LSDIR)/linux_gcc/bin/multimod
# extract LevelS, patch dipole.cc, make
$(LSDIR)/linux_gcc/bin/multimod:
	cd $(SROLLDIR)/LevelS && tar -jxvf LevelS-2015-01-09.tar.bz2
	mv $(LSDIR)/cxxmod/dipole.cc $(LSDIR)/cxxmod/dipole_levels.cc
	cp $(SROLLDIR)/LevelS/dipole_bware.cc $(LSDIR)/cxxmod/dipole.cc
	cd $(LSDIR) && export LEVELS_TARGET=linux_gcc && make


# Note that the generated parLoader source code depend on the parameter file header definition.
# This rule will generate the 2 required files: "sroll_parLoader.c" and "sroll_parLoader.h"
sroll_parLoader.c: $(PAR_LOADER_INPUT_HEADER_DEF)
	$(MAKE) -C $(PAR_LOADER_BASE_DIR) currentParContentDefHeaderFile="$(PAR_LOADER_INPUT_HEADER_DEF)" tagPrefix="$(PAR_LOADER_TAG_PREFIX)"

hprclean_parLoader.c: $(PAR_LOADER_INPUT_HEADER_HPRCLEAN_DEF)
	$(MAKE) -C $(PAR_LOADER_BASE_DIR) currentParContentDefHeaderFile="$(PAR_LOADER_INPUT_HEADER_HPRCLEAN_DEF)" tagPrefix="$(PAR_LOADER_TAG_PREFIX_HPRCLEAN)"

troll_parLoader.c: $(PAR_LOADER_INPUT_HEADER_TROLL_DEF)
	$(MAKE) -C $(PAR_LOADER_BASE_DIR) currentParContentDefHeaderFile="$(PAR_LOADER_INPUT_HEADER_TROLL_DEF)" tagPrefix="$(PAR_LOADER_TAG_PREFIX_TROLL)"

%.o: %.c alwaysCheckModules
	#$(CC) -c -o $@ $< $(CFLAGS) $(CFGSROLL) -Werror -fpic -lm
	$(CC) -c -o $@ $< $(CFLAGS) $(CFGSROLL) -fpic
#hprclean: alwaysCheckModules
#hprclean: hprclean.c $(OBJHPRCLEAN)
#	$(CC) $(CFLAGS) $(CFGSROLL) $(OPTIONAL_FLAGS) -Wl,-R$(NO_DMC_LIB_DIR):$(STIM_LIB_DIR)  -l$(NO_DMC_LIB_NAME) -l$(STIM_LIB_NAME) -o $@ $< $(OBJ) -l$(STIM_LIB_NAME) $(FFTW_LIB) $(FFTW_OMPLIB) $(CFITSIO_LIB)  -l$(NO_DMC_LIB_NAME) -lm

sroll: alwaysCheckModules
sroll: sroll.c $(OBJ)
	$(CC) $(CFLAGS) $(CFGSROLL) $(OPTIONAL_FLAGS) -Wl,-R$(NO_DMC_LIB_DIR):$(STIM_LIB_DIR) -l$(NO_DMC_LIB_NAME) -l$(STIM_LIB_NAME) -o $@ $< $(OBJ) -l$(STIM_LIB_NAME) $(FFTW_LIB) $(FFTW_OMPLIB) $(CFITSIO_LIB)  -l$(NO_DMC_LIB_NAME) -lm

srollit: alwaysCheckModules
srollit: sroll.c $(OBJ)
	$(CC) $(CFLAGS) $(CFGSROLL) $(OPTIONAL_FLAGS) -DMAXSIMU=5 -DMAXTHEOHPR=7 -Wl,-R$(NO_DMC_LIB_DIR):$(STIM_LIB_DIR) -o $@ $< $(OBJ) -l$(STIM_LIB_NAME) $(FFTW_LIB) $(FFTW_OMPLIB) $(CFITSIO_LIB) -l$(NO_DMC_LIB_NAME)  -lm

#troll: alwaysCheckModules
troll: troll.c $(OBJTROLL)
	$(CC) $(CFLAGS) $(CFGSROLL) $(OPTIONAL_FLAGS) -Wl,-R$(NO_DMC_LIB_DIR):$(STIM_LIB_DIR) -l$(NO_DMC_LIB_NAME) -l$(STIM_LIB_NAME) -o $@ $< $(OBJTROLL) -l$(STIM_LIB_NAME) $(FFTW_LIB) $(FFTW_OMPLIB) $(CFITSIO_LIB) -l$(NO_DMC_LIB_NAME) -l$(LIBPYTHONLIB) -lm -DPYTHON3

#troll_14tf: alwaysCheckModules
troll_14tf: troll.c $(OBJTROLL)
	$(CC) $(CFLAGS) $(CFGSROLL) $(OPTIONAL_FLAGS)  -DMAXTHEOHPR=14 -Wl,-R$(NO_DMC_LIB_DIR):$(STIM_LIB_DIR) -l$(NO_DMC_LIB_NAME) -l$(STIM_LIB_NAME) -o $@ $< $(OBJTROLL) -l$(STIM_LIB_NAME) $(FFTW_LIB) $(FFTW_OMPLIB) $(CFITSIO_LIB) -l$(NO_DMC_LIB_NAME) -l$(LIBPYTHONLIB) -lm -DPYTHON3

srollit545: alwaysCheckModules
srollit545: sroll.c $(OBJ)
	$(CC) $(CFLAGS) $(CFGSROLL) $(OPTIONAL_FLAGS) -DMAXSIMU=10 -DMAXTHEOHPR=9 -Wl,-R$(NO_DMC_LIB_DIR):$(STIM_LIB_DIR) -l$(NO_DMC_LIB_NAME) -l$(STIM_LIB_NAME) -o $@ $< $(OBJ) -l$(STIM_LIB_NAME) $(FFTW_LIB) $(FFTW_OMPLIB) $(CFITSIO_LIB) -l$(NO_DMC_LIB_NAME) -lm

projsroll: alwaysCheckModules
projsroll: projsroll.c $(OBJ)
	$(CC) $(CFLAGS) $(CFGSROLL) $(OPTIONAL_FLAGS) -Wl,-R$(NO_DMC_LIB_DIR):$(STIM_LIB_DIR) -l$(NO_DMC_LIB_NAME) -l$(STIM_LIB_NAME) -o $@ $< $(OBJ) $(FFTW_LIB) $(FFTW_OMPLIB) $(CFITSIO_LIB)

merge: alwaysCheckModules
merge: merge.c $(OBJ)
	$(CC) $(CFLAGS) $(CFGSROLL) $(OPTIONAL_FLAGS) -Wl,-R$(NO_DMC_LIB_DIR) -l$(NO_DMC_LIB_NAME) -o $@ $< $(OBJ) $(CFITSIO_LIB) -l$(NO_DMC_LIB_NAME) -lm

dist:
	$(RM) $(ARCHIVE_FILENAME)
	$(TAR) --exclude=".*" --exclude="CVS" --exclude="sroll*.tgz" --dereference -zcf $(ARCHIVE_FILENAME) *.c *.h Makefile

clean:
	$(RM) *.o projsroll sroll merge "$(PAR_LOADER_TAG_PREFIX)parLoader.c" "$(PAR_LOADER_TAG_PREFIX)parLoader.h" "$(PAR_LOADER_TAG_PREFIX_TROLL)parLoader.c" "$(PAR_LOADER_TAG_PREFIX_HPRCLEAN)parLoader.c" "$(PAR_LOADER_TAG_PREFIX_TROLL)parLoader.h"
	$(MAKE) -C $(STIM_BASE_DIR) clean

cleanAll: clean
#	$(MAKE) -C $(FFTW_BASE_DIR)/$(FFTWDIR) clean
	$(MAKE) -C $(CFITSIO_DIR) clean
	$(MAKE) -C $(NO_DMC_LIB_BASE_DIR) clean
#	$(MAKE) -C $(PAR_LOADER_BASE_DIR) clean
	$(MAKE) -C $(STIM_BASE_DIR) clean
	$(RM)   -rf $(FFTW_BASE_DIR)/$(FFTWDIR)
	$(RM)   -rf $(LSDIR)
