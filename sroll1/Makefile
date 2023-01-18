SHELL = /bin/sh

#CFLAGS := -g -O0 -Wall $(CFG) -std=c99
CFLAGS := -O3 -Wall $(CFG) -std=c99 -fgnu89-inline -fopenmp

# Ensure 'mpicc' (ex on magique4) is available otherwise use 'cc' (ex on edison at NERSC)
ifeq "$(shell which mpicc 2> /dev/null)" ""
CC := cc
OPTIONAL_FLAGS := -dynamic
else
CC := mpicc
OPTIONAL_FLAGS :=
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

# --- STIM package ---
STIM_LIB_NAME := stim
STIM_BASE_DIR := ../stim
STIM_DIR_INCLUDE := $(STIM_BASE_DIR)
STIM_LIB_DIR := $(STIM_BASE_DIR)
STIM_LIB_FILE := $(STIM_LIB_DIR)/stim.so


ODIR := $(BUILDDIR)
_OBJ := chealpix.o $(PAR_LOADER_TAG_PREFIX)parLoader.o
OBJ := $(patsubst %,$(ODIR)/%,$(_OBJ))

CFGSROLL := -I$(FFTWDIR_INCLUDE) -I$(CFITSIO_DIR) -I$(NO_DMC_LIB_DIR_INCLUDE) -L$(NO_DMC_LIB_DIR) -I$(STIM_DIR_INCLUDE) -L$(STIM_LIB_DIR) -I../RngStreams -I.


.PHONY: all clean dist fftw cfitsio cleanAll alwaysCheckModules


all: sroll


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
	$(MAKE) -C $(CFITSIO_DIR)

# Note that the generated parLoader source code depend on the parameter file header definition.
# This rule will generate the 2 required files: "sroll_parLoader.c" and "sroll_parLoader.h"
sroll_parLoader.c: $(PAR_LOADER_INPUT_HEADER_DEF)
	$(MAKE) -C $(PAR_LOADER_BASE_DIR) currentParContentDefHeaderFile="$(PAR_LOADER_INPUT_HEADER_DEF)" tagPrefix="$(PAR_LOADER_TAG_PREFIX)"

%.o: %.c alwaysCheckModules
	$(CC) -c -o $@ $< $(CFLAGS) $(CFGSROLL) -Werror

sroll: alwaysCheckModules
sroll: sroll.c $(OBJ)
	$(CC) $(CFLAGS) $(CFGSROLL) $(OPTIONAL_FLAGS) -Wl,-R$(NO_DMC_LIB_DIR):$(STIM_LIB_DIR) -l$(NO_DMC_LIB_NAME) -l$(STIM_LIB_NAME) -o $@ $< $(OBJ) $(FFTW_LIB) $(FFTW_OMPLIB) $(CFITSIO_LIB)

srollit: alwaysCheckModules
srollit: sroll.c $(OBJ)
	$(CC) $(CFLAGS) $(CFGSROLL) $(OPTIONAL_FLAGS) -DMAXSIMU=5 -DMAXTHEOHPR=7 -Wl,-R$(NO_DMC_LIB_DIR):$(STIM_LIB_DIR) -l$(NO_DMC_LIB_NAME) -l$(STIM_LIB_NAME) -o $@ $< $(OBJ) $(FFTW_LIB) $(FFTW_OMPLIB) $(CFITSIO_LIB)

srollit545: alwaysCheckModules
srollit545: sroll.c $(OBJ)
	$(CC) $(CFLAGS) $(CFGSROLL) $(OPTIONAL_FLAGS) -DMAXSIMU=10 -DMAXTHEOHPR=9 -Wl,-R$(NO_DMC_LIB_DIR):$(STIM_LIB_DIR) -l$(NO_DMC_LIB_NAME) -l$(STIM_LIB_NAME) -o $@ $< $(OBJ) $(FFTW_LIB) $(FFTW_OMPLIB) $(CFITSIO_LIB) 

projsroll: alwaysCheckModules
projsroll: projsroll.c $(OBJ)
	$(CC) $(CFLAGS) $(CFGSROLL) $(OPTIONAL_FLAGS) -Wl,-R$(NO_DMC_LIB_DIR):$(STIM_LIB_DIR) -l$(NO_DMC_LIB_NAME) -l$(STIM_LIB_NAME) -o $@ $< $(OBJ) $(FFTW_LIB) $(FFTW_OMPLIB) $(CFITSIO_LIB)

merge: alwaysCheckModules
merge: merge.c $(OBJ)
	$(CC) $(CFLAGS) $(CFGSROLL) $(OPTIONAL_FLAGS) -Wl,-R$(NO_DMC_LIB_DIR) -l$(NO_DMC_LIB_NAME) -o $@ $< $(OBJ) $(CFITSIO_LIB)

dist:
	$(RM) $(ARCHIVE_FILENAME)
	$(TAR) --exclude=".*" --exclude="CVS" --exclude="sroll*.tgz" --dereference -zcf $(ARCHIVE_FILENAME) *.c *.h Makefile

clean:
	$(RM) *.o projsroll sroll merge "$(PAR_LOADER_TAG_PREFIX)parLoader.c" "$(PAR_LOADER_TAG_PREFIX)parLoader.h"
	$(MAKE) -C $(STIM_BASE_DIR) clean

cleanAll: clean
	$(MAKE) -C $(FFTW_BASE_DIR)/$(FFTWDIR) clean
	$(MAKE) -C $(CFITSIO_DIR) clean
	$(MAKE) -C $(NO_DMC_LIB_BASE_DIR) clean
#	$(MAKE) -C $(PAR_LOADER_BASE_DIR) clean
	$(MAKE) -C $(STIM_BASE_DIR) clean
