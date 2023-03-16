SHELL = /bin/sh

#CFLAGS := -g -O0 -Wall $(CFG) -std=c99
CFLAGS := -O3 -lm -ldl -Wall $(CFG) -std=c99 -fgnu89-inline -fopenmp

# Ensure 'mpicc' (ex on magique4) is available otherwise use 'cc' (ex on edison at NERSC)
ifeq "$(shell which mpicc 2> /dev/null)" ""
CC := cc
OPTIONAL_FLAGS := -dynamic
else
CC := mpicc
OPTIONAL_FLAGS := -lm
endif

# --- FFTW package ---
FFTW_BASE_DIR := ../fftw
FFTW_NAME := fftw-3.3.4
FFTWDIR := $(FFTW_NAME)
FFTWDIR_INCLUDE := $(FFTW_BASE_DIR)/$(FFTWDIR)/api
FFTW_TAR_FILE := $(FFTW_NAME).tgz
FFTW_LIB := $(FFTW_BASE_DIR)/$(FFTWDIR)/.libs/libfftw3.a
FFTW_OMPLIB := $(FFTW_BASE_DIR)/$(FFTWDIR)/threads/.libs/libfftw3_omp.a -lm



DIRPYTHONPATH := /home1/datahome/jmdeloui/conda-env/healpy
DIRPYTHONINC := /home1/datahome/jmdeloui/conda-env/healpy/include/python3.8
DIRNUMPYINC  := /home1/datahome/jmdeloui/conda-env/healpy/lib/python3.8/site-packages/numpy/core/include
DIRPYTHONLIB := $(DIRPYTHONPATH)/lib
LIBPYTHONLIB := python3.8
PYTHONCONF   := $(DIRPYTHONPATH)/bin/python3.8-config
OPTIONPYTHON := -DPYTHON3

# --- CFITSIO package ---
CFITSIO_ROOT := $(SROLLDIR)/cfitsio
CFITSIO_TAR  := $(CFITSIO_ROOT)/cfitsio3410.tar.gz
CFITSIO_DIR  := $(CFITSIO_ROOT)/cfitsio
CFITSIO_LIB  := $(CFITSIO_DIR)/libcfitsio.a

# --- NO_DMC_LIB package ---
NO_DMC_LIB_BASE_DIR := ../nodmc/nodmclib
NO_DMC_LIB_DIR_INCLUDE := $(NO_DMC_LIB_BASE_DIR)/NO_DMC_LIB
NO_DMC_LIB_LIB_DIR := $(NO_DMC_LIB_BASE_DIR)/Linux-x86_64
NO_DMC_LIB_LIB_FILE := $(NO_DMC_LIB_LIB_DIR)/libno_dmc.so

# --- PAR_LOADER package ---
PAR_LOADER_BASE_DIR := ../nodmc/parloader
# Config for stim!
PAR_LOADER_INPUT_HEADER_DEF := $(shell pwd)/stim_param.h
PAR_LOADER_TAG_PREFIX := stim_

# HPR_Cleaner.o should have been built by sroll/hprclean before building libstim
STIMOBJ := $(PAR_LOADER_TAG_PREFIX)parLoader.o \
           ../RngStreams/RngStream.o \
           add_oof_and_white_noise.o \
           add_oof_noise.o \
           sim_adu.o \
           compress_decompress.o \
           correct_adc.o \
           adu_to_volts.o \
           despike_flag.o \
           bl_demod.o \
           gain_decorr.o \
           taudeconv_lib.o \
           real_filter.o \
           convolve.o \
           gapfiller.o \
           stim_tools.o \
           ../sroll/chealpix.o \
           ../sroll/spline.o \
           ../sroll/HPR_Cleaner.o \
           stim.o
#           subtract_adu.o

STIMINCLUDES := -I$(FFTWDIR_INCLUDE) -I$(CFITSIO_DIR) -I$(NO_DMC_LIB_DIR_INCLUDE) -I. -L$(NO_DMC_LIB_LIB_DIR) -I../sroll -I../RngStreams -L.

# Some usefull var for guessing the dir where to find the shared lib (libstim.so)
mkfile_path := $(dir $(realpath $(lastword $(MAKEFILE_LIST))))



all: always_check #stimexe bin2map bin2polmap

fftw: $(FFTW_OMPLIB)

$(FFTW_OMPLIB):
	cd $(FFTW_BASE_DIR) && tar zxf $(FFTW_TAR_FILE)
	cd $(FFTW_BASE_DIR)/$(FFTWDIR) && ./configure --enable-threads --enable-openmp
	$(MAKE) -C $(FFTW_BASE_DIR)/$(FFTWDIR)

cfitsio: $(CFITSIO_LIB)

$(CFITSIO_LIB) :
	cd $(CFITSIO_ROOT) && $(TAR) zxf $(CFITSIO_TAR)
	cd $(CFITSIO_DIR) && ./configure
	$(MAKE) -C $(CFITSIO_DIR)

always_check: fftw cfitsio $(PAR_LOADER_TAG_PREFIX)parLoader.c
	$(MAKE) -C $(NO_DMC_LIB_BASE_DIR)

$(PAR_LOADER_TAG_PREFIX)parLoader.c: $(PAR_LOADER_INPUT_HEADER_DEF)
	$(MAKE) -C $(PAR_LOADER_BASE_DIR) currentParContentDefHeaderFile="$(PAR_LOADER_INPUT_HEADER_DEF)" tagPrefix="$(PAR_LOADER_TAG_PREFIX)" 

%.o: %.c $(PAR_LOADER_TAG_PREFIX)parLoader.c always_check
	$(CC) $(STIMINCLUDES) -c -o $@ $< $(CFLAGS) -fpic


libstim.so: $(STIMOBJ)
	$(CC) -shared -o $@ $^  $(CFITSIO_LIB) 

#stimexe: stimexe.c libstim.so $(NO_DMC_LIB_LIB_FILE)
	$(CC) $(CFLAGS) $(STIMINCLUDES) $(OPTIONAL_FLAGS) -Wl,-R$(NO_DMC_LIB_LIB_DIR):$(mkfile_path) -o $@ $< -lno_dmc -lstim $(FFTW_LIB) $(FFTW_OMPLIB) 

#bin2map: bin2map.c $(NO_DMC_LIB_LIB_FILE)
	$(CC) $(CFLAGS) $(STIMINCLUDES) $(OPTIONAL_FLAGS) -Wl,-R$(NO_DMC_LIB_LIB_DIR):$(mkfile_path) -o $@ $< -lno_dmc -lstim $(FFTW_LIB) $(FFTW_OMPLIB) `$(PYTHONCONF) --cflag

#bin2polmap: bin2polmap.c $(NO_DMC_LIB_LIB_FILE)
	$(CC) $(CFLAGS) $(STIMINCLUDES) $(OPTIONAL_FLAGS) -Wl,-R$(NO_DMC_LIB_LIB_DIR):$(mkfile_path) -o $@ $< -lno_dmc -lstim $(FFTW_LIB) $(FFTW_OMPLIB) `$(PYTHONCONF) --cflags` -lpython3.8 `$(PYTHONCONF) --ldflags`

clean:
	$(RM) *.o *.so stimexe bin2map bin2polmap "$(PAR_LOADER_TAG_PREFIX)parLoader.c" "$(PAR_LOADER_TAG_PREFIX)parLoader.h"
