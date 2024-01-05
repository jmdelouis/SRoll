SHELL = /bin/sh
#CFLAGS := -g -O0 -Wall $(CFG) -std=c99 -fgnu89-inline -fopenmp -ffunction-sections -fdata-sections -fprofile-arcs -ftest-coverage
CFLAGS := -g -O3 -Wall $(CFG) -std=c99 -fgnu89-inline -fopenmp -ffunction-sections -fdata-sections 

CC := mpicc
OPTIONAL_FLAGS := -Wl,--print-gc-sections -Wunused-result

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
DIRPYTHONPATH := /export/home/jmdeloui/SROLL/py_sroll/
DIRPYTHONINC := /usr/include/python3.10
DIRNUMPYINC := /export/home/jmdeloui/SROLL/py_sroll/lib/python3.6/site-packages/numpy/core/include
DIRPYTHONLIB := /usr/lib/python3.10/config-3.10-x86_64-linux-gnu/
LIBPYTHONLIB := -lpython3.10
OPTIONPYTHON := -DPYTHON3
PYTHONCONF := `python-config --ldflags` `python-config --cflags`

PYTHONLIB = -L${DIRPYTHONLIB} ${LIBPYTHONLIB} -lpthread -ldl  -lutil -lm  -Xlinker -export-dynamic -Wl,-O1 -Wl,-Bsymbolic-functions

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
STIM_LIB_NAME := m
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

#CFGSROLL := -I$(FFTWDIR_INCLUDE) -I$(CFITSIO_DIR) -I$(NO_DMC_LIB_DIR_INCLUDE) -I$(DIRPYTHONINC) -I$(DIRNUMPYINC) -L$(DIRPYTHONLIB) -L$(NO_DMC_LIB_DIR) -I$(STIM_DIR_INCLUDE) -L$(STIM_LIB_DIR) -I../RngStreams -I.  $(PYTHONLIB) $(PYTHONCONF)
CFGSROLL := -I$(FFTWDIR_INCLUDE) -I$(CFITSIO_DIR) -I$(NO_DMC_LIB_DIR_INCLUDE) -I$(DIRPYTHONINC) -I$(DIRNUMPYINC) -L$(DIRPYTHONLIB) -L$(NO_DMC_LIB_DIR) -I$(STIM_DIR_INCLUDE) -I../RngStreams -I.  $(PYTHONLIB) $(PYTHONCONF)

.PHONY: all clean dist fftw cfitsio levels cleanAll alwaysCheckModules


# build hprclean before sroll, because sroll uses stimlib which needs HPR_Clean.o
all: alwaysCheckModules troll_857 troll_cfosat


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
	$(CC) -c -o $@ $< $(CFLAGS) $(CFGSROLL) -fpic -lm 
#hprclean: alwaysCheckModules
#hprclean: hprclean.c $(OBJHPRCLEAN)
#	$(CC) $(CFLAGS) $(CFGSROLL) $(OPTIONAL_FLAGS) -Wl,-R$(NO_DMC_LIB_DIR):$(STIM_LIB_DIR)  -l$(NO_DMC_LIB_NAME) -l$(STIM_LIB_NAME) -o $@ $< $(OBJ) -l$(STIM_LIB_NAME) $(FFTW_LIB) $(FFTW_OMPLIB) $(CFITSIO_LIB)  -l$(NO_DMC_LIB_NAME) -lm

sroll: alwaysCheckModules
sroll: sroll.c $(OBJ)
	$(CC) $(CFLAGS) $(CFGSROLL) $(OPTIONAL_FLAGS) -Wl,-R$(NO_DMC_LIB_DIR):$(STIM_LIB_DIR) -l$(NO_DMC_LIB_NAME) -l$(STIM_LIB_NAME) -o $@ $< $(OBJ) -l$(STIM_LIB_NAME) $(FFTW_LIB) $(FFTW_OMPLIB) $(CFITSIO_LIB)  -l$(NO_DMC_LIB_NAME) $(PYTHONLIB) -lm  

srollit: alwaysCheckModules
srollit: sroll.c $(OBJ)
	$(CC) $(CFLAGS) $(CFGSROLL) $(OPTIONAL_FLAGS) -DMAXSIMU=5 -DMAXTHEOHPR=7 -Wl,-R$(NO_DMC_LIB_DIR):$(STIM_LIB_DIR) -o $@ $< $(OBJ) -l$(STIM_LIB_NAME) $(FFTW_LIB) $(FFTW_OMPLIB) $(CFITSIO_LIB) -l$(NO_DMC_LIB_NAME) $(PYTHONLIB) -lm  

#troll: alwaysCheckModules
troll: troll.c $(OBJTROLL)
	$(CC) $(CFLAGS) $(CFGSROLL) $(OPTIONAL_FLAGS) -Wl,-R$(NO_DMC_LIB_DIR):$(STIM_LIB_DIR) -l$(NO_DMC_LIB_NAME) -l$(STIM_LIB_NAME) -o $@ $< $(OBJTROLL) -l$(STIM_LIB_NAME) $(FFTW_LIB) $(FFTW_OMPLIB) $(CFITSIO_LIB) -l$(NO_DMC_LIB_NAME) $(LIBPYTHONLIB) -DPYTHON3 $(PYTHONLIB) -lm 

#troll_14tf: alwaysCheckModules
troll_14tf: troll.c $(OBJTROLL)
	$(CC) $(CFLAGS) $(CFGSROLL) $(OPTIONAL_FLAGS)  -DMAXTHEOHPR=14 -Wl,-R$(NO_DMC_LIB_DIR):$(STIM_LIB_DIR) -l$(NO_DMC_LIB_NAME) -l$(STIM_LIB_NAME) -o $@ $< $(OBJTROLL) -l$(STIM_LIB_NAME) $(FFTW_LIB) $(FFTW_OMPLIB) $(CFITSIO_LIB) -l$(NO_DMC_LIB_NAME) $(LIBPYTHONLIB) -DPYTHON3 ${PYTHONCONF} $(PYTHONLIB) -lm 

#troll_857: alwaysCheckModules
troll_857: troll.c $(OBJTROLL)
	$(CC) $(CFLAGS) $(CFGSROLL) $(OPTIONAL_FLAGS) -DMAXTHEOHPR=7 -DMAXTHEOMAP=3 -DMAXCHAN=3 -Wl,-R$(NO_DMC_LIB_DIR):$(STIM_LIB_DIR) -l$(NO_DMC_LIB_NAME) -o $@ $< $(OBJTROLL)  $(FFTW_LIB) $(FFTW_OMPLIB) $(CFITSIO_LIB) -l$(NO_DMC_LIB_NAME) $(LIBPYTHONLIB) -DPYTHON3 ${PYTHONCONF} $(PYTHONLIB) -lm 

#troll_857: alwaysCheckModules
troll_cfosat: troll.c $(OBJTROLL)
	$(CC) $(CFLAGS) $(CFGSROLL) $(OPTIONAL_FLAGS) -DRINGSIZE=3000000 -DMAXTHEOHPR=16 -DMAXTHEOMAP=0 -DMAXCHAN=4 -Wl,-R$(NO_DMC_LIB_DIR):$(STIM_LIB_DIR) -l$(NO_DMC_LIB_NAME) -l$(STIM_LIB_NAME) -o $@ $< $(OBJTROLL) -l$(STIM_LIB_NAME) $(FFTW_LIB) $(FFTW_OMPLIB) $(CFITSIO_LIB) -l$(NO_DMC_LIB_NAME) $(LIBPYTHONLIB)  -DPYTHON3 ${PYTHONCONF} $(PYTHONLIB) -lm -fno-finite-math-only
#troll_20tf: alwaysCheckModules
troll_20tf: troll.c $(OBJTROLL)
	$(CC) $(CFLAGS) $(CFGSROLL) $(OPTIONAL_FLAGS)  -DMAXTHEOHPR=20 -Wl,-R$(NO_DMC_LIB_DIR):$(STIM_LIB_DIR) -l$(NO_DMC_LIB_NAME) -l$(STIM_LIB_NAME) -o $@ $< $(OBJTROLL) -l$(STIM_LIB_NAME) $(FFTW_LIB) $(FFTW_OMPLIB) $(CFITSIO_LIB) -l$(NO_DMC_LIB_NAME) $(LIBPYTHONLIB) -DPYTHON3 ${PYTHONCONF} $(PYTHONLIB) -lm

srollit545: alwaysCheckModules
srollit545: sroll.c $(OBJ)
	$(CC) $(CFLAGS) $(CFGSROLL) $(OPTIONAL_FLAGS) -DMAXSIMU=10 -DMAXTHEOHPR=9 -Wl,-R$(NO_DMC_LIB_DIR):$(STIM_LIB_DIR) -l$(NO_DMC_LIB_NAME) -l$(STIM_LIB_NAME) -o $@ $< $(OBJ) -l$(STIM_LIB_NAME) $(FFTW_LIB) $(FFTW_OMPLIB) $(CFITSIO_LIB) -l$(NO_DMC_LIB_NAME) ${PYTHONCONF} $(PYTHONLIB) -lm

projsroll: alwaysCheckModules
projsroll: projsroll.c $(OBJ)
	$(CC) $(CFLAGS) $(CFGSROLL) $(OPTIONAL_FLAGS) -Wl,-R$(NO_DMC_LIB_DIR):$(STIM_LIB_DIR) -l$(NO_DMC_LIB_NAME) -l$(STIM_LIB_NAME) -o $@ $< $(OBJ) $(FFTW_LIB) $(FFTW_OMPLIB) $(CFITSIO_LIB) ${PYTHONCONF} $(PYTHONLIB) -lm

#merge: alwaysCheckModules
#merge: merge.c $(OBJ)
#	$(CC) $(CFLAGS) $(CFGSROLL) $(OPTIONAL_FLAGS) -Wl,-R$(NO_DMC_LIB_DIR) -l$(NO_DMC_LIB_NAME) -o $@ $< $(OBJ) $(CFITSIO_LIB) -l$(NO_DMC_LIB_NAME) -lm ${PYTHONCONF} $(PYTHONLIB)

dist:
	$(RM) $(ARCHIVE_FILENAME)
	$(TAR) --exclude=".*" --exclude="CVS" --exclude="sroll*.tgz" --dereference -zcf $(ARCHIVE_FILENAME) *.c *.h Makefile

clean:
	$(RM) *.o projsroll sroll merge "$(PAR_LOADER_TAG_PREFIX)parLoader.c" "$(PAR_LOADER_TAG_PREFIX)parLoader.h" "$(PAR_LOADER_TAG_PREFIX_TROLL)parLoader.c" "$(PAR_LOADER_TAG_PREFIX_HPRCLEAN)parLoader.c" "$(PAR_LOADER_TAG_PREFIX_TROLL)parLoader.h"
#	$(MAKE) -C $(STIM_BASE_DIR) clean

cleanAll: clean
#	$(MAKE) -C $(FFTW_BASE_DIR)/$(FFTWDIR) clean
	$(MAKE) -C $(CFITSIO_DIR) clean
	$(MAKE) -C $(NO_DMC_LIB_BASE_DIR) clean
#	$(MAKE) -C $(PAR_LOADER_BASE_DIR) clean
#	$(MAKE) -C $(STIM_BASE_DIR) clean
	$(RM)   -rf $(FFTW_BASE_DIR)/$(FFTWDIR)
	$(RM)   -rf $(LSDIR)
