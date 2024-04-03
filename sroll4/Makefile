SHELL = /bin/sh
#CFLAGS := -g -O0 -Wall $(CFG) -std=c99 -fgnu89-inline -fopenmp -ffunction-sections -fdata-sections -fprofile-arcs -ftest-coverage
CFLAGS := -O3 -Wall $(CFG) -std=c99 -fgnu89-inline -fopenmp -ffunction-sections -fdata-sections 

CC := mpicc
OPTIONAL_FLAGS := -Wl,--print-gc-sections -Wunused-result

TAR := tar

ARCHIVE_FILENAME := sroll_`date +"%Y-%m-%d"`.tgz

BUILDDIR := .

SROLLDIR := ..

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

ifdef PYTHON
PYTHON := $(PYTHON)
else
PYTHON := python
endif

# --- Python wrapper ---
DIRPYTHONINC := $(shell $(PYTHON) -c "from distutils.sysconfig import get_python_inc; print(get_python_inc())")
DIRPYTHONLIB := $(shell $(PYTHON) -c "from distutils.sysconfig import get_config_var;import sysconfig; print(get_config_var('LIBDIR')+' -Wl,-rpath,'+get_config_var('LIBDIR'))")
LIBPYTHONLIB := $(shell $(PYTHON) -c "import sysconfig; print('-lpython'+sysconfig.get_config_var('LDVERSION'))")
OPTIONPYTHON := -DPYTHON3
#PYTHONCONF := $(shell $(PYTHON) -c "from distutils.sysconfig import get_python_inc; print(get_python_inc())")

PYTHONLIB = -L${DIRPYTHONLIB} ${LIBPYTHONLIB} -lpthread -ldl  -lutil -Xlinker -export-dynamic -Wl,-O1 -Wl,-Bsymbolic-functions


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

CFGSROLL := -I$(FFTWDIR_INCLUDE) -I$(CFITSIO_DIR) -I$(NO_DMC_LIB_DIR_INCLUDE) -I$(DIRPYTHONINC) -L$(DIRPYTHONLIB) -L$(NO_DMC_LIB_DIR) -I$(STIM_DIR_INCLUDE) -I../RngStreams -I.  $(PYTHONLIB) $(PYTHONCONF)

.PHONY: all clean dist fftw cfitsio levels cleanAll alwaysCheckModules


# build hprclean before sroll, because sroll uses stimlib which needs HPR_Clean.o
all: alwaysCheckModules troll_857 troll_cfosat


# DEBUG rule, to print a variable computed in this Makefile, for example OBJ, do: make print-OBJ
print-%  : ; @echo $* = $($*)


# This rule allow to ensure the proper (re)compilation of modules that sroll depend on.
# For now the relevant modules are: no_dmc_lib and stim.
alwaysCheckModules: fftw cfitsio troll_parLoader.c
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

troll_parLoader.c: $(PAR_LOADER_INPUT_HEADER_TROLL_DEF)
	$(MAKE) -C $(PAR_LOADER_BASE_DIR) currentParContentDefHeaderFile="$(PAR_LOADER_INPUT_HEADER_TROLL_DEF)" tagPrefix="$(PAR_LOADER_TAG_PREFIX_TROLL)"

%.o: %.c alwaysCheckModules
	$(CC) -c -o $@ $< $(CFLAGS) $(CFGSROLL) -fpic

#troll_857: alwaysCheckModules
troll_857: troll.c $(OBJTROLL)
	$(CC) $(CFLAGS) $(CFGSROLL) $(OPTIONAL_FLAGS) -DMAXEXTERNALMAP=3 -DMAXCHAN=3 -Wl,-R$(NO_DMC_LIB_DIR):$(STIM_LIB_DIR) -l$(NO_DMC_LIB_NAME) -o $@ $< $(OBJTROLL)  $(FFTW_LIB) $(FFTW_OMPLIB) $(CFITSIO_LIB) -l$(NO_DMC_LIB_NAME) $(LIBPYTHONLIB) -DPYTHON3 $(PYTHONLIB) -lm 

#troll_857: alwaysCheckModules
troll_cfosat_large: troll.c $(OBJTROLL)
	$(CC) $(CFLAGS) $(CFGSROLL) $(OPTIONAL_FLAGS) -DRINGSIZE=3000000 -DMAXSIMU=1 -DMAXEXTERNALSHPR=48  -DMAXEXTERNALMAP=1 -DMAXCHAN=2 -Wl,-R$(NO_DMC_LIB_DIR):$(STIM_LIB_DIR) -l$(NO_DMC_LIB_NAME) -l$(STIM_LIB_NAME) -o $@ $< $(OBJTROLL) -l$(STIM_LIB_NAME) $(FFTW_LIB) $(FFTW_OMPLIB) $(CFITSIO_LIB) -l$(NO_DMC_LIB_NAME) $(LIBPYTHONLIB)  -DPYTHON3 $(PYTHONLIB) -fno-finite-math-only -lm

#troll_857: alwaysCheckModules
troll_cfosat: troll.c $(OBJTROLL)
	$(CC) $(CFLAGS) $(CFGSROLL) $(OPTIONAL_FLAGS) -DRINGSIZE=3000000 -DMAXSIMU=1 -DMAXEXTERNALSHPR=16  -DMAXEXTERNALMAP=1 -DMAXCHAN=2 -Wl,-R$(NO_DMC_LIB_DIR):$(STIM_LIB_DIR) -l$(NO_DMC_LIB_NAME) -l$(STIM_LIB_NAME) -o $@ $< $(OBJTROLL) -l$(STIM_LIB_NAME) $(FFTW_LIB) $(FFTW_OMPLIB) $(CFITSIO_LIB) -l$(NO_DMC_LIB_NAME) $(LIBPYTHONLIB)  -DPYTHON3 $(PYTHONLIB) -fno-finite-math-only -lm

#troll_857: alwaysCheckModules
troll_planck: troll.c $(OBJTROLL)
	$(CC) $(CFLAGS) $(CFGSROLL) $(OPTIONAL_FLAGS) -DRINGSIZE=27664 -DMAXSIMU=1 -DMAXEXTERNALSHPR=16  -DMAXEXTERNALMAP=1 -DMAXCHAN=3 -Wl,-R$(NO_DMC_LIB_DIR):$(STIM_LIB_DIR) -l$(NO_DMC_LIB_NAME) -l$(STIM_LIB_NAME) -o $@ $< $(OBJTROLL) -l$(STIM_LIB_NAME) $(FFTW_LIB) $(FFTW_OMPLIB) $(CFITSIO_LIB) -l$(NO_DMC_LIB_NAME) $(LIBPYTHONLIB)  -DPYTHON3 $(PYTHONLIB) -fno-finite-math-only -lm
#troll_20tf: alwaysCheckModules

#troll_857: alwaysCheckModules
troll_LB: troll.c $(OBJTROLL)
	$(CC) $(CFLAGS) $(CFGSROLL) $(OPTIONAL_FLAGS) -DRINGSIZE=37200 -DMAXSIMU=1 -DMAXEXTERNALSHPR=6  -DMAXEXTERNALMAP=1 -DMAXCHAN=3 -Wl,-R$(NO_DMC_LIB_DIR):$(STIM_LIB_DIR) -o $@ $< $(OBJTROLL) -l$(NO_DMC_LIB_NAME) -l$(STIM_LIB_NAME) -l$(STIM_LIB_NAME) $(FFTW_LIB) $(FFTW_OMPLIB) $(CFITSIO_LIB) -l$(NO_DMC_LIB_NAME) $(LIBPYTHONLIB)  -DPYTHON3 $(PYTHONLIB) -fno-finite-math-only -lm

dist:
	$(RM) $(ARCHIVE_FILENAME)
	$(TAR) --exclude=".*" --exclude="CVS" --exclude="sroll*.tgz" --dereference -zcf $(ARCHIVE_FILENAME) *.c *.h Makefile

clean:
	$(RM) *.o "$(PAR_LOADER_TAG_PREFIX)parLoader.c" "$(PAR_LOADER_TAG_PREFIX)parLoader.h" "$(PAR_LOADER_TAG_PREFIX_TROLL)parLoader.c" "$(PAR_LOADER_TAG_PREFIX_HPRCLEAN)parLoader.c" "$(PAR_LOADER_TAG_PREFIX_TROLL)parLoader.h"
#	$(MAKE) -C $(STIM_BASE_DIR) clean

cleanAll: clean
	$(MAKE) -C $(FFTW_BASE_DIR)/$(FFTWDIR) clean
	$(MAKE) -C $(CFITSIO_DIR) clean
	$(MAKE) -C $(NO_DMC_LIB_BASE_DIR) clean
	$(MAKE) -C $(PAR_LOADER_BASE_DIR) clean
	$(MAKE) -C $(STIM_BASE_DIR) clean
	$(RM)   -rf $(FFTW_BASE_DIR)/$(FFTWDIR)
	$(RM)   -rf $(LSDIR)
