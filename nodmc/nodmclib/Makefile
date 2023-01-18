SHELL = /bin/sh

# Ensure 'mpicc' (ex on magique4) is available otherwise use 'cc' (ex on edison at NERSC)
ifeq "$(shell which mpicc 2> /dev/null)" ""
CC := cc
OPTIONAL_FLAGS :=
else
CC := mpicc
OPTIONAL_FLAGS :=
endif

INSTALL := install
MKDIR_P := mkdir -p
dir_guard=@mkdir -p $(@D)

NAME := no_dmc

BUILDDIR := Linux-x86_64
DIR_LIB := $(BUILDDIR)
DIR_TEST_SRC := test
DIR_EXAMPLES_SRC := examples
DIR_HEADER := NO_DMC_LIB
DIR_SRC := src
DIR_BIN := $(BUILDDIR)
DIR_TOOLS := tools
PY_DIR := pynodmclib
VPATH := $(DIR_SRC) $(DIR_HEADER)

#CFLAGS := -O3 -Wall -Werror -fpic -I$(DIR_HEADER) -I$(SROLLDIR)/cfitsio/cfitsio $(CFG) -std=c99
CFLAGS := -O3 -Wall -fpic -I$(DIR_HEADER) -I$(SROLLDIR)/cfitsio/cfitsio $(CFG) -std=c99

PREFIX?=/usr
exec_prefix=$(PREFIX)
LIBDIR ?= $(exec_prefix)/lib64
INCLUDEDIR ?= $(PREFIX)/include
BINDIR ?= $(exec_prefix)/bin

# This directory is only used when user define a DESTDIR
pythonlibdir := $(exec_prefix)/python/site-packages


# MAJOR and MINOR are retrieve from C source file so that there is a unique place for the values.
# These variable are used in Makefile for distribution and installation rules only.
version_source_file := $(DIR_HEADER)/no_dmc_version.h
MAJOR_Identifier := __NODMCLIB_MAJOR__
MINOR_Identifier := __NODMCLIB_MINOR__
MAJOR := `cat $(version_source_file) | tr -s ' ' | grep "\#define $(MAJOR_Identifier)" | cut -d" " -f 3`
MINOR := `cat $(version_source_file) | tr -s ' ' | grep "\#define $(MINOR_Identifier)" | cut -d" " -f 3`

VERSION := $(MAJOR).$(MINOR)

RELEASE_DIR := $(shell echo "$(NAME)-$(VERSION)")
RELEASE_TAR_FILE := $(RELEASE_DIR).tgz

IDIR := $(DIR_HEADER)
_DEPS := no_dmc_piolib_type_def.h \
         no_dmc_util.h \
         no_dmc_metadata.h \
         no_dmc_data_access.h \
         no_dmc_debug.h \
         no_dmc_version.h \
         ring_index.h \
         brimo.h \
         binary_file_toi.h
DEPS := $(patsubst %,$(IDIR)/%,$(_DEPS))

ODIR := $(BUILDDIR)
_OBJ := no_dmc_metadata.o \
        no_dmc_data_access.o \
        no_dmc_util.o \
        no_dmc_version.o \
        ring_index.o \
        brimo.o \
        binary_file_toi.o
OBJ := $(patsubst %,$(ODIR)/%,$(_OBJ))

# --- CFITSIO package ---
CFITSIO_ROOT := $(SROLLDIR)/cfitsio
CFITSIO_TAR  := $(CFITSIO_ROOT)/cfitsio3410.tar.gz
CFITSIO_DIR  := $(CFITSIO_ROOT)/cfitsio
CFITSIO_LIB  := $(CFITSIO_DIR)/libcfitsio.a

### PYTHON version check and update env var ###

PYTHON := "$(shell which python 2> /dev/null)"
export PYTHON

ifdef DESTDIR
	destdir_examples := "$(DESTDIR)/$(DIR_EXAMPLES_SRC)"
else
	destdir_examples := "/tmp/libno_dmc/$(DIR_EXAMPLES_SRC)"
endif
export destdir_examples



.PHONY: all lib test pylib dist clean examples

all: lib pylib

check_cmd:
ifndef PYTHON
	@echo "Makefile requires 'python' version 2.7.5 or 2.7.9"
	@exit 1
else
	@echo "python version: OK"
endif

cfitsio: $(CFITSIO_LIB)

$(CFITSIO_LIB) :
	cd $(CFITSIO_ROOT) && $(TAR) zxf $(CFITSIO_TAR)
	cd $(CFITSIO_DIR) && ./configure
	$(MAKE) -C $(CFITSIO_DIR) all shared


#lib: $(DIR_LIB)/lib$(NAME).so  $(DIR_LIB)/lib$(NAME).a
lib: $(DIR_LIB)/lib$(NAME).so

test: lib $(DIR_BIN)/test_libno_dmc

$(DIR_BIN)/test_libno_dmc: $(DIR_TEST_SRC)/test_no_dmc_lib.c
	$(dir_guard)
	$(CC) -I$(IDIR) -L$(DIR_LIB) -Wall $(CFG) -o $@ $< -l$(NAME)

$(BUILDDIR)/lib$(NAME).so: $(OBJ)
	$(dir_guard)
	$(CC) -shared -Wl,--whole-archive $(CFITSIO_LIB) -Wl,--no-whole-archive -o $@ $^

#$(BUILDDIR)/lib$(NAME).a: $(OBJ)
#	$(dir_guard)
#	cd $(BUILDDIR) && $(AR) -cvq lib$(NAME).a $(OBJ)

$(ODIR)/%.o: %.c $(DEPS)
	$(dir_guard)
	$(CC) -c -o $@ $< $(CFLAGS)

# This is the rule for generating the python wrapper library
pylib: check_cmd lib $(BUILDDIR)/nodmclib.so

$(BUILDDIR)/nodmclib.so: $(PY_DIR)/setup.py $(PY_DIR)/nodmclib/cno_dmc_data_access.pxd $(PY_DIR)/nodmclib/nodmclib.pyx
	@echo "PYTHON = $(PYTHON)"
	@echo "PYTHONPATH = $(PYTHONPATH)"
	$(dir_guard)
	$(RM) $(BUILDDIR)/nodmclib.so*
	cd $(PY_DIR) && python setup.py build_ext --build-lib ../$(BUILDDIR)

# Allow to generate the example prog
examples: lib pylib $(DIR_BIN)/example_libno_dmc

$(DIR_BIN)/example_libno_dmc: $(DIR_EXAMPLES_SRC)/example_libno_dmc.c
	$(dir_guard)
	$(CC) -I$(IDIR) -L$(DIR_LIB) -Wall $(CFG) -o $@ $< -l$(NAME)

pylib-install: check_cmd lib
ifdef DESTDIR
	cd $(PY_DIR) && python setup.py install --install-lib $(DESTDIR)$(pythonlibdir)
else
	cd $(PY_DIR) && python setup.py install
endif


# Allow to install the libraries and the tools in appropriate directories
install: lib pylib-install
	$(MKDIR_P) "$(DESTDIR)$(INCLUDEDIR)" || exit 1;
	$(INSTALL) -m 644 $(wildcard $(DIR_HEADER)/*.h) "$(DESTDIR)$(INCLUDEDIR)"
	$(MKDIR_P) "$(DESTDIR)$(LIBDIR)" || exit 1;
	$(INSTALL) -m 755 "$(DIR_LIB)/lib$(NAME).so" "$(DESTDIR)$(LIBDIR)"
	mv "$(DESTDIR)$(LIBDIR)/lib$(NAME).so" "$(DESTDIR)$(LIBDIR)/lib$(NAME).so.$(VERSION)"
	ln -s lib$(NAME).so.$(VERSION) "$(DESTDIR)$(LIBDIR)/lib$(NAME).so"
	$(MKDIR_P) "$(DESTDIR)$(BINDIR)" || exit 1;
	$(INSTALL) -m 755 "$(DIR_TOOLS)/no_dmc_metadata.py" "$(DESTDIR)$(BINDIR)"
	$(INSTALL) -m 755 "$(DIR_TOOLS)/no_dmc_metadata_fix_backendname.py" "$(DESTDIR)$(BINDIR)"
	$(MKDIR_P) $(destdir_examples) || exit 1;
	$(INSTALL) -m 644 "$(DIR_EXAMPLES_SRC)/example_libno_dmc.c" "$(DIR_EXAMPLES_SRC)/example_no_dmc_python.py" "$(DIR_EXAMPLES_SRC)/Makefile" "$(DIR_EXAMPLES_SRC)/README" $(destdir_examples)
	sed -i "s:__IDIR__:$(DESTDIR)$(INCLUDEDIR):g" $(destdir_examples)/Makefile
	sed -i "s:__DIR_LIB__:$(DESTDIR)$(LIBDIR):g" $(destdir_examples)/Makefile
	@echo -e "In order to be able to execute the examples you may set env variables as folowing:\n  export LD_LIBRARY_PATH=$(DESTDIR)$(LIBDIR)\n  export PYTHONPATH=$(DESTDIR)$(pythonlibdir)/nodmclib\n\nRequired python version: $(PYTHON)" > $(destdir_examples)/README.FIRST
#	Note: the folowing command print info about further linking to the installed library (ex: ldconfig)
	libtool --mode=finish "$(DESTDIR)$(LIBDIR)"


# This rule is used for generting the deliverable
# Test if there is already a directory named RELEASE_DIR: if true print warn message and do nothing, otherwise construct the dist
dist:
ifneq ($(wildcard $(RELEASE_DIR)), )
	@echo "*************************************************************"
	@echo "WARNING! Directory \"$(RELEASE_DIR)/\" already exist, please 'make clean' before attempting to construct dist!"
	@echo "*************************************************************"
else
	@echo "*************************************************************"
	@echo "OK we will construct the dist in \"$(RELEASE_DIR)/\"..."
	@echo "--> dist tar file will be \"$(RELEASE_TAR_FILE)\""
	@echo "*************************************************************"
	$(MKDIR_P) $(RELEASE_DIR)
	cp Makefile AUTHORS COPYING INSTALL README FAQ RELEASE-NOTES $(RELEASE_DIR)
	$(MKDIR_P) $(RELEASE_DIR)/doc
	cp doc/README $(RELEASE_DIR)/doc
	$(MKDIR_P) $(RELEASE_DIR)/$(DIR_HEADER)
	cp $(wildcard $(DIR_HEADER)/*.h) $(RELEASE_DIR)/$(DIR_HEADER)
	$(MKDIR_P) $(RELEASE_DIR)/$(DIR_SRC)
	cp $(wildcard $(DIR_SRC)/*.c) $(RELEASE_DIR)/$(DIR_SRC)
	$(MKDIR_P) $(RELEASE_DIR)/$(DIR_TOOLS)
	cp $(DIR_TOOLS)/no_dmc_metadata.py $(DIR_TOOLS)/no_dmc_metadata_fix_backendname.py $(DIR_TOOLS)/README $(RELEASE_DIR)/$(DIR_TOOLS)
	$(MKDIR_P) $(RELEASE_DIR)/$(PY_DIR)
	cp $(PY_DIR)/setup.py $(PY_DIR)/README $(RELEASE_DIR)/$(PY_DIR)
	$(MKDIR_P) $(RELEASE_DIR)/$(PY_DIR)/nodmclib
	cp $(PY_DIR)/nodmclib/cno_dmc_data_access.pxd $(PY_DIR)/nodmclib/nodmclib.pyx $(PY_DIR)/nodmclib/__init__.py $(RELEASE_DIR)/$(PY_DIR)/nodmclib
	$(MKDIR_P) $(RELEASE_DIR)/$(DIR_EXAMPLES_SRC)
	cp $(wildcard $(DIR_EXAMPLES_SRC)/*.c $(DIR_EXAMPLES_SRC)/*.py) $(DIR_EXAMPLES_SRC)/Makefile $(DIR_EXAMPLES_SRC)/README $(RELEASE_DIR)/$(DIR_EXAMPLES_SRC)
	$(MKDIR_P) $(RELEASE_DIR)/$(DIR_TEST_SRC)
	cp $(DIR_TEST_SRC)/test_no_dmc_lib.c $(DIR_TEST_SRC)/compare_dmc.py $(DIR_TEST_SRC)/README $(RELEASE_DIR)/$(DIR_TEST_SRC)
	tar zcf $(RELEASE_TAR_FILE) $(RELEASE_DIR)
endif

clean:
	$(RM) $(BUILDDIR)/*.o $(BUILDDIR)/lib$(NAME).so* $(BUILDDIR)/nodmclib.so* $(BUILDDIR)/lib$(NAME).a
	$(RM) -r $(BUILDDIR)/nodmclib
	$(RM) -r $(PY_DIR)/build $(PY_DIR)/nodmclib.so $(PY_DIR)/nodmclib/nodmclib.c
	$(RM) $(DIR_BIN)/test_libno_dmc $(DIR_BIN)/example_libno_dmc
	$(RM) -r $(RELEASE_DIR)
	$(RM) $(DIR_EXAMPLES_SRC)/example_libno_dmc
	$(MAKE) -C $(CFITSIO_DIR) clean
