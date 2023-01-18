SHELL = /bin/sh

# Ensure 'mpicc' (ex on magique4) is available otherwise use 'cc' (ex on edison at NERSC)
ifeq "$(shell which mpicc 2> /dev/null)" ""
CC := cc
OPTIONAL_FLAGS :=
else
CC := mpicc
OPTIONAL_FLAGS := -lm
endif

# Define some default commands
MKDIR_P := mkdir -p
PYTHON := python

# Here specify the default paths corresponding to NO_DMC_LIB installation
NODMC_IDIR := $(shell pwd)/../nodmclib/NO_DMC_LIB
DIR_LIB := $(shell pwd)/../nodmclib/Linux-x86_64/
NO_DMC_LIB_NAME := no_dmc

# Release info
NAME := parLoader
MAJOR := 1
MINOR := 3
VERSION := $(MAJOR).$(MINOR)
RELEASE_DIR := "$(NAME)-$(VERSION)"
RELEASE_TAR_FILE := $(RELEASE_DIR).tgz

# TEST env
TEST_DIR := test
TEST_PARLOADER_TAGPREFIX := test_
OBJ_TEST := $(TEST_DIR)/$(TEST_PARLOADER_TAGPREFIX)parLoader.o


CFLAGS := -O3 -Wall -Werror -I$(NODMC_IDIR) -L$(DIR_LIB) $(CFG) -std=c99


.PHONY: all parLoader test clean dist check-env

all: parLoader

parLoader: check-env fillParLoader.py parLoader.c_GEN parLoader.h_GEN
	$(PYTHON) fillParLoader.py "$(currentParContentDefHeaderFile)" "$(tagPrefix)"

%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS)

test: $(TEST_DIR)/testParLoader

$(TEST_DIR)/$(TEST_PARLOADER_TAGPREFIX)parLoader.c: fillParLoader.py parLoader.c_GEN parLoader.h_GEN
	$(PYTHON) fillParLoader.py "$(TEST_DIR)/test_param.h" "$(TEST_PARLOADER_TAGPREFIX)"

$(TEST_DIR)/testParLoader: $(TEST_DIR)/testParLoader.c $(TEST_DIR)/$(TEST_PARLOADER_TAGPREFIX)parLoader.c $(OBJ_TEST)
	$(CC) $(CFLAGS) -Wl,-R$(DIR_LIB) -l$(NO_DMC_LIB_NAME) -o $@ $< $(OBJ_TEST)

#testMPI: $(TEST_DIR)/testParLoader_MPI
#
#$(TEST_DIR)/testParLoader_MPI: $(TEST_DIR)/testParLoader_MPI.c $(OBJ)
#	$(CC) $(CFLAGS) -Wl,-R$(DIR_LIB) -o $@ $< $(OBJ) -l$(NO_DMC_LIB_NAME)

### Ensure appropriate vars are set! ###
# currentParContentDefHeaderFile: the header describing the content structure.
# prefix: the name to be used in front of file and function name! Allow to
#         avoid collision between several 'instance' of the parloader.
check-env:
ifndef currentParContentDefHeaderFile
	$(error ERROR: 'currentParContentDefHeaderFile' must be defined!)
endif
ifndef tagPrefix
	$(error ERROR: 'tagPrefix' must be defined!)
endif
	@echo "ParLoader check-env --> OK"

dist:
ifneq ($(wildcard $(RELEASE_DIR)), )
	@echo "WARNING! Directory \"$(RELEASE_DIR)/\" already exist, please 'make clean' before attempting to construct dist!"
else
	@echo "OK we will construct the dist in \"$(RELEASE_DIR)/\"..."
	@echo "--> dist tar file will be \"$(RELEASE_TAR_FILE)\""
	$(MKDIR_P) $(RELEASE_DIR)
	cp README RELEASE-NOTES $(RELEASE_DIR)
	cp fillParLoader.py parLoader.h_GEN parLoader.c_GEN Makefile $(RELEASE_DIR)
	$(MKDIR_P) $(RELEASE_DIR)/$(TEST_DIR)
	cp $(TEST_DIR)/test_param.h $(TEST_DIR)/testParLoader.c $(TEST_DIR)/testParLoader_MPI.c $(RELEASE_DIR)/$(TEST_DIR)/
	$(MKDIR_P) $(RELEASE_DIR)/$(TEST_DIR)/samples/
	cp $(TEST_DIR)/samples/Param*.txt $(RELEASE_DIR)/$(TEST_DIR)/samples/
	tar zcf $(RELEASE_TAR_FILE) $(RELEASE_DIR)
endif

clean:
	$(RM) *.o $(TEST_DIR)/*.o $(TEST_DIR)/testParLoader $(TEST_DIR)/testParLoader_MPI $(TEST_DIR)/$(TEST_PARLOADER_TAGPREFIX)parLoader.c $(TEST_DIR)/$(TEST_PARLOADER_TAGPREFIX)parLoader.h
	$(RM) -r $(RELEASE_DIR)
