SHELL = /bin/sh


CC := gcc

NAME := no_dmc


# PATHs are automaticaly computed at installation time.

# If you do not have install the "no dmc lib" then uncomment following lines for replacement.
# Note that you may have compiled the lib before running the examples!
#IDIR := ../NO_DMC_LIB
#DIR_LIB := ../Linux-x86_64

IDIR := __IDIR__ # This path is automaticaly computed at installation time.
DIR_LIB := __DIR_LIB__ # This path is automaticaly computed at installation time.


examples: example_libno_dmc

example_libno_dmc: example_libno_dmc.c
	$(CC) -I$(IDIR) -L$(DIR_LIB) -Wl,-R$(DIR_LIB) -Wall $(CFG) -o $@ $< -l$(NAME)


clean:
	$(RM) example_libno_dmc
