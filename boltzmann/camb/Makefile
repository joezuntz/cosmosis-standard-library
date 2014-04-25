include ${COSMOSIS_DIR}/config/compilers.mk

USER_LDFLAGS = -lcosmosis_fortran

F90C=$(FC)

#Set FISHER=Y to compile bispectrum fisher matrix code
FISHER ?=

ifneq ($(FISHER),)
FFLAGS += -DFISHER
EXTCAMBFILES = Matrix_utils.o
else
EXTCAMBFILES =
endif
EQUATIONS=equations_ppf
FITSLIB=cfitsio

include ./Makefile_main
