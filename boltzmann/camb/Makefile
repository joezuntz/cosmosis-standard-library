include ${COSMOSIS_SRC_DIR}/config/compilers.mk

USER_FFLAGS=-I.
USER_LDFLAGS=-lcosmosis_fortran
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

interface: $(CAMBLIB) camb_interface.F90
	$(F90C) $(F90FLAGS) -shared camb_interface.F90 camb_mode.F90 -fPIC -o camb.so $(LDFLAGS)   -L. -lcamb_$(RECOMBINATION)

include ./Makefile_main

test:
	@echo "Alas, camb has no tests"
