include ../../../config/compilers.mk
include ../../../config/subdirs.mk

SUBDIRS = 

.PHONY: all

all: halofit_module.so
	$(process_subdirs)

test:
	$(process_subdirs)

clean:
	@echo Cleaning in $$(PWD)
	rm -f halofilt_module.so *.o *.mod
	$(process_subdirs)

halofit_module.so: halofit.F90 halofit_module.F90 ../../../cosmosis/datablock/libcosmosis_fortran.so
	$(FC) $(FFLAGS) -J../../../cosmosis/datablock -shared -o $@ $+
