include ${COSMOSIS_DIR}/config/compilers.mk

SUBDIRS = 

all: halofit_module.so

clean:
	rm -f halofilt_module.so *.o *.mod
	$(process_subdirs)

halofit_module.so: halofit.F90 halofit_module.F90 ../../../cosmosis/datablock/libcosmosis_fortran.so
	$(FC) $(FFLAGS) -J../../../cosmosis/datablock -shared -o $@ $+
