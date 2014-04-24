include ${COSMOSIS_DIR}/config/compilers.mk

all: halofit_module.so

clean:
	rm -f halofit_module.so *.o *.mod

halofit_module.so: halofit.F90 halofit_module.F90 ${COSMOSIS_DIR}/cosmosis/datablock/libcosmosis_fortran.so
	$(FC) $(FFLAGS) -shared -o $@ $+
