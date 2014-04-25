include ${COSMOSIS_SRC_DIR}/config/compilers.mk

all: halofit_module.so

test:
	@echo "Alas, halofit has no tests"

clean:
	rm -f halofit_module.so *.o *.mod

halofit_module.so: halofit.F90 halofit_module.F90 ${COSMOSIS_SRC_DIR}/cosmosis/datablock/libcosmosis_fortran.so
	$(FC) $(FFLAGS) -shared -o $@ $+
