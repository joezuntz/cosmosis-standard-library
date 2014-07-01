include ${COSMOSIS_SRC_DIR}/config/compilers.mk

all: nowiggle_module.so

test:
	@echo "Alas, CRL ST has no tests"

clean:
	rm -f nowiggle_module.so *.o *.mod
	rm -rf nowiggle_module.so.dSYM/

nowiggle_module.so: cosmo.f90  spline.f hunt.f interface_tools.f90 compute_pk_nowiggle.f90 eisensteinhu.f90 nowiggle_module.F90 ${COSMOSIS_SRC_DIR}/cosmosis/datablock/libcosmosis_fortran.so
	$(FC) $(FFLAGS) -shared -o $@ $+
