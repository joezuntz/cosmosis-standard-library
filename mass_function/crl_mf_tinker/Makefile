include ${COSMOSIS_SRC_DIR}/config/compilers.mk

all: tinker_mf_module.so

test:
	@echo "Alas, CRL Tinker has no tests"

clean:
	rm -f tinker_mf_module.so *.o *.mod
	rm -rf tinker_mf_module.so.dSYM/

tinker_mf_module.so: cosmo.f90 linearpk.f90 sigma.f90 spline.f hunt.f qromb.f chebyshev.f growth.f90 dverk.f mf_tinker.f90 interface_tools.f90 compute_mf_tinker.f90  tinker_module.F90 ${COSMOSIS_SRC_DIR}/cosmosis/datablock/libcosmosis_fortran.so
	$(FC) $(FFLAGS) -shared -o $@ $+
