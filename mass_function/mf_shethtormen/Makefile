include ${COSMOSIS_SRC_DIR}/config/compilers.mk

all: st_mf_module.so

test:
	@echo "Alas,  ST has no tests"

clean:
	rm -f st_mf_module.so *.o *.mod
	rm -rf st_mf_module.so.dSYM/

st_mf_module.so: cosmo.f90 linearpk.f90 sigma.f90 spline_cubic_set.f toms446.f90  quadpack.f90  chebyshev_series.f90  chebyshev.f90  mf_shethtormen.f90 interface_tools.f90 compute_mf_shethtormen.f90  shethtormen_module.F90 ${COSMOSIS_SRC_DIR}/cosmosis/datablock/libcosmosis_fortran.so
	$(FC) $(FFLAGS) -shared -o $@ $+
