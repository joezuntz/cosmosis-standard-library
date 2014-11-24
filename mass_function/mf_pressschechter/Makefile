include ${COSMOSIS_SRC_DIR}/config/compilers.mk

all: ps_mf_module.so

test:
	@echo "Alas,  PS has no tests"

clean:
	rm -f ps_mf_module.so *.o *.mod
	rm -rf ps_mf_module.so.dSYM/

ps_mf_module.so: cosmo.f90 linearpk.f90 sigma.f90 spline_cubic_set.f toms446.f90  quadpack.f90  chebyshev_series.f90  chebyshev.f90  mf_pressschechter.f90 interface_tools.f90 compute_mf_pressschechter.f90  pressschechter_module.F90 ${COSMOSIS_SRC_DIR}/cosmosis/datablock/libcosmosis_fortran.so
	$(FC) $(FFLAGS) -shared -o $@ $+
