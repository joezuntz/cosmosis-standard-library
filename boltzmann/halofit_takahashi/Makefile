include ${COSMOSIS_SRC_DIR}/config/compilers.mk

all: halofit_interface.so

USER_LDFLAGS=-L${COSMOSIS_SRC_DIR}/cosmosis/datablock -lcosmosis_fortran

halofit_interface.so: camb_params_stub.f90 matter_power_stub.f90 nonlinear.f90 halofit_interface.F90 ${COSMOSIS_SRC_DIR}/cosmosis/datablock/libcosmosis_fortran.so
	$(FC) $(FFLAGS) -fPIC -c -o camb_params_stub.o camb_params_stub.f90
	$(FC) $(FFLAGS) -fPIC -c -o matter_power_stub.o matter_power_stub.f90
	$(FC) $(FFLAGS) -fPIC -c -o nonlinear.o nonlinear.f90
	$(FC) $(FFLAGS)  -fPIC -shared -o halofit_interface.so halofit_interface.F90 camb_params_stub.o matter_power_stub.o nonlinear.o $(LDFLAGS)



clean:
	rm -rf *.o *~ *.flc smith2demo fsmith2  halofit_interface.so halofit_interface.so.dSYM *.mod

test:
	@echo "Alas, halofit_takahashi has no tests"