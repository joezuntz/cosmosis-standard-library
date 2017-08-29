include ${COSMOSIS_SRC_DIR}/config/compilers.mk

USER_FFLAGS=-I.
USER_LDFLAGS=-lcosmosis_fortran


INTERFACE=planck2015_lowl.so

all: $(INTERFACE)

$(INTERFACE): comm_gauss_br_mod.f90 planck2015_lowl.f90
	$(FC) $(FFLAGS) -c comm_gauss_br_mod.f90 
	$(FC) $(FFLAGS) -shared planck2015_lowl.f90 comm_gauss_br_mod.o -o $(INTERFACE) -L${CFITSIO_LIB} -lcfitsio ${LAPACK_LINK} $(LDFLAGS) 

clean:
	rm -rf *.so *.o *.mod *.dSYM

test:
	@echo "Alas, planck_2015_lowl has no tests"
