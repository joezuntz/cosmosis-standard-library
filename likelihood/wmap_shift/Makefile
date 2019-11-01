include ${COSMOSIS_SRC_DIR}/config/compilers.mk

USER_LDFLAGS=-lcosmosis_fortran

all: wmap_shift.so

wmap_shift.so: wmap_shift.f90
	$(FC) $(FFLAGS) -shared -o wmap_shift.so wmap_shift.f90 $(LDFLAGS) 

clean:
	rm -f *.o *.mod wmap_shift.so
	rm -rf wmap_shift.so.dSYM

test:
	@echo Alas, wmap_shift has no tests

