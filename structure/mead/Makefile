include ${COSMOSIS_SRC_DIR}/config/compilers.mk

all: hmcode.so mead_interface.so

mead_interface.so: hmcode.so mead_interface.f90
	$(FC) $(FFLAGS) -shared -o $@ $+ $(LDFLAGS) -lcosmosis_fortran -lcosmosis

test: test.f90 hmcode.so
	$(FC) $(FFLAGS) -o $@ $+ $(LDFLAGS)
	

clean:
	rm -f hmcode
	rm -f hmcode.so
	rm -f mead_interface.so
	rm -rf test
	rm -rf power.dat
	rm -rf *.dSYM/
	rm -rf *.mod

hmcode.so: hmcode.f90
	$(FC) $(FFLAGS) -shared -o $@ $+ $(LDFLAGS)
