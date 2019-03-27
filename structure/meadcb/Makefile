include ${COSMOSIS_SRC_DIR}/config/compilers.mk

all: libmhmcamb.a mead_interface.so

mead_interface.so: libmhmcamb.a mead_interface.f90
	$(FC) $(FFLAGS) -shared -o $@ $+ -L. -lmhmcamb $(LDFLAGS) -lcosmosis_fortran -lcosmosis



test: test.f90 libmhmcamb.a
	$(FC) $(FFLAGS) -o $@ $< -L. -lmhmcamb $(LDFLAGS)


clean:
	rm -f mhmcamb
	rm -f libmhmcamb.a
	rm -f mhmcamb.o
	rm -f mead_interface.so
	rm -rf test
	rm -rf power.dat
	rm -rf *.dSYM/
	rm -rf *.mod

libmhmcamb.a: mhmcamb.f90
	$(FC) $(FFLAGS) -c  $+ $(LDFLAGS)
	$(AR) rc $@ mhmcamb.o
