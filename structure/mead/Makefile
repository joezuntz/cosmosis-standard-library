include ${COSMOSIS_SRC_DIR}/config/compilers.mk

all: libhmcode.so mead_interface.so

mead_interface.so: libhmcode.so mead_interface.f90
	$(FC) $(FFLAGS) -shared -o $@ $+ -L. -lhmcode $(LDFLAGS) -lcosmosis_fortran -lcosmosis



test: test.f90 libhmcode.so
	$(FC) $(FFLAGS) -o $@ $< -L. -lhmcode $(LDFLAGS)


clean:
	rm -f hmcode
	rm -f libhmcode.so
	rm -f hmcode.o
	rm -f mead_interface.so
	rm -rf test
	rm -rf power.dat
	rm -rf *.dSYM/
	rm -rf *.mod

libhmcode.so: HMcode.f90
	$(FC) $(FFLAGS) -c  $+ $(LDFLAGS)
	$(FC) -shared -o $(CURDIR)/$@ HMcode.o
