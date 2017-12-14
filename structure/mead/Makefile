include ${COSMOSIS_SRC_DIR}/cosmosis/compilers.mk

all: libhmcode.a mead_interface.so

mead_interface.so: libhmcode.a mead_interface.f90
	$(FC) $(FFLAGS) -shared -o $@ $+ -L. -lhmcode $(LDFLAGS)  -lcosmosis



test: test.f90 libhmcode.a
	$(FC) $(FFLAGS) -o $@ $< -L. -lhmcode $(LDFLAGS)


clean:
	rm -f hmcode
	rm -f libhmcode.a
	rm -f hmcode.o
	rm -f mead_interface.so
	rm -rf test
	rm -rf power.dat
	rm -rf *.dSYM/
	rm -rf *.mod

libhmcode.a: HMcode.f90
	$(FC) $(FFLAGS) -c  $+ $(LDFLAGS)
	$(AR) rc $@ HMcode.o
