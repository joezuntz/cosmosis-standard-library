include ${COSMOSIS_SRC_DIR}/config/compilers.mk

all: hmcode.so test

test: test.f90 hmcode.so
	$(FC) $(FFLAGS) -o $@ $+ $(LDFLAGS)
	

clean:
	rm -f hmcode
	rm -rf *.dSYM/

hmcode.so: hmcode.f90
	$(FC) $(FFLAGS) -shared -o $@ $+ $(LDFLAGS)
