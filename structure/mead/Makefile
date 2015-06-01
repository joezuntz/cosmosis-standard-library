include ${COSMOSIS_SRC_DIR}/config/compilers.mk

all: hmcode

test:
	@echo "Alas, CRL ST has no tests"

clean:
	rm -f hmcode
	rm -rf *.dSYM/

hmcode: hmcode.f90
	$(FC) $(FFLAGS) -o $@ $+
