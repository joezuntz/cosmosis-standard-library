include ${COSMOSIS_SRC_DIR}/config/compilers.mk

USER_CFLAGS= -I  ${GSL_INC}
USER_LDFLAGS=  -lcosmosis -L ${GSL_LIB} -lgsl -lgslcblas

all: interface.so

test:
	@echo Alas, growth_factor has no tests

interface.so: growthfactor.c qromb.c polint.c trapzd.c nrutil.c qromo.c midpnt.c  interface.c
	$(CC) $(CFLAGS)  -shared -o interface.so growthfactor.c qromb.c polint.c trapzd.c nrutil.c qromo.c midpnt.c interface.c $(LDFLAGS)

clean:
	rm -rf interface.so interface.so.dSYM  *.o
