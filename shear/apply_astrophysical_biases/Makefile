include ${COSMOSIS_SRC_DIR}/config/compilers.mk

USER_CFLAGS= -I ../limber -I ${GSL_INC}
USER_LDFLAGS= -L../limber -llimber -lcosmosis -L ${GSL_LIB} -lgsl -lgslcblas 

all: interface.so

interface.so: interface.c 
	$(CC) $(CFLAGS)  -shared -o interface.so interface.c $(LDFLAGS)

clean:
	rm -rf interface.so interface.so.dSYM

test:
	@echo Alas, shear/apply_astrophysical_biases has no tests

