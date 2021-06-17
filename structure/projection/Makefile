include ${COSMOSIS_SRC_DIR}/config/compilers.mk

#USER_CFLAGS= -I ../limber -I ${GSL_INC}
#USER_LDFLAGS= -L../limber -llimber -lcosmosis -L ${GSL_LIB} -lgsl -lgslcblas 
USER_CFLAGS= -I ${GSL_INC}  -O2
USER_LDFLAGS= -lcosmosis -L ${GSL_LIB} -lgsl -lgslcblas 

all: src/spec_tools.so

src/spec_tools.so: src/kernel.c src/limber.c src/non_limber.c src/utils.c src/interp2d.c 
	$(CC) $(CFLAGS)  -shared -o $(CURDIR)/src/spec_tools.so src/kernel.c src/non_limber.c src/limber.c src/utils.c src/interp2d.c $(LDFLAGS)

clean:
	rm -rf src/spec_tools.so src/spec_tools.so.dSYM

test:
	@echo Alas, shear/spectra has no tests

