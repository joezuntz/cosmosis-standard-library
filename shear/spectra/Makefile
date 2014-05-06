include ${COSMOSIS_SRC_DIR}/config/compilers.mk

USER_CFLAGS= -I ../limber -I ${GSL_INC}
USER_LDFLAGS=-lcosmosis -L ${GSL_LIB} -lgsl -L../limber -llimber

all: interface.so

interface.so: utils.c shear_shear.c interface.c ../limber/interp2d.o  ../limber/limber.o
	$(CC) $(CFLAGS) $(LDFLAGS) -shared -o interface.so utils.c shear_shear.c interface.c ../limber/interp2d.o  ../limber/limber.o

clean:
	rm -rf interface.so interface.so.dSYM
