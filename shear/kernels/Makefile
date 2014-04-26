include ${COSMOSIS_SRC_DIR}/config/compilers.mk

USER_CFLAGS= -I ../limber -I ${GSL_INC}
USER_LDFLAGS=-lcosmosis -L ${GSL_LIB} -lgsl -L../limber -llimber

all: interface.so

interface.so: utils.c shear_shear.c interface.c ../limber/liblimber.a
	$(CC) $(CFLAGS) $(LDFLAGS) -shared -o interface.so utils.c shear_shear.c interface.c

clean:
	rm -rf interface.so interface.so.dSYM
