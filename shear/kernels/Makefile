include ${COSMOSIS_SRC_DIR}/config/compilers.mk

CFLAGS+= -I ../limber
LDFLAGS+= -lcosmosis -lgsl -L../limber -llimber

all: interface.so

interface.so:utils.c shear_shear.c interface.c
	echo $(CFLAGS)
	$(CC) $(CFLAGS) $(LDFLAGS) -shared -o interface.so utils.c shear_shear.c interface.c

clean:
	rm -rf interface.so interface.so.dSYM
