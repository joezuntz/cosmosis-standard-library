include ../../../Makefile.mine

CFLAGS+=-I../../../datablock -std=c99 -I ../limber
LDFLAGS+=-L../../../datablock -lcosmosis -lgsl -L../limber -llimber

interface.so:utils.c shear_shear.c interface.c
	$(CC) $(CFLAGS) $(LDFLAGS) -shared -o interface.so utils.c shear_shear.c interface.c

clean:
	rm -rf interface.so interface.so.dSYM