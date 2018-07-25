include ${COSMOSIS_SRC_DIR}/config/compilers.mk

USER_LDFLAGS=-L$(GSL_LIB)  -lgsl -lgslcblas -lm -lcosmosis
USER_CFLAGS=-I$(GSL_INC)

emu.so: emu.c emu_interface.c params.h
	${CC}  ${CFLAGS} emu.c emu_interface.c -shared  -o  interface.so ${LDFLAGS}

clean:
	rm -f *.o