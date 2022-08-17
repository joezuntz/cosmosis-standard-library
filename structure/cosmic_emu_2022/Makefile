include ${COSMOSIS_SRC_DIR}/config/compilers.mk

USER_LDFLAGS=-L$(GSL_LIB)  -lgsl -lgslcblas -lm -lcosmosis
USER_CFLAGS=-I$(GSL_INC) -D EMU22_DATA_DIR=\"${PWD}\"

all: emu.so

emu.so: emu.c io.c emu_interface.c
	${CC}  ${CFLAGS} emu.c io.c emu_interface.c -shared  -o  emu.so ${LDFLAGS}

test:
	@echo "Alas, cosmic_emu has no tests"

clean:
	rm -f *.o emu.so
	rm -rf emu.so.dSYM
