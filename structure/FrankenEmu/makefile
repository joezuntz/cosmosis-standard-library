include ${COSMOSIS_SRC_DIR}/config/compilers.mk

USER_CFLAGS= -I ${GSL_INC}
USER_LDFLAGS= -lcosmosis -L ${GSL_LIB} -lgsl -lgslcblas 


LIBRARY_FILES = emu.o hubble.o emu_noh.o
P = main.o 

all: interface.so

emu.exe: $(P) $(LIBRARY_FILES)
	$(CC) -o $@ $(CFLAGS) $(P) $(LIBRARY_FILES)  $(LDFLAGS)

interface.so: interface.c $(LIBRARY_FILES)
	$(CC) $(CFLAGS) -shared -o $@  $^ $(LDFLAGS)


clean:
	/bin/rm -rf *.o *.so *.dSYM emu.exe

