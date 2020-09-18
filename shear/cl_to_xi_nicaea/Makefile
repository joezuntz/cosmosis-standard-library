include ${COSMOSIS_SRC_DIR}/config/compilers.mk

USER_LDFLAGS=-L${GSL_LIB} -lgsl -lgslcblas -lcosmosis -L${FFTW_LIBRARY} -lfftw3
USER_CFLAGS=-I${GSL_INC} -Wno-deprecated-declarations -I${FFTW_INCLUDE_DIR} -std=c99

OBJS = io.c errorlist.c maths.c cl_to_xi_interface.c

all: nicaea_interface.so

nicaea_interface.so: $(OBJS)
	$(CC) $(CFLAGS) -shared -o $@ $^ $(LDFLAGS)

test:
	@echo "Alas, cl_to_xi_nicaea has no tests"

clean:
	rm -f nicaea_interface.so
	rm -rf nicaea_interface.so.dSYM
