include ${COSMOSIS_SRC_DIR}/config/compilers.mk

USER_CXXFLAGS+= -I ${GSL_INC}

ifneq (Darwin, $(OS))
  USER_CXXFLAGS+=-Wl,-rpath,${PWD}
endif

USER_LDFLAGS+=  -lcosmosis ${LAPACK_LINK} -L ${GSL_LIB} -lgslcblas -lgsl

ifneq (Darwin, $(OS))
  LDFLAGS+=-ldl
endif

all: libjla.so jla.so

libjla.so: src/jla.o
	$(CXX) -shared $(LDFLAGS) -o $(CURDIR)/$@ $^

jla.so: interface.cc libjla.so
	$(CXX) $(CXXFLAGS) -shared -o $@ interface.cc -L. -ljla $(LDFLAGS)

test_jla: test_jla.c jla.so
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

clean:
	rm -f src/*.o
	rm -f libjla.so jla.so
	rm -rf jla.so.dSYM
	rm -f test_jla
	rm -rf test_jla.dSYM

test: test_jla
	@echo Running test_jla
	@./test_jla && echo "passed" || echo "  failed !"
