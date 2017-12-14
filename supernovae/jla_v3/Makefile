include ${COSMOSIS_SRC_DIR}/cosmosis/compilers.mk


USER_CXXFLAGS+= -I ${GSL_INC}
USER_LDFLAGS+=  -lcosmosis ${LAPACK_LINK}   -lgfortran -L ${GSL_LIB} -lgslcblas -lgsl 


all: libjla.a jla.so

libjla.a: src/jla.o
	$(AR) $(ARFLAGS) $@ $^

jla.so: interface.cc libjla.a
	$(CXX) $(CXXFLAGS) -shared -o $@ interface.cc -L. -ljla $(LDFLAGS)
#test_jla: src/test.o libjla.a
#	$(CXX) -o $@ $^ $(CXXFLAGS) $(LFLAGS) $(CLASS_LFLAGS)

#src/test.o: src/test.cc
#	$(CXX) -c -o $@ $^  $(CXXFLAGS) $(CLASS_CXXFLAGS)

clean:
	rm -f src/*.o
	rm -f libjla.a

test:
	@echo Alas, supernovae/jla_v3 has no tests

