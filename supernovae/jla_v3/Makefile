include ${COSMOSIS_SRC_DIR}/config/compilers.mk


USER_CXXFLAGS+= -I ${GSL_INC}
USER_LDFLAGS+=  -lcosmosis ${LAPACK_LINK} -L ${GSL_LIB} -lgslcblas -lgsl


all: libjla.so jla.so

libjla.so: src/jla.o
	$(CXX) -shared $(LDFLAGS) -o $(CURDIR)/$@ $^

jla.so: interface.cc libjla.so
	$(CXX) $(CXXFLAGS) -shared -o $@ interface.cc -L. -ljla $(LDFLAGS)
#test_jla: src/test.o libjla.a
#	$(CXX) -o $@ $^ $(CXXFLAGS) $(LFLAGS) $(CLASS_LFLAGS)

#src/test.o: src/test.cc
#	$(CXX) -c -o $@ $^  $(CXXFLAGS) $(CLASS_CXXFLAGS)

clean:
	rm -f src/*.o
	rm -f libjla.so jla.so
	rm -rf jla.so.dSYM

test:
	@echo Alas, supernovae/jla_v3 has no tests

