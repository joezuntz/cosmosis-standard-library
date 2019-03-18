include ${COSMOSIS_SRC_DIR}/config/compilers.mk
USER_CXXFLAGS+=-I${GSL_INC}
USER_LDFLAGS+=-L${GSL_LIB} -lgsl -lgslcblas


all: sigma.so

sigma.so : sigma.cpp
	${CXX} ${CXXFLAGS} sigma.cpp -shared -o $(CURDIR)/sigma.so ${LDFLAGS}

clean :
	rm -f sigma.so
	rm -rf sigma.so.dSYM


