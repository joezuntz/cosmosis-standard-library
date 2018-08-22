include ${COSMOSIS_SRC_DIR}/config/compilers.mk

all: sigma.so

sigma.so : sigma.cpp
	${CXX} ${CXXFLAGS} sigma.cpp -shared -o sigma.so -lgsl -lgslcblas

clean :
	rm -f sigma.so


