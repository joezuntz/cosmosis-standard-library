include ${COSMOSIS_SRC_DIR}/config/compilers.mk

all: sigma.so

sigma.so : sigma.cpp
	${CXX} ${CXXFLAGS} sigma.cpp -shared -o sigma.so
	#g++ -fopenmp -lm -fPIC -O3 -std=c++11 sigma.cpp -shared -o sigma.so

clean :
	rm sigma.so


