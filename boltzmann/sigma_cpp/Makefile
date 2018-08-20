all: sigma.so

sigma.so : sigma.cpp
	g++ -fopenmp -lm -fPIC -O3 -std=c++11 sigma.cpp -shared -o sigma.so

clean :
	rm sigma.so
