all: sigma.so

sigma.so : sigma.cpp files/*
	g++ -fopenmp -lm -fPIC -O3 sigma.cpp -shared -o sigma.so

clean :
	rm sigma.so
