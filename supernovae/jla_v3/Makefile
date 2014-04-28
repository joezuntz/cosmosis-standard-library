include ../../../Makefile.mine
# CLASS is required to compile the test code
#CLASSROOT=/home/marc/soft/class_v2.0.4
LDFLAGS+= -lblas -llapack -lgsl
CXXFLAGS+= -O3 -Wall -framework veclib -I../../../datablock -std=c++11

#CLASS_CXXFLAGS= -I$(CLASSROOT)/include/ -Isrc/
#CLASS_LFLAGS= -fopenmp $(CLASSROOT)/libclass.a

all: libjla.a jla.so

libjla.a: src/jla.o src/ini.o
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
