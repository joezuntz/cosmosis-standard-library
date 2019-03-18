# modified for compatibility with COSMOSIS build

OBJECTS =  clusters.o fgas.o lensing.o wrapper.o xray.o ConfigFile.o
SEPOBJ = fwrapper.o

MYLIB = libclusters.so

default: objects $(MYLIB)

listobjects:
	@echo $(OBJECTS)

objects: $(OBJECTS) $(SEPOBJ)

$(MYLIB): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -shared -o $(CURDIR)/$(MYLIB) -L$(GSL_LIB) -lgsl $(OBJECTS)

clusters.o: clusters.hpp util.hpp
fgas.o: clusters.hpp fgas.hpp lensing.hpp util.hpp xray.hpp
lensing.o: clusters.hpp lensing.hpp  util.hpp ConfigFile.h
wrapper.o: clusters.hpp fgas.hpp lensing.hpp wrapper.hpp xray.hpp
xray.o: clusters.hpp util.hpp xray.hpp

%.o: %.cpp
	$(CXX) -o $@ -c $(CXXFLAGS) $*.cpp

clean:
	rm -f $(OBJECTS) $(SEPOBJ) $(MYLIB)
