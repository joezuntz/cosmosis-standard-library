# modified for compatibility with COSMOSIS build

USER_CXXFLAGS += -O3 -fPIC

OBJECTS =  clusters.o fgas.o lensing.o wrapper.o xray.o ConfigFile.o
SEPOBJ = fwrapper.o

ARCHIVES = libclusters.a

default: objects archives

listobjects:
	@echo $(OBJECTS)

objects: $(OBJECTS) $(SEPOBJ)

archives: $(ARCHIVES)

$(ARCHIVES): $(OBJECTS)

clusters.o: clusters.hpp util.hpp
fgas.o: clusters.hpp fgas.hpp lensing.hpp util.hpp xray.hpp
lensing.o: clusters.hpp lensing.hpp  util.hpp ConfigFile.h
wrapper.o: clusters.hpp fgas.hpp lensing.hpp wrapper.hpp xray.hpp
xray.o: clusters.hpp util.hpp xray.hpp



%.o: %.cpp
	$(CXX) -o $@ -c $(CXXFLAGS) $(INCLUDE) $*.cpp

%.a:
	rm -f $@; ar -r $@ $(OBJECTS)

clean:
	rm -f $(OBJECTS) $(SEPOBJ) $(ARCHIVES)
