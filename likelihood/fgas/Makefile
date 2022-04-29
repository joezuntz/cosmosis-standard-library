
# If you already have your own Makefile you can 
# replace all of this, but you need to keep this line
# at the top:
include ${COSMOSIS_SRC_DIR}/config/compilers.mk

USER_CXXFLAGS += -I${GSL_INC}
USER_LDFLAGS += -lcosmosis -L${GSL_LIB}  -lgsl -lgslcblas -Wl,-rpath,$(PWD)/src

#For python, leave this as is.  For C/C++/Fortran, 
#comment the first line and uncomment the second:
#all: 
all: fgas.so

# ... and then provide the command to compile it, perhaps
# like this for a simple example in C:
src/libclusters.so:
	$(MAKE) -C src

fgas.so: fgas_cosmosis.cpp src/libclusters.so
	$(CXX) $(CXXFLAGS) fgas_cosmosis.cpp -shared -o fgas.so  -Lsrc -lclusters $(LDFLAGS)

#Replace this to put in some kind of test command for
#your module, if you have one.
test:
	@echo "Alas, cluster_fgas/fgas module has no tests"

#Add anything else to clean here
clean:
	rm -fr fgas.so fgas.so.dSYM
	$(MAKE) -C src clean
