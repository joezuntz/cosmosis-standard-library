include ${COSMOSIS_SRC_DIR}/config/compilers.mk

COLIB     = modules/libcolib.a

ALL_LIBS = cl_to_cosebis/cl_to_cosebis.so

all: $(ALL_LIBS)
     
$(COLIB): modules/*.cc modules/*.h
	cd modules && $(MAKE)

cl_to_cosebis/cl_to_cosebis.so: $(COLIB) cl_to_cosebis/cl_to_cosebis_interface.cc
	cd cl_to_cosebis && $(MAKE)

clean:
	rm -rf */*.so */*.so.dSYM
	cd modules && $(MAKE) clean

	
	
