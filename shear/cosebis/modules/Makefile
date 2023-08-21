include ${COSMOSIS_SRC_DIR}/config/compilers.mk
#-----------------------------------
#--- objects to build
#-----------------------------------
USER_LDFLAGS += -L${GSL_LIB}  -lgsl -lgslcblas
USER_CXXFLAGS += -I${GSL_INC} -D COSEBIS_DIR="\"${PWD}/../\""

ifneq (Darwin, $(OS))
  USER_CXXFLAGS+=-Wl,-rpath,${PWD}
endif


OBJS     = matrix.o function_cosebis.o ticker.o tostring.o tolower.o timestamp.o \
		   globaldef.o Integrate.o COSEBIs.o WnLog.o TpnLogRadian.o errors.o 

LIBNAME  = libcolib.a

HEADER   = $(OBJS:.o=.h) 

#-----------------------------------
#--- libraries
#-----------------------------------

# make objects (components) of library
%.o : %.cc %.h
	$(CXX)  $(CXXFLAGS) ${USER_CXXFLAGS} -c $<

# make library
all: besselzero04.h $(HEADER) $(OBJS)
	rm -f $(LIBNAME)
	ar rv $(LIBNAME) $(OBJS)
	ranlib $(LIBNAME)

besselzero04.h:
	python generate-bessel-zeros.py

.PHONY: clean
clean:
	rm -f $(LIBNAME)
	rm -f *.o besselzero04.h
