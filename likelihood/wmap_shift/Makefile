include ${COSMOSIS_DIR}/config/compilers.mk
LDFLAGS+= -lcosmosis_fortran -lcosmosis 

all: wmap_shift.so

wmap_shift.so: wmap_shift.f90
	$(F90) $(FFLAGS) -shared -o wmap_shift.so wmap_shift.f90 $(LDFLAGS) 

clean:
	rm -f *.o *.mod wmap_shift.so