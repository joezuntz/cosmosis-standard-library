include ${COSMOSIS_SRC_DIR}/config/compilers.mk

USER_FFLAGS=-I.
#USER_LDFLAGS=-lcosmosis_fortran

export EQUATIONS ?= equations_ppf
export RECOMBINATION ?= recfast

INTERFACE=camb.so
CAMBLIB=libcamb_$(RECOMBINATION).so
CAMBDIR=${PWD}/camb_Jan15

all: $(INTERFACE)

# purposefully chosen to be non-file target to
# always decend into CAMBDIR
$(CAMBDIR)/$(CAMBLIB)::
	cd $(CAMBDIR) && $(MAKE) $(CAMBLIB)

$(INTERFACE): $(CAMBDIR)/$(CAMBLIB) camb_interface.F90 camb_module.F90 $(CAMBDIR)/$(CAMBLIB)
	$(FC) $(FFLAGS) -shared -Wl,-rpath,@loader_path/../../../cosmosis/datablock  camb_interface.F90 camb_module.F90 -o $(INTERFACE) -I$(CAMBDIR) -lcamb_$(RECOMBINATION) -lcosmosis_fortran

clean:
	cd $(CAMBDIR) && $(MAKE) clean
	rm -rf *.so *.o *.mod camb.dSYM

test:
	@echo "Alas, camb has no tests"

install: $(CAMBDIR)/$(CAMBLIB)
	cp $(CAMBDIR)/$(CAMBLIB) /usr/local/lib
