include ${COSMOSIS_SRC_DIR}/config/compilers.mk

USER_FFLAGS=-I.
USER_LDFLAGS=-lcosmosis_fortran

export EQUATIONS ?= equations_ppf
export RECOMBINATION ?= recfast

INTERFACE=camb.so
CAMBLIB=libcamb_$(RECOMBINATION).a
CAMBDIR=camb_Nov13

all: $(INTERFACE)

# purposefully chosen to be non-file target to
# always decend into CAMBDIR
$(CAMBDIR)/$(CAMBLIB)::
	cd $(CAMBDIR) && $(MAKE) $(CAMBLIB)

$(INTERFACE): $(CAMBDIR)/$(CAMBLIB) camb_interface.F90 camb_module.F90
	$(FC) $(FFLAGS) -shared camb_interface.F90 camb_module.F90 -o $(INTERFACE) $(LDFLAGS) -I$(CAMBDIR)  -L$(CAMBDIR) -lcamb_$(RECOMBINATION)

clean:
	cd $(CAMBDIR) && $(MAKE) clean
	rm -rf *.so *.o *.mod *.dSYM

test:
	@echo "Alas, camb has no tests"
