include ${COSMOSIS_SRC_DIR}/config/compilers.mk

PLANCK_DIR=plc-2.0

DEFINES=-D HAS_LAPACK -D LAPACK_CLIK -D NOHEALPIX -D CLIK_LENSING -D 'CLIKSVNVERSION="$(VERSION)"'
USER_CFLAGS=-I$(PLANCK_DIR)/include -I$(PLANCK_DIR)/src -D_GNU_SOURCE $(DEFINES)
USER_LDFLAGS=-L$(PLANCK_DIR)/lib -lcosmosis_fortran -lcosmosis $(RPATH_FLAGS) -L$(CFITSIO_LIB) -lcfitsio $(LAPACK_LINK)

RPATH_FLAGS=-Wl,-rpath,$(PWD)/$(PLANCK_DIR)/lib

ifeq ($(CFITSIO_INC)$(CFITSIO_LIB),)
$(error To compile the Planck likelihood you must activate the UPS package cfitsio or set the environment variables CFITSIO_INC and CFITSIO_LIB)
endif

ifeq ($(LAPACK_LINK),)
$(error To compile the Planck likelihood you must have the environment variable LAPACK_LINK set to tell me how to link to lapack on your system e.g. -L/usr/local/lib -llapack. If you are using the config/setup-cosmosis script then make sure everything is up to date and open a new terminal)
endif

all: planck planck_interface.so

planck: $(PLANCK_DIR)
	cd $(PLANCK_DIR) && $(MAKE) cosmosis && $(MAKE) install_cosmosis

planck_interface.so: planck planck_interface.c
	$(CC) -shared $(CFLAGS) -o $@ planck_interface.c -lclik $(LDFLAGS) $(RPATH_FLAGS)

clean:
	rm -rf *.o *.so *.dSYM
	rm -rf $(PLANCK_DIR)/lib
	rm -rf $(PLANCK_DIR)/bin
	rm -rf $(PLANCK_DIR)/share
	rm -rf $(PLANCK_DIR)/include
	cd $(PLANCK_DIR) && $(MAKE) clean

test:
	@echo Alas, planck has no tests

