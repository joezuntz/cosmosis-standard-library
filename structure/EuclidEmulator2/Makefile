include ${COSMOSIS_SRC_DIR}/config/compilers.mk

PYTHON ?= python
PYTHONVERSION = $(shell $(PYTHON) -c 'import sys;print(f"{sys.version_info.major}.{sys.version_info.minor}")')
PYTHONINSTALL = $(PWD)/ee_install/lib/python$(PYTHONVERSION)/site-packages

DATA = ee2_bindata.dat
HEADERS = src/cosmo.h src/emulator.h src/units_and_constants.h
SRC = src/euclidemu2.pyx src/cosmo.cxx src/emulator.cxx

all: ee_install/done

ee_install/done: $(SRC) $(HEADERS) $(DATA)
	export CC=$(CC) CXX=$(CXX); export PYTHONPATH=${PYTHONPATH}:$(PYTHONINSTALL);  $(PYTHON) setup.py install   --prefix=$(PWD)/ee_install;  $(PYTHON)  setup.py install   --prefix=$(PWD)/ee_install
	touch ee_install/done

clean:
	rm -rf build dist ee_install

.PHONY: clean all
