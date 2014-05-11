include ${COSMOSIS_SRC_DIR}/config/compilers.mk
include ${COSMOSIS_SRC_DIR}/config/subdirs.mk

SUBDIRS = planck wmap_shift

ifeq ($(WMAPDATA_VERSION),v5_00)
SUBDIRS += wmap9
endif
ifeq ($(WMAPDATA_VERSION),v4_00)
SUBDIRS += wmap7
endif
