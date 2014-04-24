include ../config/compilers.mk
include ../config/subdirs.mk

SUBDIRS = boltzmann

.PHONY: all

all:
	$(process_subdirs)

test:
	$(process_subdirs)

clean:
	@echo Cleaning in $$(PWD)
	$(process_subdirs)
