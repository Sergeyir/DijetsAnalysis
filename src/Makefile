# Set the shell.
SHELL=/usr/bin/env bash

# Include the configuration.
-include Makefile.inc

# Check distribution (use local version first, then installed version).
ifneq ("$(wildcard ../lib/libpythia8.*)","")
  PREFIX_LIB=../lib
  PREFIX_INCLUDE=../include
endif
CXX_COMMON:=-I$(PREFIX_INCLUDE) $(CXX_COMMON) $(GZIP_LIB)
CXX_COMMON+= -L$(PREFIX_LIB) -Wl,-rpath,$(PREFIX_LIB) -lpythia8 -ldl
#PYTHIA=$(PREFIX_LIB)/libpythia8$(LIB_SUFFIX)

# Rules without physical targets (secondary expansion for specific rules).
.SECONDEXPANSION:
.PHONY: all clean

generate: generate.cpp
ifeq ($(ROOT_USE)${FASTJET3_USE},11)
	$(CXX) $@.cpp -o $@.exe -w $(CXX_COMMON) \
	$(FASTJET3_INCLUDE) $(FASTJET3_LIB) \
	$(ROOT_LIB) `$(ROOT_CONFIG) --cflags --glibs`
else
	$(error Error: $@ requires ROOT and FASTJET)
endif

# Clean.
clean:
	rm generate.exe \
	rm -f *~; rm -f \
