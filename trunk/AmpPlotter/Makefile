
INC_DIR :=  -I.. -I$(shell root-config --incdir) -I$(AMPTOOLS)

SRCDIRS := AmpPlotter

CXX := gcc
CXX_FLAGS := -O3 $(shell root-config --cflags)

TARGET_LIBS := $(addsuffix .a, $(addprefix lib, $(SRCDIRS)))

DEFAULT := libAmpPlotter.a

export

.PHONY: default clean

default: lib $(DEFAULT)

lib:
	mkdir lib

lib%.a: 
	@$(MAKE) -C $(subst lib,, $(subst .a,, $@ )) LIB=$@
	@cp $(subst lib,, $(subst .a,, $@))/$@ lib/


clean: $(addprefix clean_, $(SRCDIRS))
	-rm -f lib/*.a

clean_%:
	@-cd $(subst clean_,, $@) && $(MAKE) clean
