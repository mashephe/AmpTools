
INC_DIR :=  -I.. -I$(shell root-config --incdir) -I$(AMPTOOLS)

SRCDIRS := AmpPlotter

CXX := gcc
CXX_FLAGS := -O3 $(shell root-config --cflags)

# for more detailed info, use VERBOSE=1
ifdef VERBOSE
	Q :=
	vecho = @true
else
	Q := @
	vecho = @echo
endif

TARGET_LIBS := $(addsuffix .a, $(addprefix lib, $(SRCDIRS)))

DEFAULT := libAmpPlotter.a

export

.PHONY: default clean

default: lib $(DEFAULT)

lib:
	$(Q)mkdir lib

lib%.a: 
	$(Q)$(MAKE) -C $(subst lib,, $(subst .a,, $@ )) LIB=$@
	$(Q)cp $(subst lib,, $(subst .a,, $@))/$@ lib/


clean: $(addprefix clean_, $(SRCDIRS))
	$(Q)-rm -f lib/*.a

clean_%:
	$(Q)-cd $(subst clean_,, $@) && $(MAKE) clean
