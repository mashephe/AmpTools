
ifndef AMPTOOLS_ROOT
$(error Please set AMPTOOLS_ROOT to point to the root of the AmpTools source tree.)
endif

include $(AMPTOOLS_ROOT)/Makefile.settings

###########################################################
# Most things below here probably don't need to be changed
#
# (expert option for profile and tracing)
# to build code instrumented for Vampir Trace use VTRACE=1

ifdef VTRACE
CXX := vtcxx -vt:inst manual
CXX_FLAGS += -DVTRACE
endif

INC_DIR :=  -I.. -I$(shell root-config --incdir) -I$(AMPTOOLS)
SRCDIRS := AmpPlotter

# check verbose output
ifdef VERBOSE
	Q :=
	vecho = @true
else
	Q := @
	vecho = @echo
endif

CXX_FLAGS += $(shell root-config --cflags)

DEFAULT := libAmpPlotter.a

export

.PHONY: default clean
.PRECIOUS: %.o

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
