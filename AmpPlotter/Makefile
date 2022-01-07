
############################################################
# Using make to build AmpPlotter
#
# The user should be sure that root-config is in the path prior to building.
# The enviornment variable AMPTOOLS must be set to the AmpTools library
# top level.
#
# options:
#
#     VERBOSE=1   (shows make actions)
#
#  Note that the GUI AmpPlotter library has no GPU or MPI code.
#
# Below are the common user-configurable options for building AmpPlotter

CXX := gcc
CXX_FLAGS := -O3 $(shell root-config --cflags)


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
