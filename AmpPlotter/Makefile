
ifndef AMPTOOLS_HOME
$(error Please set AMPTOOLS_HOME to point to the root of the AmpTools source tree.)
endif

include $(AMPTOOLS_HOME)/Makefile.settings

###########################################################
# Most things below here probably don't need to be changed
#
# (expert option for profile and tracing)
# to build code instrumented for Score-P use SCOREP=1

ifdef SCOREP
CXX := scorep-g++ --user
CXX_FLAGS += -DSCOREP_USER_ENABLE
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

.PHONY: default clean force
.PRECIOUS: %.o

default: lib $(DEFAULT)

lib:
	$(Q)mkdir lib

lib%.a: compiler_flags
	$(Q)$(MAKE) -C $* LIB=$@
	$(Q)cp $*/$@ lib/


clean: $(addprefix clean_, $(SRCDIRS))
	$(Q)-rm -f lib/*.a

clean_%:
	$(Q)-cd $* && $(MAKE) clean

compiler_flags: force
	@echo '$(CXX_FLAGS)' | cmp -s - $@ || echo '$(CXX_FLAGS)' > $@
