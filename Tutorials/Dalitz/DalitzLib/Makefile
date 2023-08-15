
ifndef PROJECT
PROJECT := Dalitz
endif

ifndef PROJECT_HOME
PROJECT_HOME := ..
endif

ifndef AMPTOOLS_HOME
$(error Please set AMPTOOLS_HOME to point to the root of the AmpTools source tree. )
endif

include $(AMPTOOLS_HOME)/Makefile.settings

############################################################
# Most things below here probably don't need to be changed
#
# (expert option for profile and tracing)
# to build code instrumented for Score-P use SCOREP=1

INC_DIR :=  -I$(PROJECT_HOME)/$(PROJECT)Lib -I$(shell root-config --incdir) -I$(AMPTOOLS)
CXX_FLAGS += $(shell root-config --cflags)

ifdef SCOREP
CXX := scorep-g++ --user
CXX_FLAGS += -DSCOREP_USER_ENABLE
INC_DIR += -I$(SCOREP_INC)
endif

ifndef INSTALL_LIB
INSTALL_LIB = $(PROJECT_HOME)/lib
endif

SRCDIRS := $(PROJECT)DataIO $(PROJECT)Amp $(PROJECT)Plot

# check if GPU build is requested
ifdef GPU
INC_DIR += -I$(CUDA_INSTALL_PATH)/include
CXX_FLAGS += -DGPU_ACCELERATION
SUFFIX := _GPU
endif

# setup the targets
DEFAULT := lib$(PROJECT)$(SUFFIX).a
TARGET_LIBS := $(SRCDIRS:%=lib%$(SUFFIX).a)

export

.PHONY: test default clean lib

default: lib $(DEFAULT)

lib:
	$(Q)-mkdir $(INSTALL_LIB) >& /dev/null ||:

$(DEFAULT): $(TARGET_LIBS)
	$(foreach lib, $(TARGET_LIBS), $(shell cd $(INSTALL_LIB); ar -x $(lib) ) )
	$(Q)(cd $(INSTALL_LIB) && ar -rsc $@ .*$(SUFFIX).o && ranlib $@)
	$(Q)(cd $(INSTALL_LIB) && rm -f .*$(SUFFIX).o)
	$(Q)(cd $(INSTALL_LIB) && rm -f $(TARGET_LIBS))
	$(Q)-rm -f $(INSTALL_LIB)/__*
	$(vecho) "=== Build of $@ is complete. ==="

# descend into the directory obtained by stripping lib and
# $(SUFFIX).a and execute make with LIB set to the target
lib%$(SUFFIX).a:
	$(Q)$(MAKE) -C $* LIB=$@
	$(Q)cp $*/$@ $(INSTALL_LIB)/

clean:
	$(Q)-cd $(PROJECT)DataIO && $(MAKE) clean
	$(Q)-cd $(PROJECT)Amp && $(MAKE) clean
	$(Q)-cd $(PROJECT)Plot && $(MAKE) clean
