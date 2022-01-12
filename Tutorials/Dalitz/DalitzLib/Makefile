
ifndef PROJECT
PROJECT := Dalitz
endif

ifndef PROJECT_ROOT
PROJECT_ROOT := ..
endif

ifndef AMPTOOLS_ROOT
$(error Please set AMPTOOLS_ROOT to point to the root of the AmpTools source tree. )
endif

include $(AMPTOOLS_ROOT)/Makefile.settings

############################################################
# Most things below here probably don't need to be changed
#
# (expert option for profile and tracing)
# to build code instrumented for Vampir Trace use VTRACE=1

ifndef INSTALL_LIB
INSTALL_LIB = $(PROJECT_ROOT)/lib
endif

ifdef VTRACE
CXX := vtcxx -vt:inst manual

NVCC := vtnvcc
CXX_FLAGS += -DVTRACE
endif

INC_DIR :=  -I$(PROJECT_ROOT)/$(PROJECT)Lib -I$(shell root-config --incdir) -I$(AMPTOOLS)
CXX_FLAGS += $(shell root-config --cflags)

SRCDIRS := $(PROJECT)DataIO $(PROJECT)Amp $(PROJECT)Plot

# check verbose output
ifdef VERBOSE
	Q :=
	vecho = @true
else
	Q := @
	vecho = @echo
endif

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
#$(addprefix clean_,$(SRCDIRS))

default: lib $(DEFAULT)

lib:
	$(Q)-mkdir $(INSTALL_LIB) >& /dev/null ||:
	
#test:
#	$(vecho) $(foreach exe,$(wildcard $(PROJECT)Exe/.*_bin),$(subst _bin,,$(subst $(PROJECT)Exe/.,bin/,$(exe))))

$(DEFAULT): $(TARGET_LIBS)
	$(foreach lib, $(TARGET_LIBS), $(shell cd $(INSTALL_LIB); ar -x $(lib) ) )
	$(Q)(cd $(INSTALL_LIB) && ar -rsc $@ .*$(SUFFIX).o && ranlib $@)
	$(Q)(cd $(INSTALL_LIB) && rm -f .*$(SUFFIX).o)
	$(Q)(cd $(INSTALL_LIB) && rm -f $(TARGET_LIBS))
	$(Q)-rm -f $(INSTALL_LIB)/__*
	$(vecho) "=== Build of $@ is complete. ==="

# descend into the directory obtained by stripping lib and
# $(SUFFIX).a and execute make with LIB set to the target
lib%.a:
	$(Q)$(MAKE) -C $(subst lib,, $(subst $(SUFFIX).a,, $@ )) LIB=$@
	$(Q)cp $(subst lib,, $(subst $(SUFFIX).a,, $@))/$@ $(INSTALL_LIB)/
	
clean:
	$(Q)-cd $(PROJECT)DataIO && $(MAKE) clean
	$(Q)-cd $(PROJECT)Amp && $(MAKE) clean
	$(Q)-cd $(PROJECT)Plot && $(MAKE) clean
