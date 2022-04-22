
ifndef AMPTOOLS_HOME
$(error Please set AMPTOOLS_HOME to point to the root of the AmpTools source tree. )
endif

include $(AMPTOOLS_HOME)/Makefile.settings

############################################################
# Most things below here probably don't need to be changed
#
# (expert option for profile and tracing)
# to build code instrumented for Score-P use SCOREP=1

ifdef SCOREP
CXX := scorep-g++ --user

CXX_FLAGS += -DSCOREP_USER_ENABLE
endif

INC_DIR :=  -I.. -I$(shell root-config --incdir) -I$(SCOREP_INC)
CXX_FLAGS += $(shell root-config --cflags)

SRCDIRS := UpRootMinuit MinuitInterface IUAmpTools
ALLDIRS := $(SRCDIRS) IUAmpToolsMPI GPUManager

# check if GPU build is requested
ifdef GPU
SRCDIRS += GPUManager
INC_DIR += -I$(CUDA_INSTALL_PATH)/include
CXX_FLAGS += -DGPU_ACCELERATION
SUFFIX := _GPU
endif

# check if MPI build is requested
ifdef MPI
CXX := $(MPICXX)
ifdef SCOREP
CXX := scorep --mpp=mpi --user mpicxx
endif
SRCDIRS += IUAmpToolsMPI
CXX_FLAGS += -DUSE_MPI
SUFFIX := $(SUFFIX)_MPI
endif

# setup the targets
DEFAULT := libAmpTools$(SUFFIX).a
TARGET_LIBS := $(SRCDIRS:%=lib%$(SUFFIX).a)

export

.PHONY: default clean

default: lib $(DEFAULT)

lib:
	$(Q)mkdir lib

$(DEFAULT): $(TARGET_LIBS)
	$(foreach lib, $(TARGET_LIBS), $(shell cd lib; ar -x $(lib) ) )
	$(Q)(cd lib && ar -rsc $@ .*$(SUFFIX).o && ranlib $@)
	$(Q)(cd lib && rm -f .*$(SUFFIX).o)
	$(Q)(cd lib && rm -f $(TARGET_LIBS))
	$(Q)-rm -f lib/__*
	$(vecho) "=== Build of $@ is complete. ==="

# descend into the directory obtained by stripping lib and
# $(SUFFIX).a and execute make with LIB set to the target
lib%$(SUFFIX).a:
	$(Q)$(MAKE) -C $* LIB=$@
	$(Q)cp $*/$@ lib/

# cleaning
clean: $(addprefix clean_, $(ALLDIRS))
	$(Q)-rm -f lib/*

clean_%:
	$(Q)-cd $* && $(MAKE) clean
