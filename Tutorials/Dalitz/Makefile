############################################################
# Using make to build the Dalitz Tutorial
#
# The user should be sure that root-config is in the path prior to building.
# The user should also define the following environment variables:
#
#     AMPTOOLS
#     AMPPLOTTER
#
# options:
#     MPI=1   (compile with MPI support)
#     GPU=1   (compile GPU-accelerated code, works with MPI also)
#     GPU_ARCH=arch_flag  (set the GPU architecture flag to arch_flag, e.g, sm_20)
#
#     VERBOSE=1   (shows make actions)
#
# Below are the common user-configurable options for building AmpTools.

PROJECT := Dalitz

CXX := gcc
CXX_FLAGS := -O3

# to use MPI, an MPI compiler should be in the user path -- this will be used
# to compile MPI executables.  In the project libraries (amplitudes, etc.)
# there should be no MPI specific code.

ifdef MPI
MPICXX := mpicxx
endif

# to use CUDA the environment variable CUDA_INSTALL_PATH must be defined
# and the proper architecture specified for the GPU.  In addition
# the nvcc compiler (likely in $CUDA_INSTALL_PATH/bin) needs to be in
# the user's path

NVCC :=	nvcc
NVCC_FLAGS := -m64

ifndef GPU_ARCH
GPU_ARCH := sm_20
endif

NVCC_FLAGS += -arch=$(GPU_ARCH)

ifdef GPU  # otherwise these variables will be empty to avoid problems
# the path to cuda specific libraries can vary with architecture
CUDA_LIBS := -L$(CUDA_INSTALL_PATH)/lib64
CUDA_LINK := -lcudart
endif

############################################################
# Most things below here probably don't need to be changed
#
# (expert option for profile and tracing)
# to build code instrumented for Vampir Trace use VTRACE=1

ifdef VTRACE
CXX := vtcxx -vt:inst manual

ifdef MPI
CXX += -vt:mpi
endif

NVCC := vtnvcc
CXX_FLAGS += -DVTRACE
endif

INC_DIR :=  -I.. -I$(shell root-config --incdir) -I$(AMPTOOLS)
CXX_FLAGS += $(shell root-config --cflags)

SRCDIRS := $(PROJECT)DataIO $(PROJECT)Amp $(PROJECT)Plot
ALLDIRS := $(SRCDIRS) $(PROJECT)Exe

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
ATSUFFIX := _GPU
endif

# check if MPI build is requested
#   * unlike AmpTools Makefile, don't modify SUFFIX because the only thing
#     this affects is the executables which will be handled separately
#   * however, we do need the right suffix to link MPI enabled
#     AmpTools library
ifdef MPI
CXX_FLAGS += -DUSE_MPI
ATSUFFIX += _MPI
endif

# setup the targets
DEFAULT := lib$(PROJECT)$(SUFFIX).a
TARGET_LIBS := $(SRCDIRS:%=lib%$(SUFFIX).a)

export

.PHONY: test default exe clean

default: lib bin $(DEFAULT) exe

lib:
	$(Q)mkdir lib
	
bin:
	$(Q)mkdir bin

test:
	$(vecho) $(foreach exe,$(wildcard $(PROJECT)Exe/.*_bin),$(subst _bin,,$(subst $(PROJECT)Exe/.,bin/,$(exe))))

$(DEFAULT): $(TARGET_LIBS)
	$(foreach lib, $(TARGET_LIBS), $(shell cd lib; ar -x $(lib) ) )
	$(Q)(cd lib && ar -rsc $@ .*$(SUFFIX).o && ranlib $@)
	$(Q)(cd lib && rm -f .*$(SUFFIX).o)
	$(Q)(cd lib && rm -f $(TARGET_LIBS))
	$(Q)-rm -f lib/__*
	$(vecho) "=== Build of $@ is complete. ==="

# descend into the directory obtained by stripping lib and
# $(SUFFIX).a and execute make with LIB set to the target
lib%.a:
	$(Q)$(MAKE) -C $(subst lib,, $(subst $(SUFFIX).a,, $@ )) LIB=$@
	$(Q)cp $(subst lib,, $(subst $(SUFFIX).a,, $@))/$@ lib/

exe:
	$(Q)$(MAKE) -C $(PROJECT)Exe
	$(foreach exe,$(wildcard $(PROJECT)Exe/.*_bin), \
	$(shell cp $(exe) $(subst _bin,,$(subst $(PROJECT)Exe/.,bin/,$(exe)))))
	$(vecho) "=== Build of $(PROJECT) executables is complete. ==="

# cleaning
clean: $(addprefix clean_, $(ALLDIRS))
	$(Q)-rm -f lib/*
	$(Q)-rm -f bin/*

clean_%:
	$(Q)-cd $(subst clean_,, $@) && $(MAKE) clean
