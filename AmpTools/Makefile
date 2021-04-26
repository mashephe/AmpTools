
INC_DIR :=  -I.. -I$(shell root-config --incdir)

SRCDIRS := UpRootMinuit MinuitInterface IUAmpTools
ALLDIRS := $(SRCDIRS) IUAmpToolsMPI GPUManager

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

#To build in MPI-exentsions type: make MPI=1
#(this can be used with GPU, e.g.:  make MPI=1 GPU=1)

ifdef MPI

SRCDIRS += IUAmpToolsMPI

MPICC := mpicxx 
CXX_FLAGS += -DUSE_MPI

endif

TARGET_LIBS := $(addsuffix .a, $(addprefix lib, $(SRCDIRS)))

#To build GPU-accelerated code type: make GPU=1

ifdef GPU

SRCDIRS += GPUManager
TARGET_LIBS_GPU :=  $(addsuffix _GPU.a, $(addprefix lib, $(SRCDIRS)))

NVCC :=	nvcc
CUDA_FLAGS := -m64 -arch=sm_75
INC_DIR += -I$(CUDA_INSTALL_PATH)/include 

CXX_FLAGS += -DGPU_ACCELERATION
DEFAULT := libAmpTools_GPU.a

else

DEFAULT := libAmpTools.a

endif

# To build code instrumented for Vampir Trace  VTRACE=1

ifdef VTRACE

CXX := vtcxx -vt:inst manual
MPICC := vtcxx -vt:mpi -vt:inst manual
NVCC := vtnvcc

CXX_FLAGS += -DVTRACE

endif

export

.PHONY: default clean

default: lib $(DEFAULT)

lib:
	$(Q)mkdir lib

libAmpTools.a: $(TARGET_LIBS)
	$(foreach lib, $(TARGET_LIBS), $(shell cd lib; ar -x $(lib) ) )
	$(Q)(cd lib && ar -rsc $@ *.o && ranlib $@)
	$(Q)(cd lib && rm -f *.o)
	$(vecho) "=== Build of $@ is complete. ==="

libAmpTools_GPU.a: $(TARGET_LIBS) $(TARGET_LIBS_GPU)
	$(foreach lib_GPU, $(TARGET_LIBS_GPU), $(shell cd lib; ar -x $(lib_GPU) ) )
	$(Q)(cd lib && ar -rsc $@ *.o && ranlib $@)
	$(Q)(cd lib && rm -f *.o)
	$(vecho) "=== Build of $@ is complete. ==="

lib%_GPU.a: 
	$(Q)$(MAKE) -C $(subst lib,, $(subst _GPU.a,, $@ )) LIB=$@
	$(Q)cp $(subst lib,, $(subst _GPU.a,, $@))/$@ lib/

lib%.a: 
	$(Q)$(MAKE) -C $(subst lib,, $(subst .a,, $@ )) LIB=$@
	$(Q)cp $(subst lib,, $(subst .a,, $@))/$@ lib/


clean: $(addprefix clean_, $(ALLDIRS))
	$(Q)-rm -f lib/*.a

clean_%:
	$(Q)-cd $(subst clean_,, $@) && $(MAKE) clean
