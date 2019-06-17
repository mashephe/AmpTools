
INC_DIR :=  -I.. -I$(shell root-config --incdir)

SRCDIRS := UpRootMinuit MinuitInterface IUAmpTools

CXX := gcc
CXX_FLAGS := -O3 $(shell root-config --cflags)
#CXX_FLAGS := $(shell root-config --cflags) -g

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
CUDA_FLAGS := -m64 -arch=sm_20
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
	mkdir lib

libAmpTools.a: $(TARGET_LIBS)
	$(foreach lib, $(TARGET_LIBS), $(shell cd lib; ar -x $(lib) ) )
	@cd lib && ar -rsv $@ *.o && ranlib $@
	@cd lib && rm -f *.o

libAmpTools_GPU.a: $(TARGET_LIBS) $(TARGET_LIBS_GPU)
	$(foreach lib_GPU, $(TARGET_LIBS_GPU), $(shell cd lib; ar -x $(lib_GPU) ) )
	@cd lib && ar -rsv $@ *.o && ranlib $@
	@cd lib && rm -f *.o

lib%_GPU.a: 
	@$(MAKE) -C $(subst lib,, $(subst _GPU.a,, $@ )) LIB=$@
	@cp $(subst lib,, $(subst _GPU.a,, $@))/$@ lib/

lib%.a: 
	@$(MAKE) -C $(subst lib,, $(subst .a,, $@ )) LIB=$@
	@cp $(subst lib,, $(subst .a,, $@))/$@ lib/


clean: $(addprefix clean_, $(SRCDIRS))
	-rm -f lib/*.a

clean_%:
	@-cd $(subst clean_,, $@) && $(MAKE) clean
