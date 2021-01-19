
INC_DIR :=  -I.. -I$(shell root-config --incdir) -I$(AMPTOOLS)

CXX := gcc
CXX_FLAGS := -O3 $(shell root-config --cflags)

SRCDIRS := DalitzDataIO DalitzAmp DalitzPlot
TARGET_LIBS := $(addsuffix .a, $(addprefix lib, $(SRCDIRS)))

SRCDIRS_GPU := $(SRCDIRS)
TARGET_LIBS_GPU :=  $(addsuffix _GPU.a, $(addprefix lib, $(SRCDIRS_GPU)))

# for more detailed info, use VERBOSE=1
ifdef VERBOSE
	Q :=
	vecho = @true
else
	Q := @
	vecho = @echo
endif

#To build GPU-accelerated code type: make GPU=1
ifdef GPU

NVCC := nvcc
CUDA_FLAGS := -arch=sm_20
INC_DIR += -I. -I$(CUDA_INSTALL_PATH)/include

CXX_FLAGS += -DGPU_ACCELERATION
DEFAULT := libDalitz_GPU.a

# To build code instrumented for Vampir Trace  VTRACE=1

ifdef VTRACE

CXX := vtcxx -vt:inst manual
NVCC := vtnvcc

CXX_FLAGS += -DVTRACE

endif


else

DEFAULT := libDalitz.a

endif

export

.PHONY: default clean
.PRECIOUS: %.o

default: lib $(DEFAULT)

lib:
	$(Q)mkdir lib

libDalitz.a: $(TARGET_LIBS)
	$(foreach lib, $(TARGET_LIBS), $(shell cd lib; ar -x $(lib) ) )
	$(Q)(cd lib && ar -rsc $@ *.o && ranlib $@)
	$(Q)(cd lib && rm -f *.o)
	$(vecho) "=== Build of $@ is complete. ==="

libDalitz_GPU.a: $(TARGET_LIBS_GPU)
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

clean: $(addprefix clean_, $(SRCDIRS))
	$(Q)-rm -f lib/*.a

clean_%:
	$(Q)-cd $(subst clean_,, $@) && $(MAKE) clean
