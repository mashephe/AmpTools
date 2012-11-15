
INC_DIR :=  -I.. -I$(CLHEP_INCLUDE_DIR) -I$(ROOTSYS)/include -I$(AMPTOOLS)

SRCDIRS := AmpPlotter

CXX := gcc
CXX_FLAGS := -O3

TARGET_LIBS := $(addsuffix .a, $(addprefix lib, $(SRCDIRS)))

#To build GPU-accelerated code type: make GPU=1

ifdef GPU

TARGET_LIBS_GPU :=  $(addsuffix _GPU.a, $(addprefix lib, $(SRCDIRS)))

NVCC :=	nvcc
CUDA_FLAGS := -arch=compute_11  -code=compute_11
INC_DIR += -I$(CUDA_INSTALL_PATH)/include 

CXX_FLAGS += -DGPU_ACCELERATION
DEFAULT := libAmpPlotter_GPU.a

else

DEFAULT := libAmpPlotter.a

endif

export

.PHONY: default clean

default: lib $(DEFAULT)

lib:
	mkdir lib

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
