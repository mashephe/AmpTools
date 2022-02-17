
###### BUILDING PROJECT EXECUTABLES

ifndef PROJECT
PROJECT := Dalitz
endif

ifndef PROJECT_HOME
PROJECT_HOME := ..
endif

ifndef PROJECT_LIB
PROJECT_LIB := $(PROJECT_HOME)/$(PROJECT)Lib
endif

# if there are other libraries to include, the include path can be added here
INC_DIR := -I$(PROJECT_LIB)
LIB_DIR := -L$(PROJECT_HOME)/lib

# add extra library linking flags here if they are needed -- be sure the locations
# are specified in the LIB_DIR arugments above
#EXTRA_LIBS := -lOtherLibrary

######  most things below probably don't need adjustment

ifndef AMPTOOLS_HOME
$(error Please set AMPTOOLS_HOME to point to the root of the AmpTools source tree. )
endif

include $(AMPTOOLS_HOME)/Makefile.settings

# leave the binaries in this directory if the install directory hasn't been set
ifndef INSTALL_BIN
INSTALL_BIN = .
endif

INC_DIR += -I$(AMPPLOTTER) -I$(AMPTOOLS) -I$(shell root-config --incdir)
CXX_FLAGS += $(shell root-config --cflags)

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
ATSUFFIX := $(ATSUFFIX)_MPI
endif

# need some linking information available
LIB_DIR += -L$(shell root-config --libdir) $(CUDA_LIBS) \
	-L$(AMPTOOLS)/lib -L$(AMPPLOTTER)/lib

LIBS :=  -l$(PROJECT)$(SUFFIX) -lAmpPlotter -lAmpTools$(ATSUFFIX) \
	$(EXTRA_LIBS) $(shell root-config --glibs) $(CUDA_LINK) -lstdc++

# any source file falls into one of these two categories -- we'll
# assume that anything that has MPI in the name is an MPI executable
MPI_SRC := $(shell ls *.cc | grep -i mpi)
NOMPI_SRC := $(shell ls *.cc | grep -vi mpi)

TARGET_EXE := $(NOMPI_SRC:%.cc=%$(SUFFIX))

# optionally if MPI is called for, we will build the MPI executables
ifdef MPI
TARGET_EXE += $(MPI_SRC:%.cc=%$(SUFFIX))
endif

.PHONY: default bin
.PRECIOUS: .%$(SUFFIX).o .%$(SUFFIX)_mpi.o

default: bin $(TARGET_EXE)

bin:
	$(Q)-mkdir $(INSTALL_BIN) >& /dev/null ||:

# for the non-MPI executables we want to strip out the -DUSEMPI flag
# from the compiler arguments and also be sure we link against the
# non-MPI version of AmpTools

$(NOMPI_SRC:%.cc=%$(SUFFIX)): %$(SUFFIX) : .%$(SUFFIX).o
	$(vecho) "--> Linking $*$(SUFFIX)"
	$(Q)$(CXX) $(subst -DUSEMPI,,$(CXX_FLAGS)) -o $(INSTALL_BIN)/$*$(SUFFIX) $< \
	$(INC_DIR) $(LIB_DIR) $(subst _MPI,,$(LIBS))

$(MPI_SRC:%.cc=%$(SUFFIX)): %$(SUFFIX) : .%$(SUFFIX)_mpi.o
	$(vecho) "--> Linking $*$(SUFFIX)"
	$(Q)$(MPICXX) $(CXX_FLAGS) -o $(INSTALL_BIN)/$*$(SUFFIX) $< \
	$(INC_DIR) $(LIB_DIR) $(LIBS)

.%$(SUFFIX).o: %.cc .%$(SUFFIX).d
	$(vecho) "-> Compiling $<"
	$(Q)$(CXX) $(subst -DUSEMPI,,$(CXX_FLAGS)) -M -MP -MT $@ -o .$*$(SUFFIX).d $< $(INC_DIR)
	$(Q)$(CXX) $(subst -DUSEMPI,,$(CXX_FLAGS)) -c -o $@ $< $(INC_DIR)

.%$(SUFFIX)_mpi.o: %.cc .%$(SUFFIX).d
	$(vecho) "-> Compiling $<"
	$(Q)$(MPICXX) $(CXX_FLAGS) -M -MP -MT $@ -o .$*$(SUFFIX).d $< $(INC_DIR)
	$(Q)$(MPICXX) $(CXX_FLAGS) -c -o $@ $< $(INC_DIR)

DEPFILES := $(NOMPI_SRC:%.cc=.%$(SUFFIX).d)
ifdef MPI
DEPFILES += $(MPI_SRC:%.cc=.%$(SUFFIX).d)
endif
$(DEPFILES):

clean:
	$(Q)rm -f .*.o .*.d $(TARGET_EXE)

-include $(DEPFILES)
