
###### BUILDING PROJECT EXECUTABLES

ifndef PROJECT
PROJECT := Dalitz
endif

ifndef PROJECT_ROOT
PROJECT_ROOT := ..
endif

ifndef PROJECT_LIB
PROJECT_LIB := $(PROJECT_ROOT)/$(PROJECT)Lib
endif

# if there are other libraries to include, the include path can be added here
INC_DIR := -I$(PROJECT_LIB)
LIB_DIR := -L$(PROJECT_ROOT)/lib

# add extra library linking flags here if they are needed -- be sure the locations
# are specified in the LIB_DIR arugments above
#EXTRA_LIBS := -lOtherLibrary

######  most things below probably don't need adjustment

ifndef AMPTOOLS_ROOT
$(error Please set AMPTOOLS_ROOT to point to the root of the AmpTools source tree. )
endif

include $(AMPTOOLS_ROOT)/Makefile.settings

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
ATSUFFIX += _MPI
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

TARGET_EXE := $(NOMPI_SRC:%.cc=.%$(SUFFIX)_nompiexe)

# optionally if MPI is called for, we will build the MPI executables
ifdef MPI
TARGET_EXE += $(MPI_SRC:%.cc=.%$(SUFFIX)_mpiexe)
endif

.PHONY: default bin
.PRECIOUS: .%$(SUFFIX).o .%$(SUFFIX)_mpi.o

default: bin $(TARGET_EXE)

bin:
	$(Q)-mkdir $(INSTALL_BIN) >& /dev/null ||:

.%$(SUFFIX)_nompiexe: .%$(SUFFIX).o
	$(vecho) "--> Linking $*"
	$(Q)$(CXX) $(CXX_FLAGS) -o $(INSTALL_BIN)/$* $< \
	$(INC_DIR) $(LIB_DIR) $(LIBS)

.%$(SUFFIX)_mpiexe: .%$(SUFFIX)_mpi.o
	$(vecho) "--> Linking $*"
	$(Q)$(MPICXX) $(CXX_FLAGS) -o $(INSTALL_BIN)/$* $< \
	$(INC_DIR) $(LIB_DIR) $(LIBS)

.%$(SUFFIX).o: %.cc .%$(SUFFIX).d
	$(vecho) "-> Compiling $<"
	$(Q)$(CXX) $(CXX_FLAGS) -M -MP -MT $@ -o .$*$(SUFFIX).d $< $(INC_DIR)
	$(Q)$(CXX) $(CXX_FLAGS) -c -o $@ $< $(INC_DIR)

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
	$(Q)rm -f .*.o .*.d .*bin

-include $(SRCFILES:.cc=.dep)
