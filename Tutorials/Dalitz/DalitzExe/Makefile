
INC_DIR += -I$(AMPPLOTTER)

# need some linking information available
LIB_DIR := -L$(shell root-config --libdir) $(CUDA_LIBS) \
	-L$(AMPTOOLS)/lib -L$(AMPPLOTTER)/lib -L../lib
	
LIBS :=  -l$(PROJECT)$(SUFFIX) -lAmpPlotter -lAmpTools$(ATSUFFIX) \
	$(shell root-config --glibs) $(CUDA_LINK) -lstdc++

# any source file falls into one of these two categories -- we'll
# assume that anything that has MPI in the name is an MPI executable
MPI_SRC := $(shell ls *.cc | grep -i mpi)
NOMPI_SRC := $(shell ls *.cc | grep -vi mpi)

TARGET_EXE := $(NOMPI_SRC:%.cc=.%$(SUFFIX)_nompiexe)

# optionally if MPI is called for, we will build the MPI executables
ifdef MPI
TARGET_EXE += $(MPI_SRC:%.cc=.%$(SUFFIX)_mpiexe)
endif

.PHONY: default
.PRECIOUS: .%$(SUFFIX).o .%$(SUFFIX)_mpi.o

default: $(TARGET_EXE)

test:
	$(vecho) "exesuffix: $(SUFFIX)"
	$(vecho) "suffix: $(SUFFIX)"
	$(vecho) $(TARGET_EXE)

.%$(SUFFIX)_nompiexe: .%$(SUFFIX).o
	$(vecho) "--> Linking $*"
	$(Q)$(CXX) $(CXX_FLAGS) -o $(subst _nompiexe,_bin,$@) $< \
	$(INC_DIR) $(LIB_DIR) $(LIBS)

.%$(SUFFIX)_mpiexe: .%$(SUFFIX)_mpi.o
	$(vecho) "--> Linking $*"
	$(Q)$(MPICXX) $(CXX_FLAGS) -o $(subst _mpiexe,_bin,$@) $< \
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
