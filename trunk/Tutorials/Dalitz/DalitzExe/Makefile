
INC_DIR :=  -I.. -I$(shell root-config --incdir) \
  -I$(AMPTOOLS) -I$(AMPPLOTTER) -I$(DALITZ)

LIB_DIR := -L$(shell root-config --libdir) \
	-L$(AMPTOOLS)/lib -L$(AMPPLOTTER)/lib -L$(DALITZ)/lib

ROOTLIBS := $(shell root-config --glibs)

CXX := gcc
#CXX_FLAGS := -O3 -m32
CXX_FLAGS := -O3 $(shell root-config --cflags)


ifdef MPI

CXX := mpicxx -cc=$(CXX)
CXX_FLAGS += -DUSEMPI

endif


# To build code instrumented for Vampir Trace  VTRACE=1

ifdef VTRACE

CXX := vtcxx -vt:inst manual
MPICC := vtcxx -vt:mpi -vt:inst manual

CXX_FLAGS += -DVTRACE

endif


ifdef GPU

INC_DIR += -I$(CUDA_INSTALL_PATH)/include
LIB_DIR += -L$(CUDA_INSTALL_PATH)/lib

LIBS :=  -lDalitz_GPU -lAmpPlotter -lAmpTools_GPU $(ROOTLIBS) -lcudart -lstdc++
CXX_FLAGS += -DGPU_ACCELERATION

else

LIBS :=  -lDalitz -lAmpPlotter -lAmpTools $(ROOTLIBS) -lstdc++

endif


SRCFILES := $(wildcard *.cc)
TARGET_EXE := $(subst .cc,_exe,$(SRCFILES))

.PHONY: default

default: $(TARGET_EXE)

# generic build rule (requires all libraries)

%_exe: %.o
	$(CXX) $(CXX_FLAGS) -o $(subst _exe,,$@) $< -O3 $(INC_DIR) $(LIB_DIR) $(LIBS)

%.o : %.cc
	$(CXX) $(CXX_FLAGS) -M -o $*.d $< $(INC_DIR); \
	cp $*.d $*.dep; \
	sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
		-e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $*.dep; \
	rm -f $*.d $*.d.bak
	$(CXX) $(CXX_FLAGS) -c -o $@ $< $(INC_DIR)

clean:
	rm -f *.o *.dep $(subst _exe,,$(TARGET_EXE))

-include $(SRCFILES:.cc=.dep)
