SRCFILES	:= $(wildcard *.cc)
CUDAFILES	:= $(wildcard *.cu)

.PHONY: default clean
.PRECIOUS: %.o

default: $(LIB)

%.a:  $(SRCFILES:%.cc=.%$(SUFFIX).o) $(CUDAFILES:%.cu=.%$(SUFFIX).o)
	$(vecho) "--> Archiving $@"
	$(Q)ar -rsc $@ $^

.%$(SUFFIX).o: %.cc .%$(SUFFIX).d
	$(vecho) "-> Compiling $<"
	$(Q)$(CXX) $(CXX_FLAGS) -M -MP -MT $@ -o .$*$(SUFFIX).d $< $(INC_DIR)
	$(Q)$(CXX) $(CXX_FLAGS) -c -o $@ $< $(INC_DIR)

.%$(SUFFIX).o: %.cu .%$(SUFFIX).d
	$(vecho) "-> Compiling $<"
	$(Q)$(NVCC) $(NVCC_FLAGS) -M -MP -MT $@ -o .$*$(SUFFIX).d $< $(INC_DIR)
	$(Q)$(NVCC) $(NVCC_FLAGS) -c -o $@ $< $(INC_DIR)

DEPFILES := $(SRCFILES:%.cc=.%$(SUFFIX).d) $(CUDAFILES:%.cu=.%$(SUFFIX).d)
$(DEPFILES):

clean:
	$(Q)-rm -f .*.o .*.d *.a

-include $(DEPFILES)
