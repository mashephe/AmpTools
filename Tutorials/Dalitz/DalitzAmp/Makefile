
SRCFILES := $(wildcard *.cc)
FORTRAMPS := $(wildcard *.f)

ifdef GPU
CUDAFILES := $(wildcard *cu)
endif

.PHONY: default
.PRECIOUS: %.o

default: $(LIB)

%.a: $(SRCFILES:.cc=.o) $(FORTRAMPS:.f=.o) $(CUDAFILES:.cu=.o)
	$(vecho) "--> Archiving $@"
	$(Q)ar -rsc $@ $^

%.o : %.f
	$(vecho) "-> Compiling $<"
	$(Q)gfortran -c -O3 -o $@ $<

%.o : %.cc
	$(vecho) "-> Compiling $<"
	$(Q)$(CXX) $(CXX_FLAGS) -M -o $*.d $< $(INC_DIR) ; \
	cp $*.d $*.dep; \
	sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
		-e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $*.dep; \
	rm -f $*.d $*.d.bak
	$(Q)$(CXX) $(CXX_FLAGS) -c -o $@ $< $(INC_DIR)
	
%.o : %.cu
	$(vecho) "-> Compiling $<"
	$(Q)$(NVCC) $(CUDA_FLAGS) -M -o $*.d $< $(INC_DIR) ; \
	cp $*.d $*.dep; \
	sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
		-e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $*.dep; \
	rm -f $*.d $*.d.bak
	$(Q)$(NVCC) $(CUDA_FLAGS) -c -o $@ $< $(INC_DIR)

clean:
	$(Q)-rm -f *.o *.dep *.a

-include $(SRCFILES:.cc=.dep)
