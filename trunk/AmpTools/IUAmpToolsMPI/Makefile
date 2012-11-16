SRCFILES	:= $(wildcard *.cc)

.PHONY: default

default: $(LIB)

%.a: $(SRCFILES:.cc=.o) $(CUDAFILES:.cu=.o)
	ar -rsv $@ $^

%.o : %.cc
	$(MPICC) $(CXX_FLAGS) -M -o $*.d $< $(INC_DIR) ; \
	cp $*.d $*.dep; \
	sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
		-e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $*.dep; \
	rm -f $*.d $*.d.bak
	$(MPICC) $(CXX_FLAGS) -c -o $@ $< $(INC_DIR) 

clean:
	rm *.o *.dep *.a

-include $(SRCFILES:.cc=.dep)
