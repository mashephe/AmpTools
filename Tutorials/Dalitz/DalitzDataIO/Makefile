
SRCFILES := $(wildcard *.cc)

.PHONY: default

default: $(LIB)

%.a: $(SRCFILES:.cc=.o)
	ar -rsv $@ $^

%.o : %.cc
	$(CXX) $(CXX_FLAGS) -M -o $*.d $< $(INC_DIR); \
	cp $*.d $*.dep; \
	sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
		-e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $*.dep; \
	rm -f $*.d $*.d.bak
	$(CXX) $(CXX_FLAGS) -c -o $@ $< $(INC_DIR)

clean:
	rm *.o *.dep *.a

-include $(SRCFILES:.cc=.dep)
