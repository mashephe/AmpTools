SRCFILES	:= $(wildcard *.cc)

.PHONY: default
.PRECIOUS: %.o

default: $(LIB)

%.a: $(SRCFILES:.cc=.o)
	$(vecho) "--> Archiving $@"
	$(Q)ar -rsc $@ $^

%.o : %.cc
	$(vecho) "-> Compiling $<"
	$(Q)$(CXX) $(CXX_FLAGS) -M -o $*.d $< $(INC_DIR) ; \
	cp $*.d $*.dep; \
	sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
		-e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $*.dep; \
	rm -f $*.d $*.d.bak
	$(Q)$(CXX) $(CXX_FLAGS) -c -o $@ $< $(INC_DIR)

clean:
	$(Q)-rm -f *.o *.dep *.a

-include $(SRCFILES:.cc=.dep)
