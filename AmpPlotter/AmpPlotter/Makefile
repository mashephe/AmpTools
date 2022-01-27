SRCFILES := $(wildcard *.cc)

.PHONY: default clean
.PRECIOUS: .%.o

default: $(LIB)

%.a:  $(SRCFILES:%.cc=.%.o)
	$(vecho) "--> Archiving $@"
	$(Q)ar -rsc $@ $^

.%.o: %.cc .%.d
	$(vecho) "-> Compiling $<"
	$(Q)$(CXX) $(CXX_FLAGS) -M -MP -MT $@ -o .$*.d $< $(INC_DIR)
	$(Q)$(CXX) $(CXX_FLAGS) -c -o $@ $< $(INC_DIR)

DEPFILES := $(SRCFILES:%.cc=.%.d)
$(DEPFILES):

clean:
	$(Q)-rm -f .*.o .*.d *.a

-include $(SRCFILES:.cc=.dep)
