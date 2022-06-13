
PROJECT := Dalitz
PROJECT_HOME := $(shell pwd)

INSTALL_BIN := $(PROJECT_HOME)/bin
INSTALL_LIB := $(PROJECT_HOME)/lib

export

.PHONY: default clean library force

default: library exe

library: compiler_flags
	@echo "=== Building $(PROJECT) libraries ==="
	@$(MAKE) -C $(PROJECT)Lib
	
exe: library
	@echo "=== Building $(PROJECT) executables ==="
	@$(MAKE) -C $(PROJECT)Exe

clean:
	@$(MAKE) -C $(PROJECT)Lib clean
	@$(MAKE) -C $(PROJECT)Exe clean

compiler_flags: force
	@echo '$(CXX_FLAGS)' | cmp -s - $@ || echo '$(CXX_FLAGS)' > $@
