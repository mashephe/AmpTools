
PROJECT := Dalitz
PROJECT_ROOT := $(shell pwd)

INSTALL_BIN := $(PROJECT_ROOT)/bin
INSTALL_LIB := $(PROJECT_ROOT)/lib

export

.PHONY: default clean library

default: library exe

library:
	@echo "=== Building $(PROJECT) libraries ==="
	@$(MAKE) -C $(PROJECT)Lib
	
exe: library
	@echo "=== Building $(PROJECT) executables ==="
	@$(MAKE) -C $(PROJECT)Exe

clean:
	@$(MAKE) -C $(PROJECT)Lib clean
	@$(MAKE) -C $(PROJECT)Exe clean
