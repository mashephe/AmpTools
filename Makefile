
# try to catch common previous uses that will bungle things up...

ifdef MPI
$(error Use the mpi target 'make mpi' rather than setting MPI=1.)
endif

ifdef GPU
$(error Use the mpi target 'make gpu' rather than setting GPU=1.)
endif

AMPTOOLS_HOME := $(shell pwd)
AMPTOOLS := $(AMPTOOLS_HOME)/AmpTools
AMPPLOTTER := $(AMPTOOLS_HOME)/AmpPlotter

export

.PHONY: default, all, mpi, gpu, mpigpu, gpumpi, clean

default:
	@echo "=== Building AmpTools ==="
	@$(MAKE) -C AmpTools
	@echo "=== Building AmpPlotter ==="
	@$(MAKE) -C AmpPlotter
	@echo "=== Building Dalitz tutorial ==="
	@$(MAKE) -C Tutorials/Dalitz

mpi:  default
	@echo "=== Building AmpTools with MPI ==="
	@$(MAKE) -C AmpTools MPI=1
	@echo "=== Building Dalitz tutorial with MPI ==="
	@$(MAKE) -C Tutorials/Dalitz MPI=1

gpu:
	@echo "=== Building AmpTools with GPU acceleration ==="
	@$(MAKE) -C AmpTools GPU=1
	@echo "=== Building AmpPlotter ==="
	@$(MAKE) -C AmpPlotter
	@echo "=== Building Dalitz tutorial with GPU acceleration ==="
	@$(MAKE) -C Tutorials/Dalitz GPU=1

mpigpu: gpu
	@echo "=== Building AmpTools with MPI and GPU acceleration ==="
	@$(MAKE) -C AmpTools GPU=1 MPI=1
	@echo "=== Building AmpPlotter ==="
	@$(MAKE) -C AmpPlotter
	@echo "=== Building Dalitz tutorial with MPI and GPU acceleration ==="
	@$(MAKE) -C Tutorials/Dalitz MPI=1 GPU=1

gpumpi: mpigpu

clean:
	@$(MAKE) -C AmpTools clean
	@$(MAKE) -C AmpPlotter clean
	@$(MAKE) -C Tutorials/Dalitz clean

