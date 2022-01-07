
# try to catch common previous uses that will bungle things up...

ifdef MPI
$(error Use the mpi target 'make mpi' rather than setting MPI=1.)
endif

ifdef GPU
$(error Use the mpi target 'make gpu' rather than setting GPU=1.)
endif

# check for verbose output
ifdef VERBOSE
	Q :=
	vecho = @true
else
	Q := @
	vecho = @echo
endif

export

.PHONY: default, all, mpi, gpu, mpigpu, clean

default:
	$(vecho) "=== Building AmpTools ==="
	$(Q)$(MAKE) -C AmpTools
	$(vecho) "=== Building AmpPlotter ==="
	$(Q)$(MAKE) -C AmpPlotter
	$(vecho) "=== Building Dalitz tutorial ==="
	$(Q)$(MAKE) -C Tutorials/Dalitz

mpi:  default
	$(vecho) "=== Building AmpTools with MPI ==="
	$(Q)$(MAKE) -C AmpTools MPI=1
	$(vecho) "=== Building Dalitz tutorial with MPI ==="
	$(Q)$(MAKE) -C Tutorials/Dalitz MPI=1
	
gpu:
	$(vecho) "=== Building AmpTools with GPU acceleration ==="
	$(Q)$(MAKE) -C AmpTools GPU=1
	$(vecho) "=== Building AmpPlotter ==="
	$(Q)$(MAKE) -C AmpPlotter
	$(vecho) "=== Building Dalitz tutorial with GPU acceleration ==="
	$(Q)$(MAKE) -C Tutorials/Dalitz GPU=1

mpigpu: gpu
	$(vecho) "=== Building AmpTools with MPI and GPU acceleration ==="
	$(Q)$(MAKE) -C AmpTools GPU=1 MPI=1
	$(vecho) "=== Building AmpPlotter ==="
	$(Q)$(MAKE) -C AmpPlotter
	$(vecho) "=== Building Dalitz tutorial with MPI and GPU acceleration ==="
	$(Q)$(MAKE) -C Tutorials/Dalitz MPI=1 GPU=1

clean:
	$(Q)$(MAKE) -C AmpTools clean
	$(Q)$(MAKE) -C AmpPlotter clean
	$(Q)$(MAKE) -C Tutorials/Dalitz clean

