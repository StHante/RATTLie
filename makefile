# Which targets aren't actual filenames
.PHONY: default cleanup clean

# Fortran compiler:
FC = gfortran
# Flags for Fortran compiler
FFLAGS = -cpp -ffree-line-length-none -DNOZ
#-DUSE_INDEX_1_APPROX

ifdef EXTRAFFLAGS
FFLAGS := $(FFLAGS) $(EXTRAFFLAGS)
else
FFLAGS := $(FFLAGS) -O3
endif

#FFLAGS := $(FFLAGS) -Dpure='' -g -fbounds-check -fimplicit-none -fbacktrace -fcheck=all -finit-real=-inf -finit-integer=-77 -ffpe-trap=zero,invalid
#FFLAGS := $(FFLAGS) -pg
#FFLAGS := $(FFLAGS) -DNOBANDEDDER
#FFLAGS := $(FFLAGS) -DSTNUM
#FFLAGS := $(FFLAGS) -DRELAXATION

# Set flags used by Fortran
ifeq "$(FC)" "gfortran"
	INCLUDE = -I
	MODULEP = -J
else
ifeq "$(FC)" "ifort"
	INCLUDE = -I
	MODULEP = -module
endif
endif


############################################################# TARGETS ##
# default: This target will be built if make is run without arguments
default: obj/RATTLie.o makefile


# Target for all objects built from preprocessed Fortran90 files (.F90)
obj/%.o: src/%.F90 makefile
	$(FC) -c $(MODULEP) obj $(FFLAGS) -o $@ $<


# Target to clean up the object folder
cleanup:
	-rm obj/*


# Target that cleans a bit more
clean: cleanup
