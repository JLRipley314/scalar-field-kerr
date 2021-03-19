#==========================================================================
## Usage: compiles fortran code and links to FFTW
##
## Configured for two compiler options: ifort and gfortran
#==========================================================================
TOP= ./
SRC= $(TOP)Source/
BIN= $(TOP)Bin/

vpath %.mod $(BIN)
vpath %.f90 $(SRC)
#==========================================================================
## for FFTW 

INCFFTW= /usr/include
LIBFFTW= /lib/x86_64-linux-gnu
#==========================================================================
FC = gfortran#ifort#

FFLAGS= -g -fmax-errors=5 -O3 

SYSLIB= -lfftw3 
#==========================================================================
ifeq ($(FC),gfortran)
	FFLAGS+= -std=f2008 -Wall -Wextra -fimplicit-none -fopenmp \
		 -ftree-vectorize -march=native \
		 -J$(BIN) 
		#-fcheck=all  
endif

ifeq ($(FC),ifort)
	FFLAGS+= -std08 -ipo -warn declarations -warn all -qopenmp \
		-traceback \
		-module $(BIN) 
		#-check-bounds 
endif
#==========================================================================
MAIN = main.f90

OBJ = $(addprefix $(BIN), \
	mod_prec.o \
	mod_params.o \
	mod_field.o \
	mod_fields_list.o \
	mod_io.o \
	mod_cheb_fftw.o \
	mod_swal.o \
	mod_bkgrd.o \
	mod_ghp.o \
	mod_initial_data.o \
	mod_teuk.o \
	mod_write_level.o \
	)
#==========================================================================
all: default.run

%.run: $(MAIN) $(OBJ)
	$(FC) -o $(BIN)$@ $^ -I$(INCFFTW) -L$(LIBFFTW) $(SYSLIB) $(FFLAGS)   

$(BIN)%.o: %.f90
	$(FC) $(FFLAGS) -I$(INCFFTW) -L$(LIBFFTW) $(SYSLIB) -c -o $@ $^ 
#==========================================================================
clean_obj:
	$(RM) $(BIN)*.o
	$(RM) $(BIN)*.mod
clean_bin:
	$(RM) $(BIN)*.run
clean_out:
	$(RM) -r output/*
clean_all:
	$(RM) $(BIN)*.o
	$(RM) $(BIN)*.mod
	$(RM) $(BIN)*.run
	$(RM) -r output/*
