CC = g++ -std=c++11

BIN = ../Bin

SHTNS_INC = -I$(HOME)/shtns/Library/include
SHTNS_LIB = -L$(HOME)/shtns/Library/lib

CFLAGS = 

ifeq ($(debug),true) 
CFLAGS = -pg -fbounds-check
else
CFLAGS = -g
endif

CFLAGS += -fmax-errors=5 -O3 -Wall -Wextra -fopenmp -lshtns -lfftw3 -lm

ifeq ($(use_cheb),true) 
CFLAGS += -DUSE_CHEB 
BIN_NAME = default_cheb.run
else
BIN_NAME = default_fd.run
endif

ifeq ($(use_contiguous_longitudes),true) 
CFLAGS += -DSHTNS_CONTIGUOUS_LONGITUDES
endif
#=============================================================================
.PHONY: all
all: $(BIN_NAME) 
#=============================================================================
OBJ = params.o \
      cheb.o \
      finite_diff.o \
      finite_diff.o \
      sphere.o \
      grid.o \
      io.o \
      field.o \
      initial_data.o \
      scalar_eom.o \
      unit_manager.o \
      main.o

DEPS = params.hpp \
       cheb.hpp \
       finite_diff.hpp \
       finite_diff.hpp \
       sphere.hpp \
       grid.hpp \
       io.hpp \
       field.hpp \
       initial_data.hpp \
       scalar_eom.hpp \
       unit_manager.hpp
#=============================================================================
%.run: $(OBJ) 
	$(CC) -o $(BIN)/$@ $^ $(SHTNS_INC) $(SHTNS_LIB) $(CFLAGS)

%.o: %.cpp $(DEPS)
	$(CC) $(CFLAGS) $(SHTNS_INC) $(SHTNS_LIB) -c -o $@ $<

#=============================================================================
.PHONY: clean
clean:
	$(RM) *.o
