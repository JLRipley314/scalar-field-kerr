CC = g++

BIN = ../Bin

SHTNS_INC = -I/usr/include -I/home/jripley/shtns/Library/include
SHTNS_LIB = -L/lib/x86_64-linux-gnu -L/home/jripley/shtns/Library/lib

CFLAGS = -g -fmax-errors=5 -O2 -Wall -Wextra -fopenmp \
	 -lshtns -lfftw3 -lm
#CFLAGS = -g -fmax-errors=5 -O2 -Wall -Wextra \
#	 -lshtns -lfftw3 -lm #-D_GLIBCXX_DEBUG
#=============================================================================
OBJ = params.o \
      cheb.o \
      sphere.o \
      grid.o \
      io.o \
      field.o \
      initial_data.o \
      scalar_eom.o \
      main.o

DEPS = params.hpp \
       cheb.hpp \
       sphere.hpp \
       grid.hpp \
       io.hpp \
       field.hpp \
       initial_data.hpp \
       scalar_eom.hpp
#=============================================================================
.PHONY: all
all: default.run

%.run: $(OBJ) 
	$(CC) -o $(BIN)/$@ $^ $(SHTNS_INC) $(SHTNS_LIB) $(CFLAGS)

%.o: %.cpp $(DEPS)
	$(CC) $(CFLAGS) $(SHTNS_INC) $(SHTNS_LIB) -c -o $@ $<
