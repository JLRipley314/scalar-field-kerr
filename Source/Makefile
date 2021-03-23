CC = g++

BIN = ../Bin

SHTNS_INC = -I/usr/include -I/home/jripley/shtns/Library/include
SHTNS_LIB = -L/lib/x86_64-linux-gnu -L/home/jripley/shtns/Library/lib

#CFLAGS = -g -fmax-errors=5 -O2 -Wall -Wextra -fopenmp \
#	  -lshtns_omp -lfftw3_omp -lfftw3 -lm
CFLAGS = -g -fmax-errors=5 -O2 -Wall -Wextra \
	 -lshtns -lfftw3 -lm
#=============================================================================
OBJ = params.o io.o cheb.o sphere.o arr.o initial_data.o main.o

DEPS = params.hpp io.hpp cheb.hpp sphere.hpp initial_data.hpp arr.hpp
#=============================================================================
.PHONY: all
all: default.run

%.run: $(OBJ) 
	$(CC) -o $(BIN)/$@ $^ $(SHTNS_INC) $(SHTNS_LIB) $(CFLAGS)

%.o: %.cpp $(DEPS)
	$(CC) $(CFLAGS) $(SHTNS_INC) $(SHTNS_LIB) -c -o $@ $<
