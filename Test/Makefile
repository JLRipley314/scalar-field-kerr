CC = g++

BIN = ../Bin
SRC = ../Source

# googletest
GT_LIB = -lgtest -lpthread
GT_INCLUDE = -I$(HOME)/googletest/build/lib

# shtns 
SHTNS_INC = -I$(HOME)/shtns/Library/include
SHTNS_LIB = -L$(HOME)/shtns/Library/lib

CFLAGS = -g -fmax-errors=5 -O2 -Wall -Wextra -fopenmp \
	 -lshtns -lfftw3 -lm #-D_GLIBCXX_DEBUG
#=============================================================================
OBJ = main.o \
      test_finite_diff.o $(SRC)/finite_diff.o \
      test_cheb.o $(SRC)/cheb.o \
      test_sphere.o $(SRC)/sphere.o \
      test_grid.o $(SRC)/cheb.o $(SRC)/sphere.o $(SRC)/grid.o
#=============================================================================
.PHONY: all
all: test.run

.PHONY: run_test
run:
	../Bin/test.run | tee test.log

%.run: $(OBJ) 
	$(CC) -o $(BIN)/$@ $^ $(GT_INCLUDE) $(GT_LIB) $(SHTNS_INC) $(SHTNS_LIB) $(CFLAGS)

%.o: %.cpp
	$(CC) $(CFLAGS) $(GT_INCLUDE) $(GT_LIB) $(SHTNS_INC) $(SHTNS_LIB) -c -o $@ $<

$(SRC)/%.o: %.cpp
	$(CC) $(CFLAGS) $(GT_INCLUDE) $(GT_LIB) $(SHTNS_INC) $(SHTNS_LIB) -c -o $@ $<
