CC = g++

BIN = ../Bin
SRC = ../Source

GT_LIB = -lgtest -lpthread
GT_INCLUDE = -I/usr/local/include

SHTNS_INC = -I/usr/include -I/home/jripley/shtns/Library/include
SHTNS_LIB = -L/lib/x86_64-linux-gnu -L/home/jripley/shtns/Library/lib

CFLAGS = -g -fmax-errors=5 -O2 -Wall -Wextra \
	 -lshtns -lfftw3 -lm \

vpath $(SRC)
#=============================================================================
OBJ = main.o \
      test_cheb.o $(SRC)/cheb.o \
      test_sphere.o $(SRC)/sphere.o
#=============================================================================
.PHONY: all
all: test.run

%.run: $(OBJ) 
	$(CC) -o $(BIN)/$@ $^ $(GT_INCLUDE) $(GT_LIB) $(SHTNS_INC) $(SHTNS_LIB) $(CFLAGS)

%.o: %.cpp
	$(CC) $(CFLAGS) $(GT_INCLUDE) $(GT_LIB) $(SHTNS_INC) $(SHTNS_LIB) -c -o $@ $<

$(SRC)/%.o: %.cpp
	$(CC) $(CFLAGS) $(GT_INCLUDE) $(GT_LIB) $(SHTNS_INC) $(SHTNS_LIB) -c -o $@ $<
