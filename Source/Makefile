#==========================================================================
## Usage: compiles executable 
##
## Configured for gnu and intel c++ compilers 
#==========================================================================
TOP= ./
BIN= ../Bin/
#==========================================================================
## for FFTW and shtns libraries

SHTNS_INC = -I/usr/include -I/home/jripley/shtns/library/include
SHTNS_LIB = -L/lib/x86_64-linux-gnu -L/home/jripley/shtns/library/lib 

SYS_LIB = -lshtns -lfftw3 -lm
#==========================================================================
CC = g++

CFLAGS= -g -fmax-errors=5 -O2 -std=c++14 -Wall -Wextra -fopenmp
#==========================================================================
MAIN = main.cpp

OBJS= $(addprefix $(BIN), \
	io.o \
	spherical_harmonics.o \
	)
#==========================================================================
all: $(BIN)default.run

$(BIN)%.run: $(MAIN) $(OBJS)
	$(CC) -o $@ $^ $(SHTNS_INC) $(SHTNS_LIB) $(SYS_LIB) $(CFLAGS)

$(BIN)%.o: %.cpp
	$(CC) $(CFLAGS) $(SHTNS_INC) $(SHTNS_LIB) $(SYS_LIB) -c -o $@ $^ 
