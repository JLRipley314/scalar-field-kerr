CC = g++

BIN = ../Bin

SHTNS_INC = -I/usr/include -I/home/jripley/shtns/Library/include
SHTNS_LIB = -L/lib/x86_64-linux-gnu -L/home/jripley/shtns/Library/lib

CFLAGS= -g -fmax-errors=5 -O2 -Wall -Wextra -fopenmp \
	-lshtns_omp -lfftw3_omp -lfftw3 -lm
#=============================================================================
OBJ= $(addprefix $(OBJDIR), \
     field.o \
     main.o \
     )

DEPS=	
#=============================================================================
all: $(BIN)/default.run

%.run: $(OBJ) 
	$(CC) -o $@ $^ $(SHTNS_INC) $(SHTNS_LIB) $(CFLAGS)

$(BIN)%.o: %.cpp $(DEPS)
	$(CC) $(CFLAGS) -c -o $@ $<
