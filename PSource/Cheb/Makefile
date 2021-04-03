##
## Make a shared library for python bindings
##
CC = gcc

CFLAGS = -g -fmax-errors=5 -O2 -Wall -Wextra \
	 -lfftw3 -lm

SO = -shared -Wl,-soname,lib_cheb
FPIC = -fPIC

all: lib_cheb.so

%.o: %.c 
	$(CC) $(CFLAGS) -c -o $@ $<

%.so: %.o
	$(CC) $(SO) -o $@ $(FPIC) $< 
