# Installation notes

Please feel free to contact me if you are having any problems:
lloydripley[at]gmail[dot]com 

I have only tested everything on *nix machines, using the GNU C++
compiler (version 4.8.5 onwards, C++11 standard).

## Installing FFTW

FFTW is an open source library that computes fast Fourier transforms.
It can be found at http://www.fftw.org/

I assume that you have have FFTW3 installed.
If not, and if you are using a home computer, then you should be able
to install FFTW using a package manager.
If you are on a computer cluster, then it should be available
using `module load [...]`

## Installing shtns

shtns is an open source library that can perform scalar and vector
spherical harmonic transforms. It is very fast, and 
has been extensively used and tested in the context
of geophysical fluid simulations.

shtns can be found at https://nschaeff.bitbucket.io/shtns/

I assume that when you install shtns, that you will make
a new directory within shtns called `Library`, and
install shtns using use
`./configure --prefix=$(HOME)/shtns/Library`

## Installing googletest 

googletest is an open source unit test framework for C++.

googletest cane be found at https://github.com/google/googletest

The internal googletest README.md file should contain all
the necessary installation instructions.

The following stackoverflow Q&A may be useful as well.

https://stackoverflow.com/questions/13513905/how-to-set-up-googletest-as-a-shared-library-on-linux

## Making scalar-field-kerr binary

Once all the above libraries are installed, you should be able to simply
type `make` to compile a binary `default.run` in `Bin/`.	
