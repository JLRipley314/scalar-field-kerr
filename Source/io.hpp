#include <stdio.h>
#include <stdlib.h>
#include <string>

/* Basic routines for writing/reading data to file */

namespace IO {
   void write_vect(std::string fn, double *vec, int N); 
   void write_mx(std::string fn,   double *mx,  int N1, int N2); 
}
