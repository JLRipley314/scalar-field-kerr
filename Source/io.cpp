#include "io.hpp"

namespace IO {
/*==========================================================================*/
void write_vect(std::string fn, double *vec, int N)
{
   FILE *fp;
   int i;

   fp = fopen(fn.c_str(),"w");
   for (i=0;i<N;i++) {
      fprintf(fp,"%.6g ",vec[i]);
   }
   fclose(fp);
}
/*==========================================================================*/
void write_mx(std::string fn, double *mx, int N1, int N2)
{
   FILE *fp;
   int i,j;

   fp = fopen(fn.c_str(),"w");
   for (i=0;i<N1;i++) {
      for(j=0;j<N2;j++) {
         fprintf(fp,"%.6g ",mx[i*N2+j]);
      }
      fprintf(fp,"\n");
   }
   fclose(fp);
}

/*==========================================================================*/
} // IO
