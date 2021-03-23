#include <iostream>
#include <fstream>
#include <iomanip>
#include <cassert>
#include <string>

#include "io.hpp"

namespace Csv {
/*===========================================================================*/
inline static bool exists(const std::string file_name)
{
   std::ifstream f(file_name);
   return f.good();
}
/*===========================================================================*/
void write(const std::string outdir, const double time, const Field &f)
{
   std::string file_name= outdir+"/"+f.name+".csv";
   std::ofstream out;
   out.open(file_name,std::ios::app);

   if (out.is_open()) {
      const int n= f.np1.size();
      out<<time<<","<<n<<",";
      for (int i=0; i<n-1; ++i) {
         out<<std::setprecision(16)<<f.np1[i]<<",";
      }		
      out<<f.np1[n-1]<<endl;
   }
   else {
      std::cout
         <<"ERROR(Csv::write): "+file_name+" does not exist"
         <<endl;
   }
   out.close();
}
/*===========================================================================*/
}; /* Csv */
