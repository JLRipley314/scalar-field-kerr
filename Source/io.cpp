#include <iostream>
#include <fstream>
#include <iomanip>
#include <cassert>

#include "io.hpp"

namespace Csv 
{
/*===========================================================================*/
inline static bool exists(const std::string file_name)
{
   std::ifstream f(file_name);
   return f.good();
}
/*===========================================================================*/
void write(
      const std::string name, 
      const size_t itm, 
      const std::vector<double> &vals 
   )  
{
   std::string file_name= name+"_"+std::to_string(itm)+".csv";
   std::ofstream out;
   out.open(file_name,std::ios::app);

   if (out.is_open()) {
      const int n= vals.size();
      for (int i=0; i<n; ++i) {
         out<<std::setprecision(16)<<vals[i]<<std::endl;
      }
   }
   else {
      std::cout
         <<"ERROR(Csv::write): "+file_name+" does not exist"
         <<std::endl;
   }
   out.close();
}
/*===========================================================================*/
void write_unstructured(
      const std::string name, 
      const int itm, 
      const std::vector<std::vector<double>> &grid,
      const std::vector<double> &vals)
{
   std::string file_name= name+"_"+std::to_string(itm)+".csv";
   std::ofstream out;
   out.open(file_name,std::ios::app);

   if (out.is_open()) {
      const int indxs = grid[0].size();
      const int n= vals.size();
      for (int i=0; i<n; ++i) {
         for (int j=0; j<indxs; ++j) { /* grid point location */
            out<<std::setprecision(16)<<grid[i][j]<<",";
         }		
         out<<std::setprecision(16)<<vals[i]<<std::endl; /* grid point value */
      }		
   }
   else {
      std::cout
         <<"ERROR(Csv::write): "+file_name+" does not exist"
         <<std::endl;
   }
   out.close();
}
/*===========================================================================*/
}; /* Csv */
