#include <iostream>
#include <fstream>
#include <iomanip>
#include <cassert>

#include "io.hpp"
#include "arr.hpp"
#include "grid.hpp"
#include "params.hpp"

namespace Csv 
{
/*===========================================================================*/
namespace 
{
   std::vector<std::vector<double>> _grid_sphere_2d;
   std::vector<std::vector<double>> _grid_cart_3d;

   std::vector<std::string> _labels_sphere_2d = {"phi", "theta", "value"};
   std::vector<std::string> _labels_cart_3d   = {"x", "y", "z",  "value"};
}   
/*===========================================================================*/
void init()
{
   _grid_sphere_2d.resize(Params::nphi()*Params::nlat(),std::vector<double>(2,0));
   _grid_cart_3d.resize(Params::nx_nphi_nlat(),         std::vector<double>(3,0)); 

   for (size_t i=0; i<Params::nphi()*Params::nlat(); i++) {
      _grid_sphere_2d[i] = Grid::pt_sphere(i);
   }
   for (size_t i=0; i<Params::nx_nphi_nlat(); i++) {
      _grid_cart_3d[i] = Grid::pt_cart(i);
   }
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
      for (int i=0; i<n; i++) {
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
void write_sphere_2d(
      const std::string name, 
      const int itm,
      const size_t ix,
      const std::vector<double> &vals)
{
   std::string file_name= name+"_"+std::to_string(itm)+".csv";
   std::ofstream out;
   out.open(file_name,std::ios::app);

   std::vector<double> sphere_vals(Params::nphi()*Params::nlat(),0);
   Arr3d::get_row23(ix, vals, sphere_vals);

   if (out.is_open()) {

      const size_t indxs = _labels_sphere_2d.size();
      for (size_t i=0; i<indxs-1; i++) { 
         out<<_labels_sphere_2d[i]<<",";
      }		
      out<<_labels_sphere_2d[indxs-1]<<std::endl;

      const size_t n= sphere_vals.size();
      for (size_t i=0; i<n; i++) {
         for (size_t j=0; j<indxs-1; ++j) { /* grid point location */
            out<<std::setprecision(16)<<_grid_sphere_2d[i][j]<<",";
         }		
         out<<std::setprecision(16)<<sphere_vals[i]<<std::endl; /* grid point value */
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
void write_cart_3d(
      const std::string name, 
      const int itm, 
      const std::vector<double> &vals)
{
   std::string file_name= name+"_"+std::to_string(itm)+".csv";
   std::ofstream out;
   out.open(file_name,std::ios::app);

   if (out.is_open()) {

      const size_t indxs = _labels_cart_3d.size();
      for (size_t i=0; i<indxs-1; i++) { 
         out<<_labels_cart_3d[i]<<",";
      }		
      out<<_labels_cart_3d[indxs-1]<<std::endl;

      const size_t n= vals.size();
      for (size_t i=0; i<n; ++i) {
         for (size_t j=0; j<indxs-1; ++j) { /* grid point location */
            out<<std::setprecision(16)<<_grid_cart_3d[i][j]<<",";
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
