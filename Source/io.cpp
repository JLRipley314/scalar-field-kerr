#include <iostream>
#include <fstream>
#include <iomanip>
#include <cassert>

#include "io.hpp"
#include "grid.hpp"
#include "params.hpp"

namespace Csv 
{
/*===========================================================================*/
namespace 
{
   std::vector<std::vector<double>> _grid_R_th;
   std::vector<std::vector<double>> _grid_th_ph;
   std::vector<std::vector<double>> _grid_x_y_z;

   std::vector<std::string> _labels_R_th  = {"R",   "theta", "value"};
   std::vector<std::string> _labels_th_ph = {"theta", "phi", "value"};
   std::vector<std::string> _labels_x_y_z = {"x", "y", "z",  "value"};
}   
/*===========================================================================*/
void init()
{
   _grid_R_th.resize( Params::nx()  *Params::nlat(),std::vector<double>(2,0));
   _grid_th_ph.resize(Params::nphi()*Params::nlat(),std::vector<double>(2,0));
   _grid_x_y_z.resize(Params::nx_nlat_nphi(),       std::vector<double>(3,0)); 

   for (size_t i=0; i<Params::nx()*Params::nlat(); i++) {
      _grid_R_th[i] = Grid::R_th(i);
   }
   for (size_t i=0; i<Params::nphi()*Params::nlat(); i++) {
      _grid_th_ph[i] = Grid::th_ph(i);
   }
   for (size_t i=0; i<Params::nx_nlat_nphi(); i++) {
      _grid_x_y_z[i] = Grid::x_y_z(i);
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
void write_th_ph(
      const std::string name, 
      const int itm,
      const size_t ix,
      const std::vector<double> &vals)
{
   std::string file_name= name+"_"+std::to_string(itm)+".csv";
   std::ofstream out;
   out.open(file_name,std::ios::app);

   std::vector<double> th_ph_vals(Params::nphi()*Params::nlat(),0);
   Grid::get_row_th_ph(ix, vals, th_ph_vals);

   if (out.is_open()) {

      const size_t indxs = _labels_th_ph.size();
      for (size_t i=0; i<indxs-1; i++) { 
         out<<_labels_th_ph[i]<<",";
      }		
      out<<_labels_th_ph[indxs-1]<<std::endl;

      const size_t n= th_ph_vals.size();
      for (size_t i=0; i<n; i++) {
         for (size_t j=0; j<indxs-1; ++j) { /* grid point location */
            out<<std::setprecision(16)<<_grid_th_ph[i][j]<<",";
         }		
         out<<std::setprecision(16)<<th_ph_vals[i]<<std::endl; /* grid point value */
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
void write_x_y_z(
      const std::string name, 
      const int itm, 
      const std::vector<double> &vals)
{
   std::string file_name= name+"_"+std::to_string(itm)+".csv";
   std::ofstream out;
   out.open(file_name,std::ios::app);

   if (out.is_open()) {

      const size_t indxs = _labels_x_y_z.size();
      for (size_t i=0; i<indxs-1; i++) { 
         out<<_labels_x_y_z[i]<<",";
      }		
      out<<_labels_x_y_z[indxs-1]<<std::endl;

      const size_t n= vals.size();
      for (size_t i=0; i<n; ++i) {
         if (fabs(vals[i])>1e-3) {
            for (size_t j=0; j<indxs-1; ++j) { /* grid point location */
               out<<std::setprecision(16)<<_grid_x_y_z[i][j]<<",";
            }
            out<<std::setprecision(16)<<vals[i]<<std::endl; /* grid point value */
         }
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
