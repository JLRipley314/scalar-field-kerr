#include <iostream>
#include <fstream>
#include <iomanip>
#include <cassert>

#include "io.hpp"
#include "grid.hpp"

namespace Csv 
{
/*===========================================================================*/
namespace 
{
   std::vector<std::string> _labels_R_l   = {"R",   "l",     "log(value)"};
   std::vector<std::string> _labels_n_l   = {"n",   "l",     "log(value)"};
   std::vector<std::string> _labels_x_y_z = {"x", "y", "z",  "value"};
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
/* f_{n,l}, with l = sqrt(Sum_m f_{l,m}) */
/*===========================================================================*/
void write_n_psl(
      const Grid &grid,
      const std::string name, 
      const bool save_coords, 
      const int itm,
      const std::vector<double> &vals)
{
   std::string file_name= name+"_NPsl_"+std::to_string(itm)+".csv";
   std::ofstream out;
   out.open(file_name,std::ios::app);

   std::vector<double> power_spectrum(grid.nx()*grid.nl(),0);
   grid.set_n_l_coef(vals, power_spectrum);

   if (out.is_open()) {

      const size_t indxs = _labels_n_l.size();
      if (save_coords) {
         for (size_t i=0; i<indxs-1; i++) { 
            out<<_labels_n_l[i]<<",";
         }		
         out<<_labels_n_l[indxs-1]<<std::endl;
      }
      const size_t n= power_spectrum.size();
      for (size_t i=0; i<n; i++) {
         for (size_t j=0; j<indxs-1; ++j) { /* grid point location */
            out<<std::setprecision(16)<<grid.n_l(i)[j]<<",";
         }		
         out<<std::setprecision(16)
            <<((power_spectrum[i]<1e-16) ? -16 : log(power_spectrum[i])/log(10))
            <<std::endl; 
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
void write_R_psl(
      const Grid &grid,
      const std::string name, 
      const bool save_coords, 
      const int itm,
      const std::vector<double> &vals)
{
   std::string file_name= name+"_RPsl_"+std::to_string(itm)+".csv";
   std::ofstream out;
   out.open(file_name,std::ios::app);

   std::vector<double> power_spectrum(grid.nx()*grid.nl(),0);
   grid.set_angular_power_spectrum(vals, power_spectrum);

   if (out.is_open()) {

      const size_t indxs = _labels_R_l.size();
      if (save_coords) {
         for (size_t i=0; i<indxs-1; i++) { 
            out<<_labels_R_l[i]<<",";
         }		
         out<<_labels_R_l[indxs-1]<<std::endl;
      }
      const size_t n= power_spectrum.size();
      for (size_t i=0; i<n; i++) {
         for (size_t j=0; j<indxs-1; ++j) { /* grid point location */
            out<<std::setprecision(16)<<grid.R_l(i)[j]<<",";
         }		
         out<<std::setprecision(16)
            <<((power_spectrum[i]<1e-16) ? -16 : log(power_spectrum[i])/log(10))
            <<std::endl; 
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
inline double magnitude(const std::vector<double> &v)
{
   double mag = 0;
   for (const double &val: v) {
      mag += pow(val,2);
   }
   return pow(mag,0.5);
}
/*===========================================================================*/
void write_x_y_z(
      const Grid &grid, 
      const std::string name, 
      const bool save_coords, 
      const int itm,
      const double vmin,
      const double vmax,
      const double rmin,
      const double rmax, 
      const std::vector<double> &vals)
{
   std::string file_name= name+"_XYZ_"+std::to_string(itm)+".csv";
   std::ofstream out;
   out.open(file_name,std::ios::app);

   if (out.is_open()) {

      const size_t indxs = _labels_x_y_z.size();
      if (save_coords) {
         for (size_t i=0; i<indxs-1; i++) { 
            out<<_labels_x_y_z[i]<<",";
         }		
         out<<_labels_x_y_z[indxs-1]<<std::endl;
      }
      const size_t n= vals.size();
      for (size_t i=0; i<n; ++i) {
         if ((fabs(vals[i])>vmin)
         &&  (fabs(vals[i])<vmax)
         &&  (magnitude(grid.x_y_z(i))>rmin)
         &&  (magnitude(grid.x_y_z(i))<rmax)
         ) {
            for (size_t j=0; j<indxs-1; ++j) { /* grid point location */
               out<<std::setprecision(16)<<grid.x_y_z(i)[j]<<",";
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
