#ifndef _IO_HPP_
#define _IO_HPP_
/*
 * Utility functions for writing data to file.
 */

#include <string>
#include <vector>
/*=========================================================================*/
namespace Csv 
/*=========================================================================*/
{
   void init();

   void write(
         const std::string name, 
         const size_t itm, 
         const std::vector<double> &vals 
      );	
   void write_th_ph(
         const std::string name, 
         const int itm, 
         const size_t ix, 
         const std::vector<double> &vals
      );
   void write_x_y_z(
         const std::string name, 
         const int itm,
         const double vmin,
         const double vmax,
         const double rmin,
         const double rmax, 
         const std::vector<double> &vals
      );
/*=========================================================================*/
}; /* Csv */
/*=========================================================================*/
#endif /* _IO_HPP_ */
