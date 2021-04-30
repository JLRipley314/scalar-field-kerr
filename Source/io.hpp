#ifndef _IO_HPP_
#define _IO_HPP_
/*
 * Utility functions for writing data to file.
 */

#include <string>
#include <vector>

#include "grid.hpp"
/*=========================================================================*/
namespace Csv 
/*=========================================================================*/
{
   void write(
         const std::string name, 
         const size_t itm, 
         const std::vector<double> &vals 
      );	
   void write_R_psl(
         const Grid &grid, 
         const std::string name, 
         const bool save_coords, 
         const int itm, 
         const std::vector<double> &vals
      );
#if USE_CHEB
   void write_n_psl(
         const Grid &grid, 
         const std::string name, 
         const bool save_coords, 
         const int itm, 
         const std::vector<double> &vals
      );
#endif
   void write_x_y_z(
         const Grid &grid, 
         const std::string name, 
         const bool save_coords, 
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
