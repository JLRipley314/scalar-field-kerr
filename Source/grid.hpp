#ifndef _GRID_HPP_
#define _GRID_HPP_
/*
 * Header file with 
 * utility functions for accessing array elements over
 * the sphere and the full three dimensional space
 */
#include <vector>

#include "cheb.hpp"
#include "params.hpp"
#include "sphere.hpp"

namespace Grid 
{
   void init();

   size_t indx(const size_t i_x, const size_t i_th, const size_t i_ph); 
   /*
    * returns spherical polar coordinate {R, phi, theta}
    */
   std::vector<double> R_th_ph(const size_t i_x, const size_t i_th, const size_t i_ph); 
   std::vector<double> R_th_ph(const size_t i); 
   /*
    * returns cartesian coordinate {x, y, z}
    */
   std::vector<double> x_y_z(const size_t i_x, const size_t i_th, const size_t i_ph); 
   std::vector<double> x_y_z(const size_t i); 
   /*
    * returns spherical coordinate {phi, theta}
    */
   std::vector<double> th_ph(const size_t i_th, const size_t i_ph); 
   std::vector<double> th_ph(const size_t i); 
   /*
    * returns spherical coordinate {R, theta}
    */
   std::vector<double> R_th(const size_t i_x, const size_t i_th); 
   std::vector<double> R_th(const size_t i); 
   /*
    * Get row vals
    */
   void get_row_R(const size_t j, const size_t k, 
         const std::vector<double> &in,
         std::vector<double> &out
      ); 
   void get_row_th(const size_t i, const size_t k, 
         const std::vector<double> &in,
         std::vector<double> &out
      ); 
   void get_row_ph(const size_t i, const size_t j, 
         const std::vector<double> &in,
         std::vector<double> &out
      ); 
   void get_row_R_th(const size_t k, 
         const std::vector<double> &in,
         std::vector<double> &out
      );
   void get_row_R_ph(const size_t j, 
         const std::vector<double> &in,
         std::vector<double> &out
      );
   void get_row_th_ph(const size_t i, 
         const std::vector<double> &in,
         std::vector<double> &out
      );
   /*
    * Get row vals
    */
   void set_row_R(const size_t j, const size_t k, 
         const std::vector<double> &in,
         std::vector<double> &out
      ); 
   void set_row_th(const size_t i, const size_t k, 
         const std::vector<double> &in,
         std::vector<double> &out
      ); 
   void set_row_ph(const size_t i, const size_t j, 
         const std::vector<double> &in,
         std::vector<double> &out
      ); 
   void set_row_R_th(const size_t k, 
         const std::vector<double> &in,
         std::vector<double> &out
      );
   void set_row_R_ph(const size_t j, 
         const std::vector<double> &in,
         std::vector<double> &out
      );
   void set_row_th_ph(const size_t i, 
         const std::vector<double> &in,
         std::vector<double> &out
      );
}
#endif /* _GRID_HPP */
