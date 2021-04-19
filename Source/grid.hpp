#ifndef _GRID_HPP_
#define _GRID_HPP_
/*
 * Utility functions for accessing array elements over
 * the sphere and the full three dimensional space,
 * and taking derivatives over the grid.
 */
#include <vector>

#include "cheb.hpp"
#include "sphere.hpp"

namespace Grid 
{
   void init(
         const double Rmax,
         const size_t nx,
         const size_t nlat,
         const size_t nphi);

   size_t indx(const size_t i_x, const size_t i_th, const size_t i_ph); 
   /*
    * returns spherical polar coordinate {r, phi, theta}
    */
   std::vector<double> r_th_ph(const size_t i_x, const size_t i_th, const size_t i_ph); 
   std::vector<double> r_th_ph(const size_t i); 
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
/*==========================================================================*/
/* Functions acting on full 3d grid functions */
/*==========================================================================*/
void set_partial_phi(const std::vector<double> &v, std::vector<double> &dv);
/*==========================================================================*/
void set_spherical_lap(const std::vector<double> &v, std::vector<double> &ddv);
/*==========================================================================*/
void set_sphereX(const std::vector<double> &v, std::vector<double> &vX);
/*==========================================================================*/
void set_partial_r(const std::vector<double> &v, std::vector<double> &dv);
/*==========================================================================*/
/* \partial_t f - p */
/*==========================================================================*/
double norm_indep_res(
      const double dt, 
      const std::vector<double> &f_n, 
      const std::vector<double> &f_np1, 
      const std::vector<double> &p); 
/*==========================================================================*/
double total_variation(const std::vector<double> &f); 
/*==========================================================================*/
/* Low pass filter in spectral space */
/*==========================================================================*/
void filter(std::vector<double> &v);
/*==========================================================================*/
}
#endif /* _GRID_HPP */
