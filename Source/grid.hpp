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
   void cleanup();
   /*
    * returns spherical polar coordinate {R, phi, theta}
    */
   std::vector<double> pt_polar(const size_t ix, const size_t i_ph, const size_t i_th); 
   std::vector<double> pt_polar(const size_t i); 
   /*
    * returns cartesian coordinate {x, y, z}
    */
   std::vector<double> pt_cart(const size_t ix, const size_t i_ph, const size_t i_th); 
   std::vector<double> pt_cart(const size_t i); 
}
#endif /* _GRID_HPP */
