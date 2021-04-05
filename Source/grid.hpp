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
   namespace {
      size_t _nx;
      size_t _nphi;
      size_t _nlat; 
   }
   std::vector<std::vector<double>> g3d;

   void init()
   {
      const size_t size = Params::nx_nphi_nlat();

      _nx   = Params::nx();
      _nphi = Params::nphi(); 
      _nlat = Params::nlat(); 

      g3d.resize(size, std::vector<double>(3,0));

      for (size_t ix=0; ix<_nx-1; ix++) { /* do not include ix=nx-1 as r=infty there */
      for (size_t ip=0; ip<_nphi; ip++) {
      for (size_t it=0; it<_nlat; it++) {
         double r     = pow(Params::cl(),2)/Cheb::pt(ix);
         double phi   = Sphere::phi(  ip);
         double theta = Sphere::theta(it);  
         g3d[Params::nphi()*Params::nlat()*ix + Sphere::indx_Sph(ip, it)] = {r, phi, theta};
      }
      }
      }
   }
   /*
    * returns {r, phi, theta}
    */
   inline std::vector<double> pt(const size_t ix, const size_t i_ph, const size_t i_th) 
   {
      return g3d[_nphi*_nlat*ix + _nlat*i_ph + i_th];
   }

};
#endif /* _GRID_HPP */
