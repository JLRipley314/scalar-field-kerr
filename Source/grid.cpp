#include "grid.hpp"

namespace Grid 
{
namespace {
   size_t _nx;
   size_t _nphi;
   size_t _nlat; 

   std::vector<std::vector<double>> _spherical_coords;
   std::vector<std::vector<double>> _polar_coords;
   std::vector<std::vector<double>> _cartesian_coords;
}
/*===========================================================================*/
void init()
{
   const size_t size = Params::nx_nphi_nlat();

   _nx   = Params::nx();
   _nphi = Params::nphi(); 
   _nlat = Params::nlat(); 

   _spherical_coords.resize(_nphi*_nlat, std::vector<double>(3,0));
   _polar_coords.resize(size, std::vector<double>(3,0));
   _cartesian_coords.resize(size, std::vector<double>(3,0));

   for (size_t ix=0; ix<_nx-1; ix++) { /* do not include ix=nx-1 as r=infty there */
   for (size_t ip=0; ip<_nphi; ip++) {
   for (size_t it=0; it<_nlat; it++) {
      double R     = Cheb::pt(ix);
      double phi   = Sphere::phi(  ip);
      double theta = Sphere::theta(it);  
      _polar_coords[Params::nphi()*Params::nlat()*ix + Sphere::indx_Sph(ip, it)] = {R, phi, theta};
      _cartesian_coords[Params::nphi()*Params::nlat()*ix + Sphere::indx_Sph(ip, it)] = {
         R*cos(phi)*sin(theta), 
         R*sin(phi)*sin(theta),
         R*cos(theta)
      };
   }
   }
   }
   for (size_t ip=0; ip<_nphi; ip++) {
   for (size_t it=0; it<_nlat; it++) {
      double phi   = Sphere::phi(  ip);
      double theta = Sphere::theta(it);  
      _spherical_coords[Sphere::indx_Sph(ip, it)] = {phi, theta};
   }
   }
}
/*===========================================================================*/
void cleanup()
{
}
/*===========================================================================*/
std::vector<double> pt_sphere(const size_t i_ph, const size_t i_th) 
{
   return _spherical_coords[_nlat*i_ph + i_th];
}
/*===========================================================================*/
std::vector<double> pt_sphere(const size_t i) 
{
   return _spherical_coords[i];
}
/*===========================================================================*/
std::vector<double> pt_polar(const size_t ix, const size_t i_ph, const size_t i_th) 
{
   return _polar_coords[_nphi*_nlat*ix + _nlat*i_ph + i_th];
}
/*===========================================================================*/
std::vector<double> pt_polar(const size_t i) 
{
   return _polar_coords[i];
}
/*===========================================================================*/
std::vector<double> pt_cart(const size_t ix, const size_t i_ph, const size_t i_th) 
{
   return _cartesian_coords[_nphi*_nlat*ix + _nlat*i_ph + i_th];
}
/*===========================================================================*/
std::vector<double> pt_cart(const size_t i) 
{
   return _cartesian_coords[i];
}
/*===========================================================================*/
} /* Grid */
