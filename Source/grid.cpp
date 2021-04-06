#include "grid.hpp"

namespace Grid 
{
namespace {
   size_t _nx;
   size_t _nphi;
   size_t _nlat; 

   std::vector<std::vector<double>> _g3d;
}
/*===========================================================================*/
void init()
{
   const size_t size = Params::nx_nphi_nlat();

   _nx   = Params::nx();
   _nphi = Params::nphi(); 
   _nlat = Params::nlat(); 

   _g3d.resize(size, std::vector<double>(3,0));

   for (size_t ix=0; ix<_nx-1; ix++) { /* do not include ix=nx-1 as r=infty there */
   for (size_t ip=0; ip<_nphi; ip++) {
   for (size_t it=0; it<_nlat; it++) {
      double R     = Cheb::pt(ix);
      double phi   = Sphere::phi(  ip);
      double theta = Sphere::theta(it);  
      _g3d[Params::nphi()*Params::nlat()*ix + Sphere::indx_Sph(ip, it)] = {R, phi, theta};
   }
   }
   }
}
/*===========================================================================*/
void cleanup()
{
}
/*===========================================================================*/
std::vector<double> pt(const size_t ix, const size_t i_ph, const size_t i_th) 
{
   return _g3d[_nphi*_nlat*ix + _nlat*i_ph + i_th];
}
/*===========================================================================*/
std::vector<double> pt(const size_t i) 
{
   return _g3d[i];
}
/*===========================================================================*/
} /* Grid */
