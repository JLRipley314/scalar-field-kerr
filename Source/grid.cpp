#include <cassert>
#include <iostream>

/*
 * NOTE: we use a COMPACTIFIED version of the radial coordinate
 * different from what we we use to evovle; i.e. NOT R,
 * but instead rr = r/(1+r/L) = 1/(1 + R/L)
 */

#include "grid.hpp"

#define INDX_R_TH(ix,it) ((_nlat)*(ix) + (it))
#define INDX_R_PH(ix,it) ((_nphi)*(ix) + (it))

#if SHTNS_CONTIGUOUS_LONGITUDES
#define INDX_R_TH_PH(ix,it,ip) ((_nphi)*(_nlat)*(ix) + (_nlat)*(ip) + (it))
#define INDX_TH_PH(it,ip) ((_nlat)*(ip) + (it))
#else 
#define INDX_R_TH_PH(ix,it,ip) ((_nlat)*(_nphi)*(ix) + (_nphi)*(it) + (ip))
#define INDX_TH_PH(it,ip) ((_nphi)*(it) + (ip))
#endif

/*===========================================================================*/
namespace Grid 
{
namespace {
   size_t _nx;
   size_t _nphi;
   size_t _nlat; 

   std::vector<std::vector<double>> _R_th_ph;
   std::vector<std::vector<double>> _x_y_z;
   std::vector<std::vector<double>> _th_ph;
   std::vector<std::vector<double>> _R_th;
}
/*===========================================================================*/
double compactification(const double cl, const double R)
{
   return cl/(1.0 + (R/cl));
}
/*===========================================================================*/
void init()
{
   _nx   = Params::nx();
   _nlat = Params::nlat(); 
   _nphi = Params::nphi(); 

   _R_th_ph.resize(_nx*_nphi*_nlat, std::vector<double>(3,0));
   _x_y_z.resize(  _nx*_nphi*_nlat, std::vector<double>(3,0));
   for (size_t ix=0; ix<_nx;   ix++) {
   for (size_t it=0; it<_nlat; it++) {
   for (size_t ip=0; ip<_nphi; ip++) {
      double R     = compactification(Params::cl(),Cheb::pt(ix));
      double theta = Sphere::theta(it);  
      double phi   = Sphere::phi(  ip);
      _R_th_ph[INDX_R_TH_PH(ix,it,ip)] = {
         R, 
         theta,
         phi 
      };
      _x_y_z[INDX_R_TH_PH(ix,it,ip)] = {
         R*cos(phi)*sin(theta), 
         R*sin(phi)*sin(theta),
         R*cos(theta)
      };
   }
   }
   }
   _th_ph.resize(_nphi*_nlat, std::vector<double>(2,0));

   for (size_t it=0; it<_nlat; it++) {
   for (size_t ip=0; ip<_nphi; ip++) {
      double theta = Sphere::theta(it);  
      double phi   = Sphere::phi(  ip);
      _th_ph[INDX_TH_PH(it, ip)] = {
         theta,
         phi
      };
   }
   }
   _R_th.resize(_nx*_nlat, std::vector<double>(2,0));

   for (size_t ix=0; ix<_nx;   ix++) {
   for (size_t it=0; it<_nlat; it++) {
      double R     = compactification(Params::cl(),Cheb::pt(ix));
      double theta = Sphere::theta(it);  
      _R_th[INDX_R_TH(ix,it)] = {
         R, 
         theta
      };
   }
   }
}
/*=========================================================================*/
size_t indx(const size_t i_x, const size_t i_th, const size_t i_ph) 
{
   return INDX_R_TH_PH(i_x, i_th, i_ph);
}
/*===========================================================================*/
std::vector<double> R_th_ph(const size_t i_x, const size_t i_th, const size_t i_ph) 
{
   return _R_th_ph[INDX_R_TH_PH(i_x, i_th, i_ph)];
}
/*===========================================================================*/
std::vector<double> R_th_ph(const size_t i) 
{
   return _R_th_ph[i];
}
/*===========================================================================*/
std::vector<double> x_y_z(const size_t i_x, const size_t i_th, const size_t i_ph) 
{
   return _x_y_z[INDX_R_TH_PH(i_x, i_th, i_ph)];
}
/*===========================================================================*/
std::vector<double> x_y_z(const size_t i) 
{
   return _x_y_z[i];
}
/*===========================================================================*/
std::vector<double> th_ph(const size_t i_th, const size_t i_ph) 
{
   return _th_ph[INDX_TH_PH(i_th, i_ph)];
}
/*===========================================================================*/
std::vector<double> th_ph(const size_t i) 
{
   return _th_ph[i];
}
/*===========================================================================*/
std::vector<double> R_th(const size_t i_x, const size_t i_th) 
{
   return _R_th[INDX_R_TH(i_x, i_th)];
}
/*===========================================================================*/
std::vector<double> R_th(const size_t i) 
{
   return _R_th[i];
}
/*=========================================================================*/
void get_row_R(const size_t it, const size_t ip, 
      const std::vector<double> &in,
      std::vector<double> &out) 
{
   assert(in.size() ==_nx*_nphi*_nlat);
   assert(out.size()==_nx);
   for (size_t ix=0; ix<_nx; ix++) {
      out[ix] = in[INDX_R_TH_PH(ix,it,ip)];
   }
} 
/*=========================================================================*/
void get_row_th(const size_t ix, const size_t ip, 
      const std::vector<double> &in,
      std::vector<double> &out) 
{
   assert(in.size() ==_nx*_nphi*_nlat);
   assert(out.size()==_nlat);
   for (size_t it=0; it<_nlat; it++) {
      out[it] = in[INDX_R_TH_PH(ix,it,ip)];
   }
} 
/*=========================================================================*/
void get_row_ph(const size_t ix, const size_t it, 
      const std::vector<double> &in,
      std::vector<double> &out) 
{
   assert(in.size() ==_nx*_nphi*_nlat);
   assert(out.size()==_nphi);
   for (size_t ip=0; ip<_nphi; ip++) {
      out[ip] = in[INDX_R_TH_PH(ix,it,ip)];
   }
} 
/*=========================================================================*/
void get_row_R_th(const size_t ip, 
      const std::vector<double> &in,
      std::vector<double> &out)
{
   assert(in.size() ==_nx*_nphi*_nlat);
   assert(out.size()==_nx*_nlat);
   for (size_t ix=0; ix<_nx;   ix++) {
   for (size_t it=0; it<_nlat; it++) {
      out[INDX_R_TH(ix,it)] = in[INDX_R_TH_PH(ix,it,ip)];
   }
   }
}
/*=========================================================================*/
void get_row_R_ph(const size_t it, 
      const std::vector<double> &in,
      std::vector<double> &out)
{
   assert(in.size() ==_nx*_nphi*_nlat);
   assert(out.size()==_nx*_nlat);
   for (size_t ix=0; ix<_nx;   ix++) {
   for (size_t ip=0; ip<_nphi; ip++) {
      out[INDX_R_PH(ix,ip)] = in[INDX_R_TH_PH(ix,it,ip)];
   }
   }
}
/*=========================================================================*/
void get_row_th_ph(const size_t ix, 
      const std::vector<double> &in,
      std::vector<double> &out)
{
   assert(in.size() ==_nx*_nphi*_nlat);
   assert(out.size()==_nphi*_nlat);
   for (size_t it=0; it<_nlat; it++) {
   for (size_t ip=0; ip<_nphi; ip++) {
      out[INDX_TH_PH(it,ip)] = in[INDX_R_TH_PH(ix,it,ip)];
   }
   }
}
/*=========================================================================*/
void set_row_R(const size_t it, const size_t ip, 
      const std::vector<double> &in,
      std::vector<double> &out) 
{
   assert(in.size() ==_nx);
   assert(out.size()==_nx*_nphi*_nlat);
   for (size_t ix=0; ix<_nx; ix++) {
      out[INDX_R_TH_PH(ix,it,ip)] = in[ix];
   }
} 
/*=========================================================================*/
void set_row_th(const size_t ix, const size_t ip, 
      const std::vector<double> &in,
      std::vector<double> &out) 
{
   assert(in.size() ==_nlat);
   assert(out.size()==_nx*_nphi*_nlat);
   for (size_t it=0; it<_nlat; it++) {
      out[INDX_R_TH_PH(ix,it,ip)] = in[it];
   }
} 
/*=========================================================================*/
void set_row_ph(const size_t ix, const size_t it, 
      const std::vector<double> &in,
      std::vector<double> &out) 
{
   assert(in.size() ==_nphi);
   assert(out.size()==_nx*_nphi*_nlat);
   for (size_t ip=0; ip<_nlat; ip++) {
      out[INDX_R_TH_PH(ix,it,ip)] = in[ip];
   }
} 
/*=========================================================================*/
void set_row_R_th(const size_t ip, 
      const std::vector<double> &in,
      std::vector<double> &out)
{
   assert(in.size() ==_nx*_nlat);
   assert(out.size()==_nx*_nphi*_nlat);
   for (size_t ix=0; ix<_nx;   ix++) {
   for (size_t it=0; it<_nlat; it++) {
      out[INDX_R_TH_PH(ix,it,ip)] = in[INDX_R_TH(ix,it)];
   }
   }
}
/*=========================================================================*/
void set_row_R_ph(const size_t it, 
      const std::vector<double> &in,
      std::vector<double> &out)
{
   assert(in.size() ==_nx*_nphi);
   assert(out.size()==_nx*_nphi*_nlat);
   for (size_t ix=0; ix<_nx;   ix++) {
   for (size_t ip=0; ip<_nphi; ip++) {
      out[INDX_R_TH_PH(ix,it,ip)] = in[INDX_R_PH(ix,ip)];
   }
   }
}
/*=========================================================================*/
void set_row_th_ph(const size_t ix, 
      const std::vector<double> &in,
      std::vector<double> &out)
{
   assert(in.size() ==_nphi*_nlat);
   assert(out.size()==_nx*_nphi*_nlat);
   for (size_t it=0; it<_nlat; it++) {
   for (size_t ip=0; ip<_nphi; ip++) {
      out[INDX_R_TH_PH(ix,it,ip)] = in[INDX_TH_PH(it,ip)];
   }
   }
}
/*===========================================================================*/
} /* Grid */

#undef INDX_R_TH
#undef INDX_R_PH
#undef INDX_R_TH_PH
#undef INDX_TH_PH
