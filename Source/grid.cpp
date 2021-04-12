#include <cassert>
#include <iostream>
#include <iomanip>

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
   size_t _nlat; 
   size_t _nphi;

   std::vector<std::vector<double>> _r_th_ph;
   std::vector<std::vector<double>> _R_th_ph;
   std::vector<std::vector<double>> _x_y_z;
   std::vector<std::vector<double>> _th_ph;
   std::vector<std::vector<double>> _R_th;

   std::vector<double> _partial_R_to_partial_r;
}
/*===========================================================================*/
void init(
      const double cl,
      const size_t nx,
      const size_t nlat,
      const size_t nphi)
{
   _nx   = nx;
   _nlat = nlat; 
   _nphi = nphi; 

   _r_th_ph.resize(_nx*_nphi*_nlat, std::vector<double>(3,0));
   _R_th_ph.resize(_nx*_nphi*_nlat, std::vector<double>(3,0));
   _x_y_z.resize(  _nx*_nphi*_nlat, std::vector<double>(3,0));
   for (size_t ix=0; ix<_nx;   ix++) {
   for (size_t it=0; it<_nlat; it++) {
   for (size_t ip=0; ip<_nphi; ip++) {
      double R     = Cheb::pt(ix);
      double theta = Sphere::theta(it);  
      double phi   = Sphere::phi(  ip);

      _r_th_ph[INDX_R_TH_PH(ix,it,ip)][0] = (fabs(R-cl)>1e-16) ? (R/(1.0 - (R/cl))) : 1e12; 
      _r_th_ph[INDX_R_TH_PH(ix,it,ip)][1] = theta; 
      _r_th_ph[INDX_R_TH_PH(ix,it,ip)][2] = phi; 

      _R_th_ph[INDX_R_TH_PH(ix,it,ip)][0] = R;
      _R_th_ph[INDX_R_TH_PH(ix,it,ip)][1] = theta;
      _R_th_ph[INDX_R_TH_PH(ix,it,ip)][2] = phi;

      _x_y_z[INDX_R_TH_PH(ix,it,ip)][0] = R*cos(phi)*sin(theta); 
      _x_y_z[INDX_R_TH_PH(ix,it,ip)][1] = R*sin(phi)*sin(theta);
      _x_y_z[INDX_R_TH_PH(ix,it,ip)][2] = R*cos(theta);
   }
   }
   }
   _th_ph.resize(_nphi*_nlat, std::vector<double>(2,0));

   for (size_t it=0; it<_nlat; it++) {
   for (size_t ip=0; ip<_nphi; ip++) {
      double theta = Sphere::theta(it);  
      double phi   = Sphere::phi(  ip);
      _th_ph[INDX_TH_PH(it, ip)][0] = theta;
      _th_ph[INDX_TH_PH(it, ip)][1] = phi;
   }
   }
   _R_th.resize(_nx*_nlat, std::vector<double>(2,0));

   for (size_t ix=0; ix<_nx;   ix++) {
   for (size_t it=0; it<_nlat; it++) {
      double R     = Cheb::pt(ix);
      double theta = Sphere::theta(it);  
      _R_th[INDX_R_TH(ix,it)][0] = R;
      _R_th[INDX_R_TH(ix,it)][0] = theta;
   }
   }
   _partial_R_to_partial_r.resize(_nx,0);
   
   for (size_t ix=0; ix<_nx; ix++) {
      _partial_R_to_partial_r[ix] = pow(1.0 - (Cheb::pt(ix)/cl),2);
   }
}
/*=========================================================================*/
size_t indx(const size_t i_x, const size_t i_th, const size_t i_ph) 
{
   return INDX_R_TH_PH(i_x, i_th, i_ph);
}
/*===========================================================================*/
std::vector<double> r_th_ph(const size_t i_x, const size_t i_th, const size_t i_ph) 
{
   return _r_th_ph[INDX_R_TH_PH(i_x, i_th, i_ph)];
}
/*===========================================================================*/
std::vector<double> r_th_ph(const size_t i) 
{
   return _r_th_ph[i];
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
/*==========================================================================*/
void set_partial_phi(const std::vector<double> &v, std::vector<double> &dv)
{
   assert(v.size() ==_nx*_nlat*_nphi);
   assert(dv.size()==_nx*_nlat*_nphi);

#pragma omp parallel for
   for (size_t ix=0; ix<_nx; ix++) {
      std::vector<double> inter(   _nlat*_nphi);
      std::vector<double> inter_dv(_nlat*_nphi);

      get_row_th_ph(ix, v, inter); 

      Sphere::partial_phi(inter, inter_dv);

      set_row_th_ph(ix, inter_dv, dv); 
   }
}
/*==========================================================================*/
void set_spherical_lap(const std::vector<double> &v, std::vector<double> &ddv)
{
   assert(v.size()  ==_nx*_nlat*_nphi);
   assert(ddv.size()==_nx*_nlat*_nphi);

#pragma omp parallel for
   for (size_t ix=0; ix<_nx; ix++) {
      std::vector<double> inter(    _nlat*_nphi);
      std::vector<double> inter_ddv(_nlat*_nphi);

      get_row_th_ph(ix, v, inter); 

      Sphere::laplace_beltrami(inter, inter_ddv);

      set_row_th_ph(ix, inter_ddv, ddv); 
   }
}
/*==========================================================================*/
void set_partial_r(const std::vector<double> &v, std::vector<double> &dv)
{
   assert(v.size() ==_nx*_nlat*_nphi);
   assert(dv.size()==_nx*_nlat*_nphi);

#pragma omp parallel for
   for (size_t it=0; it<_nlat; it++) {
   for (size_t ip=0; ip<_nphi; ip++) {
      std::vector<double> inter(   _nx);
      std::vector<double> inter_dv(_nx);

      get_row_R(it, ip, v, inter); 

      Cheb::der(inter, inter_dv);

      for (size_t ix=0; ix<_nx; ix++) {
         inter_dv[ix] *= _partial_R_to_partial_r[ix];
      }

      set_row_R(it, ip, inter_dv, dv); 
   }
   }
}
/*==========================================================================*/
/* Low pass filter in spectral space */
/*==========================================================================*/
void filter(std::vector<double> &v)
{
   assert(v.size()==_nx*_nlat*_nphi);

#pragma omp parallel for
   for (size_t ix=0; ix<_nx; ix++) {
      std::vector<double> inter_sphere(_nlat*_nphi);
      get_row_th_ph(ix, v, inter_sphere); 
      Sphere::filter(inter_sphere);
      set_row_th_ph(ix, inter_sphere, v); 
   }
#pragma omp parallel for
   for (size_t it=0; it<_nlat; it++) {
   for (size_t ip=0; ip<_nphi; ip++) {
      std::vector<double> inter_radial(_nx);
      get_row_R(it, ip, v, inter_radial); 
      Cheb::filter(inter_radial);
      set_row_R(it, ip, inter_radial, v); 
   }
   }
}
/*==========================================================================*/
double norm_indep_res(
      const std::vector<double> &f, 
      const std::vector<double> &q) 
{
   assert(f.size()==_nx*_nlat*_nphi);
   assert(q.size()==_nx*_nlat*_nphi);

   std::vector<double> df(_nx*_nlat*_nphi,0);

   set_partial_r(f, df);

   double res = 0;

   for (size_t i=0; i<_nx*_nlat*_nphi; i++) {
      res += pow(df[i]-q[i],2);
   }
   res = pow(res/(_nx*_nlat*_nphi),0.5);

   return res;
}
/*==========================================================================*/
double total_variation(const std::vector<double> &v) 
{
   assert(v.size()==_nx*_nlat*_nphi);

   double tv = 0;
   for (size_t ix=0; ix<_nx-1;   ix++) {
   for (size_t it=0; it<_nlat-1; it++) {
   for (size_t ip=0; ip<_nphi-1; ip++) {
      tv += pow(
               pow(
                  v[indx(ix+1, it, ip)] 
               -  v[indx(ix,   it, ip)]
               ,2)
            +  pow(
                  v[indx(ix, it+1, ip)] 
               -  v[indx(ix, it,   ip)]
               ,2)
            +  pow(
                  v[indx(ix, it, ip+1)] 
               -  v[indx(ix, it, ip)]
               ,2)
            ,
            0.5); 
   }
   }
   }
   return tv;
}
/*===========================================================================*/
} /* Grid */

#undef INDX_R_TH
#undef INDX_R_PH
#undef INDX_R_TH_PH
#undef INDX_TH_PH
