#ifndef GRID_HPP__
#define GRID_HPP__
/*
 * Utility functions for accessing array elements over
 * the sphere and the full three dimensional space,
 * and taking derivatives over the grid.
 */
#include <cassert>
#include <vector>

#include "cheb.hpp"
#include "finite_diff.hpp"
#include "sphere.hpp"

class Grid
{
public:
   Grid(
      const double cl,
      const double Rmin,
      const double Rmax,
      const size_t nx,
      const size_t nl,
      const size_t nm,
      const size_t nlat,
      const size_t nphi
      );
   ~Grid();

   std::vector<double> compute_ylm(const size_t l_ang, const size_t m_ang) const;
   /* Functions acting on full 3d grid functions 
    */
   void set_partial_phi(  const std::vector<double> &v, std::vector<double> &dv) const;
   void set_spherical_lap(const std::vector<double> &v, std::vector<double> &ddv) const;
   void set_sphereX(      const std::vector<double> &v, std::vector<double> &vX) const;
   void set_partial_R(    const std::vector<double> &v, std::vector<double> &dv) const;
   void set_partial_r(    const std::vector<double> &v, std::vector<double> &dv) const;
   void set_partial2_r(   const std::vector<double> &v, std::vector<double> &ddv) const;
   void set_angular_power_spectrum(const std::vector<double> &v, std::vector<double> &p) const;
#if USE_CHEB
   void set_n_l_coef(const std::vector<double> &v, std::vector<double> &p) const;
#endif
   /* \partial_t f - p 
    * */
   double norm_indep_res(
         const double dt, 
         const std::vector<double> &f_n, 
         const std::vector<double> &f_np1, 
         const std::vector<double> &p) const; 

   double total_variation(const std::vector<double> &f) const; 
   /* Low pass filter in spectral space 
    */
   void filter(std::vector<double> &v) const;
   /*=========================================================================*/
   inline void get_row_R(const size_t it, const size_t ip, 
         const std::vector<double> &in,
         std::vector<double> &out) const 
   {
      assert(in.size() ==_nx*_nphi*_nlat);
      assert(out.size()==_nx);
      for (size_t ix=0; ix<_nx; ix++) {
         out[ix] = in[indx_R_th_ph(ix,it,ip)];
      }
   } 
   inline void set_row_R(const size_t it, const size_t ip, 
         const std::vector<double> &in,
         std::vector<double> &out) const 
   {
      assert(in.size() ==_nx);
      assert(out.size()==_nx*_nphi*_nlat);
      for (size_t ix=0; ix<_nx; ix++) {
         out[indx_R_th_ph(ix,it,ip)] = in[ix];
      }
   } 
   /*=========================================================================*/
   inline void get_row_th(const size_t ix, const size_t ip, 
         const std::vector<double> &in,
         std::vector<double> &out) const 
   {
      assert(in.size() ==_nx*_nphi*_nlat);
      assert(out.size()==_nlat);
      for (size_t it=0; it<_nlat; it++) {
         out[it] = in[indx_R_th_ph(ix,it,ip)];
      }
   } 
   inline void set_row_th(const size_t ix, const size_t ip, 
         const std::vector<double> &in,
         std::vector<double> &out) const 
   {
      assert(in.size() ==_nlat);
      assert(out.size()==_nx*_nphi*_nlat);
      for (size_t it=0; it<_nlat; it++) {
         out[indx_R_th_ph(ix,it,ip)] = in[it];
      }
   } 
   /*=========================================================================*/
   inline void get_row_ph(const size_t ix, const size_t it, 
         const std::vector<double> &in,
         std::vector<double> &out) const 
   {
      assert(in.size() ==_nx*_nphi*_nlat);
      assert(out.size()==_nphi);
      for (size_t ip=0; ip<_nphi; ip++) {
         out[ip] = in[indx_R_th_ph(ix,it,ip)];
      }
   } 
   inline void set_row_ph(const size_t ix, const size_t it, 
         const std::vector<double> &in,
         std::vector<double> &out) const 
   {
      assert(in.size() ==_nphi);
      assert(out.size()==_nx*_nphi*_nlat);
      for (size_t ip=0; ip<_nlat; ip++) {
         out[indx_R_th_ph(ix,it,ip)] = in[ip];
      }
   } 
   /*=========================================================================*/
   inline void get_row_n(const size_t il, 
         const std::vector<double> &in, 
         std::vector<double> &out) const
   {
      assert(in.size() ==_nx*_nl);
      assert(out.size()==_nx);
      for (size_t ix=0; ix<_nx; ix++) {
         out[ix] = in[indx_r_l(ix,il)];
      }
   } 
   inline void set_row_n(const size_t il, 
         const std::vector<double> &in, 
         std::vector<double> &out) const
   {
      assert(in.size() ==_nx);
      assert(out.size()==_nx*_nl);
      for (size_t ix=0; ix<_nx; ix++) {
         out[indx_r_l(ix,il)] = in[ix];
      }
   } 
   /*=========================================================================*/
   inline void get_row_R_th(const size_t ip, 
         const std::vector<double> &in,
         std::vector<double> &out) const
   {
      assert(in.size() ==_nx*_nphi*_nlat);
      assert(out.size()==_nx*_nlat);
      for (size_t ix=0; ix<_nx;   ix++) {
      for (size_t it=0; it<_nlat; it++) {
         out[indx_r_th(ix,it)] = in[indx_R_th_ph(ix,it,ip)];
      }
      }
   }
   inline void set_row_R_th(const size_t ip, 
         const std::vector<double> &in,
         std::vector<double> &out) const
   {
      assert(in.size() ==_nx*_nlat);
      assert(out.size()==_nx*_nphi*_nlat);
      for (size_t ix=0; ix<_nx;   ix++) {
      for (size_t it=0; it<_nlat; it++) {
         out[indx_R_th_ph(ix,it,ip)] = in[indx_r_th(ix,it)];
      }
      }
   }
   /*=========================================================================*/
   inline void get_row_R_ph(const size_t it, 
         const std::vector<double> &in,
         std::vector<double> &out) const
   {
      assert(in.size() ==_nx*_nphi*_nlat);
      assert(out.size()==_nx*_nlat);
      for (size_t ix=0; ix<_nx;   ix++) {
      for (size_t ip=0; ip<_nphi; ip++) {
         out[indx_r_ph(ix,ip)] = in[indx_R_th_ph(ix,it,ip)];
      }
      }
   }
   inline void set_row_R_ph(const size_t it, 
         const std::vector<double> &in,
         std::vector<double> &out) const
   {
      assert(in.size() ==_nx*_nphi);
      assert(out.size()==_nx*_nphi*_nlat);
      for (size_t ix=0; ix<_nx;   ix++) {
      for (size_t ip=0; ip<_nphi; ip++) {
         out[indx_R_th_ph(ix,it,ip)] = in[indx_r_ph(ix,ip)];
      }
      }
   }
   /*=========================================================================*/
   inline void get_row_th_ph(const size_t ix, 
         const std::vector<double> &in,
         std::vector<double> &out) const
   {
      assert(in.size() ==_nx*_nphi*_nlat);
      assert(out.size()==_nphi*_nlat);
      for (size_t it=0; it<_nlat; it++) {
      for (size_t ip=0; ip<_nphi; ip++) {
         out[indx_th_ph(it,ip)] = in[indx_R_th_ph(ix,it,ip)];
      }
      }
   }
   inline void set_row_th_ph(const size_t ix, 
         const std::vector<double> &in,
         std::vector<double> &out) const
   {
      assert(in.size() ==_nphi*_nlat);
      assert(out.size()==_nx*_nphi*_nlat);

      for (size_t it=0; it<_nlat; it++) {
      for (size_t ip=0; ip<_nphi; ip++) {
         out[indx_R_th_ph(ix,it,ip)] = in[indx_th_ph(it,ip)];
      }
      }
   }
   /*======================================================================*/
   inline void set_row_l(const size_t ix, 
         const std::vector<double> &lvals, 
         std::vector<double> &Rlvals) const
   {
      assert(lvals.size() ==_nl);
      assert(Rlvals.size()==_nx*_nl);
      for (size_t il=0; il<_nl; il++) {
         Rlvals[indx_r_l(ix,il)] = lvals[il];
      }
   } 
/*=========================================================================*/
private:
   const double _cl;
   const double _Rmin;
   const double _Rmax;
   const size_t _nx;
   const size_t _nl;
   const size_t _nm;
   const size_t _nlat;
   const size_t _nphi;

   Sphere _sphere;
#if USE_CHEB
   Cheb _radial;
#else
   FD _radial;
#endif
   std::vector<std::vector<double>> _R_th_ph;
   std::vector<std::vector<double>> _x_y_z;
   std::vector<std::vector<double>> _x_z;
   std::vector<std::vector<double>> _th_ph;
   std::vector<std::vector<double>> _R_th;
   std::vector<std::vector<double>> _R_l;
#if USE_CHEB
   std::vector<std::vector<double>> _n_l;
#endif
   std::vector<double> _dR_over_dr;
   std::vector<double> _d2R_over_dr2;

public:
   inline double cl() const {return _cl;}
   inline double Rmin() const {return _Rmin;}
   inline double Rmax() const {return _Rmax;}
   inline size_t nx() const {return _nx;}
   inline size_t nl() const {return _nl;}
   inline size_t nm() const {return _nm;}
   inline size_t nlat() const {return _nlat;}
   inline size_t nphi() const {return _nphi;}

   /* 
    * accesing indices of 1,2, and 3d arrays 
    */
   inline size_t indx_r_th(const size_t ix, const size_t it) const
   {
      return (_nlat)*(ix) + (it);
   }
   inline size_t indx_r_ph(const size_t ix, const size_t ip) const
   {
      return (_nphi)*(ix) + (ip);
   }
   inline size_t indx_r_l(const size_t ix, const size_t il) const
   {
         return (_nl)*(ix) + (il);
   } 
   inline size_t indx_R_th_ph(const size_t ix, const size_t it, const size_t ip) const 
   {
      #if SHTNS_CONTIGUOUS_LONGITUDES
      return (_nphi)*(_nlat)*(ix) + (_nlat)*(ip) + (it);
      #else
      return (_nlat)*(_nphi)*(ix) + (_nphi)*(it) + (ip);
      #endif
   }
   inline size_t indx_th_ph(const size_t it, const size_t ip) const 
   {
      #if SHTNS_CONTIGUOUS_LONGITUDES
      return (_nlat)*(ip) + (it);
      #else 
      return (_nphi)*(it) + (ip);
      #endif
   }
   inline double R_to_r(const double R) const
   {
      return (fabs(R-_cl)>1e-16) ? (R/(1.0 - (R/_cl))) : 1e12; 
   }
   inline std::vector<double> R_th_ph(const size_t ix, const size_t it, const size_t ip) const 
   {
      return _R_th_ph[indx_R_th_ph(ix, it, ip)];
   }
   inline std::vector<double> R_th_ph(const size_t i) const 
   {
      return _R_th_ph[i];
   }
   inline std::vector<double> x_y_z(const size_t ix, const size_t it, const size_t ip) const 
   {
      return _x_y_z[indx_R_th_ph(ix, it, ip)];
   }
   inline std::vector<double> x_y_z(const size_t i) const 
   {
      return _x_y_z[i];
   }
   inline std::vector<double> th_ph(const size_t it, const size_t ip) const 
   {
      return _th_ph[indx_th_ph(it, ip)];
   }
   inline std::vector<double> th_ph(const size_t i) const 
   {
      return _th_ph[i];
   }
   inline std::vector<double> R_l(const size_t ix, const size_t il) const 
   {
      return _R_l[indx_r_l(ix, il)];
   }
   inline std::vector<double> R_l(const size_t i) const 
   {
      return _R_l[i];
   }
#if USE_CHEB
   inline std::vector<double> n_l(const size_t ix, const size_t il) const 
   {
      return _n_l[indx_r_l(ix, il)];
   }
   inline std::vector<double> n_l(const size_t i) const 
   {
      return _n_l[i];
   }
#endif
};
#endif /* _GRID_HPP */
