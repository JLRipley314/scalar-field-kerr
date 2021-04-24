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

   size_t indx(const size_t i_x, const size_t i_th, const size_t i_ph) const; 
   size_t sphere_indx(const size_t i_th, const size_t i_ph) const; 

   std::vector<double> r_th_ph(const size_t i_x, const size_t i_th, const size_t i_ph) const; 
   std::vector<double> R_th_ph(const size_t i_x, const size_t i_th, const size_t i_ph) const; 
   std::vector<double> x_y_z(  const size_t i_x, const size_t i_th, const size_t i_ph) const; 
   std::vector<double> r_th_ph(const size_t i) const; 
   std::vector<double> R_th_ph(const size_t i) const; 
   std::vector<double> x_y_z(  const size_t i) const; 

   std::vector<double> th_ph(const size_t i_th, const size_t i_ph) const; 
   std::vector<double> R_l(  const size_t i_x,  const size_t i_l) const; 
   std::vector<double> n_l(  const size_t i_x, const size_t i_l) const; 
   std::vector<double> th_ph(const size_t i) const; 
   std::vector<double> R_l(  const size_t i) const; 
   std::vector<double> n_l(  const size_t i) const; 

   std::vector<double> compute_ylm(const size_t l_ang, const size_t m_ang) const;
   /* Functions acting on full 3d grid functions 
    */
   void set_partial_phi(  const std::vector<double> &v, std::vector<double> &dv) const;
   void set_spherical_lap(const std::vector<double> &v, std::vector<double> &ddv) const;
   void set_sphereX(      const std::vector<double> &v, std::vector<double> &vX) const;
   void set_partial_r(    const std::vector<double> &v, std::vector<double> &dv) const;
   void set_angular_power_spectrum(const std::vector<double> &v, std::vector<double> &p) const;
   void set_n_l_coef(              const std::vector<double> &v, std::vector<double> &p) const;
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

   void get_row_R(const size_t j, const size_t k, 
         const std::vector<double> &in,
         std::vector<double> &out
      ) const; 
   void set_row_R(const size_t j, const size_t k, 
         const std::vector<double> &in,
         std::vector<double> &out
      ) const; 

   void get_row_th(const size_t i, const size_t k, 
         const std::vector<double> &in,
         std::vector<double> &out
      ) const; 
   void set_row_th(const size_t i, const size_t k, 
         const std::vector<double> &in,
         std::vector<double> &out
      ) const; 

   void get_row_ph(const size_t i, const size_t j, 
         const std::vector<double> &in,
         std::vector<double> &out
      ) const; 
   void set_row_ph(const size_t i, const size_t j, 
         const std::vector<double> &in,
         std::vector<double> &out
      ) const; 

   void get_row_n(const size_t il, 
         const std::vector<double> &in, 
         std::vector<double> &out
      ) const;
   void set_row_n(const size_t il, 
         const std::vector<double> &in, 
         std::vector<double> &out
      ) const;

   void get_row_R_th(const size_t k, 
         const std::vector<double> &in,
         std::vector<double> &out
      ) const;
   void set_row_R_th(const size_t k, 
         const std::vector<double> &in,
         std::vector<double> &out
      ) const;

   void get_row_R_ph(const size_t j, 
         const std::vector<double> &in,
         std::vector<double> &out
      ) const;
   void set_row_R_ph(const size_t j, 
         const std::vector<double> &in,
         std::vector<double> &out
      ) const;

   void get_row_th_ph(const size_t i, 
         const std::vector<double> &in,
         std::vector<double> &out
      ) const;
   void set_row_th_ph(const size_t i, 
         const std::vector<double> &in,
         std::vector<double> &out
      ) const;

   void set_row_l(const size_t ix, 
         const std::vector<double> &lvals, 
         std::vector<double> &Rlvals
      ) const;


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
   Cheb _cheb;

   std::vector<std::vector<double>> _r_th_ph;
   std::vector<std::vector<double>> _R_th_ph;
   std::vector<std::vector<double>> _x_y_z;
   std::vector<std::vector<double>> _x_z;
   std::vector<std::vector<double>> _th_ph;
   std::vector<std::vector<double>> _R_th;
   std::vector<std::vector<double>> _R_l;
   std::vector<std::vector<double>> _n_l;

   std::vector<double> _partial_R_to_partial_r;

public:
   double cl() const {return _cl;}
   double Rmin() const {return _Rmin;}
   double Rmax() const {return _Rmax;}
   size_t nx() const {return _nx;}
   size_t nl() const {return _nl;}
   size_t nm() const {return _nm;}
   size_t nlat() const {return _nlat;}
   size_t nphi() const {return _nphi;}
};
#endif /* _GRID_HPP */
