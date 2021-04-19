#include <cassert>
#include <iostream>
#include <iomanip>

#include "scalar_eom.hpp"
#include "params.hpp"
#include "grid.hpp"
/*==========================================================================*/
namespace Eom {
/*==========================================================================*/
namespace {
   size_t _n;
   double _dt;

   double _km1;
   double _k0;
   double _k1;
   double _k2;

   double _v2;
   double _v3;
   double _v4;
   std::vector<double> _pre;
   std::vector<double> _p_p;
   std::vector<double> _p_dr_f;
   std::vector<double> _p_dr_p;
   std::vector<double> _p_dr_dr_f;
   std::vector<double> _p_dphi_dr_f;
   std::vector<double> _p_lap_f;

   std::vector<double> _p_p_p;
   std::vector<double> _p_p_dr_f;
   std::vector<double> _p_dr_f_dr_f;
   std::vector<double> _p_dr_f_dphi_f;
   std::vector<double> _p_sphereX_f;

   std::vector<double> _dr_f;
   std::vector<double> _lap_f;
   std::vector<double> _dr_p;
   std::vector<double> _dr_dr_f;
   std::vector<double> _dphi_f;
   std::vector<double> _dphi_dr_f;

   std::vector<double> _sphereX_f;
}
/*==========================================================================*/
void init()
{
   _n  = Params::nx_nlat_nphi();
   _dt = Params::dt();

   _pre.resize(_n,0);

   _p_p.resize(        _n,0);
   _p_dr_f.resize(     _n,0);
   _p_dr_p.resize(     _n,0);
   _p_dr_dr_f.resize(  _n,0);
   _p_dphi_dr_f.resize(_n,0);
   _p_lap_f.resize(    _n,0);

   _p_p_p.resize(        _n,0);
   _p_p_dr_f.resize(     _n,0);
   _p_dr_f_dr_f.resize(  _n,0);
   _p_dr_f_dphi_f.resize(_n,0);
   _p_sphereX_f.resize(  _n,0);

   _dr_f.resize(     _n);
   _lap_f.resize(    _n);
   _dr_p.resize(     _n);
   _dr_dr_f.resize(  _n);
   _dphi_f.resize(   _n);
   _dphi_dr_f.resize(_n);

   const size_t nx   = Params::nx();
   const size_t nlat = Params::nlat();
   const size_t nphi = Params::nphi();

   const double cl  = Params::cl();

   const double m = Params::bh_mass();
   const double a = Params::bh_spin();

   _km1 = Params::km1();
   _k0  = Params::k0();
   _k1  = Params::k1();
   _k2  = Params::k2();

   _v2  = Params::V2();
   _v3  = Params::V3();
   _v4  = Params::V4();

   for (size_t ix=0; ix<nx;   ix++) {
   for (size_t ip=0; ip<nphi; ip++) {
   for (size_t it=0; it<nlat; it++) {
      const size_t indx = Grid::indx(ix,it,ip);
      
      std::vector<double> R_th_ph = Grid::R_th_ph(ix, it, ip); 
      std::vector<double> r_th_ph = Grid::r_th_ph(ix, it, ip); 

      const double inv_r = (fabs(R_th_ph[0]-cl)<1e-16) ? 0 : (1.0/r_th_ph[0]);
      const double th = r_th_ph[1];

      /* Sigma and Delta:
       * Divide by r^2 to reduce infty/infy type errors 
       * in computing coefficients. */
      const double Sigma = 1.0 + pow(inv_r,2)*pow(a,2)*pow(cos(th),2); 

      const double rp = 
         (fabs(m-fabs(a))<1e-16) ? (m)   : (
         (fabs(a)        <1e-16) ? (2*m) : (m + pow((m-a)*(m+a),0.5))
         );
      const double rm = 
         (fabs(m-a)<1e-16) ? (m)   : (pow(a,2)/rp)
         ;
      const double Delta = (1.0 - rp*inv_r)*(1.0 - rm*inv_r); 

      _p_p[indx]         = (2.0*m/Sigma)*pow(inv_r,2);
      _p_dr_f[indx]      = 2.0*(inv_r - (m*pow(inv_r,2)))/Sigma;
      _p_dr_p[indx]      = 4.0*m*inv_r/Sigma;
      _p_dr_dr_f[indx]   = (Delta/Sigma);
      _p_dphi_dr_f[indx] = (2.0*a/Sigma)*pow(inv_r,2);
      _p_lap_f[indx]     = (1.0/Sigma)*pow(inv_r,2);

      _p_p_p[indx]         = -(1.0 + (2.0*m*inv_r/Sigma));
      _p_p_dr_f[indx]      = 4.0*m*inv_r/Sigma;
      _p_dr_f_dr_f[indx]   = Delta/Sigma;
      _p_dr_f_dphi_f[indx] = 2.0*a*pow(inv_r,2)/Sigma;
      _p_sphereX_f[indx]   = pow(inv_r,2)/Sigma;

      _pre[indx] = 1.0 + (2.0*m*inv_r/Sigma);

      _p_p[indx]         /= _pre[indx];
      _p_dr_f[indx]      /= _pre[indx];
      _p_dr_p[indx]      /= _pre[indx];
      _p_dr_dr_f[indx]   /= _pre[indx];
      _p_dphi_dr_f[indx] /= _pre[indx];
      _p_lap_f[indx]     /= _pre[indx];

      _p_p_p[indx]         /= _pre[indx];
      _p_p_dr_f[indx]      /= _pre[indx];
      _p_dr_f_dr_f[indx]   /= _pre[indx];
      _p_dr_f_dphi_f[indx] /= _pre[indx]; 
      _p_sphereX_f[indx]   /= _pre[indx];
   }
   }
   }
}
/*==========================================================================*/
void set_k(
      const std::vector<double> &f,
      const std::vector<double> &p,
      std::vector<double> &f_k,
      std::vector<double> &p_k
      )
{
   Grid::set_partial_r(    f,    _dr_f);
   Grid::set_partial_r(    p,    _dr_p);
   Grid::set_partial_r(_dr_f, _dr_dr_f);

   Grid::set_partial_phi(    f, _dphi_f);
   Grid::set_partial_phi(_dr_f, _dphi_dr_f);

   Grid::set_sphereX(f, _sphereX_f);

   Grid::set_spherical_lap(f, _lap_f);

   for (size_t i=0; i<_n; i++) {
      const double inverse_k = 1.0/( 
            _k0 
         +  _km1/(f[i]+1)
         +  _k1*f[i]
         +  (_k2/2.0)*pow(f[i],2)
         );
      const double kprime = 
      -  _km1/pow(f[i]+1,2)
      +  _k1
      +  _k2*f[i]
      ;
      const double vprime = inverse_k*(
            _v2*f[i]
         +  (_v3/2.0)*pow(f[i],2)
         +  (_v4/6.0)*pow(f[i],3)
         );
         
      f_k[i] = 
         p[i]
      ;
      p_k[i] = 
      +  _p_p[i]*p[i] 

      +  _p_dr_f[i]*_dr_f[i]
      +  _p_dr_p[i]*_dr_p[i] 

      +  _p_dr_dr_f[i]*_dr_dr_f[i]

      +  _p_dphi_dr_f[i]*_dphi_dr_f[i]

      +  _p_lap_f[i]*_lap_f[i]

      +  0.5*(kprime*inverse_k)*(
            _p_p_p[i]*pow(p[i],2)
         +  _p_p_dr_f[i]*p[i]*_dr_f[i]
         +  _p_dr_f_dr_f[i]*pow(_dr_f[i],2)
         +  _p_dr_f_dphi_f[i]*_dr_f[i]*_dphi_f[i]
         +  _p_sphereX_f[i]*_sphereX_f[i]
         )

      -  (1.0/_pre[i])*vprime
      ;
   }
}
/*==========================================================================*/
/* Fourth order Runge-Kutta integrator */
/*==========================================================================*/
void time_step(Field &f, Field &p)
{
   assert(f.size==_n);
   assert(p.size==_n);

   Grid::filter(f.n);
   Grid::filter(p.n);
   /*--------------------------------------------*/
   set_k(f.n, p.n, 
         f.k1, p.k1
      );
   /*--------------------------------------------*/
   for (size_t i=0; i<_n; i++) {
      f.l2[i] = f.n[i] + 0.5*_dt*f.k1[i];
      p.l2[i] = p.n[i] + 0.5*_dt*p.k1[i];
   }
   set_k(f.l2, p.l2, 
         f.k2, p.k2
      );
   /*--------------------------------------------*/
   for (size_t i=0; i<_n; i++) {
      f.l3[i] = f.n[i] + 0.5*_dt*f.k2[i];
      p.l3[i] = p.n[i] + 0.5*_dt*p.k2[i];
   }
   set_k(f.l3, p.l3, 
         f.k3, p.k3
      );
   /*--------------------------------------------*/
   for (size_t i=0; i<_n; i++) {
      f.l4[i] = f.n[i] + _dt*f.k3[i];
      p.l4[i] = p.n[i] + _dt*p.k3[i];
   }
   set_k(f.l4, p.l4,
         f.k4, p.k4
      );
   /*--------------------------------------------*/
   for (size_t i=0; i<_n; i++) {
      f.np1[i] = f.n[i] + (_dt/6.0)*(f.k1[i] + 2.0*f.k2[i] + 2.0*f.k3[i] + f.k4[i]);
      p.np1[i] = p.n[i] + (_dt/6.0)*(p.k1[i] + 2.0*p.k2[i] + 2.0*p.k3[i] + p.k4[i]);
   }
}
/*==========================================================================*/
} /* Eom */
