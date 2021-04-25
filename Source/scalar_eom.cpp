#include <cassert>
#include <iostream>
#include <iomanip>

#include "scalar_eom.hpp"
/*==========================================================================*/
Scalar_eom::Scalar_eom(
      const Grid &grid,
      const Params &params
      )
:  _n{grid.nx()*grid.nlat()*grid.nphi()},
   _dt{params.dt()},
   _km1{params.km1()},
   _k0{params.k0()},
   _k1{params.k1()},
   _k2{params.k2()},
   _v2{params.V2()},
   _v3{params.V3()},
   _v4{params.V4()},
  
   _pre(_n,0),

   _p_p(        _n,0),
   _p_dr_f(     _n,0),
   _p_dr_p(     _n,0),
   _p_dr_dr_f(  _n,0),
   _p_dphi_dr_f(_n,0),
   _p_lap_f(    _n,0),

   _p_p_p(        _n,0),
   _p_p_dr_f(     _n,0),
   _p_dr_f_dr_f(  _n,0),
   _p_dr_f_dphi_f(_n,0),
   _p_sphereX_f(  _n,0),

   _rho_vv(     _n,0),
   _rho_vr(     _n,0),
   _rho_rr(     _n,0),
   _rho_rphi(   _n,0),
   _rho_sphereX(_n,0),

   _not_at_spat_infty(_n,0)
{
   std::cout<<"Initializing Scalar_eom"<<std::endl;
   const size_t nx   = grid.nx();
   const size_t nlat = grid.nlat();
   const size_t nphi = grid.nphi();

   const double cl = params.cl();

   const double m = params.bh_mass();
   const double a = params.bh_spin();

   for (size_t ix=0; ix<nx;   ix++) {
   for (size_t ip=0; ip<nphi; ip++) {
   for (size_t it=0; it<nlat; it++) {
      const size_t indx = grid.indx_r_th_ph(ix,it,ip);
      
      std::vector<double> R_th_ph = grid.R_th_ph(ix, it, ip); 
      std::vector<double> r_th_ph = grid.r_th_ph(ix, it, ip); 

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
         (fabs(m-fabs(a))<1e-16) ? (m)   : (pow(a,2)/rp)
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

      _rho_vv[indx]      = 0.5*(1.0 + 2.0*m*inv_r/Sigma);
      _rho_vr[indx]      = 2.0*m*inv_r/Sigma;
      _rho_rr[indx]      = 0.5*( 
            ((1.0+(2*m*inv_r)+pow(a*inv_r,2))/Sigma)
         -  (4.0*m*inv_r/(2.0*m*inv_r + Sigma))     
         );
      _rho_rphi[indx]    = a*pow(inv_r,2)/Sigma;
      _rho_sphereX[indx] = 0.5*pow(inv_r,2)/Sigma;
   }
   }
   }
   for (size_t i=0; i<_n; i++) {
      std::vector<double> R_th_ph = grid.R_th_ph(i); 
      if (fabs(R_th_ph[0]-cl)<1e-16) {
         _not_at_spat_infty[i] = 0.0;
      } else {
         _not_at_spat_infty[i] = 1.0;
      }
   }
   std::cout<<"Finished initializing Scalar_eom"<<std::endl;
}
/*==========================================================================*/
Scalar_eom::~Scalar_eom()
{
}
/*==========================================================================*/
void Scalar_eom::set_k(
      const Grid &grid,
      const std::vector<double> &f,
      const std::vector<double> &p,
      std::vector<double> &f_k,
      std::vector<double> &p_k
      ) const
{
   assert(f.size()  ==_n);
   assert(p.size()  ==_n);
   assert(f_k.size()==_n);
   assert(p_k.size()==_n);

   std::vector<double> _dr_f(     _n);
   std::vector<double> _lap_f(    _n);
   std::vector<double> _dr_p(     _n);
   std::vector<double> _dr_dr_f(  _n);
   std::vector<double> _dphi_f(   _n);
   std::vector<double> _dphi_dr_f(_n);
   std::vector<double> _sphereX_f(_n);

   grid.set_partial_r(    f,    _dr_f);
   grid.set_partial_r(    p,    _dr_p);
   grid.set_partial_r(_dr_f, _dr_dr_f);

   grid.set_partial_phi(    f, _dphi_f);
   grid.set_partial_phi(_dr_f, _dphi_dr_f);

   grid.set_sphereX(f, _sphereX_f);

   grid.set_spherical_lap(f, _lap_f);

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
      f_k[i] *= _not_at_spat_infty[i];
      p_k[i] *= _not_at_spat_infty[i];
   }
}
/*==========================================================================*/
void Scalar_eom::set_level(
      const int level,
      Field &f,
      Field &p) const
{
   if (level==2) {
      for (size_t i=0; i<_n; i++) {
         f.l2[i] = f.n[i] + 0.5*_dt*f.k1[i];
         p.l2[i] = p.n[i] + 0.5*_dt*p.k1[i];
      }
   } else
   if (level==3) {
      for (size_t i=0; i<_n; i++) {
         f.l3[i] = f.n[i] + 0.5*_dt*f.k2[i];
         p.l3[i] = p.n[i] + 0.5*_dt*p.k2[i];
      }
   } else
   if (level==4) {
      for (size_t i=0; i<_n; i++) {
         f.l4[i] = f.n[i] + _dt*f.k3[i];
         p.l4[i] = p.n[i] + _dt*p.k3[i];
      }
   } else
   if (level==5) {
      for (size_t i=0; i<_n; i++) {
         f.np1[i] = f.n[i] + (_dt/6.0)*(f.k1[i] + 2.0*f.k2[i] + 2.0*f.k3[i] + f.k4[i]);
         p.np1[i] = p.n[i] + (_dt/6.0)*(p.k1[i] + 2.0*p.k2[i] + 2.0*p.k3[i] + p.k4[i]);
      }
   } else {
      std::cout<<"ERROR(set_level): level = "<<level<<std::endl;
   }
}
/*==========================================================================*/
/* Fourth order Runge-Kutta integrator */
/*==========================================================================*/
void Scalar_eom::time_step(const Grid &grid, Field &f, Field &p) const
{
   assert(f.size==_n);
   assert(p.size==_n);

   grid.filter(f.n);
   grid.filter(p.n);

   set_k(grid, f.n, p.n, f.k1, p.k1);
   set_level(2, f, p);

   set_k(grid, f.l2, p.l2, f.k2, p.k2);
   set_level(3, f, p);

   set_k(grid, f.l3, p.l3, f.k3, p.k3);
   set_level(4, f, p);

   set_k(grid, f.l4, p.l4, f.k4, p.k4);
   set_level(5, f, p);
}
/*==========================================================================*/
void Scalar_eom::set_rho(
      const Grid &grid, 
      const std::vector<double> &f,
      const std::vector<double> &p,
      std::vector<double> &rho
      ) const
{
   assert(f.size()  ==_n);
   assert(p.size()  ==_n);
   assert(rho.size()==_n);

   std::vector<double> _dr_f(_n);
   std::vector<double> _dphi_f(_n);
   std::vector<double> _sphereX_f(_n);

   grid.set_partial_r(  f,      _dr_f);
   grid.set_partial_phi(f,    _dphi_f);
   grid.set_sphereX(    f, _sphereX_f);

   for (size_t i=0; i<_n; i++) {
      const double k = ( 
            _k0 
         +  _km1/(f[i]+1)
         +  _k1*f[i]
         +  (_k2/2.0)*pow(f[i],2)
         );
      const double v = (
            (_v2/2.0) *pow(f[i],2)
         +  (_v3/6.0 )*pow(f[i],3)
         +  (_v4/24.0)*pow(f[i],4)
         );
         
      rho[i] = 
      +  k*(
            _rho_vv[i]*pow(p[i],2) 
         +  _rho_vr[i]*p[i]*_dr_f[i]
         +  _rho_rr[i]*pow(_dr_f[i],2)
         +  _rho_rphi[i]*_dr_f[i]*_dphi_f[i]
         +  _rho_sphereX[i]*_sphereX_f[i]
      )
      +  v
      ;
   }
}
