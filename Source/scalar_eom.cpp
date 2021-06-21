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

   _dr_f(     _n,0),
   _lap_f(    _n,0),
   _dr_p(     _n,0),
   _dr_dr_f(  _n,0),
   _dphi_f(   _n,0),
   _dphi_dr_f(_n,0),
   _sphereX_f(_n,0)
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
      const size_t indx = grid.indx_R_th_ph(ix,it,ip);
      
      std::vector<double> R_th_ph = grid.R_th_ph(ix, it, ip); 
  
      const double r = grid.R_to_r(R_th_ph[0]);

      const double inv_r = (fabs(R_th_ph[0]-cl)<1e-16) ? 0 : (1.0/r);
      const double th = R_th_ph[1];

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

      _p_p[indx]         = 2.*m*pow(inv_r,2)/Sigma;
      _p_dr_f[indx]      = 2.*(inv_r - (m*pow(inv_r,2)))/Sigma;
      _p_dr_p[indx]      = 4.*m*inv_r/Sigma;
      _p_dr_dr_f[indx]   = Delta/Sigma;
      _p_dphi_dr_f[indx] = 2.*a*pow(inv_r,2)/Sigma;
      _p_lap_f[indx]     = pow(inv_r,2)/Sigma;

      _p_p_p[indx]         = -(1. + (2.*m*inv_r/Sigma));
      _p_p_dr_f[indx]      = 4.*m*inv_r/Sigma;
      _p_dr_f_dr_f[indx]   = Delta/Sigma;
      _p_dr_f_dphi_f[indx] = 2.*a*pow(inv_r,2)/Sigma;
      _p_sphereX_f[indx]   = pow(inv_r,2)/Sigma;

      _pre[indx] = 1./(1. + (2.*m*inv_r/Sigma));

      _p_p[indx]         *= _pre[indx];
      _p_dr_f[indx]      *= _pre[indx];
      _p_dr_p[indx]      *= _pre[indx];
      _p_dr_dr_f[indx]   *= _pre[indx];
      _p_dphi_dr_f[indx] *= _pre[indx];
      _p_lap_f[indx]     *= _pre[indx];

      _p_p_p[indx]         *= _pre[indx];
      _p_p_dr_f[indx]      *= _pre[indx];
      _p_dr_f_dr_f[indx]   *= _pre[indx];
      _p_dr_f_dphi_f[indx] *= _pre[indx]; 
      _p_sphereX_f[indx]   *= _pre[indx];

      _rho_vv[indx]      = 0.5*(1. + 2.*m*inv_r/Sigma);
      _rho_vr[indx]      = 2. *m*inv_r/Sigma;
      _rho_rr[indx]      = 0.5*( 
            ((1. + (2.*m*inv_r)+pow(a*inv_r,2))/Sigma)
         -  (4.*m*inv_r/(2.*m*inv_r + Sigma))     
         );
      _rho_rphi[indx]    =   a*pow(inv_r,2)/Sigma;
      _rho_sphereX[indx] = 0.5*pow(inv_r,2)/Sigma;
   }
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
      )
{
   assert(f.size()==_n);
   assert(p.size()==_n);
   assert(f_k.size()==_n);
   assert(p_k.size()==_n);

   #pragma omp parallel sections
   {
      #pragma omp section
      {
         grid.set_partial_r(f, _dr_f);
      }
      #pragma omp section
      {
         grid.set_partial_r(p, _dr_p);
      }
      #pragma omp section
      {
         grid.set_partial2_r(f, _dr_dr_f);
      }
      #pragma omp section
      {
         grid.set_partial_phi(    f, _dphi_f);
         grid.set_partial_r(_dphi_f, _dphi_dr_f);
      }
      #pragma omp section
      {
         grid.set_sphereX(f, _sphereX_f);
      }
      #pragma omp section
      {
         grid.set_spherical_lap(f, _lap_f);
      }
   }
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

      -  _pre[i]*vprime
      ;
   }
}
/*==========================================================================*/
void Scalar_eom::set_level(
      const int level,
      Field &f,
      Field &p) const
{
   if (level==1) {
#pragma omp parallel for
      for (size_t i=0; i<_n; i++) {
         f.l[i] = f.n[i];
         p.l[i] = p.n[i];
      }
   } else
   if (level==2) {
#pragma omp parallel for
      for (size_t i=0; i<_n; i++) {
         f.l[i] = f.n[i] + 0.5*_dt*f.k[i];
         p.l[i] = p.n[i] + 0.5*_dt*p.k[i];

         f.np1[i] = f.n[i] + (_dt/6.0)*(f.k[i]);
         p.np1[i] = p.n[i] + (_dt/6.0)*(p.k[i]);
      }
   } else
   if (level==3) {
#pragma omp parallel for
      for (size_t i=0; i<_n; i++) {
         f.l[i] = f.n[i] + 0.5*_dt*f.k[i];
         p.l[i] = p.n[i] + 0.5*_dt*p.k[i];

         f.np1[i] += (_dt/3.0)*(f.k[i]);
         p.np1[i] += (_dt/3.0)*(p.k[i]);
      }
   } else
   if (level==4) {
#pragma omp parallel for
      for (size_t i=0; i<_n; i++) {
         f.l[i] = f.n[i] + _dt*f.k[i];
         p.l[i] = p.n[i] + _dt*p.k[i];

         f.np1[i] += (_dt/3.0)*(f.k[i]);
         p.np1[i] += (_dt/3.0)*(p.k[i]);
      }
   } else
   if (level==5) {
#pragma omp parallel for
      for (size_t i=0; i<_n; i++) {
         f.np1[i] += (_dt/6.0)*(f.k[i]);
         p.np1[i] += (_dt/6.0)*(p.k[i]);
      }
   } else {
      std::cout<<"ERROR(set_level): level = "<<level<<std::endl;
   }
}
/*==========================================================================*/
/* Fourth order Runge-Kutta integrator; 
 * we reuse the levels and time derivatives (k) to save memory */
/*==========================================================================*/
void Scalar_eom::time_step(const Grid &grid, Field &f, Field &p)
{
   assert(f.size==_n);
   assert(p.size==_n);

   #pragma omp parallel sections
   {
      #pragma omp section
      {
         grid.filter(f.n);
      }
      #pragma omp section
      {
         grid.filter(p.n);
      }
   }
   set_level(1, f, p);

   set_k(grid, f.l, p.l, f.k, p.k);
   set_level(2, f, p);

   set_k(grid, f.l, p.l, f.k, p.k);
   set_level(3, f, p);

   set_k(grid, f.l, p.l, f.k, p.k);
   set_level(4, f, p);

   set_k(grid, f.l, p.l, f.k, p.k);
   set_level(5, f, p);
}
/*==========================================================================*/
void Scalar_eom::set_rho(
      const Grid &grid, 
      const std::vector<double> &f,
      const std::vector<double> &p,
      std::vector<double> &rho
      )
{
   assert(f.size()  ==_n);
   assert(p.size()  ==_n);
   assert(rho.size()==_n);

   #pragma omp parallel sections
   {
      #pragma omp section
      {
         grid.set_partial_r(f, _dr_f);
      }
      #pragma omp section
      {
         grid.set_partial_phi(f, _dphi_f);
      }
      #pragma omp section
      {
         grid.set_sphereX(f, _sphereX_f);
      }
   }
   for (size_t i=0; i<_n; i++) {
      const double k = ( 
            _k0 
         +  _km1/(f[i]+1)
         +  _k1*f[i]
         +  (_k2/2.)*pow(f[i],2)
         );
      const double v = (
            (_v2/2.) *pow(f[i],2)
         +  (_v3/6. )*pow(f[i],3)
         +  (_v4/24.)*pow(f[i],4)
         );
         
      rho[i] = 
      +  k*(1./2.)*(
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
