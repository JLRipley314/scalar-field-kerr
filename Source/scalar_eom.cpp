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

   cs_M_R(grid.nlat()*grid.nphi(),0),
   cs_P_L(grid.nlat()*grid.nphi(),0)
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

      const size_t indx_sphere = grid.indx_th_ph(it,ip);

      /* Right boundary */
      if (ix==0) {
         const double inter = 2.0*m*inv_r/Sigma;

         cs_M_R[indx_sphere] = 
            pow(1.0 + inter, -1)*(
            -  inter
            -  pow(
                  pow(inter,2)
               +  (1.0 + inter)*(Delta/Sigma)
            , 0.5)
            );
      }
      /* Left boundary */
      if (ix==nx-1) {
         const double inter = 2.0*m*inv_r/Sigma;

         cs_P_L[indx_sphere] = 
            pow(1.0 + inter, -1)*(
            -  inter
            +  pow(
                  pow(inter,2)
               +  (1.0 + inter)*(Delta/Sigma)
            , 0.5)
            );
      }
   }
   }
   }
#if USE_HYPERBOLOIDAL 
   for (size_t ix=0; ix<nx;   ix++) {
   for (size_t ip=0; ip<nphi; ip++) {
   for (size_t it=0; it<nlat; it++) {
      const size_t indx = grid.indx_R_th_ph(ix,it,ip);
      
      std::vector<double> R_th_ph = grid.R_th_ph(ix, it, ip); 

      const double R  = R_th_ph[0];
      const double Th = R_th_ph[1];
      const double Om[i] = 1. - (R/cl);
      _inv_r[i] = Om/R;

      const double Sigma = pow(R,2) + pow(a*cos(Th)*Om,2); 

      _p_p[indx] = (
            (2*pow(a,2)*Om*(
               -  2*m*Om*R 
               + cl*(2*m*Om*(2 + Om) + R)
               ) 
             - 4*m*R*(
                  R*(-2*m*Om + R) 
                + cl*(2*m*Om*(1 + Om) + R - Om*R)
                )
             )/cl
         )/Sigma 
         ;
      _p_dr_f[indx] = (
            (2*Om*R*(
                  pow(a,2)*pow(Om,2)*(-cl + R) 
               +  R*(cl*m*Om + R*(-2*m*Om + R))
               ))/cl
         )/Sigma 
         ;
      _p_dr_p[indx] = ( 
            R*(
                 (4*m*Om - R)*R*(2*m*Om + R) 
               - pow(a,2)*pow(Om,2)*(4*m*Om + R)
            )
         )/Sigma
         ;
      _p_dr_dr_f[indx] = ( 
            pow(Om,2)*pow(R,2)*(pow(a,2)*pow(Om,2) + R*(-2*m*Om + R))  
         )/Sigma
         ;
      _p_dphi_dr_f[indx] = (
            2*a*pow(Om,2)*pow(R,2)
         )/Sigma
         ;
      _p_dphi_p[indx] = ( 
            -2*a*R*(4*m*Om + R)
         )/Sigma
         ;
      _p_dphi_f[indx] = (
            -2*a*Om*R
         )/Sigma  
         ;
      _p_lap_f[indx] = 
            pow(R,2)/Sigma 
         ;
      _p_f[indx] = ( 
            2*Om*(pow(a,2)*Om - m*R)
         )/Sigma 
         ;
      /*---------------------------------*/
      _p_p_p[indx] = (
            (pow(Om,2)*(
               pow(a,2) - 32*pow(m,2) 
            +  (16*m*Om*(-4*pow(m,2)*R + pow(a,2)*(2*m*Om + R)))/pow(R,2) 
            -  pow(a,2)*cos(2*th)
            ))/2.
         )/Sigma         
         ;
      _p_p_dr_f[indx] = (
            (-2*pow(Om,2)*(
            -  8*pow(m,2)*pow(Om,2)*R 
            +  pow(R,3) 
            +  pow(a,2)*pow(Om,2)*(4*m*Om + R))
            )/R 
         )/Sigma 
         ;
      _p_p_dphi_f[indx] = (
            (-2*a*pow(Om,2)*(4*m*Om + R))/R 
         )/Sigma
         ;
      _p_dr_f_dr_f[indx] = (
            pow(a,2)*pow(Om,6) + pow(Om,4)*R*(-2*m*Om + R)
         )/Sigma 
         ;
      _p_dr_f_dphi_f[indx] = (
            2*a*pow(Om,4)
         )/Sigma 
         ;
      _p_sphereX_f[indx] =
            pow(Om,2)/Sigma 
         ;
      _p_dr_f_f[indx] = (
            (
            -  2*pow(Om,3)*(pow(a,2)*pow(Om,2) 
            +  R*(-2*m*Om + R))
            )/R
         )/Sigma
         ;
      _p_dphi_f_f[indx] = (
            (-2*a*pow(Om,3))/R
         )/Sigma
         ;
      _p_f_f[indx] = (
            (
               pow(a,2)*pow(Om,4) 
            +  pow(Om,2)*R*(-2*m*Om + R)
            )/pow(R,2)
         )/Sigma
         ;
      _p_p_f[indx] = (
            (2*Om*(
            -  8*pow(m,2)*pow(Om,2)*R 
            +  pow(R,3) 
            +  pow(a,2)*pow(Om,2)*(4*m*Om + R)
            ))/pow(R,2)
         )/Sigma
         ;
      /*---------------------------------*/
      _pre[indx] = (-1.)*( 
            (-32*pow(m,2)*R*(2*m*Om + R) 
             + pow(a,2)*(32*pow(m,2)*pow(Om,2) + 16*m*Om*R + pow(R,2)) 
             - pow(a,2)*pow(R,2)*cos(2*th)
             )/2.
         )/Sigma
         ;
      /*---------------------------------*/
      _p_p[indx]         *= _pre[indx];
      _p_dr_f[indx]      *= _pre[indx];
      _p_dr_p[indx]      *= _pre[indx];
      _p_dphi_p[indx]    *= _pre[indx];
      _p_dphi_f[indx]    *= _pre[indx];
      _p_dr_dr_f[indx]   *= _pre[indx];
      _p_dphi_dr_f[indx] *= _pre[indx];
      _p_lap_f[indx]     *= _pre[indx];
      _p_f[indx]         *= _pre[indx];

      _p_p_p[indx]         *= _pre[indx]; 
      _p_p_dr_f[indx]      *= _pre[indx]; 
      _p_p_dphi_f[indx]    *= _pre[indx]; 
      _p_dr_f_dr_f[indx]   *= _pre[indx]; 
      _p_dr_f_dphi_f[indx] *= _pre[indx]; 
      _p_sphereX_f[indx]   *= _pre[indx];
      _p_dr_f_f[indx]      *= _pre[indx];
      _p_dphi_f_f[indx]    *= _pre[indx];
      _p_f_f[indx]         *= _pre[indx]; 
      _p_p_f[indx]         *= _pre[indx];
      /*---------------------------------*/
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
#endif

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
#if USE_HYPERBOLOIDAL
   std::vector<double> _dR_f(     _n);
   std::vector<double> _lap_f(    _n);
   std::vector<double> _dR_p(     _n);
   std::vector<double> _dR_dR_f(  _n);
   std::vector<double> _dphi_f(   _n);
   std::vector<double> _dphi_dR_f(_n);
   std::vector<double> _sphereX_f(_n);
   #pragma omp parallel sections
   {
      #pragma omp section
      {
         grid.set_partial_R(f, _dR_f);
      }
      #pragma omp section
      {
         grid.set_partial_R(p, _dR_p);
      }
      #pragma omp section
      {
         grid.set_partial2_R(f, _dR_dR_f);
      }
      #pragma omp section
      {
         grid.set_partial_phi(    f, _dphi_f);
         grid.set_partial_R(_dphi_f, _dphi_dR_f);
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
      if (fabs(_inv_r[i])<1e-16) {
         f_k[i] = 0;
         p_k[i] = 0;
      } else {
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
               (_v2    /_inv_r[i])   *f[i]
            +  (_v3/2.0)         *pow(f[i],2)
            +  (_v4/6.0)*_inv_r[i]*pow(f[i],3)
            );

         f_k[i] = 
            p[i]
         ;
         p_k[i] = 
            _p_p[i]*p[i] 

         +  _p_dR_f[i]*_dR_f[i]

         +  _p_dR_p[i]*_dR_p[i]

         +  _p_dR_dR_f[i]*_dR_dR_f[i]

         +  _p_dphi_dR_f[i]*_dphi_dR_f[i]

         +  _p_dphi_p[i]*_dphi[i]

         +  _p_dphi_f[indx]*_dphi_f[i]

         +  _p_lap_f[i]*_lap_f[i]

         +  _p_f[i]*f[i]
 
         +  0.5*(kprime*inverse_k)*(
               _p_p_p[i]*p[i]*p[i]

               _p_p_dR_f[i]*p[i]*_dR_f[i]

               _p_p_dphi_f[i]*p[i]*_dphi_f[i]

               _p_dR_f_dR_f[i]*_dR_f[i]*_dR_f[i] 

               _p_dR_f_dphi_f[i]*_dR_f[i]*_dphi_f[i]  

               _p_sphereX_f[i]*_sphereX_f[i]

               _p_dR_f_f[i]*_dR_f[i]*f[i]

               _p_dphi_f_f[i]*_dphi_f[i]*f[i] 

               _p_f_f[i]*f[i]*f[i] 

               _p_p_f[i]*p[i]*f[i]
            )

         +  _pre[i]*vprime
         ;
      }
   }
#endif
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
