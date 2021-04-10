#include <cassert>
#include <iostream>

#include "scalar_eom.hpp"
#include "params.hpp"
#include "grid.hpp"
/*==========================================================================*/
namespace Eom {
/*==========================================================================*/
namespace {
   size_t _n;
   double _constraint_damping;

   std::vector<double> _p_f1;
   std::vector<double> _p_f2;
   std::vector<double> _p_f3;
   std::vector<double> _p_p;
   std::vector<double> _p_q;
   std::vector<double> _p_dr_p;
   std::vector<double> _p_dr_q;
   std::vector<double> _p_dphi_q;
   std::vector<double> _p_lap_f;
}
/*==========================================================================*/
void init()
{
   _n = Params::nx_nlat_nphi();
   _constraint_damping = Params::constraint_damping();

   _p_f1.resize(_n,0);
   _p_f2.resize(_n,0);
   _p_f3.resize(_n,0);
   _p_p.resize(_n,0);
   _p_q.resize(_n,0);
   _p_dr_p.resize(_n,0);
   _p_dr_q.resize(_n,0);
   _p_dphi_q.resize(_n,0);
   _p_lap_f.resize(_n,0);

   const size_t nx   = Params::nx();
   const size_t nlat = Params::nlat();
   const size_t nphi = Params::nphi();

   const double cl  = Params::cl();

   const double V2 = Params::V_2();
   const double V3 = Params::V_3();
   const double V4 = Params::V_4();


   for (size_t ix=0; ix<nx;   ix++) {
   for (size_t ip=0; ip<nphi; ip++) {
   for (size_t it=0; it<nlat; it++) {
      const size_t indx = Grid::indx(ix,it,ip);
      
      const double inv_r = Cheb::pt(ix)/pow(cl,2);

      const double m = Params::bh_mass();
      const double a = Params::bh_spin()/m;

      /* divide by r^2 to reduce infty/infy type errors 
       * in computing coefficients. */
      const double Sigma = 1.0 + pow(inv_r,2)*pow(a,2)*pow(cos(Sphere::theta(it)),2); 
      const double Delta = 1.0 + pow(inv_r,2)*pow(a,2) - inv_r*(2.0*m); 

      _p_f1[indx] = V2;
      _p_f2[indx] = V3/2.0;
      _p_f3[indx] = V4/6.0;

      _p_p[indx] = (2.0*m/Sigma) * pow(inv_r,2);
      _p_q[indx] = 2.0*(inv_r - (m*pow(inv_r,2)))/Sigma;

      _p_dr_p[indx] = 4.0*m*inv_r/Sigma;
      _p_dr_q[indx] = (Delta/Sigma);

      _p_dphi_q[indx] = (2.0*a/Sigma) * pow(inv_r,2);

      _p_lap_f[indx] = (1.0/Sigma) * pow(inv_r,2);

      const double pre = 1.0 + (2.0*m*inv_r/Sigma);

      _p_f1[indx] /= pre;
      _p_f2[indx] /= pre;
      _p_f3[indx] /= pre;

      _p_p[indx] /= pre;
      _p_q[indx] /= pre;

      _p_dr_p[indx] /= pre;
      _p_dr_q[indx] /= pre;

      _p_dphi_q[indx] /= pre;

      _p_lap_f[indx] /= pre;
   }
   }
   }
}
/*==========================================================================*/
void set_k(
      const std::vector<double> &f,
      const std::vector<double> &p,
      const std::vector<double> &q,
      std::vector<double> &dr_f,
      std::vector<double> &lap_f,
      std::vector<double> &dr_p,
      std::vector<double> &dr_q,
      std::vector<double> &dphi_q,
      std::vector<double> &f_k,
      std::vector<double> &p_k,
      std::vector<double> &q_k
      )
{
   Grid::set_partial_r(f, dr_f);
   Grid::set_partial_r(p, dr_p);
   Grid::set_partial_r(q, dr_q);

   Grid::set_partial_phi(f, dphi_q);

   Grid::set_spherical_lap(f, lap_f);

   for (size_t i=0; i<_n; i++) {
      f_k[i] = 
         p[i]
      ;
      p_k[i] = 
      +  _p_f1[i]*f[i] 
      +  _p_f2[i]*pow(f[i],2) 
      +  _p_f3[i]*pow(f[i],3) 

      +  _p_p[i]*p[i] 
      +  _p_q[i]*q[i]

      +  _p_dr_p[i]*dr_p[i] 
      +  _p_dr_q[i]*dr_q[i]

      +  _p_dphi_q[i]*dphi_q[i]

      +  _p_lap_f[i]*lap_f[i]
      ;
      q_k[i] = 
         dr_p[i] 
      +  _constraint_damping*(q[i] - dr_f[i])
      ;
   }
}
/*==========================================================================*/
/* Fourth order Runge-Kutta integrator */
/*==========================================================================*/
void time_step(Field &f, Field &p, Field &q)
{
   const double dt = Params::dt();
   const size_t n = Params::nx_nlat_nphi();

//   Grid::filter(f.n);
//   Grid::filter(p.n);
//   Grid::filter(q.n);

   std::vector<double> dr_f(  n);
   std::vector<double> lap_f( n);
   std::vector<double> dr_p(  n);
   std::vector<double> dr_q(  n);
   std::vector<double> dphi_q(n);
   /*--------------------------------------------*/
   set_k(
         f.n, p.n, q.n, 
         dr_f, lap_f, dr_p, dr_q, dphi_q, 
         f.k1, p.k1, q.k1
      );
   for (size_t i=0; i<n; i++) {
      f.l2[i] = f.n[i] + 0.5*dt*f.k1[i];
      p.l2[i] = p.n[i] + 0.5*dt*p.k1[i];
      q.l2[i] = q.n[i] + 0.5*dt*q.k1[i];
   }
   /*--------------------------------------------*/
   set_k(
         f.l2, p.l2, q.l2, 
         dr_f, lap_f, dr_p, dr_q, dphi_q, 
         f.k2, p.k2, q.k2
      );
   for (size_t i=0; i<n; i++) {
      f.l3[i] = f.n[i] + 0.5*dt*f.k2[i];
      p.l3[i] = p.n[i] + 0.5*dt*p.k2[i];
      q.l3[i] = q.n[i] + 0.5*dt*q.k2[i];
   }
   /*--------------------------------------------*/
   set_k(
         f.l3, p.l3, q.l3, 
         dr_f, lap_f, dr_p, dr_q, dphi_q, 
         f.k3, p.k3, q.k3
      );
   for (size_t i=0; i<n; i++) {
      f.l4[i] = f.n[i] + dt*f.k3[i];
      p.l4[i] = p.n[i] + dt*p.k3[i];
      q.l4[i] = q.n[i] + dt*q.k3[i];
   }
   /*--------------------------------------------*/
   set_k(
         f.l4, p.l4, q.l4,
         dr_f, lap_f, dr_p, dr_q, dphi_q, 
         f.k4, p.k4, q.k4
      );
   for (size_t i=0; i<n; i++) {
      f.np1[i] = f.n[i] + (dt/6.0)*(f.k1[i] + 2.0*f.k2[i] + 2.0*f.k3[i] + f.k4[i]);
      p.np1[i] = p.n[i] + (dt/6.0)*(p.k1[i] + 2.0*p.k2[i] + 2.0*p.k3[i] + p.k4[i]);
      q.np1[i] = q.n[i] + (dt/6.0)*(q.k1[i] + 2.0*q.k2[i] + 2.0*q.k3[i] + q.k4[i]);
   }
}
/*==========================================================================*/
} /* Eom */
