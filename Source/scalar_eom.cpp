#include "scalar_eom.hpp"
#include "arr.hpp"
#include "cheb.hpp"
#include "sphere.hpp"
#include "params.hpp"
/*==========================================================================*/
namespace Eom {
/*==========================================================================*/
namespace {
   std::vector<double> p_f;
   std::vector<double> p_p;
   std::vector<double> p_q;
   std::vector<double> p_dr_p;
   std::vector<double> p_dr_q;
   std::vector<double> p_dphi_q;
   std::vector<double> p_lap_f;
}
/*==========================================================================*/
void init()
{
   const size_t ntotal = Params::nx_nphi_nlat();
   p_f.resize(ntotal,0);
   p_p.resize(ntotal,0);
   p_q.resize(ntotal,0);
   p_dr_p.resize(ntotal,0);
   p_dr_q.resize(ntotal,0);
   p_dphi_q.resize(ntotal,0);
   p_lap_f.resize(ntotal,0);

   const size_t nx   = Params::nx();
   const size_t nlat = Params::nlat();
   const size_t nphi = Params::nphi();

   const double cl  = Params::cl();

   for (size_t ix=0; ix<nx-1; ix++) { /* do not include ix=nx-1 as r=infty there */
   for (size_t ip=0; ip<nphi; ip++) {
   for (size_t it=0; it<nlat; it++) {
      const size_t indx     = Arr3d::indx(  ix,ip,it);
//      const size_t indx_Sph = Sphere::indx_Sph(ip,it);

      double r = pow(cl,2)/Cheb::pt(ix);

      double m = Params::bh_mass();
      double a = Params::bh_spin()/m;

      /* divide by r^2 to reduce infty/infy type errors 
       * in computing coefficients. */
      double Sigma = 1.0 + pow(a/r,2)*pow(cos(Sphere::theta(it)),2); 
      double Delta = 1.0 + pow(a/r,2) - (2.0*m/r); 

      p_f[indx] = 0.0;
      p_p[indx] = (2.0*m/Sigma) / pow(r,2);
      p_q[indx] = 2.0*((1.0/r) - (m/pow(r,2)))/Sigma;

      p_dr_p[indx] = 2.0*m*(1.0/r)/Sigma;
      p_dr_q[indx] = (Delta/Sigma) / pow(r,2);

      p_dphi_q[indx] = (2.0*a/Sigma) / pow(r,2);

      p_lap_f[indx] = (1.0/Sigma) / pow(r,2);

      double pre = 1.0 + (2.0*m*(1.0/r)/Sigma);

      p_f[indx] /= pre;
      p_p[indx] /= pre;
      p_q[indx] /= pre;

      p_dr_p[indx] /= pre;
      p_dr_q[indx] /= pre;

      p_dphi_q[indx] /= pre;

      p_lap_f[indx] /= pre;
   }
   }
   }
}
/*==========================================================================*/
void cleanup()
{
return;
}
/*==========================================================================*/
void set_partial_r(const std::vector<double> &v, std::vector<double> dv)
{
   std::vector<double> inter(   Params::nx());
   std::vector<double> inter_dv(Params::nx());

   for (size_t ip=0; ip<Params::nphi(); ip++) {
   for (size_t it=0; it<Params::nlat(); it++) {
      Arr3d::get_row1(ip, it, v, inter); 

      Cheb::der(inter, inter_dv);

      Arr3d::set_row1(ip, it, inter_dv, dv); 
   }
   }
}
/*==========================================================================*/
void set_spherical_lap(const std::vector<double> &v, std::vector<double> ddv)
{
   std::vector<double> inter(    Params::nlat()*Params::nphi());
   std::vector<double> inter_ddv(Params::nlat()*Params::nphi());

   for (size_t ix=0; ix<Params::nx(); ix++) {
      Arr3d::get_row23(ix, v, inter); 

      Sphere::laplace_beltrami(inter, inter_ddv);

      Arr3d::set_row23(ix, inter_ddv, ddv); 
   }
}
/*==========================================================================*/
void set_partial_phi(const std::vector<double> &v, std::vector<double> dv)
{
   std::vector<double> inter(   Params::nlat()*Params::nphi());
   std::vector<double> inter_dv(Params::nlat()*Params::nphi());

   for (size_t ix=0; ix<Params::nx(); ix++) {
      Arr3d::get_row23(ix, v, inter); 

      Sphere::partial_phi(inter, inter_dv);

      Arr3d::set_row23(ix, inter_dv, dv); 
   }
}
/*==========================================================================*/
/* Low pass filter in spectral space */
/*==========================================================================*/
void filter(std::vector<double> &v)
{
   std::vector<double> inter_radial(Params::nx());
   std::vector<double> inter_sphere(Params::nlat()*Params::nphi());

   for (size_t ix=0; ix<Params::nx(); ix++) {
      Arr3d::get_row23(ix, v, inter_sphere); 
      Sphere::filter(inter_sphere);
      Arr3d::set_row23(ix, inter_sphere, v); 
   }

   for (size_t ip=0; ip<Params::nphi(); ip++) {
   for (size_t it=0; it<Params::nlat(); it++) {
      Arr3d::get_row1(ip, it, v, inter_radial); 
      Cheb::filter(inter_radial);
      Arr3d::set_row1(ip, it, inter_radial, v); 
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
   set_partial_r(f, dr_f);
   set_partial_r(p, dr_p);
   set_partial_r(q, dr_q);

   set_partial_phi(f, dphi_q);

   set_spherical_lap(f, lap_f);

   for (size_t i=0; i<Params::nx_nphi_nlat(); i++) {
      f_k[i] = p[i];
      p_k[i] = 
         p_f[i]*f[i] 
      +  p_p[i]*p[i] 
      +  p_q[i]*q[i]

      +  p_dr_p[i]*dr_p[i] 
      +  p_dr_q[i]*dr_q[i]

      +  p_dphi_q[i]*dphi_q[i]

      +  p_lap_f[i]*lap_f[i]
      ;
      q_k[i] = dr_p[i];
   }
}
/*==========================================================================*/
/* Fourth order Runge-Kutta integrator */
/*==========================================================================*/
void time_step(Field &f, Field &p, Field &q)
{
   const double dt = Params::dt();
   const size_t n = Params::nx_nphi_nlat();

   filter(f.n);
   filter(p.n);
   filter(q.n);

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
