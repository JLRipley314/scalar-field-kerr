#include "initial_data.hpp"
#include "arr.hpp"
#include "params.hpp"
#include "cheb.hpp"
#include "sphere.hpp"

using Arr3d::indx;
using Sphere::indx_Sph;
/*==========================================================================*/
namespace ID 
{
/*==========================================================================*/
/* time symmetric compact bump for phi of angular dependence Y_{lm} */
/*==========================================================================*/
void ingoing_pulse(
      std::vector<double> &f,
      std::vector<double> &p,
      std::vector<double> &q)
{
   const size_t nx   = Params::nx();
   const size_t nlat = Params::nlat();
   const size_t nphi = Params::nphi();

   const double cl  = Params::cl();
   const double amp = Params::amp();
   const double rl  = Params::rl();
   const double ru  = Params::ru();

   const int l_ang  = Params::l_ang();
   const int m_ang  = Params::m_ang();

   double max_val = 0.0;

   double width = abs(rl-ru);
   if (abs(width)<1e-14) width = 1e-14; /* don't want to divide by zero */

   std::vector<double> ylm = Sphere::compute_ylm(l_ang, m_ang);

   for (size_t ix=0; ix<nx-1; ix++) { /* do not include ix=nx-1 as r=infty there */
   for (size_t ip=0; ip<nphi; ip++) {
   for (size_t it=0; it<nlat; it++) {
      double r    = pow(cl,2)/Cheb::pt(ix);
      double bump = 0.0;

      if ((r<ru) && (r>rl)) {
         bump = exp(-1.0*width/(r-rl))*exp(-2.0*width/(ru-r));
      }
      f[indx(ix,ip,it)] = pow((r-rl)/width,2)*pow((ru-r)/width,2)*bump; 

      q[indx(ix,ip,it)] = (
         (2.0*(((r-rl)/width)   )*pow(((ru-r)/width),2))
      -  (2.0*pow((r-rl)/width,2)*((ru-r)/width       ))
      +  (1.0*(1.0              )*pow(((ru-r)/width),2))
      -  (2.0*pow((r-rl)/width,2)*(1.0                ))
      )*bump/width;

      q[indx(ix,ip,it)] *= -pow(r/cl,2);
      /*
       * time symmetric for now
       */ 
      p[indx(ix,ip,it)] = 0.0;
      /*
       * give angular structure Y_{lm} 
       */
      f[indx(ix,ip,it)] *= ylm[indx_Sph(ip,it)];
      q[indx(ix,ip,it)] *= ylm[indx_Sph(ip,it)];
      p[indx(ix,ip,it)] *= ylm[indx_Sph(ip,it)];

      if (abs(f[indx(ix,ip,it)]) > max_val) {
         max_val = abs(f[indx(ix,ip,it)]);
      }
   }
   }
   }
   /* 
    * rescale so initial amplitude is amp
    */
   for (size_t ix=0; ix<nx;   ix++) {
   for (size_t it=0; it<nlat; it++) {
   for (size_t ip=0; ip<nphi; ip++) {
      p[indx(ix,ip,it)] *= amp / max_val;
      p[indx(ix,ip,it)] *= amp / max_val;
      p[indx(ix,ip,it)] *= amp / max_val; 
   }
   }
   }
}
/*==========================================================================*/
} /* ID */
/*==========================================================================*/
