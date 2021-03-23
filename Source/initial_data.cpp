#include "initial_data.hpp"
#include "arr.hpp"
#include "params.hpp"
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
      const std::vector<double> &rv,
      std::vector<double> &f,
      std::vector<double> &p,
      std::vector<double> &q)
{
   const size_t nx   = Params::nx();
   const size_t nlat = Params::nlat();
   const size_t nphi = Params::nphi();

   const double cl  = Params::cl();
   const double amp = Params::amp();
   const double rl  = Params::r_l();
   const double ru  = Params::r_u();

   const int l_ang  = Params::l_ang();
   const int m_ang  = Params::m_ang();

   double max_val = 0.0;

   double width = abs(rl-ru);
   if (abs(width)<1e-14) width = 1e-14; /* don't want to divide by zero */

   std::vector<double> ylm = Sphere::compute_ylm(l_ang, m_ang);

   for (size_t ix=0; ix<nx;   ix++) {
   for (size_t it=0; it<nlat; it++) {
   for (size_t ip=0; ip<nphi; ip++) {
      double r    = rv[ix];
      double bump = 0.0;

      if ((r<ru) && (r>rl)) {
         bump = exp(-1.0*width/(r-rl))*exp(-2.0*width/(ru-r));
      }
      f[indx(ix,it,ip)] = pow((r-rl)/width,2)*pow((ru-r)/width,2)*bump; 

      q[indx(ix,it,ip)] = (
         (2.0*(((r-rl)/width)   )*pow(((ru-r)/width),2))
      -  (2.0*pow((r-rl)/width,2)*((ru-r)/width       ))
      +  (1.0*(1.0              )*pow(((ru-r)/width),2))
      -  (2.0*pow((r-rl)/width,2)*(1.0                ))
      )*bump/width;

      q[indx(ix,it,ip)] *= -pow(r/cl,2);
      /*
       * time symmetric for now
       */ 
      p[indx(ix,it,ip)] = 0.0;
      /*
       * give angular structure Y_{lm} 
       */
      f[indx(ix,it,ip)] *= ylm[indx_Sph(it,ip)];
      q[indx(ix,it,ip)] *= ylm[indx_Sph(it,ip)];
      p[indx(ix,it,ip)] *= ylm[indx_Sph(it,ip)];

      if (abs(f[indx(ix,it,ip)]) > max_val) {
         max_val = abs(f[indx(ix,it,ip)]);
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
      p[indx(ix,it,ip)] *= amp / max_val;
      p[indx(ix,it,ip)] *= amp / max_val;
      p[indx(ix,it,ip)] *= amp / max_val; 
   }
   }
   }
}
/*==========================================================================*/
} /* ID */
/*==========================================================================*/
