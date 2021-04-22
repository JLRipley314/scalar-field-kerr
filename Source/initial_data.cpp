#include <iostream>

#include "initial_data.hpp"
#include "grid.hpp"
#include "params.hpp"
#include "cheb.hpp"
#include "sphere.hpp"
/*==========================================================================*/
namespace ID 
{
/*==========================================================================*/
/* time symmetric compact bump for phi of angular dependence Y_{lm} */
/*==========================================================================*/
void compact_pulse(
      std::string initial_data_direction,
      std::vector<double> &f,
      std::vector<double> &p)
{
   const size_t nx   = Params::nx();
   const size_t nphi = Params::nphi();
   const size_t nlat = Params::nlat();

   const double amp = Params::amp();
   const double rl  = Params::rl();
   const double ru  = Params::ru();

   const int l_ang  = Params::l_ang();
   const int m_ang  = Params::m_ang();

   double max_val = 0.0;

   double width = fabs(ru-rl);
   if (fabs(width)<1e-14) width = 1e-14; /* don't want to divide by zero */

   std::vector<double> ylm = Sphere::compute_ylm(l_ang, m_ang);

   for (size_t ix=0; ix<nx-1; ix++) { /* do not include ix=nx-1 as r=infty there */
   for (size_t ip=0; ip<nphi; ip++) {
   for (size_t it=0; it<nlat; it++) {
      const std::vector<double> loc = Grid::r_th_ph(ix,it,ip);
      const double r = loc[0];
      double bump = 0.0;

      if ((r<ru) && (r>rl)) {
         bump = exp(-1.0*width/(r-rl))*exp(-2.0*width/(ru-r));
      }
      f[Grid::indx(ix,it,ip)] = pow((r-rl)/width,2)*pow((ru-r)/width,2)*bump; 

      const double partial_r_f = (
            2.0*   ((r-rl)/width  )*pow((ru-r)/width,2) 
         -  2.0*pow((r-rl)/width,2)*   ((ru-r)/width  )
         +                          pow((ru-r)/width,2)
         -  2.0*pow((r-rl)/width,2) 
         )*(bump/width)
         ;
      if (initial_data_direction=="ingoing") {
         p[Grid::indx(ix,it,ip)] =  partial_r_f;
      } else
      if (initial_data_direction=="outgoing") { 
         p[Grid::indx(ix,it,ip)] = -partial_r_f;
      } else
      if (initial_data_direction=="time_symmetric") { 
         p[Grid::indx(ix,it,ip)] = 0; 
      } else {
         std::cout
            <<"Using default initial data direction: time_symmetric."
            <<std::endl;
         p[Grid::indx(ix,it,ip)] = 0; 
      } 
      /*
       * give angular structure Y_{lm} 
       */
      f[Grid::indx(ix,it,ip)] *= ylm[Sphere::indx(it,ip)];
      p[Grid::indx(ix,it,ip)] *= ylm[Sphere::indx(it,ip)];

      if (fabs(f[Grid::indx(ix,it,ip)]) > max_val) {
         max_val = fabs(f[Grid::indx(ix,it,ip)]);
      }
   }
   }
   }
   /* 
    * rescale so initial amplitude is amp
    */
   for (size_t ix=0; ix<nx;   ix++) {
   for (size_t ip=0; ip<nphi; ip++) {
   for (size_t it=0; it<nlat; it++) {
      f[Grid::indx(ix,it,ip)] *= amp / max_val;
      p[Grid::indx(ix,it,ip)] *= amp / max_val;
   }
   }
   }
}
/*==========================================================================*/
} /* ID */
/*==========================================================================*/
