#include <iostream>
#include <iomanip>

#include "initial_data.hpp"
/*==========================================================================*/
namespace ID 
{
/*==========================================================================*/
/* time symmetric compact bump for phi of angular dependence Y_{lm} */
/*==========================================================================*/
void compact_pulse(
      const Params &params,
      const Grid &grid,
      std::vector<double> &f,
      std::vector<double> &p)
{
   std::cout<<"Setting compact_pulse initial data"<<std::endl;

   const size_t nx   = grid.nx();
   const size_t nphi = grid.nphi();
   const size_t nlat = grid.nlat();

   double max_val = 0.0;

   const std::string initial_data_direction = params.initial_data_direction(); 

   const double amp = params.amp();
   const double rl  = params.rl();
   const double ru  = params.ru();

   const int l_ang = params.l_ang(); 
   const int m_ang = params.m_ang();

   double width = fabs(ru-rl);
   if (fabs(width)<1e-14) width = 1e-14; /* don't want to divide by zero */

   std::vector<double> ylm = grid.compute_ylm(l_ang, m_ang);

   for (size_t ix=0; ix<nx-1; ix++) { /* do not include ix=nx-1 as r=infty there */
   for (size_t ip=0; ip<nphi; ip++) {
   for (size_t it=0; it<nlat; it++) {
      const std::vector<double> loc = grid.r_th_ph(ix,it,ip);

      const double r = loc[0];
      double bump = 0.0;

      if ((r<ru) && (r>rl)) {
         bump = exp(-1.0*width/(r-rl))*exp(-2.0*width/(ru-r));
      }
      f[grid.indx_r_th_ph(ix,it,ip)] = pow((r-rl)/width,4)*pow((ru-r)/width,4)*bump; 

      const double partial_r_f = (
            4.0*pow((r-rl)/width,3)*pow((ru-r)/width,4) 
         -  4.0*pow((r-rl)/width,4)*pow((ru-r)/width,3)
         +      pow((r-rl)/width,2)*pow((ru-r)/width,4)
         -  2.0*pow((r-rl)/width,4)*pow((ru-r)/width,2)  
         )*(bump/width)
         ;
      if (initial_data_direction=="ingoing") {
         p[grid.indx_r_th_ph(ix,it,ip)] =  partial_r_f;
      } else
      if (initial_data_direction=="outgoing") { 
         p[grid.indx_r_th_ph(ix,it,ip)] = -partial_r_f;
      } else
      if (initial_data_direction=="time_symmetric") { 
         p[grid.indx_r_th_ph(ix,it,ip)] = 0; 
      } else {
         std::cout
            <<"Using default initial data direction: time_symmetric."
            <<std::endl;
         p[grid.indx_r_th_ph(ix,it,ip)] = 0; 
      } 
      /*
       * give angular structure Y_{lm} 
       */
      f[grid.indx_r_th_ph(ix,it,ip)] *= ylm[grid.indx_th_ph(it,ip)];
      p[grid.indx_r_th_ph(ix,it,ip)] *= ylm[grid.indx_th_ph(it,ip)];

      if (fabs(f[grid.indx_r_th_ph(ix,it,ip)]) > max_val) {
         max_val = fabs(f[grid.indx_r_th_ph(ix,it,ip)]);
      }
   }
   }
   }
   /* 
    * rescale so initial amplitude is amp
    */
   if (fabs(max_val)>1e-16) {
      for (size_t ix=0; ix<nx;   ix++) {
      for (size_t ip=0; ip<nphi; ip++) {
      for (size_t it=0; it<nlat; it++) {
         f[grid.indx_r_th_ph(ix,it,ip)] *= amp / max_val;
         p[grid.indx_r_th_ph(ix,it,ip)] *= amp / max_val;
      }
      }
      }
   }
   std::cout<<"Finished setting compact_pulse initial data"<<std::endl;
}
/*==========================================================================*/
} /* ID */
/*==========================================================================*/
