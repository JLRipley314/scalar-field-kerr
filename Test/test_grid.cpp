#include "gtest/gtest.h"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <vector>

#include "../Source/cheb.hpp"
#include "../Source/sphere.hpp"
#include "../Source/grid.hpp"

/*==========================================================================*/
/* computes derivatives using Grid functions and analytically,
 * and takes one norm of difference of values. 
 */
std::vector<double> get_norm_diff(
      const size_t nx,
      const size_t nl,
      const size_t nm,
      const size_t nlat,
      const size_t nphi,
      const double Rmin,
      const double Rmax)
{
   const size_t n = nx*nlat*nphi;

   const double cl = Rmax;

   const double rl    = 0.1*(Rmax-Rmin) + Rmin;
   const double ru    = 0.5*(Rmax-Rmin) + Rmin;
   const double width = ru-rl;

   Grid grid(cl, Rmin, Rmax, nx, nl, nm, nlat, nphi);

   std::vector<double> v(n,0);
   std::vector<double> dr_v1( n,0);
   std::vector<double> d2r_v1(n,0);
   std::vector<double> dphi_v1(n,0);
   std::vector<double> X_v1(n,0);
   std::vector<double> lap_v1(n,0);
   std::vector<double> dr_v2(n,0);
   std::vector<double> d2r_v2(n,0);
   std::vector<double> dphi_v2(n,0); 
   std::vector<double> X_v2(n,0);
   std::vector<double> lap_v2(n,0); 

   for (size_t i=0; i<n; i++) {
      const std::vector<double> R_th_ph = grid.R_th_ph(i);
      const double r     = grid.R_to_r(R_th_ph[0]);
      const double theta = R_th_ph[1];
      const double phi   = R_th_ph[2];

      double bump = 0.0;

      if ((r<ru) && (r>rl)) {
         bump = 10*exp(-1.0*width/(r-rl))*exp(-2.0*width/(ru-r));
      }

      double rval = pow((r-rl)/width,4)*pow((ru-r)/width,4)*bump;

      const double der_rval = (
            4.0*pow((r-rl)/width,3)*pow((ru-r)/width,4) 
         -  4.0*pow((r-rl)/width,4)*pow((ru-r)/width,3)
         +      pow((r-rl)/width,2)*pow((ru-r)/width,4)
         -  2.0*pow((r-rl)/width,4)*pow((ru-r)/width,2)  
         )*(bump/width);

      v[i] = rval*(
            1.0 
         +  pow(sin(theta),2)/(1.0 + pow(cos(phi),2))
         );
      dr_v1[i] = der_rval*(
            1.0 
         +  pow(sin(theta),2)/(1.0 + pow(cos(phi),2))
         );
      dphi_v1[i] = rval*(
            pow(sin(theta),2)*sin(2*phi)/pow(1.0+pow(cos(phi),2),2)
         );
      X_v1[i] = rval*pow(sin(theta),2)*(
            pow(cos(theta),2)*pow(3.0+cos(2*phi),2)
         +  pow(sin(2*phi),2)
         )/pow(1.0 + pow(cos(phi),2),4)
         ;
      lap_v1[i] = rval*(
            31.0
         +  57.0*cos(2*theta)
         +  72.0*pow(cos(theta),2)*cos(2*phi)
         -  6.0*cos(4*phi)*pow(sin(theta),2)
         )/pow(3.0 + cos(2*phi),3)
         ;
   }
   grid.set_partial_r(     v,    dr_v2);
   grid.set_partial_r( dr_v1,   d2r_v1);
   grid.set_partial2_r(    v,   d2r_v2);
   grid.set_partial_phi(   v,  dphi_v2);
   grid.set_sphereX(       v,     X_v2);
   grid.set_spherical_lap( v,   lap_v2);

   double norm_dr   = 0;
   double norm_d2r  = 0;
   double norm_dphi = 0;
   double norm_X    = 0;
   double norm_lap  = 0;
   for (size_t i=0; i<n; i++) {
      norm_dr   += fabs(dr_v1[i]   - dr_v2[i]);
      norm_d2r  += fabs(d2r_v1[i]  - d2r_v2[i]);
      norm_dphi += fabs(dphi_v1[i] - dphi_v2[i]);
      norm_X    += fabs(X_v1[i]    - X_v2[i]);
      norm_lap  += fabs(lap_v1[i]  - lap_v2[i]);

   }
   norm_dr   /= n;
   norm_d2r  /= n;
   norm_dphi /= n;
   norm_X    /= n;
   norm_lap  /= n;

   std::vector<double> norms(5,0);

   norms[0] = norm_dr;
   norms[1] = norm_d2r;
   norms[2] = norm_dphi;
   norms[3] = norm_X;
   norms[4] = norm_lap;

   return norms;
}
/*==========================================================================*/
/* Adds random noise to function, filters, and returns total variation
 * of both 
 */
std::vector<double> get_total_variation(
      const size_t nx,
      const size_t nl,
      const size_t nm,
      const size_t nlat,
      const size_t nphi,
      const double Rmin,
      const double Rmax)
{
   const size_t n = nx*nlat*nphi;

   const double cl = Rmax;

   const double rl    = 0.1*(Rmax-Rmin) + Rmin;
   const double ru    = 0.5*(Rmax-Rmin) + Rmin;
   const double width = ru-rl;

   Grid grid(cl, Rmin, Rmax, nx, nl, nm, nlat, nphi);
   /*
    * noisy grid + filter of it
    */
   std::vector<double> v(n,0);
   std::vector<double> filter_v(n,0);

   for (size_t i=0; i<n; i++) {
      const std::vector<double> loc = grid.R_th_ph(i);
      const double r     = grid.R_to_r(loc[0]);
      const double theta = loc[1];
      const double phi   = loc[2];

      double bump = 0.0;

      if ((r<ru) && (r>rl)) {
         bump = exp(-1.0*width/(r-rl))*exp(-2.0*width/(ru-r));
      }

      double rval = pow((r-rl)/width,2)*pow((ru-r)/width,2)*bump;

      v[i] = rval*(
            1.0 
         +  pow(sin(theta),2)/(1.0 + pow(cos(phi),2))
         );
      v[i] *= (1.0 + (0.1*rand()/RAND_MAX));
      filter_v[i] = v[i];
   }
   grid.filter(filter_v);
   /*
    * total variation of grid and filtered version of it 
    */
   double tv_v        = 0;
   double tv_filter_v = 0;

   for (size_t ix=0; ix<nx-1;   ix++) {
   for (size_t it=0; it<nlat-1; it++) {
   for (size_t ip=0; ip<nphi-1; ip++) {
      tv_v += pow(
               pow(
                  v[grid.indx_R_th_ph(ix+1, it, ip)] 
               -  v[grid.indx_R_th_ph(ix,   it, ip)]
               ,2)
            +  pow(
                  v[grid.indx_R_th_ph(ix, it+1, ip)] 
               -  v[grid.indx_R_th_ph(ix, it,   ip)]
               ,2)
            +  pow(
                  v[grid.indx_R_th_ph(ix, it, ip+1)] 
               -  v[grid.indx_R_th_ph(ix, it, ip)]
               ,2)
            ,
            0.5); 

      tv_filter_v += pow(
               pow(
                  filter_v[grid.indx_R_th_ph(ix+1, it, ip)] 
               -  filter_v[grid.indx_R_th_ph(ix,   it, ip)]
               ,2)
            +  pow(
                  filter_v[grid.indx_R_th_ph(ix, it+1, ip)] 
               -  filter_v[grid.indx_R_th_ph(ix, it,   ip)]
               ,2)
            +  pow(
                  filter_v[grid.indx_R_th_ph(ix, it, ip+1)] 
               -  filter_v[grid.indx_R_th_ph(ix, it, ip)]
               ,2)
            ,
            0.5); 
   }
   }
   }
   std::vector<double> tvs(2,0);

   tvs[0] = tv_v;
   tvs[1] = tv_filter_v;

   return tvs;
}
/*==========================================================================*/
/* Testing application of cheb and sphere to 3d works 
 */
TEST(test_grid, test_dr_dphi_lap) {

   const size_t Rmin = 1;
   const size_t Rmax = 10;

   const size_t nx1   = 48;
   const size_t nl1   = 20;
   const size_t nm1   = 16;
   const size_t nlat1 = 2*nl1+2;
   const size_t nphi1 = nlat1;
   std::vector<double> norms1 = get_norm_diff(nx1, nl1, nm1, nlat1, nphi1, Rmin, Rmax);
   /* 
    * test convergence of the radial 
    */
   const size_t nx2   = 64;
   std::vector<double> norms2 = get_norm_diff(nx2, nl1, nm1, nlat1, nphi1, Rmin, Rmax);
   EXPECT_LT(norms2[0], norms1[0]);
   EXPECT_LT(norms2[1], norms1[1]);
   /* 
    * test convergence of the angular derivatives 
    */
   const size_t nl2   = 40;
   const size_t nm2   = 32;
   const size_t nlat2 = 2*nl2 + 4;
   const size_t nphi2 = nlat2;
   std::vector<double> norms3 = get_norm_diff(nx1, nl2, nm2, nlat2, nphi2, Rmin, Rmax);
   EXPECT_LT(norms3[2], norms1[2]);
   EXPECT_LT(norms3[3], norms1[3]);
   EXPECT_LT(norms3[4], norms1[4]);
}
/*==========================================================================*/
/* Testing grid filter is TVD 
 */
TEST(test_grid, test_tvd) {

   const size_t Rmin = 1;
   const size_t Rmax = 100;

   const size_t nx1   = 96;
   const size_t nl1   = 16;
   const size_t nm1   = 14;
   const size_t nlat1 = 2*nl1 + 2;
   const size_t nphi1 = nlat1;
   std::vector<double> tv1 = get_total_variation(nx1, nl1, nm1, nlat1, nphi1, Rmin, Rmax);

   const size_t nx2   = 128;
   const size_t nl2   = 24;
   const size_t nm2   = 20;
   const size_t nlat2 = 2*nl2 + 2;
   const size_t nphi2 = nlat2;
   std::vector<double> tv2 = get_total_variation(nx2, nl2, nm2, nlat2, nphi2, Rmin, Rmax);

   EXPECT_LT(tv1[1], tv1[0]);
   EXPECT_LT(tv2[1], tv2[0]);
}
