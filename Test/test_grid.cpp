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
      const size_t nl,
      const size_t nx,
      const size_t nlat,
      const size_t nphi,
      const double Rmin,
      const double Rmax)
{
   const size_t n = nx*nlat*nphi;

   const double cl = Rmax;

   const double rl    = 0.1*(Rmax-Rmin) + Rmin;
   const double ru    = 0.2*(Rmax-Rmin) + Rmin;
   const double width = ru-rl;

   Sphere::init(nl, nlat, nphi);
   Cheb::init(nx, Rmin, Rmax);
   Grid::init(cl, nx, nlat, nphi);

   std::vector<double> v(n,0);
   std::vector<double> dr_v1(n,0);
   std::vector<double> dphi_v1(n,0);
   std::vector<double> lap_v1(n,0);
   std::vector<double> dr_v2(n,0);
   std::vector<double> dphi_v2(n,0); 
   std::vector<double> lap_v2(n,0); 

   for (size_t i=0; i<n; i++) {
      const std::vector<double> loc_r = Grid::r_th_ph(i);
      const double r     = loc_r[0];
      const double theta = loc_r[1];
      const double phi   = loc_r[2];

      double bump = 0.0;

      if ((r<ru) && (r>rl)) {
         bump = exp(-1.0*width/(r-rl))*exp(-2.0*width/(ru-r));
      }

      double rval = pow((r-rl)/width,2)*pow((ru-r)/width,2)*bump;

      double der_rval = (
         (2.0*((  (r-rl)/width) )*pow(((ru-r)/width),2))
      -  (2.0*pow((r-rl)/width,2)*(    (ru-r)/width   ))
      +  (1.0*(1.0              )*pow(((ru-r)/width),2))
      -  (2.0*pow((r-rl)/width,2)*(1.0                ))
      )*bump/width;


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
      lap_v1[i] = rval*(
            31.0
         +  57.0*cos(2*theta)
         +  72.0*pow(cos(theta),2)*cos(2*phi)
         -  6.0*cos(4*phi)*pow(sin(theta),2)
         )/pow(3.0 + cos(2*phi),3)
         ;
   }
   Grid::set_partial_r(    v,   dr_v2);
   Grid::set_partial_phi(  v, dphi_v2);
   Grid::set_spherical_lap(v,  lap_v2);

   double norm_dr   = 0;
   double norm_dphi = 0;
   double norm_lap  = 0;
   for (size_t i=0; i<n; i++) {
      norm_dr   += fabs(dr_v1[i]   - dr_v2[i]);
      norm_dphi += fabs(dphi_v1[i] - dphi_v2[i]);
      norm_lap  += fabs(lap_v1[i]  - lap_v2[i]);
   }
   norm_dr   /= n;
   norm_dphi /= n;
   norm_lap  /= n;

   std::vector<double> norms(3,0);

   norms[0] = norm_dr;
   norms[1] = norm_dphi;
   norms[2] = norm_lap;

   return norms;
}
/*==========================================================================*/
/* Testing application of cheb and sphere to 3d works 
 */
TEST(grid_test, test_dr_dphi_lap) {

   const size_t Rmin = 1;
   const size_t Rmax = 10;

   const size_t nx1   = 48;
   const size_t nl1   = 16;
   const size_t nlat1 = 2*nl1 + 2;
   const size_t nphi1 = nlat1;
   std::vector<double> norms1 = get_norm_diff(nl1, nx1, nlat1, nphi1, Rmin, Rmax);

   const size_t nx2   = 64;
   const size_t nl2   = 24;
   const size_t nlat2 = 2*nl2 + 2;
   const size_t nphi2 = nlat2;
   std::vector<double> norms2 = get_norm_diff(nl2, nx2, nlat2, nphi2, Rmin, Rmax);

   EXPECT_LT(norms2[0], norms1[0]);
   EXPECT_LT(norms2[1], norms1[1]);
   EXPECT_LT(norms2[2], norms1[2]);
}
