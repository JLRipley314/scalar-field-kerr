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
/* Testing application of cheb and sphere to 3d works 
 */
TEST(grid_test, test_dr) {
   const size_t nl   = pow(2,5);
   const size_t nlat = nl + 2;
   const size_t nphi = 3*nl;

   const size_t nx   = pow(2,5);
   const double Rmin = 1;
   const double Rmax = 10;
   const double cl   = Rmax;

   const size_t n = nx*nlat*nphi;

   Sphere::init(nl, nlat, nphi);
   Cheb::init(nx, Rmin, Rmax);
   Grid::init(cl, nx, nlat, nphi);

   std::vector<double> v(  n,0);

   std::vector<double> dr_v1(n,0);
   std::vector<double> dr_v2(n,0);

   std::vector<double> dphi_v1(n,0);
   std::vector<double> dphi_v2(n,0);

   std::vector<double> lap_v1(n,0);
   std::vector<double> lap_v2(n,0);

   for (size_t i=0; i<nx*nlat*nphi; i++) {
      const std::vector<double> loc_r = Grid::r_th_ph(i);
      const double r = loc_r[0];
      const double theta = loc_r[1];
      const double phi   = loc_r[2];

      v[i] = (1.0/r)*(
            1.0 
         +  pow(sin(theta),2)*pow(cos(phi),2)
         );
      dr_v1[i] = (-1.0/pow(r,2))*(
            1.0 
         +  pow(sin(theta),2)*pow(cos(phi),2)
         );
      dphi_v1[i] = (1.0/r)*(
         -  2.0*pow(sin(theta),2)*cos(phi)*sin(phi)
         );
      lap_v1[i] = (1.0/r)*(
         0.5*(
               1.0
            +  3.0*cos(2.0*theta)
            -  6.0*cos(2.0*phi)*pow(sin(theta),2)
            ) 
         );
   }
   Grid::set_partial_r(    v,   dr_v2);
   Grid::set_partial_phi(  v, dphi_v2);
   Grid::set_spherical_lap(v,  lap_v2);

   double norm_dr   = 0;
   double norm_dphi = 0;
   double norm_lap  = 0;
   for (size_t i=0; i<nx*nlat*nphi; i++) {
      norm_dr   += fabs(dr_v1[i]   - dr_v2[i]);
      norm_dphi += fabs(dphi_v1[i] - dphi_v2[i]);
      norm_lap  += fabs(lap_v1[i]  - lap_v2[i]);
   }
   norm_dr   /= nx*nlat*nphi;
   norm_dphi /= nx*nlat*nphi;
   norm_lap  /= nx*nlat*nphi;

   std::cout<<"norm_dr:   "<<norm_dr<<std::endl;
   std::cout<<"norm_dphi: "<<norm_dphi<<std::endl;
   std::cout<<"norm_lap:  "<<norm_lap<<std::endl;
}
