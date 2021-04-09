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
   std::vector<double> dv1(n,0);
   std::vector<double> dv2(n,0);

   for (size_t i=0; i<nx*nlat*nphi; i++) {
      const std::vector<double> loc_r = Grid::r_th_ph(i);
      const double r = loc_r[0];

      v[i]   =  1.0/r;
      dv1[i] = 0;
      dv2[i] = -1.0/pow(r,2);
   }
   Grid::set_partial_r(v, dv1);

   double norm = 0;
   for (size_t i=0; i<nx*nlat*nphi; i++) {
      norm += fabs(dv1[i]-dv2[i]);
   }
   norm /= nx*nlat*nphi;

   std::cout<<"norm "<<norm<<std::endl;
}
