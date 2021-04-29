#include "gtest/gtest.h"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <vector>

#include "../Source/finite_diff.hpp"

/*==========================================================================*/
/* Testing we can compute derivatives of simple functions 
 */
TEST(test_finite_diff, derivatives) {
   const double eps = 7e-14;
   const size_t nr = pow(2,7);
   FD fd(nr, 0.0, 1.0);

   std::vector<double> v(  nr,0);
   std::vector<double> dv1(nr,0);
   std::vector<double> dv2(nr,0);
   /* 
    * fill in values
    */
   for (size_t i=0; i<nr; i++) {
      v[i]   = fd.pt(i)*sin(fd.pt(i));
      dv1[i] = fd.pt(i)*cos(fd.pt(i)) + sin(fd.pt(i));
   }
   fd.der(v, dv2);

   for (size_t i=0; i<nr; i++) {
      EXPECT_LT(fabs(dv1[i] - dv2[i]), eps);
   }
}
/*==========================================================================*/
/* Testing Kreiss-Oliger filter is total variation diminishing
 * when acting on ``rough'' data. 
 */
TEST(test_finite_diff, filter_is_TVD) {
   const size_t nr = pow(2,7);
   FD fd(nr, 0.0, 1.0);

   std::vector<double> po1(nr,0);
   std::vector<double> po2(nr,0);
   /* 
    * fill in values
    */
   for (size_t i=0; i<nr; i++) {
      po1[i] = tan(fd.pt(i))*(1.0 + (0.1*rand()/RAND_MAX));
      po2[i] = po1[i];
   }
   fd.filter(po2);
   /* 
    * compute total variations 
    */
   double tv1 = 0;
   double tv2 = 0;
   for (size_t i=0; i<nr-1; i++) {
      tv1 += fabs(po1[i+1] - po1[i]);
      tv2 += fabs(po2[i+1] - po2[i]);
   }
   EXPECT_LT(tv2, tv1);
}
