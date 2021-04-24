#include "gtest/gtest.h"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <vector>

#include "../Source/cheb.hpp"

/*==========================================================================*/
/* Testing we can go to/from Chebyshev space,
 * without changing the value of the array (to within truncation error). 
 */
TEST(test_cheb, to_and_from) {
   const double eps = 1e-14;
   const size_t nr = pow(2,5);
   Cheb cheb(nr, 0.0, 1.0);

   std::vector<double> po1(nr,0);
   std::vector<double> po2(nr,0);
   std::vector<double> ch( nr,0);

   /* 
    * fill in values
    */
   for (size_t i=0; i<nr; i++) {
      po1[i] = tan(cheb.pt(i));
   }

   cheb.to_ch(po1, ch);
   cheb.to_po(ch, po2);

   for (size_t i=0; i<nr; i++) {
      EXPECT_LT(fabs(po1[i] - po2[i]), eps);
   }
}
/*==========================================================================*/
/* Testing we can compute derivatives of simple functions 
 */
TEST(test_cheb, derivatives) {
   const double eps = 5e-14;
   const size_t nr = pow(2,5);
   Cheb cheb(nr, 0.0, 1.0);

   std::vector<double> v(  nr,0);
   std::vector<double> dv1(nr,0);
   std::vector<double> dv2(nr,0);
   /* 
    * fill in values
    */
   for (size_t i=0; i<nr; i++) {
      v[i]   = cheb.pt(i)*sin(cheb.pt(i));
      dv1[i] = cheb.pt(i)*cos(cheb.pt(i)) + sin(cheb.pt(i));
   }
   cheb.der(v, dv2);

   for (size_t i=0; i<nr; i++) {
      EXPECT_LT(fabs(dv1[i] - dv2[i]), eps);
   }
}
/*==========================================================================*/
/* Testing spectral filter is total variation diminishing
 * when acting on ``rough'' data. 
 */
TEST(test_cheb, filter_is_TVD) {
   const size_t nr = pow(2,5);
   Cheb cheb(nr, 0.0, 1.0);

   std::vector<double> po1(nr,0);
   std::vector<double> po2(nr,0);
   /* 
    * fill in values
    */
   for (size_t i=0; i<nr; i++) {
      po1[i] = tan(cheb.pt(i))*(1.0 + (0.1*rand()/RAND_MAX));
      po2[i] = po1[i];
   }
   cheb.filter(po2);
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
