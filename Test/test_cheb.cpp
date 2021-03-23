#include "gtest/gtest.h"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <vector>

#include "../Source/cheb.hpp"

const double eps = 1e-14;
/*==========================================================================*/
/* Testing we can go to/from Chebyshev space,
 * without changing the value of the array (to within truncation error). 
 */
TEST(cheb_test, to_and_from) {
   const size_t nr = pow(2,5);
   Cheb::init(nr, 0.0, 1.0);

   std::vector<double> po1(nr,0);
   std::vector<double> po2(nr,0);
   std::vector<double> ch( nr,0);

   /* 
    * fill in values
    */
   for (size_t i=0; i<nr; i++) {
      po1[i] = tan(Cheb::pt(i));
   }

   Cheb::to_ch(po1, ch);
   Cheb::to_po(ch, po2);

   for (size_t i=0; i<nr; i++) {
      EXPECT_TRUE(abs(po1[i] - po2[i]) < eps);
   }
   Cheb::cleanup();
}
/*==========================================================================*/
/* Testing spectral filter is total variation diminishing
 * when acting on ``rough'' data. 
 */
TEST(cheb_test, filter_is_TVD) {
   const size_t nr = pow(2,5);
   Cheb::init(nr, 0.0, 1.0);

   std::vector<double> po1(nr,0);
   std::vector<double> po2(nr,0);
   /* 
    * fill in values
    */
   for (size_t i=0; i<nr; i++) {
      po1[i] = tan(Cheb::pt(i))*(1.0 + (0.1*rand()/RAND_MAX));
      po2[i] = po1[i];
   }
   Cheb::filter(po2);
   /* 
    * compute total variations 
    */
   double tv1 = 0;
   double tv2 = 0;
   for (size_t i=1; i<nr-1; i++) {
      tv1 += abs(po1[i+1] - po1[i]);
      tv2 += abs(po2[i+1] - po2[i]);
   }
   EXPECT_TRUE(tv2 < tv1);
   Cheb::cleanup();
}
