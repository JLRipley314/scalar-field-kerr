#include "gtest/gtest.h"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

#include "../Source/cheb.hpp"

const double eps = 1e-14;
/*
 * Testing we can go to/from Chebyshev space,
 * without changing the value of the array (to within truncation error). 
 */
TEST(ChebTest, toAndFrom) {
   const size_t nr = pow(2,5);
   Cheb::init(nr, 0.0, 1.0);

   std::vector<double> po1(nr,0);
   std::vector<double> po2(nr,0);
   std::vector<double> ch( nr,0);

   /* 
    * fill in values
    */
   for (size_t i=0; i<nr; i++) {
      po1[i] = Cheb::pt(i);
   }

   Cheb::to_ch(po1, ch);
   Cheb::to_po(ch, po2);

   for (size_t i=0; i<nr; i++) {
      std::cout
         <<std::setw(16)<<i
         <<std::setw(16)<<ch[i]
         <<std::setw(16)<<po1[i]-po2[i]
         <<std::endl;
//      EXPECT_TRUE(abs(po1[i] - po2[i]) < eps);
   }
   Cheb::cleanup();
}
