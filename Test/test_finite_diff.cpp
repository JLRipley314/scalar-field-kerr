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
std::vector<double> norm_derivatives(const size_t n) {
   FD fd(n, 0.0, 1.0);

   std::vector<double> v(n,0);
   std::vector<double> dv1(n,0);
   std::vector<double> dv2(n,0);
   std::vector<double> ddv1(n,0);
   std::vector<double> ddv2(n,0);
   /* 
    * fill in values
    */
   for (size_t i=0; i<n; i++) {
      v[i]    =   fd.pt(i)*sin(fd.pt(i));
      dv1[i]  =   fd.pt(i)*cos(fd.pt(i)) +     sin(fd.pt(i));
      ddv1[i] = - fd.pt(i)*sin(fd.pt(i)) + 2.0*cos(fd.pt(i));
   }
   fd.der( v,  dv2);
   fd.der2(v, ddv2);

   double norm_der1 = 0;
   double norm_der2 = 0;
   for (size_t i=0; i<n; i++) {
      norm_der1 += fabs(dv1[i]  -  dv2[i]);
      norm_der2 += fabs(ddv1[i] - ddv2[i]);
   }
   norm_der1 /= n;
   norm_der2 /= n;
   std::vector<double> out = {norm_der1, norm_der2};
   return out;
}
/*==========================================================================*/
/* Testing convergence of derivatives to true answer (6th order convergence) 
 */
TEST(test_finite_diff, derivatives) {
   std::vector<double> norms_n1 = norm_derivatives(33);
   std::vector<double> norms_n2 = norm_derivatives(65);

   std::cout
      <<std::setw(16)<<norms_n1[0]
      <<std::setw(16)<<norms_n2[0]
      <<std::setw(16)<<norms_n1[1]
      <<std::setw(16)<<norms_n2[1]
      <<std::endl;

   EXPECT_LT(norms_n2[0], (1./64.)*norms_n1[0]);
   EXPECT_LT(norms_n2[1], (1./64.)*norms_n1[1]);
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
