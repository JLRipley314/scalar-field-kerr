#include "gtest/gtest.h"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

#include "../Source/sphere.hpp"

const double eps = 1e-14;
/*
 * Testing we can go to/from spherical harmonic space, 
 * without changing the value of the array (to within truncation error). 
 */
TEST(SphereTest, toAndFrom) {
   const size_t nl   = pow(2,5);
   const size_t nlat = nl + 2;
   const size_t nphi = 3*nl;
   Sphere::init(nl, nlat, nphi);

   std::vector<double> po1(Sphere::nSph(),0);
   std::vector<double> po2(Sphere::nSph(),0);
   std::vector<cplx>   ylm(Sphere::nYlm(),0);
   /* 
    * fill in values
    */
   for (size_t ip=0; ip<Sphere::nphi(); ip++) {
   for (size_t it=0; it<Sphere::nlat(); it++) {
      po1[Sphere::indx_Sph(it,ip)] = 
            1.0 
         +  pow(sin(Sphere::theta(it)),2)*pow(cos(Sphere::phi(ip)),2)
         ;
   }
   }
   Sphere::to_Ylm(po1, ylm);
   Sphere::to_Sph(ylm, po2);

   for (size_t ip=0; ip<Sphere::nphi(); ip++) {
   for (size_t it=0; it<Sphere::nlat(); it++) {
//      std::cout
//         <<std::setw(16)<<Sphere::theta(it)
//         <<std::setw(16)<<Sphere::phi(ip)
//         <<std::setw(16)<<po1[Sphere::indx_Sph(it,ip)]-po2[Sphere::indx_Sph(it,ip)]
//         <<std::endl;
      EXPECT_TRUE(
            abs(po1[Sphere::indx_Sph(it,ip)]-po2[Sphere::indx_Sph(it,ip)])
            < eps
         );
   }
   }
   Sphere::cleanup();
}
