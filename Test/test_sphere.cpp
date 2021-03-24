#include "gtest/gtest.h"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <vector>

#include "../Source/sphere.hpp"

const double eps = 1e-14;
/*==========================================================================*/
/* Testing we can go to/from spherical harmonic space, 
 * without changing the value of the array (to within truncation error). 
 */
TEST(sphere_test, to_and_from) {
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
      po1[Sphere::indx_Sph(ip,it)] = 
            1.0 
         +  pow(sin(Sphere::theta(it)),2)*pow(cos(Sphere::phi(ip)),2)
         ;
   }
   }
   Sphere::to_Ylm(po1, ylm);
   Sphere::to_Sph(ylm, po2);

   for (size_t ip=0; ip<Sphere::nphi(); ip++) {
   for (size_t it=0; it<Sphere::nlat(); it++) {
      EXPECT_TRUE(
            abs(po1[Sphere::indx_Sph(ip,it)]-po2[Sphere::indx_Sph(ip,it)])
            < eps
         );
   }
   }
   Sphere::cleanup();
}
/*==========================================================================*/
/* Testing if partial_{\phi} operator works 
 */
TEST(sphere_test, filter_is_TVD) {
   const size_t nl   = pow(2,5);
   const size_t nlat = nl + 2;
   const size_t nphi = 3*nl;
   Sphere::init(nl, nlat, nphi);

   std::vector<double> po(  Sphere::nSph(),0);
   std::vector<double> dpo1(Sphere::nSph(),0);
   std::vector<double> dpo2(Sphere::nSph(),0);
   std::vector<cplx>   ylm(Sphere::nYlm(),0);
   /* 
    * fill in values
    */
   for (size_t ip=0; ip<Sphere::nphi(); ip++) {
   for (size_t it=0; it<Sphere::nlat(); it++) {
      po1[Sphere::indx_Sph(ip,it)] = 
            pow(sin(Sphere::theta(it)),2)
            *pow(cos(Sphere::phi(ip)),2)
            *(1.0 + (0.1*rand()/RAND_MAX))
         ;
      po2[Sphere::indx_Sph(ip,it)] = po1[Sphere::indx_Sph(ip,it)];
   }
   }
   Sphere::filter(po2);
   /* 
    * compute total variations 
    */
   double tv1 = 0;
   double tv2 = 0;
   for (size_t ip=1; ip<Sphere::nphi()-1; ip++) {
   for (size_t it=1; it<Sphere::nlat()-1; it++) {
      tv1 += 
         abs(po1[Sphere::indx_Sph(ip  ,it)]
         -   po1[Sphere::indx_Sph(ip+1,it)])
         +  
         abs(po1[Sphere::indx_Sph(ip,it  )]
          -  po1[Sphere::indx_Sph(ip,it+1)]);

      tv2 += 
         abs(po2[Sphere::indx_Sph(ip  ,it)]
         -   po2[Sphere::indx_Sph(ip+1,it)])
         +  
         abs(po2[Sphere::indx_Sph(ip,it  )]
          -  po2[Sphere::indx_Sph(ip,it+1)]);
   }
   }
   EXPECT_TRUE(abs(tv2 < tv1));

   Sphere::cleanup();
}
/*==========================================================================*/
/* Testing that the spectral filter is total variation diminishing
 * when acting on ``rough'' data. 
 */
TEST(sphere_test, filter_is_TVD) {
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
      po1[Sphere::indx_Sph(ip,it)] = 
            pow(sin(Sphere::theta(it)),2)
            *pow(cos(Sphere::phi(ip)),2)
            *(1.0 + (0.1*rand()/RAND_MAX))
         ;
      po2[Sphere::indx_Sph(ip,it)] = po1[Sphere::indx_Sph(ip,it)];
   }
   }
   Sphere::filter(po2);
   /* 
    * compute total variations 
    */
   double tv1 = 0;
   double tv2 = 0;
   for (size_t ip=1; ip<Sphere::nphi()-1; ip++) {
   for (size_t it=1; it<Sphere::nlat()-1; it++) {
      tv1 += 
         abs(po1[Sphere::indx_Sph(ip  ,it)]
         -   po1[Sphere::indx_Sph(ip+1,it)])
         +  
         abs(po1[Sphere::indx_Sph(ip,it  )]
          -  po1[Sphere::indx_Sph(ip,it+1)]);

      tv2 += 
         abs(po2[Sphere::indx_Sph(ip  ,it)]
         -   po2[Sphere::indx_Sph(ip+1,it)])
         +  
         abs(po2[Sphere::indx_Sph(ip,it  )]
          -  po2[Sphere::indx_Sph(ip,it+1)]);
   }
   }
   EXPECT_TRUE(abs(tv2 < tv1));

   Sphere::cleanup();
}
