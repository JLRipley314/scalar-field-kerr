#include "gtest/gtest.h"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <vector>

#include "../Source/sphere.hpp"

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
      po1[Sphere::indx(it,ip)] = 
            1.0 
         +  pow(sin(Sphere::theta(it)),2)*pow(cos(Sphere::phi(ip)),2)
         ;
   }
   }
   /* 
    * Default transform 
    */
   Sphere::to_Ylm(po1, ylm);
   Sphere::to_Sph(ylm, po2);

   for (size_t ip=0; ip<Sphere::nphi(); ip++) {
   for (size_t it=0; it<Sphere::nlat(); it++) {
      EXPECT_LT(
            fabs(po1[Sphere::indx(it,ip)]-po2[Sphere::indx(it,ip)]),
            1e-14
         );
   }
   }
   Sphere::cleanup();
}
/*==========================================================================*/
/* Testing partial_{\phi} operator acts correctly
 */
TEST(sphere_test, partial_phi) {
   const size_t nl   = pow(2,5);
   const size_t nlat = nl + 2;
   const size_t nphi = 3*nl;
   Sphere::init(nl, nlat, nphi);

   std::vector<double> v(  Sphere::nSph(),0);
   std::vector<double> dv1(Sphere::nSph(),0);
   std::vector<double> dv2(Sphere::nSph(),0);
   /* 
    * fill in values
    */
   for (size_t ip=0; ip<Sphere::nphi(); ip++) {
   for (size_t it=0; it<Sphere::nlat(); it++) {
      v[Sphere::indx(it,ip)] = 
            1.0 
         +  pow(sin(Sphere::theta(it)),2)*pow(cos(Sphere::phi(ip)),2)
         ;
      dv2[Sphere::indx(it,ip)] = 
            pow(sin(Sphere::theta(it)),2)*(
                  -  2.0*sin(Sphere::phi(ip))*cos(Sphere::phi(ip))
                  )
         ;
   }
   }
   /* 
    * Default transform 
    */
   Sphere::partial_phi(v, dv1);

   for (size_t ip=0; ip<Sphere::nphi(); ip++) {
   for (size_t it=0; it<Sphere::nlat(); it++) {
      EXPECT_LT(
            fabs(dv1[Sphere::indx(it,ip)]-dv2[Sphere::indx(it,ip)]),
            5e-14
         );
   }
   }
   Sphere::cleanup();
}
/*==========================================================================*/
/* Testing Spherical Laplace-Beltrami operator acts correctly
 */
TEST(sphere_test, laplace_beltrami) {
   const size_t nl   = pow(2,5);
   const size_t nlat = nl + 2;
   const size_t nphi = 3*nl;
   Sphere::init(nl, nlat, nphi);

   std::vector<double> v(   Sphere::nSph(),0);
   std::vector<double> ddv1(Sphere::nSph(),0);
   std::vector<double> ddv2(Sphere::nSph(),0);
   /* 
    * fill in values
    */
   for (size_t ip=0; ip<Sphere::nphi(); ip++) {
   for (size_t it=0; it<Sphere::nlat(); it++) {
      v[Sphere::indx(it,ip)] = 
            1.0 
         +  pow(sin(Sphere::theta(it)),2)*pow(cos(Sphere::phi(ip)),2)
         ;
      ddv2[Sphere::indx(it,ip)] = 
         0.5*(
               1.0
            +  3.0*cos(2.0*Sphere::theta(it))
            -  6.0*cos(2.0*Sphere::phi(ip))*pow(sin(Sphere::theta(it)),2)
            ) 
         ;
   }
   }
   /* 
    * Default transform 
    */
   Sphere::laplace_beltrami(v, ddv1);

   for (size_t ip=0; ip<Sphere::nphi(); ip++) {
   for (size_t it=0; it<Sphere::nlat(); it++) {
      EXPECT_LT(
            fabs(ddv1[Sphere::indx(it,ip)]-ddv2[Sphere::indx(it,ip)]),
            5e-12
         );
   }
   }
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
      po1[Sphere::indx(it,ip)] = 
            pow(sin(Sphere::theta(it)),2)
            *pow(cos(Sphere::phi(ip)),2)
            *(1.0 + (0.1*rand()/RAND_MAX))
         ;
      po2[Sphere::indx(it,ip)] = po1[Sphere::indx(it,ip)];
   }
   }
   /* 
    * Default filter 
    */
   Sphere::filter(po2);
   /* 
    * compute total variations 
    */
   double tv1 = 0;
   double tv2 = 0;
   for (size_t ip=1; ip<Sphere::nphi()-1; ip++) {
   for (size_t it=1; it<Sphere::nlat()-1; it++) {
      tv1 += 
         fabs(po1[Sphere::indx(it,ip)]
         -    po1[Sphere::indx(it,ip+1)])
         +  
         fabs(po1[Sphere::indx(it,  ip)]
          -   po1[Sphere::indx(it+1,ip)]);

      tv2 += 
         fabs(po2[Sphere::indx(it,ip  )]
         -    po2[Sphere::indx(it,ip+1)])
         +  
         fabs(po2[Sphere::indx(it,  ip)]
          -   po2[Sphere::indx(it+1,ip)]);
   }
   }
   EXPECT_LT(tv2, tv1);

   Sphere::cleanup();
}
