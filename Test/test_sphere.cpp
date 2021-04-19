#include "gtest/gtest.h"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <vector>

#include "../Source/sphere.hpp"

/*==========================================================================*/
/* Testing that the we can make spherical harmonics with all
 * the correct signs and normalizations.
 * Namely: we use the Condon-Shortly phase convention,
 * and normalize the spherical harmonics so they are orthonormal on
 * the unit sphere. 
 */
TEST(test_sphere, make_Ylm) {
   const size_t nl   = pow(2,5);
   const size_t nm   = pow(2,4);
   const size_t nlat = nl + 2;
   const size_t nphi = 3*nl;
   Sphere::init(nl, nm, nlat, nphi);

   std::vector<double> y00 = Sphere::compute_ylm(0, 0);
   std::vector<double> y10 = Sphere::compute_ylm(1, 0);
   std::vector<double> y11 = Sphere::compute_ylm(1, 1);

   for (size_t ip=0; ip<Sphere::nphi(); ip++) {
   for (size_t it=0; it<Sphere::nlat(); it++) {
      const double theta = Sphere::theta(it);
      const double phi   = Sphere::phi(ip);
      EXPECT_LT(fabs(
            pow(1.0/(4.0*M_PI),0.5)
         -  y00[Sphere::indx(it,ip)]
         ),
            1e-15
         );
      EXPECT_LT(fabs(
            pow(3.0/(4.0*M_PI),0.5)*cos(theta) 
         - y10[Sphere::indx(it,ip)]
         ),
            1e-15
         );
      EXPECT_LT(fabs(
         -  pow(3.0/(2.0*M_PI),0.5)*sin(theta)*cos(phi) 
         -  y11[Sphere::indx(it,ip)]
         ), 
            1e-15
         );
   }
   }

   Sphere::cleanup();
}
/*==========================================================================*/
/* Testing we can go to/from spherical harmonic space,
 * (for both real and complex arrays) 
 * without changing the value of the array (to within truncation error). 
 */
TEST(test_sphere, to_and_from) {
   const size_t nl   = pow(2,5);
   const size_t nm   = pow(2,4);
   const size_t nlat = nl + 2;
   const size_t nphi = 3*nl;
   Sphere::init(nl, nm, nlat, nphi);

   std::vector<double> r_1(Sphere::nSph(),0);
   std::vector<double> r_2(Sphere::nSph(),0);
   std::vector<cplx>   c_1(Sphere::nSph(),0);
   std::vector<cplx>   c_2(Sphere::nSph(),0);
   std::vector<cplx>   r_ylm(Sphere::nYlm(),0);
   std::vector<cplx>   c_ylm(pow(nl+1,2),0);
   /* 
    * fill in values
    */
   for (size_t ip=0; ip<Sphere::nphi(); ip++) {
   for (size_t it=0; it<Sphere::nlat(); it++) {
      r_1[Sphere::indx(it,ip)] = 
            1.0 
         +  pow(sin(Sphere::theta(it)),2)*pow(cos(Sphere::phi(ip)),2)
         ;
      c_1[Sphere::indx(it,ip)] = 
            1.0 
         +  pow(sin(Sphere::theta(it)),2)*pow(cos(Sphere::phi(ip)),2)
         ;
   }
   }
   /* 
    * Real and complex transforms 
    */
   Sphere::to_Ylm(r_1, r_ylm);
   Sphere::to_Sph(r_ylm, r_2);

   Sphere::to_Ylm(c_1, c_ylm);
   Sphere::to_Sph(c_ylm, c_2);

   for (size_t ip=0; ip<Sphere::nphi(); ip++) {
   for (size_t it=0; it<Sphere::nlat(); it++) {
      EXPECT_LT(
            fabs(r_1[Sphere::indx(it,ip)]-r_2[Sphere::indx(it,ip)]),
            1e-14
         );
   }
   }
   for (size_t ip=0; ip<Sphere::nphi(); ip++) {
   for (size_t it=0; it<Sphere::nlat(); it++) {
      EXPECT_LT(
            0.5*fabs(c_1[Sphere::indx(it,ip)].real()-c_2[Sphere::indx(it,ip)].real())
         +  0.5*fabs(c_1[Sphere::indx(it,ip)].imag()-c_2[Sphere::indx(it,ip)].imag()),
            1e-14
         );
   }
   }
   Sphere::cleanup();
}
/*==========================================================================*/
/* Testing partial_{\phi} operator acts correctly
 */
TEST(test_sphere, partial_phi) {
   const size_t nl   = pow(2,5);
   const size_t nm   = pow(2,4);
   const size_t nlat = nl + 2;
   const size_t nphi = 3*nl;
   Sphere::init(nl, nm, nlat, nphi);

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
/* Testing raising and lowering operators 
 */
TEST(test_sphere, raise_and_lower) {
   const size_t nl   = 20;
   const size_t nm   = 14;
   const size_t nlat = 40;
   const size_t nphi = 32;
   Sphere::init(nl, nm, nlat, nphi);

   std::vector<double> v(Sphere::nSph(),0);
   std::vector<cplx>  r1(Sphere::nSph(),0);
   std::vector<cplx>  r2(Sphere::nSph(),0);
   std::vector<cplx>  l1(Sphere::nSph(),0);
   std::vector<cplx>  l2(Sphere::nSph(),0);
   /* 
    * fill in values
    */
   for (size_t ip=0; ip<Sphere::nphi(); ip++) {
   for (size_t it=0; it<Sphere::nlat(); it++) {
      /* Y_{l=1,m=0} */
      v[Sphere::indx(it,ip)] = 
            0.5*pow(3/M_PI,0.5)*cos(Sphere::theta(it))
         ;
      /* sqrt(2) * Y_{l=1,m=1} */
      r2[Sphere::indx(it,ip)] = 
         pow(1*(1+1),0.5)*(
            -  0.5*pow(3/(2.0*M_PI),0.5)*sin(Sphere::theta(it))*(
                  cos(Sphere::phi(ip))
               +  cplx(0,1)*sin(Sphere::phi(ip))
               )
            )
         ;
      /* sqrt(2) * Y_{l=1,m=-1} */
      l2[Sphere::indx(it,ip)] = 
         pow(1*(1+1),0.5)*(
               0.5*pow(3/(2.0*M_PI),0.5)*sin(Sphere::theta(it))*(
                  cos(Sphere::phi(ip))
               -  cplx(0,1)*sin(Sphere::phi(ip))
               )
            )
         ;
   }
   }
   /* 
    * Default transform 
    */
   Sphere::raise(v, r1);
   Sphere::lower(v, l1);

   for (size_t ip=0; ip<Sphere::nphi(); ip++) {
   for (size_t it=0; it<Sphere::nlat(); it++) {
      EXPECT_LT(
            0.5*fabs(r1[Sphere::indx(it,ip)].real()-r2[Sphere::indx(it,ip)].real())
         +  0.5*fabs(r1[Sphere::indx(it,ip)].imag()-r2[Sphere::indx(it,ip)].imag()),
            5e-12
         );
   }
   }
   Sphere::cleanup();
}
/*==========================================================================*/
/* Testing sphereX operator acts correctly 
 */
TEST(test_sphere, sphereX) {
   const size_t nl   = 24;
   const size_t nm   = 12;
   const size_t nlat = 48;
   const size_t nphi = 28;
   Sphere::init(nl, nm, nlat, nphi);

   std::vector<double> v(  Sphere::nSph(),0);
   std::vector<double> vX1(Sphere::nSph(),0);
   std::vector<double> vX2(Sphere::nSph(),0);
   /* 
    * fill in values
    */
   for (size_t ip=0; ip<Sphere::nphi(); ip++) {
   for (size_t it=0; it<Sphere::nlat(); it++) {
      const double theta = Sphere::theta(it);
      const double phi   = Sphere::phi(ip);

    v[Sphere::indx(it,ip)] =
            1.0 
         +  pow(sin(theta),2)*pow(cos(phi),2)
         ;
      vX2[Sphere::indx(it,ip)] =
            4.0*pow(sin(theta)*cos(phi),2)*(
                  pow(cos(theta)*cos(phi),2)
               +  pow(sin(phi),2)
               )
         ;
   }
   }
   /* 
    * Default transform 
    */
   Sphere::sphereX(v, vX1);

   for (size_t ip=0; ip<Sphere::nphi(); ip++) {
   for (size_t it=0; it<Sphere::nlat(); it++) {
      EXPECT_LT(
            fabs(vX1[Sphere::indx(it,ip)]-vX2[Sphere::indx(it,ip)]),
            5e-12
         );
   }
   }
   Sphere::cleanup();
}
/*==========================================================================*/
/* Testing Spherical Laplace-Beltrami operator acts correctly
 */
TEST(test_sphere, laplace_beltrami) {
   const size_t nl   = 32;
   const size_t nm   = 24;
   const size_t nlat = 34;
   const size_t nphi = 60;
   Sphere::init(nl, nm, nlat, nphi);

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
TEST(test_sphere, filter_is_TVD) {
   const size_t nl   = pow(2,5);
   const size_t nm   = pow(2,4);
   const size_t nlat = nl + 2;
   const size_t nphi = 3*nl;
   Sphere::init(nl, nm, nlat, nphi);

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
   for (size_t ip=0; ip<Sphere::nphi()-1; ip++) {
   for (size_t it=0; it<Sphere::nlat()-1; it++) {
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
