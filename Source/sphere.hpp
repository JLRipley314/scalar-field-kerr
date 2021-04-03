/*
 * Spherical harmonics and derivatives 
 */
#ifndef _SPHERE_HPP_
#define _SPHERE_HPP_

#include <vector>
#include <complex>
#include <cmath>
#include <fftw3.h>
#include <shtns.h>
/*===========================================================================*/
namespace Sphere {
/*===========================================================================*/
   void init(const size_t nl, const size_t nphi, const size_t nlat);
   void cleanup();

   /*
    * NOTE the indexing! we are using "SHT_NATIVE_LAYOUT", so in
    * fact theta varies the fastest; i.e. we index as (phi, theta) 
    */
   size_t indx_Sph(const size_t i_ph, const size_t i_th);

   size_t nlat();
   size_t nphi();
   size_t nSph();
   size_t nYlm();

   double theta(const size_t i_th);  
   double phi(  const size_t i_ph);

   /*
    *  These are all threadsafe. 
    */
   void to_Sph(const std::vector<cplx> &ylm, std::vector<double> &sph);
   void to_Ylm(const std::vector<double> &sph, std::vector<cplx> &ylm);

   void laplace_beltrami(const std::vector<double> &v, std::vector<double> &ddv);
   void partial_phi(const std::vector<double> &v, std::vector<double> &dv);
   void filter(std::vector<double> &v);

   std::vector<double> compute_ylm(const int l_ang, const int m_ang);
/*===========================================================================*/
} /* Sphere */
/*===========================================================================*/
#endif /* _SPHERE_HPP_ */
