/*
 * Spherical harmonics and derivatives 
 */
#ifndef _SPHERE_HPP_
#define _SPHERE_HPP_

#include <stdbool.h>
/*
 * SHTNS_CONTIGUOUS_LONGITUDES true:
 * means theta is the ``fastest'' index,
 * i.e. the indexing is 
 * _nlat*ip + it.
 * Setting it to false means we index as
 * _nphi*it + ip
 *
 * For the library shtns, usually the fastest transforms use
 * contiguous longitudinal indexing.
 *
 * If you use the function indx, you ALWAYS access 
 * the angular coordinates as
 * indx(it,ip),
 * which matches the usual notation (theta,phi) in spherical
 * polar coordinates. 
 */
#define SHTNS_CONTIGUOUS_LONGITUDES true

#include <vector>
#include <complex>
#include <cmath>
#include <fftw3.h>
#include <shtns.h>
/*===========================================================================*/
namespace Sphere {
/*===========================================================================*/
   void init(
         const size_t nl, 
         const size_t nm, 
         const size_t nlat, 
         const size_t nphi
         );
   void cleanup();

   size_t indx(const size_t i_th, const size_t i_ph);

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

   void to_Sph(const std::vector<cplx> &ylm, std::vector<cplx> &sph);
   void to_Ylm(const std::vector<cplx> &sph, std::vector<cplx> &ylm);

   void laplace_beltrami(const std::vector<double> &v, std::vector<double> &ddv);
   void partial_phi(const std::vector<double> &v, std::vector<double> &dv);
   void raise(const std::vector<double> &v, std::vector<cplx> &rv);
   void lower(const std::vector<double> &v, std::vector<cplx> &lv);

   /* sphereX
    * =
    * (\partial_{\theta} f)^2 + (1/sin^2(theta)) (\partial_{\phi}f)^2 
    * =
    * - (L_+f)(L_-f) + (\partial_{\phi}f)^2 */
   void sphereX(const std::vector<double> &v, std::vector<double> &vX);

   void filter(std::vector<double> &v);

   std::vector<double> compute_ylm(const int l_ang, const int m_ang);
/*===========================================================================*/
} /* Sphere */
/*===========================================================================*/
#endif /* _SPHERE_HPP_ */
