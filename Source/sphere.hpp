/*
 * Spherical harmonics and derivatives 
 */
#ifndef _SPHERE_HPP_
#define _SPHERE_HPP_

#define SHTNS_CONTIGUOUS_LONGITUDES true

#include <vector>
#include <complex>
#include <cmath>
#include <fftw3.h>
#include <shtns.h>
/*===========================================================================*/
namespace Sphere {
/*===========================================================================*/
   void init(const size_t nl, const size_t nlat, const size_t nphi);
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

   void laplace_beltrami(const std::vector<double> &v, std::vector<double> &ddv);
   void partial_phi(const std::vector<double> &v, std::vector<double> &dv);
   void filter(std::vector<double> &v);

   std::vector<double> compute_ylm(const int l_ang, const int m_ang);
/*===========================================================================*/
} /* Sphere */
/*===========================================================================*/
#endif /* _SPHERE_HPP_ */
