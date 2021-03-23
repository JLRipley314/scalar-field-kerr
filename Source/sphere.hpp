#ifndef _SPHERE_HPP_
#define _SPHERE_HPP_
/*===========================================================================*/
/* Spherical harmonics and derivatives */
/*===========================================================================*/
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

   void to_Sph(const std::vector<cplx> &ylm, std::vector<double> &sph);
   void to_Ylm(const std::vector<double> &sph, std::vector<cplx> &ylm);

   size_t indx_Sph(const size_t i_th, const size_t i_ph);

   size_t nlat();
   size_t nphi();
   size_t nSph();
   size_t nYlm();

   double theta(const size_t i_th);  
   double phi(  const size_t i_ph);

   std::vector<double> compute_ylm(const int l_ang, const int m_ang);
/*===========================================================================*/
} /* Sphere */
/*===========================================================================*/
#endif /* _SPHERE_HPP_ */
