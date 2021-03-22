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
   void init(const size_t nl);
   void cleanup();

   void to_Sph(const std::vector<cplx> &ylm, std::vector<double> &sph);
   void to_Ylm(const std::vector<double> &sph, std::vector<cplx> &ylm);

   size_t indx_Sph(const size_t i_ph, const size_t i_th);

   double nlat();
   double nphi();
   double nSph();
   double nYlm();

   double theta(const size_t i_th);  
   double phi(  const size_t i_ph);
/*===========================================================================*/
} /* Sphere */
/*===========================================================================*/
#endif /* _SPHERE_HPP_ */
