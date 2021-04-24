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

class Sphere
{
public:
   Sphere(
         const size_t nl, 
         const size_t nm, 
         const size_t nlat, 
         const size_t nphi
         );
   ~Sphere();

   size_t nlat() const;
   size_t nphi() const;
   size_t nSph() const;
   size_t nYlm() const;

   double theta(const size_t i_th) const; 
   double phi(  const size_t i_ph) const;
   /*
    *  These are all threadsafe. 
    */
   void to_Sph(const std::vector<cplx> &ylm, std::vector<double> &sph) const;
   void to_Ylm(const std::vector<double> &sph, std::vector<cplx> &ylm) const;

   void to_Sph(const std::vector<cplx> &ylm, std::vector<cplx> &sph) const;
   void to_Ylm(const std::vector<cplx> &sph, std::vector<cplx> &ylm) const;

   void laplace_beltrami(const std::vector<double> &v, std::vector<double> &ddv) const;
   void partial_phi(     const std::vector<double> &v, std::vector<double> &dv) const;
   void sphereX(         const std::vector<double> &v, std::vector<double> &vX) const;
   void filter(std::vector<double> &v) const;

   void raise(const std::vector<double> &v, std::vector<cplx> &rv) const;
   void lower(const std::vector<double> &v, std::vector<cplx> &lv) const;

   void power_spectrum(const std::vector<double> &v, std::vector<double> &p) const;

   std::vector<double> compute_ylm(const size_t l_ang, const size_t m_ang) const;

private:

   const size_t _mres = 1; /* angular periodicity (2*pi/mres) */
   const double _eps  = 0; /* polar optimization threshold */

   const size_t _lmax; /* highest l mode resolved */  
   const size_t _mmax; /* highest azimuthal wavenumber described */

   const size_t _nlat; /* nlat Gauss-Legendre nodes (0,pi) 
                          from sampling theorem need nlat > lmax */
   const size_t _nphi; /* nphi equally spaced nodes; from [0,2*pi/mres)
                          from sampling theorem need nphi > 2*mmax */
   
   std::vector<double> _lap; /* Laplace-Beltrami operator on the unit sphere:
                               \Delta Y_{lm} = -l(l+1) Y_{lm} */

   std::vector<double> _low_pass; /* low pass filter in 
                                     Spherical harmonic space */

   size_t _nSph; /* number of spatial points */
   size_t _nYlm; /* number of spherical harmonic coefficients */

   shtns_cfg _shtns; /* spherical harmonic struct */

   double* allocate_real(const size_t size) const;
   cplx*   allocate_cplx(const size_t size) const;
   void    free(double *tmp) const;
   void    free(cplx   *tmp) const;

public:

   inline size_t indx(const size_t i_th, const size_t i_ph) const
   {
      #if SHTNS_CONTIGUOUS_LONGITUDES
      return _nlat*i_ph + i_th;
      #else 
      return _nphi*i_th + i_ph;
      #endif
   }
};
#endif /* _SPHERE_HPP_ */
