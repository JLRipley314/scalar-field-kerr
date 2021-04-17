#include <cassert>
#include "sphere.hpp"

#include <iostream>
/*==========================================================================*/
namespace Sphere {
/*==========================================================================*/
namespace {
   const size_t _mres = 1; /* angular periodicity (2*pi/mres) */
   const double _eps  = 0; /* polar optimization threshold */

   size_t _lmax; /* highest l mode resolved */  
   size_t _mmax; /* highest azimuthal wavenumber described */

   size_t _nYlm; /* number of spherical harmonic coefficients */
   size_t _nSph; /* number of spatial points */

   size_t _nlat; /* nlat Gauss-Legendre nodes (0,pi) 
                    from sampling theorem need nlat > lmax */
   size_t _nphi; /* nphi equally spaced nodes; from [0,2*pi/mres)
                    from sampling theorem need nphi > 2*mmax */
   
   shtns_cfg _shtns; /* spherical harmonic struct */

   std::vector<double> _lap; /* Laplace-Beltrami operator on the unit sphere:
                               \Delta Y_{lm} = -l(l+1) Y_{lm} */

   std::vector<double> _low_pass; /* low pass filter in 
                                     Spherical harmonic space */
}
/*==========================================================================*/
size_t nlat() { return _nlat; }
size_t nphi() { return _nphi; }
size_t nSph() { return _nSph; }
size_t nYlm() { return _nYlm; }
/*==========================================================================*/
double theta(const size_t i_th)
{
   return acos(_shtns->ct[i_th]);
}
double phi(const size_t i_ph)
{
   return i_ph*2.0*M_PI/((_shtns->nphi)*(_shtns->mres));
}
/*==========================================================================*/
void init(
      const size_t nl, const size_t nm, 
      const size_t nlat, const size_t nphi)
{
   _lmax = nl;
   _mmax = nm;
   _nlat = nlat; 
   _nphi = nphi; 

   shtns_verbose(0);     /* displays informations during initialization. */
   shtns_use_threads(1); /* enable multi-threaded transforms (if supported). */

   _shtns = shtns_create(_lmax, _mmax, _mres, sht_orthonormal);
#if SHTNS_CONTIGUOUS_LONGITUDES
   const int num_spat_allocate =
      shtns_set_grid(
            _shtns, 
            shtns_type(int(sht_gauss) | int(SHT_NATIVE_LAYOUT)), 
            _eps, 
            _nlat, _nphi
         );
#else
   const int num_spat_allocate =
      shtns_set_grid(
            _shtns, 
            shtns_type(int(sht_gauss) | int(SHT_PHI_CONTIGUOUS)), 
            _eps, 
            _nlat, _nphi
         );
#endif
   /* 
    * check set shtns and grid 
    */
   assert(NSPAT_ALLOC(_shtns)==size_t(num_spat_allocate)); 
   /* 
    * Set Laplacian in spherical harmonic space: \Delta Y_{lm} = -l(l+1) Y_{lm}
    */
   _lap.resize(_lmax,0);
   for (size_t i=0; i<_lmax; i++) {
      _lap[i] = - double(i)*double(i+1.0);
   }
   /* 
    * Low pass filter in spherical harmonic space 
    */
   _low_pass.resize(_lmax,0);
   for (size_t i=0; i<_lmax; i++) {
      _low_pass[i] = exp(-40.0*pow(double(i)/_lmax,16));
   }
   _nSph = NSPAT_ALLOC(_shtns);
   _nYlm = _shtns->nlm;
}
/*==========================================================================*/
void cleanup()
{
   shtns_destroy(_shtns);
}
/*==========================================================================*/
/* Need to use fftw_malloc instead of vectors
 * when using shtns. */
/*==========================================================================*/
double* allocate_real(const size_t size)
{
   double *tmp = (double *) fftw_malloc(size * sizeof(double));
   assert(tmp!=nullptr);
   return tmp;
}
/*==========================================================================*/
cplx* allocate_cplx(const size_t size)
{
   cplx *tmp = (cplx *) fftw_malloc(size * sizeof(cplx));
   assert(tmp!=nullptr);
   return tmp;
}
/*==========================================================================*/
void free(double *tmp)
{
   fftw_free(tmp);
   tmp = nullptr;
}
/*==========================================================================*/
void free(cplx *tmp)
{
   fftw_free(tmp);
   tmp = nullptr;
}
/*==========================================================================*/
size_t indx(const size_t i_th, const size_t i_ph)
{
#if SHTNS_CONTIGUOUS_LONGITUDES
   return _nlat*i_ph + i_th;
#else 
   return _nphi*i_th + i_ph;
#endif
}
/*==========================================================================*/
void to_Ylm(const std::vector<double> &sph, std::vector<cplx> &ylm)
{
   assert(sph.size()==NSPAT_ALLOC(_shtns));
   assert(ylm.size()==_shtns->nlm);

   double *Sph_tmp = allocate_real(NSPAT_ALLOC(_shtns));
   cplx   *Ylm_tmp = allocate_cplx(_shtns->nlm);

   for (size_t i=0; i<_nSph; i++) {
      Sph_tmp[i] = sph[i];
   }
   spat_to_SH(_shtns, Sph_tmp, Ylm_tmp); 

   for (size_t lm=0; lm<_nYlm; lm++) {
      ylm[lm] = Ylm_tmp[lm]; 
   }
   free(Sph_tmp);
   free(Ylm_tmp);
}
/*==========================================================================*/
void to_Sph(const std::vector<cplx> &ylm, std::vector<double> &sph)
{
   assert(sph.size()==NSPAT_ALLOC(_shtns));
   assert(ylm.size()==_shtns->nlm);

   double *Sph_tmp = allocate_real(NSPAT_ALLOC(_shtns));
   cplx   *Ylm_tmp = allocate_cplx(_shtns->nlm);

   for (size_t lm=0; lm<_nYlm; lm++) {
      Ylm_tmp[lm] = ylm[lm];
   }
   SH_to_spat(_shtns, Ylm_tmp, Sph_tmp); 

   for (size_t i=0; i<_nSph; i++) {
      sph[i] = Sph_tmp[i];
   }
   free(Sph_tmp);
   free(Ylm_tmp);
}
/*==========================================================================*/
void to_Ylm(const std::vector<cplx> &sph, std::vector<cplx> &ylm)
{
   assert(sph.size()==NSPAT_ALLOC(_shtns));
   assert(ylm.size()==pow(_shtns->nlm+1,2));

   cplx *Sph_tmp = allocate_cplx(NSPAT_ALLOC(_shtns));
   cplx *Ylm_tmp = allocate_cplx(pow(_shtns->nlm+1,2));

   for (size_t i=0; i<_nSph; i++) {
      Sph_tmp[i] = sph[i];
   }
   spat_cplx_to_SH(_shtns, Sph_tmp, Ylm_tmp); 

   for (size_t lm=0; lm<_nYlm; lm++) {
      ylm[lm] = Ylm_tmp[lm]; 
   }
   free(Sph_tmp);
   free(Ylm_tmp);
}
/*==========================================================================*/
void to_Sph(const std::vector<cplx> &ylm, std::vector<cplx> &sph)
{
   assert(sph.size()==NSPAT_ALLOC(_shtns));
   assert(ylm.size()==pow(_shtns->nlm+1,2));

   cplx *Sph_tmp = allocate_cplx(NSPAT_ALLOC(_shtns));
   cplx *Ylm_tmp = allocate_cplx(pow(_shtns->nlm+1,2));

   for (size_t lm=0; lm<_nYlm; lm++) {
      Ylm_tmp[lm] = ylm[lm];
   }
   SH_to_spat_cplx(_shtns, Ylm_tmp, Sph_tmp); 

   for (size_t i=0; i<_nSph; i++) {
      sph[i] = Sph_tmp[i];
   }
   free(Sph_tmp);
   free(Ylm_tmp);
}
/*==========================================================================*/
void laplace_beltrami(
      const std::vector<double> &v, 
      std::vector<double> &ddv)
{
   double *Sph_tmp = allocate_real(NSPAT_ALLOC(_shtns));
   cplx   *Ylm_tmp = allocate_cplx(_shtns->nlm);

   assert(v.size()  ==_nlat*_nphi);
   assert(ddv.size()==_nlat*_nphi);

   for (size_t i=0; i<_nSph; i++) {
      Sph_tmp[i] = v[i];
   }
   spat_to_SH(_shtns, Sph_tmp, Ylm_tmp);

   for (size_t l=0; l<_lmax;               l++) {
   for (size_t m=0; m<std::min(l+1,_mmax); m++) {
      Ylm_tmp[LM(_shtns,l,m)] *= _lap[l];
   }
   }
   SH_to_spat(_shtns, Ylm_tmp, Sph_tmp);

   for (size_t i=0; i<_nSph; i++) {
      ddv[i] = Sph_tmp[i];
   }
   free(Sph_tmp);
   free(Ylm_tmp);
}
/*==========================================================================*/
/* partial_{\phi} operator */
/*==========================================================================*/
void partial_phi(const std::vector<double> &v, std::vector<double> &dv)
{
   double *Sph_tmp = allocate_real(NSPAT_ALLOC(_shtns));
   cplx   *Ylm_tmp = allocate_cplx(_shtns->nlm);

   assert(v.size() ==_nlat*_nphi);
   assert(dv.size()==_nlat*_nphi);

   for (size_t i=0; i<_nSph; i++) {
      Sph_tmp[i] = v[i];
   }
   spat_to_SH(_shtns, Sph_tmp, Ylm_tmp);

   const cplx img(0,1);

   for (size_t l=0; l<_lmax;               l++) {
   for (size_t m=0; m<std::min(l+1,_mmax); m++) {
      Ylm_tmp[LM(_shtns,l,m)] *= img*cplx(m);
   }
   }
   SH_to_spat(_shtns, Ylm_tmp, Sph_tmp);

   for (size_t i=0; i<_nSph; i++) {
      dv[i] = Sph_tmp[i];
   }
   free(Sph_tmp);
   free(Ylm_tmp);
}
/*==========================================================================*/
/* low pass filter in spherical harmonic coefficient space
 * Note that only positive m are stored in spherical harmonic space,
 * as we only deal with real scalar fields. */
/*==========================================================================*/
void filter(std::vector<double> &v)
{
   double *Sph_tmp = allocate_real(NSPAT_ALLOC(_shtns));
   cplx   *Ylm_tmp = allocate_cplx(_shtns->nlm);

   for (size_t i=0; i<_nSph; i++) {
      Sph_tmp[i] = v[i];
   }
   spat_to_SH(_shtns, Sph_tmp, Ylm_tmp);

   for (size_t l=0; l<_lmax;               l++) {
   for (size_t m=0; m<std::min(l+1,_mmax); m++) {
      Ylm_tmp[LM(_shtns,l,m)] *= _low_pass[l]*_low_pass[m];
   }
   }
   SH_to_spat(_shtns, Ylm_tmp, Sph_tmp);

   for (size_t i=0; i<_nSph; i++) {
      v[i] = Sph_tmp[i];
   }
   free(Sph_tmp);
   free(Ylm_tmp);
}
/*==========================================================================*/
/* returns values for spherical harmonic in real space Y_{l_ang,m_ang} */
/*==========================================================================*/
std::vector<double> compute_ylm(const int l_ang, const int m_ang)
{
   std::vector<cplx>   in( _nYlm,0.0);
   std::vector<double> out(_nSph,0.0);

   in[LM(_shtns,l_ang,m_ang)] = 1;

   to_Sph(in, out);

   return out;
}
/*==========================================================================*/
} /* Sphere */
