#include <cassert>
#include "sphere.hpp"

#include <iostream>
/* 
 * for thread safe versions
 */
#define ALLOCATE_TMP \
   double *Sph_tmp = (double *) fftw_malloc( NSPAT_ALLOC(shtns_) * sizeof(double)); \
   cplx   *Ylm_tmp = (cplx *)   fftw_malloc(         shtns_->nlm * sizeof(cplx));   \
   assert(Sph_tmp!=nullptr); \
   assert(Ylm_tmp!=nullptr);
#define FREE_TMP \
   fftw_free(Sph_tmp); \
   fftw_free(Ylm_tmp);
/*==========================================================================*/
namespace Sphere {
/*==========================================================================*/
namespace {
   const size_t mres_ = 1;     /* angular periodicity (2*pi/mres) */
   const double eps_  = 1e-12; /* polar optimization threshold */

   size_t lmax_; /* highest l mode resolved */  
   size_t mmax_; /* highest azimuthal wavenumber described */
   
   size_t nYlm_; /* number of spherical harmonic coefficients */
   size_t nSph_; /* number of spatial points */


   size_t nlat_; /* nlat Gauss-Legendre nodes (0,pi) 
                    from sampling theorem need nlat > lmax */
   size_t nphi_; /* nphi equally spaced nodes; from [0,2*pi/mres)
                    from sampling theorem need nphi > 2*mmax */

   shtns_cfg shtns_; /* spherical harmonic struct */
   double *Sph_;     /* spherical space */
   cplx   *Ylm_;     /* spherical harmonic coefficients */

   std::vector<double> lap_; /* Laplace-Beltrami operator on the unit sphere:
                               \Delta Y_{lm} = -l(l+1) Y_{lm} */

   std::vector<double> low_pass_; /* low pass filter in 
                                     Spherical harmonic space */
}
/*==========================================================================*/
size_t nlat() { return nlat_; }
size_t nphi() { return nphi_; }
size_t nSph() { return nSph_; }
size_t nYlm() { return nYlm_; }
/*==========================================================================*/
double theta(const size_t i_th)
{
   return acos(shtns_->ct[i_th]);
}
double phi(const size_t i_ph)
{
   return i_ph*2.0*M_PI/((shtns_->nphi)*(shtns_->mres));
}
/*==========================================================================*/
void init(const size_t nl, const size_t nphi, const size_t nlat)
{
   lmax_ = nl;
   mmax_ = lmax_;
   nlat_ = nlat; 
   nphi_ = nphi; 

   shtns_verbose(0);     /* displays informations during initialization. */
   shtns_use_threads(1); /* enable multi-threaded transforms (if supported). */
   
   shtns_ = shtns_init( 
         shtns_type(int(sht_gauss) | int(SHT_NATIVE_LAYOUT)), 
         lmax_, mmax_, mres_, 
         nlat_, nphi_
      );
   /* 
    * Set Laplacian in spherical harmonic space: \Delta Y_{lm} = -l(l+1) Y_{lm}
    */
   lap_.resize(lmax_,0);
   for (size_t i=0; i<lmax_; i++) {
      lap_[i] = - double(i)*double(i+1.0);
   }
   /* 
    * Low pass filter in spherical harmonic space 
    */
   low_pass_.resize(lmax_,0);
   for (size_t i=0; i<lmax_; i++) {
      low_pass_[i] = exp(-36*pow(double(i)/lmax_,10));
   }
   /* 
    * Memory allocation : the use of fftw_malloc is required because we need proper 16-byte alignement.
    */
   Sph_ = (double *) fftw_malloc( NSPAT_ALLOC(shtns_) * sizeof(double));
   Ylm_ = (cplx *)   fftw_malloc(         shtns_->nlm * sizeof(cplx));
   assert(Sph_!=nullptr);
   assert(Ylm_!=nullptr);

   nSph_ = NSPAT_ALLOC(shtns_);
   nYlm_ = shtns_->nlm;
}
/*==========================================================================*/
void cleanup()
{
   assert(Sph_!=nullptr);
   assert(Ylm_!=nullptr);
   fftw_free(Sph_);
   fftw_free(Ylm_);

   shtns_destroy(shtns_);
}
/*==========================================================================*/
size_t indx_Sph(const size_t i_ph, const size_t i_th)
{
   return nlat_*i_ph + i_th;
}
/*==========================================================================*/
void to_Ylm(const std::vector<double> &sph, std::vector<cplx> &ylm)
{
   for (size_t i=0; i<nSph_; i++) {
      Sph_[i] = sph[i];
   }
   spat_to_SH(shtns_, Sph_, Ylm_); 

   for (size_t lm=0; lm<nYlm_; lm++) {
      ylm[lm] = Ylm_[lm]; 
   }
}
/*==========================================================================*/
void to_Ylm_ts(const std::vector<double> &sph, std::vector<cplx> &ylm)
{
   ALLOCATE_TMP;

   for (size_t i=0; i<nSph_; i++) {
      Sph_tmp[i] = sph[i];
   }
   spat_to_SH(shtns_, Sph_tmp, Ylm_tmp); 

   for (size_t lm=0; lm<nYlm_; lm++) {
      ylm[lm] = Ylm_tmp[lm]; 
   }
   FREE_TMP;
}
/*==========================================================================*/
void to_Sph(const std::vector<cplx> &ylm, std::vector<double> &sph)
{
   for (size_t lm=0; lm<nYlm_; lm++) {
      Ylm_[lm] = ylm[lm];
   }
   SH_to_spat(shtns_, Ylm_, Sph_); 

   for (size_t i=0; i<nSph_; i++) {
      sph[i] = Sph_[i];
   }
}
/*==========================================================================*/
void to_Sph_ts(const std::vector<cplx> &ylm, std::vector<double> &sph)
{
   ALLOCATE_TMP;

   for (size_t lm=0; lm<nYlm_; lm++) {
      Ylm_tmp[lm] = ylm[lm];
   }
   SH_to_spat(shtns_, Ylm_tmp, Sph_tmp); 

   for (size_t i=0; i<nSph_; i++) {
      sph[i] = Sph_tmp[i];
   }
   FREE_TMP;
}
/*==========================================================================*/
/* Laplace-Beltrami operator on the unit sphere:
 * \Delta Y_{lm} = -l(l+1) Y_{lm} */
/*==========================================================================*/
void laplace_beltrami(
      const std::vector<double> &v, 
      std::vector<double> &ddv)
{
   assert(v.size()  ==nlat_*nphi_);
   assert(ddv.size()==nlat_*nphi_);

   for (size_t i=0; i<nSph_; i++) {
      Sph_[i] = v[i];
   }
   spat_to_SH(shtns_, Sph_, Ylm_);

   for (size_t l=0; l<lmax_; l++) {
   for (size_t m=0; m<=l;    m++) {
      Ylm_[LM(shtns_,l,m)] *= lap_[l];
   }
   }
   SH_to_spat(shtns_, Ylm_, Sph_);

   for (size_t i=0; i<nSph_; i++) {
      ddv[i] = Sph_[i];
   }
}
/*==========================================================================*/
void laplace_beltrami_ts(
      const std::vector<double> &v, 
      std::vector<double> &ddv)
{
   ALLOCATE_TMP;

   assert(v.size()  ==nlat_*nphi_);
   assert(ddv.size()==nlat_*nphi_);

   for (size_t i=0; i<nSph_; i++) {
      Sph_tmp[i] = v[i];
   }
   spat_to_SH(shtns_, Sph_tmp, Ylm_tmp);

   for (size_t l=0; l<lmax_; l++) {
   for (size_t m=0; m<=l;    m++) {
      Ylm_tmp[LM(shtns_,l,m)] *= lap_[l];
   }
   }
   SH_to_spat(shtns_, Ylm_tmp, Sph_tmp);

   for (size_t i=0; i<nSph_; i++) {
      ddv[i] = Sph_tmp[i];
   }
   FREE_TMP;
}
/*==========================================================================*/
/* partial_{\phi} operator */
/*==========================================================================*/
void partial_phi(const std::vector<double> &v, std::vector<double> &dv)
{
   assert(v.size() ==nlat_*nphi_);
   assert(dv.size()==nlat_*nphi_);

   for (size_t i=0; i<nSph_; i++) {
      Sph_[i] = v[i];
   }
   spat_to_SH(shtns_, Sph_, Ylm_);

   const cplx img(0,1);

   for (size_t l=0; l<lmax_; l++) {
   for (size_t m=0; m<=l;    m++) {
      Ylm_[LM(shtns_,l,m)] *= img*cplx(m);
   }
   }
   SH_to_spat(shtns_, Ylm_, Sph_);

   for (size_t i=0; i<nSph_; i++) {
      dv[i] = Sph_[i];
   }
}
/*==========================================================================*/
void partial_phi_ts(const std::vector<double> &v, std::vector<double> &dv)
{
   ALLOCATE_TMP;

   assert(v.size() ==nlat_*nphi_);
   assert(dv.size()==nlat_*nphi_);

   for (size_t i=0; i<nSph_; i++) {
      Sph_tmp[i] = v[i];
   }
   spat_to_SH(shtns_, Sph_tmp, Ylm_tmp);

   const cplx img(0,1);

   for (size_t l=0; l<lmax_; l++) {
   for (size_t m=0; m<=l;    m++) {
      Ylm_tmp[LM(shtns_,l,m)] *= img*cplx(m);
   }
   }
   SH_to_spat(shtns_, Ylm_tmp, Sph_tmp);

   for (size_t i=0; i<nSph_; i++) {
      dv[i] = Sph_tmp[i];
   }
   FREE_TMP;
}
/*==========================================================================*/
/* low pass filter in spherical harmonic coefficient space
 * Note that only positive m are stored in spherical harmonic space,
 * as we only deal with real scalar fields. */
/*==========================================================================*/
void filter(std::vector<double> &v)
{
   for (size_t i=0; i<nSph_; i++) {
      Sph_[i] = v[i];
   }
   spat_to_SH(shtns_, Sph_, Ylm_);

   for (size_t l=0; l<lmax_; l++) {
   for (size_t m=0; m<l;     m++) {
      Ylm_[LM(shtns_,l,m)] *= low_pass_[l]*low_pass_[m];
   }
   }
   SH_to_spat(shtns_, Ylm_, Sph_);

   for (size_t i=0; i<nSph_; i++) {
      v[i] = Sph_[i];
   }
}
/*==========================================================================*/
void filter_ts(std::vector<double> &v)
{
   ALLOCATE_TMP;
   for (size_t i=0; i<nSph_; i++) {
      Sph_tmp[i] = v[i];
   }
   spat_to_SH(shtns_, Sph_tmp, Ylm_tmp);

   for (size_t l=0; l<lmax_; l++) {
   for (size_t m=0; m<l;     m++) {
      Ylm_tmp[LM(shtns_,l,m)] *= low_pass_[l]*low_pass_[m];
   }
   }
   SH_to_spat(shtns_, Ylm_tmp, Sph_tmp);

   for (size_t i=0; i<nSph_; i++) {
      v[i] = Sph_tmp[i];
   }
   FREE_TMP;
}
/*==========================================================================*/
/* returns values for spherical harmonic in real space Y_{l_ang,m_ang} */
/*==========================================================================*/
std::vector<double> compute_ylm(const int l_ang, const int m_ang)
{
   std::vector<cplx>   in(nYlm_, 0.0);
   std::vector<double> out(nSph_,0.0);

   in[LM(shtns_,l_ang,m_ang)] = 1;

   to_Sph(in, out);

   return out;
}
/*==========================================================================*/
std::vector<double> compute_ylm_ts(const int l_ang, const int m_ang)
{
   std::vector<cplx>   in(nYlm_, 0.0);
   std::vector<double> out(nSph_,0.0);

   in[LM(shtns_,l_ang,m_ang)] = 1;

   to_Sph_ts(in, out);

   return out;
}
/*==========================================================================*/
} /* Sphere */

#undef ALLOCATE_TMP
#undef FREE_TMP
