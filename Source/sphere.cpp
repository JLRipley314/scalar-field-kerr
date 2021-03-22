#include <cassert>
#include "sphere.hpp"

/*==========================================================================*/
namespace Sphere {
/*==========================================================================*/
namespace {
   shtns_cfg shtns_; /* spherical harmonic struct */

   size_t lmax_; /* highest l mode resolved */  
   size_t mmax_; /* highest azimuthal wavenumber described */
   
   size_t nYlm_; /* number of spherical harmonic coefficients */
   size_t nSph_; /* number of spatial points */

   const size_t mres_ = 1;  /* angular periodicity (2*pi/mres) */

   const double eps_ = 1e-12; /* polar optimization threshold */

   size_t nlat_; /* nlat Gauss-Legendre nodes (0,pi) 
                    from sampling theorem need nlat > lmax */
   size_t nphi_; /* nphi equally spaced nodes; from [0,2*pi/mres)
                    from sampling theorem need nphi > 2*mmax */

   double *Sph_; /* spherical space */
   cplx   *Ylm_; /* spherical harmonic coefficients */

   std::vector<double> lap; /* Laplace-Beltrami operator on the unit sphere:
                               \Delta Y_{lm} = -l(l+1) Y_{lm} */
}
/*==========================================================================*/
void init(const size_t nl)
{
   lmax_ = nl;
   mmax_ = lmax_;
   nlat_ = lmax_ + 2; 
   if (nl%2==0) { nphi_ = 3*mmax_; } 
   else {         nphi_ = 2*mmax_; }

   shtns_verbose(1);     /* displays informations during initialization. */
   shtns_use_threads(0); /* enable multi-threaded transforms (if supported). */
   
   shtns_ = shtns_init( 
         shtns_type(int(sht_gauss) | int(SHT_NATIVE_LAYOUT)), 
         lmax_, mmax_, mres_, 
         nlat_, nphi_
      );
   /* 
    * Set Laplacian in spherical harmonic space: \Delta Y_{lm} = -l(l+1) Y_{lm}
    */
   lap.resize(lmax_,0);
   for (size_t i=0; i<lmax_; i++) {
      lap[i] = - i*(i+1);
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
/* spatial index */
/*==========================================================================*/
inline size_t I_SPH(const size_t i_ph, const size_t i_th) {
   return i_ph*nlat_ + i_th;
}
/*==========================================================================*/
size_t indx_Sph(const size_t i_ph, const size_t i_th) {
   return i_ph*nlat_ + i_th;
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
/* Laplace-Beltrami operator on the unit sphere:
 * \Delta Y_{lm} = -l(l+1) Y_{lm} */
/*==========================================================================*/
void laplace_beltrami(
      const std::vector<double> v, 
      std::vector<double> ddv)
{
   assert(v.size()  ==nlat_*nphi_);
   assert(ddv.size()==nlat_*nphi_);

   for (size_t i=0; i<nSph_; i++) {
      Sph_[i] = v[i];
   }
   spat_to_SH(shtns_, Sph_, Ylm_);

   for (size_t lm=0; lm<nYlm_; lm++) {
      const size_t l = shtns_->li[lm];
      Ylm_[lm] *= lap[l];
   }
   SH_to_spat(shtns_, Ylm_, Sph_);

   for (size_t i=0; i<nSph_; i++) {
      ddv[i] = Sph_[i];
   }
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
/* returns values for spherical harmonic in real space Y_{l_ang,m_ang} */
/*==========================================================================*/
std::vector<double> compute_ylm(const int l_ang, const int m_ang)
{
   std::vector<cplx>   in(nYlm_,{0.0,0.0});
   std::vector<double> out(nSph_,0.0);

   in[LM(shtns_,l_ang,m_ang)] = 1;

   to_Sph(in, out);

   return out;
}
/*==========================================================================*/
} /* Sphere */
