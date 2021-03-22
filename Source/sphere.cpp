#include <cassert>
#include "sphere.hpp"

/*==========================================================================*/
namespace Sphere {
/*==========================================================================*/
namespace {
   shtns_cfg shtns_; /* spherical harmonic struct */

   size_t lmax_; /* highest l mode resolved */  
   size_t mmax_; /* highest azimuthal wavenumber described */

   const size_t mres_ = 1;  /* angular periodicity (2*pi/mres) */

   const double eps_ = 1e-12; /* polar optimization threshold */

   size_t nlat_; /* nlat gauss nodes (0,pi) 
                    from sampling theorem need nlat > 2*lmax */
   size_t nphi_; /* nphi equally spaced nodes; from [0,2*pi/mres)
                    from sampling theorem need nphi > 2*mmax */

   double *Sph_; /* spherical space */
   cplx   *Ylm_; /* spherical harmonic coefficients */

   std::vector<double> lap;
}
/*==========================================================================*/
void init(const size_t nl)
{
   lmax_ = nl;
   mmax_ = lmax_;
   if (nl%2==0) {
      nlat_ = 3*lmax_; nphi_ = 3*mmax_;
   } else {
      nlat_ = 2*lmax_; nphi_ = 2*mmax_;
   }
   shtns_verbose(1);     /* displays informations during initialization. */
   shtns_use_threads(0); /* enable multi-threaded transforms (if supported). */
   
   shtns_ = shtns_init( 
         shtns_type(int(sht_gauss) | int(SHT_NATIVE_LAYOUT)), 
         lmax_, mmax_, mres_, 
         nlat_, nphi_
      );
   /* 
    * Set Laplacian in spherical harmonic space 
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
}
/*==========================================================================*/
/* spatial index */
/*==========================================================================*/
inline size_t I_SPH(const size_t i_ph, const size_t i_th) {
   return i_ph*nlat_ + i_th;
}
size_t indx_Sph(const size_t i_ph, const size_t i_th) {
   return i_ph*nlat_ + i_th;
}
/*==========================================================================*/
void to_Ylm(const std::vector<double> &sph, std::vector<cplx> &ylm)
{
   for (size_t ip=0; ip<nphi_; ip++) {
   for (size_t it=0; it<nlat_; it++) {
      Sph_[I_SPH(ip,it)] = sph[I_SPH(ip,it)];
   }
   }
   spat_to_SH(shtns_, Sph_, Ylm_); 

   for (size_t il=0;   il<=lmax_;  il++) {
   for (size_t im=-il; im<=il;     im++) {
      ylm[LM(shtns_,il,im)] = Ylm_[LM(shtns_,il,im)]; 
   }
   }
}
/*==========================================================================*/
void to_Sph(const std::vector<cplx> &ylm, std::vector<double> &sph)
{
   for (size_t il=0;   il<=lmax_;  il++) {
   for (size_t im=-il; im<=il;     im++) {
      Ylm_[LM(shtns_,il,im)] = ylm[LM(shtns_,il,im)]; 
   }
   }
   SH_to_spat(shtns_, Ylm_, Sph_); 

   for (size_t ip=0; ip<nphi_; ip++) {
   for (size_t it=0; it<nlat_; it++) {
      sph[I_SPH(ip,it)] = Sph_[I_SPH(ip,it)];
   }
   }
}
/*==========================================================================*/
void laplace_beltrami(
      const std::vector<double> v, 
      std::vector<double> ddv)
{
   assert(v.size()  ==nlat_*nphi_);
   assert(ddv.size()==nlat_*nphi_);

   for (size_t ip=0; ip<nphi_; ip++) {
   for (size_t it=0; it<nlat_; it++) {
      Sph_[I_SPH(ip,it)] = v[I_SPH(ip,it)];
   }
   }
   spat_to_SH(shtns_, Sph_, Ylm_);
   for (size_t il=0;   il<=lmax_; il++) {
   for (size_t im=-il; im<=il;    im++) {
      Ylm_[LM(shtns_,il,im)] *= lap[il];
   }
   }
   SH_to_spat(shtns_, Ylm_, Sph_);
   for (size_t ip=0; ip<nphi_; ip++) {
   for (size_t it=0; it<nlat_; it++) {
      ddv[I_SPH(ip,it)] = Sph_[I_SPH(ip,it)]; 
   }
   }
}
/*==========================================================================*/
void cleanup()
{
   assert(Sph_ !=nullptr);
   assert(Ylm_!=nullptr);
   fftw_free(Sph_);
   fftw_free(Ylm_);

   shtns_destroy(shtns_);
}
/*==========================================================================*/
double nlat() { return nlat_; }
double nphi() { return nphi_; }
double nSph() { return nlat_*nphi_; }
double nYlm() { return shtns_->nlm;  }
/*==========================================================================*/
double theta(const size_t i_th)
{
   return acos(shtns_->ct[i_th]);
}
double phi(const size_t i_ph)
{
   return PHI_RAD(shtns_,i_ph);
}
/*==========================================================================*/
} /* Sphere */
