#include <cassert>
#include "sphere.hpp"

#include <iostream>
#include <iomanip>
/*==========================================================================*/
Sphere::Sphere(
      const size_t nl, 
      const size_t nm, 
      const size_t nlat, 
      const size_t nphi)
:  _lmax{nl},
   _mmax{nm},
   _nlat{nlat}, 
   _nphi{nphi}, 
   _lap(_lmax,0),
   _low_pass(_lmax,0)
{
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
   for (size_t i=0; i<_lmax; i++) {
      _lap[i] = - double(i)*double(i+1.0);
   }
   /* 
    * Low pass filter in spherical harmonic space 
    */
   for (size_t i=0; i<_lmax; i++) {
      _low_pass[i] = exp(-40.0*pow(double(i)/_lmax,16));
   }
   _nSph = NSPAT_ALLOC(_shtns);
   _nYlm = _shtns->nlm;
}
/*==========================================================================*/
Sphere::~Sphere()
{
   shtns_destroy(_shtns);
}
/*==========================================================================*/
size_t Sphere::nlat() const { return _nlat; }
size_t Sphere::nphi() const { return _nphi; }
size_t Sphere::nSph() const { return _nSph; }
size_t Sphere::nYlm() const { return _nYlm; }
/*==========================================================================*/
double Sphere::theta(const size_t i_th) const
{
   return acos(_shtns->ct[i_th]);
}
double Sphere::phi(const size_t i_ph) const
{
   return i_ph*2.0*M_PI/((_shtns->nphi)*(_shtns->mres));
}
/*==========================================================================*/
/* Need to use fftw_malloc instead of vectors
 * when using shtns. */
/*==========================================================================*/
double* Sphere::allocate_real(const size_t size) const
{
   double *tmp = (double *) fftw_malloc(size * sizeof(double));
   assert(tmp!=nullptr);
   return tmp;
}
/*==========================================================================*/
cplx* Sphere::allocate_cplx(const size_t size) const
{
   cplx *tmp = (cplx *) fftw_malloc(size * sizeof(cplx));
   assert(tmp!=nullptr);
   return tmp;
}
/*==========================================================================*/
void Sphere::free(double *tmp) const
{
   fftw_free(tmp);
   tmp = nullptr;
}
/*==========================================================================*/
void Sphere::free(cplx *tmp) const
{
   fftw_free(tmp);
   tmp = nullptr;
}
/*==========================================================================*/
void Sphere::to_Ylm(const std::vector<double> &sph, std::vector<cplx> &ylm) const
{
   assert(sph.size()==_nSph);
   assert(ylm.size()==_nYlm);

   double *Sph_tmp = allocate_real(_nSph);
   cplx   *Ylm_tmp = allocate_cplx(_nYlm);

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
void Sphere::to_Sph(const std::vector<cplx> &ylm, std::vector<double> &sph) const
{
   assert(sph.size()==_nSph);
   assert(ylm.size()==_nYlm);

   double *Sph_tmp = allocate_real(_nSph);
   cplx   *Ylm_tmp = allocate_cplx(_nYlm);

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
void Sphere::to_Ylm(const std::vector<cplx> &sph, std::vector<cplx> &ylm) const
{
   assert(sph.size()==_nSph);
   assert(ylm.size()==pow(_lmax+1,2));

   cplx *Sph_tmp = allocate_cplx(_nSph);
   cplx *Ylm_tmp = allocate_cplx(pow(_lmax+1,2));

   for (size_t i=0; i<_nSph; i++) {
      Sph_tmp[i] = sph[i];
   }
   spat_cplx_to_SH(_shtns, Sph_tmp, Ylm_tmp); 

   for (size_t i=0; i<pow(_lmax+1,2); i++) {
      ylm[i] = Ylm_tmp[i]; 
   }
   free(Sph_tmp);
   free(Ylm_tmp);
}
/*==========================================================================*/
void Sphere::to_Sph(const std::vector<cplx> &ylm, std::vector<cplx> &sph) const
{
   assert(sph.size()==_nSph);
   assert(ylm.size()==pow(_lmax+1,2));

   cplx *Sph_tmp = allocate_cplx(_nSph);
   cplx *Ylm_tmp = allocate_cplx(pow(_lmax+1,2));

   for (size_t i=0; i<pow(_lmax+1,2); i++) {
      Ylm_tmp[i] = ylm[i];
   }
   SH_to_spat_cplx(_shtns, Ylm_tmp, Sph_tmp); 

   for (size_t i=0; i<_nSph; i++) {
      sph[i] = Sph_tmp[i];
   }
   free(Sph_tmp);
   free(Ylm_tmp);
}
/*==========================================================================*/
void Sphere::laplace_beltrami(
      const std::vector<double> &v, 
      std::vector<double> &ddv) const
{
   double *Sph_tmp = allocate_real(_nSph);
   cplx   *Ylm_tmp = allocate_cplx(_nYlm);

   assert(v.size()  ==_nSph);
   assert(ddv.size()==_nSph);

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
void Sphere::partial_phi(const std::vector<double> &v, std::vector<double> &dv) const
{
   double *Sph_tmp = allocate_real(_nSph);
   cplx   *Ylm_tmp = allocate_cplx(_nYlm);

   assert(v.size() ==_nSph);
   assert(dv.size()==_nSph);

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
/* Raising operator: notice that you obtain a complex function
 * L_+ Y_{l,m} = sqrt((l-m)(l+m+1)) Y_{l,m+1}. 
 * For the coefficients, that means
 * (L_+ f)_{l,m} = sqrt((l-m+1)(l+m)) f_{l,m-1},
 * so that (L_+ f)_{l,-l} = 0 */
/*==========================================================================*/
void Sphere::raise(const std::vector<double> &v, std::vector<cplx> &rv) const
{
   cplx *Sph_tmp = allocate_cplx(_nSph);
   cplx *Ylm_tmp = allocate_cplx(pow(_lmax+1,2));
   assert(v.size() ==_nSph);
   assert(rv.size()==_nSph);

   for (size_t i=0; i<_nSph; i++) {
      Sph_tmp[i] = cplx(v[i],0);
   }
   spat_cplx_to_SH(_shtns, Sph_tmp, Ylm_tmp);

   for (size_t l=0; l<_lmax; l++) {
      for (int m=int(l); m>-int(l); m--) {
         Ylm_tmp[l*(l+1)+m] = pow((l-m+1)*(l+m),0.5)*Ylm_tmp[l*(l+1)+(m-1)];
      }
      Ylm_tmp[l*(l+1)+(-l)] = 0;
   }
   SH_to_spat_cplx(_shtns, Ylm_tmp, Sph_tmp);

   for (size_t i=0; i<_nSph; i++) {
      rv[i] = Sph_tmp[i];
   }
   free(Sph_tmp);
   free(Ylm_tmp);
}
/*==========================================================================*/
/* Lowering operator: notice that you obtain a complex function
 * L_- Y_{l,m} = sqrt((l+m)(l-m+1)) Y_{l,m-1}. 
 * For the coefficients, that means
 * (L_- f)_{l,m} = sqrt((l+m+1)(l-m)) f_{l,m+1},
 * so that (L_- f)_{l,l} = 0 */
/*==========================================================================*/
void Sphere::lower(const std::vector<double> &v, std::vector<cplx> &lv) const
{
   cplx *Sph_tmp = allocate_cplx(_nSph);
   cplx *Ylm_tmp = allocate_cplx(pow(_lmax+1,2));

   assert(v.size() ==_nSph);
   assert(lv.size()==_nSph);

   for (size_t i=0; i<_nSph; i++) {
      Sph_tmp[i] = cplx(v[i],0);
   }
   spat_cplx_to_SH(_shtns, Sph_tmp, Ylm_tmp);

   for (size_t l=0; l<_lmax; l++) {
      for (int m=-int(l); m<int(l); m++) {
         Ylm_tmp[l*(l+1)+m] = pow((l+m+1)*(l-m),0.5)*Ylm_tmp[l*(l+1)+(m+1)];
      }
      Ylm_tmp[l*(l+1)+(l)] = 0;
   }
   SH_to_spat_cplx(_shtns, Ylm_tmp, Sph_tmp);

   for (size_t i=0; i<_nSph; i++) {
      lv[i] = Sph_tmp[i];
   }
   free(Sph_tmp);
   free(Ylm_tmp);
}
/*==========================================================================*/
/* SphereX operator: this is
 * (\partial_{\theta} f)^2 + (1/sin^2(theta)) (\partial_{\phi}f)^2 
 * =
 * - (L_+f)(L_-f) + (\partial_{\phi}f)^2 */
/*==========================================================================*/
void Sphere::sphereX(const std::vector<double> &v, std::vector<double> &vX) const
{
   assert(v.size() ==_nSph);
   assert(vX.size()==_nSph);

   std::vector<cplx>   raise_v(_nSph,0);
   std::vector<cplx>   lower_v(_nSph,0);
   std::vector<double> partial_phi_v(_nSph,0);

   raise(v,             raise_v);
   lower(v,             lower_v);
   partial_phi(v, partial_phi_v);

   for (size_t i=0; i<_nSph; i++) {
      cplx vrvl = raise_v[i]*lower_v[i]; /* This should be real */
      vX[i] = -vrvl.real() + pow(partial_phi_v[i],2);
   }
}
/*==========================================================================*/
/* low pass filter in spherical harmonic coefficient space
 * Note that only positive m are stored in spherical harmonic space,
 * as we only deal with real scalar fields. */
/*==========================================================================*/
void Sphere::filter(std::vector<double> &v) const
{
   assert(v.size()==_nSph);

   double *Sph_tmp = allocate_real(_nSph);
   cplx   *Ylm_tmp = allocate_cplx(_nYlm);

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
void Sphere::power_spectrum(const std::vector<double> &v, std::vector<double> &p) const
{
   assert(v.size()==_nSph);
   assert(p.size()==_lmax);

   double *Sph_tmp = allocate_real(_nSph);
   cplx   *Ylm_tmp = allocate_cplx(_nYlm);

   for (size_t i=0; i<_nSph; i++) {
      Sph_tmp[i] = v[i];
   }
   spat_to_SH(_shtns, Sph_tmp, Ylm_tmp);

   for (size_t l=0; l<_lmax; l++) {

      const double r_Ylm = Ylm_tmp[LM(_shtns,l,0)].real();
      const double i_Ylm = Ylm_tmp[LM(_shtns,l,0)].imag();
      p[l] = (pow(r_Ylm,2) + pow(i_Ylm,2));

      for (size_t m=1; m<std::min(l+1,_mmax); m++) {
         /*
          * We double count positive m as we do not
          * count negative m.
          */
         const double r_Ylm = Ylm_tmp[LM(_shtns,l,m)].real();
         const double i_Ylm = Ylm_tmp[LM(_shtns,l,m)].imag();

         p[l] += 2.0*(pow(r_Ylm,2) + pow(i_Ylm,2));
      }
   }
   free(Sph_tmp);
   free(Ylm_tmp);
}
/*==========================================================================*/
/* returns values for spherical harmonic in real space Y_{l_ang,m_ang} */
/*==========================================================================*/
std::vector<double> Sphere::compute_ylm(const size_t l_ang, const size_t m_ang) const
{
   std::vector<cplx>   in( _nYlm,0.0);
   std::vector<double> out(_nSph,0.0);

   in[LM(_shtns,l_ang,m_ang)] = 1;

   to_Sph(in, out);

   return out;
}
