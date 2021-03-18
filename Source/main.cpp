#include "io.hpp"
#include "spherical_harmonics.hpp"

int main()
{
   shtns_cfg shtns;              // handle to a sht transform configuration
   long int lmax,mmax,nlat,nphi,mres, NLM;
   cmplx  *Slm, *Tlm;    // spherical harmonics coefficients (l,m space): complex numbers.
   double *Sh, *Th;              // real space : theta,phi
   long int i,im;
   double t;

   lmax = 5;       nlat = 32;
   mmax = 3;       nphi = 10;
   mres = 1;
   shtns_verbose(1);       // displays informations during initialization.
   shtns_use_threads(0);   // enable multi-threaded transforms (if supported).
   shtns = shtns_init( sht_gauss, lmax, mmax, mres, nlat, nphi );
   // shtns = shtns_create(lmax, mmax, mres, sht_orthonormal | SHT_REAL_NORM);
   // shtns_set_grid(shtns, sht_gauss, 0.0, nlat, nphi);
   NLM = shtns->nlm;

   // Memory allocation : the use of fftw_malloc is required because we need proper 16-byte alignement.
   // allocate spatial fields.
   Sh = (double *) fftw_malloc( NSPAT_ALLOC(shtns) * sizeof(double));
   Th = (double *) fftw_malloc( NSPAT_ALLOC(shtns) * sizeof(double));

   // allocate SH representations.
   Slm = (cmplx *) fftw_malloc( NLM * sizeof(cmplx));
   Tlm = (cmplx *) fftw_malloc( NLM * sizeof(cmplx));

   // SH_to_spat
   for (unsigned int lm=0; lm<shtns->nlm; lm++) {
      Slm[lm] = 0.0;  
      Tlm[lm] = 0.0;
   }

   //      Slm[LM(shtns, 1,1)] = sh11_st(shtns); // access to SH coefficient
   Slm[LM(shtns, 2,0)] = 1.0;
   //      Slm[LiM(shtns, 1,0)] = sh10_ct(shtns);
   //      Slm[LiM(shtns, 0,0)] = 0.5*sh00_1(shtns);
   SH_to_spat(shtns, Slm,Sh);
   SHtor_to_spat(shtns, Slm,Th,Sh);
   IO::write_vect("ylm",(double *) Slm,NLM*2);
   IO::write_mx("spat",Sh,nphi,nlat);

   // compute value of SH expansion at a given physical point.
   double t2;
   SHqst_to_point(shtns, Tlm, Tlm, Slm, shtns->ct[nlat/3], 2.*M_PI/(mres*nphi),&t2,&t2,&t);
   printf("check if SH_to_point coincides with SH_to_spat : %f = %f\n",t,Sh[nlat/3]);
   printf("ct*st = %f\n",shtns->ct[nlat/3]*shtns->st[nlat/3]);

   // check non-linear behaviour
   for (im=0;im<nphi;im++) {
      for (i=0;i<nlat;i++) {
         Sh[im*nlat+i] *= Sh[im*nlat+i];
      }
   }
   spat_to_SH(shtns, Sh, Tlm);
   IO::write_vect("ylm_nl",(double *) Tlm, NLM*2);

   // vector transform
   for (unsigned int lm=0; lm<shtns->nlm; lm++) {
      Slm[lm] = 0.0;  
      Tlm[lm] = 0.0;
   }
   SHsphtor_to_spat(shtns, Slm,Tlm, Sh,Th);  // vector transform
   IO::write_mx("spatt",Sh,nphi,nlat);
   IO::write_mx("spatp",Th,nphi,nlat);

   // spat_to_SH
   for (im=0;im<nphi;im++) {
      for (i=0;i<nlat;i++) {
         Sh[im*nlat+i] = shtns->ct[i]; // use cos(theta) array
      }
   }
   spat_to_SH(shtns, Sh,Slm);
   IO::write_vect("ylm_v",(double *) Slm,NLM*2);
}
