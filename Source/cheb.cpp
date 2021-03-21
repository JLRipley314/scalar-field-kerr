#include <cmath>
#include <cassert>
#include "cheb.hpp"

/*==========================================================================*/
Cheb::Cheb(const int n, const double jacobian) 
:  n{n},
   jacobian{jacobian},
   low_pass(n,0)
{
   /* 
    * for low pass filter 
    */
   for (int i=0; i<n; i++) {
      low_pass[i] = exp(-40.0*pow(double(i)/n,10));
   }   
   /* 
    * initialize fftw discrete cosine transform plan 
    */
   in  = (double *)fftw_malloc(n*sizeof(double));
   out = (double *)fftw_malloc(n*sizeof(double));
   assert(in!=nullptr);
   assert(in!=nullptr);

   for (int i=0; i<n; i++) {
      in[i]  = 0;
      out[i] = 0;
   }

   plan_dct = fftw_plan_r2r_1d(n,in,out,FFTW_REDFT00,FFTW_PATIENT);
}
/*==========================================================================*/
/* position space to Chebyshev space */
/*==========================================================================*/
void Cheb::to_ch(const std::vector<double> &po, std::vector<double> &ch) 
{
   /*
    * compute Fourier transform
    */
   assert(po.size()==n);
   assert(ch.size()==n);
   for (int i=0; i<n; i++) {
      in[i]  = po[i];
      out[i] = ch[i];
   }
   fftw_execute(plan_dct);
   /*
    * normalize Chebyshev coefficients
    */
   ch[0]   = out[0]  /2.0*(n-1);
   ch[n-1] = out[n-1]/2.0*(n-1);
   for (int i=1; i<n-1; i++) {
      ch[i] = out[i]/(n-1);
   }
}
/*==========================================================================*/
/* Chebyshev space to position space */
/*==========================================================================*/
void Cheb::to_po(const std::vector<double> &ch, std::vector<double> &po) 
{
   /*
    * compute Fourier transform
    */
   assert(ch.size()==n);
   assert(po.size()==n);
   for (int i=0; i<n; i++) {
      in[i]  = ch[i];
      out[i] = po[i];
   }
   fftw_execute(plan_dct);
   /*
    * normalize transform to position space 
    */
   po[0]   = out[0];
   po[n-1] = out[n-1];
   for (int i=1; i<n-1; i++) {
      po[i] = out[i]/2.0;
   }
}
/*==========================================================================*/
/* Compute derivative over interval */
/*==========================================================================*/
void Cheb::der(
      const std::vector<double> &v, 
      std::vector<double> &ch, std::vector<double> &dv)
{
   assert(v.size()==n);
   assert(ch.size()==n);
   assert(dv.size()==n);

   to_ch(v,ch);
   /*
    * to start Chebyshev derivative recurrence relation 
    */
   ch[n-1] = 0;
   ch[n-2] = 0;
   /* 
    * use dv as a temporary array
    */
   for (int i=0; i<n; i++) {dv[i] = ch[i];}
   /* 
    * apply Chebyshev derivative recurrence relation 
    */
   for (int i=n-2; i>=1; i--) { 
      ch[i-1] = 2.0*(i-1)*dv[i] + ch[i+1];
   } 
   ch[0] /= 2.0;
   /* 
    * Normalize derivative to inverval 
    */
   to_po(ch,dv);
   for (int i=0; i<n; i++) {
      dv[i] *= jacobian;
   }
}
/*==========================================================================*/
/* Low pass filter of Chebyshev coefficients */
/*==========================================================================*/
void Cheb::filter(
      std::vector<double> &ch, std::vector<double> &v)
{
   assert(ch.size()==n);
   assert(v.size()==n);

   to_ch(v,ch);
   for (int i=0; i<n; i++) { 
      ch[i] *= low_pass[i];
   } 
   to_po(ch,v);
}
/*==========================================================================*/
Cheb::~Cheb()
{
   fftw_free(in);
   fftw_free(out);
   fftw_destroy_plan(plan_dct);
}
