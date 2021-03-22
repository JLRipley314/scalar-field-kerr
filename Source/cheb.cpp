#include <cmath>
#include <cassert>
#include "cheb.hpp"

/*==========================================================================*/
namespace Cheb {
/*==========================================================================*/
namespace {
      size_t n_;

      double lower_;
      double upper_;
      double jacobian_; 
      /* 
       * Chebyshev points over interval [lower,upper] 
       */
      std::vector<double> pts_;

      std::vector<double> low_pass_;
      /* 
       * for fftw Fourier transform 
       */
      double *in_;
      double *out_;
      fftw_plan plan_dct_;
}
/*==========================================================================*/
/* setup fftw and vectors */
/*==========================================================================*/
void init(const int n, const double lower, const double upper) 
{
   n_ = n;
   lower_ = lower;
   upper_ = upper;
   jacobian_ = (upper_-lower_) / 2.0;

   pts_.resize(n_,0.0);
   low_pass_.resize(n_,0.0);
   /* 
    * Chebyshev points on interval [lower, upper]
    */
   for (size_t i=0; i<n_; i++) {
      pts_[i] = jacobian_*cos(M_PI*double(i)/(n_-1)) + ((upper_+lower_)/2.0);
   }
   /* 
    * for low pass filter in Chebyshev space 
    */
   for (size_t i=0; i<n_; i++) {
      low_pass_[i] = exp(-40.0*pow(double(i)/n_,10));
   }   
   /* 
    * initialize fftw discrete cosine transform plan 
    */
   in_  = (double *)fftw_malloc(n_*sizeof(double));
   out_ = (double *)fftw_malloc(n_*sizeof(double));
   assert(in_ !=nullptr);
   assert(out_!=nullptr);

   for (size_t i=0; i<n_; i++) {
      in_[i]  = 0;
      out_[i] = 0;
   }
   plan_dct_ = fftw_plan_r2r_1d(n_,in_,out_,FFTW_REDFT00,FFTW_PATIENT);

   assert(plan_dct_!=nullptr);
}
/*==========================================================================*/
void cleanup()
{
   assert(in_      !=nullptr);
   assert(out_     !=nullptr);
   assert(plan_dct_!=nullptr);

   fftw_free(in_);
   fftw_free(out_);
   fftw_destroy_plan(plan_dct_);
   fftw_cleanup();
}
/*==========================================================================*/
/* position space to Chebyshev space */
/*==========================================================================*/
void to_ch(const std::vector<double> &po, std::vector<double> &ch) 
{
   /*
    * compute Fourier transform
    */
   assert(po.size()==n_);
   assert(ch.size()==n_);
   for (size_t i=0; i<n_; i++) {
      in_[i]  = po[i];
      out_[i] = ch[i];
   }
   fftw_execute(plan_dct_);
   /*
    * normalize Chebyshev coefficients
    */
   ch[0]    = out_[0]   /(2.0*(n_-1));
   ch[n_-1] = out_[n_-1]/(2.0*(n_-1));
   for (size_t i=1; i<n_-1; i++) {
      ch[i] = out_[i]/(n_-1);
   }
}
/*==========================================================================*/
/* Chebyshev space to position space */
/*==========================================================================*/
void to_po(const std::vector<double> &ch, std::vector<double> &po) 
{
   /*
    * compute Fourier transform
    */
   assert(ch.size()==n_);
   assert(po.size()==n_);
   for (size_t i=0; i<n_; i++) {
      in_[i]  = ch[i];
      out_[i] = po[i];
   }
   /*
    * normalize coefficients for Fourier transform 
    */
   for (size_t i=1; i<n_-1; i++) {
      in_[i] /= 2.0;
   }
   fftw_execute(plan_dct_);

   for (size_t i=0; i<n_; i++) {
      po[i] = out_[i];
   }
}
/*==========================================================================*/
/* Compute derivative over interval */
/*==========================================================================*/
void der(const std::vector<double> &v, 
      std::vector<double> &ch, std::vector<double> &dv)
{
   assert(v.size()==n_);
   assert(ch.size()==n_);
   assert(dv.size()==n_);

   to_ch(v,ch);
   /*
    * to start Chebyshev derivative recurrence relation 
    */
   ch[n_-1] = 0;
   ch[n_-2] = 0;
   /* 
    * use dv as a temporary array
    */
   for (size_t i=0; i<n_; i++) {dv[i] = ch[i];}
   /* 
    * apply Chebyshev derivative recurrence relation 
    */
   for (size_t i=n_-2; i>=1; i--) { 
      ch[i-1] = 2.0*(i-1)*dv[i] + ch[i+1];
   } 
   ch[0] /= 2.0;
   /* 
    * Normalize derivative to inverval 
    */
   to_po(ch,dv);
   for (size_t i=0; i<n_; i++) {
      dv[i] *= jacobian_;
   }
}
/*==========================================================================*/
/* Low pass filter of Chebyshev coefficients */
/*==========================================================================*/
void filter(std::vector<double> &ch, std::vector<double> &v)
{
   assert(ch.size()==n_);
   assert(v.size()==n_);

   to_ch(v,ch);
   for (size_t i=0; i<n_; i++) { 
      ch[i] *= low_pass_[i];
   } 
   to_po(ch,v);
}
/*==========================================================================*/
size_t n() { return n_; }
double lower() { return lower_; }
double upper() { return upper_; }
double pt(const size_t i)    { return pts_[i]; }
/*==========================================================================*/
} /* Cheb */
