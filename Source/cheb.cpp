#include <cmath>
#include <cassert>
#include "cheb.hpp"

/* 
 * for thread safe versions
 */
#define ALLOCATE_TMP \
   double *in_tmp  = (double *)fftw_malloc(n_*sizeof(double)); \
   double *out_tmp = (double *)fftw_malloc(n_*sizeof(double)); \
   assert(in_tmp !=nullptr); \
   assert(out_tmp!=nullptr);
#define FREE_TMP \
   fftw_free(in_tmp); \
   fftw_free(out_tmp);
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
       * stand-in for chebyshev coefficients and position space
       * for some routines 
       */
      std::vector<double> internal_po_;
      std::vector<double> internal_ch_;
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
   internal_po_.resize(n_,0.0);
   internal_ch_.resize(n_,0.0);
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
   fftw_execute_r2r(plan_dct_, in_, out_);
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
void to_ch_ts(const std::vector<double> &po, std::vector<double> &ch) 
{
   ALLOCATE_TMP;
   /*
    * compute Fourier transform
    */
   assert(po.size()==n_);
   assert(ch.size()==n_);
   for (size_t i=0; i<n_; i++) {
      in_tmp[i]  = po[i];
      out_tmp[i] = ch[i];
   }
   fftw_execute_r2r(plan_dct_, in_tmp, out_tmp);
   /*
    * normalize Chebyshev coefficients
    */
   ch[0]    = out_tmp[0]   /(2.0*(n_-1));
   ch[n_-1] = out_tmp[n_-1]/(2.0*(n_-1));
   for (size_t i=1; i<n_-1; i++) {
      ch[i] = out_tmp[i]/(n_-1);
   }
   FREE_TMP;
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
   fftw_execute_r2r(plan_dct_, in_, out_);

   for (size_t i=0; i<n_; i++) {
      po[i] = out_[i];
   }
}
/*==========================================================================*/
void to_po_ts(const std::vector<double> &ch, std::vector<double> &po) 
{
   ALLOCATE_TMP;
   /*
    * compute Fourier transform
    */
   assert(ch.size()==n_);
   assert(po.size()==n_);
   for (size_t i=0; i<n_; i++) {
      in_tmp[i]  = ch[i];
      out_tmp[i] = po[i];
   }
   /*
    * normalize coefficients for Fourier transform 
    */
   for (size_t i=1; i<n_-1; i++) {
      in_tmp[i] /= 2.0;
   }
   fftw_execute_r2r(plan_dct_, in_tmp, out_tmp);

   for (size_t i=0; i<n_; i++) {
      po[i] = out_tmp[i];
   }
   FREE_TMP;
}
/*==========================================================================*/
/* Compute derivative over interval */
/*==========================================================================*/
void der(const std::vector<double> &v, std::vector<double> &dv)
{
   assert(v.size( )==n_);
   assert(dv.size()==n_);

   to_ch(v,internal_ch_);
   /*
    * to start Chebyshev derivative recurrence relation 
    */
   internal_ch_[n_-1] = 0;
   internal_ch_[n_-2] = 0;
   /* 
    * use dv as a temporary array
    */
   for (size_t i=0; i<n_; i++) {
      dv[i] = internal_ch_[i];
   }
   /* 
    * apply Chebyshev derivative recurrence relation 
    */
   for (size_t i=n_-2; i>=1; i--) { 
      internal_ch_[i-1] = 2.0*i*dv[i] + internal_ch_[i+1];
   } 
   internal_ch_[0] /= 2.0;
   /* 
    * Normalize derivative to inverval 
    */
   to_po(internal_ch_,dv);
   for (size_t i=0; i<n_; i++) {
      dv[i] /= jacobian_;
   }
}
/*==========================================================================*/
void der_ts(const std::vector<double> &v, std::vector<double> &dv)
{
   std::vector<double> internal_ch_tmp(n_,0);

   assert(v.size( )==n_);
   assert(dv.size()==n_);

   to_ch_ts(v,internal_ch_tmp);
   /*
    * to start Chebyshev derivative recurrence relation 
    */
   internal_ch_tmp[n_-1] = 0;
   internal_ch_tmp[n_-2] = 0;
   /* 
    * use dv as a temporary array
    */
   for (size_t i=0; i<n_; i++) {
      dv[i] = internal_ch_tmp[i];
   }
   /* 
    * apply Chebyshev derivative recurrence relation 
    */
   for (size_t i=n_-2; i>=1; i--) { 
      internal_ch_tmp[i-1] = 2.0*i*dv[i] + internal_ch_tmp[i+1];
   } 
   internal_ch_tmp[0] /= 2.0;
   /* 
    * Normalize derivative to inverval 
    */
   to_po_ts(internal_ch_tmp,dv);
   for (size_t i=0; i<n_; i++) {
      dv[i] /= jacobian_;
   }
}
/*==========================================================================*/
/* Low pass filter of Chebyshev coefficients */
/*==========================================================================*/
void filter(std::vector<double> &v)
{
   assert(v.size()==n_);

   to_ch(v,internal_ch_);
   for (size_t i=0; i<n_; i++) { 
      internal_ch_[i] *= low_pass_[i];
   } 
   to_po(internal_ch_,v);
}
/*==========================================================================*/
void filter_ts(std::vector<double> &v)
{
   std::vector<double> internal_ch_tmp(n_,0);
   assert(v.size()==n_);

   to_ch_ts(v,internal_ch_tmp);
   for (size_t i=0; i<n_; i++) { 
      internal_ch_tmp[i] *= low_pass_[i];
   } 
   to_po_ts(internal_ch_tmp,v);
}
/*==========================================================================*/
size_t n() { return n_; }
double lower() { return lower_; }
double upper() { return upper_; }
double pt(const size_t i) { return pts_[i]; }
/*==========================================================================*/
} /* Cheb */

#undef ALLOCATE_TMP
#undef FREE_TMP
