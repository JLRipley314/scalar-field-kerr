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
   std::vector<double> in_( n_,0);
   std::vector<double> out_(n_,0);

   plan_dct_ = fftw_plan_r2r_1d(n_,in_.data(),out_.data(),FFTW_REDFT00,FFTW_PATIENT);

   assert(plan_dct_!=nullptr);
}
/*==========================================================================*/
void cleanup()
{
   assert(plan_dct_!=nullptr);

   fftw_destroy_plan(plan_dct_);
   fftw_cleanup();
}
/*==========================================================================*/
/* position space to Chebyshev space */
/*==========================================================================*/
void to_ch(std::vector<double> &po, std::vector<double> &ch) 
{
   /*
    * compute Fourier transform
    */
   assert(po.size()==n_);
   assert(ch.size()==n_);
   fftw_execute_r2r(plan_dct_, po.data(), ch.data());
   /*
    * normalize Chebyshev coefficients
    */
   ch[0]    /= (2.0*(n_-1));
   ch[n_-1] /= (2.0*(n_-1));
   for (size_t i=1; i<n_-1; i++) {
      ch[i] /= (n_-1);
   }
}
/*==========================================================================*/
/* Chebyshev space to position space */
/*==========================================================================*/
void to_po(std::vector<double> &ch, std::vector<double> &po) 
{
   /*
    * compute Fourier transform
    */
   assert(ch.size()==n_);
   assert(po.size()==n_);
   /*
    * normalize coefficients for Fourier transform 
    */
   for (size_t i=1; i<n_-1; i++) {
      ch[i] /= 2.0;
   }
   fftw_execute_r2r(plan_dct_, ch.data(), po.data());
}
/*==========================================================================*/
/* Compute derivative over interval */
/*==========================================================================*/
void der(std::vector<double> &v, std::vector<double> &dv)
{
   std::vector<double> ch_tmp(n_,0);
   assert(v.size( )==n_);
   assert(dv.size()==n_);

   to_ch(v,ch_tmp);
   /*
    * to start Chebyshev derivative recurrence relation 
    */
   ch_tmp[n_-1] = 0;
   ch_tmp[n_-2] = 0;
   /* 
    * use dv as a temporary array
    */
   for (size_t i=0; i<n_; i++) {
      dv[i] = ch_tmp[i];
   }
   /* 
    * apply Chebyshev derivative recurrence relation 
    */
   for (size_t i=n_-2; i>=1; i--) { 
      ch_tmp[i-1] = 2.0*i*dv[i] + ch_tmp[i+1];
   } 
   ch_tmp[0] /= 2.0;
   /* 
    * Normalize derivative to inverval 
    */
   to_po(ch_tmp,dv);
   for (size_t i=0; i<n_; i++) {
      dv[i] /= jacobian_;
   }
}
/*==========================================================================*/
/* Low pass filter of Chebyshev coefficients */
/*==========================================================================*/
void filter(std::vector<double> &v)
{
   std::vector<double> ch_tmp(n_,0);
   assert(v.size()==n_);

   to_ch(v,ch_tmp);
   for (size_t i=0; i<n_; i++) { 
      ch_tmp[i] *= low_pass_[i];
   } 
   to_po(ch_tmp,v);
}
/*==========================================================================*/
size_t n() { return n_; }
double lower() { return lower_; }
double upper() { return upper_; }
double pt(const size_t i) { return pts_[i]; }
/*==========================================================================*/
} /* Cheb */
