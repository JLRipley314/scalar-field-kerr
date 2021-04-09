#include <cmath>
#include <cassert>
#include <iostream>
#include "cheb.hpp"
/*==========================================================================*/
namespace Cheb {
/*==========================================================================*/
namespace {
      size_t _n;

      double _lower;
      double _upper;
      double _jacobian; 
      /* 
       * Chebyshev points over interval [lower,upper] 
       */
      std::vector<double> _pts;
      std::vector<double> _low_pass;
      /* 
       * for fftw Fourier transform 
       */
      fftw_plan _plan_dct;
}
/*==========================================================================*/
/* setup fftw and vectors */
/*==========================================================================*/
void init(const int n, const double lower, const double upper) 
{
   _n = n;
   _lower = lower;
   _upper = upper;
   _jacobian = (_upper-_lower) / 2.0;

   _pts.resize(_n,0.0);
   _low_pass.resize(_n,0.0);
   /* 
    * Chebyshev points on interval [lower, upper]
    */
   for (size_t i=0; i<_n; i++) {
      _pts[i] = _jacobian*cos(M_PI*double(i)/(_n-1.0)) + ((_upper+_lower)/2.0);
   }
   /* 
    * for low pass filter in Chebyshev space 
    */
   for (size_t i=0; i<_n; i++) {
      _low_pass[i] = exp(-40.0*pow(double(i)/_n,10));
   }   
   /* 
    * initialize fftw discrete cosine transform plan 
    */
   std::vector<double> i_n( _n,0);
   std::vector<double> out_(_n,0);

   _plan_dct = fftw_plan_r2r_1d(_n,i_n.data(),out_.data(),FFTW_REDFT00,FFTW_PATIENT);

   assert(_plan_dct!=nullptr);
}
/*==========================================================================*/
void cleanup()
{
   assert(_plan_dct!=nullptr);

   fftw_destroy_plan(_plan_dct);
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
   assert(po.size()==_n);
   assert(ch.size()==_n);
   fftw_execute_r2r(_plan_dct, po.data(), ch.data());
   /*
    * normalize Chebyshev coefficients
    */
   ch[0]    /= (2.0*(_n-1));
   ch[_n-1] /= (2.0*(_n-1));
   for (size_t i=1; i<_n-1; i++) {
      ch[i] /= (_n-1);
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
   assert(ch.size()==_n);
   assert(po.size()==_n);
   /*
    * normalize coefficients for Fourier transform 
    */
   for (size_t i=1; i<_n-1; i++) {
      ch[i] /= 2.0;
   }
   fftw_execute_r2r(_plan_dct, ch.data(), po.data());
}
/*==========================================================================*/
/* Compute derivative over interval */
/*==========================================================================*/
void der(std::vector<double> &v, std::vector<double> &dv)
{
   std::vector<double> ch_tmp(_n,0);
   assert(v.size( )==_n);
   assert(dv.size()==_n);

   to_ch(v,ch_tmp);
   /*
    * to start Chebyshev derivative recurrence relation 
    */
   ch_tmp[_n-1] = 0;
   ch_tmp[_n-2] = 0;
   /* 
    * use dv as a temporary array
    */
   for (size_t i=0; i<_n; i++) {
      dv[i] = ch_tmp[i];
   }
   /* 
    * apply Chebyshev derivative recurrence relation 
    */
   for (size_t i=_n-2; i>=1; i--) { 
      ch_tmp[i-1] = 2.0*i*dv[i] + ch_tmp[i+1];
   } 
   ch_tmp[0] /= 2.0;
   /* 
    * Normalize derivative to inverval 
    */
   to_po(ch_tmp,dv);
   for (size_t i=0; i<_n; i++) {
      dv[i] /= _jacobian;
   }
}
/*==========================================================================*/
/* Low pass filter of Chebyshev coefficients */
/*==========================================================================*/
void filter(std::vector<double> &v)
{
   std::vector<double> ch_tmp(_n,0);
   assert(v.size()==_n);

   to_ch(v,ch_tmp);
   for (size_t i=0; i<_n; i++) { 
      ch_tmp[i] *= _low_pass[i];
   } 
   to_po(ch_tmp,v);
}
/*==========================================================================*/
size_t n() { return _n; }
double lower() { return _lower; }
double upper() { return _upper; }
double pt(const size_t i) { return _pts[i]; }
/*==========================================================================*/
} /* Cheb */
