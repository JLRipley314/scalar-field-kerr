#ifndef _CHEB_HPP_
#define _CHEB_HPP_

#include <vector>
#include <fftw3.h>
/*===========================================================================*/
/* Class to store Chebyshev points over interval [lower,upper]
 * and compute Chebyshev derivatives */
/*===========================================================================*/
namespace Cheb {
/*---------------------------------------------------------------------------*/
   namespace {
         int n_;

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
/*---------------------------------------------------------------------------*/
   void init(const int n, const double lower, const double upper);
   void cleanup();

   void to_ch(const std::vector<double> &po, std::vector<double> &ch);
   void to_po(const std::vector<double> &ch, std::vector<double> &po);

   void der(const std::vector<double> &v, 
         std::vector<double> &ch, std::vector<double> &dv
         );

   void filter(std::vector<double> &ch, std::vector<double> &v);

   inline double pt(const size_t i) { return pts_[i]; }
};
/*===========================================================================*/
#endif /* _CHEB_HPP_ */
