#ifndef _CHEB_HPP_
#define _CHEB_HPP_

#include <vector>
#include <fftw3.h>
/*===========================================================================*/
/* Class to store Chebyshev points over interval [lower,upper]
 * and compute Chebyshev derivatives */
/*===========================================================================*/
class Cheb {
   private:
      const int n;

      const double lower;
      const double upper;
      const double jacobian; /* (upper-lower)/2.0 */
      /* 
       * Chebyshev points over interval [lower,upper] 
       */
      std::vector<double> pts;

      std::vector<double> low_pass;
      /* 
       * for fftw Fourier transform 
       */
      double *in;
      double *out;
      fftw_plan plan_dct;

   public:
      Cheb(const int n, 
            const double lower, 
            const double upper
         );
      ~Cheb();

      void to_ch(const std::vector<double> &po, std::vector<double> &ch);
      void to_po(const std::vector<double> &ch, std::vector<double> &po);

      void der(const std::vector<double> &v, 
            std::vector<double> &ch, std::vector<double> &dv
            );

      void filter(std::vector<double> &ch, std::vector<double> &v);

      inline double pt(const size_t i) { return pts[i]; }
};
/*===========================================================================*/
#endif /* _CHEB_HPP_ */
