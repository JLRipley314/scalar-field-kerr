#ifndef _CHEB_HPP_
#define _CHEB_HPP_
/*===========================================================================*/
/* Stores Chebyshev points over interval [lower,upper]
 * and computes Chebyshev derivatives */
/*===========================================================================*/

#include <vector>
#include <fftw3.h>
/*===========================================================================*/
namespace Cheb {
/*---------------------------------------------------------------------------*/
   void init(const int n, const double lower, const double upper);
   void cleanup();

   void to_ch(const std::vector<double> &po, std::vector<double> &ch);
   void to_po(const std::vector<double> &ch, std::vector<double> &po);

   void der(const std::vector<double> &v, 
         std::vector<double> &ch, std::vector<double> &dv
         );

   void filter(std::vector<double> &ch, std::vector<double> &v);

   size_t n();   
   double lower();   
   double upper();   
/*===========================================================================*/
} /* Cheb */
/*===========================================================================*/
#endif /* _CHEB_HPP_ */
