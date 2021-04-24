/*
 * Stores Chebyshev points over interval [lower,upper]
 * and computes Chebyshev derivatives 
 */
#ifndef _CHEB_HPP_
#define _CHEB_HPP_

#include <vector>
#include <fftw3.h>

class Cheb {
public:
   Cheb(const size_t n, const double lower, const double upper);
   ~Cheb();

   size_t n() const;
   double lower() const; 
   double upper() const; 
   double pt(const size_t i) const;

   void to_ch( std::vector<double> &po, std::vector<double> &ch) const;
   void to_po( std::vector<double> &ch, std::vector<double> &po) const;
   void der(   std::vector<double> &v, std::vector<double> &dv) const;
   void filter(std::vector<double> &v) const;

private:
   const size_t _n;

   const double _lower;
   const double _upper;
   const double _jacobian; 
   /* 
    * Chebyshev points over interval [lower,upper] 
    */
   std::vector<double> _pts;
   std::vector<double> _low_pass;
   /* 
    * for fftw Fourier transform 
    */
   fftw_plan _plan_dct;
};
/*===========================================================================*/
#endif /* _CHEB_HPP_ */
