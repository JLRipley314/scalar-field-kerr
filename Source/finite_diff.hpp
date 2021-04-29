/*
 * Computes finite differences 
 */
#ifndef _FD_HPP_
#define _FD_HPP_

#include <vector>

class FD {
public:
   FD(const size_t n, const double lower, const double upper);
   ~FD();

   size_t n() const;
   double lower() const; 
   double upper() const; 
   double pt(const size_t i) const;

   void der(   std::vector<double> &v, std::vector<double> &dv ) const;
   void der2(  std::vector<double> &v, std::vector<double> &ddv) const;
   void filter(std::vector<double> &v) const;

private:
   const size_t _n;

   const double _lower;
   const double _upper;
   const double _dx; 
   /* 
    * Points over interval [lower,upper] 
    */
   std::vector<double> _pts;
};
/*===========================================================================*/
#endif /* _FD_HPP_ */
