/*
 * Computes finite differences 
 */
#ifndef FD_HPP__
#define FD_HPP__

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
};
/*===========================================================================*/
#endif /* _FD_HPP_ */
