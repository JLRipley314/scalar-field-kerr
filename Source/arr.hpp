/*
 * Routines to manipulate 3d arrays 
 */
#ifndef _ARR_HPP_
#define _ARR_HPP_

#include <vector>
/*=========================================================================*/
namespace Arr3d {
/*=========================================================================*/
namespace {
   size_t nx_;
   size_t ny_;
   size_t nz_;
}
/*=========================================================================*/
void init(const size_t nx, const size_t ny, const size_t nz);
void cleanup();

std::vector<double> arr3d(const double val);

inline size_t indx(const size_t i, const size_t j, const size_t k) 
{
   return ny_*nz_*i + nz_*j + k;
}

void row1(const size_t j, const size_t k, 
      const std::vector<double> &in,
      std::vector<double> &out
   ); 
void row2(const size_t i, const size_t k, 
      const std::vector<double> &in,
      std::vector<double> &out
   ); 
void row3(const size_t i, const size_t j, 
      const std::vector<double> &in,
      std::vector<double> &out
   ); 
void row12(const size_t k, 
      const std::vector<double> &in,
      std::vector<double> &out
   );
void row13(const size_t j, 
      const std::vector<double> &in,
      std::vector<double> &out
   );
void row23(const size_t i, 
      const std::vector<double> &in,
      std::vector<double> &out
   );
/*=========================================================================*/
}; /* Arr3d */
/*=========================================================================*/
#endif /* _ARR_HPP_ */
