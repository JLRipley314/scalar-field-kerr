/*
 * Routines to manipulate 3d arrays 
 */
#ifndef _ARR_HPP_
#define _ARR_HPP_

#include <vector>
/*=========================================================================*/
namespace Arr3d 
{
/*=========================================================================*/
namespace {
   size_t _ny;
   size_t _nz;
}
/*=========================================================================*/
void init(const size_t nx, const size_t ny, const size_t nz);
void cleanup();

std::vector<double> arr3d(const double val);

inline size_t indx(const size_t i, const size_t j, const size_t k) 
{
   return _ny*_nz*i + _nz*j + k;
}

void get_row1(const size_t j, const size_t k, 
      const std::vector<double> &in,
      std::vector<double> &out
   ); 
void get_row2(const size_t i, const size_t k, 
      const std::vector<double> &in,
      std::vector<double> &out
   ); 
void get_row3(const size_t i, const size_t j, 
      const std::vector<double> &in,
      std::vector<double> &out
   ); 
void get_row12(const size_t k, 
      const std::vector<double> &in,
      std::vector<double> &out
   );
void get_row13(const size_t j, 
      const std::vector<double> &in,
      std::vector<double> &out
   );
void get_row23(const size_t i, 
      const std::vector<double> &in,
      std::vector<double> &out
   );

void set_row1(const size_t j, const size_t k, 
      const std::vector<double> &in,
      std::vector<double> &out
   ); 
void set_row2(const size_t i, const size_t k, 
      const std::vector<double> &in,
      std::vector<double> &out
   ); 
void set_row3(const size_t i, const size_t j, 
      const std::vector<double> &in,
      std::vector<double> &out
   ); 
void set_row12(const size_t k, 
      const std::vector<double> &in,
      std::vector<double> &out
   );
void set_row13(const size_t j, 
      const std::vector<double> &in,
      std::vector<double> &out
   );
void set_row23(const size_t i, 
      const std::vector<double> &in,
      std::vector<double> &out
   );
/*=========================================================================*/
}; /* Arr3d */
/*=========================================================================*/
#endif /* _ARR_HPP_ */
