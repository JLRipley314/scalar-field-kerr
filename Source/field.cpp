/*
 * Field class to store 3D data
 */
#include <cassert>
#include "field.hpp"

/*=========================================================================*/
Field3d::Field3d(
      const std::string name, 
      const size_t nx,
      const size_t ny,
      const size_t nz)
:  name{name},
   vals(nx*nx*nz,0),
   nx{nx},
   ny{ny},
   nz{nz}
{
}
/*=========================================================================*/
void Field3d::row1(const size_t j, const size_t k, std::vector<double> &v)
{
   assert(v.size()==nx);
   for (size_t i=0; i<nx; i++) {
      v[i] = vals[indx(i, j, k)];
   }
} 
/*=========================================================================*/
void Field3d::row2(const size_t i, const size_t k, std::vector<double> &v)
{
   assert(v.size()==ny);
   for (size_t j=0; j<ny; j++) {
      v[j] = vals[indx(i, j, k)];
   }
} 
/*=========================================================================*/
void Field3d::row3(const size_t i, const size_t j, std::vector<double> &v) 
{
   assert(v.size()==nz);
   for (size_t k=0; k<nz; k++) {
      v[k] = vals[indx(i, j, k)];
   }
} 
/*=========================================================================*/
void Field3d::row12(const size_t k, std::vector<double> &v) 
{
   assert(v.size()==nx*ny);
   for (size_t i=0; i<nx; i++) {
   for (size_t j=0; j<ny; j++) {
      v[ny*i + j] = vals[indx(i, j, k)];
   }
   }
}
/*=========================================================================*/
void Field3d::row13(const size_t j, std::vector<double> &v) 
{
   assert(v.size()==nx*nz);
   for (size_t i=0; i<nx; i++) {
   for (size_t k=0; k<nz; k++) {
      v[nz*i + k] = vals[indx(i, j, k)];
   }
   }
}
/*=========================================================================*/
void Field3d::row23(const size_t i, std::vector<double> &v) 
{
   assert(v.size()==ny*nz);
   for (size_t j=0; j<ny; j++) {
   for (size_t k=0; k<nz; k++) {
      v[nz*j + k] = indx(i, j, k);
   }
   }
}
/*=========================================================================*/
Field3d::~Field3d()
{
}
