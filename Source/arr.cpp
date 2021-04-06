/*
 * Array classes 
 */
#include <cassert>
#include "arr.hpp"

#define INDX2D(n2,i1,i2) ((n2)*(i1) + (i2))
#define INDX3D(n2,n3,i1,i2,i3) ((n2)*(n3)*(i1) + (n3)*(i2) + (i3))
/*=========================================================================*/
namespace Arr3d 
{
   namespace 
   {
      size_t _nx;
      size_t _ny;
      size_t _nz;
   }
/*=========================================================================*/
void init(const size_t nx, const size_t ny, const size_t nz)
{
   _nx = nx;
   _ny = ny;
   _nz = nz;
}
/*=========================================================================*/
void cleanup()
{
};
/*=========================================================================*/
size_t indx(const size_t i, const size_t j, const size_t k) 
{
   return _ny*_nz*i + _nz*j + k;
}
/*=========================================================================*/
std::vector<double> arr3d(const double val)
{
   std::vector<double> out(_nx*_ny*_nz,val);
   return out;
}
/*=========================================================================*/
void get_row1(const size_t j, const size_t k, 
      const std::vector<double> &in,
      std::vector<double> &out) 
{
   assert(in.size() ==_nx*_ny*_nz);
   assert(out.size()==_nx);
   for (size_t i=0; i<_nx; i++) {
      out[i] = in[INDX3D(_ny,_nz,i,j,k)];
   }
} 
/*=========================================================================*/
void get_row2(const size_t i, const size_t k, 
      const std::vector<double> &in,
      std::vector<double> &out) 
{
   assert(in.size() ==_nx*_ny*_nz);
   assert(out.size()==_ny);
   for (size_t j=0; j<_ny; j++) {
      out[j] = in[INDX3D(_ny,_nz,i,j,k)];
   }
} 
/*=========================================================================*/
void get_row3(const size_t i, const size_t j, 
      const std::vector<double> &in,
      std::vector<double> &out) 
{
   assert(in.size() ==_nx*_ny*_nz);
   assert(out.size()==_nz);
   for (size_t k=0; k<_nz; k++) {
      out[k] = in[INDX3D(_ny,_nz,i,j,k)];
   }
} 
/*=========================================================================*/
void get_row12(const size_t k, 
      const std::vector<double> &in,
      std::vector<double> &out)
{
   assert(in.size() ==_nx*_ny*_nz);
   assert(out.size()==_nx*_ny);
   for (size_t i=0; i<_nx; i++) {
   for (size_t j=0; j<_ny; j++) {
      out[INDX2D(_ny,i,j)] = in[INDX3D(_ny,_nz,i,j,k)];
   }
   }
}
/*=========================================================================*/
void get_row13(const size_t j, 
      const std::vector<double> &in,
      std::vector<double> &out)
{
   assert(in.size() ==_nx*_ny*_nz);
   assert(out.size()==_nx*_nz);
   for (size_t i=0; i<_nx; i++) {
   for (size_t k=0; k<_nz; k++) {
      out[INDX2D(_nz,i,k)] = in[INDX3D(_ny,_nz,i,j,k)];
   }
   }
}
/*=========================================================================*/
void get_row23(const size_t i, 
      const std::vector<double> &in,
      std::vector<double> &out)
{
   assert(in.size() ==_nx*_ny*_nz);
   assert(out.size()==_ny*_nz);
   for (size_t j=0; j<_ny; j++) {
   for (size_t k=0; k<_nz; k++) {
      out[INDX2D(_nz,j,k)] = in[INDX3D(_ny,_nz,i,j,k)];
   }
   }
}
/*=========================================================================*/
void set_row1(const size_t j, const size_t k, 
      const std::vector<double> &in,
      std::vector<double> &out) 
{
   assert(in.size() ==_nx);
   assert(out.size()==_nx*_ny*_nz);
   for (size_t i=0; i<_nx; i++) {
      out[INDX3D(_ny,_nz,i,j,k)] = in[i];
   }
} 
/*=========================================================================*/
void set_row2(const size_t i, const size_t k, 
      const std::vector<double> &in,
      std::vector<double> &out) 
{
   assert(in.size() ==_ny);
   assert(out.size()==_nx*_ny*_nz);
   for (size_t j=0; j<_ny; j++) {
      out[INDX3D(_ny,_nz,i,j,k)] = in[j];
   }
} 
/*=========================================================================*/
void set_row3(const size_t i, const size_t j, 
      const std::vector<double> &in,
      std::vector<double> &out) 
{
   assert(in.size() ==_nz);
   assert(out.size()==_nx*_ny*_nz);
   for (size_t k=0; k<_nz; k++) {
      out[INDX3D(_ny,_nz,i,j,k)] = in[k];
   }
} 
/*=========================================================================*/
void set_row12(const size_t k, 
      const std::vector<double> &in,
      std::vector<double> &out)
{
   assert(in.size() ==_nx*_ny);
   assert(out.size()==_nx*_ny*_nz);
   for (size_t i=0; i<_nx; i++) {
   for (size_t j=0; j<_ny; j++) {
      out[INDX3D(_ny,_nz,i,j,k)] = in[INDX2D(_ny,i,j)];
   }
   }
}
/*=========================================================================*/
void set_row13(const size_t j, 
      const std::vector<double> &in,
      std::vector<double> &out)
{
   assert(in.size() ==_nx*_nz);
   assert(out.size()==_nx*_ny*_nz);
   for (size_t i=0; i<_nx; i++) {
   for (size_t k=0; k<_nz; k++) {
      out[INDX3D(_ny,_nz,i,j,k)] = in[INDX2D(_nz,i,k)];
   }
   }
}
/*=========================================================================*/
void set_row23(const size_t i, 
      const std::vector<double> &in,
      std::vector<double> &out)
{
   assert(in.size() ==_ny*_nz);
   assert(out.size()==_nx*_ny*_nz);
   for (size_t j=0; j<_ny; j++) {
   for (size_t k=0; k<_nz; k++) {
      out[INDX3D(_ny,_nz,i,j,k)] = in[INDX2D(_nz,j,k)];
   }
   }
}
/*=========================================================================*/
} /* Arr3d */
#undef INDX2D
#undef INDX3D
