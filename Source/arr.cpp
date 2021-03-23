/*
 * Array classes 
 */
#include <cassert>
#include "arr.hpp"

#define INDX2D(n2,i1,i2) ((n2)*(i1) + (i2))
#define INDX3D(n2,n3,i1,i2,i3) ((n2)*(n3)*(i1) + (n3)*(i2) + (i3))
/*=========================================================================*/
namespace Arr3d {
/*=========================================================================*/
void init(const size_t nx, const size_t ny, const size_t nz)
{
   nx_ = nx;
   ny_ = ny;
   nz_ = nz;
}
/*=========================================================================*/
void cleanup()
{
};
/*=========================================================================*/
std::vector<double> arr3d(const double val)
{
   std::vector<double> out(nx_*ny_*nz_,val);
   return out;
}
/*=========================================================================*/
void row1(const size_t j, const size_t k, 
      const std::vector<double> &in,
      std::vector<double> &out) 
{
   assert(in.size() ==nx_*ny_*nz_);
   assert(out.size()==nx_);
   for (size_t i=0; i<nx_; i++) {
      out[i] = in[INDX3D(ny_,nz_,i,j,k)];
   }
} 
/*=========================================================================*/
void row2(const size_t i, const size_t k, 
      const std::vector<double> &in,
      std::vector<double> &out) 
{
   assert(in.size() ==nx_*ny_*nz_);
   assert(out.size()==ny_);
   for (size_t j=0; j<ny_; j++) {
      out[j] = in[INDX3D(ny_,nz_,i,j,k)];
   }
} 
/*=========================================================================*/
void row3(const size_t i, const size_t j, 
      const std::vector<double> &in,
      std::vector<double> &out) 
{
   assert(in.size() ==nx_*ny_*nz_);
   assert(out.size()==nz_);
   for (size_t k=0; k<nz_; k++) {
      out[k] = in[INDX3D(ny_,nz_,i,j,k)];
   }
} 
/*=========================================================================*/
void row12(const size_t k, 
      const std::vector<double> &in,
      std::vector<double> &out)
{
   assert(in.size() ==nx_*ny_*nz_);
   assert(out.size()==nx_*ny_);
   for (size_t i=0; i<nx_; i++) {
   for (size_t j=0; j<ny_; j++) {
      out[INDX2D(ny_,i,j)] = in[INDX3D(ny_,nz_,i,j,k)];
   }
   }
}
/*=========================================================================*/
void row13(const size_t j, 
      const std::vector<double> &in,
      std::vector<double> &out)
{
   assert(in.size() ==nx_*ny_*nz_);
   assert(out.size()==nx_*nz_);
   for (size_t i=0; i<nx_; i++) {
   for (size_t k=0; k<nz_; k++) {
      out[INDX2D(nz_,i,k)] = in[INDX3D(ny_,nz_,i,j,k)];
   }
   }
}
/*=========================================================================*/
void row23(const size_t i, 
      const std::vector<double> &in,
      std::vector<double> &out)
{
   assert(in.size() ==nx_*ny_*nz_);
   assert(out.size()==ny_*nz_);
   for (size_t j=0; j<ny_; j++) {
   for (size_t k=0; k<nz_; k++) {
      out[INDX2D(nz_,j,k)] = in[INDX3D(ny_,nz_,i,j,k)];
   }
   }
}
/*=========================================================================*/
} /* Arr3d */
#undef INDX2D
#undef INDX3D
