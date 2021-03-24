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
void get_row1(const size_t j, const size_t k, 
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
void get_row2(const size_t i, const size_t k, 
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
void get_row3(const size_t i, const size_t j, 
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
void get_row12(const size_t k, 
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
void get_row13(const size_t j, 
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
void get_row23(const size_t i, 
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
void set_row1(const size_t j, const size_t k, 
      const std::vector<double> &in,
      std::vector<double> &out) 
{
   assert(in.size() ==nx_);
   assert(out.size()==nx_*ny_*nz_);
   for (size_t i=0; i<nx_; i++) {
      out[INDX3D(ny_,nz_,i,j,k)] = in[i];
   }
} 
/*=========================================================================*/
void set_row2(const size_t i, const size_t k, 
      const std::vector<double> &in,
      std::vector<double> &out) 
{
   assert(in.size() ==ny_);
   assert(out.size()==nx_*ny_*nz_);
   for (size_t j=0; j<ny_; j++) {
      out[INDX3D(ny_,nz_,i,j,k)] = in[j];
   }
} 
/*=========================================================================*/
void set_row3(const size_t i, const size_t j, 
      const std::vector<double> &in,
      std::vector<double> &out) 
{
   assert(in.size() ==nz_);
   assert(out.size()==nx_*ny_*nz_);
   for (size_t k=0; k<nz_; k++) {
      out[INDX3D(ny_,nz_,i,j,k)] = in[k];
   }
} 
/*=========================================================================*/
void set_row12(const size_t k, 
      const std::vector<double> &in,
      std::vector<double> &out)
{
   assert(in.size() ==nx_*ny_);
   assert(out.size()==nx_*ny_*nz_);
   for (size_t i=0; i<nx_; i++) {
   for (size_t j=0; j<ny_; j++) {
      out[INDX3D(ny_,nz_,i,j,k)] = in[INDX2D(ny_,i,j)];
   }
   }
}
/*=========================================================================*/
void set_row13(const size_t j, 
      const std::vector<double> &in,
      std::vector<double> &out)
{
   assert(in.size() ==nx_*nz_);
   assert(out.size()==nx_*ny_*nz_);
   for (size_t i=0; i<nx_; i++) {
   for (size_t k=0; k<nz_; k++) {
      out[INDX3D(ny_,nz_,i,j,k)] = in[INDX2D(nz_,i,k)];
   }
   }
}
/*=========================================================================*/
void set_row23(const size_t i, 
      const std::vector<double> &in,
      std::vector<double> &out)
{
   assert(in.size() ==ny_*nz_);
   assert(out.size()==nx_*ny_*nz_);
   for (size_t j=0; j<ny_; j++) {
   for (size_t k=0; k<nz_; k++) {
      out[INDX3D(ny_,nz_,i,j,k)] = in[INDX2D(nz_,j,k)];
   }
   }
}
/*=========================================================================*/
} /* Arr3d */
#undef INDX2D
#undef INDX3D
