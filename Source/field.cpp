/*
 * Field class to store 3D data
 */
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
Field3d::~Field3d()
{
}
