#ifndef _FIELD_HPP_
#define _FIELD_HPP_
/*
 * Field class to store 3D data
 */
#include <string>
#include <vector>

/*=========================================================================*/
class Field {
   private:
      const size_t nx;
      const size_t ny;
      const size_t nz;
   public:
      std::string name;
      std::vector<double> vals;

      Field(const std::string name, 
            const size_t nx, 
            const size_t ny, 
            const size_t nz
         );
      ~Field();

      inline size_t indx(const size_t i, const size_t j, const size_t k) 
      {
         return ny*nz*(i) + nz*(j) + (k);
      }
};
/*=========================================================================*/
#endif /* _FIELD_HPP_ */
