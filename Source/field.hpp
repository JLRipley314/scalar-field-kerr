#ifndef _FIELD_HPP_
#define _FIELD_HPP_
/*
 * Field class to store 3D data
 */
#include <string>
#include <vector>

/*=========================================================================*/
class Field3d {
   public:
      std::string name;
      std::vector<double> vals;

      Field3d(const std::string name, 
            const size_t nx, 
            const size_t ny, 
            const size_t nz
         );
      ~Field3d();

      inline size_t indx(const size_t i, const size_t j, const size_t k) 
      {
         return ny*nz*i + nz*j + k;
      }

      void row1(const size_t j, const size_t k, std::vector<double> &v); 
      void row2(const size_t i, const size_t k, std::vector<double> &v); 
      void row3(const size_t i, const size_t j, std::vector<double> &v); 
      void row12(const size_t k, std::vector<double> &v); 
      void row13(const size_t j, std::vector<double> &v); 
      void row23(const size_t i, std::vector<double> &v); 

      const size_t nx;
      const size_t ny;
      const size_t nz;
};
/*=========================================================================*/
#endif /* _FIELD_HPP_ */
