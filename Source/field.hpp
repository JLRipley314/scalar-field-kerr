/*
 * contains all levels for time integration
 */
#ifndef FIELD_HPP__
#define FIELD_HPP__

#include <string>
#include <vector>

class Field 
{
   public:
      Field(
            const std::string name,
            const size_t size,
            const double init_val
         );

      ~Field();

      void shift();

      const std::string name;

      std::vector<double> n;
      std::vector<double> l;
      std::vector<double> np1;

      std::vector<double> k;

      const size_t size;
};

#endif /* _FIELD_HPP_ */
