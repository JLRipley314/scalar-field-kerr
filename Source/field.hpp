/*
 * contains all levels for time integration
 */
#ifndef _FIELD_HPP_
#define _FIELD_HPP_

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
      std::vector<double> l2;
      std::vector<double> l3;
      std::vector<double> l4;
      std::vector<double> np1;

      std::vector<double> k1;
      std::vector<double> k2;
      std::vector<double> k3;
      std::vector<double> k4;
      std::vector<double> k5;

   private:
      const size_t size;
};

#endif /* _FIELD_HPP_ */
