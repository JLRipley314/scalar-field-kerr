#include "field.hpp"

/*==========================================================================*/
Field::Field(
         std::string name,
         const size_t num,
         const double init_val
      )
:  name{name}
{
   n.resize(num,init_val);
   l2.resize(num,init_val);
   l3.resize(num,init_val);
   l4.resize(num,init_val);
   np1.resize(num,init_val);

   k1.resize(num,init_val);
   k2.resize(num,init_val);
   k3.resize(num,init_val);
   k4.resize(num,init_val);
   k5.resize(num,init_val);
}
/*==========================================================================*/
Field::~Field()
{
}
/*==========================================================================*/
