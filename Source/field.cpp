#include "field.hpp"

/*==========================================================================*/
Field::Field(
         std::string name,
         const size_t num,
         const double init_val
      )
:  name{name},
   n(num,init_val),
   l2(num,init_val),
   l3(num,init_val),
   l4(num,init_val),
   np1(num,init_val),

   k1(num,init_val),
   k2(num,init_val),
   k3(num,init_val),
   k4(num,init_val),
   k5(num,init_val)
{
}
/*==========================================================================*/
Field::~Field()
{
}
/*==========================================================================*/
