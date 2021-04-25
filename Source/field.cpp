#include "field.hpp"
#include <iostream>
/*==========================================================================*/
Field::Field(
         std::string name,
         const size_t size,
         const double init_val
      )
:  name{name},
   n(size,init_val),
   l2(size,init_val),
   l3(size,init_val),
   l4(size,init_val),
   np1(size,init_val),

   k1(size,init_val),
   k2(size,init_val),
   k3(size,init_val),
   k4(size,init_val),
   k5(size,init_val),
   size{size}
{
   std::cout<<"Initializing Field "<<name<<std::endl;
   std::cout<<"Finished initializing Field "<<name<<std::endl;
}
/*==========================================================================*/
Field::~Field()
{
}
/*==========================================================================*/
void Field::shift()
{
   for (size_t i=0; i<size; i++) {
      n[i] = np1[i];
      k1[i] = k5[i];
   }
}


