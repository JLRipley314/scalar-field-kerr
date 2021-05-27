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
   l(size,init_val),
   np1(size,init_val),

   k(size,init_val),
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
#pragma omp parallel for
   for (size_t i=0; i<size; i++) {
      n[i] = np1[i];
   }
}


