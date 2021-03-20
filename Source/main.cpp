#include <cmath>
#include <iostream>

#include "field.hpp"

int main() {
   const size_t nr = pow(2,6);
   const size_t nt = pow(2,5);
   const size_t np = pow(2,5);

   Field f("scalar",nr,nt,np);
      
   return EXIT_SUCCESS;
}
