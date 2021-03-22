#include <cmath>
#include <iostream>

#include "cheb.hpp"
#include "field.hpp"

int main() {
   {
      const size_t nr = pow(2,6);
   const size_t nt = pow(2,5);
   const size_t np = pow(2,5);

   Cheb cheb(nr, 0.0, 1.0);

   Field3d f("scalar",nr,nt,np);
 
   return EXIT_SUCCESS;
   }
}
