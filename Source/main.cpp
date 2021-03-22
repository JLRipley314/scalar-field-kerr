#include <cmath>
#include <iostream>

#include "cheb.hpp"
#include "sphere.hpp"
#include "field.hpp"

int main() {
   const size_t nr = pow(2,6);
   const size_t nl = pow(2,4);

   Cheb::init(nr, 0.0, 1.0);
   Sphere::init(nl);

   Field3d f("scalar",nr,Sphere::nlat(),Sphere::nphi());

   Cheb::cleanup();
   Sphere::cleanup();
 
   return EXIT_SUCCESS;
}
