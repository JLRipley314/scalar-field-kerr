#include <cmath>
#include <iostream>

#include "cheb.hpp"
#include "sphere.hpp"
#include "arr.hpp"
#include "initial_data.hpp"

int main() {
   const size_t nr = pow(2,6);
   const size_t nl = pow(2,4);
   const size_t nlat = nl + 2;
   const size_t nphi = 3*nl;

   Cheb::init(nr, 0.0, 1.0);
   Sphere::init(nl, nlat, nphi);

   std::vector<double> r = Arr3d::arr3d(1.0);
   std::vector<double> f = Arr3d::arr3d(1.0);
   std::vector<double> p = Arr3d::arr3d(1.0);
   std::vector<double> q = Arr3d::arr3d(1.0);

   ID::ingoing_pulse(r,f,p,q);

   Cheb::cleanup();
   Sphere::cleanup();
 
   return EXIT_SUCCESS;
}
