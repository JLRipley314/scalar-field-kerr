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
   
   const double lower = 0.0;
   const double upper = 1.0;

   Cheb::init(  nr, lower, upper);
   Sphere::init(nl, nlat, nphi);
   Arr3d::init( nr, nlat, nphi);

   std::vector<double> f = Arr3d::arr3d(1.0);
   std::vector<double> p = Arr3d::arr3d(1.0);
   std::vector<double> q = Arr3d::arr3d(1.0);

   ID::ingoing_pulse(f,p,q);

   Cheb::  cleanup();
   Sphere::cleanup();
   Arr3d:: cleanup();
 
   return EXIT_SUCCESS;
}
