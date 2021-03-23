#include <cassert>
#include <cmath>
#include <iostream>
#include <string>

#include "params.hpp"
#include "cheb.hpp"
#include "sphere.hpp"
#include "arr.hpp"
#include "initial_data.hpp"

int main(int argc, char **argv)
{
   assert(argc==2);
   const std::string output_dir= argv[1];

   Params::init(output_dir);

   Arr3d::init( Params::nx(), Params::nlat(), Params::nphi());
   Sphere::init(Params::nl(), Params::nlat(), Params::nphi());
   Cheb::init(  Params::nx(), Params::rbl(),  Params::rbu());

   std::vector<double> f = Arr3d::arr3d(0.0);
   std::vector<double> p = Arr3d::arr3d(0.0);
   std::vector<double> q = Arr3d::arr3d(0.0);
   /* 
    * initial data 
    */
   ID::ingoing_pulse(f,p,q);

   Cheb::  cleanup();
   Sphere::cleanup();
   Arr3d:: cleanup();
 
   return EXIT_SUCCESS;
}
