#include <cassert>
#include <cmath>
#include <iostream>
#include <string>

#include "params.hpp"
#include "cheb.hpp"
#include "sphere.hpp"
#include "arr.hpp"
#include "field.hpp"
#include "initial_data.hpp"
#include "scalar_eom.hpp"

int main(int argc, char **argv)
{
   assert(argc==2);
   const std::string output_dir= argv[1];

   Params::init(output_dir);

   Arr3d::init( Params::nx(), Params::nphi(), Params::nlat());
   Sphere::init(Params::nl(), Params::nphi(), Params::nlat());
   Cheb::init(  Params::nx(), Params::rbl(),  Params::rbu());
   /* 
    * equations of motion 
    */
   Eom::init();

   Field f("f", Params::nx_nlat_nphi(), 0.0);
   Field p("p", Params::nx_nlat_nphi(), 0.0);
   Field q("q", Params::nx_nlat_nphi(), 0.0);
   /* 
    * initial data 
    */
   ID::ingoing_pulse(f.n, p.n, q.n);
   /* 
    * evolve in time 
    */
   for (size_t itm=0; itm<Params::nt(); itm++) {
      Eom::time_step(f, p, q);

      if (itm%Params::t_step_save()==0) {
         /* save to file */
      }
   }

   Cheb::  cleanup();
   Sphere::cleanup();
   Arr3d:: cleanup();
 
   return EXIT_SUCCESS;
}
