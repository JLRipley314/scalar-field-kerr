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
   assert(argc==3);
   const std::string param_file= argv[1];
   const std::string output_dir= argv[2];

   std::cout<<"Initializing Parameters"<<std::endl;
   Params::init(param_file);

   std::cout<<"Initializing Arr3d"<<std::endl;
   Arr3d::init( Params::nx(), Params::nphi(), Params::nlat());

   std::cout<<"Initializing Sphere"<<std::endl;
   Sphere::init(Params::nl(), Params::nphi(), Params::nlat());

   std::cout<<"Initializing Cheb"<<std::endl;
   Cheb::init(  Params::nx(), Params::rmin(), Params::rmax());

   std::cout<<"Initializing Equations of motion"<<std::endl;
   Eom::init();

   std::cout<<"Initializing Fields"<<std::endl;
   Field f("f", Params::nx_nlat_nphi(), 0.0);
   Field p("p", Params::nx_nlat_nphi(), 0.0);
   Field q("q", Params::nx_nlat_nphi(), 0.0);

   std::cout<<"Setting initial data"<<std::endl;
   ID::ingoing_pulse(f.n, p.n, q.n);

   std::cout<<"Beginning evolution"<<std::endl;
   for (size_t itm=0; itm<Params::nt(); itm++) {
      Eom::time_step(f, p, q);

      if (itm%Params::t_step_save()==0) {
         std::cout<<itm*Params::dt()/Params::bh_mass()<<std::endl;
         /* save to file */
      }
   }

   std::cout<<"Cleaning up"<<std::endl;
   Cheb::  cleanup();
   Sphere::cleanup();
   Arr3d:: cleanup();
 
   return EXIT_SUCCESS;
}
