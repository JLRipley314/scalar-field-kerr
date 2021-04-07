#include <cassert>
#include <cmath>
#include <iostream>
#include <string>

#include "params.hpp"
#include "cheb.hpp"
#include "sphere.hpp"
#include "grid.hpp"
#include "field.hpp"
#include "initial_data.hpp"
#include "scalar_eom.hpp"
#include "io.hpp"

int main(int argc, char **argv)
{
   assert(argc==3);
   const std::string param_file= argv[1];
   const std::string output_dir= argv[2];

   std::cout<<"Initializing Parameters"<<std::endl;
   Params::init(param_file);

   std::cout<<"Initializing Sphere"<<std::endl;
   Sphere::init(Params::nl(), Params::nlat(), Params::nphi());

   std::cout<<"Initializing Cheb"<<std::endl;
   Cheb::init(Params::nx(), Params::Rmin(), Params::Rmax());

   std::cout<<"Initializing Grid"<<std::endl;
   Grid::init();

   std::cout<<"Initializing Equations of motion"<<std::endl;
   Eom::init();

   std::cout<<"Initializing IO routines"<<std::endl;
   Csv::init();

   std::cout<<"Initializing Fields"<<std::endl;
   Field f("f", Params::nx_nlat_nphi(), 0.0);
   Field p("p", Params::nx_nlat_nphi(), 0.0);
   Field q("q", Params::nx_nlat_nphi(), 0.0);

   std::cout<<"Setting initial data"<<std::endl;
   ID::ingoing_pulse(f.n, p.n, q.n);

   size_t save_indx = 0;
   Csv::write_x_y_z(output_dir+"/"+f.name, save_indx, f.n);
   Csv::write_x_y_z(output_dir+"/"+p.name, save_indx, p.n);
   Csv::write_x_y_z(output_dir+"/"+q.name, save_indx, q.n);

   std::cout<<"Beginning evolution"<<std::endl;
   for (size_t itm=0; itm<Params::nt(); itm++) {

      Eom::time_step(f, p, q);

      /* save to file */
      if (itm%Params::t_step_save()==0) {
         std::cout<<itm*Params::dt()/Params::bh_mass()<<std::endl;

         save_indx += 1;
         Csv::write_x_y_z(output_dir+"/"+f.name, save_indx, f.np1);
         Csv::write_x_y_z(output_dir+"/"+p.name, save_indx, p.np1);
         Csv::write_x_y_z(output_dir+"/"+q.name, save_indx, q.np1);
      }
      f.shift();
      p.shift();
      q.shift();
   }

   std::cout<<"Cleaning up"<<std::endl;
   Cheb::  cleanup();
   Sphere::cleanup();
 
   return EXIT_SUCCESS;
}
