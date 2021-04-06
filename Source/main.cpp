#include <cassert>
#include <cmath>
#include <iostream>
#include <string>

#include "params.hpp"
#include "cheb.hpp"
#include "sphere.hpp"
#include "arr.hpp"
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

   std::cout<<"Initializing Arr3d"<<std::endl;
   Arr3d::init(Params::nx(), Params::nphi(), Params::nlat());

   std::cout<<"Initializing Sphere"<<std::endl;
   Sphere::init(Params::nl(), Params::nphi(), Params::nlat());

   std::cout<<"Initializing Cheb"<<std::endl;
   Cheb::init(  Params::nx(), Params::Rmin(), Params::Rmax());

   std::cout<<"Initializing Grid"<<std::endl;
   Grid::init();
   /* For writing grid points to file 
    */
   std::vector<std::vector<double>> grid_cart(Params::nx_nphi_nlat(),std::vector<double>(3,0));
   std::vector<std::string> labels_cart = {"x", "y", "z", "value"};
   for (size_t i=0; i<Params::nx_nphi_nlat(); i++) {
      grid_cart[i] = Grid::pt_cart(i);
   }

   std::cout<<"Initializing Equations of motion"<<std::endl;
   Eom::init();

   std::cout<<"Initializing Fields"<<std::endl;
   Field f("f", Params::nx_nphi_nlat(), 0.0);
   Field p("p", Params::nx_nphi_nlat(), 0.0);
   Field q("q", Params::nx_nphi_nlat(), 0.0);

   std::cout<<"Setting initial data"<<std::endl;
   ID::ingoing_pulse(f.n, p.n, q.n);

   size_t save_indx = 0;
   Csv::write_unstructured(output_dir+"/"+f.name, save_indx, labels_cart, grid_cart, f.n);
   Csv::write_unstructured(output_dir+"/"+p.name, save_indx, labels_cart, grid_cart, p.n);
   Csv::write_unstructured(output_dir+"/"+q.name, save_indx, labels_cart, grid_cart, q.n);

   std::cout<<"Beginning evolution"<<std::endl;
   for (size_t itm=0; itm<Params::nt(); itm++) {

      Eom::time_step(f, p, q);

      /* save to file */
      if (itm%Params::t_step_save()==0) {
         save_indx += 1;
         std::cout<<itm*Params::dt()/Params::bh_mass()<<std::endl;
         Csv::write_unstructured(output_dir+"/"+f.name, save_indx, labels_cart, grid_cart, f.np1);
         Csv::write_unstructured(output_dir+"/"+p.name, save_indx, labels_cart, grid_cart, p.np1);
         Csv::write_unstructured(output_dir+"/"+q.name, save_indx, labels_cart, grid_cart, q.np1);
      }
      f.shift();
      p.shift();
      q.shift();
   }

   std::cout<<"Cleaning up"<<std::endl;
   Eom::   cleanup();
   Grid::  cleanup();
   Cheb::  cleanup();
   Sphere::cleanup();
   Arr3d:: cleanup();
 
   return EXIT_SUCCESS;
}
