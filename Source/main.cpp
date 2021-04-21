#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>
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
   Sphere::init(Params::nl(), Params::nm(), Params::nlat(), Params::nphi());

   std::cout<<"Initializing Cheb"<<std::endl;
   Cheb::init(Params::nx(), Params::Rmin(), Params::Rmax());

   std::cout<<"Initializing Grid"<<std::endl;
   Grid::init(Params::cl(), Params::nx(), Params::nlat(), Params::nphi(), Params::nl());

   std::cout<<"Initializing Equations of motion"<<std::endl;
   Eom::init();

   std::cout<<"Initializing IO routines"<<std::endl;
   Csv::init();

   std::cout<<"Initializing Fields"<<std::endl;
   Field f("f", Params::nx_nlat_nphi(), 0.0);
   Field p("p", Params::nx_nlat_nphi(), 0.0);

   std::vector<double> rho(Params::nx_nlat_nphi(), 0.0);

   std::cout<<"Setting initial data"<<std::endl;

   ID::ingoing_pulse(f.n, p.n);

   Eom::set_rho(f.n, p.n, rho);

   size_t save_indx = 0;

//   const double save_rmin = Params::Rmin();
//   const double save_rmax = Params::Rmin() + ((Params::Rmax() - Params::Rmin())/4.0);

//   Csv::write_x_y_z(output_dir+"/"+f.name, save_indx, 1e-8, 1e6, save_rmin, save_rmax, f.n);
//   Csv::write_x_y_z(output_dir+"/rho",     save_indx, 1e-8, 1e6, save_rmin, save_rmax, rho);

   Csv::write_th_ph(output_dir+"/horizon_"+f.name, save_indx, Params::nx()-1, f.n);
   Csv::write_th_ph(output_dir+"/horizon_rho",     save_indx, Params::nx()-1, rho);

   Csv::write_x_z(output_dir+"/"+f.name, save_indx, 0, f.n);
   Csv::write_x_z(output_dir+"/rho",     save_indx, 0, rho);

   Csv::write_n_psl(output_dir+"/"+f.name, save_indx, f.n);
   Csv::write_n_psl(output_dir+"/rho",     save_indx, rho);

   const double res = Grid::norm_indep_res(Params::dt(), f.n, f.np1, p.n);
   const double tv  = Grid::total_variation(f.n);
   std::cout
      <<std::setw(16)<<0
      <<std::setw(16)<<res
      <<std::setw(16)<<tv
      <<std::endl;

   std::cout<<"Beginning evolution"<<std::endl;
   for (size_t itm=1; itm<Params::nt(); itm++) {

      Eom::time_step(f, p);

      /* save to file */
      if (itm%Params::t_step_save()==0) {
         const double res = Grid::norm_indep_res(Params::dt(), f.n, f.np1, p.np1);
         const double tv  = Grid::total_variation(f.np1); 
         std::cout
            <<std::setw(16)<<itm*Params::dt()/Params::bh_mass()
            <<std::setw(16)<<res
            <<std::setw(16)<<tv
            <<std::endl;
         if (std::isnan(res)) return EXIT_SUCCESS;

         Eom::set_rho(f.np1, p.np1, rho);

         save_indx += 1;

//         Csv::write_x_y_z(output_dir+"/"+f.name, save_indx, 1e-8, 1e6, save_rmin, save_rmax, f.np1);
//         Csv::write_x_y_z(output_dir+"/rho",     save_indx, 1e-8, 1e6, save_rmin, save_rmax, rho);

         Csv::write_th_ph(output_dir+"/horizon_"+f.name, save_indx, Params::nx()-1, f.np1);
         Csv::write_th_ph(output_dir+"/horizon_rho",     save_indx, Params::nx()-1, rho);

         Csv::write_x_z(output_dir+"/"+f.name, save_indx, 0, f.np1);
         Csv::write_x_z(output_dir+"/rho",     save_indx, 0, rho);

         Csv::write_n_psl(output_dir+"/"+f.name, save_indx, f.np1);
         Csv::write_n_psl(output_dir+"/rho",     save_indx, rho);
      }
      f.shift();
      p.shift();
   }

   std::cout<<"Cleaning up"<<std::endl;
   Cheb::  cleanup();
   Sphere::cleanup();
 
   return EXIT_SUCCESS;
}
