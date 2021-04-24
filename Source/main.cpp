#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>

#include "params.hpp"
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
   Params params(param_file);

   std::cout<<"Initializing Grid"<<std::endl;
   Grid grid(
      params.cl(),
      params.Rmin(),
      params.Rmax(),
      params.nx(),
      params.nl(),
      params.nm(),
      params.nlat(),
      params.nphi() 
      );

   std::cout<<"Initializing Equations of motion"<<std::endl;
   Eom eom(grid, params);

   std::cout<<"Initializing Fields"<<std::endl;
   Field f("f", params.nx_nlat_nphi(), 0.0);
   Field p("p", params.nx_nlat_nphi(), 0.0);

   std::vector<double> rho(params.nx_nlat_nphi(), 0.0);

   std::cout<<"Setting initial data"<<std::endl;

   ID::compact_pulse(
         grid,
         params,
         f.n, p.n
      );

   eom.set_rho(grid, f.n, p.n, rho);

   size_t save_indx = 0;

   const double save_rmin = params.Rmin();
   const double save_rmax = params.Rmin() + ((params.Rmax() - params.Rmin())/8.0);

   Csv::write_x_y_z(grid, output_dir+"/"+f.name, save_indx, 1e-8, 1e6, save_rmin, save_rmax, f.n);
   Csv::write_x_y_z(grid, output_dir+"/rho",     save_indx, 1e-8, 1e6, save_rmin, save_rmax, rho);

   Csv::write_R_psl(grid, output_dir+"/"+f.name, save_indx, f.n);
   Csv::write_R_psl(grid, output_dir+"/rho",     save_indx, rho);

   Csv::write_n_psl(grid, output_dir+"/"+f.name, save_indx, f.n);
   Csv::write_n_psl(grid, output_dir+"/rho",     save_indx, rho);

   const double res = grid.norm_indep_res(params.dt(), f.n, f.np1, p.n);
   const double tv  = grid.total_variation(f.n);
   std::cout
      <<std::setw(16)<<0
      <<std::setw(16)<<res
      <<std::setw(16)<<tv
      <<std::endl;

   std::cout<<"Beginning evolution"<<std::endl;
   for (size_t itm=1; itm<params.nt(); itm++) {

      eom.time_step(grid, f, p);

      /* save to file */
      if (itm%params.t_step_save()==0) {
         const double res = grid.norm_indep_res(params.dt(), f.n, f.np1, p.np1);
         const double tv  = grid.total_variation(f.np1); 
         std::cout
            <<std::setw(16)<<itm*params.dt()/params.bh_mass()
            <<std::setw(16)<<res
            <<std::setw(16)<<tv
            <<std::endl;
         if (std::isnan(res)) return EXIT_SUCCESS;

         eom.set_rho(grid, f.np1, p.np1, rho);

         save_indx += 1;

         Csv::write_x_y_z(grid, output_dir+"/"+f.name, save_indx, 1e-8, 1e6, save_rmin, save_rmax, f.np1);
         Csv::write_x_y_z(grid, output_dir+"/rho",     save_indx, 1e-8, 1e6, save_rmin, save_rmax, rho);

         Csv::write_R_psl(grid, output_dir+"/"+f.name, save_indx, f.np1);
         Csv::write_R_psl(grid, output_dir+"/rho",     save_indx, rho);

         Csv::write_n_psl(grid, output_dir+"/"+f.name, save_indx, f.np1);
         Csv::write_n_psl(grid, output_dir+"/rho",     save_indx, rho);
      }
      f.shift();
      p.shift();
   } 
   return EXIT_SUCCESS;
}
