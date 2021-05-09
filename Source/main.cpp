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
#include "unit_manager.hpp"

int main(int argc, char **argv)
{
   assert(argc==3);
   const std::string param_file= argv[1];
   const std::string output_dir= argv[2];

   Params params(param_file);
   size_t save_indx = 0;

   Field f("f", params.nx()*params.nlat()*params.nphi(), 0.0);
   Field p("p", params.nx()*params.nlat()*params.nphi(), 0.0);
   std::vector<double> rho(params.nx()*params.nlat()*params.nphi(), 0.0);
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
   Scalar_eom scalar_eom(grid, params);

   ID::compact_pulse(params, grid, f.n,   p.n);
   ID::compact_pulse(params, grid, f.np1, p.np1);

   std::cout
      <<std::setw(16)<<0;
   size_t i=0;
   double res = grid.norm_indep_res(params.dt(), f.n, f.np1, p.n);
   double tv  = grid.total_variation(f.n);
   std::cout
      <<std::setw(16)<<"Grid "<<i
      <<std::setw(16)<<res
      <<std::setw(16)<<tv
      <<std::endl;

   const double save_rmin = params.Rmin();
   const double save_rmax = params.Rmax();
   const double val_min = 1e-4;
   const double val_max = 1e4;

   bool save_coords = true;
   scalar_eom.set_rho(grid, f.np1, p.np1, rho);

   Csv::write_x_y_z(grid, output_dir+"/"+f.name, save_coords, save_indx, val_min, val_max, save_rmin, save_rmax, f.np1);
   Csv::write_x_y_z(grid, output_dir+"/rho",     save_coords, save_indx, val_min, val_max, save_rmin, save_rmax, rho);

   Csv::write_R_psl(grid, output_dir+"/"+f.name, save_coords, save_indx, f.np1);
   Csv::write_R_psl(grid, output_dir+"/rho",     save_coords, save_indx, rho);

   std::cout<<"Beginning evolution"<<std::endl;

   for (size_t itm=1; itm<params.nt(); itm++) {
      scalar_eom.time_step(grid, f, p);
      /* 
       * save to file 
       * */
      if (itm%params.t_step_save()==0) {
         save_indx += 1;
         std::cout
            <<std::setw(16)<<params.dt()*itm;
         size_t i=0;
         res = grid.norm_indep_res(params.dt(), f.n, f.np1, p.n);
         tv  = grid.total_variation(f.np1);
         std::cout
            <<std::setw(16)<<"Grid "<<i
            <<std::setw(16)<<res
            <<std::setw(16)<<tv
            <<std::endl;

            scalar_eom.set_rho(grid, f.np1, p.np1, rho);

            Csv::write_x_y_z(grid, output_dir+"/"+f.name, save_coords, save_indx, val_min, val_max, save_rmin, save_rmax, f.np1);
            Csv::write_x_y_z(grid, output_dir+"/rho",     save_coords, save_indx, val_min, val_max, save_rmin, save_rmax, rho);

            Csv::write_R_psl(grid, output_dir+"/"+f.name, save_coords, save_indx, f.np1);
            Csv::write_R_psl(grid, output_dir+"/rho",     save_coords, save_indx, rho);
      }
   f.shift();
   p.shift();
   } 
   return EXIT_SUCCESS;
}
