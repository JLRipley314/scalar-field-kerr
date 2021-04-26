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

   Unit_manager::init(params);
   Unit_manager::set_initial_data(params);
   Unit_manager::write_time_stdout(0, params);
   Unit_manager::write_to_file(save_indx, output_dir, params); 

   std::cout<<"Beginning evolution"<<std::endl;

   for (size_t itm=1; itm<params.nt(); itm++) {
      Unit_manager::time_step();
      /* 
       * save to file 
       * */
      if (itm%params.t_step_save()==0) {
         save_indx += 1;
         Unit_manager::write_time_stdout(itm*params.dt(), params);
         Unit_manager::write_to_file(save_indx, output_dir, params); 
      }
      Unit_manager::shift();
   } 
   Unit_manager::cleanup();

   return EXIT_SUCCESS;
}
