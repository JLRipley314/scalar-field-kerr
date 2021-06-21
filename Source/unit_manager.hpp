#ifndef UNIT_MANAGER__
#define UNIT_MANAGER__

#include <vector>
#include <string>

#include "params.hpp"
#include "grid.hpp"
#include "scalar_eom.hpp"
#include "field.hpp"

/*===========================================================================*/
class Unit
{
public:
   Unit(
      const Params &params,
      const string loc,
      const double Rmin,
      const double Rmax,
      const size_t nx
      );
   ~Unit();

   void set_k(const int level);
   void set_level(const int level);
   void shift();

   const string loc;

   const size_t local_nx;

   Field f;
   Field p;

   std::vector<double> ghost_cells_left;
   std::vector<double> ghost_cells_right;
   std::vector<double> rho;

   Grid grid;
   Scalar_eom scalar_eom;
};
/*===========================================================================*/
/* Holds all evolution units and evolves them.
 * We use a namespace so we do not have to worry about
 * default constructors for the Unit classes, etc. */
/*===========================================================================*/
namespace Unit_manager
{
   void init(const Params &params);
   void cleanup();

   extern std::vector<Unit*> units;

   void time_step(/*const size_t itm*/);
   void shift();

   void set_initial_data(const Params &params);

   void write_time_stdout(const double time, const Params &params);

   void write_to_file(
         const size_t save_indx, 
         const std::string output_dir, 
         const Params &params
      );
}
#endif /* _UNIT_MANAGER_ */
