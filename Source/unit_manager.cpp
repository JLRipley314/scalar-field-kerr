#include <cassert>
#include <iostream>
#include <iomanip>
#include "unit_manager.hpp"

#include "io.hpp"
#include "initial_data.hpp"
/*=========================================================================*/
Unit::Unit(
      const Params &params,
      const string loc,
      const double Rmin,
      const double Rmax,
      const size_t nx)
:  loc{loc},
   f("f", (6+nx)*params.nlat()*params.nphi(), 0.0), // including ghost cells
   p("p", (6+nx)*params.nlat()*params.nphi(), 0.0),
   rho(nx*params.nlat()*params.nphi(), 0.0),
   grid(
      params.cl(),
      Rmin,
      Rmax,
      nx,
      params.nl(),
      params.nm(),
      params.nlat(),
      params.nphi() 
   ),
   scalar_eom(grid, params)
{
}
/*=========================================================================*/
Unit::~Unit()
{
}
/*=========================================================================*/
void Unit::set_k(const int level)
{
   if (level==1) {
      scalar_eom.set_k(grid, f.n, p.n, f.k1, p.k1);
   } else
   if (level==2) {
      scalar_eom.set_k(grid, f.l2, p.l2, f.k2, p.k2);
   } else
   if (level==3) {
      scalar_eom.set_k(grid, f.l3, p.l3, f.k3, p.k3);
   } else
   if (level==4) {
      scalar_eom.set_k(grid, f.l4, p.l4, f.k4, p.k4);
   } else {
      std::cout<<"ERROR(Unit::set_k): level = "<<level<<std::endl;
   }
}
/*=========================================================================*/
void Unit::set_level(const int level)
{
   scalar_eom.set_level(level, f, p);
}
/*=========================================================================*/
void Unit::shift()
{
   f.shift();
   p.shift(); 
}
/*=========================================================================*/
/*=========================================================================*/
namespace Unit_manager
{
/*=========================================================================*/
namespace {
   size_t _ngrids;
   size_t _nlat;
   size_t _nphi;
}
std::vector<Unit*> units = {};
/*=========================================================================*/
void init(const Params &params)
{
   std::cout<<"Initializing Unit_manager"<<std::endl;
   _ngrids = params.ngrids();
   _nlat   = params.nlat();
   _nphi   = params.nphi();

   assert(_ngrids>0);
   assert(_nlat  >0);
   assert(_nphi  >0);

   for (size_t i=0; i<_ngrids; i++) {
      std::cout<<"Unit "<<i<<std::endl;

      Unit *unit = new Unit(
            params,
            params.Rvals(i),
            params.Rvals(i+1),
            params.nxs(i)
            );
      assert(unit!=nullptr);
      units.push_back(unit);
   }
   std::cout<<"Finished initializing Unit_manager"<<std::endl;
}
/*=========================================================================*/
void cleanup()
{
   std::cout<<"Cleanup Unit_manager"<<std::endl;
   for (Unit *u: units) {
      delete u;
      u = nullptr;
   }
   std::cout<<"Finished Cleanup Unit_manager"<<std::endl;
}
/*===========================================================================*/
void shift()
{
   for (Unit *u: units) {
      u->shift();
   }
}
/*===========================================================================*/
static void _fill_ghost_cells(const int level, Unit *left, Unit *right)
{
   assert(left !=nullptr);
   assert(right!=nullptr);

   for (size_t i=0; i<_nlat*_nphi; i++) {
      left->f.l[(_nlat*_nphi)*(_nx-1) + i] = right->f.l[(_nlat*_nphi)*(3) + i];
      left->f.l[(_nlat*_nphi)*(_nx-2) + i] = right->f.l[(_nlat*_nphi)*(4) + i];
      left->f.l[(_nlat*_nphi)*(_nx-3) + i] = right->f.l[(_nlat*_nphi)*(5) + i];

      right->f.l[(_nlat*_nphi)*(0) + i] = left->f.l[(_nlat*_nphi)*(_nx-7) + i];
      right->f.l[(_nlat*_nphi)*(1) + i] = left->f.l[(_nlat*_nphi)*(_nx-6) + i];
      right->f.l[(_nlat*_nphi)*(2) + i] = left->f.l[(_nlat*_nphi)*(_nx-5) + i];
   }
}
/*===========================================================================*/
/* RK4 time integration */
/*===========================================================================*/
void time_step(/*const size_t itm*/)
{
//   #pragma omp parallel for
   for (Unit *u: units) {
      #pragma omp parallel sections
      {
         #pragma omp section
         {
            u->grid.filter(u->f.n);
         }
         #pragma omp section
         {
            u->grid.filter(u->p.n);
         }
      }
   }
   /*--------------------------------------*/
//   #pragma omp parallel for
   for (Unit *u: units) {
      u->set_k(1);
   }
   for (size_t i=0; i<_ngrids-1; i++) {
      _fill_ghost_cells(1, units[i], units[i+1]);
   }
   /*--------------------------------------*/
//   #pragma omp parallel for
   for (Unit *u: units) {
      u->set_level(2);
      u->set_k(2);
   }
   for (size_t i=0; i<_ngrids-1; i++) {
      _fill_ghost_cells(2, units[i], units[i+1]);
   }
   /*--------------------------------------*/
//   #pragma omp parallel for
   for (Unit *u: units) {
      u->set_level(3);
      u->set_k(3);
   }
   for (size_t i=0; i<_ngrids-1; i++) {
      _fill_ghost_cells(3, units[i], units[i+1]);
   }
   /*--------------------------------------*/
//   #pragma omp parallel for
   for (Unit *u: units) {
      u->set_level(4);
      u->set_k(4);
   }
   for (size_t i=0; i<_ngrids-1; i++) {
      _fill_ghost_cells(4, units[i], units[i+1]);
   }
   /*--------------------------------------*/
//   #pragma omp parallel for
   for (Unit *u: units) {
      u->set_level(5);
   }
}
/*===========================================================================*/
/* Sets initial data*/ 
/*===========================================================================*/
void set_initial_data(const Params &params)
{
   for (Unit *u: units) {
      ID::compact_pulse(params, u->grid, u->f.n,   u->p.n);
      ID::compact_pulse(params, u->grid, u->f.np1, u->p.np1);
   }
}
/*===========================================================================*/
void write_time_stdout(const double time, const Params &params)
{
   std::cout
      <<std::setw(16)<<time;
   size_t i=0;

   for (Unit *u: units) {
      i += 1;
      const double res = u->grid.norm_indep_res(params.dt(), u->f.n, u->f.np1, u->p.n);
      const double tv  = u->grid.total_variation(u->f.n);
      std::cout
         <<std::setw(16)<<"Grid "<<i
         <<std::setw(16)<<res
         <<std::setw(16)<<tv;
   }
   std::cout<<std::endl;
}
/*===========================================================================*/
/* Saves all the unit fields to file */ 
/*===========================================================================*/
void write_to_file(
      const size_t save_indx, 
      const std::string output_dir, 
      const Params &params)
{
   const double save_rmin = params.Rmin();
   const double save_rmax = params.Rmax();
   const double val_min = 1e-4;
   const double val_max = 1e4;

   bool save_coords = true;

   size_t indx=1;
   for (Unit *u: units) {
      u->scalar_eom.set_rho(u->grid, u->f.np1, u->p.np1, u->rho);

      Csv::write_x_y_z(u->grid, output_dir+"/"+u->f.name, save_coords, save_indx, val_min, val_max, save_rmin, save_rmax, u->f.np1);
      Csv::write_x_y_z(u->grid, output_dir+"/rho",        save_coords, save_indx, val_min, val_max, save_rmin, save_rmax, u->rho);

      Csv::write_R_psl(u->grid, output_dir+"/"+u->f.name, save_coords, save_indx, u->f.np1);
      Csv::write_R_psl(u->grid, output_dir+"/rho",        save_coords, save_indx, u->rho);
      indx += 1;

      save_coords = false;
   }
}
/*===========================================================================*/
}
