#include <cassert>
#include <iostream>
#include <iomanip>
#include "unit_manager.hpp"

#include "io.hpp"
#include "initial_data.hpp"
/*=========================================================================*/
Unit::Unit(
      const Params &params,
      const double Rmin,
      const double Rmax,
      const size_t nx)
:  f("f", nx*params.nlat()*params.nphi(), 0.0),
   p("p", nx*params.nlat()*params.nphi(), 0.0),
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
inline double average(const double a, const double b)
{
   return (a+b)/2.0;
}
/*===========================================================================*/
/* At the boundary of the grid cells we need to set boundary
 * conditions for the fields */
/*===========================================================================*/
void _set_boundaries(const int level, Unit *left, Unit *right)
{
   assert(left !=nullptr);
   assert(right!=nullptr);

   const size_t nx = right->grid.nx();

   if (level==1) { 
      for (size_t ip=0; ip<_nphi; ip++) {
      for (size_t it=0; it<_nlat; it++) {
         const double p = average(
               left->p.k1[  left->grid.indx_r_th_ph(   0, it, ip)], 
            +  right->p.k1[right->grid.indx_r_th_ph(nx-1, it, ip)]
            );

         left->p.k1[  left->grid.indx_r_th_ph(   0, it, ip)] = p; 
         right->p.k1[right->grid.indx_r_th_ph(nx-1, it, ip)] = p;
      }
      }
   } else
   if (level==2) { 
      for (size_t ip=0; ip<_nphi; ip++) {
      for (size_t it=0; it<_nlat; it++) {
         const double p = average(
               left->p.k2[  left->grid.indx_r_th_ph(   0, it, ip)], 
            +  right->p.k2[right->grid.indx_r_th_ph(nx-1, it, ip)]
            );

         left->p.k2[  left->grid.indx_r_th_ph(   0, it, ip)] = p; 
         right->p.k2[right->grid.indx_r_th_ph(nx-1, it, ip)] = p;
      }
      }
   } else
   if (level==3) { 
      for (size_t ip=0; ip<_nphi; ip++) {
      for (size_t it=0; it<_nlat; it++) {
         const double p = average(
               left->p.k3[  left->grid.indx_r_th_ph(   0, it, ip)], 
            +  right->p.k3[right->grid.indx_r_th_ph(nx-1, it, ip)]
            );

         left->p.k3[  left->grid.indx_r_th_ph(   0, it, ip)] = p; 
         right->p.k3[right->grid.indx_r_th_ph(nx-1, it, ip)] = p;
      }
      }
   } else
   if (level==4) { 
      for (size_t ip=0; ip<_nphi; ip++) {
      for (size_t it=0; it<_nlat; it++) {
         const double p = average(
               left->p.k4[  left->grid.indx_r_th_ph(   0, it, ip)], 
            +  right->p.k4[right->grid.indx_r_th_ph(nx-1, it, ip)]
            );

         left->p.k4[  left->grid.indx_r_th_ph(   0, it, ip)] = p; 
         right->p.k4[right->grid.indx_r_th_ph(nx-1, it, ip)] = p;
      }
      }
   } else {
      std::cout<<"ERROR(_set_boundaries): level = "<<level<<std::endl;
   }
}
/*===========================================================================*/
/* RK4 time integration */
/*===========================================================================*/
void time_step(const size_t itm)
{
   #pragma omp parallel for
   for (Unit *u: units) {
      u->grid.filter(u->f.n);
      u->grid.filter(u->p.n);
   }
   /*--------------------------------------*/
   #pragma omp parallel for
   for (Unit *u: units) {
      u->set_k(1);
   }
   for (size_t i=0; i<_ngrids-1; i++) {
      _set_boundaries(1, units[i], units[i+1]);
   }
   /*--------------------------------------*/
   #pragma omp parallel for
   for (Unit *u: units) {
      u->set_level(2);
      u->set_k(2);
   }
   for (size_t i=0; i<_ngrids-1; i++) {
      _set_boundaries(2, units[i], units[i+1]);
   }
   /*--------------------------------------*/
   #pragma omp parallel for
   for (Unit *u: units) {
      u->set_level(3);
      u->set_k(3);
   }
   for (size_t i=0; i<_ngrids-1; i++) {
      _set_boundaries(3, units[i], units[i+1]);
   }
   /*--------------------------------------*/
   #pragma omp parallel for
   for (Unit *u: units) {
      u->set_level(4);
      u->set_k(4);
   }
   for (size_t i=0; i<_ngrids-1; i++) {
      _set_boundaries(4, units[i], units[i+1]);
   }
   /*--------------------------------------*/
   #pragma omp parallel for
   for (Unit *u: units) {
      u->set_level(5);
   }
/*   for (size_t i=0; i<_ngrids-1; i++) {
      const size_t nx = units[i+1]->grid.nx();
      for (size_t ip=0; ip<_nphi; ip++) {
      for (size_t it=0; it<_nlat; it++) {
         if (fabs(
            units[i  ]->p.np1[units[i  ]->grid.indx_r_th_ph(   0, it, ip)] 
         -  units[i+1]->p.np1[units[i+1]->grid.indx_r_th_ph(nx-1, it, ip)]
         )>1e-16 ||
         fabs(
            units[i  ]->f.np1[units[i  ]->grid.indx_r_th_ph(   0, it, ip)] 
         -  units[i+1]->f.np1[units[i+1]->grid.indx_r_th_ph(nx-1, it, ip)]
         )>1e-16) {
            std::cout
               <<itm<<"\t"
               <<ip<<"\t"
               <<it<<"\t" 
               <<
                  units[i  ]->p.np1[units[i  ]->grid.indx_r_th_ph(   0, it, ip)] 
               -  units[i+1]->p.np1[units[i+1]->grid.indx_r_th_ph(nx-1, it, ip)]
               <<"\t"
               <<
                  units[i  ]->f.np1[units[i  ]->grid.indx_r_th_ph(   0, it, ip)] 
               -  units[i+1]->f.np1[units[i+1]->grid.indx_r_th_ph(nx-1, it, ip)]
               <<std::endl;
            std::quick_exit(EXIT_FAILURE); 
         }
      }
      }
   }
*/
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

   for (Unit *u: units) {
      u->scalar_eom.set_rho(u->grid, u->f.np1, u->p.np1, u->rho);

      Csv::write_x_y_z(u->grid, output_dir+"/"+u->f.name, save_coords, save_indx, val_min, val_max, save_rmin, save_rmax, u->f.np1);
      Csv::write_x_y_z(u->grid, output_dir+"/rho",        save_coords, save_indx, val_min, val_max, save_rmin, save_rmax, u->rho);

      Csv::write_R_psl(u->grid, output_dir+"/"+u->f.name, save_coords, save_indx, u->f.np1);
      Csv::write_R_psl(u->grid, output_dir+"/rho",        save_coords, save_indx, u->rho);

      Csv::write_n_psl(u->grid, output_dir+"/"+u->f.name, save_coords, save_indx, u->f.np1);
      Csv::write_n_psl(u->grid, output_dir+"/rho",        save_coords, save_indx, u->rho);

      save_coords = false;
   }
}
/*===========================================================================*/
}
