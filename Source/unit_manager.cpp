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
inline double _average(const double a, const double b)
{
   return (a+b)/2.0;
}
/*===========================================================================*/
/* At the boundary of the grid cells we need to set boundary
 * conditions for the fields */
/*===========================================================================*/
void _partial_t_p( 
      const Grid &grid_L,
      const Grid &grid_R,
      const std::vector<double> &p_L,
      const std::vector<double> &p_R,
      const std::vector<double> &cs_M_L,
      const std::vector<double> &cs_P_R,
      std::vector<double> &dt_p_L,
      std::vector<double> &dt_p_R
      )
{
   assert(p_L.size()==grid_L.nx()*_nphi*_nlat);
   assert(p_R.size()==grid_R.nx()*_nphi*_nlat);

   assert(cs_M_L.size()==_nphi*_nlat);
   assert(cs_P_R.size()==_nphi*_nlat);

   assert(dt_p_L.size()==grid_L.nx()*_nphi*_nlat);
   assert(dt_p_R.size()==grid_R.nx()*_nphi*_nlat);

   std::vector<double> dr_p_L(p_L.size(),0);
   std::vector<double> dr_p_R(p_R.size(),0);

   grid_L.set_partial_r(p_L, dr_p_L);
   grid_R.set_partial_r(p_R, dr_p_R);

   const size_t ix_L = 0;
   const size_t ix_R = grid_R.nx()-1;

   for (size_t ip=0; ip<_nphi; ip++) {
   for (size_t it=0; it<_nlat; it++) {
      const size_t indx_L = grid_L.indx_r_th_ph(ix_L, it, ip);
      const size_t indx_R = grid_R.indx_r_th_ph(ix_R, it, ip);

      const size_t indx_S = grid_R.indx_th_ph(it, ip);

      const double cs_M = cs_M_L[indx_S];
      const double cs_P = cs_P_R[indx_S];

      const double dt_psi_P = dt_p_L[indx_L] + cs_M*dr_p_L[indx_L];
      const double dt_psi_M = dt_p_R[indx_R] + cs_P*dr_p_R[indx_R];

      const double dt_p = (
            cs_P*dt_psi_P - cs_M*dt_psi_M
         )/(
            cs_P - cs_M
         );
/*      std::cout
         <<std::setw(16)<<cs_M
         <<std::setw(16)<<cs_P
         <<std::setw(16)<<dt_p_L[indx_L] 
         <<std::setw(16)<<dt_p_R[indx_R] 
         <<std::setw(16)<<dt_p 
         <<std::endl;
*/
      dt_p_L[indx_L] = dt_p;
      dt_p_R[indx_R] = dt_p;
   }
   }
}
/*===========================================================================*/
void _set_boundaries(const int level, Unit *left, Unit *right)
{
   assert(left !=nullptr);
   assert(right!=nullptr);

   if (level==1) { 
      _partial_t_p( 
            left->grid,              right->grid,
            left->p.n,               right->p.n,
            left->scalar_eom.cs_M_R, right->scalar_eom.cs_P_L,
            left->p.k1,              right->p.k1
            );
   } else
   if (level==2) { 
      _partial_t_p( 
            left->grid,              right->grid,
            left->p.l2,              right->p.l2,
            left->scalar_eom.cs_M_R, right->scalar_eom.cs_P_L,
            left->p.k2,              right->p.k2
            );
   } else
   if (level==3) { 
      _partial_t_p( 
            left->grid,              right->grid,
            left->p.l3,              right->p.l3,
            left->scalar_eom.cs_M_R, right->scalar_eom.cs_P_L,
            left->p.k3,              right->p.k3
            );
   } else
   if (level==4) { 
      _partial_t_p( 
            left->grid,              right->grid,
            left->p.l4,              right->p.l4,
            left->scalar_eom.cs_M_R, right->scalar_eom.cs_P_L,
            left->p.k4,              right->p.k4
            );
   } else {
      std::cout<<"ERROR(_set_boundaries): level = "<<level<<std::endl;
   }
}
/*===========================================================================*/
/* RK4 time integration */
/*===========================================================================*/
void time_step(/*const size_t itm*/)
{
//   #pragma omp parallel for
   for (Unit *u: units) {
      u->grid.filter(u->f.n);
      u->grid.filter(u->p.n);
   }
   /*--------------------------------------*/
//   #pragma omp parallel for
   for (Unit *u: units) {
      u->set_k(1);
   }
   for (size_t i=0; i<_ngrids-1; i++) {
      _set_boundaries(1, units[i], units[i+1]);
   }
   /*--------------------------------------*/
//   #pragma omp parallel for
   for (Unit *u: units) {
      u->set_level(2);
      u->set_k(2);
   }
   for (size_t i=0; i<_ngrids-1; i++) {
      _set_boundaries(2, units[i], units[i+1]);
   }
   /*--------------------------------------*/
//   #pragma omp parallel for
   for (Unit *u: units) {
      u->set_level(3);
      u->set_k(3);
   }
   for (size_t i=0; i<_ngrids-1; i++) {
      _set_boundaries(3, units[i], units[i+1]);
   }
   /*--------------------------------------*/
//   #pragma omp parallel for
   for (Unit *u: units) {
      u->set_level(4);
      u->set_k(4);
   }
   for (size_t i=0; i<_ngrids-1; i++) {
      _set_boundaries(4, units[i], units[i+1]);
   }
   /*--------------------------------------*/
//   #pragma omp parallel for
   for (Unit *u: units) {
      u->set_level(5);
   }
/*
   for (size_t i=0; i<_ngrids-1; i++) {
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
