/*
 * Field equations for scalar field about a Kerr black hole.
 */
#ifndef SCALAR_EOM_HPP__ 
#define SCALAR_EOM_HPP__

#include <vector>

#include "grid.hpp"
#include "params.hpp"
#include "field.hpp"
/*==========================================================================*/
class Scalar_eom 
{
public:
   Scalar_eom(  const Grid &grid,
         const Params &params
      );
   ~Scalar_eom();

   void time_step(const Grid &grid, Field &f, Field &p);

   void set_k(
         const Grid &grid,
         const std::vector<double> &f,
         const std::vector<double> &p,
         std::vector<double> &f_k,
         std::vector<double> &p_k
         );

   void set_level(
         const int level,
         Field &f,
         Field &p) const;

   /* local energy density */
   void set_rho(
         const Grid &grid, 
         const std::vector<double> &f,
         const std::vector<double> &p,
         std::vector<double> &rho
         );

private:
   const size_t _n;
   const double _dt;

   const double _km1;
   const double _k0;
   const double _k1;
   const double _k2;

   const double _v2;
   const double _v3;
   const double _v4;

   std::vector<double> _pre;

   std::vector<double> _p_p;
   std::vector<double> _p_dr_f;
   std::vector<double> _p_dr_p;
   std::vector<double> _p_dr_dr_f;
   std::vector<double> _p_dphi_dr_f;
   std::vector<double> _p_lap_f;

   std::vector<double> _p_p_p;
   std::vector<double> _p_p_dr_f;
   std::vector<double> _p_dr_f_dr_f;
   std::vector<double> _p_dr_f_dphi_f;
   std::vector<double> _p_sphereX_f;

   /* 
    * for computing energy density rho 
    */
   std::vector<double> _rho_vv;
   std::vector<double> _rho_vr;
   std::vector<double> _rho_rr;
   std::vector<double> _rho_rphi;
   std::vector<double> _rho_sphereX;

   /* 
    * for computing time evolution of the
    * scalar field 
    */
   std::vector<double> _dr_f;
   std::vector<double> _lap_f;
   std::vector<double> _dr_p;
   std::vector<double> _dr_dr_f;
   std::vector<double> _dphi_f;
   std::vector<double> _dphi_dr_f;
   std::vector<double> _sphereX_f;
};
#endif /* _SCALAR_EOM_HPP_ */
