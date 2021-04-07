/*
 * utility routines for reading from parameter file
 * and accessing parameters in code
 */
#ifndef _PARAMS_HPP_
#define _PARAMS_HPP_

#include <string>

/*===========================================================================*/
namespace Params
{
  /* 
   * read in params from file
   */ 
   void init(const std::string output_dir);
   /*
    * parameters 
    */
   size_t nt();
   size_t nx();
   size_t nl();

   size_t nx_nlat_nphi();

   size_t nlat();
   size_t nphi();
   size_t t_step_save();

   double dt();

   double cl();

   double Rmax();
   double Rmin();

   double bh_mass();
   double bh_spin();
   /* 
    * for the potentials 
    */
   double V_2();
   double V_3();
   double V_4();
   /* 
    * for the initial data 
    */
   std::string id();

   std::string initial_data_direction();

   double amp();
   double rl();
   double ru();

   int l_ang();
   int m_ang();
}
#endif // _PARAMS_HPP_
