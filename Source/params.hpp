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
   size_t nlat();
   size_t nphi();
   size_t t_step_save();

   int direction();

   double dt();

   double curv();
   double cl();

   double rbl();
   double rbu();

   int initial_exc_i();

   double bh_mass();
   /* 
    * for the potentials 
    */
   double V_0();
   double V_1();
   double V_2();
   double V_3();
   double V_4();
   /* 
    * for the initial data 
    */
   std::string id();

   double amp();
   double r_l();
   double r_u();

   int l_ang();
   int m_ang();
}
#endif // _PARAMS_HPP_
