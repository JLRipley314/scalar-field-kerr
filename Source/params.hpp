#ifndef _PARAMS_HPP_
#define _PARAMS_HPP_

#include <string>

/*===========================================================================*/
namespace Params
{
   size_t nt();
   size_t nx();
   size_t nlat();
   size_t nphi();
   int t_step_save();

   int direction();

   double dt();

   double curv();
   double cl();

   int initial_exc_i();

   double bh_mass();
/*---------------------------------------------------------------------------*/
/* for the potentials */
   double V_0();
   double V_1();
   double V_2();
   double V_3();
   double V_4();
/*---------------------------------------------------------------------------*/
/* for the initial data */
   std::string id();

   double amp();
   double r_l();
   double r_u();

   int l_ang();
   int m_ang();
   
   void init(const std::string output_dir);
}
#endif // _PARAMS_HPP_
