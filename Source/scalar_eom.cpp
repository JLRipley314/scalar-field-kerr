#include <vector>

#include "scalar_eom.hpp"
#include "arr.hpp"
#include "cheb.hpp"
#include "sphere.hpp"
#include "params.hpp"
/*==========================================================================*/
namespace Eom {
/*==========================================================================*/
namespace {
   std::vector<double> a_ff;
   std::vector<double> a_fp;
   std::vector<double> a_fq;

   std::vector<double> a_pf;
   std::vector<double> a_pp;
   std::vector<double> a_pq;

   std::vector<double> a_qf;
   std::vector<double> a_qp;
   std::vector<double> a_qq;

   std::vector<double> b_ff;
   std::vector<double> b_fp;
   std::vector<double> b_fq;

   std::vector<double> b_pf;
   std::vector<double> b_pp;
   std::vector<double> b_pq;

   std::vector<double> b_qf;
   std::vector<double> b_qp;
   std::vector<double> b_qq;
}
/*==========================================================================*/
void init()
{
};
void cleanup()
{
};
/*==========================================================================*/
void set_der_r(const std::vector<double> &v, std::vector<double> dv)
{
   std::vector<double> inter(   Params::nx());
   std::vector<double> inter_dv(Params::nx());

   for (size_t ip=0; ip<Params::nphi(); ip++) {
   for (size_t it=0; it<Params::nlat(); it++) {
      Arr3d::get_row1(ip, it, v, inter); 

      Cheb::der(inter, inter_dv);

      Arr3d::set_row1(ip, it, inter_dv, dv); 
   }
   }
}
/*==========================================================================*/
void set_spherical_lap(const std::vector<double> &v, std::vector<double> ddv)
{
   std::vector<double> inter(    Params::nlat()*Params::nphi());
   std::vector<double> inter_ddv(Params::nlat()*Params::nphi());

   for (size_t ix=0; ix<Params::nx(); ix++) {
      Arr3d::get_row23(ix, v, inter); 

      Sphere::laplace_beltrami(inter, inter_ddv);

      Arr3d::set_row23(ix, inter_ddv, ddv); 
   }
}
/*==========================================================================*/
void set_partial_phi(const std::vector<double> &v, std::vector<double> dv)
{
   std::vector<double> inter(   Params::nlat()*Params::nphi());
   std::vector<double> inter_dv(Params::nlat()*Params::nphi());

   for (size_t ix=0; ix<Params::nx(); ix++) {
      Arr3d::get_row23(ix, v, inter); 

      Sphere::partial_phi(inter, inter_dv);

      Arr3d::set_row23(ix, inter_dv, dv); 
   }
}
/*==========================================================================*/
void time_step(Field &f, Field &p, Field &q)
{

}
/*==========================================================================*/
} /* Eom */
