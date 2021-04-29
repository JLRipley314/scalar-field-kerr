#include <cmath>
#include <cassert>
#include <iostream>
#include <iomanip>
#include "finite_diff.hpp"
/*==========================================================================*/
FD::FD(const size_t n, const double lower, const double upper)
:  _n{n},
   _lower{lower},
   _upper{upper},
   _dx{(_upper-_lower)/(_n-1)},
   _pts(_n,0)
{
   std::cout<<"Initializing Finite Differences"<<std::endl;
   /* 
    * Chebyshev points on interval [lower, upper]
    */
   for (size_t i=0; i<_n; i++) {
      _pts[i] = _dx*i + _lower;
   }
   std::cout<<"Finished initializing Finite Differences"<<std::endl;
}
/*==========================================================================*/
FD::~FD()
{
}
/*==========================================================================*/
/* Compute derivative over interval */
/*==========================================================================*/
void FD::der(std::vector<double> &v, std::vector<double> &dv) const
{
   assert( v.size()==_n);
   assert(dv.size()==_n);

   dv[0] = (
      -  (49./20.)*v[0]
      +  (6.     )*v[1]
      -  (15./2. )*v[2]
      +  (20./3. )*v[3]
      -  (15./4. )*v[4]
      +  (6./5.  )*v[5]
      -  (1./6.  )*v[6]
      )/(_dx)
      ;
   dv[1] = (
      -  (1./6.  )*v[0]
      -  (77./60.)*v[1]   
      +  (5./2.  )*v[2]   
      -  (5./3.  )*v[3]   
      +  (5./6.  )*v[4]   
      -  (1./4.  )*v[5]   
      +  (1./30. )*v[6]   
      )/(_dx)
      ;
   dv[2] = (
         (1./30.)*v[0]
      -  (2./5. )*v[1]   
      -  (7./12.)*v[2]   
      +  (4./3. )*v[3]   
      -  (1./2. )*v[4]   
      +  (2./15.)*v[5]   
      -  (1./60.)*v[6]   
      )/(_dx)
      ;

   dv[_n-1] = -(
      -  (49./20.)*v[_n-1-0]
      +  (6.     )*v[_n-1-1]
      -  (15./2. )*v[_n-1-2]
      +  (20./3. )*v[_n-1-3]
      -  (15./4. )*v[_n-1-4]
      +  (6./5.  )*v[_n-1-5]
      -  (1./6.  )*v[_n-1-6]
      )/(_dx)
      ;
   dv[_n-2] = -(
      -  (1./6.  )*v[_n-1-0]
      -  (77./60.)*v[_n-1-1]   
      +  (5./2.  )*v[_n-1-2]   
      -  (5./3.  )*v[_n-1-3]   
      +  (5./6.  )*v[_n-1-4]   
      -  (1./4.  )*v[_n-1-5]   
      +  (1./30. )*v[_n-1-6]   
      )/(_dx)
      ;
   dv[_n-3] = -(
         (1./30.)*v[_n-1-0]
      -  (2./5. )*v[_n-1-1]   
      -  (7./12.)*v[_n-1-2]   
      +  (4./3. )*v[_n-1-3]   
      -  (1./2. )*v[_n-1-4]   
      +  (2./15.)*v[_n-1-5]   
      -  (1./60.)*v[_n-1-6]   
      )/(_dx)
      ;
   for (size_t i=3; i<_n-3; i++) {
      dv[i] = (
         (1./60.)*v[i+3]
      -  (3./20.)*v[i+2]
      +  (3./4. )*v[i+1]
      -  (3./4. )*v[i-1]
      +  (3./20.)*v[i-2]
      -  (1./60.)*v[i-3] 
      )/_dx;       
   }
}
/*==========================================================================*/
void FD::der2(std::vector<double> &v, std::vector<double> &ddv) const
{
   assert(  v.size()==_n);
   assert(ddv.size()==_n);

   const double dx2 = pow(_dx,2);

   for (size_t i=3; i<_n-3; i++) {
      ddv[i] = (
         (1./90. )*v[i+3]
      -  (3./20. )*v[i+2]
      +  (3./2.  )*v[i+1]
      -  (49./18.)*v[i  ]
      +  (3./2.  )*v[i-1]
      -  (3./20. )*v[i-2]
      +  (1./90. )*v[i-3] 
      )/dx2;
   }
}
/*==========================================================================*/
/* Kreiss-Oliger filter */
/*==========================================================================*/
void FD::filter(std::vector<double> &v) const
{
   assert(v.size()==_n);

   const double eps = 0.5;

   std::vector<double> tmp(_n,0);

   for (size_t i=4; i<_n-4; i++) {
      tmp[i]=
                (eps/128.)*(  1.0) *v[i-4]
      +         (eps/128.)*( -8.0) *v[i-3]
      +         (eps/128.)*( 28.0) *v[i-2]
      +         (eps/128.)*(-56.0) *v[i-1]
      +  (1.0 + (eps/128.)*( 70.0))*v[i]
      +         (eps/128.)*(-56.0) *v[i+1]
      +         (eps/128.)*( 28.0) *v[i+2]
      +         (eps/128.)*( -8.0) *v[i+3]
      +         (eps/128.)*(  1.0) *v[i+4]
      ;
   }
   for (size_t i=0; i<_n; i++) { 
      v[i] = tmp[i];
   } 

}
/*==========================================================================*/
size_t FD::n() const { return _n; }
double FD::lower() const { return _lower; }
double FD::upper() const { return _upper; }
double FD::pt(const size_t i) const { return _pts[i]; }
/*==========================================================================*/


