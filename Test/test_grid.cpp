#include "gtest/gtest.h"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <vector>

#include "../Source/cheb.hpp"
#include "../Source/sphere.hpp"
#include "../Source/grid.hpp"

/*==========================================================================*/
/* computes derivatives using Grid functions and analytically,
 * and takes one norm of difference of values. 
 */
std::vector<double> get_norm_diff(
      const size_t nx,
      const size_t nl,
      const size_t nlat,
      const size_t nphi,
      const double Rmin,
      const double Rmax)
{
   const size_t n = nx*nlat*nphi;

   const double cl = Rmax;

   const double rl    = 0.1*(Rmax-Rmin) + Rmin;
   const double ru    = 0.2*(Rmax-Rmin) + Rmin;
   const double width = ru-rl;

   Sphere::init(nl, nlat, nphi);
   Cheb::init(nx, Rmin, Rmax);
   Grid::init(cl, nx, nlat, nphi);

   std::vector<double> v(n,0);
   std::vector<double> dr_v1(n,0);
   std::vector<double> dphi_v1(n,0);
   std::vector<double> lap_v1(n,0);
   std::vector<double> dr_v2(n,0);
   std::vector<double> dphi_v2(n,0); 
   std::vector<double> lap_v2(n,0); 

   for (size_t i=0; i<n; i++) {
      const std::vector<double> loc_r = Grid::r_th_ph(i);
      const double r     = loc_r[0];
      const double theta = loc_r[1];
      const double phi   = loc_r[2];

      double bump = 0.0;

      if ((r<ru) && (r>rl)) {
         bump = exp(-1.0*width/(r-rl))*exp(-2.0*width/(ru-r));
      }

      double rval = pow((r-rl)/width,2)*pow((ru-r)/width,2)*bump;

      double der_rval = (
         (2.0*((  (r-rl)/width) )*pow(((ru-r)/width),2))
      -  (2.0*pow((r-rl)/width,2)*(    (ru-r)/width   ))
      +  (1.0*(1.0              )*pow(((ru-r)/width),2))
      -  (2.0*pow((r-rl)/width,2)*(1.0                ))
      )*bump/width;


      v[i] = rval*(
            1.0 
         +  pow(sin(theta),2)/(1.0 + pow(cos(phi),2))
         );
      dr_v1[i] = der_rval*(
            1.0 
         +  pow(sin(theta),2)/(1.0 + pow(cos(phi),2))
         );
      dphi_v1[i] = rval*(
            pow(sin(theta),2)*sin(2*phi)/pow(1.0+pow(cos(phi),2),2)
         );
      lap_v1[i] = rval*(
            31.0
         +  57.0*cos(2*theta)
         +  72.0*pow(cos(theta),2)*cos(2*phi)
         -  6.0*cos(4*phi)*pow(sin(theta),2)
         )/pow(3.0 + cos(2*phi),3)
         ;
   }
   Grid::set_partial_r(    v,   dr_v2);
   Grid::set_partial_phi(  v, dphi_v2);
   Grid::set_spherical_lap(v,  lap_v2);

   double norm_dr   = 0;
   double norm_dphi = 0;
   double norm_lap  = 0;
   for (size_t i=0; i<n; i++) {
      norm_dr   += fabs(dr_v1[i]   - dr_v2[i]);
      norm_dphi += fabs(dphi_v1[i] - dphi_v2[i]);
      norm_lap  += fabs(lap_v1[i]  - lap_v2[i]);
   }
   norm_dr   /= n;
   norm_dphi /= n;
   norm_lap  /= n;

   std::vector<double> norms(3,0);

   norms[0] = norm_dr;
   norms[1] = norm_dphi;
   norms[2] = norm_lap;

   return norms;
}
/*==========================================================================*/
/* Adds random noise to function, filters, and returns total variation
 * of both 
 */
std::vector<double> get_total_variation(
      const size_t nx,
      const size_t nl,
      const size_t nlat,
      const size_t nphi,
      const double Rmin,
      const double Rmax)
{
   const size_t n = nx*nlat*nphi;

   const double cl = Rmax;

   const double rl    = 0.1*(Rmax-Rmin) + Rmin;
   const double ru    = 0.5*(Rmax-Rmin) + Rmin;
   const double width = ru-rl;

   Sphere::init(nl, nlat, nphi);
   Cheb::init(nx, Rmin, Rmax);
   Grid::init(cl, nx, nlat, nphi);

   std::vector<double> v(n,0);
   std::vector<double> filter_v(n,0);

   for (size_t i=0; i<n; i++) {
      const std::vector<double> loc_r = Grid::r_th_ph(i);
      const double r     = loc_r[0];
      const double theta = loc_r[1];
      const double phi   = loc_r[2];

      double bump = 0.0;

      if ((r<ru) && (r>rl)) {
         bump = exp(-1.0*width/(r-rl))*exp(-2.0*width/(ru-r));
      }

      double rval = pow((r-rl)/width,2)*pow((ru-r)/width,2)*bump;

      v[i] = rval*(
            1.0 
         +  pow(sin(theta),2)/(1.0 + pow(cos(phi),2))
         );
      v[i] *= (1.0 + (0.1*rand()/RAND_MAX));
      filter_v[i] = v[i];
   }
   Grid::filter(filter_v);

   double tv_v        = 0;
   double tv_filter_v = 0;

   for (size_t ix=0; ix<nx-1;   ix++) {
   for (size_t it=0; it<nlat-1; it++) {
   for (size_t ip=0; ip<nphi-1; ip++) {
      tv_v += pow(
               pow(
                  v[Grid::indx(ix+1, it, ip)] 
               -  v[Grid::indx(ix,   it, ip)]
               ,2)
            +  pow(
                  v[Grid::indx(ix, it+1, ip)] 
               -  v[Grid::indx(ix, it,   ip)]
               ,2)
            +  pow(
                  v[Grid::indx(ix, it, ip+1)] 
               -  v[Grid::indx(ix, it, ip)]
               ,2)
            ,
            0.5); 

      tv_filter_v += pow(
               pow(
                  filter_v[Grid::indx(ix+1, it, ip)] 
               -  filter_v[Grid::indx(ix,   it, ip)]
               ,2)
            +  pow(
                  filter_v[Grid::indx(ix, it+1, ip)] 
               -  filter_v[Grid::indx(ix, it,   ip)]
               ,2)
            +  pow(
                  filter_v[Grid::indx(ix, it, ip+1)] 
               -  filter_v[Grid::indx(ix, it, ip)]
               ,2)
            ,
            0.5); 
   }
   }
   }
   std::vector<double> tvs(2,0);

   tvs[0] = tv_v;
   tvs[1] = tv_filter_v;

   return tvs;
}
/*==========================================================================*/
/* Testing application of cheb and sphere to 3d works 
 */
TEST(grid_test, test_dr_dphi_lap) {

   const size_t Rmin = 1;
   const size_t Rmax = 10;

   const size_t nx1   = 48;
   const size_t nl1   = 16;
   const size_t nlat1 = 2*nl1 + 2;
   const size_t nphi1 = nlat1;
   std::vector<double> norms1 = get_norm_diff(nx1, nl1, nlat1, nphi1, Rmin, Rmax);

   const size_t nx2   = 64;
   const size_t nl2   = 24;
   const size_t nlat2 = 2*nl2 + 2;
   const size_t nphi2 = nlat2;
   std::vector<double> norms2 = get_norm_diff(nx2, nl2, nlat2, nphi2, Rmin, Rmax);

   EXPECT_LT(norms2[0], norms1[0]);
   EXPECT_LT(norms2[1], norms1[1]);
   EXPECT_LT(norms2[2], norms1[2]);
}
/*==========================================================================*/
/* Testing grid filter is TVD 
 */
TEST(grid_test, test_tvd) {

   const size_t Rmin = 1;
   const size_t Rmax = 100;

   const size_t nx1   = 96;
   const size_t nl1   = 16;
   const size_t nlat1 = 2*nl1 + 2;
   const size_t nphi1 = nlat1;
   std::vector<double> tv1 = get_total_variation(nx1, nl1, nlat1, nphi1, Rmin, Rmax);

   const size_t nx2   = 128;
   const size_t nl2   = 24;
   const size_t nlat2 = 2*nl2 + 2;
   const size_t nphi2 = nlat2;
   std::vector<double> tv2 = get_total_variation(nx2, nl2, nlat2, nphi2, Rmin, Rmax);

   EXPECT_LT(tv1[1], tv1[0]);
   EXPECT_LT(tv2[1], tv2[0]);
}
