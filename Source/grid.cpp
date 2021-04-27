#include <cassert>
#include <iostream>
#include <iomanip>

#include "grid.hpp"

/*===========================================================================*/
Grid::Grid(
      const double cl,
      const double Rmin,
      const double Rmax,
      const size_t nx,
      const size_t nl,
      const size_t nm,
      const size_t nlat,
      const size_t nphi)
:  _cl{cl},
   _Rmin{Rmin},
   _Rmax{Rmax},
   _nx{nx},
   _nl{nl},
   _nm{nm},
   _nlat{nlat},
   _nphi{nphi},
   _sphere(_nl, _nm, _nlat, _nphi),
   _cheb(_nx, _Rmin, _Rmax)
{
   std::cout<<"Initializing Grid"<<std::endl;

   _r_th_ph.resize(_nx*_nphi*_nlat, std::vector<double>(3,0));
   _R_th_ph.resize(_nx*_nphi*_nlat, std::vector<double>(3,0));
   _x_y_z.resize(  _nx*_nphi*_nlat, std::vector<double>(3,0));
   for (size_t ix=0; ix<_nx;   ix++) {
   for (size_t it=0; it<_nlat; it++) {
   for (size_t ip=0; ip<_nphi; ip++) {
      double R     = _cheb.pt(ix);
      double theta = _sphere.theta(it);  
      double phi   = _sphere.phi(  ip);

      _r_th_ph[indx_r_th_ph(ix,it,ip)][0] = (fabs(R-cl)>1e-16) ? (R/(1.0 - (R/cl))) : 1e12; 
      _r_th_ph[indx_r_th_ph(ix,it,ip)][1] = theta; 
      _r_th_ph[indx_r_th_ph(ix,it,ip)][2] = phi; 

      _R_th_ph[indx_r_th_ph(ix,it,ip)][0] = R;
      _R_th_ph[indx_r_th_ph(ix,it,ip)][1] = theta;
      _R_th_ph[indx_r_th_ph(ix,it,ip)][2] = phi;

      _x_y_z[indx_r_th_ph(ix,it,ip)][0] = R*cos(phi)*sin(theta); 
      _x_y_z[indx_r_th_ph(ix,it,ip)][1] = R*sin(phi)*sin(theta);
      _x_y_z[indx_r_th_ph(ix,it,ip)][2] = R*cos(theta);
   }
   }
   }
   _th_ph.resize(_nphi*_nlat, std::vector<double>(2,0));
   for (size_t it=0; it<_nlat; it++) {
   for (size_t ip=0; ip<_nphi; ip++) {
      double theta = _sphere.theta(it);  
      double phi   = _sphere.phi(  ip);
      _th_ph[indx_th_ph(it, ip)][0] = theta;
      _th_ph[indx_th_ph(it, ip)][1] = phi;
   }
   }
   _x_z.resize( _nx*_nlat, std::vector<double>(2,0));
   _R_th.resize(_nx*_nlat, std::vector<double>(2,0));
   for (size_t ix=0; ix<_nx;   ix++) {
   for (size_t it=0; it<_nlat; it++) {
      double R     = _cheb.pt(ix);
      double theta = _sphere.theta(it);  
      _R_th[indx_r_th(ix,it)][0] = R;
      _R_th[indx_r_th(ix,it)][1] = theta;
      _x_z[ indx_r_th(ix,it)][0] = R*sin(theta); 
      _x_z[ indx_r_th(ix,it)][1] = R*cos(theta);
   }
   }
   _n_l.resize(_nx*_nl, std::vector<double>(2,0));
   _R_l.resize(_nx*_nl, std::vector<double>(2,0));
   for (size_t ix=0; ix<_nx; ix++) {
   for (size_t il=0; il<_nl; il++) {
      double R = _cheb.pt(ix);
      double l = il;  
      _R_l[indx_r_l(ix,il)][0] = R;
      _R_l[indx_r_l(ix,il)][1] = l;
      _n_l[indx_r_l(ix,il)][0] = double(ix);
      _n_l[indx_r_l(ix,il)][1] = l;
   }
   }
   _partial_R_to_partial_r.resize(_nx,0); 
   for (size_t ix=0; ix<_nx; ix++) {
      _partial_R_to_partial_r[ix] = pow(1.0 - (_cheb.pt(ix)/cl),2);
   }
   std::cout<<"Finished initializing Grid"<<std::endl;
}
/*=========================================================================*/
Grid::~Grid()
{
}
/*=========================================================================*/
std::vector<double> Grid::compute_ylm(const size_t l_ang, const size_t m_ang) const
{
   return _sphere.compute_ylm(l_ang, m_ang);
}
/*=========================================================================*/
void Grid::get_row_R(const size_t it, const size_t ip, 
      const std::vector<double> &in,
      std::vector<double> &out) const 
{
   assert(in.size() ==_nx*_nphi*_nlat);
   assert(out.size()==_nx);
   for (size_t ix=0; ix<_nx; ix++) {
      out[ix] = in[indx_r_th_ph(ix,it,ip)];
   }
} 
void Grid::set_row_R(const size_t it, const size_t ip, 
      const std::vector<double> &in,
      std::vector<double> &out) const 
{
   assert(in.size() ==_nx);
   assert(out.size()==_nx*_nphi*_nlat);
   for (size_t ix=0; ix<_nx; ix++) {
      out[indx_r_th_ph(ix,it,ip)] = in[ix];
   }
} 
/*=========================================================================*/
void Grid::get_row_th(const size_t ix, const size_t ip, 
      const std::vector<double> &in,
      std::vector<double> &out) const 
{
   assert(in.size() ==_nx*_nphi*_nlat);
   assert(out.size()==_nlat);
   for (size_t it=0; it<_nlat; it++) {
      out[it] = in[indx_r_th_ph(ix,it,ip)];
   }
} 
void Grid::set_row_th(const size_t ix, const size_t ip, 
      const std::vector<double> &in,
      std::vector<double> &out) const 
{
   assert(in.size() ==_nlat);
   assert(out.size()==_nx*_nphi*_nlat);
   for (size_t it=0; it<_nlat; it++) {
      out[indx_r_th_ph(ix,it,ip)] = in[it];
   }
} 
/*=========================================================================*/
void Grid::get_row_ph(const size_t ix, const size_t it, 
      const std::vector<double> &in,
      std::vector<double> &out) const 
{
   assert(in.size() ==_nx*_nphi*_nlat);
   assert(out.size()==_nphi);
   for (size_t ip=0; ip<_nphi; ip++) {
      out[ip] = in[indx_r_th_ph(ix,it,ip)];
   }
} 
void Grid::set_row_ph(const size_t ix, const size_t it, 
      const std::vector<double> &in,
      std::vector<double> &out) const 
{
   assert(in.size() ==_nphi);
   assert(out.size()==_nx*_nphi*_nlat);
   for (size_t ip=0; ip<_nlat; ip++) {
      out[indx_r_th_ph(ix,it,ip)] = in[ip];
   }
} 
/*=========================================================================*/
void Grid::get_row_n(const size_t il, 
      const std::vector<double> &in, 
      std::vector<double> &out) const
{
   assert(in.size() ==_nx*_nl);
   assert(out.size()==_nx);
   for (size_t ix=0; ix<_nx; ix++) {
      out[ix] = in[indx_r_l(ix,il)];
   }
} 
void Grid::set_row_n(const size_t il, 
      const std::vector<double> &in, 
      std::vector<double> &out) const
{
   assert(in.size() ==_nx);
   assert(out.size()==_nx*_nl);
   for (size_t ix=0; ix<_nx; ix++) {
      out[indx_r_l(ix,il)] = in[ix];
   }
} 
/*=========================================================================*/
void Grid::get_row_R_th(const size_t ip, 
      const std::vector<double> &in,
      std::vector<double> &out) const
{
   assert(in.size() ==_nx*_nphi*_nlat);
   assert(out.size()==_nx*_nlat);
   for (size_t ix=0; ix<_nx;   ix++) {
   for (size_t it=0; it<_nlat; it++) {
      out[indx_r_th(ix,it)] = in[indx_r_th_ph(ix,it,ip)];
   }
   }
}
void Grid::set_row_R_th(const size_t ip, 
      const std::vector<double> &in,
      std::vector<double> &out) const
{
   assert(in.size() ==_nx*_nlat);
   assert(out.size()==_nx*_nphi*_nlat);
   for (size_t ix=0; ix<_nx;   ix++) {
   for (size_t it=0; it<_nlat; it++) {
      out[indx_r_th_ph(ix,it,ip)] = in[indx_r_th(ix,it)];
   }
   }
}
/*=========================================================================*/
void Grid::get_row_R_ph(const size_t it, 
      const std::vector<double> &in,
      std::vector<double> &out) const
{
   assert(in.size() ==_nx*_nphi*_nlat);
   assert(out.size()==_nx*_nlat);
   for (size_t ix=0; ix<_nx;   ix++) {
   for (size_t ip=0; ip<_nphi; ip++) {
      out[indx_r_ph(ix,ip)] = in[indx_r_th_ph(ix,it,ip)];
   }
   }
}
void Grid::set_row_R_ph(const size_t it, 
      const std::vector<double> &in,
      std::vector<double> &out) const
{
   assert(in.size() ==_nx*_nphi);
   assert(out.size()==_nx*_nphi*_nlat);
   for (size_t ix=0; ix<_nx;   ix++) {
   for (size_t ip=0; ip<_nphi; ip++) {
      out[indx_r_th_ph(ix,it,ip)] = in[indx_r_ph(ix,ip)];
   }
   }
}
/*=========================================================================*/
void Grid::get_row_th_ph(const size_t ix, 
      const std::vector<double> &in,
      std::vector<double> &out) const
{
   assert(in.size() ==_nx*_nphi*_nlat);
   assert(out.size()==_nphi*_nlat);
   for (size_t it=0; it<_nlat; it++) {
   for (size_t ip=0; ip<_nphi; ip++) {
      out[indx_th_ph(it,ip)] = in[indx_r_th_ph(ix,it,ip)];
   }
   }
}
void Grid::set_row_th_ph(const size_t ix, 
      const std::vector<double> &in,
      std::vector<double> &out) const
{
   assert(in.size() ==_nphi*_nlat);
   assert(out.size()==_nx*_nphi*_nlat);

   for (size_t it=0; it<_nlat; it++) {
   for (size_t ip=0; ip<_nphi; ip++) {
      out[indx_r_th_ph(ix,it,ip)] = in[indx_th_ph(it,ip)];
   }
   }
}
/*=========================================================================*/
void Grid::set_row_l(const size_t ix, 
      const std::vector<double> &lvals, 
      std::vector<double> &Rlvals) const
{
   assert(lvals.size() ==_nl);
   assert(Rlvals.size()==_nx*_nl);
   for (size_t il=0; il<_nl; il++) {
      Rlvals[indx_r_l(ix,il)] = lvals[il];
   }
} 
/*==========================================================================*/
void Grid::set_partial_phi(const std::vector<double> &v, std::vector<double> &dv) const
{
   assert(v.size() ==_nx*_nlat*_nphi);
   assert(dv.size()==_nx*_nlat*_nphi);

#pragma omp parallel for
   for (size_t ix=0; ix<_nx; ix++) {
      std::vector<double> inter(   _nlat*_nphi);
      std::vector<double> inter_dv(_nlat*_nphi);

      get_row_th_ph(ix, v, inter); 

      _sphere.partial_phi(inter, inter_dv);

      set_row_th_ph(ix, inter_dv, dv); 
   }
}
/*==========================================================================*/
void Grid::set_spherical_lap(const std::vector<double> &v, std::vector<double> &ddv) const
{
   assert(v.size()  ==_nx*_nlat*_nphi);
   assert(ddv.size()==_nx*_nlat*_nphi);

#pragma omp parallel for
   for (size_t ix=0; ix<_nx; ix++) {
      std::vector<double> inter(    _nlat*_nphi);
      std::vector<double> inter_ddv(_nlat*_nphi);

      get_row_th_ph(ix, v, inter); 

      _sphere.laplace_beltrami(inter, inter_ddv);

      set_row_th_ph(ix, inter_ddv, ddv); 
   }
}
/*==========================================================================*/
void Grid::set_sphereX(const std::vector<double> &v, std::vector<double> &vX) const
{
   assert(v.size() ==_nx*_nlat*_nphi);
   assert(vX.size()==_nx*_nlat*_nphi);

#pragma omp parallel for
   for (size_t ix=0; ix<_nx; ix++) {
      std::vector<double> inter(   _nlat*_nphi);
      std::vector<double> inter_vX(_nlat*_nphi);

      get_row_th_ph(ix, v, inter); 

      _sphere.sphereX(inter, inter_vX);

      set_row_th_ph(ix, inter_vX, vX); 
   }
}
/*==========================================================================*/
void Grid::set_partial_r(const std::vector<double> &v, std::vector<double> &dv) const
{
   assert(v.size() ==_nx*_nlat*_nphi);
   assert(dv.size()==_nx*_nlat*_nphi);

#pragma omp parallel for
   for (size_t it=0; it<_nlat; it++) {
   for (size_t ip=0; ip<_nphi; ip++) {
      std::vector<double> inter(   _nx);
      std::vector<double> inter_dv(_nx);

      get_row_R(it, ip, v, inter); 

      _cheb.der(inter, inter_dv);

      for (size_t ix=0; ix<_nx; ix++) {
         inter_dv[ix] *= _partial_R_to_partial_r[ix];
      }

      set_row_R(it, ip, inter_dv, dv); 
   }
   }
}
/*==========================================================================*/
void Grid::set_angular_power_spectrum(const std::vector<double> &v, std::vector<double> &p) const
{
   assert(v.size()==_nx*_nlat*_nphi);
   assert(p.size()==_nx*_nl);

#pragma omp parallel for
   for (size_t ix=0; ix<_nx; ix++) {
      std::vector<double> inter__sphere(_nlat*_nphi, 0);
      std::vector<double> inter__sphere_p(_nl, 0);

      get_row_th_ph(ix, v, inter__sphere); 

      _sphere.power_spectrum(inter__sphere, inter__sphere_p);

      set_row_l(ix, inter__sphere_p, p); 
   }
}
/*==========================================================================*/
void Grid::set_n_l_coef(const std::vector<double> &v, std::vector<double> &p) const
{
   assert(v.size()==_nx*_nlat*_nphi);
   assert(p.size()==_nx*_nl);

   set_angular_power_spectrum(v, p);

#pragma omp parallel for
   for (size_t il=0; il<_nl; il++) {
      std::vector<double> inter_v(_nx, 0);
      std::vector<double> inter_p(_nx, 0);

      get_row_n(il, p, inter_v); 

      /* undo square root of Psl power spectrum */
      for (size_t ix=0; ix<_nx; ix++) {
         inter_v[ix] = pow(inter_v[ix],0.5);
      } 
      /* compute power spectrum */
      _cheb.to_ch(inter_v, inter_p);

      set_row_n(il, inter_p, p); 
   }
}
/*==========================================================================*/
/* Low pass filter in spectral space, but keep the boundaries
 * of the radial points fixed. */
/*==========================================================================*/
void Grid::filter(std::vector<double> &v) const
{
   assert(v.size()==_nx*_nlat*_nphi);

#pragma omp parallel for
   for (size_t ix=0; ix<_nx; ix++) {
      std::vector<double> inter__sphere(_nlat*_nphi);
      get_row_th_ph(ix, v, inter__sphere); 
      _sphere.filter(inter__sphere);
      set_row_th_ph(ix, inter__sphere, v); 
   }
#pragma omp parallel for
   for (size_t it=0; it<_nlat; it++) {
   for (size_t ip=0; ip<_nphi; ip++) {
      std::vector<double> inter_radial(_nx);
      get_row_R(it, ip, v, inter_radial); 
      const double lwr = inter_radial[0    ];
      const double upr = inter_radial[_nx-1];
      _cheb.filter(inter_radial);
      inter_radial[0    ] = lwr;
      inter_radial[_nx-1] = upr;
      set_row_R(it, ip, inter_radial, v); 
   }
   }
}
/*==========================================================================*/
double Grid::norm_indep_res(
      const double dt, 
      const std::vector<double> &f_n, 
      const std::vector<double> &f_np1, 
      const std::vector<double> &p) const 
{
   assert(f_n.size()  ==_nx*_nlat*_nphi);
   assert(f_np1.size()==_nx*_nlat*_nphi);
   assert(p.size()    ==_nx*_nlat*_nphi);

   double res = 0;

   for (size_t i=0; i<_nx*_nlat*_nphi; i++) {
      res += pow(p[i] - ((f_np1[i]-f_n[i])/dt),2);
   }
   res = pow(res/(_nx*_nlat*_nphi),0.5);

   return res;
}
/*==========================================================================*/
double Grid::total_variation(const std::vector<double> &v) const 
{
   assert(v.size()==_nx*_nlat*_nphi);

   double tv = 0;
   for (size_t ix=0; ix<_nx-1;   ix++) {
   for (size_t it=0; it<_nlat-1; it++) {
   for (size_t ip=0; ip<_nphi-1; ip++) {
      tv += pow(
               pow(
                  v[indx_r_th_ph(ix+1, it, ip)] 
               -  v[indx_r_th_ph(ix,   it, ip)]
               ,2)
            +  pow(
                  v[indx_r_th_ph(ix, it+1, ip)] 
               -  v[indx_r_th_ph(ix, it,   ip)]
               ,2)
            +  pow(
                  v[indx_r_th_ph(ix, it, ip+1)] 
               -  v[indx_r_th_ph(ix, it, ip)]
               ,2)
            ,
            0.5); 
   }
   }
   }
   return tv/(_nx*_nlat*_nphi);
}
/*===========================================================================*/
