/*
 * utility routines for reading from parameter file
 * and accessing parameters in code
 */
#ifndef _PARAMS_HPP_
#define _PARAMS_HPP_

#include <string>
#include <vector>

/*===========================================================================*/
class Params
{
public:
   Params() {}; /* default: does nothing */
   Params(const std::string output_dir);
   ~Params();

private: 
   /* access the parameter variables with getter functions; 
    * see below under `public'
    */
   size_t _nt;
   size_t _nl;
   size_t _nm;
   size_t _nlat;
   size_t _nphi;
   size_t _t_step_save;

   size_t _ngrids;
   std::vector<double> _Rvals;
   std::vector<size_t> _nxs;

   double _dt;

   double _cl;

   double _Rmin;
   double _Rmax;

   double _bh_mass;
   double _bh_spin;

   double _km1;
   double _k0;
   double _k1;
   double _k2;

   double _V2;
   double _V3;
   double _V4;

   std::string _id;
   std::string _initial_data_direction;

   double _amp;
   double _rl;
   double _ru;

   int _l_ang;
   int _m_ang;

   /* reads line from parameter file */
   std::string _read(
         const std::string param_file, 
         const std::string find_this_var);

   void _read_array(
         const std::string param_file, 
         const std::string find_this_var,
         std::vector<size_t> &arr);

   void _read_array(
         const std::string param_file, 
         const std::string find_this_var,
         std::vector<double> &arr);

public:

   inline size_t nt() const {return _nt;}
   inline size_t nl() const {return _nl;}
   inline size_t nm() const {return _nm;}

   inline size_t ngrids() const {return _ngrids;}
   inline size_t nxs(const size_t i) const {return _nxs[i];}

   inline size_t nxn() const {return _nl;}

   inline size_t nlat() const {return _nlat;}
   inline size_t nphi() const {return _nphi;}
   inline size_t t_step_save() const {return _t_step_save;}

   inline double dt() const {return _dt;}

   inline double cl() const {return _cl;}

   inline double Rmax() const {return _Rmax;}
   inline double Rmin() const {return _Rmin;}
   inline double Rvals(const size_t i) const {return _Rvals[i];}

   inline double bh_mass() const {return _bh_mass;}
   inline double bh_spin() const {return _bh_spin;}

   inline double km1() const {return _km1;}
   inline double k0()  const {return _k0;}
   inline double k1()  const {return _k1;}
   inline double k2()  const {return _k2;}

   inline double V2() const {return _V2;}
   inline double V3() const {return _V3;}
   inline double V4() const {return _V4;}

   inline std::string id() const {return _id;}

   inline std::string initial_data_direction() const {return _initial_data_direction;}

   inline double amp() const {return _amp;}
   inline double rl() const {return _rl;}
   inline double ru() const {return _ru;}

   inline int l_ang() const {return _l_ang;}
   inline int m_ang() const {return _m_ang;}

};
#endif // _PARAMS_HPP_
