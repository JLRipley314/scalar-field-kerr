#include <cassert>
#include <fstream>
#include <iostream>
#include <string>
#include <iomanip>
#include <fstream>

#include "params.hpp"

namespace Params
{
   /*========================================================================*/
   namespace 
   {
      size_t _nt;
      size_t _nx;
      size_t _nl;
      size_t _nlat;
      size_t _nphi;
      size_t _t_step_save;

      double _dt;

      double _cl;

      double _Rmin;
      double _Rmax;

      double _bh_mass;
      double _bh_spin;

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

      std::string read(
            const std::string param_file, 
            const std::string find_this_var)
      {
         std::ifstream infile(param_file);
         assert(infile.good());

         std::string name, val;

         while (infile >> name >> val) {
            if (name==find_this_var) {
               std::cout<<name<<" "<<val<<std::endl;
               return val;
            }
         }
         std::cout<<"ERROR: did not find "<<find_this_var<<std::endl;
         std::quick_exit(0);
      }
   }
   /*========================================================================*/
   size_t nt() {return _nt;};
   size_t nx() {return _nx;};
   size_t nl() {return _nl;};

   size_t nx_nlat_nphi() { return _nx*_nlat*_nphi; }
   size_t nxn() {return _nl;};

   size_t nlat() {return _nlat;};
   size_t nphi() {return _nphi;};
   size_t t_step_save() {return _t_step_save;};

   double dt() {return _dt;};

   double cl() {return _cl;};

   double Rmax() {return _Rmax;};
   double Rmin() {return _Rmin;};

   double bh_mass() {return _bh_mass;};
   double bh_spin() {return _bh_spin;};

   double V2() {return _V2;};
   double V3() {return _V3;};
   double V4() {return _V4;};

   std::string id() {return _id;};

   std::string initial_data_direction() {return _initial_data_direction;};

   double amp() {return _amp;};
   double rl() {return _rl;};
   double ru() {return _ru;};

   int l_ang() {return _l_ang;};
   int m_ang() {return _m_ang;};
   /*========================================================================*/
   void init(const std::string param_file)
   {
      _nt   = std::stoi(read(param_file,"nt"));
      _nx   = std::stoi(read(param_file,"nx"));
      _nl   = std::stoi(read(param_file,"nl"));
      _nlat = std::stoi(read(param_file,"nlat"));
      _nphi = std::stoi(read(param_file,"nphi"));

      _t_step_save = std::stoi(read(param_file,"t_step_save"));

      _dt = std::stod(read(param_file,"dt"));

      _cl =   std::stod(read(param_file,"compactification_length"));

      _Rmax =   std::stod(read(param_file,"Rmax"));
      _Rmin =   std::stod(read(param_file,"Rmin"));

      _bh_mass = std::stod(read(param_file,"black_hole_mass"));
      _bh_spin = std::stod(read(param_file,"black_hole_spin"));
      /*---------------------------------------------------------------------*/
      /* for the potentials */
      _V2 = std::stod(read(param_file,"V2"));
      _V3 = std::stod(read(param_file,"V3"));
      _V4 = std::stod(read(param_file,"V4"));
      /*---------------------------------------------------------------------*/
      /* for the initial data */
      _id = read(param_file,"initial_data_type");

      _initial_data_direction = read(param_file,"initial_data_direction");

      _amp = std::stod(read(param_file,"amp"));
      _rl = std::stod(read(param_file,"rl"));
      _ru = std::stod(read(param_file,"ru"));

      _l_ang = std::stoi(read(param_file,"l_ang"));
      _m_ang = std::stoi(read(param_file,"m_ang"));
   }
   /*========================================================================*/
}
