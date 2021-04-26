#include <cassert>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <fstream>
#include <sstream>

#include "params.hpp"

/*========================================================================*/
Params::Params(const std::string param_file)
{
   std::cout<<"Initializing Params"<<std::endl;
   _nt   = std::stoi(_read(param_file,"nt"));
   _nl   = std::stoi(_read(param_file,"nl"));
   _nm   = std::stoi(_read(param_file,"nm"));
   _nlat = std::stoi(_read(param_file,"nlat"));
   _nphi = std::stoi(_read(param_file,"nphi"));

   _ngrids = std::stoi(_read(param_file,"ngrids"));
   _read_array(param_file,"nxs",     _nxs);
   _read_array(param_file,"Rvals", _Rvals);

   assert(_ngrids  ==  _nxs.size());
   assert(_ngrids+1==_Rvals.size());

   _t_step_save = std::stoi(_read(param_file,"t_step_save"));

   _dt = std::stod(_read(param_file,"dt"));

   _cl =   std::stod(_read(param_file,"compactification_length"));

   _Rmax =   std::stod(_read(param_file,"Rmax"));
   _Rmin =   std::stod(_read(param_file,"Rmin"));

   _bh_mass = std::stod(_read(param_file,"black_hole_mass"));
   _bh_spin = std::stod(_read(param_file,"black_hole_spin"));
   /*---------------------------------------------------------------------*/
   /* for the potentials */
   _km1 = std::stod(_read(param_file,"km1"));
   _k0  = std::stod(_read(param_file,"k0"));
   _k1  = std::stod(_read(param_file,"k1"));
   _k2  = std::stod(_read(param_file,"k2"));

   _V2 = std::stod(_read(param_file,"V2"));
   _V3 = std::stod(_read(param_file,"V3"));
   _V4 = std::stod(_read(param_file,"V4"));
   /*---------------------------------------------------------------------*/
   /* for the initial data */
   _id = _read(param_file,"initial_data_type");

   _initial_data_direction = _read(param_file,"initial_data_direction");

   _amp = std::stod(_read(param_file,"amp"));
   _rl = std::stod(_read(param_file,"rl"));
   _ru = std::stod(_read(param_file,"ru"));

   _l_ang = std::stoi(_read(param_file,"l_ang"));
   _m_ang = std::stoi(_read(param_file,"m_ang"));
   std::cout<<"Finished initializing Params"<<std::endl;
}
/*========================================================================*/
Params::~Params()
{
}
/*========================================================================*/
std::string Params::_read(
      const std::string param_file, 
      const std::string find_this_var)
{
   std::ifstream infile(param_file);
   assert(infile.good());

   std::string line;

   while (std::getline(infile,line)) {
      std::stringstream linestream(line);

      std::string name;
      std::getline(linestream, name, ' ');
      if (name==find_this_var) {

         std::string val;
         linestream >> val;
         
         std::cout<<name<<"\t"<<val<<std::endl;

         return val;
      }
   }
   std::cout<<"ERROR: did not find "<<find_this_var<<std::endl;
   std::quick_exit(0);
}
/*========================================================================*/
void Params::_read_array(
      const std::string param_file, 
      const std::string find_this_var,
      std::vector<size_t> &arr)
{
   std::ifstream infile(param_file);
   assert(infile.good());

   std::string line;

   while (std::getline(infile,line)) {
      std::stringstream linestream(line);

      std::string name;
      std::getline(linestream, name, ' ');

      if (name==find_this_var) {
         std::cout<<name;
         size_t val;
         while (linestream >> val) {
            arr.push_back(size_t(val));
            std::cout<<" "<<val;
         }
         std::cout<<std::endl;
         return;
      }
   }
   std::cout<<"ERROR: did not find "<<find_this_var<<std::endl;
   std::quick_exit(0);
}
/*========================================================================*/
void Params::_read_array(
      const std::string param_file, 
      const std::string find_this_var,
      std::vector<double> &arr)
{
   std::ifstream infile(param_file);
   assert(infile.good());

   std::string line;

   while (std::getline(infile,line)) {
      std::stringstream linestream(line);

      std::string name;
      std::getline(linestream, name, ' ');
      if (name==find_this_var) {
         std::cout<<name;
         double val;
         while (linestream >> val) {
            arr.push_back(val);
            std::cout<<" "<<val;
         }
         std::cout<<std::endl;
         return;
      }
   }
   std::cout<<"ERROR: did not find "<<find_this_var<<std::endl;
   std::quick_exit(0);
}
