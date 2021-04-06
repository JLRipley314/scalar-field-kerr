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
      size_t nt_;
      size_t nx_;
      size_t nl_;
      size_t nlat_;
      size_t nphi_;
      size_t t_step_save_;

      double dt_;

      double cl_;

      double Rmin_;
      double Rmax_;

      double bh_mass_;
      double bh_spin_;

      double V_1_;
      double V_2_;
      double V_3_;
      double V_4_;

      std::string id_;
      std::string initial_data_direction_;

      double amp_;
      double rl_;
      double ru_;

      int l_ang_;
      int m_ang_;

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
   size_t nt() {return nt_;};
   size_t nx() {return nx_;};
   size_t nl() {return nl_;};

   size_t nx_nphi_nlat() { return nx_*nphi_*nlat_; }
   size_t nxn() {return nl_;};

   size_t nlat() {return nlat_;};
   size_t nphi() {return nphi_;};
   size_t t_step_save() {return t_step_save_;};

   double dt() {return dt_;};

   double cl() {return cl_;};

   double Rmax() {return Rmax_;};
   double Rmin() {return Rmin_;};

   double bh_mass() {return bh_mass_;};
   double bh_spin() {return bh_spin_;};

   double V_1() {return V_1_;};
   double V_2() {return V_2_;};
   double V_3() {return V_3_;};
   double V_4() {return V_4_;};

   std::string id() {return id_;};

   std::string initial_data_direction() {return initial_data_direction_;};

   double amp() {return amp_;};
   double rl() {return rl_;};
   double ru() {return ru_;};

   int l_ang() {return l_ang_;};
   int m_ang() {return m_ang_;};
   /*========================================================================*/
   void init(const std::string param_file)
   {
      nt_   = std::stoi(read(param_file,"nt"));
      nx_   = std::stoi(read(param_file,"nx"));
      nl_   = std::stoi(read(param_file,"nl"));
      nlat_ = std::stoi(read(param_file,"nlat"));
      nphi_ = std::stoi(read(param_file,"nphi"));

      t_step_save_ = std::stoi(read(param_file,"t_step_save"));

      dt_ = std::stod(read(param_file,"dt"));

      cl_ =   std::stod(read(param_file,"compactification_length"));

      Rmax_ =   std::stod(read(param_file,"R_max"));
      Rmin_ =   std::stod(read(param_file,"R_min"));

      bh_mass_ = std::stod(read(param_file,"black_hole_mass"));
      bh_spin_ = std::stod(read(param_file,"black_hole_spin"));
      /*---------------------------------------------------------------------*/
      /* for the potentials */
      V_1_ = std::stod(read(param_file,"V_1"));
      V_2_ = std::stod(read(param_file,"V_2"));
      V_3_ = std::stod(read(param_file,"V_3"));
      V_4_ = std::stod(read(param_file,"V_4"));
      /*---------------------------------------------------------------------*/
      /* for the initial data */
      id_ = read(param_file,"initial_data_type");

      initial_data_direction_ = read(param_file,"initial_data_direction");

      amp_ = std::stod(read(param_file,"amp"));
      rl_ = std::stod(read(param_file,"rl"));
      ru_ = std::stod(read(param_file,"ru"));

      l_ang_ = std::stoi(read(param_file,"l_ang"));
      m_ang_ = std::stoi(read(param_file,"m_ang"));
   }
   /*========================================================================*/
}
