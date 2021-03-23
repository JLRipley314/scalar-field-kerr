#include <cassert>
#include <fstream>
using std::ifstream;
#include <iostream>
using std::cout; 
using std::endl;
#include <string>
using std::string;
using std::stoi;
using std::stod;
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
      size_t nlat_;
      size_t nphi_;
      int t_step_save_;
      int direction_;

      double dt_;

      double curv_;
      double cl_;

      int initial_exc_i_;

      double bh_mass_;

      double V_0_;
      double V_1_;
      double V_2_;
      double V_3_;
      double V_4_;

      std::string id_;

      double amp_;
      double r_l_;
      double r_u_;

      int l_ang_;
      int m_ang_;

      string read(const string output_dir, const string find_this_var)
      {
         ifstream infile(output_dir+"/params.txt");
         assert(infile.good());

         string name, val;

         while (infile >> name >> val) {
            if (name==find_this_var) {
               return val;
            }
         }
         cout << "ERROR: did not find " << find_this_var << endl;
         std::quick_exit(0);
      }
   }
   /*========================================================================*/
    size_t nt() {return nt_;};
    size_t nx() {return nx_;};
    size_t nlat() {return nlat_;};
    size_t nphi() {return nphi_;};
    int t_step_save() {return t_step_save_;};
    int direction() {return direction_;};

    double dt() {return dt_;};

    double curv() {return curv_;};
    double cl()   {return cl_;};

    int initial_exc_i() {return initial_exc_i_;};

    double bh_mass() {return bh_mass_;};

    double V_0() {return V_0_;};
    double V_1() {return V_1_;};
    double V_2() {return V_2_;};
    double V_3() {return V_3_;};
    double V_4() {return V_4_;};

    std::string id() {return id_;};

    double amp() {return amp_;};
    double r_l() {return r_l_;};
    double r_u() {return r_u_;};

    int l_ang() {return l_ang_;};
    int m_ang() {return m_ang_;};
   /*========================================================================*/
   void init(const string output_dir)
   {
      nt_   = stoi(read(output_dir,"nt"));
      nx_   = stoi(read(output_dir,"nx"));
      nlat_ = stoi(read(output_dir,"nlat"));
      nphi_ = stoi(read(output_dir,"nphi"));

      t_step_save_ = stoi(read(output_dir,"t_step_save"));

      direction_ = stoi(read(output_dir,"direction"));

      dt_ = stod(read(output_dir,"dt"));

      curv_ = stod(read(output_dir,"curvature"));
      cl_ =   stod(read(output_dir,"compactification_scale"));

      initial_exc_i_ = stoi(read(output_dir,"initial_exc_i"));

      bh_mass_ = stod(read(output_dir,"bh_mass"));
      /*---------------------------------------------------------------------*/
      /* for the potentials */
      V_0_ = stod(read(output_dir,"V_0"));
      V_1_ = stod(read(output_dir,"V_1"));
      V_2_ = stod(read(output_dir,"V_2"));
      V_3_ = stod(read(output_dir,"V_3"));
      V_4_ = stod(read(output_dir,"V_4"));
      /*---------------------------------------------------------------------*/
      /* for the initial data */
      id_ = read(output_dir,"initial_data_type");

      amp_ = stod(read(output_dir,"amp"));
      r_l_ = stod(read(output_dir,"r_l"));
      r_u_ = stod(read(output_dir,"r_u"));

      l_ang_ = stoi(read(output_dir,"l_ang"));
      m_ang_ = stoi(read(output_dir,"m_ang"));
   }
   /*========================================================================*/
}
