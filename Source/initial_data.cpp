#include "initial_data.hpp"
#include "sphere.hpp"

/*==========================================================================*/
namespace ID 
{
/*==========================================================================*/
/* time symmetric compact bump for phi of angular dependence Y_{lm} */
/*==========================================================================*/
void ingoing_pulse(
      const double amp,
      const double width,
      const double rl,
      const double ru,
      const int l_ang,
      const in m_ang,
      const std::vector<double> rv,
      Field3d &f,
      Field3d &q,
      Field3d &f
   )
{
   assert(r.size()==f.nx());

   double max_val = 0.0;

   double width = abs(rl-ru);
   if (abs(width)<1e-14) width = 1e-14; /* don't want to divide by zero */

   std::vector<double> ylm = Sphere::compute_ylm(l_ang, m_ang);

   for (size_t i=0; i<f.nx; i++) {
   for (size_t j=0; j<f.ny; j++) {
   for (size_t k=0; k<f.nz; k++) {
      double r    = rv[i];
      double bump = 0.0;

      if ((r<ru) && (r>rl)) {
         bump = exp(-1.0*width/(r-rl))*exp(-2.0*width/(ru-r));
      }
      f.vals[f.indx(i,j,k)] = (((r-rl)/width)**2) * (((ru-r)/width)**2) * bump; 

      q.vals[q.indx(i,j,j)] = (
         (2.0*(((r-rl)/width)   )*pow(((ru-r)/width),2))
      -  (2.0*pow((r-rl)/width,2)*((ru-r)/width       ))
      +  (1.0*(1.0_             )*pow(((ru-r)/width),2))
      -  (2.0*pow((r-rl)/width,2)*(1.0                ))
      )*bump/width;

      q.vals[q.indx(i,j,j)] *= -pow(r/cl,2);
      /*
       * time symmetric for now
       */ 
      p.vals[p.indx(i,j,j)] = 0.0;
      /*
       * give angular structure Y_{lm} 
       */
      f.vals[f.indx(i,j,k)] *= ylm[Sphere::indx_Sph(j,k)];
      q.vals[q.indx(i,j,k)] *= ylm[Sphere::indx_Sph(j,k)];
      p.vals[p.indx(i,j,k)] *= ylm[Sphere::indx_Sph(j,k)];

      if (abs(f.vals[f.indx(i,j,k)]) > max_val) {
         max_val = abs(f.vals[f.indx(i,j,k)])
      }
   }
   }
   }
   /* 
    * rescale so initial amplitude is amp
    */
   for (size_t i=0; i<f.nx; i++) {
   for (size_t j=0; j<f.ny; j++) {
   for (size_t k=0; k<f.nz; k++) {
      p.vals[f.indx(i,j,k)] *= amp / max_val;
      p.vals[q.indx(i,j,k)] *= amp / max_val;
      p.vals[q.indx(i,j,k)] *= amp / max_val; 
   }
   }
   }
}
/*==========================================================================*/
} /* ID */
/*==========================================================================*/
#endif /* _ID_HPP_ */

