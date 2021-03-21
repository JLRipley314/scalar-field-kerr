#include <vector>
#include <fftw3.h>

class Cheb {
   private:
      const int n;

      const double jacobian;

      std::vector<double> low_pass;
      /* for fftw Fourier transform */
      double *in;
      double *out;
      fftw_plan plan_dct;
   public:
      Cheb(const int n, const double jacobian);
      ~Cheb();

      void to_ch(const std::vector<double> &po, std::vector<double> &ch);
      void to_po(const std::vector<double> &ch, std::vector<double> &po);

      void der(const std::vector<double> &v, 
            std::vector<double> &ch, std::vector<double> &dv
            );

      void filter(std::vector<double> &ch, std::vector<double> &v);
};
