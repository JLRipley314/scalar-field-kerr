#ifndef _ID_HPP_
#define _ID_HPP_

#include <vector>

/*==========================================================================*/
namespace ID 
{
/*==========================================================================*/
void ingoing_pulse(
      const double amp,
      const double width,
      const double rl,
      const double ru,
      const std::vector<double> r,
      Field3d &f,
      Field3d &p,
      Field3d &q
   );
/*==========================================================================*/
} /* ID */
/*==========================================================================*/
#endif /* _ID_HPP_ */
