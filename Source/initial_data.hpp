/*
 * writes initial data for scalar field variables 
 */
#ifndef _ID_HPP_
#define _ID_HPP_

#include <vector>

/*==========================================================================*/
namespace ID 
{
/*==========================================================================*/
void ingoing_pulse(
      std::vector<double> &f,
      std::vector<double> &p,
      std::vector<double> &q
   );
/*==========================================================================*/
} /* ID */
/*==========================================================================*/
#endif /* _ID_HPP_ */