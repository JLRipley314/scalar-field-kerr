/*
 * writes initial data for scalar field variables 
 */
#ifndef _ID_HPP_
#define _ID_HPP_

#include <string>
#include <vector>

/*==========================================================================*/
namespace ID 
{
/*==========================================================================*/
void compact_pulse(
      std::string initial_data_direction,
      std::vector<double> &f,
      std::vector<double> &p
   );
/*==========================================================================*/
} /* ID */
/*==========================================================================*/
#endif /* _ID_HPP_ */
