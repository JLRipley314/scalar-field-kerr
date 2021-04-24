/*
 * writes initial data for scalar field variables 
 */
#ifndef _ID_HPP_
#define _ID_HPP_

#include <string>
#include <vector>

#include "grid.hpp"
#include "params.hpp"

/*==========================================================================*/
namespace ID 
{
/*==========================================================================*/
void compact_pulse(
      const Grid &grid,
      const Params &params,
      std::vector<double> &f,
      std::vector<double> &p
   );
/*==========================================================================*/
} /* ID */
/*==========================================================================*/
#endif /* _ID_HPP_ */
