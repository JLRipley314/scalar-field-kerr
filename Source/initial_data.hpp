/*
 * writes initial data for scalar field variables 
 */
#ifndef ID_HPP__
#define ID_HPP__

#include <string>
#include <vector>

#include "grid.hpp"
#include "params.hpp"

/*==========================================================================*/
namespace ID 
{
/*==========================================================================*/
void compact_pulse(
      const Params &params,
      const Grid &grid,
      std::vector<double> &f,
      std::vector<double> &p
   );
/*==========================================================================*/
} /* ID */
/*==========================================================================*/
#endif /* _ID_HPP_ */
