/*
 * Field equations for scalar field about a Kerr black hole.
 */
#ifndef _SCALAR_EOM_HPP_ 
#define _SCALAR_EOM_HPP_

#include <vector>
#include "field.hpp"

/*==========================================================================*/
namespace Eom {
/*==========================================================================*/
void init();

void time_step(Field &f, Field &p);

/* local energy density */
void set_rho(
      const std::vector<double> &f,
      const std::vector<double> &p,
      std::vector<double> &rho
      );

/*==========================================================================*/
} /* Eom */
#endif /* _SCALAR_EOM_HPP_ */
