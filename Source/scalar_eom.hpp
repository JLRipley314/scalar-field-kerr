/*
 * Field equations for scalar field about a Kerr black hole.
 */
#ifndef _SCALAR_EOM_HPP_ 
#define _SCALAR_EOM_HPP_

#include "field.hpp"

/*==========================================================================*/
namespace Eom {
/*==========================================================================*/
void init();
void cleanup();

void time_step(Field &f, Field &p, Field &q);

/*==========================================================================*/
} /* Eom */
#endif /* _SCALAR_EOM_HPP_ */
