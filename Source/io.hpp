/*
 * Utility functions for writing data to file.
 */
#ifndef _IO_HPP_
#define _IO_HPP_

#include <string>

#include "field.hpp"
/*=========================================================================*/
namespace Csv 
/*=========================================================================*/
{
   void write(const std::string outdir, const double time, 
         const class Field &field
      );	
/*=========================================================================*/
}; /* Csv */
/*=========================================================================*/
#endif /* _IO_HPP_ */
