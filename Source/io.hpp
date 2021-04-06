#ifndef _IO_HPP_
#define _IO_HPP_
/*
 * Utility functions for writing data to file.
 */

#include <string>
#include <vector>
/*=========================================================================*/
namespace Csv 
/*=========================================================================*/
{
   void write(
         const std::string name, 
         const size_t itm, 
         const std::vector<double> &vals 
      );	
   void write_unstructured(
         const std::string name, 
         const int itm, 
         const std::vector<std::vector<double>> &grid,
         const std::vector<double> &vals
      );
/*=========================================================================*/
}; /* Csv */
/*=========================================================================*/
#endif /* _IO_HPP_ */
