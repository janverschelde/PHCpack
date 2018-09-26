// polysyshom.h contains prototypes of templated data types for homotopies
// of polynomial systems with complex coefficients of different precisions.

#ifndef __POLYSYSHOM_H__
#define __POLYSYSHOM_H__

#include <iostream>
#include <fstream>
#include "stdlib.h"
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>

template <class ComplexType, class RealType>
class PolySysHom
{
   public:

      PolySys<ComplexType,RealType>* start_sys;
      PolySys<ComplexType,RealType>* target_sys;
      int dim;

      PolySysHom ( PolySys<ComplexType,RealType>* start_sys,
                   PolySys<ComplexType,RealType>* target_sys )
      {
         if(start_sys->dim != target_sys->dim)
         {
            std::cout << "start system and end system " << std::endl;
         }
         else
         {
            this->start_sys = start_sys;
            this->target_sys = target_sys;
            dim = start_sys->dim;
         }
      }

      void print()
      {
         std::cout << "Start System : " << std::endl;
         start_sys->print();
         std::cout << "Target System : " << std::endl;
         target_sys->print();
      }
};

#endif
