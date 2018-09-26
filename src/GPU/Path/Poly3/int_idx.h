// int_idx.h contains the class definition for an index type

#ifndef __INT_IDX_H__
#define __INT_IDX_H__

#include <iostream>
#include <fstream>
#include "stdlib.h"
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>

class int_idx
{
   public:

      int eq_idx;
      int mon_idx;

      int_idx ( int i, int j )
      {
         eq_idx = i;
         mon_idx = j;
      }
};

#endif
