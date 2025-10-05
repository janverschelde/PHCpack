/*
    This file is a component of DEMiCs
    Copyright (C) 2007 Tomohiko Mizutani, Masakazu Kojima and Akiko Takeda

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#ifndef __INPUTDATA_H
#define __INPUTDATA_H

#include "global.h"

class dataSet
{
   public:

      dataSet();       // constructor sets everythin to 0 or NULL
      ~dataSet();      // destructor deallocates memory

      int Dim;         // dimension of the points
      int supN;        // number of distinct supports
      int termSumNum;  // total number of points in all supports
      int termMax;     // largest number in the array termSet
      int typeMax;     // largest number in the array type

      int* termSet;    // array of supN integers with number of points
                       // in each support set
      int* termStart;  // array of supN+1 integers with the index to the 
                       // first point in each support set, the last element 
                       // in termStart is the total number of points
      int* type;       // array of supN integers with the number
                       // of occurrences of each distinct support

      double* support; // coordinates of the points in each support
      double* coef;

      char* outFile;

      void support_in ( int rowIdx, int colIdx, double elem )
      /*
       * DESCRIPTION :
       *   Assigns the elem to the coordinate of position colIdx
       *   of support with index rowIdx. */
      {
         support[colIdx + Dim * rowIdx] = elem; 
      };

      double support_out ( int rowIdx, int colIdx )
      /*
       * DESCRIPTION :
       *   Returns the coordinate of index colIdx
       *   of the point with index rowIdx. */
      {
         return (support[colIdx + Dim * rowIdx]);
      };
  
      void getInputFile ( char* inputFile );
      /*
       * DESCRIPTION :
       *   Opens the file with name defined by inputFile,
       *   parses the information on the file into the data
       *   on the support sets. */

      void info_preamble();
      /*
       * DESCRIPTION :
       *   Writes the dimension, number of distinct supports,
       *   the number of points in each support and
       *   the number of occurrences of each support to screen. */
  
      void info_supports();
      /*
       * DESCRIPTION :
       *    Writes each support set to screen. */
};

#endif
