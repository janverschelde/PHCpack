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

      dataSet();  // constructor sets everythin to 0 or NULL
      ~dataSet(); // destructor deallocates memory

      int Dim, supN, termSumNum;
      int termMax, typeMax;

      int* termSet, * termStart, * type;

      double* support;
      double* coef;

      char* outFile;

      void support_in ( int rowIdx, int colIdx, double elem )
      {
         support[colIdx + Dim * rowIdx] = elem; 
      };

      double support_out ( int rowIdx, int colIdx )
      {
         return (support[colIdx + Dim * rowIdx]);
      };
  
  ///// functions  /////

      void getInputFile ( char* inputFile );

  ///// output infomation on display /////

      void info_pre();
      void info_support();
};

#endif
