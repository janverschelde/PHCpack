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

#ifndef __RELTAB_H
#define __RELTAB_H

#include "global.h"
#include "inputData.h"
#include "simplex.h"

class reltab
{
   private:

      int Dim;
      int supN;
      int maxConst;
      int termSumNum;
  
      int* termSet;
      int* termStart;
      int* re_termStart;
      int* firIdx;

      double unbLP;
      double totalLP;

      int row;
      int col;

      int nbN;
      int nfN;

      double *invB;

      double *p_sol;
      double *d_sol;

      int *basisIdx;
      int *nbIdx; 

      int *nf_pos;

      int *negIdx;
      double *val;

      int *feasIdx_a;
      int *feasIdx_b;

      simplex *Simplex;

      void get_init_triData ( int lab, int idx );
      void get_init_squData ( int lab_a, int lab_b, int idx_a, int idx_b, 
                              int colPos, int rowPos );

      void init_data();
      void init_tri ( int lab, int idx );
      void init_squ ( int lab_a, int lab_b, int idx_a, int idx_b );

      void put_data();
      void put_frIdx ( int frIdx );

      void makeTri();
      void makeSqu();

      void findAllFeasLPs_tri ( int lab, int idx, int frIdx );
      void findAllFeasLPs_squ ( int lab_a, int lab_b, int idx_a, int idx_b, 
                                int colPos, int rowPos );

      void table_in ( int row, int col, int elem )
      {
         table[row + col * termSumNum] = elem;
      };

      int table_out ( int row, int col )
      {
         return (table[row + col * termSumNum]);
      };

      void info_invB();

      void info_p_sol();
      void info_d_sol();

      void info_basisIdx();
      void info_nbIdx(); 
  
      void info_nf_pos();
      void info_feasIdx_tri(int num);
      void info_feasIdx_squ(int num_a, int num_b);

      void info_allTable();
      void info_table();

      double invB_out ( int rowIdx, int colIdx )
      {
         return (invB[colIdx + Dim * rowIdx]);
      };

   public:

      reltab();
      ~reltab();

      int* table;

      void allocateAndIni
       ( simplex& ori_Simplex, int** ori_firIdx, int ori_Dim, int ori_supN,
         int ori_termSumNum, int* ori_termSet, int* ori_termStart,
         int* ori_re_termStart );

      void makeTable ( double& total_unbLP_tab );
};

#endif
