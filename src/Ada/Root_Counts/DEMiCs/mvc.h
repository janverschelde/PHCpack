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

#ifndef __MVC_H
#define __MVC_H

#include "global.h"
#include "inputData.h"
#include "iTest.h"
#include "fTest.h"
#include "simplex.h"
#include "reltab.h"

class mvc
{
   int Dim;
   int supN;

   int row;
   int col;

   int termSumNum;
   int termMax;
   int maxLength;

   double total_iter;
   double total_feasLP;
   double total_LPs;
   double total_1PT;
   double total_2PT;
   double total_triLPs_mLP;

   double total_unbLP_tab;

   double* lvl_1PT;
   double* lvl_2PT;
   double* actNode;

   int* mfNum;

   int* firIdx;
   int* repN;

   int** negIdx;

   int* termSet; 
   int* termStart;
   int* re_termStart;

   int* type; 

   int* mRepN;
   int** mFeaIdx;
   int* mFea;

   int** trNeg;
   int* sp;

   int* candIdx;
   double* trMat;

   int* table;
 
   reltab Reltab;
   simplex Simplex;
   lvData* lv;
   iLvData* iLv;

   void getMemory ( int depth, int lvl, int length );
  
   void initMemoryCheck ( ftData& Data, int depth );
   void memoryCheck ( ftData* Data, int depth );

   void get_candIdx ( inifData& curInif );

   int chooseSup ( int depth, theData* curNode, 
                   inifData* curInif, inifData* nextInif );

   void fUpdateDirRed ( inifData* curInif, inifData* nextInif, 
                        theData* curNode, int* curRsp, int depth );

   void updateDirRed ( inifData* curInif, inifData* nextInif, 
                       theData* curNode, int* curRsp, int depth );

  // flag -- CONTINUE (Go to feasTest or findMixedCell)
  //      -- STOP (Break loop)
   int findUnbDir ( inifData* nextInif, theData* curNode, 
                    int* nextRsp, int* curRsp, int depth );

   int findUnbDir_art ( inifData* nextInif, theData* curNode, 
                        int* nextRsp, int*curRsp, int depth );

   int checkDir ( uData** corPtr, uData* tarPtr,
                  double* tar_dir, double tar_red,
                  int* nf_pos, int* basisIdx, int nfN );

   int checkDir_art ( uData** corPtr, uData* tarPtr,
                      double* tar_dir, double tar_red,
                      int* nf_pos, int* basisIdx, int nfN );

   void skipPtr ( uData** curr, uData** fHead );
  
   void get_tuple_index ( ftData* Node, ftData* Data, int length );

   void dbg_init_transMat ( theData* curNode );

   void dbg_transMat ( theData* preNode, theData* curNode );
   void check_transMat ( theData* preNode, theData* curNode );

   void check_init_transRed ( theData* curNode );

   int checkSign_red ( double curRed, double tarRed )
   {
      int flag;

      if(curRed > tarRed + PLUSZERO) flag = POSITIVE;
      else if(curRed < tarRed + MINUSZERO) flag = NEGATIVE;

      return (flag);
   };

   int checkNonNeg_dir ( double curDirElem, double tarDirElem )
   {
      if(curDirElem < tarDirElem + PLUSZERO) return(TRUE);
      else return(FALSE);
   };

   int checkNonPos_dir ( double curDirElem, double tarDirElem )
   {
      if(curDirElem > tarDirElem + MINUSZERO) return(TRUE);
      else return(FALSE);
   };

   int checkZero_dir ( double curDirElem, double tarDirElem )
   {
      double val;

      val = curDirElem - tarDirElem;

      if(MINUSZERO < val && val < PLUSZERO)
         return(TRUE);
      else
         return(FALSE);
   };

   int table_out ( int row, int col )
   {
      return (table[row + col * termSumNum]);
   };

   void info_neg ( int termSet, int** negIdx );
  
   void info_sp ( int depth );

   void info_parent_node ( int depth );
   void info_tuple ( int lvl, int depth );

   void info_all_dirRed ( int depth, ftData* Node, inifData* nextInif );

   void info_mFea ( int length );
   void info_firIdx ( int length );
   void info_fIdx ( ftData* Data );
   void info_candIdx();

   void info_elemNum ( int length, ftData* Data, ftData Node );
   void info_prop_elemNum ( int length, ftData* Data, ftData Node );
   void info_table();

   public:
  
      mvc();
      ~mvc();

  //// for mixed volume computation /////////////////////

#ifdef compile4phc
      void initialize_with_lifting
             ( dataSet& Data, double* lifvals, int seedNum, int output );
      /*
       * Extends the allocateAndIni with an array of lifting values
       * for all points in the supports. */
#endif
      void allocateAndIni ( dataSet& Data, int seedNum, int output );
      void initFeasTest ( int depth );
      void initCheck ( int depth, ftData& Data );
      void initLP ( ftData& Data, int** negIdx, int depth, 
                    int idx, int& feaNum );

  // flag  -- STOP (No node)
  //       -- Otherwise (Find a node)
      int feasTest(int depth, theData* parent);
      int upFeasTest(int& depth);

      void findMixedCell(int depth, theData* parent);
      void findAllMixedCells(int depth);

      int iCheck
            ( int depth, theData* parent, ftData& Data, inifData& inifData );

      void iLP ( theData* parent, ftData& Data, int depth, int idx_one, 
                 int fst_pivInIdx, int sub_fst_pivInIdx, 
                 int preNbN, int& feaNum );

      void iLP_Art ( theData* parent, ftData& Data, int depth, int idx_one, 
                     int fst_pivInIdx, int sub_fst_pivInIdx, 
                     int preNbN, int& feaNum );

  // flag -- FNN (Output a node and go to findeNextNode)
  //      -- STOP (Break loop)
      int findNode ( int depth, int& lvl, int& feaNum, ftData* Data );


  // flag -- CONTINUE (Go to "SLIDE" and findNode)
  //      -- FNN (Again go to findNextNode)
  //      -- STOP(Terminate)
      int findNextNode ( int depth, int& lvl, int& feaNum, ftData* Data );

      void findUpNode ( ftData* Data, ftData** pre, ftData** cur, int& lvl, 
                        int polyDim, int depth );

      int mLP ( ftData* Pre, ftData* Cur, ftData* Data,
                int* repIdx, int* feaIdx, int tarIdx, int* mRepN,
                int totalN, int depth, int& feaNum, int lvl, int length );

      int checkBasis ( theData* target, int sub_sIdx );
      int checkAnotherBasis ( int repIdx, int dist, theData** target );

      void get_firIdx ( ftData data_a, ftData data_b, int sn, int lvl );

      void info_cpuTime ( double cpuTime_start, double cpuTime_end );
      void info_final();

      void Enum();
};

#endif
