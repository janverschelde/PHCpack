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

#ifndef __FTEST_H
#define __FTEST_H

#include "global.h"
#include "inputData.h"
#include "iTest.h"

class theData
{
   private:

      int row, col;
      int termS;

   public:

      theData();
      ~theData();

      theData *next;

      int flag;
      int polyDim;
      int nbN;
      int nfN;

      int artV;
      int pivOutNum;

      int fIdx;

      int sw;

      double *invB;  
      double *transMat;
      double *transRed;

      double *p_sol; 
      double *d_sol; 

      double *redVec;

      int *basisIdx; 
      int *nbIdx; 

      int *nf_pos;
      int *rIdx; 

      int *pivOutList;
      int *pivOutCheck;

      double *invB_ptr;  
      double *transMat_ptr;
      double *transRed_ptr;

      double *p_sol_ptr; 
      double *d_sol_ptr; 

      double *redVec_ptr;

      int *basisIdx_ptr; 
      int *nbIdx_ptr; 

      int *nf_pos_ptr;

      int* nodeLabel;

      void create ( int ori_row, int ori_col, int ori_termS, int ori_polyDim );

      void joint();
      void iJoint();
      void mJoint();

      void clear();
      void clear_transMat();

      void put_info(int repIdx, int& idx2, int& lNbN, int& lNfN);

      double invB_out ( int rowIdx, int colIdx )
      {
         return (invB[colIdx + row * rowIdx]);
      };

      double transMat_out ( int rowIdx, int colIdx )
      {
         return (transMat[colIdx + row * rowIdx]);
      };

      double invB_ptr_out ( int rowIdx, int colIdx )
      {
         return (invB_ptr[colIdx + row * rowIdx]);
      };

      double transMat_ptr_out ( int rowIdx, int colIdx )
      {
         return (transMat_ptr[colIdx + row * rowIdx]);
      };

      void info_p_sol();
      void info_d_sol();
      void info_invB();
      void info_transMat();
      void info_transRed();
      void info_basisIdx();
      void info_nf_pos();
      void info_nbIdx();
      void info_redVec();
      void info_rIdx();
      void info_pivOutIdx();

      void info_p_sol_ptr();
      void info_d_sol_ptr();
      void info_invB_ptr();
      void info_transMat_ptr();
      void info_transRed_ptr();
      void info_basisIdx_ptr();
      void info_nf_pos_ptr();
      void info_nbIdx_ptr();
      void info_redVec_ptr();

      void info_fIdx();
      void info_node();

};

class ftData
{
   private:

      int Dim;

   public:

      ftData();
      ~ftData();

      int elemNum;

      theData* cur;
      theData* parent;

      theData* limit;

      theData* head;
      theData* last;

      void create_elem(int row, int col, int termS, int polyDim);
      void add_elem();

      void mark();

      void clear();
      void clear_transMat();

      void delete_cur();
      void delete_all();

      void delete_addedElem();

      void init_ptr(){ parent = (cur = head); };

      void make_init_data(int termSumNum, int supN, int termS, int reTermS);
      void next_data(); 

      void copy(int col, theData* pre_data);
      void get_ptr(theData* pre_data);

      // iCheck -->
      void create_rIdx(int nbN, int repIdx, int* candIdx);
      void init_info();

      void get_nbIdx_rIdx ( int preNbN, int repIdx, int* candIdx, 
                            int reTermS, theData* pre_data );

      void iCopy ( int nbN, int nfN, int repIdx, 
                   int termS, int reTermS, int* candIdx, theData* pre_data );

      void iGetPtr(theData* pre_data);

  // <--

  // mCheck -->
      void output(int repIdx, int* idx2, int* nbN, int* nfN);

      void decrease_nfN();
      void copy_rIdx(theData* pre_data, int termS);
      void copy_pivOutIdx(theData* pre_data);
      void get_nf_pos(theData* pre_data, int nfN, int idx2);

      void mCopy(int nbN, int nfN, int idx2, int termS, theData* pre_data);
      void mGetPtr(theData* pre_data);

  // <--

      void put_sup(int* sup);

      void info_parent_nbN_nfN();

      void info_parent_p_sol();
      void info_parent_d_sol();
      void info_parent_invB();
      void info_parent_transMat();
      void info_parent_transRed();
      void info_parent_basisIdx();
      void info_parent_nf_pos();
      void info_parent_nbIdx();
      void info_parent_redVec();
      void info_parent_rIdx();
      void info_parent_pivOutIdx();

      void info_parent_p_sol_ptr();
      void info_parent_d_sol_ptr();
      void info_parent_invB_ptr();
      void info_parent_transMat_ptr();
      void info_parent_transRed_ptr();
      void info_parent_basisIdx_ptr();
      void info_parent_nf_pos_ptr();
      void info_parent_nbIdx_ptr();
      void info_parent_redVec_ptr();
      void info_parent_pivOutIdx_ptr();

      void info_parent();
      void info_parent_ptr();
      void info_parent_node();

      void info_cur_nbN_nfN();

      void info_cur_p_sol();
      void info_cur_d_sol();
      void info_cur_invB();
      void info_cur_transMat();
      void info_cur_transRed();
      void info_cur_basisIdx();
      void info_cur_nf_pos();
      void info_cur_nbIdx();
      void info_cur_redVec();
      void info_cur_rIdx();
      void info_cur_pivOutIdx();

      void info_cur_p_sol_ptr();
      void info_cur_d_sol_ptr();
      void info_cur_invB_ptr();
      void info_cur_transMat_ptr();
      void info_cur_transRed_ptr();
      void info_cur_basisIdx_ptr();
      void info_cur_nf_pos_ptr();
      void info_cur_nbIdx_ptr();
      void info_cur_redVec_ptr();

      void info_cur();
      void info_cur_ptr();
      void info_cur_node();

      void info_all_node();
      void info_all_cur();
      void info_all_nodeNum();

      void info_numElem();
};

class lvData
{

   private:

      int Dim;

      int length;
      int termMax;

  public:

     lvData();
     ~lvData();

     int* mRepN;
     int** mFeaIdx;
     int* mFea;

     ftData* fTest;
     ftData* Node;

     void create ( int depth, int supN, int Dim, 
                   int ori_length, int ori_termMax );

     void get_info ( int** g_mRepN, int*** g_mFeaIdx, int** g_mFea );
     void init_ptr();

     void info_mFea();
};

#endif
