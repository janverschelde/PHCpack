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

#ifndef __SIMPLEX_H
#define __SIMPLEX_H

#include "global.h"
#include "inputData.h"
#include "iTest.h"
#include "fTest.h"

class supportSet
{
   int row, col;

   public:

      supportSet();
      ~supportSet();

      double* supMat;
      double* costVec;

      void allocSupp ( dataSet& Data, int level, int num, double* lifting );
      void allocAux ( dataSet& Data );

      void supMat_in ( int rowIdx, int colIdx, double elem )
      {
         supMat[rowIdx + colIdx * row] = elem;
      };

      void supMat_neg ( int rowIdx, int colIdx )
      {
         supMat[rowIdx + colIdx * row] *= -1;
      };

      double supMat_out ( int rowIdx, int colIdx )
      {
         return (supMat[rowIdx + colIdx * row]);
      };

      double redVal ( double* d_sol, int idx, int ii )
      {
         int i;
         double val = 0;

         // cout << "<< redVal >>\n\n";
         // cout << "d_sol: ";

         for(i = 0; i < row; i++)
         {
            // cout << d_sol[i] << " ";
            val += d_sol[i] * supMat[i + ii];
         }
         // cout << "\n";
         // cout << "idx: " << idx << "\n";
         // cout << "costVec: " << costVec[idx] << "\n\n";

         return(costVec[idx] - val);
      };

      void info_sup();
      void info_costVec();
};

class simplex
{
   private:

      int Dim;  
      int supN; 
      int termSumNum;  
  
      int repIdx;

      int* candIdx;
  
      int* firIdx;
      int* termSet;    
      int* termStart;
      int* re_termStart; 

      int output;

      double mixedVol;
      int mixedCell;

      int* ip;
      double* weight;
      double* vol;

      double* eye;

      int nbN;
      int nfN;

      int artV;
      int pivOutNum;

      int frIdx;

      supportSet** Supp; 

      double** oriSupp;

      double* invB; // row oriented
      double* transMat; // row oriented
      double* transRed; 

      double* p_sol;
      double* d_sol;

      double* p1_d_sol;
      double* fst_d_sol;
  
      double* aux_cvec; 
      double* dir; 

      double* fst_redVec; 
      double* redVec;

      int* basisIdx;
      int* nf_pos;
      int* nbIdx;
      int* rIdx;

      int *pivOutList;
      int *pivOutCheck;

      double* tmp_newInvB;  
      double* tmp_transMat;  
  
      int* nIdx;  
  
      double* pre_p_sol;
      double* pre_d_sol;

      double* pre_redVec;

      int* pre_basisIdx;
      int* pre_nbIdx;
      int* pre_nf_pos;
  
      double* pre_invB;
      double* pre_transMat;
      double* pre_transRed;

  // relation table
      int checkFrIdx();
      void elimFrIdx ( int sub_pivOutIdx );

  // phase1
      void reMakeNonBasisIdx ( int reTermS );
      void reMakeNonBasisIdx_tab();

      void elimArt ( int depth, int preNbN, int termS, int reTermS, int& iter );
      void calRedCost ( int pivInIdx, double& redCost );
      int isZeroDirEle ( int termS, int idx, int preNbN, int& sub_pivInIdx );
      void IP_vec_mat();

  // reduced cost
      int reducedCost_tab_p1 ( int& enterIdx, int& sub_enterIdx,
                               double& redCost );
      int reducedCost_tab ( int& enterIdx, int& sub_enterIdx, double& redCost );

      int reducedCost_p1 ( int& enterIdx, int& sub_enterIdx, double& redCost );

      int reducedCost ( int& enterIdx, int& sub_enterIdx, double& redCost );

      int reducedCost_Bland ( int& enterIdx, int& sub_enterIdx,
                              double& redCost );

      int reducedCost_mFst ( int& enterIdx, int& sub_enterIdx, 
                             int pivOutIdx, int sub_pivOutIdx,
                             double& redCost );

      int reducedCost_iFst ( int& enterIdx, int& sub_enterIdx, 
                             int pivOutIdx, int sub_pivOutIdx, double& redCost,
                             int termS, int reTermS, int preNbN );
  
      void extend_nbIdx ( int cIdx, int pre_pivInIdx, int pre_pivOutIdx, 
                          int pre_length, int reTermS, int& cnt );

      int extend_nbIdx_comp ( int& non_basisIdx, int cIdx, int pre_pivInIdx,
                              int pre_pivOutIdx, 
                              int pre_length, int reTermS, int& cnt );

      void getIdx ( int& level, int& idx, int& idx2, int& ii, int d_nbIdx );

  // ratio test
      int ratioTest ( int redFlag, int pivInIdx, int sub_pivInIdx,
                      int& pivOutIdx, int& sub_pivOutIdx, double& theta );

      int ratioTest_artFst ( int redFlag, int pivInIdx, int sub_pivInIdx,
                             int& pivOutIdx, int& sub_pivOutIdx,
                             double& theta );

      int ratioTest_art ( int redFlag, int pivInIdx, int sub_pivInIdx,
                          int& pivOutIdx, int& sub_pivOutIdx, double& theta );

      int ratioTest_art_Bland ( int redFlag, int pivInIdx, int sub_pivInIdx,
                                int& pivOutIdx, int& sub_pivOutIdx,
                                double& theta );

      int ratioTest_frIdx ( int pivInIdx );

      void IP_mat_vec ( int pivInIdx );
      void IP_mat_vec_fst ( int pivInIdx );

  //
      void update_p1_d_sol ( int pivInIdx, int sub_pivOutIdx );

      void modify_p_sol ( int pivInIdx );
      void calElem ( int idx );

  // create new basis and nonbasis
      void createNewBandN_tab ( int pivInIdx, int sub_pivInIdx, 
                                int pivOutIdx, int sub_pivOutIdx, 
                                double theta, double redCost );

      void createNewBandN_p1 ( int pivInIdx, int sub_pivInIdx, 
                               int pivOutIdx, int sub_pivOutIdx, 
                               double theta, double redCost, int termS, 
                               int reTermS );

      void createNewBandN ( int pivInIdx, int sub_pivInIdx,
                            int pivOutIdx, int sub_pivOutIdx,
                            double theta, double redCost, int termS,
                            int reTermS );

      void createNewBandN_iFst ( int pivInIdx, int sub_pivInIdx,
                                 int pivOutIdx, int sub_pivOutIdx,
                                 double theta, double redCost,
                                 int termS, int reTermS );

      void createNewBandN_mFst ( int pivInIdx, int sub_pivInIdx, 
                                 int pivOutIdx, int sub_pivOutIdx, 
                                 double theta, double redCost, int termS,
                                 int reTermS );

      void createNewBandN_art ( int pivInIdx, int sub_pivInIdx, 
                                int pivOutIdx, int sub_pivOutIdx, 
                                double redCost, int termS, int reTermS );

  //
      void invB_in ( int rowIdx, int colIdx, double elem )
      {
         invB[colIdx + Dim * rowIdx] = elem; 
      };

      double invB_out ( int rowIdx, int colIdx )
      {
         return (invB[colIdx + Dim * rowIdx]);
      };

  //
      double transMat_out ( int rowIdx, int colIdx )
      {
         return (transMat[colIdx + Dim * rowIdx]);
      };

      void supp_in ( int lvl, int rowIdx, int colIdx, double elem )
      {
         oriSupp[lvl][rowIdx + colIdx * Dim] = elem;
      };

      double supp_out ( int lvl, int rowIdx, int colIdx )
      {
         return(oriSupp[lvl][rowIdx + colIdx * Dim]);
      };

      int isZero ( double val )
      {
         if(MINUSZERO < val && val < PLUSZERO) return (TRUE);
         else return (FALSE);
      }

  //
      void info_p_sol();
      void info_d_sol();
      void info_p1_d_sol();
      void info_invB();
      void info_transMat();
      void info_transRed();
      void info_basisIdx();
      void info_nf_pos();
      void info_nbIdx();
      void info_rIdx();
      void info_redVec();
      void info_dir();
      void info_frIdx();
      void info_candIdx();
      void info_repIdx();

      void info_oriSup();

   public:

      simplex();
      ~simplex();

  //
      double* lifting; 
 
  //
      void get_iNbN_nfN ( theData** cur, int lNbN, int lNfN );
       // initCheck, iCheck
      void get_mNbN_nfN ( theData* parent, theData** cur );  // mCheck

      void get_repIdx_candIdx ( int* ori_candIdx, int ori_repIdxa ); 
      void get_parent ( theData* parent );
      void get_cur ( theData** cur ); 

      void get_res ( ftData& iData ); // fSolLP, solLP_art
      void get_pivOutNum ( theData** cur ); 

  //
      void get_nbN_nfN ( int ori_nbN, int ori_nfN );
      void get_p_sol ( double* ori_p_sol );
      void get_d_sol ( double* ori_d_sol );
      void get_basisIdx ( int* ori_basisIdx );
      void get_nf_pos ( int* ori_nf_pos );
      void get_nbIdx ( int* ori_nbIdx );
      void get_invB ( double* invB );
      void get_frIdx ( int ori_frIdx );

      void copy_p1_d_sol ( theData* cur )
      {
         memcpy(p1_d_sol, cur->d_sol, sizeof(double) * Dim); 
      };
 
      void copy_eye(theData** cur)
      {
         memcpy((*cur)->transMat, eye, sizeof(double) * Dim * Dim);
      };

#ifdef compile4phc
      void initialize_with_lifting
             ( dataSet& Data, double* lifvals,
               int* ori_firIdx, int seedNum, int ori_output );
      /*
       * Extends the definition of allocateAndIni with lifting values
       * for each point in the supports. */
#endif
      void allocateAndIni ( dataSet& Data, int* ori_firIdx,
                            int seedNum, int ori_output );

  // For relation table
      int tSolLP ( int& iter, int mode ); 

  // For Phase-1 and 2
      int fSolLP ( int termS, int reTermS, int& iter ); 

  // iCheck
      void fstRed_candIdx ( inifData& curInif, int** mCandIdx, 
                            int& pivInIdx, int& sub_pivInIdx );

      void cal_redVec
             ( int termS, int reTermS, int fst_pivInIdx, theData** cur );
      double put_redCost ( int fst_pivInIdx )
      {
         return (fst_redVec[fst_pivInIdx] - fst_redVec[repIdx]);
      };

  // iCheck_art

      int solLP_art ( int depth, int idx_one, int fst_pivIn, int preNbN, 
                      int termS, int reTermS, int& iter ); 

      int solLP_art_Bland ( int pivInIdx, int sub_pivInIdx, 
                            int pivOutIdx, int sub_pivOutIdx,
                            int redFlag, double theta, double redCost, 
                            int termS, int reTermS, int& iter );

  // For mLP
      int solLP ( int depth, int fst_pivInIdx, int fst_sub_pivInIdx,
                  double fst_redCost, int mode, int termS, int reTermS,
                  int preNbN, int& iter ); 

      int solLP_Bland ( int pivInIdx, int sub_pivInIdx,
                        int pivOutIdx, int sub_pivOutIdx,
                        int redFlag, double theta, double redCost,
                        int termS, int reTermS, int& iter );

      int initIter ( int mode, int fst_pivInIdx, int fst_sub_pivInIdx,
                     double fst_redCost, int& redFlag, int& pivInIdx,
                     int& sub_pivInIdx, int& pivOutIdx, int& sub_pivOutIdx,
                     double& theta, double& redCost, 
                     int termS, int reTermS, int preNbN );

  //
      void calMixedVol ( lvData* lv, int* sp, int supN );

      double lu ( int n, double* a );
      double matinv ( int n, double* a, double* a_inv );

  //
      double put_elem_supp ( int lvl, int idx, int row, int col )
      {
         return(Supp[lvl][idx].supMat_out(row, col));
      };

      void mult_elem_supp ( int lvl, int idx, int row, int col )
      {
         Supp[lvl][idx].supMat_neg(row, col);
      };

  //
      void check_dirRed(theData* parent, int depth);
      void dbg_dirRed(theData* parent, inifData* nextInif, int depth);

  //
      void info_mv();
      void info_allSup();
      void info_allCostVec();
      void info_lifting();
      void info_simplexData();
};

#endif
