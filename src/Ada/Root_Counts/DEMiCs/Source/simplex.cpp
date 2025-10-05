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

/*
    To compile for use with PHCpack, add
    -Dcompile4phc=1
    to the compiler statement.
    By doing so, the code in the #ifdef compile4phc
    statements is compiled.
 */

#include "simplex.h"

#ifdef compile4phc
#include "outputData.h"
#endif

// JV added 23 April 2018
#include <iomanip>
// JV addeed 11 May 2018
#include <sstream>

supportSet::supportSet()
{
   row = 0;
   col = 0;

   supMat = NULL;
   costVec = NULL;
}

supportSet::~supportSet()
{
   delete [] supMat;
   delete [] costVec;
}

void supportSet::allocSupp
 ( dataSet& Data, int level, int num, double* lifting )
{
   int i, j, cnt = 0;

   // cout << "<< allocSupp >>\n";

   row = Data.Dim;
   col = Data.termSet[level];

   supMat = new double [row * col * sizeof(double)];
   assert(supMat);
   memset(supMat, 0, row * col * sizeof(double));

   costVec = new double [col - 1];
   assert(costVec);
   memset(costVec, 0, (col - 1) * sizeof(double));

   for(j = 0; j < col; j++)
   {
      if(j != num)
      {
         for(i = 0; i < row; i++)
         {
             supMat_in(i, cnt, Data.support_out(Data.termStart[level] + j, i) 
                - Data.support_out(Data.termStart[level] + num, i));
         }
         costVec[cnt] = lifting[Data.termStart[level] + j]
                      - lifting[Data.termStart[level] + num];
         cnt++;
      }
   }
   col--;
}

void supportSet::allocAux ( dataSet& Data )
{
   int i;
  
   row = (col = Data.Dim);

   supMat = new double [row * col * sizeof(double)];
   assert(supMat);
   memset(supMat, 0, row * col * sizeof(double));

   costVec = new double [col];
   assert(costVec);

   for(i = 0; i < col; i++)
   {
      supMat_in(i, i, 1);
      costVec[i] = 1;
   }
}

void supportSet::info_sup()
{
   int i, j;

   for(j = 0; j < row; j++)
   {
      for(i = 0; i < col; i++)
      {
         cout.width(3);
         cout << supMat_out(j, i) << " ";
      }
      cout << "\n";
   }
}

void supportSet::info_costVec()
{
   int i;

   for(i = 0; i < col; i++)
   {
      cout.width(3);
      cout << costVec[i] << " ";
   }
   cout << "\n";
}

simplex::simplex ( void )
{
   Dim = 0; 
   supN = 0; 
   termSumNum = 0; 

   repIdx = 0;

   output = 0;

   mixedVol = 0;
   mixedCell = 0;

   ip = NULL;
   weight = NULL;
   vol = NULL;

   eye = NULL;

   nbN = 0; 
   nfN = 0; 

   artV = 0; 
   pivOutNum = 0;

   frIdx = 0;

   Supp = NULL;
   lifting = NULL;
   oriSupp = NULL;

   candIdx = NULL;

   firIdx = NULL;
   termSet = NULL; 
   termStart = NULL;
   re_termStart = NULL; 

   invB = NULL; 
   transMat = NULL;
   transRed = NULL;
  
   p_sol = NULL; 
   d_sol = NULL; 

   p1_d_sol = NULL;
   fst_d_sol = NULL; 
  
   aux_cvec = NULL; 
   dir = NULL; 

   fst_redVec = NULL; 
   redVec = NULL; 

   basisIdx = NULL; 
   nf_pos = NULL; 
   nbIdx = NULL;  
   rIdx = NULL; 

   pivOutList = NULL;
   pivOutCheck = NULL;

   tmp_newInvB = NULL;
   tmp_transMat = NULL;

   nIdx = NULL; 

   pre_p_sol = NULL;
   pre_d_sol = NULL;

   pre_redVec = NULL;

   pre_basisIdx = NULL;
   pre_nbIdx = NULL;
   pre_nf_pos = NULL;

   pre_invB = NULL;
   pre_transMat = NULL;
   pre_transRed = NULL;
}

simplex::~simplex ( void )
{
   int i;

   delete [] re_termStart;
   delete [] p1_d_sol;
   delete [] ip;
   delete [] weight;
   delete [] vol;
   delete [] eye;
   delete [] lifting;
   delete [] fst_d_sol;
   delete [] aux_cvec;
   delete [] dir;
   delete [] fst_redVec;
   delete [] tmp_newInvB;
   delete [] tmp_transMat;
   delete [] nIdx;

   if(Supp)
   {
      for(i = 0; i <= supN; i++) delete [] Supp[i];

      delete [] Supp;
      Supp = NULL;
   }
   if(oriSupp)
   {
      for(i = 0; i < supN; i++) delete [] oriSupp[i];

      delete [] oriSupp;
      oriSupp = NULL;
   }
}

#ifdef compile4phc
void simplex::initialize_with_lifting
 ( dataSet& Data, double* lifvals,
   int* ori_firIdx, int seedNum, int ori_output )
{
   int i, j, k, cnt;
   int tmp_nDim, tmp_mDim;

   double Rand_Max;
  
   Dim = (tmp_nDim = (tmp_mDim = Data.Dim));

   supN = Data.supN;
   termSumNum = Data.termSumNum;

   tmp_nDim += termSumNum - supN;
  
   firIdx = ori_firIdx;
   output = ori_output;

   termSet = Data.termSet;    
   termStart = Data.termStart;
  
   re_termStart = new int [supN + 1];
   assert(re_termStart);
   re_termStart[0] = 0;

   p1_d_sol = new double [Dim];
   assert(p1_d_sol);
   memset(p1_d_sol, 0, sizeof(double) * Dim);

   ip = new int [Dim];
   assert(ip);
   memset(ip, 0, sizeof(int) * Dim);

   weight = new double [Dim];
   assert(weight);
  
   vol = new double [Dim * Dim];
   assert(vol);
   memset(vol, 0, Dim * Dim * sizeof(double));

   eye = new double [Dim * Dim];
   assert(eye);
   memset(eye, 0, Dim * Dim * sizeof(double));
  
   lifting = new double [termSumNum];
   assert(lifting);
  
   oriSupp = new double* [supN];
   assert(oriSupp);

   for(i = 0; i < supN; i++)
   {
      oriSupp[i] = new double [Dim * termSet[i]];
      assert(oriSupp[i]);
      memset(oriSupp[i], 0, sizeof(double) * Dim * termSet[i]);

      for(j = 0; j < Dim; j++)
      {
         for(k = termStart[i]; k < termStart[i] + termSet[i]; k++)
         {
            supp_in(i, j, k - termStart[i], Data.support_out(k, j));
         }
      }
   }
   fst_d_sol = new double [tmp_mDim];
   assert(fst_d_sol);
   memset(fst_d_sol, 0, tmp_mDim * sizeof(double));
  
   aux_cvec = new double [tmp_nDim];
   assert(aux_cvec);
   memset(aux_cvec, 0, tmp_nDim * sizeof(double));

   dir = new double [tmp_mDim];
   assert(dir);
   memset(dir, 0, tmp_mDim * sizeof(double));
  
   fst_redVec = new double [Data.termMax];
   assert(fst_redVec);
   memset(fst_redVec, 0, Data.termMax * sizeof(double));

   tmp_newInvB = new double [Dim];
   assert(tmp_newInvB);

   tmp_transMat = new double [Dim];
   assert(tmp_transMat);
  
   nIdx = new int [2 * tmp_nDim];
   assert(nIdx);
  
   srand(seedNum); // srandom(seedNum);
   Rand_Max = 2;  

   for(i = 1; i < 30; i++) Rand_Max *= 2;
   Rand_Max = (Rand_Max - 1) * 2 + 1;

   for(i = 0; i < termSumNum; i++) // copy the lifting values
   {
      lifting[i] =  lifvals[i]; // (double) random() / Rand_Max*10;
   }
   Supp = new supportSet* [supN + 1];
   assert(Supp);
  
   for(i = 0; i < supN; i++)
   {
      re_termStart[i + 1] = Data.termStart[i + 1] - i - 1;
 
      Supp[i] = new supportSet [termSet[i]];
      assert(Supp[i]);
   }
   for(i = 0; i < supN; i++)
   {
      for(j = 0; j < termSet[i]; j++)
      {
         Supp[i][j].allocSupp(Data, i, j, lifting);
      }
   }
   Supp[supN] = new supportSet [1];
   assert(Supp[supN]);

   Supp[supN][0].allocAux(Data);

   cnt = 0;
   for(i = 0; i < supN; i++)
   {
      for(j = 0; j < termSet[i] - 1; j++)
      {
         nIdx[2 * cnt] = i;
         nIdx[2 * cnt + 1] = j;
         cnt++;
      }
   }
   for(i = 0; i < Dim; i++)
   {
      aux_cvec[termSumNum - supN + i] = 1;

      nIdx[2 * cnt] = supN;
      nIdx[2 * cnt + 1] = 0;

      eye[i * ( Dim + 1)] = 1;

      cnt++;
   }

   int fail = demics_allocate_lifting(Data.supN,Data.termSet);
   cnt = 0;
   for(i=0; i<Data.supN; i++)
      for(j=0; j<Data.termSet[i]; j++)
         fail = demics_assign_lifting(i,j,lifting[cnt++]);

   if(output)
   {
      cout << "----------------------------------\n";

      cout << "* Seed number = "  << seedNum << endl;
      cout << "* Lifting values for elements in each support set"
           << endl << endl;

      cout << fixed << setprecision(16); // JV Mon 23 Apr 2018

      cnt = 0;
      for(i = 0; i < supN; i++)
      {
         cout << "S";

         cout.width(3);
         cout.setf(ios::left);
         cout <<  i + 1 << " : ";

         for(j = termStart[i]; j < termStart[i] + termSet[i]; j++)
         {
            cout.width(10);
            cout.setf(ios::left);
            cout << lifting[cnt] << " ";
            cnt++;
         }
         cout << endl;
      }
      cout << "----------------------------------\n";
      cout << endl << endl;
   }
}
#endif

void simplex::allocateAndIni
 ( dataSet& Data, int* ori_firIdx, int seedNum, int ori_output )
{
   int i, j, k, cnt;
   int tmp_nDim, tmp_mDim;

   double Rand_Max;
  
   Dim = (tmp_nDim = (tmp_mDim = Data.Dim));

   supN = Data.supN;
   termSumNum = Data.termSumNum;

   tmp_nDim += termSumNum - supN;
  
   firIdx = ori_firIdx;
   output = ori_output;

   termSet = Data.termSet;    
   termStart = Data.termStart;
  
   re_termStart = new int [supN + 1];
   assert(re_termStart);
   re_termStart[0] = 0;

   p1_d_sol = new double [Dim];
   assert(p1_d_sol);
   memset(p1_d_sol, 0, sizeof(double) * Dim);

   ip = new int [Dim];
   assert(ip);
   memset(ip, 0, sizeof(int) * Dim);

   weight = new double [Dim];
   assert(weight);
  
   vol = new double [Dim * Dim];
   assert(vol);
   memset(vol, 0, Dim * Dim * sizeof(double));

   eye = new double [Dim * Dim];
   assert(eye);
   memset(eye, 0, Dim * Dim * sizeof(double));
  
   lifting = new double [termSumNum];
   assert(lifting);
  
   oriSupp = new double* [supN];
   assert(oriSupp);

   for(i = 0; i < supN; i++)
   {
      oriSupp[i] = new double [Dim * termSet[i]];
      assert(oriSupp[i]);
      memset(oriSupp[i], 0, sizeof(double) * Dim * termSet[i]);

      for(j = 0; j < Dim; j++)
      {
         for(k = termStart[i]; k < termStart[i] + termSet[i]; k++)
         {
            supp_in(i, j, k - termStart[i], Data.support_out(k, j));
         }
      }
   }
   fst_d_sol = new double [tmp_mDim];
   assert(fst_d_sol);
   memset(fst_d_sol, 0, tmp_mDim * sizeof(double));
  
   aux_cvec = new double [tmp_nDim];
   assert(aux_cvec);
   memset(aux_cvec, 0, tmp_nDim * sizeof(double));

   dir = new double [tmp_mDim];
   assert(dir);
   memset(dir, 0, tmp_mDim * sizeof(double));
  
   fst_redVec = new double [Data.termMax];
   assert(fst_redVec);
   memset(fst_redVec, 0, Data.termMax * sizeof(double));

   tmp_newInvB = new double [Dim];
   assert(tmp_newInvB);

   tmp_transMat = new double [Dim];
   assert(tmp_transMat);
  
   nIdx = new int [2 * tmp_nDim];
   assert(nIdx);
  
   srand(seedNum); // srandom(seedNum);
   Rand_Max = 2;  

   for(i = 1; i < 30; i++) Rand_Max *= 2;
   Rand_Max = (Rand_Max - 1) * 2 + 1;

   for(i = 0; i < termSumNum; i++)
   {
      // lifting[i] = (double) random() / Rand_Max*10;
      lifting[i] = (double) rand() / Rand_Max*10;
   }
   Supp = new supportSet* [supN + 1];
   assert(Supp);
  
   for(i = 0; i < supN; i++)
   {
      re_termStart[i + 1] = Data.termStart[i + 1] - i - 1;
 
      Supp[i] = new supportSet [termSet[i]];
      assert(Supp[i]);
   }
   for(i = 0; i < supN; i++)
   {
      for(j = 0; j < termSet[i]; j++)
      {
         Supp[i][j].allocSupp(Data, i, j, lifting);
      }
   }
   Supp[supN] = new supportSet [1];
   assert(Supp[supN]);

   Supp[supN][0].allocAux(Data);

   cnt = 0;
   for(i = 0; i < supN; i++)
   {
      for(j = 0; j < termSet[i] - 1; j++)
      {
         nIdx[2 * cnt] = i;
         nIdx[2 * cnt + 1] = j;
         cnt++;
      }
   }
   for(i = 0; i < Dim; i++)
   {
      aux_cvec[termSumNum - supN + i] = 1;

      nIdx[2 * cnt] = supN;
      nIdx[2 * cnt + 1] = 0;

      eye[i * ( Dim + 1)] = 1;

      cnt++;
   }

#ifdef compile4phc
   int fail = demics_allocate_lifting(Data.supN,Data.termSet);
   cnt = 0;
   for(i=0; i<Data.supN; i++)
      for(j=0; j<Data.termSet[i]; j++)
         fail = demics_assign_lifting(i,j,lifting[cnt++]);
#endif
   if(output)
   {
      cout << "----------------------------------\n";

      cout << "* Seed number = "  << seedNum << endl;
      cout << "* Lifting values for elements in each support set"
           << endl << endl;

      cout << fixed << setprecision(16); // JV Mon 23 Apr 2018

      cnt = 0;
      for(i = 0; i < supN; i++)
      {
         cout << "S";

         cout.width(3);
         cout.setf(ios::left);
         cout <<  i + 1 << " : ";

         for(j = termStart[i]; j < termStart[i] + termSet[i]; j++)
         {
            cout.width(10);
            cout.setf(ios::left);
            cout << lifting[cnt] << " ";
            cnt++;
         }
         cout << endl;
      }
      cout << "----------------------------------\n";
      cout << endl << endl;
   }
}

void simplex::reMakeNonBasisIdx ( int reTermS )
{
   int i, cnt = 0, artNbIdxN = 0;
   int tmp_nbIdx;

   // cout << "<< reMakeNonBasisIdx >>\n\n";

   for(i = 0; i < nbN - Dim; i++)
   {
      // cout << "nbIdx: " << nbIdx[i] << "\n";

      if(nbIdx[i] < termSumNum - supN)
      {
         tmp_nbIdx = nbIdx[i];
    
         nbIdx[cnt] = tmp_nbIdx;
         rIdx[tmp_nbIdx - reTermS] = -1 * (cnt + 1);

         cnt++;
      }
      else
      {
         nbIdx[i] = 0;
         artNbIdxN++;
      }
   }
   // cout << "\n";

   if(artNbIdxN == Dim)
      artV = 0;
   else
      artV = 1;

   // memcpy(d_sol, p1_d_sol, sizeof(double) * Dim); 

   IP_vec_mat();

   nbN = Dim + cnt;

   // info_d_sol();
   // info_p1_d_sol();
}

void simplex::reMakeNonBasisIdx_tab()
{
   int i, cnt = 0;
   int tmp_nbIdx;

   // cout << "<< reMakeNonBasisIdx_tab >>\n\n";

   for(i = 0; i < nbN - Dim; i++)
   {
      // cout << "nbIdx: " << nbIdx[i] << "\n";

      if(nbIdx[i] < termSumNum - supN)
      {
         tmp_nbIdx = nbIdx[i];
    
         nbIdx[cnt] = tmp_nbIdx;

         cnt++;
      }
      else
      {
         nbIdx[i] = 0;
      }
   }
   IP_vec_mat();
   nbN = Dim + cnt;
}

void simplex::elimArt
 ( int depth, int preNbN, int termS, int reTermS, int& iter )
{
   int i, sub_pivInIdx, flag;
   int cnt = 0, cnt_t = 0;

   double redCost;

   // cout << "<< elimArt >>\n\n";

   for(i = 0; i < Dim; i++)
   {
      if(basisIdx[i] >= termSumNum - supN)
      {
         flag = isZeroDirEle(termS, i, preNbN, sub_pivInIdx);

         cnt++;

         // True  -- val > 0
         // False -- val =0

         if(flag == TRUE)
         {
	    calRedCost(nbIdx[sub_pivInIdx], redCost);

            IP_mat_vec(nbIdx[sub_pivInIdx]);

            createNewBandN_art(nbIdx[sub_pivInIdx], sub_pivInIdx, 
                               basisIdx[i], i, redCost, termS, reTermS);

            iter++;
            cnt_t++;
         }
      }
   }
   if(cnt == cnt_t)
      artV = 0;
   else
      artV = 1;
}

void simplex::calRedCost ( int pivInIdx, double& redCost )
{
   int i, ii;
   int level, idx, idx2;
   double val = 0;

   ii = 2 * pivInIdx;
   level = nIdx[ii];
   idx = nIdx[ii + 1];

   idx2 = firIdx[level];
   ii = Dim * idx;
      
   /*
     cout << "level: " << level << "\n";
     cout << "Idx2: " << idx2 << "\n";
     cout << "ii: " << ii << "\n\n";
   */

   for(i = 0; i < Dim; i++)
   {
      val += d_sol[i] * Supp[level][idx2].supMat[i + ii];
   }
   redCost = Supp[level][idx2].costVec[idx] - val;
}

int simplex::isZeroDirEle ( int termS, int idx, int preNbN, int& sub_pivInIdx )
{
   int i, j, ii, iii, idx2, level;
   double val;

   // cout << "<< isZeroDirEle >>\n\n";

   for(j = 0; j < termS - 1; j++)
   {
      // cout << "nbIdx: " << nbIdx[preNbN - Dim + j] << "\n";

      ii = 2 * nbIdx[preNbN - Dim + j];
      level = nIdx[ii];

      idx2 = firIdx[level];
      iii = Dim * nIdx[ii + 1];

      val = 0;
      ii = Dim * idx;
    
      for(i = 0; i < Dim; i++)
      {
         val += invB[ii + i] * Supp[level][idx2].supMat[i + iii];
      }
      // cout << "val: " << val << "\n\n";
    
      if(val < MINUSZERO || val > PLUSZERO)
      {
         sub_pivInIdx = preNbN - Dim + j;
         return (TRUE);
      }
   }
   // cout << "\n";

   return (FALSE);
}

void simplex::IP_vec_mat()
{
   int i, j, ii;
   int level, idx, idx2;

   // cout << "<< IP_vec_mat >> ";

   for(j = 0; j < Dim; j++)
   {
      if(basisIdx[j] < termSumNum - supN)
      {
         ii = 2 * basisIdx[j];
         level = nIdx[ii];
         idx = nIdx[ii + 1];

         idx2 = firIdx[level];
         ii = Dim * j;

         for(i = 0; i < Dim; i++)
         {
            d_sol[i] += invB[ii + i] * Supp[level][idx2].costVec[idx];
         }
      }
   }
}

int simplex::reducedCost ( int& pivInIdx, int& sub_pivInIdx, double& redCost )
{
   int j, ii, tmp_non_basisIdx, flag = OPT;
   int level, idx, idx2;
   double val, tmp_redCost;

   // cout << "<< reducedCost >>\n\n";

   redCost = 1.0E-8;

   for(j = 0; j < nbN - Dim; j++)
   {
      val = 0;
      tmp_non_basisIdx = nbIdx[j];
      // cout << "---- nbIdx: " << tmp_non_basisIdx << " ----\n";

      getIdx(level, idx, idx2, ii, 2 * tmp_non_basisIdx);

    /*
      cout << "level: " << level << "\n";
      cout << "Idx2: " << idx2 << "\n";
      cout << "ii: " << ii << "\n\n";
     */

      (redVec[tmp_non_basisIdx] 
         = (tmp_redCost = Supp[level][idx2].redVal(d_sol, idx, ii)));

      // cout << "redVal: " << tmp_redCost << "\n\n";

      if((tmp_redCost < MINUSZERO) && (fabs(tmp_redCost) > fabs(redCost)))
      {
         redCost = tmp_redCost;

         pivInIdx = tmp_non_basisIdx;
         sub_pivInIdx = j;

         flag = POSTHETA;
      }
   }
   return (flag);
}

int simplex::reducedCost_Bland
 ( int& pivInIdx, int& sub_pivInIdx, double& redCost )
{
   int j, ii, tmp_non_basisIdx, flag = OPT;
   int level, idx, idx2;

   double val, tmp_redCost;

   // cout << "<< reducedCost_Bland >>\n\n";

   pivInIdx = BIGINT;

   for(j = 0; j < nbN - Dim; j++)
   {
      val = 0;
      tmp_non_basisIdx = nbIdx[j];

      // cout << "---- nbIdx: " << tmp_non_basisIdx << " ----\n";
    
      getIdx(level, idx, idx2, ii, 2 * tmp_non_basisIdx);

    /*
      cout << "level: " << level << "\n";
      cout << "Idx2: " << idx2 << "\n";
      cout << "ii: " << ii << "\n\n";
     */

      (redVec[tmp_non_basisIdx]
         = (tmp_redCost = Supp[level][idx2].redVal(d_sol, idx, ii)));

      // cout << "redVal: " << tmp_redCost << "\n\n";
    
      if((tmp_redCost < MINUSZERO) && (tmp_non_basisIdx < pivInIdx))
      {
         redCost = tmp_redCost;

         pivInIdx = tmp_non_basisIdx;
         sub_pivInIdx = j;

         flag = POSTHETA;
      }
   }
   return (flag);
}

int simplex::reducedCost_tab_p1
 ( int& pivInIdx, int& sub_pivInIdx, double& redCost )
{
   int i, j, ii, tmp_non_basisIdx, flag = OPT;
   int level, idx, idx2;
   double val, tmp_redCost;
  
   // cout << "<< reducedCost_tab_p1 >>\n\n";

   redCost = 1.0E-8;

   for(j = 0; j < nbN - Dim; j++)
   {
      val = 0;
      tmp_non_basisIdx = nbIdx[j];
    
      // cout << "---- nbIdx: " << tmp_non_basisIdx << " ----\n";

      getIdx(level, idx, idx2, ii, 2 * tmp_non_basisIdx);

      // cout << "level: " << level << "\n";
      // cout << "Idx2: " << idx2 << "\n";
      // cout << "ii: " << ii << "\n\n";
      // cout << "Support: ";

      for(i = 0; i < Dim; i++)
      {
         val += d_sol[i] * Supp[level][idx2].supMat[i + ii];

         // cout << Supp[level][idx2].supMat[i + ii] << " ";
      }
      // cout << "\n";
      // cout << "aux_cvec: " << aux_cvec[tmp_non_basisIdx] << "\n";
      // cout << "val: " << val << "\n";

      tmp_redCost = aux_cvec[tmp_non_basisIdx] - val;
    
      // cout << " redVal:" << tmp_redCost << "\n\n";
    
      if((tmp_redCost < MINUSZERO) && (fabs(tmp_redCost) > fabs(redCost)))
      {
         redCost = tmp_redCost;

         pivInIdx = tmp_non_basisIdx;
         sub_pivInIdx = j;
 
         flag = POSTHETA;
      }
   }
   return (flag);
}

int simplex::reducedCost_tab
 ( int& pivInIdx, int& sub_pivInIdx, double& redCost )
{
   int j, ii, tmp_non_basisIdx, flag = OPT;
   int level, idx, idx2;
   double val, tmp_redCost;
  
   // cout << "<< reducedCost_tab >>\n\n";

   redCost = 1.0E-8;

   for(j = 0; j < nbN - Dim; j++)
   {
      val = 0;
      tmp_non_basisIdx = nbIdx[j];

      // cout << "---- nbIdx: " << tmp_non_basisIdx << " ----\n";
    
      getIdx(level, idx, idx2, ii, 2 * tmp_non_basisIdx);

    /*
      cout << "level: " << level << "\n";
      cout << "Idx2: " << idx2 << "\n";
      cout << "ii: " << ii << "\n\n";
     */

      tmp_redCost = Supp[level][idx2].redVal(d_sol, idx, ii);

      // cout << "redVal: " << tmp_redCost << "\n\n";

      if((tmp_redCost < MINUSZERO) && (fabs(tmp_redCost) > fabs(redCost)))
      {
         redCost = tmp_redCost;

         pivInIdx = tmp_non_basisIdx;
         sub_pivInIdx = j;

         flag = POSTHETA;
      }
   }
   return (flag);
}

int simplex::reducedCost_p1
 ( int& pivInIdx, int& sub_pivInIdx, double& redCost )
{
   int i, j, ii, tmp_non_basisIdx, flag = OPT;
   int level, idx, idx2;
   double val, tmp_redCost;
  
   // cout << "<< reducedCost_p1 >>\n\n";

   redCost = 1.0E-8;

   for(j = 0; j < nbN - Dim; j++)
   {
      val = 0;
      tmp_non_basisIdx = nbIdx[j];
    
      // cout << "---- nbIdx: " << tmp_non_basisIdx << " ----\n";

      getIdx(level, idx, idx2, ii, 2 * tmp_non_basisIdx);

    /*
      cout << "level: " << level << "\n";
      cout << "Idx2: " << idx2 << "\n";
      cout << "ii: " << ii << "\n\n";
     */

      // cout << "Support: ";
      for(i = 0; i < Dim; i++)
      {
         val += d_sol[i] * Supp[level][idx2].supMat[i + ii];

         // cout << Supp[level][idx2].supMat[i + ii] << " ";
      }
      // cout << "\n";
      // cout << "aux_cvec: " << aux_cvec[tmp_non_basisIdx] << "\n";
      // cout << "val: " << val << "\n";

      (redVec[tmp_non_basisIdx]
         = (tmp_redCost = aux_cvec[tmp_non_basisIdx] - val));
    
      // cout << " redVal:" << tmp_redCost << "\n\n";

      if((tmp_redCost < MINUSZERO) && (fabs(tmp_redCost) > fabs(redCost)))
      {
         redCost = tmp_redCost;

         pivInIdx = tmp_non_basisIdx;
         sub_pivInIdx = j;

         flag = POSTHETA;
      }
   }
   return (flag);
}

int simplex::reducedCost_mFst
 ( int& pivInIdx, int& sub_pivInIdx, int pivOutIdx, 
   int sub_pivOutIdx, double& redCost)
{
   int j, ii, tmp_non_basisIdx, flag = OPT;
   int level, idx, idx2;
   int pre_pivOutIdx, pre_sub_pivInIdx;

   double val, tmp_redCost;
  
   // cout << "<< reducedCost_mFst >>\n\n";
  
   pre_pivOutIdx = pivOutIdx;
   pre_sub_pivInIdx = sub_pivInIdx;

   redCost = 1.0E-8;

   for(j = 0; j < nbN - Dim; j++)
   {
      if(j == pre_sub_pivInIdx)
      {
         tmp_non_basisIdx = (nbIdx[j] = pre_pivOutIdx);
      }
      else
      {
         tmp_non_basisIdx = (nbIdx[j] = pre_nbIdx[j]);
      }
      val = 0;

      getIdx(level, idx, idx2, ii, 2 * tmp_non_basisIdx);

    /*
      cout << "level: " << level << "\n";
      cout << "Idx: " << idx << "\n";
      cout << "Idx2: " << idx2 << "\n";
      cout << "ii: " << ii << "\n\n";
     */

      (redVec[tmp_non_basisIdx]
         = (tmp_redCost = Supp[level][idx2].redVal(d_sol, idx, ii)));

      // cout << "---- nbIdx: " << tmp_non_basisIdx << " ----\n";
      // cout << "redCost: " << tmp_redCost << "\n\n";

      if((tmp_redCost < MINUSZERO) && (fabs(tmp_redCost) > fabs(redCost)))
      {
         redCost = tmp_redCost;
         pivInIdx = tmp_non_basisIdx;
         sub_pivInIdx = j;
         flag = POSTHETA;
      }
   }
   return (flag);
}

int simplex::reducedCost_iFst
 ( int& pivInIdx, int& sub_pivInIdx, int pivOutIdx, 
   int sub_pivOutIdx, double& redCost, 
   int termS, int reTermS, int preNbN )
{
   int i, j, ii, non_basisIdx, flag = OPT;
   int level, idx, idx2;
   int pre_pivOutIdx, pre_pivInIdx, pre_sub_pivInIdx;
   int length, pre_length, candNum, cnt;

   double tmp_redCost;
  
   // cout << "<< reducedCost_iFst >>\n\n";

   length = nbN - Dim;
   pre_length = preNbN - Dim;

   pre_pivInIdx = pivInIdx;
   pre_sub_pivInIdx = sub_pivInIdx;

   pre_pivOutIdx = pivOutIdx;

   // cout << "pivOutIdx: " << pivOutIdx << "\n";
   // cout << "sub_pivOutIdx: " << sub_pivOutIdx << "\n\n";

   memcpy(nbIdx, pre_nbIdx, sizeof(int) * pre_length);

   candNum = nbN - preNbN + 1;

   cnt = 0;
   for(i = 0; i < candNum; i++)
   {
      extend_nbIdx(candIdx[i + 1], pre_pivInIdx,
                   pre_pivOutIdx, pre_length, reTermS, cnt);
   }
   // cout << "\n";
   // info_nbIdx();

   redCost = 1.0E-8;

   for(j = 0; j < length; j++)
   {
      non_basisIdx = nbIdx[j];

      getIdx(level, idx, idx2, ii, 2 * non_basisIdx);

    /*
      cout << "level: " << level << "\n";
      cout << "Idx: " << idx << "\n";
      cout << "Idx2: " << idx2 << "\n";
      cout << "ii: " << ii << "\n\n";
     */
      (redVec[non_basisIdx]
         = (tmp_redCost = Supp[level][idx2].redVal(d_sol, idx, ii)));

      // cout << "---- nbIdx: " << non_basisIdx << " ----\n";
      // cout << "redCost: " << tmp_redCost << "\n\n";

      if((tmp_redCost < MINUSZERO) && (fabs(tmp_redCost) > fabs(redCost)))
      {
         redCost = tmp_redCost;
         pivInIdx = non_basisIdx;
         sub_pivInIdx = j;

         flag = POSTHETA;
      }
   }

  /*
    candNum = length - pre_length;

    cnt = 0;

    for(j = 0; j < candNum + 1; j++)
    {
       eFlag = extend_nbIdx_comp(non_basisIdx, candIdx[j + 1], pre_pivInIdx,
                                 pre_pivOutIdx, pre_length, reTermS, cnt);

       // cout << "eFlag: " << eFlag << "\n";
       // cout << "non_basisIdx: " << non_basisIdx << "\n\n";

       if(eFlag == TRUE)
       {
          getIdx(level, idx, idx2, ii, 2 * non_basisIdx);

          // cout << "level: " << level << "\n";
          // cout << "Idx: " << idx << "\n";
          // cout << "Idx2: " << idx2 << "\n";
          // cout << "ii: " << ii << "\n\n";

          (redVec[non_basisIdx]
             = (tmp_redCost = Supp[level][idx2].redVal(d_sol, idx, ii)));

          // cout << "---- nbIdx: " << non_basisIdx << " ----\n";
          // cout << "redCost: " << tmp_redCost << "\n\n";

          if((tmp_redCost < MINUSZERO) && (fabs(tmp_redCost) > fabs(redCost)))
          {
             redCost = tmp_redCost;
             pivInIdx = non_basisIdx;
             sub_pivInIdx = pre_length + cnt - 1;
 
             flag = POSTHETA;
          }
       }
    }
   */
   return (flag);
}

void simplex::extend_nbIdx
 ( int cIdx, int pre_pivInIdx, int pre_pivOutIdx,
   int pre_length, int reTermS, int& cnt )
{
   int idx;

   // cout << "idx: " << idx << "\n";

   if(repIdx > cIdx)
   {
      idx = cIdx + reTermS;
    
      if(idx == pre_pivInIdx)
      {
         nbIdx[pre_length + cnt] = pre_pivOutIdx;
      }
      else
      {
         nbIdx[pre_length + cnt] = cIdx + reTermS;
      }
      // cout << "nbIdx: " << nbIdx[pre_length + cnt] << "\n";
    
      cnt++;
   }
   else if(repIdx < cIdx)
   {
      idx = cIdx + reTermS - 1;

      if(idx == pre_pivInIdx)
      {
         nbIdx[pre_length + cnt] = pre_pivOutIdx;
      }
      else
      {
         nbIdx[pre_length + cnt] = cIdx + reTermS - 1;
      }
      // cout << "nbIdx: " << nbIdx[pre_length + cnt] << "\n";
      
      cnt++;
   }
}

int simplex::extend_nbIdx_comp
 ( int& non_basisIdx, int cIdx, int pre_pivInIdx, int pre_pivOutIdx,
   int pre_length, int reTermS, int& cnt )
{
   int flag = TRUE, idx;

   // cout << "idx: " << idx << "\n";

   if(repIdx > cIdx)
   {
      idx = cIdx + reTermS;
    
      if(idx == pre_pivInIdx)
      {
         non_basisIdx = (nbIdx[pre_length + cnt] = pre_pivOutIdx);
      }
      else
      {
         non_basisIdx = (nbIdx[pre_length + cnt] = cIdx + reTermS);
      }
      // cout << "nbIdx: " << nbIdx[pre_length + cnt] << "\n";
    
      cnt++;
   }
   else if(repIdx < cIdx)
   {
      idx = cIdx + reTermS - 1;

      if(idx == pre_pivInIdx)
      {
         non_basisIdx = (nbIdx[pre_length + cnt] = pre_pivOutIdx);
      }
      else
      {
         non_basisIdx = (nbIdx[pre_length + cnt] = cIdx + reTermS - 1);
      }
      // cout << "nbIdx: " << nbIdx[pre_length + cnt] << "\n";

      cnt++;
   }
   else
      flag = FALSE;

   return (flag);
}

inline void simplex::getIdx
 ( int& level, int& idx, int& idx2, int& ii, int d_nbIdx )
{
   level = nIdx[d_nbIdx];
   d_nbIdx++;

   idx = nIdx[d_nbIdx];

   idx2 = firIdx[level];
   ii = Dim * idx;
}

int simplex::ratioTest
 ( int redFlag, int pivInIdx, int sub_pivInIdx, int& pivOutIdx, 
   int& sub_pivOutIdx, double& theta )
{
   int i, flag = 0, nfPos, checker, nonNegVarNum;
   double tmp_theta;
  
   // cout << "<< Ratio Test >>\n\n";

   IP_mat_vec(pivInIdx);

   nonNegVarNum = 0;
   checker = 0;

   switch(redFlag)
   {
      case POSTHETA:
  
      theta = (tmp_theta = SMALLDOUBLE);
      
      for(i = 0; i < nfN; i++)
      {
         nfPos = nf_pos[i];
         nonNegVarNum++;

         if(dir[nfPos] < MINUSZERO)
         {
            tmp_theta = p_sol[basisIdx[nfPos]] / dir[nfPos];
            // cout << "basisIdx: " << basisIdx[nfPos] << "\n";
            // cout << "tmp_theta: " << tmp_theta << "\n\n";
         }
         else
         {
            checker++;
            tmp_theta = SMALLDOUBLE;
         }
         // if(PLUSZERO < tmp_theta - theta){
         if(theta < tmp_theta)
         {
            theta = tmp_theta;
            pivOutIdx = basisIdx[nfPos];
            sub_pivOutIdx = nfPos;
            flag = CONTINUE;
            // cout << "theta: " << theta << "\n";
            // cout << "pivOutIdx: " << pivOutIdx << "\n\n";
         }
      }
      if(checker == nonNegVarNum) flag = UNBOUNDED;
      break;

   case NEGTHETA:

      theta = (tmp_theta = BIGDOUBLE);

      for(i = 0; i < nfN; i++)
      {
         nfPos = nf_pos[i];
         nonNegVarNum++;

         if(dir[nfPos] > PLUSZERO)
         {
            tmp_theta = p_sol[basisIdx[nfPos]] / dir[nfPos];
	    // cout << "basisIdx: " << basisIdx[nfPos] << "\n";
            // cout << "tmp_theta: " << tmp_theta << "\n\n";
         }
         else
         {
            checker++;
            tmp_theta = BIGDOUBLE;
         }
         // if(theta - tmp_theta > PLUSZERO){
         if(theta > tmp_theta)
         {
            theta = tmp_theta;
            pivOutIdx = basisIdx[nfPos];
            sub_pivOutIdx = nfPos;
            flag = CONTINUE;

            // cout << "theta: " << theta << "\n";
            // cout << "pivOutIdx: " << pivOutIdx << "\n\n";
         }
      }
      if(checker == nonNegVarNum) flag = UNBOUNDED;
      break;
   }
   theta *= -1;
   return (flag);
}

int simplex::ratioTest_artFst
 ( int redFlag, int pivInIdx, int sub_pivInIdx, int& pivOutIdx,
   int& sub_pivOutIdx, double& theta )
{
   int i, flag = 0, nfPos, checker, nonNegVarNum;
   double tmp_theta;
  
   // cout << "<< RatioTest_artFst >>\n\n";

   IP_mat_vec_fst(pivInIdx);

   nonNegVarNum = 0;
   checker = 0;

   switch(redFlag)
   {
      case POSTHETA:
  
         theta = (tmp_theta = SMALLDOUBLE);
      
         for(i = 0; i < nfN; i++)
         {
            nfPos = pre_nf_pos[i];
      
            // cout << "nfPos: " << nfPos << "\n";

            if(pre_basisIdx[nfPos] < termSumNum - supN)
            {
               nonNegVarNum++;

               if(dir[nfPos] < MINUSZERO)
               {
                  tmp_theta = pre_p_sol[pre_basisIdx[nfPos]] / dir[nfPos];
                  // cout << "tmp_theta: " << tmp_theta << "\n";
               }
               else
               {
                  checker++;
                  tmp_theta = SMALLDOUBLE;
               }
               if(theta < tmp_theta)
               {
                  theta = tmp_theta;
                  pivOutIdx = pre_basisIdx[nfPos];
                  sub_pivOutIdx = nfPos;
                  flag = CONTINUE;
                  // cout << "pivOutIdx: " << pivOutIdx << "\n";
               }
            }
            // cout << "\n";
         }
         if(checker == nonNegVarNum)
            flag = UNBOUNDED;
         else
            flag = CONTINUE;

         break;

   case NEGTHETA:

      theta = (tmp_theta = BIGDOUBLE);
    
      for(i = 0; i < nfN; i++)
      {
         nfPos = pre_nf_pos[i];
         // cout << "nfPos: " << nfPos << "\n";
      
         if(pre_basisIdx[nfPos] < termSumNum - supN)
         {
            nonNegVarNum++;

            if(dir[nfPos] > PLUSZERO)
            {
               tmp_theta = pre_p_sol[pre_basisIdx[nfPos]] / dir[nfPos];
               // cout << "tmp_theta: " << tmp_theta << "\n";
            }
            else
            {
               checker++;
               tmp_theta = BIGDOUBLE;
            }
            if(theta > tmp_theta)
            {
               theta = tmp_theta;
               pivOutIdx = pre_basisIdx[nfPos];
               sub_pivOutIdx = nfPos;
               flag = CONTINUE;
 
               // cout << "pivOutIdx: " << pivOutIdx << "\n";
            }
         }
         // cout << "\n";
      }
      if(checker == nonNegVarNum)
         flag = UNBOUNDED;
      else
         flag = CONTINUE;

      break;

   }
   theta *= -1;

   return (flag);
}

int simplex::ratioTest_art
 ( int redFlag, int pivInIdx, int sub_pivInIdx, int& pivOutIdx,
   int& sub_pivOutIdx, double& theta )
{
   int i, flag = 0, nfPos, checker, nonNegVarNum;
   double tmp_theta;

   // cout << "<< Ratio Test Art >>\n\n";

   IP_mat_vec(pivInIdx);

   nonNegVarNum = 0;
   checker = 0;

   switch(redFlag)
   {
      case POSTHETA:
  
      theta = (tmp_theta = SMALLDOUBLE);
      
      for(i = 0; i < nfN; i++)
      {
         nfPos = nf_pos[i];

         if(basisIdx[nfPos] < termSumNum - supN)
         {
            nonNegVarNum++;

            if(dir[nfPos] < MINUSZERO)
            {
               tmp_theta = p_sol[basisIdx[nfPos]] / dir[nfPos];
            }
            else
            {
               checker++;
               tmp_theta = SMALLDOUBLE;
            }
            if(theta < tmp_theta)
            {
               theta = tmp_theta;
               pivOutIdx = basisIdx[nfPos];
               sub_pivOutIdx = nfPos;
               flag = CONTINUE;
	    }
         }
      }
      if(checker == nonNegVarNum)
         flag = UNBOUNDED;
      else
         flag = CONTINUE;

      break;

   case NEGTHETA:

      theta = (tmp_theta = BIGDOUBLE);
    
      for(i = 0; i < nfN; i++)
      {
         nfPos = nf_pos[i];

         if(basisIdx[nfPos] < termSumNum - supN)
         {
            nonNegVarNum++;

            if(dir[nfPos] > PLUSZERO)
            {
               tmp_theta = p_sol[basisIdx[nfPos]] / dir[nfPos];
            }
            else
            {
               checker++;
               tmp_theta = BIGDOUBLE;
	    }
            if(theta > tmp_theta)
            {
               theta = tmp_theta;
               pivOutIdx = basisIdx[nfPos];
               sub_pivOutIdx = nfPos;
               flag = CONTINUE;
            }
         }
      }
      if(checker == nonNegVarNum)
         flag = UNBOUNDED;
      else
         flag = CONTINUE;

      break;
   }
   theta *= -1;

   return (flag);
}

int simplex::ratioTest_art_Bland
 ( int redFlag, int pivInIdx, int sub_pivInIdx, int& pivOutIdx, 
   int& sub_pivOutIdx, double& theta )
{
   int i, flag = 0, nfPos, checker, nonNegVarNum;
   int tmp_basisIdx;

   double tmp_theta;

   // cout << "<< ratioTest_art_Bland >>\n\n";

   IP_mat_vec(pivInIdx);

   nonNegVarNum = 0;
   checker = 0;

   switch(redFlag)
   {
      case POSTHETA:

         theta = 0;
         pivOutIdx = BIGINT;

         for(i = 0; i < nfN; i++)
         {
            nfPos = nf_pos[i];
            tmp_basisIdx = basisIdx[nfPos];

            if(tmp_basisIdx < termSumNum - supN)
            {
               nonNegVarNum++;

               if(dir[nfPos] < MINUSZERO)
               {
                  tmp_theta = p_sol[tmp_basisIdx] / dir[nfPos];

                  if(tmp_basisIdx < pivOutIdx)
                  {
                     theta = tmp_theta;
                     pivOutIdx = tmp_basisIdx;
                     sub_pivOutIdx = nfPos;

                     flag = CONTINUE;
                  }
               }
               else
               {
                  checker++;
               }
            }
         }
         if(checker == nonNegVarNum)
            flag = UNBOUNDED;
         else
            flag = CONTINUE;

         break;

   case NEGTHETA:

      theta = 0;
      pivOutIdx = BIGINT;

      for(i = 0; i < nfN; i++)
      {
         nfPos = nf_pos[i];
         tmp_basisIdx = basisIdx[nfPos];

         if(tmp_basisIdx < termSumNum - supN)
         {
            nonNegVarNum++;

            if(dir[nfPos] > PLUSZERO)
            {
               tmp_theta = p_sol[tmp_basisIdx] / dir[nfPos];

               if(tmp_basisIdx < pivOutIdx)
               {
                  theta = tmp_theta;
                  pivOutIdx = tmp_basisIdx;
                  sub_pivOutIdx = nfPos;

                  flag = CONTINUE;
               }
	    }
            else
            {
               checker++;
	    }
         }
      }
      if(checker == nonNegVarNum)
         flag = UNBOUNDED;
      else
         flag = CONTINUE;

      break;
   }
   theta *= -1;

   return (flag);
}

int simplex::ratioTest_frIdx ( int pivInIdx )
{
   int i, flag = 0, nfPos, checker, nonNegVarNum;
   double theta, tmp_theta;

   // cout << "<< Ratio Test frIdx >>\n\n";
  
   IP_mat_vec(pivInIdx);

   nonNegVarNum = 0;
   checker = 0;

   theta = (tmp_theta = BIGDOUBLE);
    
   for(i = 0; i < nfN; i++)
   {
      nfPos = nf_pos[i];
      
      if(basisIdx[nfPos] < termSumNum - supN)
      {
         nonNegVarNum++;

         if(dir[nfPos] > PLUSZERO)
         {
            tmp_theta = p_sol[basisIdx[nfPos]] / dir[nfPos];
         }
         else
         {
            checker++;
            tmp_theta = BIGDOUBLE;
         }
         if(theta > tmp_theta)
         {
            theta = tmp_theta;
         }
      }
   }    
   if(checker == nonNegVarNum)
      flag = UNBOUNDED;
   else
      flag = OPT;

   return (flag);
}

void simplex::IP_mat_vec ( int pivInIdx )
{
   int i, j, ii, iii, level, idx2;
   int nfPos;
   double val;
  
   // cout << "<< IP_mat_vec >> \n";

   ii = 2 * pivInIdx;
   level = nIdx[ii];

   idx2 = firIdx[level];
   iii = Dim * nIdx[ii + 1];
  
  /*
    cout << "ii: " << ii << "\n";
    cout << "level: " << level << "\n";
    cout << "idx2: " << idx2 << "\n";
    cout << "iii: " << iii << "\n\n";
   */

   for(j = 0; j < nfN; j++)
   {
       nfPos = nf_pos[j];

       val = 0;
       ii = Dim * nfPos;

       for(i = 0; i < Dim; i++)
       {
          val += invB[ii + i] * Supp[level][idx2].supMat[i + iii];
       }
       dir[nfPos] = -1 * val;
   }
}

void simplex::IP_mat_vec_fst( int pivInIdx )
{
   int i, j, ii, iii, level, idx2;
   int nfPos;
   double val;

   // cout << "<< IP_mat_vec_fst >> \n";

   ii = 2 * pivInIdx;
   level = nIdx[ii];

   idx2 = firIdx[level];
   iii = Dim * nIdx[ii + 1];

  /*
    cout << "ii: " << ii << "\n";
    cout << "iii: " << iii << "\n";
    cout << "idx2: " << idx2 << "\n";
    cout << "level: " << level << "\n\n";
   */

   for(j = 0; j < nfN; j++)
   {
      nfPos = pre_nf_pos[j];

      val = 0;
      ii = Dim * nfPos;

      for(i = 0; i < Dim; i++)
      {
         val += pre_invB[ii + i] * Supp[level][idx2].supMat[i + iii];
      }
      dir[nfPos] = -1 * val;
   }
}

void simplex::createNewBandN_tab
 ( int pivInIdx, int sub_pivInIdx, int pivOutIdx, int sub_pivOutIdx,
   double theta, double redCost )
{
   int i, j, ii;
   int nfPos;

   double elem, val, vval;

   // cout << "<< createNewBandN >>\n\n";

   ii = sub_pivOutIdx * Dim;
   elem = dir[sub_pivOutIdx];
   val = (-1 - elem) / elem;
   vval = redCost / elem;

   // make new d_sol
  
   for(i = 0; i < Dim; i++)
   {
      d_sol[i] -= vval * invB[ii + i];
   }

   // make new p_sol

   for(i = 0; i < nfN; i++)
   {
      nfPos = nf_pos[i];
      
      if(nfPos != sub_pivOutIdx)
         p_sol[basisIdx[nfPos]] +=  theta * dir[nfPos];
      else
         p_sol[basisIdx[nfPos]] = 0;
   }
   p_sol[pivInIdx] = theta;

   // cout << "theta: " << theta << "\n\n";

   // make new basisIdx and non_basisIdx

   basisIdx[sub_pivOutIdx] = pivInIdx;
   nbIdx[sub_pivInIdx] = pivOutIdx;

   // make the new basis matrix   

   for(i = 0; i < Dim; i++)
   {
      tmp_newInvB[i] = (invB[ii + i] /= -1 * elem);
   }
   for(j = 0; j < nfN; j++)
   {
      nfPos = nf_pos[j];

      if(nfPos != sub_pivOutIdx)
      {
         ii = Dim * nfPos;

         for(i = 0; i < Dim; i++)
         {
            invB[ii + i] += dir[nfPos] * tmp_newInvB[i];
         }
      }
   }
}

void simplex::createNewBandN_p1
 ( int pivInIdx, int sub_pivInIdx, int pivOutIdx, int sub_pivOutIdx,
   double theta, double redCost, int termS, int reTermS )
{
   int i, j, ii;
   int nfPos;

   double elem, val, vval;

   // cout << "<< createNewBandN >>\n\n";

   ii = sub_pivOutIdx * Dim;
   elem = dir[sub_pivOutIdx];
   val = (-1 - elem) / elem;
   vval = redCost / elem;

   // make new d_sol
  
   for(i = 0; i < Dim; i++)
   {
      d_sol[i] -= vval * invB[ii + i];
   }

   // make new p_sol

   for(i = 0; i < nfN; i++)
   {
      nfPos = nf_pos[i];
      
      if(nfPos != sub_pivOutIdx)
         p_sol[basisIdx[nfPos]] +=  theta * dir[nfPos];
      else
         p_sol[basisIdx[nfPos]] = 0;
   }
   p_sol[pivInIdx] = theta;

   // cout << "theta: " << theta << "\n\n";

   // make new basisIdx and non_basisIdx

   basisIdx[sub_pivOutIdx] = pivInIdx;
   nbIdx[sub_pivInIdx] = pivOutIdx;

   if(pivOutCheck[sub_pivOutIdx] == 0)
   {
      pivOutCheck[sub_pivOutIdx] = 1;

      pivOutList[pivOutNum] = sub_pivOutIdx;
      pivOutNum++;
   }
   if(reTermS <= pivInIdx && pivInIdx < reTermS + termS - 1)
   {
      rIdx[pivInIdx - reTermS] = sub_pivOutIdx;
   }
   if(reTermS <= pivOutIdx  && pivOutIdx < reTermS + termS - 1)
   {
      rIdx[pivOutIdx - reTermS] = -1 * (1 + sub_pivInIdx);
   }

   // make the new basis matrix   

   for(i = 0; i < Dim; i++) 
   {
      tmp_newInvB[i] = (invB[ii + i] /= -1 * elem);

      tmp_transMat[i] = transMat[ii + i];
      transMat[ii + i] +=  val * tmp_transMat[i];
   }
   for(j = 0; j < nfN; j++)
   {
      nfPos = nf_pos[j];

      if(nfPos != sub_pivOutIdx)
      {
         ii = Dim * nfPos;

         for(i = 0; i < Dim; i++)
         {
            invB[ii + i] += dir[nfPos] * tmp_newInvB[i];
            transMat[ii + i] -=  tmp_transMat[i] * (dir[nfPos] / elem);
         }
      }
   }
}

void simplex::createNewBandN
 ( int pivInIdx, int sub_pivInIdx, int pivOutIdx, int sub_pivOutIdx,
   double theta, double redCost, int termS, int reTermS )
{
   int i, j, ii;
   int nfPos;

   double elem, val, vval;

   // cout << "<< createNewBandN >>\n\n";

   ii = sub_pivOutIdx * Dim;
   elem = dir[sub_pivOutIdx];
   val = (-1 - elem) / elem;
   vval = redCost / elem;

   // make new d_sol
  
   for(i = 0; i < Dim; i++)
   {
      d_sol[i] -= vval * invB[ii + i];
      transRed[i] -= vval * transMat[ii + i];
   }

   // make new p_sol

   for(i = 0; i < nfN; i++)
   {
      nfPos = nf_pos[i];
      
      if(nfPos != sub_pivOutIdx)
         p_sol[basisIdx[nfPos]] +=  theta * dir[nfPos];
      else
         p_sol[basisIdx[nfPos]] = 0;
   }
   p_sol[pivInIdx] = theta;

   // cout << "theta: " << theta << "\n\n";

   // make new basisIdx and non_basisIdx

   basisIdx[sub_pivOutIdx] = pivInIdx;
   nbIdx[sub_pivInIdx] = pivOutIdx;

   if(pivOutCheck[sub_pivOutIdx] == 0)
   {
      pivOutCheck[sub_pivOutIdx] = 1;

      pivOutList[pivOutNum] = sub_pivOutIdx;
      pivOutNum++;
   }
   if(reTermS <= pivInIdx && pivInIdx < reTermS + termS - 1)
   {
      rIdx[pivInIdx - reTermS] = sub_pivOutIdx;
   }
   if(reTermS <= pivOutIdx  && pivOutIdx < reTermS + termS - 1)
   {
      rIdx[pivOutIdx - reTermS] = -1 * (1 + sub_pivInIdx);
   }

   // make the new basis matrix   

   for(i = 0; i < Dim; i++)
   {
      tmp_newInvB[i] = (invB[ii + i] /= -1 * elem);

      tmp_transMat[i] = transMat[ii + i];
      transMat[ii + i] +=  val * tmp_transMat[i];
   }
   for(j = 0; j < nfN; j++)
   {
      nfPos = nf_pos[j];

      if(nfPos != sub_pivOutIdx)
      {
         ii = Dim * nfPos;

         for(i = 0; i < Dim; i++)
         {
            invB[ii + i] += dir[nfPos] * tmp_newInvB[i];
            transMat[ii + i] -=  tmp_transMat[i] * (dir[nfPos] / elem);
         }
      }
   }
}

void simplex::createNewBandN_iFst
 ( int pivInIdx, int sub_pivInIdx, int pivOutIdx, int sub_pivOutIdx,
   double theta, double redCost, int termS, int reTermS )
{
   int i, j, ii;
   int cnt, nfPos;

   double elem, val, vval;

   // cout << "<< createNewBandN_iFst >>\n\n";

   ii = sub_pivOutIdx * Dim;
   elem = dir[sub_pivOutIdx];
   val = (-1 - elem) / elem;
   vval = redCost / elem;

   // make new basisIdx and non_basisIdx

   cnt = 0;

   for(i = 0; i < Dim; i++)
   {
      d_sol[i] = pre_d_sol[i] - vval * pre_invB[ii + i];
      transRed[i] -= vval * eye[ii + i];

      nf_pos[i] = pre_nf_pos[i];

      if(i != sub_pivOutIdx)
      {
         basisIdx[i] = pre_basisIdx[i];
      }
      else
      {
         basisIdx[i] = pivInIdx;
      }
   }
   if(pivOutCheck[sub_pivOutIdx] == 0)
   {
      pivOutCheck[sub_pivOutIdx] = 1;

      pivOutList[pivOutNum] = sub_pivOutIdx;
      pivOutNum++;
   }
   if(reTermS <= pivInIdx && pivInIdx < reTermS + termS - 1)
   {
      rIdx[pivInIdx - reTermS] = sub_pivOutIdx;
   }
   if(reTermS <= pivOutIdx  && pivOutIdx < reTermS + termS - 1)
   {
      rIdx[pivOutIdx - reTermS] = -1 * (1 + sub_pivInIdx);
   }

   // make new p_sol

   for(i = 0; i < nfN; i++)
   {
      nfPos = nf_pos[i];
      
      if(nfPos != sub_pivOutIdx)
      {
         p_sol[basisIdx[nfPos]]
            = pre_p_sol[basisIdx[nfPos]] + theta * dir[nfPos];
      }
      else
      {
         p_sol[basisIdx[nfPos]] = 0;
      }
   }
   p_sol[pivInIdx] = theta;
  
   // make the new basis matrix   

   for(i = 0; i < Dim; i++)
   {
      tmp_newInvB[i] = (invB[ii + i] = pre_invB[ii + i] / (-1 * elem));

      tmp_transMat[i] = eye[ii + i];
      transMat[ii + i] = eye[ii + i]  + val * tmp_transMat[i];
   }
   for(j = 0; j < nfN; j++)
   {
      nfPos = nf_pos[j];

      if(nfPos != sub_pivOutIdx)
      {
         ii = Dim * nfPos;

         for(i = 0; i < Dim; i++)
         {
            invB[ii + i] = pre_invB[ii + i ] + dir[nfPos] * tmp_newInvB[i];
            transMat[ii + i]
               = eye[ii + i] - tmp_transMat[i] * (dir[nfPos] / elem);
         }
      }
   }
   // info_invB();
}

void simplex::createNewBandN_mFst
 ( int pivInIdx, int sub_pivInIdx, int pivOutIdx, int sub_pivOutIdx,
   double theta, double redCost, int termS, int reTermS )
{
   int i, j, ii;
   int cnt, nfPos;

   double elem, vval;

   // cout << "<< createNewBandN_mFst >>\n\n";

   ii = sub_pivOutIdx * Dim;
   elem = dir[sub_pivOutIdx];
   vval = redCost / elem;

   // make new basisIdx and non_basisIdx

   cnt = 0;

   for(i = 0; i < Dim; i++)
   {
      // d_sol[i] = pre_d_sol[i] - (redCost * pre_invB[ii + i]) / elem;
      d_sol[i] = pre_d_sol[i] - vval * pre_invB[ii + i];
      transRed[i] = pre_transRed[i] - vval * pre_transMat[ii + i];

      if(pre_nf_pos[i] != sub_pivOutIdx)
      {
         nf_pos[cnt] = pre_nf_pos[i];
         cnt++;
      }
      if(i != sub_pivOutIdx)
      {
         basisIdx[i] = pre_basisIdx[i];
      }
      else
      {
         basisIdx[i] = pivInIdx;
      }
   }
   nfN--;

   if(pivOutCheck[sub_pivOutIdx] == 0)
   {
      pivOutCheck[sub_pivOutIdx] = 1;

      pivOutList[pivOutNum] = sub_pivOutIdx;
      pivOutNum++;
   }
   if(reTermS <= pivInIdx && pivInIdx < reTermS + termS - 1)
   {
      rIdx[pivInIdx - reTermS] = sub_pivOutIdx;
   }
   if(reTermS <= pivOutIdx  && pivOutIdx < reTermS + termS - 1)
   {
      rIdx[pivOutIdx - reTermS] = -1 * (1 + sub_pivInIdx);
   }
   // info_rIdx();

   // make new p_sol

   for(i = 0; i < nfN; i++)
   {
      nfPos = nf_pos[i];
      
      if(nfPos != sub_pivOutIdx)
      {
         p_sol[basisIdx[nfPos]]
            = pre_p_sol[basisIdx[nfPos]] + theta * dir[nfPos];
      }
      else
      {
         p_sol[basisIdx[nfPos]] = 0;
      }
   }
   p_sol[pivInIdx] = theta;

   // info_p_sol();
  
   // make the new basis matrix   

   for(i = 0; i < Dim; i++)
   {
      tmp_newInvB[i] = (invB[ii + i] = pre_invB[ii + i] / (-1 * elem));

      tmp_transMat[i] = pre_transMat[ii + i];
      transMat[ii + i] 
        = pre_transMat[ii + i] + (-1 - elem) / elem * tmp_transMat[i];
   }
   // transMat[ii + sub_pivOutIdx]
   // = pre_transMat[ii + sub_pivOutIdx] + (-1 - elem) / elem;

   for(j = 0; j < nfN; j++)
   {
      nfPos = nf_pos[j];

      if(nfPos != sub_pivOutIdx)
      {
         ii = Dim * nfPos;

         for(i = 0; i < Dim; i++)
         {
            invB[ii + i] = pre_invB[ii + i ] + dir[nfPos] * tmp_newInvB[i];
            transMat[ii + i] 
               = pre_transMat[ii + i] - tmp_transMat[i] * (dir[nfPos] / elem);
         }
         // transMat[ii + sub_pivOutIdx]
         // = pre_transMat[ii + sub_pivOutIdx] - (dir[nfPos] / elem);
      }
   }
   // info_invB();
}

void simplex::createNewBandN_art
 ( int pivInIdx, int sub_pivInIdx, int pivOutIdx, int sub_pivOutIdx,
   double redCost, int termS, int reTermS )
{
   int i, j, ii;
   int nfPos;

   double elem, val, vval;

   ii = sub_pivOutIdx * Dim;
   elem = dir[sub_pivOutIdx];
   val = (-1 - elem) / elem;
   vval = redCost / elem;

   // make new d_sol
  
   for(i = 0; i < Dim; i++)
   {
      d_sol[i] -= vval * invB[ii + i];
      transRed[i] -= vval * transMat[ii + i];

   }

   // make new p_sol
  
   modify_p_sol(pivInIdx);

   // make new basisIdx and non_basisIdx
  
   basisIdx[sub_pivOutIdx] = pivInIdx;

   if(pivOutCheck[sub_pivOutIdx] == 0)
   {
      pivOutCheck[sub_pivOutIdx] = 1;

      pivOutList[pivOutNum] = sub_pivOutIdx;
      pivOutNum++;
   }
   if(reTermS <= pivInIdx && pivInIdx < reTermS + termS - 1)
   {
      rIdx[pivInIdx - reTermS] = sub_pivOutIdx;
   }
   for(i = sub_pivInIdx; i < nbN - Dim - 1; i++)
   {
      nbIdx[i] = nbIdx[i + 1];

      if(reTermS <= nbIdx[i]  && nbIdx[i] < reTermS + termS - 1)
      {
         rIdx[nbIdx[i] - reTermS] = -1 * (i + 1);
      }
   }
   nbN--;

   // make the new basis matrix   

   for(i = 0; i < Dim; i++)
   {
      tmp_newInvB[i] = (invB[ii + i] /= -1 * elem);

      tmp_transMat[i] = transMat[ii + i];
      transMat[ii + i] +=  val * tmp_transMat[i];
   }
   for(j = 0; j < nfN; j++)
   {
      nfPos = nf_pos[j];

      if(nfPos != sub_pivOutIdx)
      {
         ii = Dim * nfPos;

         for(i = 0; i < Dim; i++)
         {
            invB[ii + i] += dir[nfPos] * tmp_newInvB[i];
            transMat[ii + i] -=  tmp_transMat[i] * (dir[nfPos] / elem);
         }
      }
   }
}

void simplex::elimFrIdx ( int sub_pivOutIdx )
{
   int cnt = 0;

   for(int i = 0; i < nfN; i++)
   {
      if(i != sub_pivOutIdx)
      {
         nf_pos[cnt] = nf_pos[i];
         cnt++;
      }
   }
   nfN--;
}

void simplex::update_p1_d_sol ( int pivInIdx, int sub_pivOutIdx )
{
   int i, ii, idx, idx2, level;
   double redCost, val;

   // cout << "<< update_p1_d_sol >>\n\n";

   getIdx(level, idx, idx2, ii, 2 * pivInIdx);
  
   // cout << "level: " << level << "\n";
   // cout << "idx2: " << idx2 << "\n\n";

   redCost = Supp[level][idx2].redVal(p1_d_sol, idx, ii);
  
   ii = sub_pivOutIdx * Dim;
   val = redCost / dir[sub_pivOutIdx];

   // make new d_sol

   for(i = 0; i < Dim; i++)
   {
      p1_d_sol[i] -= val * invB[ii + i];
      transRed[i] -= val * transMat[ii + i];
   }
}

void simplex::modify_p_sol ( int pivInIdx )
{
   int i;
   int level, idx, idx2;

   level = nIdx[2 * pivInIdx];
   idx = nIdx[2 * pivInIdx + 1];
   idx2 = firIdx[level];

#if DBG_MODIFY
   cout << "------- << modify_p_sol >> -------\n\n";

   info_p_sol();
   info_basisIdx();

   cout << "pivInIdx: " << pivInIdx << "\n";
   cout << "supVec: ";

   for(i = 0; i < Dim; i++)
   {
      cout << Supp[level][idx2].supMat_out(i, idx) << " ";
   }
   cout << "\n\n";
#endif

   for(i = 0; i < Dim; i++)
   {
      if(isZero(p_sol[basisIdx[i]]) == TRUE &&
         isZero(Supp[level][idx2].supMat_out(i, idx)) == FALSE)
      {
         calElem(i);
      }
   }

#if DBG_MODIFY
   cout << "\n\n";
#endif

}

void simplex::calElem ( int idx )
{
   int i;

   double val, delta, randNum;
   double upper = BIGDOUBLE, lower = SMALLDOUBLE;
   double tmp_upper, tmp_lower;

#if DBG_MODIFY
   cout << "<< calElem >>\n\n";

   cout << "idx: " << idx << "\n\n";
   info_invB();

   cout.width(15);
   cout << "p_sol: ";
   for(i = 0; i < Dim; i++)
   {
      cout << p_sol[basisIdx[i]] << " ";
   }
   cout << "\n";

   cout.width(15);
   cout << "colum-invB: ";
   for(i = 0; i < Dim; i++)
   {
      cout << invB_out(i, idx) << " ";
   }
   cout << "\n\n";
#endif

   for(i = 0; i < Dim; i++)
   {
      val = invB_out(i, idx);

      if(val > PLUSZERO)
      {
         tmp_lower = -1 * p_sol[basisIdx[i]] / val;

         if(tmp_lower > lower)
         {
            lower = tmp_lower;
         }
      }
      if(val < MINUSZERO)
      {
         tmp_upper = -1 * p_sol[basisIdx[i]] / val;

         if(tmp_upper < upper)
         {
            upper = tmp_upper;
         }
      }
   }

#if DBG_MODIFY
   cout << "lower: " << lower << "\n";
   cout << "upper: " << upper << "\n\n";
#endif

   if(lower == 0 && upper == 0)
      delta = 0;
   else
   {
      if(lower == SMALLDOUBLE && upper == BIGDOUBLE)
      {
         srand(2); // srandom(2);
         delta = (double) rand() / (double) RAND_MAX;
      }
      else if(lower == SMALLDOUBLE && upper < BIGDOUBLE)
      {
         srand(3); // srandom(3);
         randNum = (double) rand() / (double) RAND_MAX;

         delta = upper - randNum / 2;
      }
      else if(lower > SMALLDOUBLE && upper == BIGDOUBLE)
      {
         srand(4); // srandom(4);
         randNum = (double) rand() / (double) RAND_MAX;

         delta = lower + randNum / 2;
      }
      else
      {
         delta = (lower + upper) / 2;
      }
   }

#if DBG_MODIFY
   cout << "delta: " << delta << "\n\n";
#endif

   for(i = 0; i < Dim; i++)
   {
      p_sol[basisIdx[i]] += delta * invB_out(i, idx);
   }

#if DBG_MODIFY
   cout << "-- modified --\n";
   info_p_sol();
#endif

}

void simplex::info_p_sol()
{
   int i;
  
   cout << "<< p_sol >> \n";
   for(i = 0; i < termSumNum - supN + Dim; i++) cout << p_sol[i] << " ";
   cout << "\n\n";
}

void simplex::info_d_sol()
{
   int i;
  
   cout << "<< d_sol >> \n";
   for(i = 0; i < Dim; i++) cout << d_sol[i] << " ";
   cout << "\n\n";
}

void simplex::info_p1_d_sol()
{
   int i;

   cout << "<< p1_d_sol >> \n";
   for(i = 0; i < Dim; i++) cout << p1_d_sol[i] << " ";
   cout << "\n\n";
}

void simplex::info_invB()
{
   int i, j;

   cout << "<< invB >> \n";

   for(i = 0; i < Dim; i++)
   {
      for(j = 0; j < Dim; j++)
      {
         cout.width(10);
         cout << invB_out(i, j) << " ";
      }
      cout << "\n";
   }
   cout << "\n\n";
}

void simplex::info_transMat()
{
   int i, j;
   double val;

   cout << "<< transMat >> \n";

   for(i = 0; i < Dim; i++)
   {
      for(j = 0; j < Dim; j++)
      {
         cout.width(10);

         val = transMat_out(i, j);

         if(val < PLUSZERO && val > MINUSZERO)
         {
            cout << "0 ";
         }
         else
         {
	    cout << transMat_out(i, j) << " ";
         }
      }
      cout << "\n";
   }
   cout << "\n\n";
}

void simplex::info_transRed()
{
   int i;

   double val;

   cout << "<< transRed >> \n";

   for(i = 0; i < Dim; i++)
   {
      cout.width(10);

      val = transRed[i];

      if(val < PLUSZERO && val > MINUSZERO)
      {
         cout << "0 ";
      }
      else
      {
         cout << val << " ";
      }
   }
   cout << "\n\n";
}

void simplex::info_basisIdx()
{
   int i;

   cout << "<< basisIdx >> \n";
   for(i = 0; i < Dim; i++) cout << basisIdx[i] << " ";
   cout<< "\n\n";
}

void simplex::info_nf_pos()
{
   int i;

   cout << "<< nf_pos >> \n";
   for(i = 0; i < nfN; i++) cout << nf_pos[i] << " ";
   cout << "\n\n";
}

void simplex::info_nbIdx()
{
   int i;

   cout << "<< nbIdx >> \n";
   for(i = 0; i < nbN - Dim; i++) cout << nbIdx[i] << " ";
   cout << "\n\n";
}

void simplex::info_rIdx()
{
   int i;

   cout << "<< rIdx >>\n";
   for(i = 0; i < nbN; i++) cout << rIdx[i] << " ";
   cout << "\n\n";
}

void simplex::info_redVec()
{
   int i;
 
   cout << "<< redVec >>\n";
   for(i = 0; i < nbN - Dim; i++) cout << redVec[i] << " ";
   cout << "\n\n";
}

void simplex::info_dir()
{
   int i;

   cout << "<< dir >> \n";

  /*
    int nfPos;

    for(i = 0; i < nfN; i++)
    {
       nfPos = nf_pos[i];

       // cout << basisIdx[nfPos] << " : " << dir[nf_pos[i]] << "\n";
       cout << nfPos << " : " << dir[nf_pos[i]] << "\n";
    }
    cout << "\n";
   */

   for(i = 0; i < Dim; i++)
   {
      cout << dir[i] << " ";
   }
   cout << "\n\n";
}

void simplex::info_frIdx()
{
   cout << "frIdx: " << frIdx << "\n\n";
}

void simplex::info_candIdx()
{
   int i;

   cout << "<< candIdx >>\n\n";

   for(i = 0; i < candIdx[0]; i++)
   {
      cout << candIdx[i + 1] << " ";
   }
   cout << "\n\n";
}

void simplex::info_repIdx()
{
   cout << "repIdx: " << repIdx << "\n\n";
}

void simplex::info_oriSup()
{
   int i, j, k;

   cout << "<< oriSup >>\n";

   for(k = 0; k < supN; k++)
   {
      for(j = 0; j < Dim; j++)
      {
         for(i = 0; i < termSet[k]; i++)
         {
            cout << supp_out(k, j, i) << " ";
         }
         cout << "\n";
      }
      cout << "\n";
   }
   cout << "\n\n";
}

int simplex::tSolLP(int& iter, int mode)
{
   int iter_1 = 0, iter_2 = 0;
   int pivInIdx, sub_pivInIdx, pivOutIdx, sub_pivOutIdx, flag, redFlag;
   double theta, redCost;

   /////////////// phase_1 ////////////////////////////////////////////

#if DBG_TSOLLP
   info_basisIdx();
   info_nbIdx();
   info_invB();
   info_p_sol();
   info_d_sol();
   info_frIdx();
#endif

   while(1)
   {
    
#if DBG_TSOLLP
      cout << "----- Phase-1. Iter: " << iter_1 << " -----\n\n";
#endif

      redFlag = reducedCost_tab_p1(pivInIdx, sub_pivInIdx, redCost);

      if(redFlag == OPT)
      {
         reMakeNonBasisIdx_tab();

#if DBG_TSOLLP
         cout << "----- OPT.Phase-1 -----\n";

         info_basisIdx();
         info_nbIdx();
         info_invB();
         info_p_sol();
         info_d_sol();
         info_dir();
         info_nf_pos();
#endif

         break;
      }

      flag = ratioTest(redFlag, pivInIdx, sub_pivInIdx, pivOutIdx,
                       sub_pivOutIdx, theta);

      createNewBandN_tab(pivInIdx, sub_pivInIdx, pivOutIdx, sub_pivOutIdx, 
                         theta, redCost);

#if DBG_TSOLLP  
      info_basisIdx();
      info_nbIdx();
      info_invB();
      info_p_sol();
      info_d_sol();
      info_dir();
      info_nf_pos();
#endif
    
      iter_1++;

      if(iter_1 > ITER)
      {
         flag = ERROR_ITER;

         break;
      }
   }

   /////////////// phase_2 ////////////////////////////////////////////
  
   if(mode == TRIANGLE && flag == CONTINUE)
   {
      flag = checkFrIdx();

#if DBG_TSOLLP	
      cout << "----- checkFrIdx -----\n";

      info_basisIdx();
      info_nbIdx();
      info_invB();
      info_p_sol();
      info_d_sol();
      info_dir();
      info_nf_pos();
#endif

   }
   if(flag == CONTINUE)
   {
      while(1)
      {

#if DBG_TSOLLP      
         cout << "----- Phase-2. Iter: " << iter_2 << " -----\n\n";
#endif
      
         redFlag = reducedCost_tab(pivInIdx, sub_pivInIdx, redCost);

         if(redFlag == OPT)
         {
            flag = OPT;

#if DBG_TSOLLP	
            cout << "----- OPT.Phase-2 -----\n";
	
            info_basisIdx();
            info_nbIdx();
            info_invB();
            info_p_sol();
            info_d_sol();
            info_dir();
            info_nf_pos();
#endif

            break;
         }
         flag = ratioTest_art(redFlag, pivInIdx, sub_pivInIdx, pivOutIdx,
                              sub_pivOutIdx, theta);

         if(flag == UNBOUNDED)
         {

#if DBG_TSOLLP	
            cout << "----- UNB.Phase-2 -----\n";
#endif
	
            break;
         }
         createNewBandN_tab(pivInIdx, sub_pivInIdx, pivOutIdx, sub_pivOutIdx, 
                            theta, redCost);
      
         iter_2++;

#if DBG_TSOLLP
         info_basisIdx();
         info_nbIdx();
         info_invB();
         info_p_sol();
         info_d_sol();
         info_dir();
         info_nf_pos();
#endif

         if(iter_2 > ITER)
         {
            flag = ERROR_ITER;

            break;
         }
      }
   }
   iter = iter_1 + iter_2;

   return (flag);
}

int simplex::checkFrIdx()
{
   int i;
   int flag = CONTINUE, redFlag;

   int pivInIdx, sub_pivInIdx, pivOutIdx, sub_pivOutIdx;

   double redCost, theta;

   for(i = 0; i < nbN - Dim; i++)
   {
      if(nbIdx[i] == frIdx)
      {
         pivInIdx = frIdx;
         sub_pivInIdx = i;

         flag = PIVOT_IN;

         break;
      }
   }
   if(flag == PIVOT_IN)
   {
      calRedCost(pivInIdx, redCost);

      if(redCost > PLUSZERO)
         redFlag = NEGTHETA;
      else
         redFlag = POSTHETA;

      flag = ratioTest_art(redFlag, pivInIdx, sub_pivInIdx,
                           pivOutIdx, sub_pivOutIdx, theta);

      if(flag == CONTINUE)
      {
         elimFrIdx(sub_pivOutIdx);

         createNewBandN_tab(pivInIdx, sub_pivInIdx, pivOutIdx, sub_pivOutIdx, 
                            theta, redCost);
      }
      else
      {

#if DBG_TSOLLP
         cout << "----- UNB.checkFrIdx -----\n";
#endif

      }
   }
   else
   {
      for(i = 0; i < Dim; i++)
      {
         if(basisIdx[i] == frIdx)
         {
            elimFrIdx(i);

            break;
         }
      }
   }
   return (flag);
}

int simplex::fSolLP( int termS, int reTermS, int& iter )
{
   int iter_1 = 0, iter_2 = 0;

   int pivInIdx, sub_pivInIdx, pivOutIdx, sub_pivOutIdx, flag, redFlag;
   double theta, redCost;

   /////////////// phase_1 ////////////////////////////////////////////

#if DBG_SOLLP
   info_basisIdx();
   info_nbIdx();
   info_rIdx();
   info_invB();
   info_p_sol();
   info_d_sol();
   info_dir();
#endif

   while(1)
   {
    
#if DBG_SOLLP
      cout << "----- Phase-1. Iter: " << iter_1 << " -----\n\n";
#endif

      // cout << "----- Phase-1. Iter: " << iter_1 << " -----\n\n";
    
      redFlag = reducedCost_p1(pivInIdx, sub_pivInIdx, redCost);

      if(redFlag == OPT)
      {
         reMakeNonBasisIdx(reTermS);

#if DBG_SOLLP
         cout << "----- OPT.Phase-1 -----\n";

         info_basisIdx();
         info_nbIdx();
         info_rIdx();
         info_invB();
         info_p_sol();
         info_d_sol();
         info_dir();
         info_nf_pos();
#endif
         break;
      }
      flag = ratioTest(redFlag, pivInIdx, sub_pivInIdx, pivOutIdx,
                       sub_pivOutIdx, theta);
  
      update_p1_d_sol(pivInIdx, sub_pivOutIdx);

      createNewBandN_p1(pivInIdx, sub_pivInIdx, pivOutIdx, sub_pivOutIdx, 
                        theta, redCost, termS, reTermS);

#if DBG_SOLLP  
      info_basisIdx();
      info_nbIdx();
      info_rIdx();
      info_invB();
      info_p_sol();
      info_d_sol();
      info_dir();
      info_nf_pos();
#endif
    
      iter_1++;

      if(iter_1 > ITER)
      {
         flag = ERROR_ITER;

         break;
      }
   }

   /////////////// phase_2 ////////////////////////////////////////////
  
   if(flag == CONTINUE)
   {
      while(1)
      {

#if DBG_SOLLP      
         cout << "----- Phase-2. Iter: " << iter_2 << " -----\n\n";
#endif
      
         redFlag = reducedCost(pivInIdx, sub_pivInIdx, redCost);

         if(redFlag == OPT)
         {
            flag = OPT;

#if DBG_SOLLP	
            cout << "----- OPT.Phase-2 -----\n";
	
            info_basisIdx();
            info_nbIdx();
            info_rIdx();
            info_invB();
            info_p_sol();
            info_d_sol();
            info_dir();
            info_nf_pos();
#endif
	
            break;
         }
         flag = ratioTest_art(redFlag, pivInIdx, sub_pivInIdx,
                              pivOutIdx, sub_pivOutIdx, theta);

         if(flag == UNBOUNDED)
         {
            break;
         }
         createNewBandN(pivInIdx, sub_pivInIdx, pivOutIdx, sub_pivOutIdx, 
                        theta, redCost, termS, reTermS);

         iter_2++;

#if DBG_SOLLP
         info_basisIdx();
         info_nbIdx();
         info_rIdx();
         info_invB();
         info_p_sol();
         info_d_sol();
         info_dir();
         info_nf_pos();
#endif

         if(iter > ITER_BLAND)
         {
            flag = solLP_art_Bland(pivInIdx,sub_pivInIdx,pivOutIdx, 
                                   sub_pivOutIdx,redFlag,theta,redCost,
                                   termS,reTermS,iter_2);
	    break;
         }
      }
   }
   iter = iter_1 + iter_2;

   return (flag);
}

void simplex::cal_redVec
 ( int termS, int reTermS, int fst_pivInIdx, theData** cur )
{
   int i, cnt;
  
   cnt = 0;

   // cout << "<< cal_redVec >>\n\n";

   for(i = 0; i < termS; i++)
   {
      if(i != fst_pivInIdx)
      {
         (*cur)->redVec[cnt+reTermS]
            = fst_redVec[i] - fst_redVec[fst_pivInIdx];
         cnt++;
      }
   }
   // cout << "\n\n";
}

void simplex::get_iNbN_nfN ( theData** cur, int lNbN, int lNfN )
{
   (*cur)->nbN = (nbN = lNbN);
   (*cur)->nfN = (nfN = lNfN);
}

void simplex::get_res ( ftData& iData )
{
   // cout << "<< get_res >>\n\n";

   iData.cur->artV = artV;
   artV = 0;

   iData.cur->nbN = nbN;
   // iData.cur->pivOutNum = pivOutNum;
}

void simplex::get_pivOutNum ( theData** cur )
{
   (*cur)->pivOutNum = pivOutNum;
}

void simplex::get_nbN_nfN ( int ori_nbN, int ori_nfN )
{
   nbN = ori_nbN;
   nfN = ori_nfN;
}

void simplex::get_p_sol ( double* ori_p_sol )
{
   p_sol = ori_p_sol;
}

void simplex::get_d_sol ( double* ori_d_sol )
{
   d_sol = ori_d_sol;
}

void simplex::get_basisIdx ( int* ori_basisIdx )
{
  basisIdx = ori_basisIdx;
}

void simplex::get_nf_pos( int* ori_nf_pos )
{
   nf_pos = ori_nf_pos;
}

void simplex::get_nbIdx ( int* ori_nbIdx )
{
   nbIdx = ori_nbIdx;
}

void simplex::get_invB ( double* ori_invB )
{
   invB = ori_invB;
}

void simplex::get_frIdx ( int ori_frIdx )
{
   frIdx = ori_frIdx;
}

void simplex::fstRed_candIdx
 ( inifData& curInif, int** mCandIdx, int& pivInIdx, int& sub_pivInIdx )
{
   int num = 0, idx;
   double tmp_redCost, redCost = 1.0E+8;

   uData *n_curr;
  
   // cout << "<< fstRed_candIdx >> \n\n";

   n_curr = curInif.fHead;

   while(n_curr != NULL)
   {
      tmp_redCost = n_curr->red;
      idx = n_curr->supLab;

      fst_redVec[idx] = tmp_redCost;
      (*mCandIdx)[num + 1] = idx;

      // cout << "idx: " << idx  <<  " redCost: " << tmp_redCost << "\n";

      if(tmp_redCost < redCost)
      {
         redCost = tmp_redCost;

         pivInIdx = idx;
         sub_pivInIdx = num;
      }
      n_curr = n_curr->fNext;
      num++;
   }
   // cout << "\n";
   (*mCandIdx)[0] = num;
}

int simplex::solLP_art
 ( int depth, int idx_one, int fst_pivIn, int preNbN, int termS, 
   int reTermS, int& iter )
{
   int pivInIdx, sub_pivInIdx, pivOutIdx, sub_pivOutIdx;
   int flag, redFlag;

   double redCost, theta;

  /////////////// phase_2 ////////////////////////////////////////////
  
#if DBG_ISOLLP
   cout << "<< iSolLP_Art >>\n\n";

   info_basisIdx();
   info_nbIdx();
   info_invB();
   info_p_sol();
   info_d_sol();
   info_nf_pos();
   info_rIdx();
   info_dir();
#endif
  
   elimArt(depth, preNbN, termS, reTermS, iter);

   while(1)
   {

#if DBG_ISOLLP
      cout << "----- Iter: " << iter << " -----\n\n";

      info_basisIdx();
      info_nbIdx();
      info_invB();
      info_p_sol();
      info_d_sol();
      info_rIdx();
      info_dir();
#endif

      redFlag = reducedCost(pivInIdx, sub_pivInIdx, redCost);

      if(redFlag == OPT)
      {

#if DBG_ISOLLP
         cout << "----- OPT -----\n\n";

         info_basisIdx();
         info_nbIdx();
         info_invB();
         info_p_sol();
         info_d_sol();
         info_rIdx();
         info_dir();
#endif

         flag = OPT;
         break;
      }
      flag = ratioTest_art(redFlag, pivInIdx, sub_pivInIdx, pivOutIdx,
                           sub_pivOutIdx, theta);

      if(flag == UNBOUNDED)
      {
         break;
      }
      createNewBandN(pivInIdx, sub_pivInIdx, pivOutIdx, sub_pivOutIdx, 
                     theta, redCost, termS, reTermS);
      iter++;
 
      if(iter > ITER_BLAND)
      {
         flag = solLP_art_Bland(pivInIdx,sub_pivInIdx,pivOutIdx,sub_pivOutIdx,
                                redFlag,theta,redCost,termS,reTermS,iter);
         break;
      }
   }
   return (flag);
}

int simplex::solLP_art_Bland
 ( int pivInIdx, int sub_pivInIdx, int pivOutIdx, int sub_pivOutIdx,
   int redFlag, double theta, double redCost, int termS, int reTermS,
   int& iter)
{

   int flag;

#if DBG_MSOLLP
   cout << "<< solLP_art_Bland >>\n\n";
#endif

   // cout << "----- Iter: " << iter << " -----\n\n";

   while(1)
   {
    
#if DBG_MSOLLP
      cout << "----- Iter: " << iter << " -----\n\n";
    
      info_basisIdx();
      info_nbIdx();
      info_invB();
      info_p_sol();
      info_d_sol();
      info_rIdx();
      info_dir();
#endif

      // cout << "----- Iter: " << iter << " -----\n\n";

      redFlag = reducedCost_Bland(pivInIdx, sub_pivInIdx, redCost);

      if(redFlag == OPT)
      {
          
#if DBG_MSOLLP
         cout << "----- OPT -----\n\n";

         info_basisIdx();
         info_nbIdx();
         info_invB();
         info_p_sol();
         info_d_sol();
         info_rIdx();
         info_dir();
#endif

         flag = OPT;
         break;
      }
      flag = ratioTest_art_Bland(redFlag, pivInIdx, sub_pivInIdx,
                                 pivOutIdx, sub_pivOutIdx, theta);
      if(flag == UNBOUNDED)
      {
         break;
      }
      createNewBandN(pivInIdx,sub_pivInIdx,pivOutIdx,sub_pivOutIdx,theta,
                     redCost,termS,reTermS);
      iter++;

      if(iter > ITER)
      {
         flag = ERROR_ITER;

         info_redVec();
         info_dir();
         info_basisIdx();
         info_nbIdx();
         info_nf_pos();
         info_invB();
         info_p_sol();
         info_d_sol();

         break;
      }
   }
   return (flag);
}

int simplex::solLP
 ( int depth, int fst_pivInIdx, int fst_sub_pivInIdx, double fst_redCost,
   int mode, int termS, int reTermS, int preNbN, int& iter )
{
   int pivInIdx, sub_pivInIdx, pivOutIdx, sub_pivOutIdx;
   int flag, redFlag;

   double theta, redCost;
  
   /////////////// Phase_2 ////////////////////////////////////////////

#if DBG_MSOLLP
   cout << "<< mSolLP >>\n\n";
#endif

   // cout << "----- Iter: " << iter << " -----\n\n";
 
   flag = initIter(mode, fst_pivInIdx, fst_sub_pivInIdx, fst_redCost,
                   redFlag, pivInIdx, sub_pivInIdx, pivOutIdx, sub_pivOutIdx,
                   theta, redCost, termS, reTermS, preNbN);
   iter++;

   if(flag != CONTINUE)
      return (flag);
   else
   {
      while(1)
      {

#if DBG_MSOLLP
         cout << "----- Iter: " << iter << " -----\n\n";
         info_basisIdx();
         info_nbIdx();
         info_invB();
         info_p_sol();
         info_d_sol();
         info_rIdx();
         info_dir();
#endif
    
         // cout << "----- Iter: " << iter << " -----\n\n";

         flag = ratioTest_art(redFlag, pivInIdx, sub_pivInIdx, pivOutIdx,
                              sub_pivOutIdx, theta);

         if(flag == UNBOUNDED)
         {
            break;
         }
         createNewBandN(pivInIdx, sub_pivInIdx, pivOutIdx, sub_pivOutIdx,
                        theta, redCost, termS, reTermS);

         redFlag = reducedCost(pivInIdx, sub_pivInIdx, redCost);

         if(redFlag == OPT)
         {

#if DBG_MSOLLP
            cout << "----- OPT -----\n\n";
            info_basisIdx();
            info_nbIdx();
            info_invB();
            info_p_sol();
            info_d_sol();
            info_rIdx();
            info_dir();
#endif
	
            flag = OPT;

            break;
         }
         iter++;

         if(iter > ITER_BLAND)
         {
            flag = solLP_Bland(pivInIdx,sub_pivInIdx,pivOutIdx,sub_pivOutIdx,
                               redFlag,theta,redCost,termS,reTermS,iter);
            break;
         }
      }
      return (flag);
   }
}

int simplex::solLP_Bland
 ( int pivInIdx, int sub_pivInIdx, int pivOutIdx, int sub_pivOutIdx,
   int redFlag, double theta, double redCost, int termS, int reTermS,
   int& iter )
{
   int flag;

#if DBG_MSOLLP
   cout << "<< solLP_Bland >>\n\n";
#endif

   // cout << "----- Iter: " << iter << " -----\n\n";

   while(1)
   {
    
#if DBG_MSOLLP
      cout << "----- Iter: " << iter << " -----\n\n";
      info_basisIdx();
      info_nbIdx();
      info_invB();
      info_p_sol();
      info_d_sol();
      info_rIdx();
      info_dir();
#endif

      // cout << "----- Iter: " << iter << " -----\n\n";

      flag = ratioTest_art_Bland(redFlag, pivInIdx, sub_pivInIdx, pivOutIdx,
                                 sub_pivOutIdx, theta);

      if(flag == UNBOUNDED)
      {
         break;
      }
      createNewBandN(pivInIdx,sub_pivInIdx,pivOutIdx,sub_pivOutIdx,theta,
                     redCost,termS,reTermS);

      redFlag = reducedCost_Bland(pivInIdx, sub_pivInIdx, redCost);

      if(redFlag == OPT)
      {
          
#if DBG_MSOLLP
         cout << "----- OPT -----\n\n";
         info_basisIdx();
         info_nbIdx();
         info_invB();
         info_p_sol();
         info_d_sol();
         info_rIdx();
         info_dir();
#endif
	
         flag = OPT;
         break;
      }
      iter++;

      if(iter > ITER)
      {
         flag = ERROR_ITER;

         info_redVec();
         info_dir();
         info_basisIdx();
         info_nbIdx();
         info_nf_pos();
         info_invB();
         info_p_sol();
         info_d_sol();

         break;
      }
   }
   return (flag);
}

int simplex::initIter
 ( int mode, int fst_pivInIdx, int fst_sub_pivInIdx, double fst_redCost,
   int& redFlag, int& pivInIdx, int& sub_pivInIdx, int& pivOutIdx, 
   int& sub_pivOutIdx, double& theta, double& redCost, 
   int termS, int reTermS, int preNbN )
{
   int flag = -1;
  
   // cout << "<< initIter >>\n\n";
  
   switch(mode)
   {
      case ICHECK:

      flag = ratioTest_artFst(POSTHETA, fst_pivInIdx, fst_sub_pivInIdx, 
                              pivOutIdx, sub_pivOutIdx, theta);

      if(flag == UNBOUNDED)
      {
         return (flag);
      }
      else
      {
         redCost = fst_redCost;
    
         createNewBandN_iFst(fst_pivInIdx,fst_sub_pivInIdx,pivOutIdx,
                             sub_pivOutIdx,theta,redCost,termS,reTermS);

         redFlag = reducedCost_iFst(fst_pivInIdx, fst_sub_pivInIdx, 
                                    pivOutIdx, sub_pivOutIdx, redCost, 
                                    termS, reTermS, preNbN);
         if(redFlag == OPT)
         {
          
#if DBG_MSOLLP
            cout << "----- OPT -----\n\n";
            info_basisIdx();
            info_nbIdx();
            info_invB();
            info_p_sol();
            info_d_sol();
            info_dir();
#endif

            return (redFlag);
         }
         else
         {
            pivInIdx = fst_pivInIdx;
            sub_pivInIdx = fst_sub_pivInIdx;

            return (CONTINUE);
         }
      }
      break;
    
   case MCHECK:
    
      flag = ratioTest_artFst(NEGTHETA, fst_pivInIdx, fst_sub_pivInIdx, 
                              pivOutIdx, sub_pivOutIdx, theta);
    
      if(flag == UNBOUNDED)
      {
         return (flag);
      }
      else
      {
         redCost = fst_redCost;
    
         createNewBandN_mFst(fst_pivInIdx,fst_sub_pivInIdx,pivOutIdx,
                             sub_pivOutIdx,theta,redCost,termS,reTermS);

         redFlag = reducedCost_mFst(fst_pivInIdx, fst_sub_pivInIdx, 
                                    pivOutIdx, sub_pivOutIdx, redCost);
         if(redFlag == OPT)
         {
          
#if DBG_MSOLLP
            cout << "----- OPT -----\n\n";
            info_basisIdx();
            info_nbIdx();
            info_invB();
            info_p_sol();
            info_d_sol();
            info_dir();
#endif

            return (redFlag);
         }
         else
         {
            pivInIdx = fst_pivInIdx;
            sub_pivInIdx = fst_sub_pivInIdx;

            return (CONTINUE);
         }
      }
      break;
   }
   return (flag); // added because end of control of non-void function ...
}

void simplex::get_mNbN_nfN ( theData* parent, theData** cur )
{
   (*cur)->nbN = (nbN = parent->nbN);
   (*cur)->nfN = (nfN = parent->nfN);
}

void simplex::get_repIdx_candIdx ( int* ori_candIdx, int ori_repIdx )
{
   repIdx = ori_repIdx;
   candIdx = ori_candIdx;
}

void simplex::get_parent ( theData* parent )
{
   // cout << "<< get_parent >> \n\n";
  	
   pre_p_sol = parent->p_sol_ptr;
   pre_d_sol = parent->d_sol_ptr;

   pre_redVec = parent->redVec_ptr;
  
   pre_basisIdx = parent->basisIdx_ptr;
   pre_nbIdx = parent->nbIdx_ptr;
   pre_nf_pos = parent->nf_pos_ptr;

   pre_invB = parent->invB_ptr;
   pre_transMat = parent->transMat_ptr;
   pre_transRed = parent->transRed_ptr;
}

void simplex::get_cur ( theData** cur )
{
   p_sol = (*cur)->p_sol;
   d_sol = (*cur)->d_sol;

   redVec = (*cur)->redVec;

   basisIdx = (*cur)->basisIdx;
   nf_pos = (*cur)->nf_pos;
   nbIdx = (*cur)->nbIdx;

   rIdx = (*cur)->rIdx;

   invB = (*cur)->invB;
   transMat = (*cur)->transMat;
   transRed = (*cur)->transRed;

   pivOutList = (*cur)->pivOutList;
   pivOutCheck = (*cur)->pivOutCheck;

   pivOutNum = (*cur)->pivOutNum;
}

void simplex::calMixedVol ( lvData* lv, int* sp, int supN )
{
   int i, j, k, ii, jj, idx;
   int cnt;
   int polyDim;
   int fIdx;

#ifdef compile4phc
   ostringstream strcell;
#endif

   double det;

   // cout << "<< calMixedVol >>\n";

   mixedCell++;

#ifdef compile4phc
   strcell << "# " << mixedCell << " : ";
#endif
   if(output)
   {
      cout << "# " << mixedCell << " : ";
   }
   cnt = 0;
   for(i = 0; i < supN; i++)
   {
      polyDim = lv[sp[i]].Node->parent->polyDim; 

      fIdx = lv[sp[i]].Node->parent->nodeLabel[0];
      jj = fIdx * Dim;

#ifdef compile4phc
      strcell << sp[i] + 1 << " : " << "( " << fIdx + 1 << " ";
#endif
      if(output)
      {
         cout << sp[i] + 1 << " : " << "( " << fIdx + 1 << " ";
      }
      for(j = 0; j < polyDim; j++)
      {
         idx = lv[sp[i]].Node->parent->nodeLabel[j + 1];
         ii = idx * Dim;

#ifdef compile4phc
         strcell << idx + 1 << " ";
#endif
         if(output)
         {
            cout << idx + 1 << " ";
         }
         // cout << "first " << cnt << " : ";

         for(k = 0; k < Dim; k++)
         {
            vol[Dim*cnt+k] = oriSupp[sp[i]][k + ii] - oriSupp[sp[i]][k + jj];
         }
         cnt++;
      }
#ifdef compile4phc
      strcell << ") ";
#endif
      if(output)
      {
         cout << ") ";
      }
   }
   det = fabs(lu(Dim, vol));
   mixedVol += det;

#ifdef compile4phc
   strcell << "volume : " << lrint(det);
   strcell << " accumulated volume : " << lrint(mixedVol);
   // cout << "strcell: " << strcell.str() << endl;
   int fail = demics_append_cell_indices(strcell.str());
#endif
   if(output)
   {
      cout << endl;
      cout << "Volume: " << det << endl << endl;
   }
}

double simplex::lu ( int n, double* a )
{
   int i, j = 0, k, ii, ik;
   double t, u, det;

   det = 0;                  

   for(k = 0; k < n; k++)
   {
      ip[k] = k;             
      u = 0;                 

      for(j = 0; j < n; j++)
      {
         t = fabs(a[n * k + j]);  if (t > u) u = t;
      }
      if (u == 0) goto EXIT; 
      weight[k] = 1 / u;     
   }
   det = 1;                 

   for(k = 0; k < n; k++)
   {
      u = -1;

      for(i = k; i < n; i++)
      { 
         ii = ip[i];            
         t = fabs(a[n * ii + k]) * weight[ii];
         if(t > u) {  u = t;  j = i;  }
      }
      ik = ip[j];

      if(j != k)
      {
         ip[j] = ip[k];  ip[k] = ik;
         det = -det; 
      }
      u = a[n * ik + k];  det *= u; 
    
      if (u == 0) goto EXIT;    
      for(i = k + 1; i < n; i++)
      { 
         ii = ip[i];
         t = (a[n * ii + k] /= u);
         for(j = k + 1; j < n; j++)
            a[n * ii + j] -= t * a[n * ik + j];
      }
   }
   EXIT:
      return det;      
}

double simplex::matinv ( int n, double* a, double* a_inv )
{
   int i, j, k, ii;
   double t, det;

   det = lu(n, a);
  
   if(det != 0)
      for(k = 0; k < n; k++)
      {
         for(i = 0; i < n; i++)
         {
            ii = ip[i];  t = (ii == k);
            for(j = 0; j < i; j++)
               t -= a[n * ii + j] * a_inv[n * j + k];
            a_inv[n * i + k] = t;
         }
         for(i = n - 1; i >= 0; i--)
         {
            t = a_inv[n * i + k];  ii = ip[i];

            for(j = i + 1; j < n; j++)
               t -= a[n * ii + j] * a_inv[n * j + k];
            a_inv[n * i + k] = t / a[n * ii + i];
         }
      }
   return det;
}

void simplex::check_dirRed ( theData* parent, int depth )
{
   int i, j, k, ii;
   int nfPos, nfN, cnt = 0;
   int *nf_pos;
  
   double val, red; 
   double *invB, *d_sol;
    
   cout << "----- << check_dirRed >> -----\n\n";

   // parent->info_invB();
  
   invB = parent->invB;
   d_sol = parent->d_sol;
   nf_pos = parent->nf_pos;
   nfN = parent->nfN;

   // info_oriSup();

   cout << "[ Direction and Recued Cost] \n\n";

   for(ii = depth + 1; ii < supN; ii++)
   {
      cnt = 0;

      cout << "--- Support: " << ii + 1 << " ---\n";

      for(k = 0; k < termSet[ii]; k++)
      {
         cout.width(3);
         cout << cnt + 1 << " : ";

         // for(j = 0; j < Dim; j++){
         for(j = 0; j < nfN; j++)
         {
            nfPos = nf_pos[j];
            val = 0;

            for(i = 0; i < Dim; i++)
            {
               val += invB[i + nfPos * Dim] * supp_out(ii, i, k);
               // cout <<  supp_out(0, i, j) << " ";
            }
            if(val < PLUSZERO && val > MINUSZERO)
            {
               cout.setf(ios::right);
               cout.width(4);
               cout << "0 ";
            }
            else
            {
               cout.setf(ios::right);
               cout.width(4);
               cout << val << " ";
            }
         }
         val = 0;

         for(i = 0; i < Dim; i++)
         {
            val += d_sol[i] * supp_out(ii, i, k);
         }
         red = lifting[termStart[ii] + k] - val;

         cout.setf(ios::right);
         cout.width(4);
         cout << " : " << red << "\n";

         cnt++;
      }
      cout << "\n";
   }
   cout << "\n\n";
}

void simplex::dbg_dirRed ( theData* parent, inifData* nextInif, int depth )
{
   int i, j, k, ii;
   int nfPos, nfN, cnt = 0;
   int *nf_pos;
  
   double val, red, diff; 
   double *invB, *d_sol;

   uData *n_curr;

   // cout << "----- << dbg_dirRed >> -----\n\n";
   // parent->info_invB();
  
   invB = parent->invB;
   d_sol = parent->d_sol;
   nf_pos = parent->nf_pos;
   nfN = parent->nfN;

   // info_oriSup();
   // cout << "[ Direction ] \n\n";

   for(ii = depth + 1; ii < supN; ii++)
   {
      cnt = 0;

      n_curr = nextInif[ii].fHead;
    
      // cout << "--- Support: " << ii + 1 << " ---\n";
    
      for(k = 0; k < termSet[ii]; k++)
      {
         cout.width(3);
         // cout << cnt + 1 << " : ";
         // For direction 
         for(j = 0; j < nfN; j++)
         {
            nfPos = nf_pos[j];
            val = 0;

            for(i = 0; i < Dim; i++)
            {
               val += invB[i + nfPos * Dim] * supp_out(ii, i, k);
            }
            // cout << val << "=" << n_curr->dir[nfPos] << " ";

            diff = val - n_curr->dir[nfPos];

            if(diff > PLUSZERO || diff < MINUSZERO)
            {
               cout << "dbg_dirRed:  ERROR -- Direction!! \n\n";
            }
         }
         // For reduced cost

         val = 0;

         for(i = 0; i < Dim; i++)
         {
            val += d_sol[i] * supp_out(ii, i, k);
         }
         red = lifting[termStart[ii] + k] - val;

         // cout << "red: " << red << " ";
         // cout << "n_curr->red: " << n_curr->red << "\n";
      
         diff = red - n_curr->red;

         if(diff > PLUSZERO || diff < MINUSZERO)
         {
            cout << "dbg_dirRed:  ERROR -- Reduced Cost!! \n\n";
         }
         n_curr = n_curr->fNext;
         cnt++;

         // cout << "\n";
      }
      // cout << "\n";
   }
   // cout << "\n";
}

void simplex::info_mv()
{
#ifdef compile4phc
   int fail = demics_store_mixed_volume(mixedVol);
#else
   cout.precision(15);
   cout << "# Mixed Cells: " << mixedCell << "\n";
   cout << "Mixed Volume: " << mixedVol << "\n\n";
#endif
}

void simplex::info_allSup()
{
   int i, j;

   cout << "<< Support Set >>\n\n";

   for(i = 0; i < supN; i++)
   {
      cout << "---- Level: " << i << " ----\n\n";

      for(j = 0; j < termSet[i]; j++)
      {
         cout << "* FrIdx: " << j << "\n";

         Supp[i][j].info_sup();

         cout << "\n";
      }
   }
   cout << "-- AuxMat -- \n";
   Supp[supN][0].info_sup();

   cout << "\n";
}

void simplex::info_allCostVec()
{
   int i, j;

   cout << "<< Cost Vector >>\n\n";

   for(i = 0; i < supN; i++)
   {
      cout << "---- Level: " << i << " ----\n\n";

      for(j = 0; j < termSet[i]; j++)
      {

         cout << "* FrIdx: " << j << "\n";

         Supp[i][j].info_costVec();

         cout << "\n";
      }
   }
   cout << "\n";
}

void simplex::info_lifting()
{
   int i, j, k, counter = 0;

   cout << "\nLifting: \n";
   for(i = 0; i < supN; i++)
   {
      for(j = termStart[i]; j < termStart[i] + termSet[i]; j++)
      {
         cout << lifting[counter] << " ";
         counter++;
      }
      cout << "\n";
   }
   cout << "\n";
  
   for(k = 0; k < supN; k++)
   {
      cout << "level: " << k << "\n";

      for(j = 0; j < termSet[k]; j++)
      {
         cout << "free index: " << j << "\n";

         for(i = termStart[k]; i < termStart[k] + termSet[k]; i++)
         {
            if(i != termStart[k] + j)
            {
               cout << lifting[i] - lifting[j + termStart[k]] << " ";
            }
         }
         cout << "\n\n";
      }
      cout << "\n";
   }
}

void simplex::info_simplexData()
{
   info_invB();
   info_p_sol();
   info_d_sol();
   info_basisIdx();
   info_nbIdx(); 
   info_nf_pos();
}
