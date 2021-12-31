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

#include "mvc.h"

mvc::mvc(void)
{
   Dim = 0;
   supN = 0;
  
   row = 0;
   col = 0;

   termSumNum = 0;
   termMax = 0;
   maxLength = 0;

   total_iter = 0;
   total_feasLP = 0;
   total_LPs = 0;
   total_1PT = 0;
   total_2PT = 0;
   total_triLPs_mLP = 0;
   total_unbLP_tab = 0;

   lvl_1PT = NULL;
   lvl_2PT = NULL;
   actNode = NULL;

   mfNum = NULL;
  
   termSet = NULL;
   termStart = NULL;
   re_termStart = NULL;
   type = NULL;

   mRepN = NULL;
   mFeaIdx = NULL;
   mFea = NULL;

   trNeg = NULL;

   firIdx = NULL;
   repN = NULL;

   sp = NULL;

   candIdx = NULL;
   trMat = NULL;

   table = NULL;

   lv = NULL;
   iLv = NULL;
}

mvc::~mvc(void)
{
   int i;

   if(trNeg)
   {
      for(i = 0; i < termSet[sp[0]]; i++) delete [] trNeg[i];

      delete [] trNeg;
      trNeg = NULL;
   }
   delete [] re_termStart;
   delete [] mfNum;
   delete [] lvl_1PT;
   delete [] lvl_2PT;
   delete [] actNode;
   delete [] firIdx;
   delete [] repN;
   delete [] sp;
   delete [] candIdx;
   delete [] trMat;
   delete [] lv;
   delete [] iLv;
}

void mvc::getMemory ( int depth, int lvl, int length )
{
   int i;
   int elemLen, div = 1;

   // cout << "<< getMemory >> \n\n";

   if(depth != 0)
      div = 2;
   else
      div = 1;

   if(lvl != length - 1)
      elemLen = termSet[depth]  / ((lvl + 1) * div);
   else
      elemLen = 1;

   // cout << "depth: " << depth << " lvl: " << lvl << "\n";
   // cout << "elemLen: " << elemLen << "\n\n";

   for(i = 0; i < termSet[depth]; i++)
   {
      lv[depth].fTest[lvl].create_elem(row,col,termSet[depth],type[depth]);
      lv[depth].fTest[lvl].add_elem();
   }
   lv[depth].fTest[lvl].mark();
   lv[depth].fTest[lvl].cur = lv[depth].fTest[lvl].head;
}

void mvc::initMemoryCheck ( ftData& Data, int depth )
{
  
   int sn = sp[depth];
  
   if(Data.cur == NULL)
   {
      // cout << "Create Memory\n\n";
      Data.create_elem(row, col, termSet[sn], type[sn]);
      Data.add_elem();
   }
}

void mvc::memoryCheck ( ftData* Data, int depth )
{
   int sn = sp[depth];

   if((*Data).cur == NULL)
   {
      // cout << "Create Memory\n\n";
      (*Data).create_elem(row, col, termSet[sn], type[sn]);
      (*Data).add_elem();
   }
}

void mvc::get_candIdx ( inifData& curInif )
{
   int num = 0;

   uData *n_curr;
  
   // cout << "<< get_candIdx_pivInIdx >> \n\n";

   n_curr = curInif.fHead;

   while(n_curr != NULL)
   {
      // cout << n_curr->supLab + 1 << " : ";

      candIdx[num + 1] = n_curr->supLab;

      n_curr = n_curr->fNext;
      num++;
   }
   candIdx[0] = num;
}

int mvc::chooseSup
 ( int depth, theData* curNode, inifData* curInif, inifData* nextInif )
{
   int flag;

   // cout << "=== chooseSup === \n\n";

   /*
     cout << "<< Before updateDireRed >>\n";
     iLv[depth + 1].info_all_feasIdx();
    */

   switch(depth)
   {
      case 0:
         fUpdateDirRed(curInif, nextInif, curNode, iLv[depth].rsp, depth);
         break;
      default:
         updateDirRed(curInif, nextInif, curNode, iLv[depth].rsp, depth);
         break;
   }

   /*
     cout << "<< After updateDireRed >>\n";
     iLv[depth + 1].info_all_feasIdx();
     info_all_dirRed(depth, lv[sp[depth]].Node, iLv[depth + 1].inif);
    */

   switch(curNode->artV)
   {
      case 0:
         flag = findUnbDir(nextInif,curNode,iLv[depth + 1].rsp,
                           iLv[depth].rsp,depth);    
         break;
      case 1:
         flag = findUnbDir_art(nextInif,curNode,iLv[depth + 1].rsp,
                               iLv[depth].rsp,depth);    
         break;
   }

   /*
     cout << "<< After findUnbDir >>\n";
     iLv[depth + 1].info_all_feasIdx();
    */

   return (flag);
}

void mvc::fUpdateDirRed
 ( inifData* curInif, inifData* nextInif, theData* curNode, int* curRsp,
   int depth )
{
   int i, j, k;
   int num, nfPos, idx, fIdx, nfN, flag;
   int length, lvl, pivOutNum, colPos, rowPos;
   int *nf_pos, *pivOutList;

   double val, preRed;
   double *transRed;

   uData *c_curr, *n_curr;

   // cout << "<< fUpdateDirRed >>\n\n";
   // transMat = curNode->transMat_ptr;
   transRed = curNode->transRed_ptr;

   nf_pos = curNode->nf_pos_ptr;
   nfN = curNode->nfN;

   pivOutNum = curNode->pivOutNum;
   pivOutList = curNode->pivOutList;

   length = supN - depth - 1;
   fIdx = firIdx[depth];
   colPos = termStart[sp[depth]];

   memcpy(trMat, curNode->transMat_ptr, sizeof(double) * Dim * Dim);

   for(j = 0; j < Dim; j++)
   {
      trMat[j + j * Dim] -= 1;

      for(i = 0; i < Dim; i++)
      {
         trMat[i + j * Dim] *=  trNeg[fIdx][i];
      }
   }
   for(j = 0; j < length; j++)
   {
      lvl = curRsp[j];
      rowPos = termStart[lvl];

      c_curr = curInif[lvl].fHead;
      n_curr = nextInif[lvl].fHead;

      num = 0;

      while(c_curr != NULL)
      {
         flag = CONTINUE;
      
         for(i = 0; i < curNode->polyDim + 1; i++)
         {
            if(table_out(colPos + curNode->nodeLabel[i], 
                         rowPos + c_curr->supLab) == UNBOUNDED)
            {
               flag = UNBOUNDED;
               break;
            }
         }
         if(flag == CONTINUE)
         {
            // cout << "----- " << num + 1 << " -----\n\n";

            n_curr->supLab = c_curr->supLab;

            // For direction
      
            for(k = 0; k < nfN; k++)
            {
               val = 0;
               nfPos = nf_pos[k];

               for(i = 0; i < pivOutNum; i++)
               {
                  idx = pivOutList[i];

                  val += trMat[idx + nfPos * Dim] * c_curr->dir[idx];
	    
                  // cout << "idx: " << idx << "\n";
                  // cout << "transMat: " << transMat[idx + nfPos * Dim]
                  //  << "\n\n";
                  // cout << "c_curr-Elem.dir: " << c_curr->dir[idx] << "\n\n";
	       }
               n_curr->dir[nfPos] = val+trNeg[fIdx][nfPos]*c_curr->dir[nfPos];

               /*
                 cout << "val: " << val << "\n";
                 cout << "trNeg: " << trNeg[fIdx][nfPos] << "\n";
                 cout << "c_curr-Elem.dir: " << c_curr->dir[nfPos] << "\n";
                 cout << "Elem.dir: " <<  n_curr->dir[nfPos] << "\n\n";
                */
            }
            // cout << "\n";
            // For reduced cost
            val = 0;
            preRed = 0;

            for(i = 0; i < Dim; i++)
            {
               val -= trNeg[fIdx][i] * transRed[i] * c_curr->dir[i];
               preRed += trNeg[fIdx][i] * c_curr->dir[i];
            }
            n_curr->red = val - preRed + c_curr->red;
         }
         else
         {
            skipPtr(&n_curr, &nextInif[lvl].fHead);
         }
         c_curr = c_curr->fNext;
         n_curr = n_curr->fNext;
         num++;
      }
      // cout << "\n\n";
      // cout << "num: " << num << "\n\n";
      
      if(n_curr != NULL)
      {
         n_curr->prev->fNext = NULL;
      }
   }
}

void mvc::updateDirRed
 ( inifData* curInif, inifData* nextInif, 
   theData* curNode, int* curRsp, int depth )
{
   int i, j, k;
   int num, nfPos, idx, nfN, flag;
   int length, pivOutNum, lvl, colPos, rowPos;

   int *nf_pos, *pivOutList;

   double val;
   double *transRed;

   uData *c_curr, *n_curr;

   // cout << "<< updateDirRed >>\n\n";
   // transMat = curNode->transMat_ptr;
   transRed = curNode->transRed_ptr;

   nf_pos = curNode->nf_pos_ptr;
   nfN = curNode->nfN;

   pivOutNum = curNode->pivOutNum;
   pivOutList = curNode->pivOutList;

   length = supN - depth - 1;
   colPos = termStart[sp[depth]];

   memcpy(trMat, curNode->transMat_ptr, sizeof(double) * Dim * Dim);

   for(i = 0; i < Dim; i++)
   {
      trMat[i + i * Dim] -= 1;
      // transMat[i * (Dim + 1)] -= 1;
   }
   for(j = 0; j < length; j++)
   {
      lvl = curRsp[j];
      rowPos = termStart[lvl];

      c_curr = curInif[lvl].fHead;
      n_curr = nextInif[lvl].fHead;

      num = 0;

      while(c_curr != NULL)
      {
         flag = CONTINUE;

         for(i = 0; i < curNode->polyDim + 1; i++)
         {
            if(table_out(colPos + curNode->nodeLabel[i], 
		     rowPos + c_curr->supLab) == UNBOUNDED)
            {
               flag = UNBOUNDED;
               break;
            }
         }
         if(flag == CONTINUE)
         {
            // cout << "----- " << num + 1 << " -----\n\n";

	    n_curr->supLab = c_curr->supLab;
      
            // For direction

            for(k = 0; k < nfN; k++)
            {
               val = 0;
               nfPos = nf_pos[k];

               for(i = 0; i < pivOutNum; i++)
               {
                  idx = pivOutList[i];
                  val += trMat[idx + nfPos * Dim] * c_curr->dir[idx];
               }
               n_curr->dir[nfPos] = val + c_curr->dir[nfPos];
	  
               // cout << n_curr->dir[nfPos] << " ";
            }
            // cout << "\n";
            // For reduced cost
            val = 0;
            for(i = 0; i < pivOutNum; i++)
            {
               idx = pivOutList[i];
               val -= transRed[idx] * c_curr->dir[idx];
            }
            n_curr->red = val + c_curr->red;
         }
         else
         {
            skipPtr(&n_curr, &nextInif[lvl].fHead);
         }
         c_curr = c_curr->fNext;
         n_curr = n_curr->fNext;
         num++;
      }
      // cout << "\n";

      if(n_curr != NULL)
      {
         n_curr->prev->fNext = NULL;
      }
   }
}

int mvc::findUnbDir
 ( inifData* nextInif, theData* curNode,
   int* nextRsp, int* curRsp, int depth )
{
   int i;
   int nfN, flag, lvl, length, cnt;
   int feasNum, min_feasNum = BIGINT, min_lvl = 0;
   int *basisIdx, *nf_pos;

   uData *n_curr, *fHead, *cor_ptr;

#if DBG_FINDUNB
   cout << "<< findUnbDir >> \n\n";
#endif

   basisIdx = curNode->basisIdx_ptr;
   nf_pos = curNode->nf_pos_ptr;
   nfN = curNode->nfN;

   length = supN - depth - 1;
  
   for(i = 0; i < length; i++)
   {
      lvl = curRsp[i];

#if DBG_FINDUNB
      cout << "-------- Support: " << lvl + 1 << " --------\n";
#endif

      n_curr = (fHead = nextInif[lvl].fHead);

      feasNum = 0;
      while(n_curr != NULL)
      {

#if DBG_FINDUNB
         cout << "-- tarIdx " << n_curr->supLab + 1 << "-- \n";
#endif

         cor_ptr = nextInif[lvl].fHead;

         flag = checkDir(&cor_ptr, n_curr, n_curr->dir, n_curr->red, 
                         nf_pos, basisIdx, nfN);

         if(flag == UNB_TAR)
         {
	    skipPtr(&n_curr, &nextInif[lvl].fHead);

#if DBG_FINDUNB
            cout << "UNB_TAR\n\n";
#endif
	
         }
         else if(flag == UNB_COR)
         {
            skipPtr(&cor_ptr, &nextInif[lvl].fHead);

            feasNum++;

#if DBG_FINDUNB
            cout << "UNB_COR\n\n";
#endif

         }
         else
         {
            feasNum++;

#if DBG_FINDUNB
            cout << "CONTINUE\n\n";
#endif

         }
         n_curr = n_curr->fNext;
      }
      if(feasNum < min_feasNum)
      {
         min_feasNum = feasNum;
         min_lvl = lvl;
      }
   }

#if DBG_FINDUNB
   cout << "min_lvl: " << min_lvl + 1 << "\n";
   cout << "min_feasNum: " << min_feasNum << "\n\n";
#endif

   sp[depth + 1] = min_lvl;

   cnt = 0;

   for(i = 0; i < length; i++)
   {
      if(curRsp[i] != min_lvl)
      {
         nextRsp[cnt] = curRsp[i];
         cnt++;
      }
   }
   if(min_feasNum <= 1)
      return (STOP);
   else
      return (CONTINUE);
}

int mvc::findUnbDir_art
 ( inifData* nextInif, theData* curNode,
   int* nextRsp, int* curRsp, int depth )
{
   int i;
   int nfN, flag, lvl, length, cnt;
   int feasNum, min_feasNum = BIGINT, min_lvl = 0;

   int *basisIdx, *nf_pos;

   uData *n_curr, *fHead, *cor_ptr;

#if DBG_FINDUNB
   cout << "<< findUnbDir_art >> \n\n";
#endif
  
   basisIdx = curNode->basisIdx_ptr;
   nf_pos = curNode->nf_pos_ptr;
   nfN = curNode->nfN;

   length = supN - depth - 1;

   for(i = 0; i < length; i++)
   {
      lvl = curRsp[i];

#if DBG_FINDUNB
      cout << "-------- Support: " << lvl + 1 << " --------\n";
#endif

      n_curr = (fHead = nextInif[lvl].fHead);
      feasNum = 0;

      while(n_curr != NULL)
      {

#if DBG_FINDUNB
         cout << "-- tarIdx " << n_curr->supLab + 1 << "-- \n";
#endif

         cor_ptr = nextInif[lvl].fHead;

         flag = checkDir_art(&cor_ptr, n_curr, n_curr->dir, n_curr->red, 
                             nf_pos, basisIdx, nfN);

         if(flag == UNB_TAR)
         {
            skipPtr(&n_curr, &nextInif[lvl].fHead);

#if DBG_FINDUNB
            cout << "UNB\n\n";
#endif

         }
         else if(flag == UNB_COR)
         {
            skipPtr(&cor_ptr, &nextInif[lvl].fHead);

            feasNum++;
         }
         else
         { 
            feasNum++;

#if DBG_FINDUNB
            cout << "CONTINUE\n\n";
#endif

         }
         n_curr = n_curr->fNext;

         // iLv[depth + 1].info_feasIdx(j);
      }
      if(feasNum < min_feasNum)
      {
         min_feasNum = feasNum;
         min_lvl = lvl;
      }
   }

#if DBG_FINDUNB
   cout << "min_lvl: " << min_lvl + 1 << "\n";
   cout << "min_feasNum: " << min_feasNum << "\n\n";
#endif

   sp[depth + 1] = min_lvl;

   cnt = 0;
   for(i = 0; i < length; i++)
   {
      if(curRsp[i] != min_lvl)
      {
         nextRsp[cnt] = curRsp[i];
         cnt++;
      }
   }
   if(min_feasNum <= 1)
      return (STOP);
   else
      return (CONTINUE);
}

int mvc::checkDir
 ( uData** corPtr, uData* tarPtr, double* tar_dir, double tar_red, 
   int* nf_pos, int* basisIdx, int nfN )
{
   int i;
   int nfPos, ans, flag, sign;

   while(*corPtr != NULL)
   {
      if(*corPtr != tarPtr)
      {

#if DBG_FINDUNB
         cout << (*corPtr)->supLab + 1 << " : ";
         cout << (*corPtr)->red - tar_red << " : ";
#endif

         sign = checkSign_red((*corPtr)->red, tar_red);

         if(sign == NEGATIVE)
         {
            flag = UNBOUNDED;
	
            for(i = 0; i < nfN; i++)
            {
               nfPos = nf_pos[i];

#if DBG_FINDUNB
               cout << (*corPtr)->dir[nfPos] - tar_dir[nfPos] << " ";
#endif

               ans = checkNonNeg_dir((*corPtr)->dir[nfPos], tar_dir[nfPos]);

               if(ans == FALSE)
               {

#if DBG_FINDUNB
                  cout << "\n";
#endif

                  flag = CONTINUE;
                  break;
               }
            }
            if(flag == UNBOUNDED)
            {

#if DBG_FINDUNB	  
               cout << "\n";
#endif
	  
               return (UNB_TAR);
            }
         }
         else
         {
            flag = UNBOUNDED;
	
            for(i = 0; i < nfN; i++)
            {
               nfPos = nf_pos[i];

#if DBG_FINDUNB
               cout << (*corPtr)->dir[nfPos] - tar_dir[nfPos] << " ";
#endif

               ans = checkNonPos_dir((*corPtr)->dir[nfPos], tar_dir[nfPos]);

               if(ans == FALSE)
               {

#if DBG_FINDUNB
                  cout << "\n";
#endif

                  flag = CONTINUE;
                  break;
               }
            }
            if(flag == UNBOUNDED)
            {

#if DBG_FINDUNB	  
               cout << "\n";
#endif
               return (UNB_COR);
            }
         }
      }
      *corPtr = (*corPtr)->fNext;
   }

#if DBG_FINDUNB	  
   cout << "\n";
#endif

   return (CONTINUE);
}

int mvc::checkDir_art
 ( uData** corPtr, uData* tarPtr, double* tar_dir, double tar_red,
   int* nf_pos, int* basisIdx, int nfN )
{
   int i;
   int nfPos;
   int ans, flag, sign, cnt, nonNegVarNum;
  
   while(*corPtr != NULL)
   {
      if(*corPtr != tarPtr)
      {

#if DBG_FINDUNB
         cout << (*corPtr)->supLab + 1 << " : ";
#endif
      
         flag = CONTINUE;
      
         for(i = 0; i < nfN; i++)
         {
            nfPos = nf_pos[i];

            if(basisIdx[nfPos] >= termSumNum - supN)
            {
               ans = checkZero_dir((*corPtr)->dir[nfPos], tar_dir[nfPos]);

               if(ans == FALSE)
               {
                  flag = STOP;
                  break;
	       }
            }
         }
         if(flag == STOP) break;

         sign = checkSign_red((*corPtr)->red, tar_red);

         if(sign == NEGATIVE)
         {
            cnt = 0;
            nonNegVarNum = 0;
            flag = STOP;

            for(i = 0; i < nfN; i++)
            {
               nfPos = nf_pos[i];

#if DBG_FINDUNB
               cout << (*corPtr)->dir[nfPos] - tar_dir[nfPos] << " ";
#endif

               if(basisIdx[nfPos] < termSumNum - supN)
               {
                  nonNegVarNum++;

                  ans = checkNonNeg_dir((*corPtr)->dir[nfPos],tar_dir[nfPos]);

                  if(ans == FALSE)
                  {

#if DBG_FINDUNB
                     cout << "\n";
#endif
                     flag = STOP;
                     break;
                  }
                  else
                     flag = UNBOUNDED;
               }
            }
            if(flag == UNBOUNDED)
            {

#if DBG_FINDUNB	  
               cout << "\n";
#endif
               return (UNB_TAR);
            }
         }
         else
         {
            cnt = 0;
            nonNegVarNum = 0;
	
            flag = STOP;

            for(i = 0; i < nfN; i++)
            {
               nfPos = nf_pos[i];

#if DBG_FINDUNB
               cout << (*corPtr)->dir[nfPos] - tar_dir[nfPos] << " ";
#endif

               if(basisIdx[nfPos] < termSumNum - supN)
               {
                  nonNegVarNum++;

                  ans = checkNonPos_dir((*corPtr)->dir[nfPos],tar_dir[nfPos]);

                  if(ans == FALSE)
                  {

#if DBG_FINDUNB
                     cout << "\n";
#endif
                     flag = STOP;

                     break;
	          }
                  else
                     flag = UNBOUNDED;
	       }
	    }
            if(flag == UNBOUNDED)
            {

#if DBG_FINDUNB	  
               cout << "\n";
#endif
               return (UNB_COR);
            }
         }
      }
      *corPtr = (*corPtr)->fNext;
   }

#if DBG_FINDUNB	  
   cout << "\n";
#endif

   return(CONTINUE);
}

void mvc::skipPtr ( uData** curr, uData** fHead )
{
   if(*curr == *fHead)
   {
      *fHead = (*curr)->fNext;
   }
   else if((*curr)->fNext != NULL)
   {
      (*curr)->prev->fNext = (*curr)->fNext;
      (*curr)->fNext->prev = (*curr)->prev;
   }
   else
   {
      (*curr)->prev->fNext = (*curr)->fNext;
   }
}

void mvc::dbg_init_transMat ( theData* curNode )
{
  int i, j;

  double val;

  curNode->info_invB();
  curNode->info_transMat();

  for(i = 0; i < Dim; i++)
  {
    for(j = 0; j < Dim; j++)
    {
      val = curNode->transMat_out(i, j) - curNode->invB_out(i, j);

      if(val > PLUSZERO || val < MINUSZERO)
      {
	cout << "dbg_init_transMat:  ERROR !! \n\n";
        break;
      }
    }
  }
}

void mvc::dbg_transMat ( theData* preNode, theData* curNode )
{
   int i, j, k;
   int nfN, nfPos;
   int* nf_pos;
  
   double val, diff;

   nfN = curNode->nfN;
   nf_pos = curNode->nf_pos;

   for(k = 0; k < nfN; k++)
   {
      nfPos = nf_pos[k];

      for(j = 0; j < Dim; j++)
      {
         val = 0;

         for(i = 0; i < Dim; i++)
         {
            val += curNode->transMat_out(nfPos, i) * preNode->invB_out(i, j);
         }
         diff = curNode->invB_out(nfPos, j) - val;

         if(diff > PLUSZERO || diff < MINUSZERO)
         {
            cout << "dbg_transMat:  ERROR !! \n\n";
         }
      }
   }
}

void mvc::check_transMat ( theData* preNode, theData* curNode )
{
   int i, j, k;
  
   double val;

   cout << "<< check_transMat >> \n\n";

   // cout << "<< Parent >> \n";
   // preNode->info_invB();

   cout << "<< Cur >> \n";
   // curNode->info_invB();
   curNode->info_transMat();
   // curNode->info_nf_pos();

   cout << "<< From transMat >>\n";
   for(k = 0; k < Dim; k++)
   {
      for(j = 0; j < Dim; j++)
      {
         val = 0;

         for(i = 0; i < Dim; i++)
         {
            val += curNode->transMat_out(k, i) * preNode->invB_out(i, j);
         }

         if(val < PLUSZERO && val > MINUSZERO)
         {
            cout.setf(ios::right);
            cout.width(10);
            cout << "0 ";
         }
         else
         {
            cout.setf(ios::right);
            cout.width(10);
            cout << val << " ";
         }
      }
      cout << "\n";
   }
   cout << "\n";
}

void mvc::info_all_dirRed ( int depth, ftData* Node, inifData* nextInif )
{
   int j, k, num, nfPos, nfN;
   int *nf_pos;

   double val;
 
   uData *n_curr;

   cout << "<< info_all_dirRed >> \n\n";

   nf_pos = Node->parent->nf_pos_ptr;
   nfN = Node->parent->nfN;

   for(j = depth + 1; j < supN; j++)
   {
      n_curr = nextInif[j].fHead;

      num = 0;

      cout << "--- Support: " << j + 1 << " ---\n";

      while(n_curr != NULL)
      {
         // cout << num + 1 << " : ";
      
         cout << n_curr->supLab  + 1 << " : ";

         for(k = 0; k < nfN; k++)
         {
            val = 0;
            nfPos = nf_pos[k];

            val = n_curr->dir[nfPos];

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
         cout << " : " << n_curr->red;

         n_curr = n_curr->fNext;
         cout << "\n";

         // num++;
      }
      cout << "\n";
   }
}

void mvc::info_neg ( int termSet, int** negIdx )
{
   int i, j;

   cout << "<< trNeg >> \n";
   for(j = 0; j < termSet; j++)
   {
      for(i = 0; i < row; i++)
      {
         cout << trNeg[j][i] << " ";
      }
      cout << "\n";
   }
   cout << "\n\n";

   cout << "<< negIdx >> \n";
   for(j = 0; j < termSet; j++)
   {
      for(i = 0; i < negIdx[j][0]; i++)
      {
         cout << negIdx[j][i + 1] << " ";
      }
      cout << "\n";
   }
   cout << "\n\n";
}

void mvc::info_sp ( int depth )
{
   int i;

   cout << "sp: ";

   for(i = 0; i < depth; i++)
   {
      cout << sp[i] << " ";
   }
   cout << "\n\n";
}

void mvc::info_parent_node ( int depth )
{
   int i;

   cout << "Node: ";
   for(i = 0; i < depth; i++)
   {
      cout << sp[i] << " : ";

      lv[sp[i]].Node->info_parent_node();
   }
   cout << "\n\n";
}

void mvc::info_tuple ( int lvl, int depth )
{
   int i;

   cout << "( ";
    
   for(i = 0; i < lvl; i++)
   {
      cout << mFeaIdx[i][mRepN[i]] + 1 << " ";
   }
   cout << ") --\n\n";
}

void mvc::info_mFea ( int length )
{
   int i;

   cout << "mFea:  ";
   for(i = 0; i < length; i++)
   {
      cout << mFea[i] << " ";
   }
   cout << "\nmRepN: ";
   for(i = 0; i < length; i++)
   {
      cout << mRepN[i] << " ";
   }
   cout << "\n\n";
}

void mvc::info_firIdx ( int length )
{
   int i;

   cout << "<< firIdx >> \n";

   for(i = 0; i <= length; i++)
   {
      cout << firIdx[i] << " ";
   }
   cout << "\n\n";
}

void mvc::info_fIdx ( ftData* Data )
{
   cout << "First Index: ";
   cout << Data[0].parent->fIdx + 1 << "\n\n";
}

void mvc::info_candIdx()
{
   cout << "candIdx: ";

   for(int i = 0; i < candIdx[0]; i++)
   {
      cout << candIdx[i + 1] << " ";
   }
   cout << "\n\n";
}

void mvc::info_elemNum( int length, ftData* Data, ftData Node )
{
   int i;

   cout.width(15);
   cout << "numElem: ";
   for(i = 0; i < length - 1; i++)
   {
      cout << Data[i].elemNum << " ";
   }
   cout << Node.elemNum << " ";
   cout << "\n\n";
}

void mvc::info_prop_elemNum ( int length, ftData* Data, ftData Node )
{
   int i;

   cout.width(15);
   cout << "prop_numElem: ";

   for(i = 0; i < length - 1; i++) Data[i].info_numElem();
   Node.info_numElem();

   cout << "\n\n";
}

void mvc::info_table()
{
   int i, j, k, u_cnt = 0, t_cnt = 0;

   cout << "<< Relation table >>\n";

   for(i = 0; i < termSumNum; i++)
   {
      for(k = 0; k < i; k++) cout << "  ";

      for(j = i + 1; j < termSumNum; j++)
      {
         if(table_out(i, j) == UNBOUNDED) u_cnt++;

         cout << table_out(i, j) << " ";
         t_cnt++;
      }
      cout << "\n";
   }
   cout << "# Unb. LPs: " << u_cnt << "\n";
   cout << "# Elem.: " << t_cnt << "\n";
   cout << "Ratio: " << (double) u_cnt / (double) t_cnt << "\n\n";
}

void mvc::info_cpuTime ( double cpuTime_start, double cpuTime_end )
{
   double totalTime, second;
   int hour, minute;

   totalTime = (cpuTime_end - cpuTime_start) / (double) CLK_TCK;

   cout.precision(8);
   cout << "CPU time: " << totalTime << "s";

   if(totalTime >= 3600)
   {
      hour = (int)totalTime / 3600;
      minute = (int)(totalTime - hour * 3600) / 60;
      second = totalTime - hour * 3600 - minute * 60;

      cout << " ( " << hour << "h" << minute << "m" << second << "s )"
           << endl;
   }
   else if(totalTime >= 60 && totalTime < 3600)
   {
      minute = (int)(totalTime) / 60;
      second = totalTime - minute * 60;

      cout << " ( " << minute << "m" << second << "s )" << endl;
   }
   else
   {
      cout << endl;
   }
}

void mvc::info_final()
{
   int i;
   double total_LP_tab;
   double levelNode, totalNode = 0;

   cout << "----- Final Info. -----\n\n";

   total_LP_tab = termSumNum * (termSumNum - 1) * 1 / 2;

   cout.precision(4);
   cout << "(Unb. LPs / # Total LPs) at Table: " 
        <<  total_unbLP_tab / total_LP_tab << "\n\n";

   cout.precision(3);
   cout << "# LPs: " <<  total_LPs << "\n";

   cout.precision(3);
   cout << "# LPs at iLP: " <<  total_1PT << "\n";

   cout.precision(3);
   cout << "# LPs at mLP: " <<  total_2PT << "\n\n";

   cout.precision(3);
   cout << "# Feas. LPs: " <<  total_feasLP << "\n";

   cout.precision(3);
   cout << "# Tri. LPs at mLP: " <<  total_triLPs_mLP << "\n\n";

   cout.precision(4);
   cout << "Ave. Iter for Feas. LPs: " << total_iter / total_feasLP << "\n\n";

 /*
   cout.width(15);
   cout << "#LPs at 1PT:";

   for(i = 0; i < supN; i++)
   {
     cout.precision(3);
     cout.width(8);

     cout << lvl_1PT[i] << " ";
   }
   cout << "\n";

   cout.width(15);
   cout << "#LPs at 2PT:";

   for(i = 0; i < supN; i++)
   {
      cout.precision(3);
      cout.width(8);

      cout << lvl_2PT[i] << " ";
   }
   cout << "\n\n";
  */

   // cout.width(15);
   // cout << "#Nodes:";

   levelNode = actNode[0];

   // cout.precision(3);
   // cout.width(8);
   // cout << levelNode << " ";

   totalNode += levelNode;

   for(i = 0; i < supN - 1; i++)
   {
      // cout.precision(3);
      // cout.width(8);

      levelNode = (termSet[i + 1] - 1) * termSet[i + 1] * actNode[i] / 2;

      // cout << levelNode << " ";

      totalNode += levelNode;
   }
   // cout << "\n\n";

   cout.precision(3);
   cout << "Total nodes: " << totalNode << "\n\n";

   cout << "-----------------------\n\n";
}

/// for mixed volume computation ///////////////////// 

#ifdef compile4phc
void mvc::initialize_with_lifting
 ( dataSet& Data, double* lifvals, int seedNum, int output )
{
   int i, j;
   int length;

   row = (Dim = Data.Dim);
   supN = Data.supN;

   termSumNum = Data.termSumNum;
   termMax = Data.termMax;
   maxLength = Data.typeMax + 1;

   col = termSumNum - supN + Dim;

   termSet = Data.termSet;
   termStart = Data.termStart;

   type = Data.type;

   mfNum = new int [supN];
   assert(mfNum);

   lvl_1PT = new double [supN];
   assert(lvl_1PT);
   memset(lvl_1PT, 0, sizeof(double) * supN);

   lvl_2PT = new double [supN];
   assert(lvl_2PT);
   memset(lvl_2PT, 0, sizeof(double) * supN);

   actNode = new double [supN];
   assert(actNode);
   memset(actNode, 0, sizeof(double) * supN);

   firIdx = new int [supN + 1];
   assert(firIdx);
   memset(firIdx, 0, sizeof(int) * (supN + 1));

   re_termStart = new int [supN + 1];
   assert(re_termStart);
   re_termStart[0] = 0;

   repN = new int [supN];
   assert(repN);
   memset(repN, 0, sizeof(int) * supN);

   sp = new int [supN];
   assert(sp);

   candIdx = new int [termMax + 1];
   assert(candIdx);

   trMat = new double [Dim * Dim];
   assert(trMat);

   lv = new lvData [supN];
   assert(lv);

   iLv = new iLvData [supN];
   assert(iLv);

   for(i = 0; i < supN; i++)
   {
      re_termStart[i + 1] = termStart[i + 1] - i - 1;

      sp[i] = i;

      lv[i].create(i, supN, Dim, type[i] + 1, termMax);
      iLv[i].create(i, supN, Dim, termMax);

      length = type[i] + 1;
 
      for(j = 0; j < length; j++)
      {
         getMemory(i, j, length);
      }
   }
   Simplex.initialize_with_lifting(Data, lifvals, firIdx, seedNum, output);
   Reltab.allocateAndIni(Simplex, &firIdx, Dim, supN, termSumNum,
                         termSet, termStart, re_termStart);
 
   iLv[0].getInit(Data, Simplex.lifting, termSet, termStart, Dim, supN);

#if DBG_CHOOSESUP
   iLv[0].info_all_dirRed();
#endif

   // Simplex.info_allCostVec();
   // Simplex.info_lifting();
   // Simplex.info_allSup();
}
#endif

void mvc::allocateAndIni ( dataSet& Data, int seedNum, int output )
{
   int i, j;
   int length;

   row = (Dim = Data.Dim);
   supN = Data.supN;

   termSumNum = Data.termSumNum;
   termMax = Data.termMax;
   maxLength = Data.typeMax + 1;

   col = termSumNum - supN + Dim;

   termSet = Data.termSet;
   termStart = Data.termStart;

   type = Data.type;

   mfNum = new int [supN];
   assert(mfNum);

   lvl_1PT = new double [supN];
   assert(lvl_1PT);
   memset(lvl_1PT, 0, sizeof(double) * supN);

   lvl_2PT = new double [supN];
   assert(lvl_2PT);
   memset(lvl_2PT, 0, sizeof(double) * supN);

   actNode = new double [supN];
   assert(actNode);
   memset(actNode, 0, sizeof(double) * supN);

   firIdx = new int [supN + 1];
   assert(firIdx);
   memset(firIdx, 0, sizeof(int) * (supN + 1));

   re_termStart = new int [supN + 1];
   assert(re_termStart);
   re_termStart[0] = 0;

   repN = new int [supN];
   assert(repN);
   memset(repN, 0, sizeof(int) * supN);

   sp = new int [supN];
   assert(sp);

   candIdx = new int [termMax + 1];
   assert(candIdx);

   trMat = new double [Dim * Dim];
   assert(trMat);

   lv = new lvData [supN];
   assert(lv);

   iLv = new iLvData [supN];
   assert(iLv);

   for(i = 0; i < supN; i++)
   {
      re_termStart[i + 1] = termStart[i + 1] - i - 1;

      sp[i] = i;

      lv[i].create(i, supN, Dim, type[i] + 1, termMax);
      iLv[i].create(i, supN, Dim, termMax);

      length = type[i] + 1;
 
      for(j = 0; j < length; j++)
      {
         getMemory(i, j, length);
      }
   }
   Simplex.allocateAndIni(Data, firIdx, seedNum, output);
   Reltab.allocateAndIni(Simplex, &firIdx, Dim, supN, termSumNum,
                         termSet, termStart, re_termStart);
 
   iLv[0].getInit(Data, Simplex.lifting, termSet, termStart, Dim, supN);

#if DBG_CHOOSESUP
   iLv[0].info_all_dirRed();
#endif

   // Simplex.info_allCostVec();
   // Simplex.info_lifting();
   // Simplex.info_allSup();
}

void mvc::initCheck ( int depth, ftData& Data )
{
   int i, j, k;
   int idx_one, negNum, sn, feaNum = 0;

   int** negIdx;

   double elem;
   double* val;

   sn = sp[depth];

   val = new double [termSet[sn] - 1];
   assert(val);

   negIdx = new int* [termSet[sn]];
   assert(negIdx);
  
   trNeg = new int* [termSet[sn]];
   assert(trNeg);
  
   for(i = 0; i < termSet[sn]; i++)
   {
      negIdx[i] = new int [Dim + 1];
      assert(negIdx[i]);

      trNeg[i] = new int [Dim];
      assert(trNeg[i]);
   }
   srand(12); // srandom(12);
   for(i = 0; i < termSet[sn] - 1; i++)
   {
      val[i] = (double) rand() / (double) RAND_MAX;
   }
   // cout << "Val: ";
   // for(i = 0; i < termSet[sn] - 1; i++) cout << val[i] << " ";
   // cout << "\n\n";

   firIdx[supN] = 0;

   for(idx_one = 0; idx_one < termSet[sn]; idx_one++)
   {

#if DBG_NODE
    cout << "------------------- Idx: " << idx_one + 1
         << "-------------------\n\n";
#endif

      // Data.create_elem(row, col, termSet[depth], type[depth]);

      initMemoryCheck(Data, depth);

      Data.cur->clear();

      firIdx[sn] = idx_one;

      Simplex.get_iNbN_nfN(&Data.cur, termSet[sn] - 1 + Dim, Dim);

      negNum = 0;

      for(i = 0; i < Dim; i++)
      {
         elem = 0;

         for(j = 0; j < termSet[sn] - 1; j++)
         {
            elem += val[j] * Simplex.put_elem_supp(sn, idx_one, i, j);
         }
         if(elem < MINUSZERO)
         {
            Data.cur->p_sol[termSumNum - supN + i] = -1 * elem;

            negIdx[idx_one][negNum + 1] = i;
            trNeg[idx_one][i] = -1;

            negNum++;

            for(k = 0; k < termSet[sn] - 1; k++)
            {
               Simplex.mult_elem_supp(sn, idx_one, i, k);
            }
         }
         else if(elem > PLUSZERO)
         {
            Data.cur->p_sol[termSumNum - supN + i] = elem;
	    trNeg[idx_one][i] = 1;
         }
         else
         {
            Data.cur->p_sol[termSumNum - supN + i] = 0;
            trNeg[idx_one][i] = 1;
         }
      }
      negIdx[idx_one][0] = negNum;

      Data.make_init_data(termSumNum, supN, termSet[sn], re_termStart[sn]);
    
      initLP(Data, negIdx, depth, idx_one, feaNum);
   }
   // info_neg(termSet[sn], negIdx);

   delete [] val;
   val = NULL;

   if(negIdx)
   {
      for(i = 0; i < termSet[sn]; i++) delete [] negIdx[i];

      delete [] negIdx;
      negIdx = NULL;
   }
}

void mvc::initLP
 ( ftData& Data, int** negIdx, int depth, int idx, int& feaNum )
{
   int j, k;
   int flag, sn, iter;

#if DBG_NODE
  cout << "<< initLP >>\n\n";
#endif

   sn = sp[depth];
  
   Simplex.get_cur(&Data.cur);
   Simplex.copy_p1_d_sol(Data.cur);

   iter = 0;

   flag = Simplex.fSolLP(termSet[sn], re_termStart[sn], iter);

   total_LPs++;
   total_1PT++;
   lvl_1PT[depth]++;

   if(flag == OPT)
   {

#if DBG_FEA
    cout << "OPT \n";
#endif

      total_iter += (double) iter;
      total_feasLP++;

      Data.cur->joint();
      Data.cur->fIdx = idx;

      Simplex.get_res(Data);
      Simplex.get_pivOutNum(&Data.cur);
    
      mFeaIdx[depth][feaNum] = idx;
      mFea[depth]++;

      feaNum++;

      for(j = 0; j < negIdx[idx][0]; j++)
      {
         for(k = 0; k < termSet[sn] - 1; k++)
         {
            Simplex.mult_elem_supp(sn, idx, negIdx[idx][j + 1], k);
         }
         for(k = 0; k < Dim; k++)
         {
            Data.cur->invB[negIdx[idx][j + 1] + k * Dim] *= -1;
         }
         Data.cur->d_sol[negIdx[idx][j + 1]] *= -1;
      }
      // Data.add_elem();
    
#if DBG_INI_CUR_INFO 
    Data.info_cur();
#endif

      Data.cur = Data.cur->next;
   }
   else if(flag == UNBOUNDED)
   {
      // Data.delete_cur();

      Data.init_info();

#if DBG_FEA
    cout << "UNB \n";
#endif

   }
   else
   {
      cout << "Error: too many iterations at initLP\n\n";
      cout << "( " << idx << " ) \n\n";

      exit(EXIT_FAILURE);
   }
}

// flag -- FNN (Nodeを出力して, findeNextNode)
//      -- STOP (ループをぬける)
int mvc::findNode ( int depth, int& lvl, int& feaNum, ftData* Data )
{
   int polyDim, flag, sn;
 
   ftData *pre, *cur;

#if DBG_SUC
   cout << "+++++++++++++++++<< findNode >>+++++++++++++++++\n";
   info_parent_node(depth);
#endif

   sn = sp[depth];
   polyDim = type[sn];

   pre = &Data[lvl - 1];
   cur = &Data[lvl];

   while(1)
   {

#if DBG_SUC
      cout << "++++++++++ lvl: " << lvl << " ++++++++++\n\n";
      info_tuple(lvl, depth);
#endif

      // cout << "------- depth: " << depth << " -------\n\n";
      // info_mFea(polyDim + 1);
    
      flag = mLP(pre, cur, Data, mFeaIdx[lvl - 1], mFeaIdx[lvl],
                 mRepN[lvl - 1], mRepN, mFea[lvl - 1], depth, mFea[lvl], 
                 lvl, polyDim + 1);

      // info_elemNum(polyDim + 1, Data, Node[depth]);
      // info_prop_elemNum(polyDim + 1, Data, Node[depth]);
            
      if(flag == NODE)
      {
         feaNum++;
         flag = FNN;
 
         break;
      }

      if(mFea[lvl] <  polyDim - lvl + 1 || lvl == polyDim)
      {  // SLIDE

#if DBG_SUC	
         cout << "-- SLIDE --\n\n";
#endif

         pre->next_data();

         if(lvl != polyDim) (*cur).delete_addedElem();

         cur->init_ptr();

         mRepN[lvl - 1]++;
         mFea[lvl] = 0;

         findUpNode(Data, &pre, &cur, lvl, polyDim, depth); // UP
      }
      else
      {  // DOWN

#if DBG_SUC
         cout << "-- DOWN --\n\n";
#endif

         lvl++;
	  
         pre = &Data[lvl - 1];
         cur = &Data[lvl];
      }
      if(lvl == 0)
      {
         flag = STOP;
         break;
      }
   }
   return (flag);
}

int mvc::findNextNode ( int depth, int& lvl, int& feaNum, ftData* Data )
{
   int polyDim, flag, sn;
 
   ftData *pre, *cur;

#if DBG_SUC
   cout << "+++++++++++++++++<< findNextNode >>+++++++++++++++++\n";
   info_parent_node(depth);
   cout << "++++++++++ lvl: " << lvl << " ++++++++++\n\n";
  info_tuple(lvl, depth);
#endif

   sn = sp[depth];
   polyDim = type[sn];

   // cout << "------- depth: " << depth << " -------\n\n";
   // info_mFea(polyDim + 1);

   flag = mLP(&Data[lvl - 1], &Data[lvl], Data, mFeaIdx[lvl - 1], mFeaIdx[lvl],
              mRepN[polyDim] + mRepN[lvl - 1], mRepN, mFea[lvl - 1], 
              depth, mFea[lvl], lvl, polyDim + 1);

   // info_elemNum(polyDim + 1, Data, Node[depth]);
   // info_prop_elemNum(polyDim + 1, Data, Node[depth]);

   if(flag == CONTINUE)
   {

#if DBG_SUC
      cout << "-- SLIDE or UP --\n\n";
#endif

      Data[lvl - 1].next_data();

      mRepN[lvl - 1]++;
      mRepN[lvl] = 0;
      mFea[lvl] = 0;

      findUpNode(Data, &pre, &cur, lvl, polyDim, depth); // UP

      if(lvl == 0)
      {
        flag = STOP;
      }
   }
   else
   {
      feaNum++;
      flag = FNN;
   }
   return(flag);
}

void mvc::findUpNode
 ( ftData* Data, ftData** pre, ftData** cur, int& lvl,
   int polyDim, int depth )
{
   while(1)
   {
      if(mFea[lvl - 1] - mRepN[lvl - 1] - 1 < polyDim - lvl + 1)
      { // UP

#if DBG_SUC
      cout << "-- UP --\n\n";
#endif
      
         mFea[lvl] = 0;
         mRepN[lvl - 1] = 0;

         lvl--;

         Data[lvl].delete_addedElem();
         Data[lvl].init_ptr();

         // Data[lvl].delete_all();
      
         if(lvl == 0)
         {
            mFea[lvl] = 0;
            mRepN[lvl] = 0;

#if DBG_SUC
         // info_sub(polyDim);
#endif

            break;
         }
         else
         {
            mFea[lvl] = 0;
            mRepN[lvl - 1]++;

            Data[lvl - 1].next_data();
	  
            (*pre) = &Data[lvl - 1];
            (*cur) = &Data[lvl];
         }
      }
      else
         break;
   }
}

int mvc::mLP
 ( ftData* Pre, ftData* Cur, ftData* Data, int* repIdx, int* feaIdx,
   int tarIdx, int* mRepN, int totalN, int depth, int& feaNum,
   int lvl, int length )
{
   int i;
   int idx2, flag, sn, iter;
   int lNbN, lNfN, fst_pivInIdx,  fst_sub_pivInIdx;
   int sub_firIdx, sub_tarIdx, colPos, rowStartPos;

   double fst_redCost;

   theData* target;

#if DBG_SUC
   cout << "---< mLP >---\n\n";
#endif

   sn = sp[depth];
  
   sub_firIdx = mRepN[0];
   sub_tarIdx = mRepN[lvl - 1];
  
   rowStartPos = termStart[sn];
   colPos = repIdx[sub_tarIdx] + termStart[sn];

   // info_fIdx(Data);

   // cout << "<< " << repIdx[sub_tarIdx] + 1 << " >>\n";
   // cout << "rowStartPos: " << rowStartPos << "\n\n";

   for(i = tarIdx + 1; i < totalN; i++)
   {

#if DBG_SUC
      cout << "-- ( " << repIdx[i] + 1 << " ) --\n\n";
#endif

#if DBG_S_PRE_INFO
      cout << "<< Pre_ptr >> \n";
      (*Pre).info_parent_ptr();
    
      cout << "<< Pre >> \n";
      (*Pre).info_parent_rIdx();
#endif

      memoryCheck(Cur, depth);
      // (*Cur).create_elem(row, col, termSet[depth], type[depth]);

      if(table_out(colPos, rowStartPos + repIdx[i]) == OPT)
      {
         if(lvl > 0)
         {
            get_firIdx(Data[0], Data[1], sn, lvl);
         }
         target = Pre->parent;

         flag = checkBasis(target, repIdx[i]);

         Simplex.get_mNbN_nfN(target, &(*Cur).cur);
         target->put_info(repIdx[i] - 1, idx2, lNbN, lNfN);

         (*Cur).cur->sw = OFF;

         if(flag == CONTINUE && lvl == 1)
         {
            flag = checkAnotherBasis(repIdx[sub_firIdx],i-sub_firIdx,&target);
	
            if(flag == OPT)
            {
               Simplex.get_mNbN_nfN(target, &(*Cur).cur);
               target->put_info(repIdx[sub_firIdx], idx2, lNbN, lNfN);

               (*Cur).cur->sw = ON;

               firIdx[sn] = repIdx[i];
            }
         }
         // info_firIdx(depth);

         if(flag == OPT)
         {

#if DBG_FEA
            cout << "OPT-1\n";
#endif

            total_triLPs_mLP++;
            actNode[depth]++;

          /*
            (*Cur).copy(col, target);
            (*Cur).mCopy(col, lNfN, idx2, termSet[sn], target);
            (*Cur).copy_pivOutIdx(target);
           */

            feaIdx[feaNum] = repIdx[i];
            feaNum++;

            (*Cur).cur->fIdx = repIdx[i];

            (*Cur).mGetPtr(target);

            (*Cur).get_nf_pos(target, lNfN, idx2);
            (*Cur).cur->mJoint();

            (*Cur).copy_rIdx(target, termSet[sn]);
            (*Cur).copy_pivOutIdx(target);      

	// (*Cur).add_elem();

#if DBG_S_CUR_INFO
         cout << "<< Cur_ptr >> \n";
         (*Cur).info_cur_ptr();
         cout << "<< Cur >> \n";
	 (*Cur).info_cur_rIdx();
#endif

         if(lvl == length - 1)
         {
            get_tuple_index(lv[sn].Node, Data, length);

            if(depth == supN - 1)
            {
               // info_parent_node(supN - 1);
               Simplex.calMixedVol(lv, sp, supN);
            }
            mRepN[lvl] += (i - tarIdx);

            (*Cur).cur = (*Cur).cur->next;

            return (NODE);
         }
         else
         {
            (*Cur).cur = (*Cur).cur->next;
         }
      }
      else
      {
         (*Cur).copy_rIdx((*Pre).parent, termSet[sn]);
         (*Cur).copy_pivOutIdx((*Pre).parent);

         Simplex.get_parent((*Pre).parent);
         Simplex.get_cur(&(*Cur).cur);

         fst_sub_pivInIdx = -1 * idx2 - 1;
         fst_pivInIdx = Pre->parent->nbIdx_ptr[fst_sub_pivInIdx];
         fst_redCost = Pre->parent->redVec_ptr[fst_pivInIdx];

         iter = 0;

         flag = Simplex.solLP(depth,fst_pivInIdx,fst_sub_pivInIdx,fst_redCost,
                              MCHECK,termSet[sn],re_termStart[sn],lNbN,iter);

         total_LPs++;
         total_2PT++;
         lvl_2PT[depth]++;

         if(flag == OPT)
         {

#if DBG_FEA	  
            cout << "OPT-2\n";
#endif
	  
            total_iter += (double) iter;
            total_feasLP++;
            actNode[depth]++;

            Simplex.get_pivOutNum(&(*Cur).cur);
	
            (*Cur).cur->joint();	
            (*Cur).decrease_nfN();
            (*Cur).cur->fIdx = repIdx[i];

            feaIdx[feaNum] = repIdx[i];
            feaNum++;

            // (*Cur).add_elem();

#if DBG_S_CUR_INFO 
            cout << "<< Cur >> \n";
            (*Cur).info_cur();
#endif

            if(lvl == length - 1)
            {
               get_tuple_index(lv[sn].Node, lv[sn].fTest, length);

               if(depth == supN - 1)
               {
                  // info_parent_node(supN - 1);
                  Simplex.calMixedVol(lv, sp, supN);
               }
               mRepN[lvl] += (i - tarIdx);

               (*Cur).cur = (*Cur).cur->next;
	  
               return (NODE);
            }
            else
            {
               (*Cur).cur = (*Cur).cur->next;
            }
         }
         else if(flag == UNBOUNDED)
         {

#if DBG_FEA
            cout << "UNB-1\n";
#endif

            (*Cur).init_info();
	
            // (*Cur).delete_cur();
         }
         else
         {
            cout << "Error: too much iterations at solLP\n\n";
	
            info_parent_node(depth);
            info_tuple(lvl, depth);
	
            cout << "( " << repIdx[i] + 1 << " ) \n\n";
	  
            exit(EXIT_FAILURE);
         }
      }
   }
   else
   {
      
#if DBG_FEA
       cout << "UNB-table\n";
#endif

       (*Cur).init_info();

   }
  }
   // (*Pre).next_data();

   return (CONTINUE);
}

int mvc::checkBasis ( theData* target, int sub_sIdx )
{
   if(target->rIdx[sub_sIdx - 1] >= 0)
      return (OPT);
   else
      return (CONTINUE);
}

int mvc::checkAnotherBasis ( int repIdx, int dist, theData** target )
{
   // cout << "===========<< checkAnotherBasis >>===========\n\n";
  
   for(int i = 0; i < dist; i++)
   {
      *target = (*target)->next;
   }
   // (*target)->info_rIdx();

   if((*target)->rIdx[repIdx] >= 0)
      return (OPT);
   else
      return (CONTINUE);
}

void mvc::get_firIdx ( ftData data_a, ftData data_b, int sn, int lvl )
{
   int sw;

   // cout << "<< get_firIdx3 >>\n\n";

   if(lvl == 1)
   {
      firIdx[sn] = data_a.parent->fIdx;
   }
   else
   {
      sw = data_b.parent->sw;

      if(sw == OFF)
      {
         firIdx[sn] = data_a.parent->fIdx;
      }
      else
      {
         firIdx[sn] = data_b.parent->fIdx;
      }
   }
}

void mvc::get_tuple_index ( ftData* Node, ftData* Data, int length )
{
   int i;
   int tmpIdx;

   // cout << "<< get_tuple_index >>\n\n";

   for(i = 0; i < length - 1; i++)
   {
      Node->parent->nodeLabel[i] = Data[i].parent->fIdx;
   }
   Node->parent->nodeLabel[length - 1] = Data[length - 1].cur->fIdx;

   if(Data[1].parent->sw == ON)
   {
      tmpIdx = Node->parent->nodeLabel[1];

      Node->parent->nodeLabel[1] = Node->parent->nodeLabel[0];
      Node->parent->nodeLabel[0] = tmpIdx;
   }
}

void mvc::initFeasTest ( int depth )
{
   int lvl, sn, length, feaNum, flag;

   sn = sp[depth];
   length = type[sn] + 1;

   lv[sn].get_info(&mRepN, &mFeaIdx, &mFea);

   lvl = 0;
   initCheck(depth, lv[sn].fTest[lvl]);

   lvl++;
   feaNum = 0;

   flag = CONTINUE;

   if(flag == CONTINUE)
   {
      flag = findNode(depth, lvl, feaNum, lv[sn].fTest);
   }
   else
   {
      flag = findNextNode(depth, lvl, feaNum, lv[sn].fTest);
   }
   // dbg_init_transMat(lv[sn].Node->parent);
}

int mvc::feasTest ( int depth, theData* parent )
{
   int lvl, length, feaNum, sn, flag;

#if DBG_SUC
   cout << "========== feasTest ==========\n\n";
#endif

   sn = sp[depth];

   length = type[sn] + 1;

   lv[sn].get_info(&mRepN, &mFeaIdx, &mFea);
  
   lvl = 0;
   flag = iCheck(depth, parent, lv[sn].fTest[lvl], iLv[depth].inif[sn]);

   lvl++;

   feaNum = 0;
   flag = CONTINUE;

   if(flag == CONTINUE)
   {
      flag = findNode(depth, lvl, feaNum, lv[sn].fTest);
   }
   else
   {
      flag = findNextNode(depth, lvl, feaNum, lv[sn].fTest);
   }
   return (flag);
}

int mvc::upFeasTest ( int& depth )
{
   int lvl, length, feaNum, flag;

#if DBG_SUC
   cout << "========== UP FeasTest ==========\n\n";
#endif

   while(1)
   {
      // lv[depth - 1].Node->delete_all();
      // iLv[depth - 1].info_rsp();

      iLv[depth].init(supN, depth - 1, iLv[depth - 1].rsp);

      // iLv[depth].info_all_feasIdx();
      // info_all_dirRed(depth - 1, lv[sp[depth - 1]].Node, lv[depth].inif);

      lv[sp[depth - 1]].Node->delete_addedElem();
      lv[sp[depth - 1]].Node->init_ptr();

      lv[sp[depth]].init_ptr();

      lv[sp[depth - 1]].get_info(&mRepN, &mFeaIdx, &mFea);

      length = type[sp[depth - 1]] + 1;

      feaNum = 0;
      lvl = type[sp[depth - 1]];
    
      flag = findNextNode(depth - 1, lvl, feaNum, lv[sp[depth - 1]].fTest);

      // STOP or FNN

      if(flag == CONTINUE)
      {
         flag = findNode(depth - 1, lvl, feaNum, lv[sp[depth - 1]].fTest);
      }
      // STOP or FNN

      depth--;

      if(flag == FNN || depth == 0) break;
   }
   return (flag);
}

void mvc::findMixedCell ( int depth, theData* parent )
{
   int lvl, length, feaNum, sn, flag;

#if DBG_SUC
   cout << "========== findMixedCell ==========\n\n";
#endif

   sn = sp[depth];
   length = type[sn] + 1;

   lv[sn].get_info(&mRepN, &mFeaIdx, &mFea);
  
   lvl = 0;
   flag = iCheck(depth, parent, lv[sn].fTest[lvl], iLv[depth].inif[sn]);

   lvl++;

   feaNum = 0;
   flag = CONTINUE;

   while(1)
   {
      if(flag == CONTINUE)
      {
         flag = findNode(depth, lvl, feaNum, lv[sn].fTest);

         if(flag == STOP) break;
      }
      else
      {
         flag = findNextNode(depth, lvl, feaNum, lv[sn].fTest);

         if(flag == STOP) break;
      }
   }
}

void mvc::findAllMixedCells ( int depth )
{
   int lvl, length, feaNum, flag;

   length = type[depth] + 1;

   lv[depth].get_info(&mRepN, &mFeaIdx, &mFea);

   lvl = 0;
   initCheck(depth, lv[depth].fTest[lvl]);

   lvl++;
   feaNum = 0;
   flag = CONTINUE;

   while(1)
   {
      if(flag == CONTINUE)
      {
         flag = findNode(depth, lvl, feaNum, lv[depth].fTest);

         lv[depth].Node->delete_addedElem();
         lv[depth].Node->init_ptr();
         // lv[depth].Node->delete_all();

         if(flag == STOP) break;
      }
      else
      {
         flag = findNextNode(depth, lvl, feaNum, lv[depth].fTest);

         lv[depth].Node->delete_addedElem();
         lv[depth].Node->init_ptr();
         // lv[depth].Node->delete_all();

         if(flag == STOP) break;
      }
   }
}

int mvc::iCheck
 ( int depth, theData* parent, ftData& Data, inifData& curInif )
{
   int feaNum = 0;
   int preNbN, sn, flag;
   int fst_pivInIdx, sub_fst_pivInIdx;

   uData* curr;
  
#if DBG_NODE
   cout << "+++++++++++++++++<< iCheck >>+++++++++++++++++\n";
   info_parent_node(depth);
#endif

#if DBG_PRE_INFO
   cout << "<< PRE >>\n";

   (*parent).info_p_sol_ptr();
   (*parent).info_d_sol_ptr();
   (*parent).info_invB_ptr();
   (*parent).info_basisIdx_ptr();
   (*parent).info_nf_pos_ptr();
   (*parent).info_nbIdx_ptr();
   (*parent).info_redVec_ptr();
#endif

   flag =  (*parent).artV;

   sn = sp[depth];

   switch(flag)
   {
      case TRUE:

#if DBG_NODE
         cout << "<< Art >> \n";
#endif

         get_candIdx(curInif);

         preNbN = (*parent).nbN;
         curr = curInif.fHead;

         while(curr != NULL)
         {

#if DBG_NODE
            cout << "------------------- Idx: " << curr->supLab + 1
                 << " -------------------\n";
#endif

            iLP_Art(parent,Data, depth,curr->supLab,fst_pivInIdx,
                    sub_fst_pivInIdx,preNbN,feaNum);
            curr = curr->fNext;
         }

#if DBG_NODE
         cout << "\n";
#endif

         if(feaNum <= 1)
         {
            mFea[0] = 0;

            // Data.delete_all();
            Data.delete_addedElem();
            Data.init_ptr();
 
            return (SLIDE);
         }
         else
            return (CONTINUE);

      case FALSE:

#if DBG_NODE
         cout << "<< NoArt >> \n";
#endif

         Simplex.fstRed_candIdx(curInif,&candIdx,fst_pivInIdx,
                                sub_fst_pivInIdx);

         // cout << "fst_pivInIdx: " << fst_pivInIdx << "\n";
         // cout << "sub_fst_pivInIdx: " << sub_fst_pivInIdx << "\n\n";

         preNbN = (*parent).nbN;
         curr = curInif.fHead;
  
         while(curr != NULL)
         {

#if DBG_NODE
            cout << "------------------- Idx: " << curr->supLab + 1
                 << " -------------------\n";
#endif

            iLP(parent,Data,depth,curr->supLab,fst_pivInIdx,
                sub_fst_pivInIdx,preNbN, feaNum);

            curr = curr->fNext;
         }

#if DBG_NODE
         cout << "\n";
#endif

      }
      if(feaNum <= 1)
      {
         mFea[0] = 0;

         // Data.delete_all();

         Data.delete_addedElem();
         Data.init_ptr();

         return (SLIDE);
      }
      else
         return (CONTINUE);
}

void mvc::iLP
 ( theData* parent, ftData& Data, int depth, int idx_one, int fst_pivInIdx,
   int sub_fst_pivInIdx, int preNbN, int& feaNum )
{
   int flag, fst_pivInIdx2, sub_fst_pivInIdx2, sn, iter;
   int termS, reTermS, repIdx, lNfN;

   double fst_redCost;

   // Data.create_elem(row, col, termSet[depth], type[depth]);

   sn = sp[depth];

   initMemoryCheck(Data, depth);

   repIdx = idx_one;
   firIdx[sn] = repIdx;
  
   termS = termStart[sn];
   reTermS = re_termStart[sn];

   // info_firIdx(depth);

   lNfN = (*parent).nfN;

   Simplex.get_iNbN_nfN(&Data.cur, preNbN + candIdx[0] - 1, lNfN);

   if(repIdx < fst_pivInIdx)
   {
      fst_pivInIdx2 = reTermS + fst_pivInIdx - 1;
      sub_fst_pivInIdx2 = preNbN - Dim + sub_fst_pivInIdx - 1;
   }
   else if(repIdx > fst_pivInIdx)
   {
      fst_pivInIdx2 = reTermS + fst_pivInIdx;
      sub_fst_pivInIdx2 = preNbN - Dim + sub_fst_pivInIdx;
   }
   if(repIdx != fst_pivInIdx)
   {
      Data.init_info();
      Data.create_rIdx(preNbN, repIdx, candIdx);
    
      Simplex.get_repIdx_candIdx(candIdx, repIdx);
      Simplex.get_parent(parent);
      Simplex.get_cur(&Data.cur);

      fst_redCost = Simplex.put_redCost(fst_pivInIdx);
    
      iter = 0;

      flag = Simplex.solLP(depth,fst_pivInIdx2,sub_fst_pivInIdx2,fst_redCost, 
                           ICHECK,termSet[sn],reTermS,preNbN,iter);

      total_LPs++;
      total_1PT++;
      lvl_1PT[depth]++;

      if(flag == OPT)
      {

#if DBG_FEA
      cout << "OPT-1\n";
#endif

         total_iter += (double) iter;
         total_feasLP++;

         Simplex.get_pivOutNum(&Data.cur);

         Data.cur->joint();
         Data.cur->fIdx = idx_one;

#if DBG_CUR_INFO
         cout << "<< Cur >>\n";
         Data.info_cur();
#endif

         // Data.add_elem();

         mFeaIdx[0][feaNum] = repIdx;
         mFea[0]++;

         feaNum++;

         Data.cur = Data.cur->next;

         // initMemoryCheck(Data, depth);
      }
      else if(flag == UNBOUNDED)
      {

#if DBG_FEA
         cout << "UNB-1\n";
#endif

         // Data.delete_cur();
      }
      else
      {
         cout << "Error: too many iterations at iLP\n";

         info_parent_node(depth);
         cout << "( " << idx_one + 1 << " ) \n\n";

         exit(EXIT_FAILURE);
      }
   }
   else
   {

#if DBG_FEA
      cout << "OPT-2\n";
#endif
    
      // Data.get_ptr(parent);
      // Data.copy(preNbN, parent);

      mFeaIdx[0][feaNum] = repIdx;
      mFea[0]++;

      feaNum++;

      Data.cur->fIdx = idx_one;

      Simplex.copy_eye(&Data.cur); 
      Simplex.cal_redVec(termSet[sn], reTermS, fst_pivInIdx, &Data.cur); 

      Data.iGetPtr(parent);
      Data.get_nbIdx_rIdx(preNbN, repIdx, candIdx, reTermS, parent);
      Data.init_info();

      Data.cur->iJoint();

#if DBG_CUR_INFO
      cout << "<< Cur_ptr >>\n";
      Data.info_cur_ptr();
      cout << "<< Cur >>\n";
      Data.info_cur_rIdx();
#endif

      // Data.add_elem();
      // initMemoryCheck(Data, depth);
    
      Data.cur = Data.cur->next;
   }
}

void mvc::iLP_Art
 ( theData* parent, ftData& Data, int depth, int idx_one, 
   int fst_pivInIdx, int sub_fst_pivInIdx, int preNbN, int& feaNum )
{
   int flag, termS, reTermS, repIdx, sn, iter, lNfN;

   sn = sp[depth];

   // Data.create_elem(row, col, termSet[depth], type[depth]);
   initMemoryCheck(Data, depth);

   repIdx = idx_one;
   firIdx[sn] = repIdx;

   termS = termStart[sn];
   reTermS = re_termStart[sn];

   // info_firIdx(depth);

   lNfN = (*parent).nfN;
   Simplex.get_iNbN_nfN(&Data.cur, preNbN + candIdx[0] - 1, lNfN);

   // Data.copy(preNbN, parent);

   Simplex.copy_eye(&Data.cur);
   Data.copy(col, parent);
   Data.iCopy(preNbN, lNfN, repIdx, termSet[sn], reTermS, candIdx, parent);
   Data.init_info();

   Simplex.get_cur(&Data.cur);

   iter = 0;

   flag = Simplex.solLP_art(depth,repIdx,fst_pivInIdx,preNbN,termSet[sn],
                            reTermS,iter);
   total_LPs++;
   total_1PT++;
   lvl_1PT[depth]++;

   if(flag == OPT)
   {

#if DBG_FEA
      cout << "OPT-1\n";
#endif

      total_iter += (double) iter;
      total_feasLP++;

      Data.cur->joint();
      Data.cur->fIdx = idx_one;

      Simplex.get_res(Data);
      Simplex.get_pivOutNum(&Data.cur);

      // Data.add_elem();

      mFeaIdx[0][feaNum] = repIdx;
      mFea[0]++;

      feaNum++;

      // Data.info_cur_transMat();

#if DBG_CUR_INFO
      cout << "Cur\n";
      Data.info_cur_ptr();
#endif

      Data.cur = Data.cur->next;
   }
   else if(flag == UNBOUNDED)
   {
      // Data.delete_cur();

#if DBG_FEA
      cout << "UNB-1\n";
#endif

   }
   else
   {
      cout << "Error: too much iterations at iLP_art\n";

      info_parent_node(depth);
      cout << "( " << idx_one + 1 << " ) \n\n";

      exit(EXIT_FAILURE);
   }
}

void mvc::Enum()
{
   int depth = 0, flag;

   Reltab.makeTable(total_unbLP_tab);
   table = Reltab.table;

   if(supN == 1)
   {
      findAllMixedCells(depth);
   }
   else
   {
      initFeasTest(depth);

      depth++;

      while(1)
      {
         if(depth == supN - 1)
         {
	    flag = chooseSup(depth - 1, lv[sp[depth - 1]].Node->parent,
                             iLv[depth - 1].inif, iLv[depth].inif);

            // info_sp(depth + 1);
            // iLv[depth].info_all_feasIdx();
            // info_all_dirRed(depth - 1, lv[sp[depth - 1]].Node,
            //                 iLv[depth].inif);

           /*
             iLv[depth].info_all_feasIdx();
             Simplex.dbg_dirRed(lv[sp[depth - 1]].Node->parent,
                                iLv[depth].inif, depth);
             Simplex.check_dirRed(lv[sp[depth - 1]].Node->parent, depth - 1);
             info_all_dirRed(depth - 1, lv[sp[depth - 1]].Node,
                             iLv[depth].inif);
            */
            if(flag == CONTINUE)
            {
               findMixedCell(depth, lv[sp[depth - 1]].Node->parent);

               flag = STOP;
            }
         }
         else
         {
           /*
             if(depth != 0)
             {
	        dbg_transMat(lv[sp[depth - 1]].Node->parent,
                             lv[sn].Node->parent);
                check_transMat(lv[sp[depth - 1]].Node->parent,
                               lv[sn].Node->parent);
             }
            */
            flag = chooseSup(depth - 1, lv[sp[depth - 1]].Node->parent,
                             iLv[depth - 1].inif, iLv[depth].inif);
            // info_sp(depth + 1);
            // iLv[depth].info_all_feasIdx();
            // info_all_dirRed(depth - 1, lv[sp[depth - 1]].Node,
            //                 iLv[depth].inif);
           /*
             iLv[depth].info_all_feasIdx();
             Simplex.dbg_dirRed(lv[sp[depth - 1]].Node->parent,
                                iLv[depth].inif, depth);
             Simplex.check_dirRed(lv[sp[depth - 1]].Node->parent, depth - 1);
             info_all_dirRed(depth - 1, lv[sp[depth - 1]].Node,
                             iLv[depth].inif);
            */

            if(flag == CONTINUE)
            {
               flag = feasTest(depth, lv[sp[depth - 1]].Node->parent);
            }
         }
         if(flag == STOP)
         {                              // SLIDE or UP
            flag = upFeasTest(depth);
	
            if(flag == STOP) break;
         }
         depth++; // DOWN
      }
   }
   Simplex.info_mv();
}
