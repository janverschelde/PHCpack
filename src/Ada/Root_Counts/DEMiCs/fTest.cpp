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

#include "fTest.h"

theData::theData()
{
   row = 0;
   col = 0;

   termS = 0;
  
   flag = 0;

   polyDim = 0;
  
   nbN = 0;
   nfN = 0;

   artV = 0;
   pivOutNum = 0;

   fIdx = 0;

   sw = 0;

   invB = NULL;
   transMat = NULL;
   transRed = NULL;

   p_sol = NULL;
   d_sol = NULL;

   redVec = NULL;

   basisIdx = NULL;
   nbIdx = NULL;

   nf_pos = NULL;
   rIdx = NULL;

   next = NULL;

   pivOutList = NULL;
   pivOutCheck = NULL;

  ////
   invB_ptr = NULL;  
   transMat_ptr = NULL;
   transRed_ptr = NULL;

   p_sol_ptr = NULL; 
   d_sol_ptr = NULL; 

   redVec_ptr = NULL;

   basisIdx_ptr = NULL; 
   nbIdx_ptr = NULL; 

   nf_pos_ptr = NULL;

  ////
   nodeLabel = NULL;
}

theData::~theData()
{
   delete [] invB;
   delete [] transMat;
   delete [] transRed;

   delete [] p_sol;
   delete [] d_sol;

   delete [] redVec;

   delete [] basisIdx;
   delete [] nbIdx;

   delete [] nf_pos;
   delete [] rIdx;

   delete [] pivOutList;
   delete [] pivOutCheck;

   delete [] nodeLabel;
}

void theData::create
 ( int ori_row, int ori_col, int ori_termS, int ori_polyDim )
{
   // cout << "<< create >> \n";
  
   row = ori_row;
   col = ori_col;

   termS = ori_termS;
   polyDim = ori_polyDim;

   invB = new double [row * row];
   assert(invB);
   memset(invB, 0, row * row * sizeof(double));

   transMat = new double [row * row];
   assert(transMat);
   memset(transMat, 0, row * row * sizeof(double));

   transRed = new double [row];
   assert(transRed);
   memset(transRed, 0, row * sizeof(double));

   p_sol = new double [col];
   assert(p_sol);
   memset(p_sol, 0, col * sizeof(double));

   d_sol = new double [row];
   assert(d_sol);
   memset(d_sol, 0, row * sizeof(double));

   basisIdx = new int [row];
   assert(basisIdx);
   memset(basisIdx, 0, row * sizeof(int));

   nf_pos = new int [row];
   assert(nf_pos);
   memset(nf_pos, 0, row * sizeof(int));

   nbIdx = new int [col];
   assert(nbIdx);
   memset(nbIdx, 0, col * sizeof(int));

   redVec = new double [col];
   assert(redVec);
   memset(redVec, 0, col * sizeof(double));
    
   rIdx = new int [termS];
   assert(rIdx);
   memset(rIdx, 0, termS * sizeof(int));

   pivOutList = new int [row];
   assert(pivOutList);
   memset(pivOutList, 0, row * sizeof(int));

   pivOutCheck = new int [row];
   assert(pivOutCheck);
   memset(pivOutCheck, 0, row * sizeof(int));

   nodeLabel = new int [polyDim + 1];
   assert(nodeLabel);
}

void theData::joint()
{
   invB_ptr = invB;  

   transMat_ptr = transMat; 
   transRed_ptr = transRed; 

   p_sol_ptr = p_sol; 
   d_sol_ptr = d_sol; 

   redVec_ptr = redVec; 

   basisIdx_ptr = basisIdx; 
   nbIdx_ptr = nbIdx;  

   nf_pos_ptr = nf_pos; 
}

void theData::iJoint()
{
   transMat_ptr = transMat;
   transRed_ptr = transRed;

   redVec_ptr = redVec;
   nbIdx_ptr = nbIdx; 
}

void theData::mJoint()
{
   nf_pos_ptr = nf_pos;
}

void theData::clear()
{
   nbN = 0;
   nfN = 0;
   artV = 0;

   memset(invB, 0, row * row * sizeof(double));
   memset(transMat, 0, row * row * sizeof(double));
   memset(p_sol, 0, col * sizeof(double));
   memset(d_sol, 0, row * sizeof(double));
   memset(basisIdx, 0, row * sizeof(int));
   memset(nf_pos, 0, row * sizeof(int));
   memset(nbIdx, 0, col * sizeof(int));
   memset(redVec, 0, col * sizeof(double));
   memset(rIdx, 0, termS * sizeof(int));
}

void theData::clear_transMat()
{
   memset(transMat, 0, row * row * sizeof(double));
}

void theData::put_info ( int repIdx, int& idx2, int& lNbN, int& lNfN )
{
   idx2 = rIdx[repIdx];

   lNbN = nbN;
   lNfN = nfN;
}

void theData::info_p_sol()
{
   int i;

   cout << "<< p_sol >> \n";

   for(i = 0; i < col; i++) cout << p_sol[i] << " ";
   cout << "\n\n";
}

void theData::info_d_sol()
{
   int i;

   cout << "<< d_sol >> \n";

   for(i = 0; i < row; i++) cout << d_sol[i] << " ";
   cout << "\n\n";
}

void theData::info_invB()
{
   int i, j;
   double val;

   cout << "<< invB >> \n";

   for(i = 0; i < row; i++)
   {
      for(j = 0; j < row; j++)
      {
         val = invB_out(i, j);

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
   cout << "\n\n";
}

void theData::info_transMat()
{
   int i, j;

   double val;

   cout << "<< transMat >> \n";

   //for(i = 0; i < nfN; i++){
   for(i = 0; i < row; i++)
   {
      // nfPos = nf_pos[i];

      for(j = 0; j < row; j++)
      {
         val = transMat_out(i, j);

         if(MINUSZERO < val && val < PLUSZERO)
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
   cout << "\n\n";
}

void theData::info_transRed()
{
   int i;
   double val;

   cout << "<< transRed >> \n";

   for(i = 0; i < row; i++)
   {
      val = transRed[i];

      if(MINUSZERO < val && val < PLUSZERO)
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
   cout << "\n\n";
}

void theData::info_basisIdx()
{
   int i;

   cout << "<< basisIdx >> \n";

   for(i = 0; i < row; i++) cout << basisIdx[i] << " ";
   cout << "\n\n";
}

void theData::info_nf_pos()
{
   int i;

   cout << "<< nf_pos >> \n";

   for(i = 0; i < nfN; i++) cout << nf_pos[i] << " ";
   cout << "\n\n";
}

void theData::info_nbIdx()
{
   int i;

   cout << "<< nbIdx >> \n";

   for(i = 0; i < col; i++) cout << nbIdx[i] << " ";
   cout << "\n\n";
}

void theData::info_redVec()
{
   int i;

   cout << "<< redVec >> \n";

   for(i = 0; i < col; i++)
   {
      // cout << i << " : " << redVec[i] << " | ";
      cout << redVec[i] << " ";
   }
   cout << "\n\n";
}

void theData::info_rIdx()
{
   int i;

   cout << "<< rIdx >> \n";
   for(i = 0; i < termS; i++) cout << rIdx[i] << " ";
   cout << "\n\n";
}

void theData::info_pivOutIdx()
{
   int i;

   cout << "pivOutCheck: ";
   for(i = 0; i < row; i++) cout << pivOutCheck[i] << " ";

   cout << "\n\npivOutList: ";
   for(i = 0; i < pivOutNum; i++) cout << pivOutList[i] << " ";

   cout << "\n\n";
}

void theData::info_p_sol_ptr()
{
   int i;

   cout << "<< p_sol_ptr >> \n";

   for(i = 0; i < col; i++) cout << p_sol_ptr[i] << " ";
   cout << "\n\n";
}

void theData::info_d_sol_ptr()
{
   int i;

   cout << "<< d_sol_ptr >> \n";

   for(i = 0; i < row; i++) cout << d_sol_ptr[i] << " ";
   cout << "\n\n";
}

void theData::info_invB_ptr()
{
   int i, j;
   double val;

   cout << "<< invB_ptr >> \n";

   for(i = 0; i < row; i++)
   {
      for(j = 0; j < row; j++)
      {
         val = invB_ptr_out(i, j);

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
   cout << "\n\n";
}

void theData::info_transMat_ptr()
{
   int i, j;

   double val;

   cout << "<< transMat_ptr >> \n";

   // for(i = 0; i < nfN; i++){
   for(i = 0; i < row; i++)
   {
      // nfPos = nf_pos[i];

      for(j = 0; j < row; j++)
      {
         val = transMat_ptr_out(i, j);

         if(MINUSZERO < val && val < PLUSZERO)
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
   cout << "\n\n";
}

void theData::info_transRed_ptr()
{
   int i;
   double val;

   cout << "<< transRed_ptr >> \n";

   for(i = 0; i < row; i++)
   {
      val = transRed_ptr[i];

      if(MINUSZERO < val && val < PLUSZERO)
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
   cout << "\n\n";
}

void theData::info_basisIdx_ptr()
{
   int i;

   cout << "<< basisIdx_ptr >> \n";

   for(i = 0; i < row; i++) cout << basisIdx_ptr[i] << " ";
   cout << "\n\n";
}

void theData::info_nf_pos_ptr()
{
   int i;

   cout << "<< nf_pos_ptr >> \n";

   for(i = 0; i < nfN; i++) cout << nf_pos_ptr[i] << " ";
   cout << "\n\n";
}

void theData::info_nbIdx_ptr()
{
   int i;

   cout << "<< nbIdx_ptr >> \n";

   for(i = 0; i < col; i++) cout << nbIdx_ptr[i] << " ";
   cout << "\n\n";
}

void theData::info_redVec_ptr()
{
   int i;

   cout << "<< redVec_ptr >> \n";

   for(i = 0; i < col; i++)
   {
      // cout << i << " : " << redVec[i] << " | ";
      cout << redVec_ptr[i] << " ";
   }
   cout << "\n\n";
}

void theData::info_fIdx()
{
   cout << fIdx + 1 << "\n";
}

void theData::info_node()
{
   int i;

   cout << "( ";
   for(i = 0; i < polyDim + 1; i++)
   {
      cout << nodeLabel[i] + 1 << " ";
   }
   cout << ") ";
}

ftData::ftData()
{
   elemNum = 0;

   cur = NULL;
   parent = NULL;

   limit = NULL;

   head = NULL;
   last = NULL;

   Dim = 0;
}

ftData::~ftData()
{

}

void ftData::create_elem ( int row, int col, int termS, int polyDim )
{
   theData *newData;

   newData = new theData;
  
   (*newData).create(row, col, termS, polyDim);

   cur = newData;

   Dim = row;
}

void ftData::add_elem()
{
   if(last != NULL)
      (*last).next = cur;
   else
   {
      head = cur;
      parent = cur;
   }
   last = cur;
   elemNum++;
}

void ftData::mark()
{
   limit = last;
}

void ftData::clear()
{
   theData *curr;

   curr = head;

   while(curr != NULL)
   {
      curr->clear();
      curr = curr->next;
   }
}

void ftData::clear_transMat()
{
   theData *curr;

   curr = head;

   while(curr != NULL)
   {
      curr->clear_transMat();
      curr = curr->next;
   }
}

void ftData::delete_cur()
{
   delete cur;
}

void ftData::delete_all()
{
   theData *curr, *tmp;

   curr = head;

  /*
   while(curr != NULL)
   {
      tmp = curr->next;

      if(curr->flag == ICHECK)
      {
         curr->invB = NULL;  
         curr->p_sol = NULL; 
         curr->d_sol = NULL; 
         curr->basisIdx = NULL; 
         curr->nbIdx = NULL; 
         curr->nf_pos = NULL;

         delete curr;
      }
      else if(curr->flag == MCHECK)
      {
         curr->invB = NULL;  
         curr->p_sol = NULL; 
         curr->d_sol = NULL; 
         curr->redVec = NULL;
         curr->basisIdx = NULL; 
         curr->nbIdx = NULL; 
         curr->rIdx = NULL; 
   
         delete curr;
      }
      else
         delete curr;

      curr = tmp;
    
   }
  */

   while(curr != NULL)
   {
      tmp = curr->next;

      delete curr;
      curr = tmp;
   }
   cur = NULL;
   parent = NULL;

   head = NULL;
   last = NULL;

   elemNum = 0;
}

void ftData::delete_addedElem()
{
   int num = 0;
   theData *curr, *tmp;

   // cout << "<< delete_addedElem >>\n\n";

   curr = limit->next;

   while(curr != NULL)
   {
      // cout << "----- " << num + 1 << " -----\n";

      tmp = curr->next;
      delete curr;
      curr = tmp;
      num++;
   }
   limit->next = NULL;
   cur = (last = limit);
   // cout << "\n";
}

void ftData::make_init_data ( int termSumNum, int supN, 
                              int termS, int reTermS )
{
   int i;

   for(i = 0; i < Dim; i++)
   {
      (*cur).nf_pos[i] = i;
      (*cur).invB[i * (Dim + 1)] = 1;
      (*cur).transMat[i * ( Dim + 1)] = 1;
      (*cur).basisIdx[i] = termSumNum - supN + i;
      (*cur).d_sol[i] = 1; 
   }
   for(i = 0; i < termS - 1; i++)
   {
      (*cur).nbIdx[i] = reTermS + i; 
      (*cur).rIdx[i] = -1 * (i + 1); 
   }
}

void ftData::decrease_nfN()
{
   (*cur).nfN--;
}

void ftData::next_data()
{
   if(parent != NULL) parent = parent->next;
}

void ftData::output ( int repIdx, int* idx2, int* nbN, int* nfN )
{
   *idx2 = (*parent).rIdx[repIdx - 1]; 
    
   *nbN = (*parent).nbN;
   *nfN = (*parent).nfN;
}

void ftData::create_rIdx ( int preNbN, int repIdx, int* candIdx )
{
   int i, cnt;
   int idx, candNum, tmp_val;
  
   // cout << "<< create_rIdx >> \n\n";

   tmp_val =  preNbN - Dim + 1;
   candNum = candIdx[0];

   cnt = 0;

   // cout << "candNum: " << candNum << "\n\n";

   for(i = 0; i < candNum; i++)
   {
      idx = candIdx[i + 1];

      if(repIdx > idx)
      {
         (*cur).rIdx[idx] = -1 * (tmp_val + cnt);
         cnt++;
      }
      else if(repIdx < idx)
      {
         (*cur).rIdx[idx - 1] = -1 * (tmp_val + cnt);
         cnt++;
      }
   }
}

void ftData::init_info()
{
   (*cur).pivOutNum = 0;

   memset((*cur).pivOutCheck, 0, sizeof(int) * Dim);
   memset((*cur).transRed, 0, sizeof(double) * Dim);
}

void ftData::get_nbIdx_rIdx ( int preNbN, int repIdx, int* candIdx, 
                              int reTermS, theData* pre_data )
{
   int i, cnt;
   int idx, sub_idx, candNum, pre_length;

   // cout << "<< get_nbIdx_rIdx >> \n\n";

   pre_length = preNbN - Dim;

   memcpy((*cur).nbIdx, (*pre_data).nbIdx_ptr, sizeof(int) * pre_length);

   candNum = candIdx[0];
   cnt = 0;

   for(i = 0; i < candNum; i++)
   {
      idx = candIdx[i + 1];

      if(repIdx > idx)
      {
         sub_idx = pre_length + cnt;
	
         (*cur).nbIdx[sub_idx] = idx + reTermS;
         (*cur).rIdx[idx] = -1 * (sub_idx + 1);

         cnt++;
      }
      else if(repIdx < idx)
      {
         sub_idx = pre_length + cnt;
	
         (*cur).nbIdx[sub_idx] = idx + reTermS - 1;
         (*cur).rIdx[idx - 1] = -1 * (sub_idx + 1);

         cnt++;
      }
   }
}

void ftData::iCopy ( int preNbN, int nfN, int repIdx, int termS,
                     int reTermS, int* candIdx, theData* pre_data )
{
   int i, j;
   int ii, nfPos, pre_length;
   int candNum, idx, sub_idx, cnt;

   // cout << "<< iCopy >> \n\n";
  
   memcpy((*cur).nf_pos, (*pre_data).nf_pos_ptr, sizeof(int) * nfN); 

   pre_length = preNbN - Dim;
   candNum = candIdx[0];

   cnt = 0;

   for(i = 0; i < candNum; i++)
   {
      idx = candIdx[i + 1];

      if(repIdx > idx)
      {
         sub_idx = pre_length + cnt;
	
         (*cur).nbIdx[sub_idx] = idx + reTermS;
         (*cur).rIdx[idx] = -1 * (sub_idx + 1);

         cnt++;
      }
      else if(repIdx < idx)
      {
         sub_idx = pre_length + cnt;
	
         (*cur).nbIdx[sub_idx] = idx + reTermS - 1;
         (*cur).rIdx[idx - 1] = -1 * (sub_idx + 1);
   
         cnt++;
      }
   }
   for(j = 0; j < nfN; j++)
   {
      nfPos = (*pre_data).nf_pos_ptr[j];
      ii = nfPos * Dim;

      for(i = 0; i < Dim; i++)
      {
         (*cur).invB[ii + i] = (*pre_data).invB_ptr[ii + i];
      }
   }
}

void ftData::iGetPtr(theData* pre_data)
{
   cur->invB_ptr = pre_data->invB_ptr;

   cur->p_sol_ptr = pre_data->p_sol_ptr;
   cur->d_sol_ptr = pre_data->d_sol_ptr;

   cur->basisIdx_ptr = pre_data->basisIdx_ptr;
   cur->nf_pos_ptr = pre_data->nf_pos_ptr;
}

void ftData::copy_rIdx(theData* pre_data, int termS)
{
   (*cur).artV = (*pre_data).artV;

   memcpy((*cur).rIdx, (*pre_data).rIdx, sizeof(int) * termS);
}

void ftData::copy_pivOutIdx(theData* pre_data)
{
   int i; 
   int pivOutNum, idx;

   memset((*cur).pivOutCheck, 0, sizeof(int) * Dim);

   (*cur).pivOutNum = (pivOutNum = (*pre_data).pivOutNum);

   for(i = 0; i < pivOutNum; i++)
   {
      idx = (*pre_data).pivOutList[i];

      (*cur).pivOutCheck[idx] = 1;
      (*cur).pivOutList[i] = idx;
   }
}

void ftData::get_nf_pos ( theData* pre_data, int nfN, int idx2 )
{
   int i;
   int nfPos, cnt;

   cnt = 0;

   for(i = 0; i < nfN; i++)
   {
      nfPos = (*pre_data).nf_pos_ptr[i];

      if((*pre_data).nf_pos_ptr[i] != idx2)
      {
         (*cur).nf_pos[cnt] = (*pre_data).nf_pos_ptr[i];
         cnt++;
      }
   }
   (*cur).nfN--;
}

void ftData::mCopy
 ( int nbN, int nfN,  int idx2, int termS, theData* pre_data )
{
   int i, j;
   int ii, nfPos, cnt;

   // cout << "<< mCopy >> \n\n";

   (*cur).artV = (*pre_data).artV;

   memcpy((*cur).redVec, (*pre_data).redVec, sizeof(double) * nbN); 
   memcpy((*cur).rIdx, (*pre_data).rIdx, sizeof(int) * termS);

   memcpy((*cur).transRed, (*pre_data).transRed, sizeof(double) * Dim);

   cnt = 0;
   for(j = 0; j < nfN; j++)
   {
      nfPos = (*pre_data).nf_pos[j];
      ii = nfPos * Dim;

      for(i = 0; i < Dim; i++)
      {
         (*cur).invB[ii + i] = (*pre_data).invB[ii + i];
         (*cur).transMat[ii + i] = (*pre_data).transMat[ii + i];
      }
      if((*pre_data).nf_pos[j] != idx2)
      {
         (*cur).nf_pos[cnt] = (*pre_data).nf_pos[j];
         cnt++;
      }
   }
   (*cur).nfN--;
}

void ftData::mGetPtr ( theData* pre_data )
{
   // cur->artV = pre_data->artV;

   cur->invB_ptr = pre_data->invB_ptr;

   cur->transMat_ptr = pre_data->transMat_ptr;
   cur->transRed_ptr = pre_data->transRed_ptr;

   cur->p_sol_ptr = pre_data->p_sol_ptr;
   cur->d_sol_ptr = pre_data->d_sol_ptr;

   cur->redVec_ptr = pre_data->redVec_ptr;

   cur->basisIdx_ptr = pre_data->basisIdx_ptr;
   cur->nbIdx_ptr = pre_data->nbIdx_ptr;
}

void ftData::copy ( int col, theData* pre_data )
{
   // cout << "<< n_copy >> \n\n";

   memcpy((*cur).p_sol, (*pre_data).p_sol_ptr, sizeof(double) * col);
   memcpy((*cur).d_sol, (*pre_data).d_sol_ptr, sizeof(double) * Dim);
   memcpy((*cur).basisIdx, (*pre_data).basisIdx_ptr, sizeof(int) * Dim);
   memcpy((*cur).nbIdx, (*pre_data).nbIdx_ptr, sizeof(int) * (col - Dim));
}

void ftData::get_ptr(theData* pre_data)
{
   // cout << "<< n_copy >> \n\n";

   delete [] (*cur).p_sol;
   delete [] (*cur).d_sol;
   delete [] (*cur).basisIdx;
   delete [] (*cur).nbIdx;

   (*cur).p_sol = (*pre_data).p_sol;
   (*cur).d_sol = (*pre_data).d_sol;
   (*cur).basisIdx = (*pre_data).basisIdx;
   (*cur).nbIdx = (*pre_data).nbIdx;
}

void ftData::put_sup(int* sup)
{
   // cout << "<< put_sup >>\n";

   if(parent != NULL)
   {
      (*sup) = (*parent).nodeLabel[0];
   }
}

void ftData::info_parent_nbN_nfN()
{
   cout << "nbN: " << (*parent).nbN << "\n";
   cout << "nfN: " << (*parent).nfN << "\n\n";
}

///

void ftData::info_parent_p_sol()
{
   (*parent).info_p_sol();
}

void ftData::info_parent_d_sol()
{
   (*parent).info_d_sol();
}

void ftData::info_parent_invB()
{
   (*parent).info_invB();
}

void ftData::info_parent_transMat()
{
   (*parent).info_transMat();
}

void ftData::info_parent_transRed()
{
   (*parent).info_transRed();
}

void ftData::info_parent_basisIdx()
{
   (*parent).info_basisIdx();
}

void ftData::info_parent_nf_pos()
{
   (*parent).info_nf_pos();
}

void ftData::info_parent_nbIdx()
{
   (*parent).info_nbIdx();
}

void ftData::info_parent_redVec()
{
   (*parent).info_redVec();
}

void ftData::info_parent_rIdx()
{
   (*parent).info_rIdx();
}

void ftData::info_parent_pivOutIdx()
{
   (*parent).info_pivOutIdx();
}

///

void ftData::info_parent_p_sol_ptr()
{
   (*parent).info_p_sol_ptr();
}

void ftData::info_parent_d_sol_ptr()
{
   (*parent).info_d_sol_ptr();
}

void ftData::info_parent_invB_ptr()
{
   (*parent).info_invB_ptr();
}

void ftData::info_parent_transMat_ptr()
{
   (*parent).info_transMat_ptr();
}

void ftData::info_parent_transRed_ptr()
{
   (*parent).info_transRed_ptr();
}

void ftData::info_parent_basisIdx_ptr()
{
   (*parent).info_basisIdx_ptr();
}

void ftData::info_parent_nf_pos_ptr()
{
   (*parent).info_nf_pos_ptr();
}

void ftData::info_parent_nbIdx_ptr()
{
   (*parent).info_nbIdx_ptr();
}

void ftData::info_parent_redVec_ptr()
{
   (*parent).info_redVec_ptr();
}

void ftData::info_parent()
{
   info_parent_p_sol();
   info_parent_d_sol();
   info_parent_invB();
   info_parent_basisIdx();
   info_parent_nf_pos();
   info_parent_nbIdx();
   info_parent_redVec();
   info_parent_rIdx();
}

void ftData::info_parent_ptr()
{
   info_parent_p_sol_ptr();
   info_parent_d_sol_ptr();
   info_parent_invB_ptr();
   info_parent_basisIdx_ptr();
   info_parent_nf_pos_ptr();
   info_parent_nbIdx_ptr();
   info_parent_redVec_ptr();
}

void ftData::info_parent_node()
{
   (*parent).info_node();
}

void ftData::info_cur_nbN_nfN()
{
   cout << "nbN: " << (*cur).nbN << "\n";
   cout << "nfN: " << (*cur).nfN << "\n\n";
}

void ftData::info_cur_p_sol()
{
   (*cur).info_p_sol();
}

void ftData::info_cur_d_sol()
{
   (*cur).info_d_sol();
}

void ftData::info_cur_invB()
{
   (*cur).info_invB();
}

void ftData::info_cur_transMat()
{
   (*cur).info_transMat();
}

void ftData::info_cur_transRed()
{
   (*cur).info_transRed();
}

void ftData::info_cur_basisIdx()
{
   (*cur).info_basisIdx();
}

void ftData::info_cur_nf_pos()
{
   (*cur).info_nf_pos();
}

void ftData::info_cur_nbIdx()
{
   (*cur).info_nbIdx();
}

void ftData::info_cur_redVec()
{
   (*cur).info_redVec();
}

void ftData::info_cur_rIdx()
{
   (*cur).info_rIdx();
}

void ftData::info_cur_pivOutIdx()
{
   (*cur).info_pivOutIdx();
}

void ftData::info_cur_p_sol_ptr()
{
   (*cur).info_p_sol_ptr();
}

void ftData::info_cur_d_sol_ptr()
{
   (*cur).info_d_sol_ptr();
}

void ftData::info_cur_invB_ptr()
{
   (*cur).info_invB_ptr();
}

void ftData::info_cur_transMat_ptr()
{
   (*cur).info_transMat_ptr();
}

void ftData::info_cur_transRed_ptr()
{
   (*cur).info_transRed_ptr();
}

void ftData::info_cur_basisIdx_ptr()
{
   (*cur).info_basisIdx_ptr();

}

void ftData::info_cur_nf_pos_ptr()
{
   (*cur).info_nf_pos_ptr();
}

void ftData::info_cur_nbIdx_ptr()
{
   (*cur).info_nbIdx_ptr();
}

void ftData::info_cur_redVec_ptr()
{
   (*cur).info_redVec_ptr();
}

void ftData::info_cur()
{
   info_cur_p_sol();
   info_cur_d_sol();
   info_cur_invB();
   info_cur_transMat();
   info_cur_basisIdx();
   info_cur_nf_pos();
   info_cur_nbIdx();
   info_cur_redVec();
   info_cur_rIdx();
}

void ftData::info_cur_ptr()
{
   info_cur_p_sol_ptr();
   info_cur_d_sol_ptr();
   info_cur_invB_ptr();
   info_cur_transMat_ptr();
   info_cur_basisIdx_ptr();
   info_cur_nf_pos_ptr();
   info_cur_nbIdx_ptr();
   info_cur_redVec_ptr();
}

void ftData::info_cur_node()
{
   (*cur).info_node();
}

void ftData::info_all_node()
{
   int i;

   theData *curr;

   curr = head;

   cout << "<< info_all_node >>\n\n";

   i = 0;

   while(curr != NULL)
   {
      cout << "# " << i + 1 << ": ";
      (*curr).info_node();
      cout << "\n";

      curr = curr->next;
      i++;
   }
   cout << "\n\n";
}

void ftData::info_all_cur()
{
   int num;
   theData *curr;

   num = 0;
   curr = head;

   cout << "<< info_all_cur >>\n\n";

   while(curr != NULL)
   {
      cout << "----- " << num + 1 << " -----\n";
    
      (*curr).info_p_sol();
      (*curr).info_d_sol();
      (*curr).info_invB();
      (*curr).info_basisIdx();
      (*curr).info_nf_pos();
      (*curr).info_nbIdx();
      (*curr).info_redVec();
      (*curr).info_rIdx();

      curr = curr->next;

      num++;
   }
}

void ftData::info_all_nodeNum()
{
   int i;

   theData *curr;

   curr = head;

   cout << "<< info_all_nodeNum >>\n\n";

   i = 0;

   while(curr != NULL)
   {
      cout << "# " << i + 1 << "\n";

      curr = curr->next;
      i++;
   }
   cout << "\n\n";
}

void ftData::info_numElem()
{
   int num;
   theData *curr;

   curr = head;
   num = 0;

   while(curr != NULL)
   {
      curr = curr->next;
      num++;
   }
   cout << num << " ";
}

lvData::lvData()
{
   fTest = NULL;
   Node = NULL;

   mRepN = NULL;
   mFeaIdx = NULL;
   mFea = NULL;

   Dim = 0;

   length = 0;
   termMax = 0;
}

lvData::~lvData()
{
   int i;

   for(i = 0; i < length; i++) fTest[i].delete_all();
   delete [] fTest;

   delete [] mRepN;
   delete [] mFea;
  
   if(mFeaIdx)
   {
      for(i = 0; i < length; i++) delete [] mFeaIdx[i];

      delete [] mFeaIdx;
      mFeaIdx = NULL;
   }
}

void lvData::create ( int depth, int supN, int Dim, 
                      int ori_length, int ori_termMax )
{
   int i;

   length = ori_length;
   termMax = ori_termMax;
  
   fTest = new ftData [length];
   assert(fTest);

   Node = &fTest[length - 1];
  
   mRepN = new int [length];
   assert(mRepN);
   memset(mRepN, 0, length * sizeof(int));
  
   mFeaIdx = new int* [length];
   assert(mFeaIdx);

   for(i = 0; i < length; i++)
   {
      mFeaIdx[i] = new int [termMax];
      assert(mFeaIdx[i]);
   }
   mFea = new int [length];
   assert(mFea);
   memset(mFea, 0, length * sizeof(int));
}

void lvData::get_info ( int** g_mRepN, int*** g_mFeaIdx, int** g_mFea )
{
   (*g_mRepN) = mRepN;
   (*g_mFeaIdx) = mFeaIdx;
   (*g_mFea) = mFea;
}

void lvData::init_ptr()
{
   for(int i = 0; i < length; i++)
   {
      fTest[i].init_ptr();
   }
}

void lvData::info_mFea()
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
