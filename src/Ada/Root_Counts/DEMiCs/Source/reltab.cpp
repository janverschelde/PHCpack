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

#include "reltab.h"

reltab::reltab ( void )
{
   Dim = 0;
   supN = 0;
   maxConst = 0;
   termSumNum = 0;

   row = 0;
   col = 0;

   unbLP = 0;
   totalLP = 0;

   nbN = 0;
   nfN = 0;

   termSet = NULL;
   termStart = NULL;
   firIdx = NULL;

   invB = NULL;
  
   p_sol = NULL;
   d_sol = NULL;
  
   basisIdx = NULL;
   nbIdx = NULL;

   nf_pos = NULL;

   negIdx = NULL;
   val = NULL;

   feasIdx_a = NULL;
   feasIdx_b = NULL;

   table = NULL;
}

reltab::~reltab ( void )
{
   delete [] invB;
   delete [] p_sol;
   delete [] d_sol;
   delete [] basisIdx;
   delete [] nbIdx;
   delete [] nf_pos;
 
   delete [] negIdx;
   delete [] val;
 
   delete [] feasIdx_a;
   delete [] feasIdx_b;
 
   delete [] table;
}

void reltab::allocateAndIni
 ( simplex& ori_Simplex, int** ori_firIdx,
   int ori_Dim, int ori_supN, int ori_termSumNum, 
   int* ori_termSet, int* ori_termStart, int* ori_re_termStart )
{
   Dim = ori_Dim;
   supN = ori_supN;
   termSumNum = ori_termSumNum;
   maxConst = termSumNum - supN;

   termSet = ori_termSet;
   termStart = ori_termStart;
   re_termStart = ori_re_termStart;
   firIdx = *ori_firIdx;

   Simplex = &ori_Simplex;
  
   row = ori_Dim;
   col = maxConst + Dim;
  
   invB = new double [row * row];
   assert(invB);
   memset(invB, 0, row * row * sizeof(double));
  
   p_sol = new double [col];
   assert(p_sol);
   memset(p_sol, 0, col * sizeof(double));

   d_sol = new double [row];
   assert(d_sol);
   memset(d_sol, 0, row * sizeof(double));

   basisIdx = new int [row];
   assert(basisIdx);
   memset(basisIdx, 0, row * sizeof(int));

   nbIdx = new int [col];
   assert(nbIdx);
   memset(nbIdx, 0, col * sizeof(int));

   nf_pos = new int [row];
   assert(nf_pos);
   memset(nf_pos, 0, row * sizeof(int));

   val = new double [col];
   assert(val);

   negIdx = new int [row + 1];
   assert(negIdx);

   feasIdx_a = new int [row];
   assert(feasIdx_a);

   feasIdx_b = new int [row];
   assert(feasIdx_b);
  
   table = new int [termSumNum * termSumNum];
   assert(table);
   memset(table, 0, termSumNum * termSumNum * sizeof(int));
}

void reltab::get_init_triData ( int lab, int idx )
{
   int i, j;
   int constNum, negNum, reTermS;

   double elem;

   firIdx[lab] = idx;
   constNum = termSet[lab] - 1;

   reTermS = re_termStart[lab];

   nbN = constNum + Dim;
   nfN = Dim;

   srand(4); // srandom(4);
   for(i = 0; i < constNum; i++)
   {
      val[i] = (double) rand() / (double) RAND_MAX;
      nbIdx[i] = reTermS + i; 
   }
   negNum = 0;

   for(i = 0; i < Dim; i++)
   {
      elem = 0;

      for(j = 0; j < constNum; j++)
      {
         elem += val[j] * Simplex->put_elem_supp(lab, idx, i, j);
      }

      if(elem < MINUSZERO)
      {
         p_sol[maxConst + i] = -1 * elem;
         negIdx[negNum + 1] = i;

         negNum++;

         for(j = 0; j < constNum; j++)
         {
            Simplex->mult_elem_supp(lab, idx, i, j);
         }
      }
      else if(elem > PLUSZERO)
      {
         p_sol[maxConst + i] = elem;
      }
      else
      {
         p_sol[maxConst + i] = 0;
      }
   }
   negIdx[0] = negNum;

   for(i = 0; i < Dim; i++)
   {
     nf_pos[i] = i;
     invB[i * (Dim + 1)] = 1;
     basisIdx[i] = maxConst + i;
     d_sol[i] = 1; 
   }
}

void reltab::get_init_squData
 ( int lab_a, int lab_b, int idx_a, int idx_b, int colPos, int rowPos)
{
   int i, j;
   int constNum, constNum_a, constNum_b, negNum;
   int reTermS_a, reTermS_b;

   double elem;

#if DBG_REL_TABLE
   cout << "========== ( " << "sup: " << lab_a + 1<< " idx: "
        << idx_a + 1 << " )-" 
        << "( " << "sup: " << lab_b + 1 << " idx: " << idx_b + 1
        << " ) ==========\n\n";
#endif

   firIdx[lab_a] = idx_a;
   firIdx[lab_b] = idx_b;

   constNum_a = termSet[lab_a] - 1;
   constNum_b = termSet[lab_b] - 1;

   reTermS_a = re_termStart[lab_a];
   reTermS_b = re_termStart[lab_b];

   constNum = constNum_a + constNum_b;

   nbN = constNum + Dim;
   nfN = Dim;

   srand(4); // srandom(4);
   for(i = 0; i < constNum_a; i++)
   {
      nbIdx[i] = reTermS_a + i; 
      val[i] = (double) rand() / (double) RAND_MAX;
   }
   for(i = 0; i < constNum_b; i++)
   {
      nbIdx[constNum_a + i] = reTermS_b + i; 
      val[constNum_a + i] = (double) rand() / (double) RAND_MAX;
   }
   // cout << "constNum_b: " << constNum_b << "\n";
   // cout << "lab_b: " << lab_b << "\n";
   // cout << "idx_b: " << idx_b << "\n\n";

   negNum = 0;

   for(i = 0; i < Dim; i++)
   {
      elem = 0;

      for(j = 0; j < constNum_a; j++)
      {
         elem += val[j] * Simplex->put_elem_supp(lab_a, idx_a, i, j);
      }
      for(j = 0; j < constNum_b; j++)
      {
         elem += val[j+constNum_a] * Simplex->put_elem_supp(lab_b,idx_b,i,j);
      }
      if(elem < MINUSZERO)
      {
         p_sol[maxConst + i] = -1 * elem;
         negIdx[negNum + 1] = i;

         negNum++;

         for(j = 0; j < constNum_a; j++)
         {
            Simplex->mult_elem_supp(lab_a, idx_a, i, j);
         }
         for(j = 0; j < constNum_b; j++)
         {
            Simplex->mult_elem_supp(lab_b, idx_b, i, j);
         }
      }
      else if(elem > PLUSZERO)
      {
         p_sol[maxConst + i] = elem;
      }
      else
      {
         p_sol[maxConst + i] = 0;
      }
   }
   negIdx[0] = negNum;

   for(i = 0; i < Dim; i++)
   {
      nf_pos[i] = i;
      invB[i * (Dim + 1)] = 1;
      basisIdx[i] = maxConst + i;
      d_sol[i] = 1; 
   }

#if DBG_REL_TABLE

   info_invB();
   info_p_sol();
   info_d_sol();
   info_basisIdx();
   info_nbIdx(); 
   info_nf_pos();

#endif
}

void reltab::init_data()
{
   memset(invB, 0, row * row * sizeof(double));
   memset(p_sol, 0, col * sizeof(double));
}

void reltab::init_tri ( int lab, int idx )
{
   int i, j;
   int constNum;

   constNum = termSet[lab] - 1;

   for(j = 0; j < negIdx[0]; j++)
   {
      for(i = 0; i < constNum; i++)
      {
         Simplex->mult_elem_supp(lab, idx, negIdx[j + 1], i);
      }
   }
}

void reltab::init_squ ( int lab_a, int lab_b, int idx_a, int idx_b )
{
   int i, j;
   int constNum_a, constNum_b;

   constNum_a = termSet[lab_a] - 1;
   constNum_b = termSet[lab_b] - 1;

   for(j = 0; j < negIdx[0]; j++)
   {
      for(i = 0; i < constNum_a; i++)
      {
         Simplex->mult_elem_supp(lab_a, idx_a, negIdx[j + 1], i);
      }
      for(i = 0; i < constNum_b; i++)
      {
         Simplex->mult_elem_supp(lab_b, idx_b, negIdx[j + 1], i);
      }
   }
   init_data();
}

void reltab::put_data()
{
   Simplex->get_nbN_nfN(nbN, nfN);

   Simplex->get_p_sol(p_sol);
   Simplex->get_d_sol(d_sol);

   Simplex->get_basisIdx(basisIdx);
   Simplex->get_nf_pos(nf_pos);
   Simplex->get_nbIdx(nbIdx);

   Simplex->get_invB(invB);
}

void reltab::put_frIdx(int frIdx)
{
   Simplex->get_frIdx(frIdx);
}

void reltab::makeTri()
{
   int i, j, k;
   int len, termSta, reTermS, flag, iter;

   // cout << "========== makeTri ==========\n\n";
  
   for(k = 0; k < supN; k++)
   {
      len = termSet[k];
      reTermS = re_termStart[k];
      termSta = termStart[k];

      for(j = 0; j < len - 1; j++)
      {
         for(i = j + 1; i < len; i++)
         {

#if DBG_REL_TABLE
            cout << "========== ( " << j + 1 << " )-"  << "( " << i + 1
                 << " ) ==========\n\n";
#endif

            if(table_out(termSta + j, termSta + i) != OPT)
            {
               get_init_triData(k, j);

               put_data();
               put_frIdx(reTermS + i - 1);

               iter = 0;
               flag = Simplex->tSolLP(iter, TRIANGLE);

               if(flag == OPT)
               {
                  findAllFeasLPs_tri(k, j, reTermS + i - 1);
               }
               else
               {
                  table_in(termSta + j, termSta + i, UNBOUNDED);
                  table_in(termSta + i, termSta + j, UNBOUNDED);

                  unbLP++;
	       }
               init_data();
               init_tri(k, j);
            }
         }
      }
   }
}

void reltab::makeSqu()
{
   int i, j;
   int supLab_a, supLab_b;
   int len_a, len_b, flag, iter;
   int reTermS, rowPos, colPos;

   // Simplex->info_allSup();

#if DBG_REL_TABLE
   cout << "========== makeSqu ==========\n\n";
#endif

   for(supLab_a = 0; supLab_a < supN; supLab_a++)
   {
      len_a = termSet[supLab_a];
      reTermS = re_termStart[supLab_a];

      rowPos = termStart[supLab_a + 1];
      colPos = termStart[supLab_a];

      for(supLab_b = supLab_a + 1; supLab_b < supN; supLab_b++)
      {
         len_b = termSet[supLab_b];

         for(j = 0; j < len_a; j++)
         {
	    for(i = 0; i < len_b; i++)
            {
	       if(table_out(colPos + j, rowPos + i) != OPT)
               {
	          get_init_squData(supLab_a, supLab_b, j, i, colPos, rowPos);
  
                  put_data();

                  iter = 0;
                  flag = Simplex->tSolLP(iter, SQUARE);  
	    
                  if(flag == OPT)
                  {
                     findAllFeasLPs_squ(supLab_a,supLab_b,j,i,colPos,rowPos);
                  }
                  else
                  {
                     table_in(colPos + j, rowPos + i, UNBOUNDED);
                     table_in(rowPos + i, colPos + j, UNBOUNDED);
	             unbLP++;
                  }
                  init_squ(supLab_a, supLab_b, j, i);
                  init_data();
               }
            }
         }
         rowPos += len_b;
      }
   }
}

void reltab::findAllFeasLPs_tri ( int lab, int idx, int frIdx )
{
   int i, j;
   int reTermS, termSta;
   int bIdx, tmp_idx;
   int fIdx, fIdx_a, fIdx_b;
   int num = 0;

   // cout << "<< findAllFeasLPs_tri >> \n";

   // info_basisIdx();
   // info_nf_pos();

   reTermS = re_termStart[lab];
   termSta = termStart[lab];
  
   for(i = 0; i < Dim; i++)
   {
      bIdx = basisIdx[i];

      if(bIdx < termSumNum - supN)
      {
         tmp_idx = bIdx - reTermS;

         if(tmp_idx >= idx)
         {
            fIdx = tmp_idx + 1; 
         }
         else
         {
            fIdx = tmp_idx;
         }
         table_in(termSta + idx, termSta + fIdx, OPT);
         table_in(termSta + fIdx, termSta + idx, OPT);

         feasIdx_a[num] = fIdx;

         num++;
      }
   }
   for(j = 0; j < num; j++)
   {
      fIdx_b = feasIdx_a[j];

      for(i = j + 1; i < num; i++)
      {
         fIdx_a = feasIdx_a[i];

         table_in(termSta + fIdx_b, termSta + fIdx_a, OPT);
         table_in(termSta + fIdx_a, termSta + fIdx_b, OPT);
      }
   }
   // info_feasIdx_tri(num);
}

void reltab::findAllFeasLPs_squ
 ( int lab_a, int lab_b, int idx_a, int idx_b, int colPos, int rowPos )
{
   int i, j;

   int reTermS_a, reTermS_b, constNum_a, constNum_b;
   int bIdx, tmp_idx;
   int fIdx, fIdx_a, fIdx_b;
   int num_a = 0, num_b = 0;

   // cout << "<< findAllFeasLPs_squ >> \n";
   // info_basisIdx();
   // info_nf_pos();

   constNum_a = termSet[lab_a] - 1;
   constNum_b = termSet[lab_b] - 1;

   reTermS_a = re_termStart[lab_a];
   reTermS_b = re_termStart[lab_b];

   table_in(colPos + idx_a, rowPos + idx_b, OPT);
   table_in(rowPos + idx_b, colPos + idx_a, OPT);

   for(i = 0; i < Dim; i++)
   {
      bIdx = basisIdx[i];

      if(bIdx < termSumNum - supN)
      {
         if(reTermS_a <= bIdx && bIdx < reTermS_a + constNum_a)
         {
	    tmp_idx = bIdx - reTermS_a;

            if(tmp_idx >= idx_a)
            {
               fIdx = tmp_idx + 1; 
            }
            else
            {
               fIdx = tmp_idx;
	    }
	    table_in(colPos + fIdx, rowPos + idx_b, OPT);
            table_in(rowPos + idx_b, colPos + fIdx, OPT);

	    feasIdx_a[num_a] = fIdx;

	    num_a++;
         }
         else
         {
            tmp_idx = bIdx - reTermS_b;

            if(tmp_idx >= idx_b)
            {
               fIdx = tmp_idx + 1; 
            }
            else
            {
               fIdx = tmp_idx;
            }
            table_in(colPos + idx_a, rowPos + fIdx, OPT);
            table_in(rowPos + fIdx, colPos + idx_a, OPT);

            feasIdx_b[num_b] = fIdx;

            num_b++;
         }
      }
   }
   for(j = 0; j < num_b; j++)
   {
      fIdx_b = feasIdx_b[j];

      for(i = 0; i < num_a; i++)
      {
         fIdx_a = feasIdx_a[i];

         table_in(colPos + fIdx_a, rowPos + fIdx_b, OPT);
         table_in(rowPos + fIdx_b, colPos + fIdx_a, OPT);
      }
   }
   // info_feasIdx_squ(num_a, num_b);
}

void reltab::info_invB()
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

void reltab::info_p_sol()
{
   int i;
  
   cout << "<< p_sol >> \n";
   for(i = 0; i < maxConst + Dim; i++) cout << p_sol[i] << " ";
   cout << "\n\n";
}

void reltab::info_d_sol()
{
   int i;
  
   cout << "<< d_sol >> \n";
   for(i = 0; i < Dim; i++) cout << d_sol[i] << " ";
   cout << "\n\n";
}

void reltab::info_basisIdx()
{
   int i;

   cout << "<< basisIdx >> \n";
   for(i = 0; i < Dim; i++) cout << basisIdx[i] << " ";
   cout<< "\n\n";
}

void reltab::info_nf_pos()
{
   int i;

   cout << "<< nf_pos >> \n";
   for(i = 0; i < nfN; i++) cout << nf_pos[i] << " ";
   cout << "\n\n";
}

void reltab::info_nbIdx()
{
   int i;

   cout << "<< nbIdx >> \n";
   for(i = 0; i < nbN - Dim; i++) cout << nbIdx[i] << " ";
   cout << "\n\n";
}

void reltab::info_feasIdx_tri ( int num )
{
   cout << "feasIdx: ";

   for(int i = 0; i < num; i++)
   {
      cout << feasIdx_a[i] << " ";
   }
   cout << "\n\n";
}

void reltab::info_feasIdx_squ ( int num_a, int num_b )
{
   cout << "feasIdx_a: ";

   for(int i = 0; i < num_a; i++)
   {
      cout << feasIdx_a[i] << " ";
   }
   cout << "\n\n";

   cout << "feasIdx_b: ";

   for(int i = 0; i < num_b; i++)
   {
      cout << feasIdx_b[i] << " ";
   }
   cout << "\n\n";
}

void reltab::info_allTable()
{
   int i, j;
   int unbNum = 0;

   cout << "<< All elements on Relation Table >>\n\n";

   for(j = 0; j < termSumNum; j++)
   {
      for(i = 0; i < termSumNum; i++)
      {
         cout << table_out(j, i) << " ";

         if(table_out(j, i) == UNBOUNDED) unbNum++;
      }
      cout << "\n";
   }
   cout << "\n";

   cout << "# Unb. LPs: " << unbNum / 2 << "\n\n";
}

void reltab::info_table()
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

void reltab::makeTable ( double& total_unbLP_tab )
{
   makeTri();
   makeSqu();

   total_unbLP_tab = unbLP;

   // info_table();
}
