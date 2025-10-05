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

#include "iTest.h"

uData::uData()
{
   nfN = 0;
   supLab = 0;
  
   red = 0;

   dir = NULL;

   next = NULL;
   prev = NULL;

   fNext = NULL;
}

uData::~uData()
{
   delete [] dir;
}

void uData::create ( int depth, int Dim )
{
   //nfN = Dim - depth;
   nfN = Dim;

   dir = new double [nfN];
   assert(dir);
   memset(dir, 0, nfN * sizeof(double));
}

void uData::init()
{
   memset(dir, 0, nfN * sizeof(double));

   red = 0;
   supLab = 0;
}

void uData::info_dirRed()
{
   int i;
   double val;

   for(i=0; i < nfN; i++)
   {
      val = dir[i];
      
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
   cout << " : " << red << "\n";
}

inifData::inifData()
{
   head = NULL;
   fHead = NULL;
   last = NULL;
}

inifData::~inifData()
{
   uData *curr, *tmp;

   curr = head;

   while(curr != NULL)
   {
      tmp = curr->next;

      delete curr;
      curr = tmp;
   }
   head = NULL;
   last = NULL;
}

void inifData::create ( int length, int depth, int Dim )
{
   int i;
   uData* newData;

   for(i = 0; i < length; i++)
   {
      newData = new uData;

      newData->create(depth, Dim);

      if(last != NULL)
      {
         (*last).fNext = ((*last).next = newData);
         (*newData).prev = last;
      }
      else
      {
         fHead = (head = newData);
      }
      last = newData;
   }
}

void inifData::get_info ( dataSet& Data, double* lifting,
                          int* termSet, int* termStart,
                          int depth, int Dim, int supN )
{
   int i, j;
   int termSta, termS;

   double val, sum, red;

   uData* curr;

   termSta = termStart[depth];
   termS = termSet[depth];

   curr = fHead;

   for(j = termSta; j < termSta + termS; j++)
   {
      sum = 0;

      for(i = 0; i < Dim; i++)
      {
         val = Data.support_out(j, i);
         // sum += val;

         curr->getDir(val, i);
      }
      // red = lifting[k] - sum;
      red = lifting[j];

      curr->getRed(red, j);
      curr->supLab = j - termSta;

      curr = curr->fNext;
   }
   if(curr != NULL)
   {
      curr->prev->fNext = NULL;
   }
}

void inifData::info_all_dirRed()
{
   uData *curr;

   curr = fHead;

   while(curr != NULL)
   {
      cout.setf(ios::right);
      cout.width(3);

      cout << curr->supLab + 1 << " : ";
      curr->info_dirRed();
    
      curr = curr->fNext;
   }
}

void inifData::info_feasIdx()
{
   uData *curr;

   curr = fHead;
  
   while(curr != NULL)
   {
      cout.setf(ios::right);
      cout.width(3);

      cout << curr->supLab + 1 << " ";
      curr = curr->fNext;
   }
}

void inifData::info_fNext()
{
   uData *curr;

   curr = fHead;

   cout << "<< info_fNext >>\n\n";
  
   while(curr != NULL)
   {
      // cout << curr->a[0] << " ";
      curr = curr->fNext;
   }
   cout << "\n\n";
}

void inifData::info_next()
{
   uData *curr;

   curr = head;

   cout << "<< info_next >>\n\n";

   while(curr != NULL)
   {
      // cout << curr->a[0] << " ";
      curr = curr->next;
   }
   cout << "\n\n";
}

void inifData::info_prev()
{
   uData *curr;

   curr = last;

   cout << "<< info_prev >>\n\n";

   while(curr != NULL)
   {
      // cout << curr->a[0] << " ";
      curr = curr->prev;
   }
   cout << "\n\n";
}

iLvData::iLvData()
{
   rspLen = 0;
   inifLen = 0;

   inif = NULL;
   rsp = NULL;
}

iLvData::~iLvData()
{
   delete [] inif;
   delete [] rsp;
}

void iLvData::create ( int depth, int supN, int Dim, int termMax )
{
   rspLen = supN - depth - 1;
   inifLen = supN;

   inif = new inifData [supN];
   assert(inif);

   rsp = new int [supN];
   assert(rsp);

   for(int i = 0; i < inifLen; i++)
   {
      inif[i].create(termMax, depth, Dim);
   }
}

void iLvData::getInit ( dataSet& Data, double* lifting, int* termSet,
                        int* termStart, int Dim, int supN )
{
   for(int i = 0; i < supN; i++)
   {
      inif[i].get_info(Data, lifting, termSet, termStart, i, Dim, supN);

      if(i < supN - 1) rsp[i] = i + 1;

      // rsp[i] = i;
   }
}

void iLvData::init ( int supN, int depth, int* preRsp )
{
   int i, lvl;
   int length;

   uData *curr;

   // cout << "<< init >> \n\n";

   length = supN - depth - 1;

   for(i = 0; i < length; i++)
   {
      lvl = preRsp[i];

      curr = (inif[lvl].fHead = inif[lvl].head);

      // cout << "----- Support: " << lvl + 1 << " -----\n\n";

      while(curr != NULL)
      {
         // curr->init();

         curr->fNext = curr->next;

         if(curr->next != NULL)
         {
	    curr->next->prev = curr;
         }
         curr = curr->next;
      }
   }
}

void iLvData::info_rsp()
{
   cout << "rsp: ";

   for(int i = 0; i < rspLen; i++)
   {
      cout << rsp[i] << " ";
   }
   cout << "\n\n";
}

void iLvData::info_all_dirRed()
{
   cout << "<< info_all_dir_red >>\n\n";

   for(int i = 0; i < inifLen; i++)
   {
      cout << "--- Support: " << i + 1 << " ---\n";

      inif[i].info_all_dirRed();

      cout << "\n";
   }
}

void iLvData::info_feasIdx ( int depth )
{
   cout << "<< info_feasIdx >>\n";

   inif[depth].info_feasIdx();

   cout << "\n\n";
}

void iLvData::info_all_feasIdx()
{
   cout << "<< info_all_feasIdx >>\n\n";

   for(int i = 0; i < inifLen; i++)
   {
      cout << "--- Support: " << i + 1 << " ---\n";

      inif[i].info_feasIdx();

      cout << "\n\n";
   }
}
