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

#include "inputData.h"

dataSet::dataSet()
{
   Dim = 0;
   termSumNum = 0;
   termMax = 0;
   typeMax = 0;
   supN = 0;

   support = NULL;

   termSet = NULL;
   termStart = NULL;

   type = NULL;

   coef = NULL;

   outFile = NULL;
}

dataSet::~dataSet()
{
   delete [] outFile;
   delete [] support;
   delete [] termStart;
   delete [] termSet;
   delete [] type;
   delete [] coef;
}

void dataSet::getInputFile ( char* inputFile )
{
   size_t length;
   int i, counter = 0, rowCounter, colCounter, ccounter, sum;
   char *tmpName, ch, str[STRLENGTH], *p;

   length = strlen(inputFile);

   for(i = 0; i < (signed) length; i++)
   {
      if(inputFile[i] == '.') break;
      counter++;
   }
   tmpName = new char [counter + EXCESS];
   assert(tmpName);

   outFile = new char [length + EXCESS];
   assert(outFile);

   for(i = 0; i < counter; i++)
   {
      tmpName[i] = inputFile[i];
      outFile[i] = inputFile[i];
   }
   outFile[counter] = '\0';
   strcat(outFile, ".out");

   delete [] tmpName;

   ifstream in(inputFile);

   if(!in)
   {
      cout << "Cannot open -- " << inputFile << "\n";
      exit(EXIT_FAILURE);
   }
   ccounter = 0;
   colCounter = 0;

   do
   {
      p = str;
      ch = in.peek();

      if(ch == '#')
      {
         do
         {
            ch = in.get();
         }
         while((ch != '\n') && (!in.eof()));
      }
      else if(ch == 'D')
      {
         do
         {
            ch = in.get();
         }
         while((ch != '=') && (ch != '\n') && (!in.eof()));

         do
         {
	    ch = in.peek();

            if(isdigit(ch))
            {
               while(isdigit(*p = in.get())) p++;

               in.putback(*p);
               *p = '\0';
               Dim = atoi(str);
            }
            ch = in.get();
         }
         while((ch != '\n') && (!in.eof()));
      }
      else if(ch == 'S')
      {
         do
         {
            ch = in.get();
         }
         while((ch != '=') && (ch != '\n') && (!in.eof()));

         do
         {
            ch = in.peek();
	
            if(isdigit(ch))
            {
               while(isdigit(*p = in.get())) p++;

               in.putback(*p);
               *p = '\0';

               supN = atoi(str);

               termSet = new int [supN];
               assert(termSet);

               termStart = new int [supN + 1];
               assert(termStart);

               type = new int [supN];
               assert(type);
	    }
            ch = in.get();
         }
         while((ch != '\n') && (!in.eof()));
      }
      else if(ch == 'E')
      {
         do
         {
            ch = in.get();
         }
         while((ch != '=') && (ch != '\n') && (!in.eof()));

         termSumNum = 0;
         counter = 0;

         do
         {
	    ch = in.peek();
	
            if(isdigit(ch))
            {
               p = str;
               while(isdigit(*p = in.get())) p++;
	  
               in.putback(*p);
               *p = '\0';
	  
               termSet[counter] = atoi(str);

               if(termMax < termSet[counter]) termMax = termSet[counter];
               counter++;
	    }
            ch = in.get();
         }
         while((ch != '\n') && (!in.eof()));

         for(i = 0; i < supN; i++) termSumNum += termSet[i];

         termStart[0] = (sum = 0);
         for(i = 1; i < supN + 1; i++)
         {
            sum += termSet[i - 1];
            termStart[i] = sum;
         }
         support = new double [termSumNum * Dim];
         assert(support);
         memset(support, 0, termSumNum * Dim * sizeof(double));

         coef = new double [2 * termSumNum];
         assert(coef);
      }
      else if(ch == 'T')
      {
         do
         {
            ch = in.get();
         }
         while((ch != '=') && (ch != '\n') && (!in.eof()));

         counter = 0;
         do
         {
            ch = in.peek();
	
            if(isdigit(ch))
            {
               p = str;
               while(isdigit(*p = in.get())) p++;
	  
               in.putback(*p);
               *p = '\0';
	  
               type[counter] = atoi(str);

               if(typeMax < type[counter])
               {
                  typeMax = type[counter];
               }
               counter++;
            }
            ch = in.get();
         }
         while((ch != '\n') && (!in.eof()));

         cout << "\n";
      }
      else if(isdigit(ch))
      {
         rowCounter = 0;
         colCounter++;

         do
         {
            ch = in.peek();
	
            if(isdigit(ch))
            {
               p = str;

               while(isdigit(*p = in.get())) p++;
               in.putback(*p);
               *p = '\0';

               support_in(colCounter - 1, rowCounter, atoi(str));

               rowCounter++;				   
            }
            ch = in.get();
         }
         while((ch != '\n') && (!in.eof()));
      }
      else
         in.get();
   }
   while(!in.eof());

   in.close();

#if DBG_INFO

   info_preamble();
   info_supports();

#endif

}

void dataSet::info_preamble()
{
   int i;
  
   cout << "Dim = " << Dim << "\n";
   cout << "Support = " << supN << "\n\n";
  
   cout << "Elem = ";
   for(i = 0; i < supN; i++) cout << termSet[i] << " ";
   cout << "\n";
    
   cout << "Type = ";
   for(i = 0; i < supN; i++) cout << type[i] << " ";

   cout << "\n\n";
}

void dataSet::info_supports()
{
   int i, j, k, top, counter;

   counter = (top = 0);

   for(k = 0; k < supN; k++)
   {
      for(j = top; j < termSet[k] + top; j++)
      {
         for(i = 0; i < Dim; i++)
         {
	    cout << support_out(j, i) << " ";
         }
         cout << "\n";
         counter++;
      }
      cout << "\n";
      top = counter;
   }
}
