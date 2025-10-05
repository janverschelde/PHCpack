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

#ifndef __ITEST_H
#define __ITEST_H

#include "global.h"
#include "inputData.h"

class uData
{
   private:

      int nfN;
  
   public:

      uData();
      ~uData();

      uData *next;
      uData *prev;

      uData *fNext;
  
      int supLab;

      double red;
      double* dir;

      void create ( int depth, int Dim );
      void init();

      void getDir ( double& val, int idx ){ dir[idx] = val; };
      void getRed ( double& val, int idx ){ red = val; };

      void info_dirRed();
};

class inifData
{
   private:

   public:

      inifData();
      ~inifData();

      uData* head;
      uData* fHead;
      uData* last;
 
      void create ( int length, int depth, int Dim );

      void get_info ( dataSet& Data, double* lifting, int* termSet,
                      int* termStart, int depth, int Dim, int supN );

      void info_all_dirRed();
      void info_feasIdx();
      void info_fNext();
      void info_next();
      void info_prev();
};

class iLvData
{
   private:

      int rspLen;
      int inifLen;

   public:

      iLvData();
      ~iLvData();

      inifData* inif;
      int* rsp;

      void create ( int depth, int supN, int Dim, int termMax );

      void getInit ( dataSet& Data, double* lifting, int* termSet,
                     int* termStart, int Dim, int supN );
  
      void init ( int supN, int depth, int* preRsp );

      void info_rsp();
      void info_all_dirRed();
      void info_feasIdx ( int depth );
      void info_all_feasIdx();
};

#endif
