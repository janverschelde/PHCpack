/*
 The codes in DEMiCs have been written 
 by Tomohiko Mizutani (mizutan8@is.titech.ac.jp)

 Version 0.95, 7 August, 2007
*/

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

#include "global.h"
#include "inputData.h"
#include "mvc.h"

void getInputPara
 ( int argc, char** argv, int& seedNum, int& output, char** inputFileName );
/*
 * DESCRIPTION :
 *   Reads the input parameters from the command line.
 *
 * ON ENTRY :
 *   argc    the count of command line arguments;
 *   argv    an array of strings holding the names of the arguments.
 *
 * ON RETURN :
 *   seedNum is the seed for the random number generator,
 *           if 1, then the default seed is used;
 *   output  is a flag to indicate if output is wanted or not,
 *           if 0, then no extra output is written to screen;
 *   inputFileName is the name of the input file. */

int main ( int argc, char* argv[] )
{
   int seedNum, output;
   char* inputFileName;
  
   dataSet Data;
   mvc MV_Comp;

   double cpuTime_start, cpuTime_end;
   // struct tms cpuTime;

   getInputPara(argc,argv,seedNum,output,&inputFileName);

   //times(&cpuTime);
   //cpuTime_start = (double) cpuTime.tms_utime;

   Data.getInputFile(inputFileName);

   MV_Comp.allocateAndIni(Data,seedNum,output);

   MV_Comp.Enum();

   //times(&cpuTime);
   //cpuTime_end = (double) cpuTime.tms_utime;
 
   MV_Comp.info_cpuTime(cpuTime_start,cpuTime_end);
    
   return (EXIT_SUCCESS);
}

void getInputPara
 ( int argc, char** argv, int& seedNum, int& output, char** inputFileName )
{
   int length, flag = 0;
  
   if(argc == 2)
   {
      seedNum = 1;
      output = 0;

      *inputFileName = argv[1];
   }
   if(argc > 2 && argc < 5 && argv[1][0] == '-')
   {
      length = strlen(argv[1]);

      if(argc == 3 && length == 2 && 
         argv[1][0] == '-' && argv[1][1] == 'c')
      {
         seedNum = 1;
         output = 1;

         *inputFileName = argv[2];
      }
      else if(argc == 4 && length == 2 && 
              argv[1][0] == '-' && argv[1][1] == 's')
      {
         seedNum = atoi(argv[3]);
         output = 0;

         *inputFileName = argv[2];
      }
      else if(argc == 4 && length == 3 && 
              argv[1][0] == '-' && argv[1][1] == 'c' && argv[1][2] == 's')
      {
         seedNum = atoi(argv[3]);
         output = 1;

         *inputFileName = argv[2];
      }
      else
         flag = 1;
   }
   else
   {
      if(argc == 2)
      {
         seedNum = 1;
         output = 0;

         *inputFileName = argv[1];
      }
      else
         flag = 1;
   }
   if(flag)
   {
      cerr <<"-----" << endl;
      cerr <<"Usage" << endl;
      cerr <<"-----" << endl << endl;
      cerr << argv[0] << " [option] input_file [seed_number]"
           << endl << endl << endl;

      cerr << "------" << endl;
      cerr << "Option" << endl;
      cerr << "------" << endl << endl;

      cerr << "default : " <<  argv[0] << " input_file" << endl << endl;
      cerr << "-c : Output information about mixed cells in a terminal"
           << endl;
      cerr << "    " << argv[0] << " -c input_file" << endl << endl;
      cerr << "-s : Designate a seed number" << endl ;
      cerr << "    " << argv[0] << " -s input_file seed_number"
           << endl << endl;
      cerr << "-cs : Output information about mixed cells in a terminal, ";
      cerr << "and set a seed number" << endl ;
      cerr << "    " << argv[0] << " -cs input_file seed_number" << endl;

      exit(EXIT_FAILURE);
   }
}
