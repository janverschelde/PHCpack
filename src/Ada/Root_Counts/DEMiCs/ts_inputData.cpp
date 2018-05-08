// test on getting the input data for DEMiCs

#include <string>
#include <iostream>
#include <cstring>
#include "inputData.h"

using namespace std;

void read_data_from_file ( void );
/*
 * DESCRIPTION :
 *   Prompts the user for the name of the data input file for demics
 *   and then prints the data read from the file. */

void interactive_input_data ( void );
/*
 * DESCRIPTION :
 *   Prompts the user for input data, dimension, number of distinct supports,
 *   the number of occurrences of each supports, the number of points in
 *   each support, and then prompt for each point in the supports. */

int main ( int argc, char* argv[] )
{
   cout << "Reading data from file ? (y/n) ";
   char ans; cin >> ans;

   if(ans == 'y')
      read_data_from_file();
   else
      interactive_input_data();

   return 0;
}

void read_data_from_file ( void )
{
   cout << "Give the name of the input file : ";
   string name; cin >> name;

   const char* inputFileName = name.c_str();

   char cstrname[80];

   strcpy(cstrname,inputFileName);

   printf("The name of the file : %s\n", cstrname);
  
   dataSet Data;
 
   Data.getInputFile(cstrname);

   Data.info_preamble();

   Data.info_supports();
}

void interactive_input_data ( void )
{
   dataSet Data;

   cout << "Give the dimension : ";
   cin >> Data.Dim;

   cout << "Give the number of distinct supports : ";
   cin >> Data.supN;

   Data.type = new int[Data.supN];

   for(int k=0; k<Data.supN; k++)
   {
      cout << "Give the number of occurrences of support " << k+1 << " : ";
      cin >> Data.type[k];
   }

   Data.termSet = new int[Data.supN];

   for(int k=0; k<Data.supN; k++)
   {
      cout << "Give the number of points in support " << k+1 << " : ";
      cin >> Data.termSet[k];
   }

   cout << endl;
   cout << "The dimension, the number of distinct support sets," << endl;
   cout << "the number of points in each support set, and" << endl;
   cout << "the number of occurrences of each support set :" << endl;
   cout << endl;

   Data.info_preamble();

   Data.termSumNum = 0;
   for(int k=0; k<Data.supN; k++)
      Data.termSumNum = Data.termSumNum + Data.termSet[k];

   Data.support = new double[Data.termSumNum*Data.Dim];

   int offset = 0;

   for(int i=0; i<Data.supN; i++)
   {
      cout << "Reading the points of support " << i+1 << " ..." << endl;
      for(int j=0; j<Data.termSet[i]; j++)
      { 
         cout << "Give coordinates of point " << j+1 << " : ";
         for(int k=0; k<Data.Dim; k++)
         {
            double x; cin >> x;
            Data.support_in(offset+j,k,x);
         }
      }
      offset = offset + Data.termSet[i];
   }

   cout << endl;
   cout << "The points in the support sets : " << endl << endl;

   Data.info_supports();
}
