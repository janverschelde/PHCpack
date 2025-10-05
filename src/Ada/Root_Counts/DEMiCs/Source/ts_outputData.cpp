/* test on the operations in outputData */

#include <iostream>
#include "outputData.h"

using namespace std;

void print_lifting ( int nbrsup, int* crdsup );
/*
 * DESCRIPTION :
 *   Given in nbrsup the number of different supports
 *   and in crdsup the cardinality of each supports,
 *   prints the lifting value for all points. */

void test_lifting ( void );
/*
 * DESCRIPTION :
 *   Prompts the user for the number of supports,
 *   cardinalities of the supports, and then lifting values. */

void test_strings ( void );
/*
 * DESCRIPTION :
 *   Test the storage and retrieval of strings. */

void test_mixed_volume ( void );
/*
 * DESCRIPTION :
 *   Tests storing and retrieving of mixed volume. */

int main ( int argc, char* argv[] )
{
   adainit();

   test_mixed_volume();
   test_strings();
   test_lifting();

   adafinal();

   return 0;
}

void print_lifting ( int nbrsup, int* crdsup )
{
   for(int i=0; i<nbrsup; i++)
   {
      cout << "lifting values for support " << i << " : " << endl;
      for(int j=0; j<crdsup[i]; j++)
      {
         double val;
         int fail = demics_retrieve_lifting(i,j,&val);
         cout << " " << val;
      }
      cout << endl;
   }
}

void test_lifting ( void )
{
   cout << "Give the number of distinct supports : ";
   int nbrsup; cin >> nbrsup;

   int* crdsup = new int[nbrsup];

   for(int k=0; k<nbrsup; k++)
   {
      cout << "Give the number of points in support " << k+1 << " : ";
      cin >> crdsup[k];
   }
   int fail = demics_allocate_lifting(nbrsup,crdsup);

   do
   {
      cout << "Give an index to a support : "; 
      int idxsup; cin >> idxsup;

      cout << "Give an index to a point : ";
      int idxpnt; cin >> idxpnt;

      cout << "Give a lifting value : ";
      double val; cin >> val;

      fail = demics_assign_lifting(idxsup,idxpnt,val);

      print_lifting(nbrsup,crdsup);

      cout << endl;
      cout << "Continue ? (y/n) ";
      char ch; cin >> ch;
 
      if(ch != 'y') break;
   }
   while(true);

   fail = demics_clear_lifting();
}

void test_strings ( void )
{
   int fail;

   do
   {
      cout << "Give a string : ";
      string strcell; cin >> strcell;

      fail = demics_append_cell_indices(strcell);

      do
      {
         cout << "Give an index (0 to quit) : ";
         int idx; cin >> idx;

         if(idx <= 0) break;

         char cell[80];
         fail = demics_retrieve_cell_indices(idx,cell);
         cout << "The retrieved string : " << cell << endl;
      }
      while(true);

      cout << "Continue ? (y/n) ";
      char ch; cin >> ch;
 
      if(ch != 'y') break;
   }
   while(true);
}

void test_mixed_volume ( void )
{
   int mv;

   cout << "Give the mixed volume : "; cin >> mv;

   int fail = demics_store_mixed_volume(mv);

   int retmv;

   fail = demics_retrieve_mixed_volume(&retmv);

   cout << "The retrieved mixed volume : " << retmv << endl;
}
