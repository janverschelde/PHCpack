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

int main ( int argc, char* argv[] )
{
   adainit();

   cout << "Give the number of distinct supports : ";
   int nbrsup; cin >> nbrsup;

   int* crdsup = new int[nbrsup];

   for(int k=0; k<nbrsup; k++)
   {
      cout << "Give the number of points in support " << k+1 << " : ";
      cin >> crdsup[k];
   }
   int fail = allocate_lifting(nbrsup,crdsup);

   do
   {
      cout << "Give an index to a support : "; 
      int idxsup; cin >> idxsup;

      cout << "Give an index to a point : ";
      int idxpnt; cin >> idxpnt;

      cout << "Give a lifting value : ";
      double val; cin >> val;

      fail = assign_lifting(idxsup,idxpnt,val);

      print_lifting(nbrsup,crdsup);

      cout << endl;
      cout << "Continue ? (y/n) ";
      char ch; cin >> ch;
 
      if(ch != 'y') break;
   }
   while(true);

   fail = clear_lifting();

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
         int fail = retrieve_lifting(i,j,&val);
         cout << " " << val;
      }
      cout << endl;
   }
}
