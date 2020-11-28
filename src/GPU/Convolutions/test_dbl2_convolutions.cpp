/* Tests the product of two series in double double precision. */

#include <iostream>
#include "dbl2_convolutions_testers.h"

using namespace std;

int main ( void )
{
   cout << "Give the seed (0 for time) : ";
   int seed; cin >> seed;
  
   cout << "Give a degree larger than one : ";
   int deg; cin >> deg;

   cout << "Give the verbose level : ";
   int vrb; cin >> vrb;

   int fail = main_dbl2_test(seed,deg,vrb);

   if(fail == 0)
      cout << "All tests passed." << endl;
   else
      cout << "Number of failed tests : " << fail << endl;

   return 0;
}
