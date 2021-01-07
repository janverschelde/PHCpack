/* Runs polynomial evaluation and differentiation in deca double precision
 * for some specific cases to demonstrate the performance. */

#include <iostream>
#include "random_polynomials.h"
#include "dbl10_polynomials_testers.h"

using namespace std;

int main ( void )
{
   int seed,dim,nva,nbr,pwr,deg,vrb;

   cout << "Give the seed (0 for time) : "; cin >> seed;

   dim = 16; nva = 4; nbr = products_count(dim,nva); pwr = 1; deg = 31;

   vrb = 2;

   int fail = main_dbl10_test_polynomial(seed,dim,nbr,nva,pwr,deg,vrb);

   if(fail == 0)
      cout << "All tests passed." << endl;
   else
      cout << "Number of failed tests : " << fail << endl;

   return 0;
}
