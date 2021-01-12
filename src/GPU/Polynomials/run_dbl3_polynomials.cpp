/* Runs polynomial evaluation and differentiation in triple double precision
 * for some specific cases to demonstrate the performance. */

#include <iostream>
#include "random_polynomials.h"
#include "dbl3_polynomials_testers.h"

using namespace std;

int main ( void )
{
   int seed,dim,nva,nbr,pwr,deg,vrb,fail,mode;

   cout << "Give the seed (0 for time) : "; cin >> seed;
   cout << "Enter 0 (GPU only), 1 (CPU only), or 2 (GPU+CPU) : ";
   cin >> mode;

   dim = 16; nva = 4; nbr = products_count(dim,nva); pwr = 1; vrb = 2;

   fail = test_dbl3_sequence(seed,dim,nva,nbr,pwr,vrb,true,mode);

   dim = 128; nva = 64; nbr = 128; pwr = 1; vrb = 2;

   fail += test_dbl3_sequence(seed,dim,nva,nbr,pwr,vrb,true,mode);

   dim = 128; nva = 2; nbr = products_count(dim,nva); pwr = 1; vrb = 2;

   fail += test_dbl3_sequence(seed,dim,nva,nbr,pwr,vrb,true,mode);

   if(mode == 2)
   {
      if(fail == 0)
         cout << "All tests passed." << endl;
      else
         cout << "Number of failed tests : " << fail << endl;
   }
   return 0;
}
