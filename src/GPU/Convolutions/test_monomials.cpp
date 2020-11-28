/* Tests monomial evaluation and differentiation in double, double double,
 * and triple double precision. */

#include <iostream>
#include "dbl_monomials_testers.h"
#include "dbl2_monomials_testers.h"
#include "dbl3_monomials_testers.h"

using namespace std;

int main ( void )
{
   cout << "Give the seed (0 for time) : ";
   int seed; cin >> seed;

   cout << "Give the dimension : ";
   int dim;  cin >> dim;

   cout << "Give the number of variables, <= "; cout << dim;
   cout << " : "; int nvr; cin >> nvr;

   // cout << "Give the largest power of each variable : "; cin >> pwr;
   const int pwr=1;

   cout << "Give the degree of the series : ";
   int deg; cin >> deg;

   cout << "Give the verbose level : ";
   int vrb; cin >> vrb;

   int fail = main_dbl_test(seed,dim,nvr,pwr,deg,vrb);

   if(fail == 0)
      cout << "All tests in double precision passed." << endl;
   else
      cout << "Number of failed tests in double precison : "
           << fail << endl;

   int fail2 = main_dbl2_test(seed,dim,nvr,pwr,deg,vrb);

   if(fail2 == 0)
      cout << "All tests in double double precision passed." << endl;
   else
      cout << "Number of failed tests in double double precision : "
           << fail2 << endl;

   int fail3 = main_dbl3_test(seed,dim,nvr,pwr,deg,vrb);

   if(fail3 == 0)
      cout << "All tests in triple double precision passed." << endl;
   else
      cout << "Number of failed tests in triple double precision : "
           << fail3 << endl;

   return 0;
}
