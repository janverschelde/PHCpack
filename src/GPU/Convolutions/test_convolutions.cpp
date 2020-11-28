/* Tests the product of two series in double, double double, triple double
 * quad double, penta double, octo double, and deca double precision. */

#include <iostream>
#include "dbl_convolutions_testers.h"
#include "dbl2_convolutions_testers.h"
#include "dbl3_convolutions_testers.h"
#include "dbl4_convolutions_testers.h"
#include "dbl5_convolutions_testers.h"
#include "dbl8_convolutions_testers.h"
#include "dbl10_convolutions_testers.h"

using namespace std;

int main ( void )
{
   cout << "Give the seed (0 for time) : ";
   int seed; cin >> seed;
  
   cout << "Give a degree larger than one : ";
   int deg; cin >> deg;

   cout << "Give the verbose level : ";
   int vrb; cin >> vrb;

   int fail = main_dbl_test(seed,deg,vrb);

   if(fail == 0)
      cout << "All tests in double precision passed." << endl;
   else
      cout << "Number of failed tests in double precison : "
           << fail << endl;

   int fail2 = main_dbl2_test(seed,deg,vrb);

   if(fail2 == 0)
      cout << "All tests in double double precision passed." << endl;
   else
      cout << "Number of failed tests in double double precision : "
           << fail2 << endl;

   int fail3 = main_dbl3_test(seed,deg,vrb);

   if(fail3 == 0)
      cout << "All tests in triple double precision passed." << endl;
   else
      cout << "Number of failed tests in triple double precision : "
           << fail3 << endl;

   int fail4 = main_dbl4_test(seed,deg,vrb);

   if(fail4 == 0)
      cout << "All tests in quad double precision passed." << endl;
   else
      cout << "Number of failed tests in quad double precision : "
           << fail4 << endl;

   int fail5 = main_dbl5_test(seed,deg,vrb);

   if(fail5 == 0)
      cout << "All tests in penta double precision passed." << endl;
   else
      cout << "Number of failed tests in penta double precision : "
           << fail5 << endl;

   int fail8 = main_dbl8_test(seed,deg,vrb);

   if(fail8 == 0)
      cout << "All tests in octo double precision passed." << endl;
   else
      cout << "Number of failed tests in octo double precision : "
           << fail8 << endl;

   int fail10 = main_dbl10_test(seed,deg,vrb);

   if(fail10 == 0)
      cout << "All tests in deca double precision passed." << endl;
   else
      cout << "Number of failed tests in deca double precision : "
           << fail10 << endl;

   return 0;
}
