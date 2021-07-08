/* Tests operations on matrix factorizations in double double precision */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <time.h>
#include "dbl2_factors_testers.h"

using namespace std;

int main ( void )
{
   srand(time(NULL));

   cout << "testing a complex lu factorization ..." << endl;
   test_factors_cmplx2_lufac();

   cout << "testing a real lu factorization ..." << endl;
   test_factors_real2_lufac();

   return 0;
}
