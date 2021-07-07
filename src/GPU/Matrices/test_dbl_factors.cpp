/* Tests operations on matrix factorizations in double precision */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <time.h>
#include "dbl_factors_testers.h"

using namespace std;

int main ( void )
{
   srand(time(NULL));

   cout << "testing a complex lu factorization ..." << endl;
   test_factors_cmplx_lufac();

   cout << endl;

   cout << "testing a real lu factorization ..." << endl;
   test_factors_real_lufac();

   return 0;
}
