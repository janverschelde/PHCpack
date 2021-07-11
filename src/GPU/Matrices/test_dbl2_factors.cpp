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

   cout << "testing a complex QR decomposition ..." << endl;
   test_factors_cmplx2_houseqr();

   cout << endl;

   cout << "testing a real QR decomposition ..." << endl;
   test_factors_real2_houseqr();

   cout << endl;

   cout << "testing a complex LU factorization ..." << endl;
   test_factors_cmplx2_lufac();

   cout << endl;

   cout << "testing a real LU factorization ..." << endl;
   test_factors_real2_lufac();

   return 0;
}
