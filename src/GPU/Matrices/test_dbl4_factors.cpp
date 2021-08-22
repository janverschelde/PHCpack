/* Tests operations on matrix factorizations in quad double precision */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <time.h>
#include "dbl4_factors_testers.h"

using namespace std;

int main ( void )
{
   srand(time(NULL));

   cout << "testing a real LU factorization ..." << endl;
   test_factors_real4_lufac();
   
   cout << endl;

   cout << "testing a complex LU factorization ..." << endl;
   test_factors_cmplx4_lufac();

   cout << endl;

   cout << "testing a real QR decomposition ..." << endl;
   test_factors_real4_houseqr();

   cout << endl;

   cout << "testing a complex QR decomposition ..." << endl;
   test_factors_cmplx4_houseqr();

   return 0;
}
