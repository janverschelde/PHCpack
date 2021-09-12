/* Tests operations on matrix factorizations in octo double precision */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <time.h>
#include "dbl8_factors_testers.h"

using namespace std;

int main ( void )
{
   srand(time(NULL));

   cout << "testing a real LU factorization ..." << endl;
   test_factors_real8_lufac();
   
   cout << endl;

   cout << "testing a complex LU factorization ..." << endl;
   test_factors_cmplx8_lufac();

   cout << endl;

   cout << "testing a real QR decomposition ..." << endl;
   test_factors_real8_houseqr();

   cout << endl;

   cout << "testing a complex QR decomposition ..." << endl;
   test_factors_cmplx8_houseqr();

   return 0;
}
