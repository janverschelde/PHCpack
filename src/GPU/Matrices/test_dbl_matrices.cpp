/* Tests operations on series matrices in double precision */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <time.h>
#include "dbl_matrices_testers.h"

using namespace std;

int main ( void )
{
   srand(time(NULL));

   cout << "testing the real upper solver ..." << endl;
   test_real_upper_solver();
   cout << "testing the complex upper solver ..." << endl;
   test_cmplx_upper_solver();

   cout << "testing the real LU factorization ..." << endl;
   test_real_lufac();
   cout << "testing the complex LU factorization ..." << endl;
   test_cmplx_lufac();

   cout << "testing the real LU solver ..." << endl;
   test_real_lu_solver();
   cout << "testing the complex LU solver ..." << endl;
   test_cmplx_lu_solver();

   return 0;
}
