/* Tests operations on series vectors in double precision */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <time.h>
#include "dbl_matrices_testers.h"

using namespace std;

int main ( void )
{
   srand(time(NULL));

   cout << "testing the real LU solver ..." << endl;
   test_real_lu_solver();
   cout << "testing the real LU factorization ..." << endl;
   test_real_lufac();
   cout << "testing a real inner product ..." << endl;
   test_real_inner_product();
   cout << "testing a complex inner product ..." << endl;
   test_cmplx_inner_product();
   cout << "testing a real matrix-vector product ..." << endl;
   test_real_matrix_vector_product();
   cout << "testing a complex matrix-vector product ..." << endl;
   test_cmplx_matrix_vector_product();
   cout << "testing the real upper solver ..." << endl;
   test_real_upper_solver();
   cout << "testing the complex upper solver ..." << endl;
   test_cmplx_upper_solver();

   return 0;
}
