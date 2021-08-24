/* Tests the operations on tiled accelerated back substitution
 * in quad double precision */

#include <iostream>
#include <cstdlib>
#include <time.h>
#include "dbl4_tabs_testers.h"

using namespace std;

int main ( void )
{
   srand(time(NULL));

   cout << "Testing the real upper inverse of one tile ..." << endl;
   test_real4_upper_inverse();

   cout << "Testing the complex upper inverse of one tile ..." << endl;
   test_cmplx4_upper_inverse();

   cout << endl
        << "Testing tiling on a real upper triangular matrix ..."
        << endl;
   test_real4_upper_tiling();

   cout << endl
        << "Testing tiling on a complex upper triangular matrix ..."
        << endl;

   test_cmplx4_upper_tiling();

   return 0;
}
