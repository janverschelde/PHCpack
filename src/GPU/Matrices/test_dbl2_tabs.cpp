/* Tests the operations on tiled accelerated back substitution
 * in double double precision */

#include <iostream>
#include <cstdlib>
#include <time.h>
#include "dbl2_tabs_testers.h"

using namespace std;

int main ( void )
{
   srand(time(NULL));

   cout << "Testing the real upper inverse ..." << endl;
   test_real2_upper_inverse();

   cout << endl
        << "Testing tiling on a real upper triangular matrix ..."
        << endl;
   test_real2_upper_tiling();

   return 0;
}