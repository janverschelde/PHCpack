/* Tests operations on tiled accelerated back substitution
 * in double precision */

#include <iostream>
#include <cstdlib>
#include <time.h>
#include "dbl_tabs_testers.h"

using namespace std;

int main ( void )
{
   srand(time(NULL));

   cout << "testing the real upper inverse ..." << endl;
   test_real_upper_inverse();

   return 0;
}
