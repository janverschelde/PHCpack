/* Tests operations on blocked accelerated QR decomposition
 * in double double precision */

#include <iostream>
#include <cstdlib>
#include <time.h>
#include "dbl2_baqr_testers.h"

using namespace std;

int main ( void )
{
   srand(time(NULL));

   cout << "Testing the real blocked QR ..." << endl;
   test_real2_blocked_qr();

   cout << endl;

   cout << "Testing the complex blocked QR ..." << endl;
   test_cmplx2_blocked_qr();

   return 0;
}
