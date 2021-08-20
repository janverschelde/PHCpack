/* Tests operations on blocked accelerated QR decomposition
 * and tiled accelerated back substitution in double precision */

#include <iostream>
#include <cstdlib>
#include <time.h>
#include "prompt_baqr_setup.h"
#include "dbl_qrbs_testers.h"

using namespace std;

int main ( void )
{
   int seed,szt,nbt,nrows,vrb,mode;

   prompt_baqr_setup(&seed,&szt,&nbt,&nrows,&vrb,&mode);

   if(seed == 0)
      srand(time(NULL));
   else
      srand(seed);

   cout << "Testing the real blocked QR and back substitution ..." << endl;
   test_real_blocked_qrbs(seed,szt,nbt,nrows,vrb,mode);

   cout << endl;

   prompt_baqr_setup(&seed,&szt,&nbt,&nrows,&vrb,&mode);

   if(seed == 0)
      srand(time(NULL));
   else
      srand(seed);

   cout << "Testing the complex blocked QR and back substitution ..." << endl;
   test_cmplx_blocked_qrbs(seed,szt,nbt,nrows,vrb,mode);

   return 0;
}
