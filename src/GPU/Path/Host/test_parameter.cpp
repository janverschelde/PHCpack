// Test on the interactive tuning of parameters.

#include "parameter.h"

using namespace std;

int main ( void )
{
   int precision;
  
   cout << "Give the precision (16, 32, or 64) : ";
   cin >> precision;

   Parameter pars(precision);

   pars.tune();

   return 0;
}
