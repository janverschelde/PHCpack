/* Tests polynomial evaluation and differentiation
 * in octo double precision. */

#include <iostream>
#include "prompt_for_setup.h"
#include "dbl8_polynomials_testers.h"

using namespace std;

int main ( void )
{
   int seed,dim,nva,nbr,pwr,deg,vrb;

   prompt_for_setup(&seed,&dim,&nbr,&nva,&pwr,&deg,&vrb);

   int fail = main_dbl8_test_polynomial(seed,dim,nbr,nva,pwr,deg,vrb);

   if(fail == 0)
      cout << "All tests passed." << endl;
   else
      cout << "Number of failed tests : " << fail << endl;

   return 0;
}
