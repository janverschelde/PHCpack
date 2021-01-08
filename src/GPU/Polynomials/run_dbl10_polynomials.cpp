/* Runs polynomial evaluation and differentiation in deca double precision
 * for some specific cases to demonstrate the performance. */

#include <iostream>
#include "random_polynomials.h"
#include "dbl10_polynomials_testers.h"

using namespace std;

int main ( void )
{
   int seed,dim,nva,nbr,pwr,deg,vrb,fail;

   cout << "Give the seed (0 for time) : "; cin >> seed;

   dim = 16; nva = 4; nbr = products_count(dim,nva); pwr = 1; vrb = 2;

   deg = 15;
   fail = main_dbl10_test_polynomial
             (seed,dim,nbr,nva,pwr,deg,vrb,1.0e-152,true);
   deg = 31;
   cout << "---> running for degree 31 ..." << endl;
   fail += main_dbl10_test_polynomial
              (seed,dim,nbr,nva,pwr,deg,vrb,1.0e-152,false);
   deg = 63;
   cout << "---> running for degree 63 ..." << endl;
   fail += main_dbl10_test_polynomial
              (seed,dim,nbr,nva,pwr,deg,vrb,1.0e-152,false);
   deg = 95;
   cout << "---> running for degree 95 ..." << endl;
   fail += main_dbl10_test_polynomial
              (seed,dim,nbr,nva,pwr,deg,vrb,1.0e-152,false);
   deg = 127;
   cout << "---> running for degree 127 ..." << endl;
   fail += main_dbl10_test_polynomial
              (seed,dim,nbr,nva,pwr,deg,vrb,1.0e-152,false);
   deg = 152;
   cout << "---> running for degree 152 ..." << endl;
   fail += main_dbl10_test_polynomial
              (seed,dim,nbr,nva,pwr,deg,vrb,1.0e-152,false);

   if(fail == 0)
      cout << "All tests passed." << endl;
   else
      cout << "Number of failed tests : " << fail << endl;

   return 0;
}
