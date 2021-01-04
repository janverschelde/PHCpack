/* Tests polynomial evaluation and differentiation in double,
 * double double, triple double, quad double, penta double,
 * octo double, and deca double precision. */

#include <iostream>
#include "random_polynomials.h"
#include "dbl_polynomials_testers.h"
#include "dbl2_polynomials_testers.h"
#include "dbl3_polynomials_testers.h"
#include "dbl4_polynomials_testers.h"
#include "dbl5_polynomials_testers.h"
#include "dbl8_polynomials_testers.h"
#include "dbl10_polynomials_testers.h"

using namespace std;

int main ( void )
{
   cout << "Give the seed (0 for time) : ";
   int seed; cin >> seed;

   cout << "Give the dimension : ";
   int dim;  cin >> dim;

   cout << "Give the variables per monomial (0 for random polynomial) : ";
   int nva; cin >> nva;

   int nbr; // number of monomials, not counting the constant

   if(nva > 0)
   {
      cout << "Enter 0 for products, other number of cyclic : ";
      cin >> nbr;

      if(nbr == 0)
         nbr = products_count(dim,nva);
      else
         nbr = dim;

      cout << "-> number of monomials : " << nbr << endl;
   }
   else
   {
      cout << "Give the number of terms : ";
      cin >> nbr;
   }
   // cout << "Give the largest power of each variable : "; cin >> pwr;
   const int pwr=1;

   cout << "Give the degree of the series : ";
   int deg; cin >> deg;

   cout << "Give the verbose level : ";
   int vrb; cin >> vrb;

   int fail = main_dbl_test_polynomial(seed,dim,nbr,nva,pwr,deg,vrb);

   if(fail == 0)
      cout << "All tests in double precision passed." << endl;
   else
      cout << "Number of failed tests in double precision : "
           << fail << endl;

   int fail2 = main_dbl2_test_polynomial(seed,dim,nbr,nva,pwr,deg,vrb);

   if(fail2 == 0)
      cout << "All tests in double double precision passed." << endl;
   else
      cout << "Number of failed tests in double double precision : "
           << fail << endl;

   int fail3 = main_dbl3_test_polynomial(seed,dim,nbr,nva,pwr,deg,vrb);

   if(fail3 == 0)
      cout << "All tests in triple double precision passed." << endl;
   else
      cout << "Number of failed tests in triple double precision : "
           << fail << endl;

   int fail4 = main_dbl4_test_polynomial(seed,dim,nbr,nva,pwr,deg,vrb);

   if(fail4 == 0)
      cout << "All tests in quad double precision passed." << endl;
   else
      cout << "Number of failed tests in quad double precision : "
           << fail << endl;

   int fail5 = main_dbl5_test_polynomial(seed,dim,nbr,nva,pwr,deg,vrb);

   if(fail5 == 0)
      cout << "All tests in penta double precision passed." << endl;
   else
      cout << "Number of failed tests in penta double precision : "
           << fail << endl;

   int fail8 = main_dbl8_test_polynomial(seed,dim,nbr,nva,pwr,deg,vrb);

   if(fail8 == 0)
      cout << "All tests in octo double precision passed." << endl;
   else
      cout << "Number of failed tests in octo double precision : "
           << fail << endl;

   int fail10 = main_dbl10_test_polynomial(seed,dim,nbr,nva,pwr,deg,vrb);

   if(fail10 == 0)
      cout << "All tests in deca double precision passed." << endl;
   else
      cout << "Number of failed tests in deca double precision : "
           << fail << endl;

   int sumfail = fail + fail2 + fail3 + fail4 + fail5 + fail8 + fail10;

   if(sumfail == 0)
      cout << "All tests in all precisions passed." << endl;
   else
      cout << "Total number of failed tests : " << sumfail << endl;

   return 0;
}
