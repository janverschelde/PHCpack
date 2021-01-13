/* Runs polynomial evaluation and differentiation in double,
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

void run
 ( int seed, int dim, int nva, int nbr, int pwr, int vrb, int mode );
/*
 * DESCRIPTION :
 *   For increasing precision and for an increasing sequence of degrees,
 *   runs tests in double precision.
 *
 * ON ENTRY :
 *   seed     seed for the random number generator;
 *   dim      dimension, total number of variables;
 *   nva      number of variables per monomial (for products and cyclic);
 *   nbr      number of terms in the polynomial;
 *   pwr      highest power of each variable;
 *   vrb      verbose level, if zero, then no output is written,
 *            otherwise, the higher the value, the more input;
 *   mode     the mode of execution, either 0, 1, or 2, as follows:
 *            0 : GPU only; 1 : CPU only; 2 : GPU and CPU. */

int main ( void )
{
   int seed,dim,nva,nbr,pwr,vrb,mode;

   cout << "Give the seed (0 for time) : "; cin >> seed;
   cout << "Enter 0 (GPU only), 1 (CPU only), or 2 (GPU+CPU) : ";
   cin >> mode;

   dim = 16; nva = 4; nbr = products_count(dim,nva); pwr = 1; vrb = 2;

   run(seed,dim,nva,nbr,pwr,vrb,mode);

   dim = 128; nva = 64; nbr = 128; pwr = 1; vrb = 2;

   run(seed,dim,nva,nbr,pwr,vrb,mode);

   dim = 128; nva = 2; nbr = products_count(dim,nva); pwr = 1; vrb = 2;

   run(seed,dim,nva,nbr,pwr,vrb,mode);

   return 0;
}

void run
 ( int seed, int dim, int nva, int nbr, int pwr, int vrb, int mode )
{
   cout << endl << "running in double precision ..." << endl;

   int fail = test_dbl_sequence(seed,dim,nva,nbr,pwr,vrb,true,mode);

   if(mode == 2)
   {
      if(fail == 0)
         cout << "All tests in double precision passed." << endl;
      else
         cout << "Number of failed tests in double precision : "
              << fail << endl;
   }
   cout << endl << "running in double double precision ..." << endl;

   int fail2 = test_dbl2_sequence(seed,dim,nva,nbr,pwr,vrb,false,mode);

   if(mode == 2)
   {
      if(fail2 == 0)
         cout << "All tests in double double precision passed." << endl;
      else
         cout << "Number of failed tests in double double precision : "
              << fail2 << endl;
   }
   cout << endl << "running in triple double precision ..." << endl;

   int fail3 = test_dbl3_sequence(seed,dim,nva,nbr,pwr,vrb,false,mode);

   if(mode == 2)
   {
      if(fail3 == 0)
         cout << "All tests in triple double precision passed." << endl;
      else
         cout << "Number of failed tests in triple double precision : "
              << fail3 << endl;
   }
   cout << endl << "running in quad double precision ..." << endl;

   int fail4 = test_dbl4_sequence(seed,dim,nva,nbr,pwr,vrb,false,mode);

   if(mode == 2)
   {
      if(fail4 == 0)
         cout << "All tests in quad double precision passed." << endl;
      else
         cout << "Number of failed tests in quad double precision : "
              << fail4 << endl;
   }
   cout << endl << "running in penta double precision ..." << endl;

   int fail5 = test_dbl5_sequence(seed,dim,nva,nbr,pwr,vrb,false,mode);

   if(mode == 2)
   {
      if(fail5 == 0)
         cout << "All tests in penta double precision passed." << endl;
      else
         cout << "Number of failed tests in penta double precision : "
              << fail5 << endl;
   }
   cout << endl << "running in octo double precision ..." << endl;

   int fail8 = test_dbl8_sequence(seed,dim,nva,nbr,pwr,vrb,false,mode);

   if(mode == 2)
   {
      if(fail8 == 0)
         cout << "All tests in octo double precision passed." << endl;
      else
         cout << "Number of failed tests in octo double precision : "
              << fail8 << endl;
   }
   cout << endl << "running in deca double precision ..." << endl;

   int fail10 = test_dbl10_sequence(seed,dim,nva,nbr,pwr,vrb,false,mode);

   if(mode == 2)
   {
      if(fail10 == 0)
         cout << "All tests in deca double precision passed." << endl;
      else
         cout << "Number of failed tests in deca double precision : "
              << fail10 << endl;

     int sumfail = fail + fail2 + fail3 + fail4 + fail5 + fail8 + fail10;

     if(sumfail == 0)
        cout << "All tests in all precisions passed." << endl;
     else
        cout << "Total number of failed tests : " << sumfail << endl;
   }
}
