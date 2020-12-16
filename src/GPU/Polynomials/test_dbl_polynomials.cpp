/* Tests polynomial evaluation and differentiation in double precision. */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include "random_numbers.h"
#include "random_monomials.h"
#include "dbl_monomials_testers.h"

using namespace std;

double test_dbl_real_polynomial
 ( int dim, int nbr, int pwr, int deg, int verbose );
/*
 * DESCRIPTION :
 *   Tests the evaluation and differentiation for random real data.
 *   Returns the sum of all errors.
 * 
 * ON ENTRY :
 *   dim      dimension, total number of variables;
 *   nbr      number of terms in the polynomial;
 *   pwr      highest power of each variable;
 *   deg      truncation degree of the series;
 *   verbose  if zero, then no output is written. */

int main_dbl_test_polynomial
 ( int seed, int dim, int nbr, int pwr, int deg, int vrblvl );
/*
 * DESCRIPTION :
 *   Runs tests on a random polynomial in double precision.
 *   Returns 0 if all tests passed,
 *   otherwise, returns the number of failed tests.
 *
 * ON ENTRY :
 *   seed     seed for the random number generator;
 *   dim      dimension, total number of variables;
 *   nbr      number of terms in the polynomial;
 *   pwr      highest power of each variable;
 *   deg      truncation degree of the series;
 *   vrblvl   is the verbose level, if 0 then no output. */

int main ( void )
{
   cout << "Give the seed (0 for time) : ";
   int seed; cin >> seed;

   cout << "Give the dimension : ";
   int dim;  cin >> dim;

   cout << "Give the number of terms : ";
   int nbr; cin >> nbr;

   // cout << "Give the largest power of each variable : "; cin >> pwr;
   const int pwr=1;

   cout << "Give the degree of the series : ";
   int deg; cin >> deg;

   cout << "Give the verbose level : ";
   int vrb; cin >> vrb;

   int fail = main_dbl_test_polynomial(seed,dim,nbr,pwr,deg,vrb);

   if(fail == 0)
      cout << "All tests passed." << endl;
   else
      cout << "Number of failed tests : " << fail << endl;

   return 0;
}

int main_dbl_test_polynomial
 ( int seed, int dim, int nbr, int pwr, int deg, int vrblvl )
{
   int seedused;

   if(seed != 0)
   {
      srand(seed);
      seedused = seed;
   }
   else
   {
      const int timevalue = time(NULL); // for a random seed
      srand(timevalue);
      seedused = timevalue;
   }
   double realsum = test_dbl_real_polynomial(dim,nbr,pwr,deg,vrblvl-1);

   const double tol = 1.0e-12;

   int fail = int(realsum > tol);

   if(vrblvl > 0)
   {
      cout << scientific << setprecision(2);
      cout << "Sum of all errors in double precision :" << endl;
      cout << "  on real data : " << realsum;
      if(realsum < tol)
         cout << "  pass," << endl;
      else
         cout << "  fail!" << endl;

      cout << "  Seed used : " <<  seedused << endl;
   }
   return fail;
}

double test_dbl_real_polynomial
 ( int dim, int nbr, int pwr, int deg, int verbose )
{
   if(nbr < 1)
      return 0.0;
   else if(nbr == 1) // one term is the test on one monomial
   {
      const int nvr = 1 + (rand() % (dim-1));
      return test_dbl_real(dim,nvr,pwr,deg,verbose);
   }
   else
   {
      double **input = new double*[dim]; // dim series of degree deg
      for(int i=0; i<dim; i++) input[i] = new double[deg+1];
      // The output are dim+1 power series of degree deg
      // for the evaluated and differentiated polynomial.
      double **output_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output_h[i] = new double[deg+1];
      double **output_d = new double*[dim+1];
      for(int i=0; i<=dim; i++) output_d[i] = new double[deg+1];

      make_real_input(dim,deg,input);

      if(verbose > 0)
      {
         cout << scientific << setprecision(16);
         cout << "Random input series :" << endl;
         for(int i=0; i<dim; i++)
         {
            cout << "-> coefficients of series " << i << " :" << endl;
            for(int j=0; j<=deg; j++) cout << input[i][j] << endl;
         }
      }
      double *cst = new double[deg+1]; // constant coefficient series
      for(int i=0; i<=deg; i++) cst[i] = random_double();

      if(verbose > 0)
      {
         cout << "Coefficient series of the constant term :" << endl;
         for(int j=0; j<=deg; j++) cout << cst[j] << endl;
      }
      double **cff = new double*[nbr-1]; // coefficient series of terms
      for(int i=0; i<nbr-1; i++) cff[i] = new double[deg+1];
      int *nvr = new int[nbr-1]; // number of variables in each monomial
      for(int i=0; i<nbr-1; i++) nvr[i] = 1 + (rand() % (dim-1));
      int **idx = new int*[nbr-1];  // indices of variables in monomials
      for(int i=0; i<nbr-1; i++) idx[i] = new int[nvr[i]];
      int **exp = new int*[nbr-1];  // exponents of the variables
      for(int i=0; i<nbr-1; i++) exp[i] = new int[nvr[i]];

      bool fail = false;

      for(int i=0; i<nbr-1; i++)
      {
         fail = make_real_monomial(dim,nvr[i],pwr,deg,idx[i],exp[i],cff[i]);
         if(fail) return 1.0;

         if(verbose > 0)
         {
            cout << "Generated random monomial " << i << " :" << endl;
            cout << "   the indices :";
            for(int j=0; j<nvr[i]; j++) cout << " " << idx[i][j];
            cout << endl;

            cout << " the exponents :";
            for(int j=0; j<nvr[i]; j++) cout << " " << exp[i][j];
            cout << endl;
            cout << " coefficient series :" << endl;
            for(int j=0; j<=deg; j++) cout << " " << cff[i][j] << endl;
         }
      }
      return 0.0;
   }
}
