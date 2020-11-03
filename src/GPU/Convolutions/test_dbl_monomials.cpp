/* Tests the evaluation and differentiation of a monomial
 * in double precision. */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include "random_numbers.h"
#include "dbl_monomials_host.h"

using namespace std;

bool sorted_insert ( int n, int *data );
/*
 * DESCRIPTION :
 *   Inserts data[n] in the sequence of n sorted numbers in data.
 *   Returns true if data[n] was already inserted, that is:
 *   there is an index k less than n, for which data[k] == data[n].
 *   Returns false if data[n] is not a duplicate number. */

bool make_real_monomial
 ( int dim, int nvr, int pwr, int deg, int *idx, int *exp, double *cff );
/*
 * DESCRIPTION :
 *   Makes a monomial in several variables, with a power series coefficient,
 *   generating random exponents and coefficients.
 *   Writes an error message and returns true if nvr > dim.
 *
 * ON ENTRY :
 *   dim     dimension, total number of variables;
 *   nvr     number of variables with positive power in the monomial;
 *   pwr     largest power of a variable;
 *   deg     degree of the power series coefficient;
 *   exp     space allocated for nvr integers;
 *   cff     space allocated for deg doubles.
 *
 * ON RETURN :
 *   idx     nvr integers in the range from 0 to dim-1,
 *           idx(k) is the index of the k-th variable in the monomial;
 *   exp     nvr positive integers with the powers of the variables,
 *           exp(k) is the power of the variable with index idx(k);
 *   cff     deg doubles with the coefficients of the power series. */

void common_factors ( int nvr, int *exp, int *nbrfac, int *expfac );
/*
 * DESCRIPTION :
 *   Extracts all exponents strictly larger than one.
 *  
 * ON ENTRY :
 *   nvr     number of variables in exp, exp[k] >= 1,
 *           for all k from 0 to nvr-1;
 *   exp     exponents of a monomial;
 *   expfac  space for nvr integers.
 *
 * ON RETURN :
 *   nbrfac  number of exponents in exp strictly larger than one;
 *   expfac  exponents of the common factor,
 *           if exp[k] > 1, then expfac[k] = exp[k]-1.  */

void make_real_input ( int dim, int deg, double **data );
/*
 * DESCRIPTION :
 *   Generates input series, as many as dim, of degree deg.
 *
 * ON ENTRY :
 *   dim      dimension of the input;
 *   deg      degree of the power series;
 *   data     space allocated for dim series of degree deg.
 *
 * ON RETURN :
 *   data     data[i][j] is the j-th coefficient of the i-th series,
 *            for i in 0..dim-1 and j in 0..deg. */

int main ( void )
{
   int dim,nvr,pwr,deg;

   cout << "Give the dimension : "; cin >> dim;
   cout << "Give the number of variables, <= "; cout << dim;
   cout << " : "; cin >> nvr;
   cout << "Give the largest power of each variable : "; cin >> pwr;
   cout << "Give the degree of the series : "; cin >> deg;

   int *idx = new int[nvr];         // indices of variables in the monomial
   int *exp = new int[nvr];         // exponents of the variables
   int *expfac = new int[nvr];      // exponents of common factor
   int nbrfac;                      // number of common factors
   double *cff = new double[deg+1]; // series coefficient

   // The input are dim power series of degree deg,
   // the output are nvr+1 power series of degree deg,
   // for the evaluated and differentiated monomial.

   double **input = new double*[dim];
   for(int i=0; i<dim; i++) input[i] = new double[deg];
   double **output = new double*[dim+1];
   for(int i=0; i<=dim; i++) output[i] = new double[deg];

   srand(time(NULL));

   bool fail = make_real_monomial(dim,nvr,pwr,deg,idx,exp,cff);

   if(!fail)
   {
      cout << "Generated a random monomial :" << endl;
      cout << "   the indices :";
      for(int i=0; i<nvr; i++) cout << " " << idx[i];
      cout << endl;

      cout << " the exponents :";
      for(int i=0; i<nvr; i++) cout << " " << exp[i];
      cout << endl;

      common_factors(nvr,exp,&nbrfac,expfac);

      cout << "common factors :";
      for(int i=0; i<nvr; i++) cout << " " << expfac[i];
      cout << endl;
      cout << "number of common factors : " << nbrfac << endl;

      cout << scientific << setprecision(16);
      cout << "the coefficients :" << endl;
      for(int i=0; i<=deg; i++) cout << " " << cff[i] << endl;;
   }

   make_real_input(dim,deg,input);

   cout << "Random input series :" << endl;
   for(int i=0; i<dim; i++)
   {
      cout << "-> coefficients of series " << i << " :" << endl;
      for(int j=0; j<=deg; j++) cout << input[i][j] << endl;
   }

   CPU_dbl_evaldiff(dim,nvr,deg,idx,input,output);

   cout << "The value of the product :" << endl;
   for(int i=0; i<=deg; i++)
      cout << output[dim][i] << endl;

   for(int k=0; k<nvr; k++)
   {
      cout << "-> derivative for index " << idx[k] << " :" << endl;
      for(int i=0; i<=deg; i++)
         cout << output[idx[k]][i] << endl;
   }

   return 0;
}

bool sorted_insert ( int n, int *data )
{
   if(n == 0)
      return false;
   else
   {
      int nbr = data[n];
      int idx = n;

      for(int i=0; i<n; i++)
         if(data[i] >= nbr)
         {
            idx = i; break;
         }

      if(idx == n)
         return false;                  // sequence is already sorted
      else
      {
         if(data[idx] == nbr)           // found duplicate number
            return true;
         else
         {
            for(int i=n; i>idx; i--)
               data[i] = data[i-1];     // shift the numbers

            data[idx] = nbr;            // insert number
            return false;
         }
      }
   }
}

bool make_real_monomial
 ( int dim, int nvr, int pwr, int deg, int *idx, int *exp, double *cff )
{
   bool fail;

   if(nvr > dim)
   {
      cout << "ERROR: nvr = " << nvr << " > " << dim << " dim" << endl;
      return true;
   }
   else
   {
      for(int i=0; i<=deg; i++) cff[i] = random_double();

      for(int i=0; i<nvr; i++)
      {
         exp[i] = 1 + (rand() % pwr);
         do
         {
            idx[i] = rand() % dim;
            fail = sorted_insert(i,idx);
         }
         while(fail);
      }
      return false;
   }
}

void common_factors ( int nvr, int *exp, int *nbrfac, int *expfac )
{
   *nbrfac = 0;

   for(int i=0; i<nvr; i++)
   {
      if(exp[i] <= 1)
         expfac[i] = 0;
      else
      {
         expfac[i] = exp[i] - 1;
         *nbrfac = *nbrfac + 1;
      }
   }
}

void make_real_input ( int dim, int deg, double **data )
{
   for(int i=0; i<dim; i++)
      for(int j=0; j<=deg; j++) data[i][j] = random_double();
}
