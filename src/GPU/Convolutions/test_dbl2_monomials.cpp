/* Tests the evaluation and differentiation of a monomial
 * in double double precision. */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <vector_types.h>
#include "random_numbers.h"
#include "random2_vectors.h"
#include "dbl2_monomials_host.h"
#include "dbl2_monomials_kernels.h"

using namespace std;

bool sorted_insert ( int n, int *data );
/*
 * DESCRIPTION :
 *   Inserts data[n] in the sequence of n sorted numbers in data.
 *   Returns true if data[n] was already inserted, that is:
 *   there is an index k less than n, for which data[k] == data[n].
 *   Returns false if data[n] is not a duplicate number. */

bool make_real_monomial
 ( int dim, int nvr, int pwr, int deg, int *idx, int *exp,
   double *cffhi, double *cfflo );
/*
 * DESCRIPTION :
 *   Makes a monomial in several variables, with a power series coefficient,
 *   generating random exponents and real coefficients.
 *   Writes an error message and returns true if nvr > dim.
 *
 * ON ENTRY :
 *   dim     dimension, total number of variables;
 *   nvr     number of variables with positive power in the monomial;
 *   pwr     largest power of a variable;
 *   deg     degree of the power series coefficient;
 *   exp     space allocated for nvr integers;
 *   cffhi   space allocated for deg+1 doubles;
 *   cfflo   space allocated for deg+1 doubles.
 *
 * ON RETURN :
 *   idx     nvr integers in the range from 0 to dim-1,
 *           idx(k) is the index of the k-th variable in the monomial;
 *   exp     nvr positive integers with the powers of the variables,
 *           exp(k) is the power of the variable with index idx(k);
 *   cffhi   deg+1 doubles with the high doubles of the coefficients;
 *   cfflo   deg+1 doubles with the low doubles of the coefficients. */

bool make_complex_monomial
 ( int dim, int nvr, int pwr, int deg, int *idx, int *exp,
   double *cffrehi, double *cffrelo, double *cffimhi, double *cffimlo );
/*
 * DESCRIPTION :
 *   Makes a monomial in several variables, with a power series coefficient,
 *   generating random exponents and complex coefficients.
 *   Writes an error message and returns true if nvr > dim.
 *
 * ON ENTRY :
 *   dim     dimension, total number of variables;
 *   nvr     number of variables with positive power in the monomial;
 *   pwr     largest power of a variable;
 *   deg     degree of the power series coefficient;
 *   exp     space allocated for nvr integers;
 *   cffrehi has space allocated for deg+1 doubles,
 *           for the high doubles of the real parts of the coefficients;
 *   cffrelo has space allocated for deg+1 doubles,
 *           for the low doubles of the real parts of the coefficients;
 *   cffimhi has space allocated for deg+1 doubles,
 *           for the high doubles of the imaginary parts of the coefficients.
 *   cffimlo has space allocated for deg+1 doubles,
 *           for the low doubles of the imaginary parts of the coefficients.
 *
 * ON RETURN :
 *   idx     nvr integers in the range from 0 to dim-1,
 *           idx(k) is the index of the k-th variable in the monomial;
 *   exp     nvr positive integers with the powers of the variables,
 *           exp(k) is the power of the variable with index idx(k);
 *   cffrehi holds deg+1 doubles with the high doubles of the real parts
 *           of the coefficients of the power series;
 *   cffrelo holds deg+1 doubles with the low doubles of the real parts
 *           of the coefficients of the power series;
 *   cffimhi holds deg+1 doubles with the high doubles of the imaginary parts
 *           of the coefficients of the power series;
 *   cffimhi holds deg+1 doubles with the low doubles of the imaginary parts
 *           of the coefficients of the power series. */

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

void make_real_input ( int dim, int deg, double **datahi, double **datalo );
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
 *   datahi   datahi[i][j] is the high double of the j-th coefficient of 
 *            the i-th series, for i in 0..dim-1 and j in 0..deg;
 *   datalo   datalo[i][j] is the low double of the j-th coefficient of 
 *            the i-th series, for i in 0..dim-1 and j in 0..deg. */

void make_complex_input
 ( int dim, int deg, double **datarehi, double **datarelo,
   double **dataimhi, double **dataimlo );
/*
 * DESCRIPTION :
 *   Generates input series, as many as dim, of degree deg.
 *
 * ON ENTRY :
 *   dim      dimension of the input;
 *   deg      degree of the power series;
 *   datarehi has space allocated for the high doubles of the real parts
 *            of dim series of degree deg.
 *   datarelo has space allocated for the low doubles of the real parts
 *            of dim series of degree deg;
 *   dataimhi has space allocated for the high doubles of the imaginary parts
 *            of dim series of degree deg;
 *   dataimlo has space allocated for the high doubles of the imaginary parts
 *            of dim series of degree deg.
 *
 * ON RETURN :
 *   datarehi contains the high doubles of the real parts of the data,
 *   datarehi contains the low doubles of the real parts of the data,
 *   dataimhi contains the high doubles of the imaginary parts of the data,
 *   dataimlo contains the low doubles of the imaginary parts of the data,
 *            data[i][j] is the j-th coefficient of the i-th series,
 *            for i in 0..dim-1 and j in 0..deg. */

int test_real ( int dim, int nvr, int pwr, int deg );
/*
 * DESCRIPTION :
 *   Tests the evaluation and differentiation for random real data.
 * 
 * ON ENTRY :
 *   dim      dimension, total number of variables;
 *   nvr      number of variables in the product;
 *   pwr      highest power of each variable;
 *   deg      truncation degree of the series. */

int test_complex ( int dim, int nvr, int pwr, int deg );
/*
 * DESCRIPTION :
 *   Tests the evaluation and differentiation for random complex data.
 * 
 * ON ENTRY :
 *   dim      dimension, total number of variables;
 *   nvr      number of variables in the product;
 *   pwr      highest power of each variable;
 *   deg      truncation degree of the series. */

int main ( void )
{
   int dim,nvr,pwr,deg;

   cout << "Give the dimension : "; cin >> dim;
   cout << "Give the number of variables, <= "; cout << dim;
   cout << " : "; cin >> nvr;
   cout << "Give the largest power of each variable : "; cin >> pwr;
   cout << "Give the degree of the series : "; cin >> deg;

   cout << endl << "Testing for real input data ... " << endl;
   test_real(dim,nvr,pwr,deg);
   // cout << endl << "Testing for complex input data ..." << endl;
   // test_complex(dim,nvr,pwr,deg);

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
 ( int dim, int nvr, int pwr, int deg, int *idx, int *exp,
   double *cffhi, double *cfflo )
{
   bool fail;

   if(nvr > dim)
   {
      cout << "ERROR: nvr = " << nvr << " > " << dim << " dim" << endl;
      return true;
   }
   else
   {
      for(int i=0; i<=deg; i++) random_double_double(&cffhi[i],&cfflo[i]);

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

bool make_complex_monomial
 ( int dim, int nvr, int pwr, int deg, int *idx, int *exp,
   double *cffrehi, double *cffrelo, double *cffimhi, double *cffimlo )
{
   bool fail;

   if(nvr > dim)
   {
      cout << "ERROR: nvr = " << nvr << " > " << dim << " dim" << endl;
      return true;
   }
   else
   {
      double rnd;

      for(int i=0; i<=deg; i++)
      {
         rnd = random_angle();        
         cffrehi[i] = cos(rnd); cffrelo[i] = 0.0;
         cffimhi[i] = sin(rnd); cffimlo[i] = 0.0;
      }

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

void make_real_input ( int dim, int deg, double **datahi, double **datalo )
{
   for(int i=0; i<dim; i++)
      for(int j=0; j<=deg; j++)
         random_double_double(&datahi[i][j],&datalo[i][j]);
}

void make_complex_input
 ( int dim, int deg, double **datarehi, double **datarelo,
   double **dataimhi, double **dataimlo )
{
   double rnd;

   for(int i=0; i<dim; i++)
      for(int j=0; j<=deg; j++)
      {
         rnd = random_angle();
         datarehi[i][j] = cos(rnd); datarelo[i][j] = 0.0;
         dataimhi[i][j] = sin(rnd); dataimlo[i][j] = 0.0;
      }
}

int test_real ( int dim, int nvr, int pwr, int deg )
{
   int *idx = new int[nvr];           // indices of variables in the monomial
   int *exp = new int[nvr];           // exponents of the variables
   int *expfac = new int[nvr];        // exponents of common factor
   int nbrfac;                        // number of common factors
   double *cffhi = new double[deg+1]; // high doubles of series coefficient
   double *cfflo = new double[deg+1]; // low doubles of series coefficient

   // The input are dim power series of degree deg,
   // the output are nvr+1 power series of degree deg,
   // for the evaluated and differentiated monomial.

   double **inputhi = new double*[dim];
   for(int i=0; i<dim; i++) inputhi[i] = new double[deg+1];
   double **inputlo = new double*[dim];
   for(int i=0; i<dim; i++) inputlo[i] = new double[deg+1];
   double **outputhi_h = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputhi_h[i] = new double[deg+1];
   double **outputhi_d = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputhi_d[i] = new double[deg+1];
   double **outputlo_h = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputlo_h[i] = new double[deg+1];
   double **outputlo_d = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputlo_d[i] = new double[deg+1];

   srand(time(NULL));

   bool fail = make_real_monomial(dim,nvr,pwr,deg,idx,exp,cffhi,cfflo);

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
      for(int i=0; i<=deg; i++)
         cout << " " << cffhi[i] << " " << cfflo[i] << endl;;
   }

   make_real_input(dim,deg,inputhi,inputlo);

   cout << "Random input series :" << endl;
   for(int i=0; i<dim; i++)
   {
      cout << "-> coefficients of series " << i << " :" << endl;
      for(int j=0; j<=deg; j++)
         cout << inputhi[i][j] << "  " << inputlo[i][j] << endl;
   }

   CPU_dbl2_evaldiff(dim,nvr,deg,idx,cffhi,cfflo,inputhi,inputlo,
                     outputhi_h,outputlo_h);
   cout << "The value of the product :" << endl;
   for(int i=0; i<=deg; i++)
      cout << outputhi_h[dim][i] << "  " << outputlo_h[dim][i] << endl;
 
   if(nvr > 2)
   {
      GPU_dbl2_evaldiff(deg+1,dim,nvr,deg,idx,cffhi,cfflo,inputhi,inputlo,
                        outputhi_d,outputlo_d);
      cout << "The value of the product computed on the GPU :" << endl;
      for(int i=0; i<=deg; i++)
         cout << outputhi_d[dim][i] << "  "
              << outputlo_d[dim][i] << endl;
   }
   for(int k=0; k<nvr; k++)
   {
      cout << "-> derivative for index " << idx[k] << " :" << endl;
      for(int i=0; i<=deg; i++)
         cout << outputhi_h[idx[k]][i] << "  "
              << outputlo_h[idx[k]][i] << endl;

      if(nvr > 2)
      {
         cout << "-> derivative for index " << idx[k]
              << " computed on GPU :" << endl;
         for(int i=0; i<=deg; i++)
            cout << outputhi_d[idx[k]][i] << "  "
                 << outputlo_d[idx[k]][i] << endl;
      }
   }
   return 0;
}

int test_complex ( int dim, int nvr, int pwr, int deg )
{
   int *idx = new int[nvr];           // indices of variables in the monomial
   int *exp = new int[nvr];           // exponents of the variables
   int *expfac = new int[nvr];        // exponents of common factor
   int nbrfac;                        // number of common factors
   double *cffre = new double[deg+1]; // real parts of coefficients
   double *cffim = new double[deg+1]; // imaginary parts of coefficients

   // The input are dim power series of degree deg,
   // the output are nvr+1 power series of degree deg,
   // for the evaluated and differentiated monomial.
/*
   double **inputre = new double*[dim];
   double **inputim = new double*[dim];
   for(int i=0; i<dim; i++)
   {
      inputre[i] = new double[deg+1];
      inputim[i] = new double[deg+1];
   }
   double **outputre_h = new double*[dim+1];
   double **outputim_h = new double*[dim+1];
   for(int i=0; i<=dim; i++)
   {
      outputre_h[i] = new double[deg+1];
      outputim_h[i] = new double[deg+1];
   }
   double **outputre_d = new double*[dim+1];
   double **outputim_d = new double*[dim+1];
   for(int i=0; i<=dim; i++)
   {
      outputre_d[i] = new double[deg+1];
      outputim_d[i] = new double[deg+1];
   }

   srand(time(NULL));

   bool fail = make_complex_monomial(dim,nvr,pwr,deg,idx,exp,cffre,cffim);

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
      for(int i=0; i<=deg; i++)
         cout << " " << cffre[i] << "  " << cffim[i] << endl;;
   }

   make_complex_input(dim,deg,inputre,inputim);

   cout << "Random input series :" << endl;
   for(int i=0; i<dim; i++)
   {
      cout << "-> coefficients of series " << i << " :" << endl;
      for(int j=0; j<=deg; j++)
         cout << inputre[i][j] << "  " << inputim[i][j] << endl;
   }

   CPU_cmplx_evaldiff(dim,nvr,deg,idx,cffre,cffim,inputre,inputim,
                      outputre_h,outputim_h);

   cout << "The value of the product :" << endl;
   for(int i=0; i<=deg; i++)
      cout << outputre_h[dim][i] << "  " << outputim_h[dim][i] << endl;
   if(nvr > 2)
   {
      GPU_cmplx_evaldiff(deg+1,dim,nvr,deg,idx,cffre,cffim,inputre,inputim,
                         outputre_d,outputim_d);

      cout << "The value of the product computed on the GPU :" << endl;
      for(int i=0; i<=deg; i++) 
         cout << outputre_d[dim][i] << "  " << outputim_d[dim][i] << endl;
   }
   for(int k=0; k<nvr; k++)
   {
      cout << "-> derivative for index " << idx[k] << " :" << endl;
      for(int i=0; i<=deg; i++)
         cout << outputre_h[idx[k]][i] << "  "
              << outputim_h[idx[k]][i] << endl;
      if(nvr > 2)
      {
         cout << "-> derivative for index " << idx[k]
              << " computed on GPU :" << endl;
         for(int i=0; i<=deg; i++)
            cout << outputre_d[idx[k]][i] << "  "
                 << outputim_d[idx[k]][i] << endl;
      }
   }

 */
   return 0;
}
