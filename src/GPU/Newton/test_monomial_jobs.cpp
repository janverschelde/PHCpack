/* Computes the jobs to evaluate and differentiate a monomial system. */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <time.h>
#include "convolution_jobs.h"
#include "unimodular_matrices.h"

using namespace std;

void prompt_monomial_setup
 ( int *seed, int *dim, int *deg, int *nbritr, int *size, int *vrblvl );
/*
 * DESCRIPTION :
 *   Prompts for the parameters to define a monomial system.
 *
 * ON RETURN :
 *   seed      the seed for the random number generator (0 for time);
 *   dim       the dimension is the number of monomials
 *             and the maximum number of variables in each monomial;
 *   deg       degree of the series;
 *   nbritr    number of unimodular multiplications in the making of
 *             the exponent matrix;
 *   size      size of the numbers (only if nbritr > 0);
 *   vrblvl    verbose level (0 if silent). */

void list_convolution_jobs
 ( int dim, int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   int *count, int vrblvl );
/*
 * DESCRIPTION :
 *   Lists and counts the convolution jobs.
 *
 * ON ENTRY :
 *   dim       number of monomials;
 *   nvr       nvr[i] is the number of variables in the i-th monomial;
 *   idx       idx[i] are the indices of the variables in monomial i;
 *   exp       exp[i] are the exponents of the variables in monomial i;
 *   nbrfac    nbrfac[i] are the number of exponents > 1 in monomial i;
 *   expfac    expfac[i] are the exponents in the i-th polynomial
 *             that are larger than one, minus one in the factor,
 *             if exp[i][k] > 1, then expfac[i][k] = exp[i][k] - 1;
 *   vrblvl    is the verbose level.
 *
 * ON RETURN :
 *   count     the number of convolution jobs. */

void list_multiplication_jobs
 ( int dim, int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   int *count, int vrblvl );
/*
 * DESCRIPTION :
 *   Lists and counts the multiplication jobs.
 *
 * ON ENTRY :
 *   dim       number of monomials;
 *   nvr       nvr[i] is the number of variables in the i-th monomial;
 *   idx       idx[i] are the indices of the variables in monomial i;
 *   exp       exp[i] are the exponents of the variables in monomial i;
 *   nbrfac    nbrfac[i] are the number of exponents > 1 in monomial i;
 *   expfac    expfac[i] are the exponents in the i-th polynomial
 *             that are larger than one, minus one in the factor,
 *             if exp[i][k] > 1, then expfac[i][k] = exp[i][k] - 1;
 *   vrblvl    is the verbose level.
 *
 * ON RETURN :
 *   count     the number of multiplication jobs. */

int main ( void )
{
   srand(time(NULL));

   cout << "testing jobs to evaluate and differentiate monomial system ..."
        << endl;
/*
 * 0. defines the parameters of the monomial system
 */
   int seed,dim,deg,size,nbritr,vrblvl;

   prompt_monomial_setup(&seed,&dim,&deg,&nbritr,&size,&vrblvl);

   if(seed == 0) seed = time(NULL);
   srand(seed);

   cout << "making an integer matrix of dimension " << dim
        << " ..." << endl;
/*
 * 1. defining the unimodular matrix of monomials
 */
   int **rowsA = new int*[dim];  // exponents in the rows
   int *nvr = new int[dim];      // number of variables in each monomial
   int **idx = new int*[dim];    // indexes of variables in each monomial
   int **exp = new int*[dim];    // exponents of the variables
   int *nbrfac = new int[dim];   // number of exponents > 1 in each monomial
   int **expfac = new int*[dim]; // exponents of the common factors

   make_monomial_system
      (dim,size,1,nbritr,nvr,idx,exp,nbrfac,expfac,rowsA,vrblvl);
/*
 * 2. list the convolution and multiplication jobs 
 */
   int cnvjobcnt;
   bool verbose = (vrblvl > 0);

   list_convolution_jobs(dim,nvr,idx,exp,nbrfac,expfac,&cnvjobcnt,vrblvl);

   cout << "number of convolution jobs : " << cnvjobcnt << endl;

   ConvolutionJobs jobs(dim);

   jobs.make(dim,nvr,idx,verbose);

   for(int k=0; k<jobs.get_depth(); k++)
   {
      cout << "jobs at layer " << k << " :" << endl;
      for(int i=0; i<jobs.get_layer_count(k); i++)
         cout << jobs.get_job(k,i) << endl;
   }
   cout << "dimension : " << dim << endl;
   cout << "number of monomials : " << dim << endl;
   cout << "number of convolution jobs : " << jobs.get_count() << endl;
   cout << "number of layers : " << jobs.get_depth() << endl;
   cout << "frequency of layer counts :" << endl;
   int checksum = 0;
   for(int i=0; i<jobs.get_depth(); i++)
   {
      cout << i << " : " << jobs.get_layer_count(i) << endl;
      checksum = checksum + jobs.get_layer_count(i); 
   }
   cout << "layer count sum : " << checksum << endl;

   int muljobcnt;

   list_multiplication_jobs(dim,nvr,idx,exp,nbrfac,expfac,&muljobcnt,vrblvl);

   cout << "number of multiplication jobs : " << muljobcnt << endl;

   cout << "seed used : " << seed << endl;

   return 0;
}

void prompt_monomial_setup
 ( int *seed, int *dim, int *deg, int *nbritr, int *size, int *vrblvl )
{
   cout << "-> give the seed (0 for time) : "; cin >> *seed;
   cout << "-> give the dimension : "; cin >> *dim;
   cout << "-> give the degree of the series : "; cin >> *deg;
   cout << "-> give the number of unimodular multiplications (0 for input) : ";
   cin >> *nbritr;
   if(*nbritr == 0) *size = 0;
   if(*nbritr > 0)
   {
      cout << "-> give the size of the numbers : ";
      cin >> *size;
   }
   cout << "-> verbose level (0 for silent) : "; cin >> *vrblvl;
}

void list_convolution_jobs
 ( int dim, int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   int *count, int vrblvl )
{
   *count = 0;
   for(int i=0; i<dim; i++) // common factors in the coefficients
   {
      if(nbrfac[i] > 0) // there are common factors in monomial i
      {
         for(int j=0; j<nvr[i]; j++) // run over all exponents
         {
            if(expfac[i][j] > 0) // the j-th exponent with variable idx[i][j]
            {
               int idxvar = idx[i][j];

               for(int k=0; k<expfac[i][j]; k++)
               {
                  if(vrblvl > 0)
                     cout << "multiplying input " << idxvar
                          << " with coefficient " << i << endl;
                  // CPU_dbl_product(deg,input[idxvar],cff[i],acc);
                  *count = *count + 1;             
               }
            }
         }
      }
   }
}

void list_multiplication_jobs
 ( int dim, int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   int *count, int vrblvl )
{
   *count = 0;

   for(int i=0; i<dim; i++) // multiply derivatives with the powers
   {
      if(nbrfac[i] > 0) // there are common factors in monomial i
      {
         for(int j=0; j<nvr[i]; j++) // run over all exponents
         {
            if(expfac[i][j] > 0) // the j-th exponent with variable idx[i][j]
            {
               int idxvar = idx[i][j];
               // double factor = (double) exp[i][j];
               // multiply derivative w.r.t. idxvar with factor
               // for(int k=0; k<=deg; k++)
               //   output[i][idxvar][k] = factor*output[i][idxvar][k];
               if(vrblvl > 0)
                  cout << "multiplying output[" << i << "][" << idxvar
                       << "] with " << exp[i][j] << endl; 
               *count = *count + 1;
            }
         }
      }
   }
}
