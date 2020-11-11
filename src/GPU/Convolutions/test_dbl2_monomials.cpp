// Tests the evaluation and differentiation of a monomial
// in double double precision.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <vector_types.h>
#include "random_monomials.h"
#include "random2_monomials.h"
#include "dbl2_monomials_host.h"
#include "dbl2_monomials_kernels.h"

using namespace std;

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
   cout << endl << "Testing for complex input data ..." << endl;
   test_complex(dim,nvr,pwr,deg);

   return 0;
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

   bool fail = make_real2_monomial(dim,nvr,pwr,deg,idx,exp,cffhi,cfflo);

   if(fail) return 1;

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

   make_real2_input(dim,deg,inputhi,inputlo);

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

   double errsum = 0.0;
   double errtot = 0.0;
 
   if(nvr > 2)
   {

      GPU_dbl2_evaldiff(deg+1,dim,nvr,deg,idx,cffhi,cfflo,inputhi,inputlo,
                        outputhi_d,outputlo_d);
      cout << "The value of the product computed on the GPU :" << endl;
      for(int i=0; i<=deg; i++)
      {
         cout << outputhi_d[dim][i] << "  " << outputlo_d[dim][i] << endl;
         errsum = errsum
                + abs(outputhi_h[dim][i] - outputhi_d[dim][i])
                + abs(outputlo_h[dim][i] - outputlo_d[dim][i]);
      }
      cout << "Sum of errors : " << errsum << endl; errtot += errsum;
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

         errsum = 0.0;
         for(int i=0; i<=deg; i++)
         {
            cout << outputhi_d[idx[k]][i] << "  "
                 << outputlo_d[idx[k]][i] << endl;
            errsum = errsum
                   + abs(outputhi_h[idx[k]][i] - outputhi_d[idx[k]][i])
                   + abs(outputlo_h[idx[k]][i] - outputlo_d[idx[k]][i]);
            cout << "Sum of errors : " << errsum << endl; errtot += errsum;
         }
      }
   }
   cout << "Total sum of all errors : " << errtot << endl;

   return 0;
}

int test_complex ( int dim, int nvr, int pwr, int deg )
{
   int *idx = new int[nvr];             // indices of variables in the monomial
   int *exp = new int[nvr];             // exponents of the variables
   int *expfac = new int[nvr];          // exponents of common factor
   int nbrfac;                          // number of common factors
   double *cffrehi = new double[deg+1]; // high real parts of coefficients
   double *cffrelo = new double[deg+1]; // low real parts of coefficients
   double *cffimhi = new double[deg+1]; // high imaginary parts of coefficients
   double *cffimlo = new double[deg+1]; // low imaginary parts of coefficients

   // The input are dim power series of degree deg,
   // the output are nvr+1 power series of degree deg,
   // for the evaluated and differentiated monomial.

   double **inputrehi = new double*[dim];
   double **inputrelo = new double*[dim];
   double **inputimhi = new double*[dim];
   double **inputimlo = new double*[dim];
   for(int i=0; i<dim; i++)
   {
      inputrehi[i] = new double[deg+1];
      inputrelo[i] = new double[deg+1];
      inputimhi[i] = new double[deg+1];
      inputimlo[i] = new double[deg+1];
   }
   double **outputrehi_h = new double*[dim+1];
   double **outputrelo_h = new double*[dim+1];
   double **outputimhi_h = new double*[dim+1];
   double **outputimlo_h = new double*[dim+1];
   for(int i=0; i<=dim; i++)
   {
      outputrehi_h[i] = new double[deg+1];
      outputrelo_h[i] = new double[deg+1];
      outputimhi_h[i] = new double[deg+1];
      outputimlo_h[i] = new double[deg+1];
   }
   double **outputrehi_d = new double*[dim+1];
   double **outputrelo_d = new double*[dim+1];
   double **outputimhi_d = new double*[dim+1];
   double **outputimlo_d = new double*[dim+1];
   for(int i=0; i<=dim; i++)
   {
      outputrehi_d[i] = new double[deg+1];
      outputrelo_d[i] = new double[deg+1];
      outputimhi_d[i] = new double[deg+1];
      outputimlo_d[i] = new double[deg+1];
   }
   srand(time(NULL));

   bool fail = make_complex2_monomial
     (dim,nvr,pwr,deg,idx,exp,cffrehi,cffrelo,cffimhi,cffimlo);

   if(fail) return 1;

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
   {
      cout << cffrehi[i] << "  " << cffrelo[i] << endl;
      cout << cffimhi[i] << "  " << cffimlo[i] << endl;
   }

   make_complex2_input(dim,deg,inputrehi,inputrelo,inputimhi,inputimlo);

   cout << "Random input series :" << endl;
   for(int i=0; i<dim; i++)
   {
      cout << "-> coefficients of series " << i << " :" << endl;
      for(int j=0; j<=deg; j++)
      {
         cout << inputrehi[i][j] << "  " << inputrelo[i][j] << endl;
         cout << inputimhi[i][j] << "  " << inputimlo[i][j] << endl;
      }
   }
   CPU_cmplx2_evaldiff
      (dim,nvr,deg,idx,cffrehi,cffrelo,cffimhi,cffimlo,
       inputrehi,inputrelo,inputimhi,inputimlo,
       outputrehi_h,outputrelo_h,outputimhi_h,outputimlo_h);

   cout << "The value of the product :" << endl;
   for(int i=0; i<=deg; i++)
   {
      cout << outputrehi_h[dim][i] << "  " << outputrelo_h[dim][i] << endl;
      cout << outputimhi_h[dim][i] << "  " << outputimlo_h[dim][i] << endl;
   }

   double errsum = 0.0;
   double errtot = 0.0;

   if(nvr > 2)
   {
      GPU_cmplx2_evaldiff(deg+1,dim,nvr,deg,idx,
         cffrehi,cffrelo,cffimhi,cffimlo,
         inputrehi,inputrelo,inputimhi,inputimlo,
         outputrehi_d,outputrelo_d,outputimhi_d,outputimlo_d);

      cout << "The value of the product computed on the GPU :" << endl;
      for(int i=0; i<=deg; i++) 
      {
         cout << outputrehi_d[dim][i] << "  " << outputrelo_d[dim][i] << endl;
         cout << outputimhi_d[dim][i] << "  " << outputimlo_d[dim][i] << endl;
         errsum = errsum
                + abs(outputrehi_h[dim][i] - outputrehi_d[dim][i])
                + abs(outputrelo_h[dim][i] - outputrelo_d[dim][i])
                + abs(outputimhi_h[dim][i] - outputimhi_d[dim][i])
                + abs(outputimlo_h[dim][i] - outputimlo_d[dim][i]);
      }
      cout << "The sum of errors : " << errsum << endl; errtot += errsum;
   }

   for(int k=0; k<nvr; k++)
   {
      cout << "-> derivative for index " << idx[k] << " :" << endl;
      for(int i=0; i<=deg; i++)
      {
         cout << outputrehi_h[idx[k]][i] << "  "
              << outputrelo_h[idx[k]][i] << endl;
         cout << outputimhi_h[idx[k]][i] << "  "
              << outputimlo_h[idx[k]][i] << endl;
      }
      if(nvr > 2)
      {
         cout << "-> derivative for index " << idx[k]
              << " computed on GPU :" << endl;

         errsum = 0.0;
         for(int i=0; i<=deg; i++)
         {
            cout << outputrehi_d[idx[k]][i] << "  "
                 << outputrelo_d[idx[k]][i] << endl;
            cout << outputimhi_d[idx[k]][i] << "  "
                 << outputimlo_d[idx[k]][i] << endl;
            errsum = errsum
                   + abs(outputrehi_h[idx[k]][i] - outputrehi_d[idx[k]][i])
                   + abs(outputrelo_h[idx[k]][i] - outputrelo_d[idx[k]][i])
                   + abs(outputimhi_h[idx[k]][i] - outputimhi_d[idx[k]][i])
                   + abs(outputimlo_h[idx[k]][i] - outputimlo_d[idx[k]][i]);
         }
         cout << "The sum of errors : " << errsum << endl; errtot += errsum;
      }
   }
   cout << "Total sum of all errors : " << errtot << endl;

   return 0;
}
