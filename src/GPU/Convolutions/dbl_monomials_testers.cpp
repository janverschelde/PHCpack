/* The file dbl_monomials_testers.cpp contains the definitions of
 * functions with prototypes in dbl_monomials_testers.h. */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <vector_types.h>
#include "random_monomials.h"
#include "dbl_monomials_host.h"
#include "dbl_monomials_kernels.h"
#include "dbl_monomials_testers.h"

using namespace std;

int main_dbl_test
 ( int seed, int dim, int nvr, int pwr, int deg, int vrblvl )
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
   double realsum = test_dbl_real(dim,nvr,pwr,deg,vrblvl-1);
   double complexsum = test_dbl_complex(dim,nvr,pwr,deg,vrblvl-1);

   const double tol = 1.0e-12;

   int fail = int(realsum > tol) + int(complexsum > tol);

   if(vrblvl > 0)
   {
      cout << scientific << setprecision(2);
      cout << "Sum of all errors in double precision :" << endl;
      cout << "  on real data : " << realsum;
      if(realsum < tol)
         cout << "  pass," << endl;
      else
         cout << "  fail!" << endl;

      cout << "  on complex data : " << complexsum;
      if(complexsum < tol)
         cout << "  pass.";
      else
         cout << "  fail!";

      cout << "  Seed used : " <<  seedused << endl;
   }
   return fail;
}

double test_dbl_real ( int dim, int nvr, int pwr, int deg, int verbose )
{
   int *idx = new int[nvr];         // indices of variables in the monomial
   int *exp = new int[nvr];         // exponents of the variables
   int *expfac = new int[nvr];      // exponents of common factor
   int nbrfac;                      // number of common factors
   double *cff = new double[deg+1]; // series coefficient

   // The input are dim power series of degree deg,
   // the output are nvr+1 power series of degree deg,
   // for the evaluated and differentiated monomial.

   double **input = new double*[dim];
   for(int i=0; i<dim; i++) input[i] = new double[deg+1];
   double **output_h = new double*[dim+1];
   for(int i=0; i<=dim; i++) output_h[i] = new double[deg+1];
   double **output_d = new double*[dim+1];
   for(int i=0; i<=dim; i++) output_d[i] = new double[deg+1];

   bool fail = make_real_monomial(dim,nvr,pwr,deg,idx,exp,cff);

   if(fail) return 1.0;
 
   if(verbose > 0)
   {
      cout << "Generated a random monomial :" << endl;
      cout << "   the indices :";
      for(int i=0; i<nvr; i++) cout << " " << idx[i];
      cout << endl;

      cout << " the exponents :";
      for(int i=0; i<nvr; i++) cout << " " << exp[i];
      cout << endl;
   }
   common_factors(nvr,exp,&nbrfac,expfac);

   if(verbose > 0)
   {
      cout << "common factors :";
      for(int i=0; i<nvr; i++) cout << " " << expfac[i];
      cout << endl;
      cout << "number of common factors : " << nbrfac << endl;

      cout << scientific << setprecision(16);
      cout << "the coefficients :" << endl;
      for(int i=0; i<=deg; i++) cout << " " << cff[i] << endl;;
   }
   make_real_input(dim,deg,input);

   if(verbose > 0)
   {
      cout << "Random input series :" << endl;
      for(int i=0; i<dim; i++)
      {
         cout << "-> coefficients of series " << i << " :" << endl;
         for(int j=0; j<=deg; j++) cout << input[i][j] << endl;
      }
   }
   double errsum = 0.0;
   double errtot = 0.0;

   CPU_dbl_evaldiff(dim,nvr,deg,idx,cff,input,output_h);

   if(verbose > 0)
   {
      cout << "The value of the product :" << endl;
      for(int i=0; i<=deg; i++) cout << output_h[dim][i] << endl;
   }
   if(nvr == dim) // the product of all input series equals one
   {
      for(int i=0; i<=deg; i++)
         errsum = errsum + abs(output_h[dim][i] - cff[i]);

      if(verbose > 0) cout << "Coefficient error : " << errsum << endl;

      errtot += errsum;
   }
   if(nvr > 2)
   {
      GPU_dbl_evaldiff(deg+1,dim,nvr,deg,idx,cff,input,output_d);

      if(verbose > 0)
      {
         cout << "The value of the product computed on the GPU :" << endl;
      }
      errsum = 0.0;
      for(int i=0; i<=deg; i++)
      {
         if(verbose > 0) cout << output_d[dim][i] << endl;

         errsum = errsum + abs(output_h[dim][i] - output_d[dim][i]);
      }
      if(verbose > 0) cout << "Sum of errors : " << errsum << endl;

      errtot += errsum;

      if(nvr == dim) // the product of all input series equals one
      {
         errsum = 0.0;
         for(int i=0; i<=deg; i++)
            errsum = errsum + abs(output_d[dim][i] - cff[i]);

         if(verbose > 0) cout << "Coefficient error : " << errsum << endl;

         errtot += errsum;
      }
   }
   for(int k=0; k<nvr; k++)
   {
      if(verbose > 0)
      {
         cout << "-> derivative for index " << idx[k] << " :" << endl;
         for(int i=0; i<=deg; i++) cout << output_h[idx[k]][i] << endl;
      }
      if(nvr > 2)
      {
         if(verbose > 0)
         {
            cout << "-> derivative for index " << idx[k]
                 << " computed on GPU :" << endl;
         }
         errsum = 0.0;
         for(int i=0; i<=deg; i++)
         {
            if(verbose > 0) cout << output_d[idx[k]][i] << endl;

            errsum = errsum + abs(output_h[idx[k]][i] - output_d[idx[k]][i]);
         }
         if(verbose > 0) cout << "Sum of errors : " << errsum << endl;

         errtot += errsum;
      }
   }
   if(verbose > 0) cout << "Total sum of all errors : " << errtot << endl;

   return errtot;
}

double test_dbl_complex ( int dim, int nvr, int pwr, int deg, int verbose )
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
   bool fail = make_complex_monomial(dim,nvr,pwr,deg,idx,exp,cffre,cffim);

   if(fail) return 1.0;

   if(verbose > 0)
   {
      cout << "Generated a random monomial :" << endl;
      cout << "   the indices :";
      for(int i=0; i<nvr; i++) cout << " " << idx[i];
      cout << endl;

      cout << " the exponents :";
      for(int i=0; i<nvr; i++) cout << " " << exp[i];
      cout << endl;
   }
   common_factors(nvr,exp,&nbrfac,expfac);

   if(verbose > 0)
   {
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

   if(verbose > 0)
   {
      cout << "Random input series :" << endl;
      for(int i=0; i<dim; i++)
      {
         cout << "-> coefficients of series " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
            cout << inputre[i][j] << "  " << inputim[i][j] << endl;
      }
   }
   double errsum = 0.0;
   double errtot = 0.0;

   CPU_cmplx_evaldiff(dim,nvr,deg,idx,cffre,cffim,inputre,inputim,
                      outputre_h,outputim_h);
   if(verbose > 0)
   {
      cout << "The value of the product :" << endl;
      for(int i=0; i<=deg; i++)
         cout << outputre_h[dim][i] << "  " << outputim_h[dim][i] << endl;
   }
   if(nvr == dim) // the product of all input series equals one
   {
      for(int i=0; i<=deg; i++)
         errsum = errsum
                + abs(outputre_h[dim][i] - cffre[i])
                + abs(outputim_h[dim][i] - cffim[i]);

      if(verbose > 0) cout << "Coefficient error : " << errsum << endl;

      errtot += errsum;
   }
   if(nvr > 2)
   {
      GPU_cmplx_evaldiff(deg+1,dim,nvr,deg,idx,cffre,cffim,inputre,inputim,
                         outputre_d,outputim_d);
      errsum = 0.0;

      if(verbose > 0)
      {
         cout << "The value of the product computed on the GPU :" << endl;
      }
      for(int i=0; i<=deg; i++) 
      {
         if(verbose > 0)
         {
            cout << outputre_d[dim][i] << "  " << outputim_d[dim][i] << endl;
         }
         errsum = errsum
                + abs(outputre_h[dim][i] - outputre_d[dim][i])
                + abs(outputim_h[dim][i] - outputim_d[dim][i]);
      }
      if(verbose > 0) cout << "The sum of errors : " << errsum << endl;

      errtot += errsum;

      if(nvr == dim) // the product of all input series equals one
      {
         errsum = 0.0;
         for(int i=0; i<=deg; i++)
            errsum = errsum
                   + abs(outputre_d[dim][i] - cffre[i])
                   + abs(outputim_d[dim][i] - cffim[i]);

         if(verbose > 0) cout << "Coefficient error : " << errsum << endl;

         errtot += errsum;
      }
   }
   for(int k=0; k<nvr; k++)
   {
      if(verbose > 0)
      {
         cout << "-> derivative for index " << idx[k] << " :" << endl;
         for(int i=0; i<=deg; i++)
            cout << outputre_h[idx[k]][i] << "  "
                 << outputim_h[idx[k]][i] << endl;
      }
      if(nvr > 2)
      {
         if(verbose > 0)
         {
            cout << "-> derivative for index " << idx[k]
                 << " computed on GPU :" << endl;
         }
         errsum = 0.0;
         for(int i=0; i<=deg; i++)
         {
            if(verbose > 0)
            {
               cout << outputre_d[idx[k]][i] << "  "
                    << outputim_d[idx[k]][i] << endl;
            }
            errsum = errsum
                   + abs(outputre_h[idx[k]][i] - outputre_d[idx[k]][i])
                   + abs(outputim_h[idx[k]][i] - outputim_d[idx[k]][i]);
         }
         if(verbose > 0) cout << "The sum of errors : " << errsum << endl;

         errtot += errsum;
      }
   }
   if(verbose > 0) cout << "Total sum of all errors : " << errtot << endl;

   return errtot;
}
