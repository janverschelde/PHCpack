/* The file dbl3_monomials_testers.cpp contains the definitions of
 * functions with prototypes in dbl3_monomials_testers.h. */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <vector_types.h>
#include "random_monomials.h"
#include "random3_monomials.h"
#include "dbl3_monomials_host.h"
#include "dbl3_monomials_kernels.h"
#include "dbl3_monomials_testers.h"

using namespace std;

int main_dbl3_test
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
   double realsum = test_dbl3_real(dim,nvr,pwr,deg,vrblvl-1);
   double complexsum = test_dbl3_complex(dim,nvr,pwr,deg,vrblvl-1);

   const double tol = 1.0e-44;

   int fail = int(realsum > tol) + int(complexsum > tol);

   if(vrblvl > 0)
   {
      cout << "Sum of all errors in triple double precision :" << endl;
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

double test_dbl3_real ( int dim, int nvr, int pwr, int deg, int verbose )
{
   int *idx = new int[nvr];           // indices of variables in the monomial
   int *exp = new int[nvr];           // exponents of the variables
   int *expfac = new int[nvr];        // exponents of common factor
   int nbrfac;                        // number of common factors
   double *cffhi = new double[deg+1]; // high doubles of series coefficient
   double *cffmi = new double[deg+1]; // middle doubles of series coefficient
   double *cfflo = new double[deg+1]; // low doubles of series coefficient

   // The input are dim power series of degree deg,
   // the output are nvr+1 power series of degree deg,
   // for the evaluated and differentiated monomial.

   double **inputhi = new double*[dim];
   for(int i=0; i<dim; i++) inputhi[i] = new double[deg+1];
   double **inputmi = new double*[dim];
   for(int i=0; i<dim; i++) inputmi[i] = new double[deg+1];
   double **inputlo = new double*[dim];
   for(int i=0; i<dim; i++) inputlo[i] = new double[deg+1];
   double **outputhi_h = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputhi_h[i] = new double[deg+1];
   double **outputhi_d = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputhi_d[i] = new double[deg+1];
   double **outputmi_h = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputmi_h[i] = new double[deg+1];
   double **outputmi_d = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputmi_d[i] = new double[deg+1];
   double **outputlo_h = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputlo_h[i] = new double[deg+1];
   double **outputlo_d = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputlo_d[i] = new double[deg+1];

   bool fail = make_real3_monomial(dim,nvr,pwr,deg,idx,exp,cffhi,cffmi,cfflo);

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
         cout << cffhi[i] << "  " << cffmi[i] << "  " << cfflo[i] << endl;
   }
   make_real3_input(dim,deg,inputhi,inputmi,inputlo);

   if(verbose > 0)
   {
      cout << "Random input series :" << endl;
      for(int i=0; i<dim; i++)
      {
         cout << "-> coefficients of series " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
            cout << inputhi[i][j] << "  " << inputmi[i][j]
                                  << "  " << inputlo[i][j] << endl;
      }
   }
   double errsum = 0.0;
   double errtot = 0.0;

   CPU_dbl3_evaldiff(dim,nvr,deg,idx,cffhi,cffmi,cfflo,
                     inputhi,inputmi,inputlo,
                     outputhi_h,outputmi_h,outputlo_h);
   if(verbose > 0)
   {
      cout << "The value of the product :" << endl;
      for(int i=0; i<=deg; i++)
         cout << outputhi_h[dim][i] << "  " << outputmi_h[dim][i] 
                                    << "  " << outputlo_h[dim][i] << endl;
   }
   if(nvr == dim) // the product of all input series equals one
   {
      for(int i=0; i<=deg; i++)
         errsum = errsum
                + abs(outputhi_h[dim][i] - cffhi[i])
                + abs(outputmi_h[dim][i] - cffmi[i])
                + abs(outputlo_h[dim][i] - cfflo[i]);

      if(verbose > 0) cout << "Coefficient error : " << errsum << endl;

      errtot += errsum;
   }
   if(nvr > 2)
   {
      GPU_dbl3_evaldiff(deg+1,dim,nvr,deg,idx,cffhi,cffmi,cfflo,inputhi,
                        inputmi,inputlo,outputhi_d,outputmi_d,outputlo_d);

      if(verbose > 0)
         cout << "The value of the product computed on the GPU :" << endl;

      for(int i=0; i<=deg; i++)
      {
         if(verbose > 0)
         {
            cout << outputhi_d[dim][i] << "  "
                 << outputmi_d[dim][i] << "  "
                 << outputlo_d[dim][i] << endl;
         }
         errsum = errsum
                + abs(outputhi_h[dim][i] - outputhi_d[dim][i])
                + abs(outputmi_h[dim][i] - outputmi_d[dim][i])
                + abs(outputlo_h[dim][i] - outputlo_d[dim][i]);
      }
      if(verbose > 0) cout << "Sum of errors : " << errsum << endl;

      errtot += errsum;

      if(nvr == dim) // the product of all input series equals one
      {
         errsum = 0.0;
         for(int i=0; i<=deg; i++)
            errsum = errsum
                   + abs(outputhi_d[dim][i] - cffhi[i])
                   + abs(outputmi_d[dim][i] - cffmi[i])
                   + abs(outputlo_d[dim][i] - cfflo[i]);

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
            cout << outputhi_h[idx[k]][i] << "  "
                 << outputmi_h[idx[k]][i] << "  "
                 << outputlo_h[idx[k]][i] << endl;
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
               cout << outputhi_d[idx[k]][i] << "  "
                    << outputmi_d[idx[k]][i] << "  "
                    << outputlo_d[idx[k]][i] << endl;
            }
            errsum = errsum
                   + abs(outputhi_h[idx[k]][i] - outputhi_d[idx[k]][i])
                   + abs(outputmi_h[idx[k]][i] - outputmi_d[idx[k]][i])
                   + abs(outputlo_h[idx[k]][i] - outputlo_d[idx[k]][i]);
         }
         if(verbose > 0) cout << "Sum of errors : " << errsum << endl;

         errtot += errsum;
      }
      if(verbose > 0) cout << "Total sum of all errors : " << errtot << endl;
   }
   return errtot;
}

double test_dbl3_complex ( int dim, int nvr, int pwr, int deg, int verbose )
{
   int *idx = new int[nvr];       // indices of variables in the monomial
   int *exp = new int[nvr];       // exponents of the variables
   int *expfac = new int[nvr];    // exponents of common factor
   int nbrfac;                    // number of common factors
   double *cffrehi = new double[deg+1]; // high real parts of coefficients
   double *cffremi = new double[deg+1]; // middle real parts of coefficients
   double *cffrelo = new double[deg+1]; // low real parts of coefficients
   double *cffimhi = new double[deg+1]; // high imaginary coefficients
   double *cffimmi = new double[deg+1]; // middle imaginary coefficients
   double *cffimlo = new double[deg+1]; // low imaginary coefficients

   // The input are dim power series of degree deg,
   // the output are nvr+1 power series of degree deg,
   // for the evaluated and differentiated monomial.

   double **inputrehi = new double*[dim];
   double **inputremi = new double*[dim];
   double **inputrelo = new double*[dim];
   double **inputimhi = new double*[dim];
   double **inputimmi = new double*[dim];
   double **inputimlo = new double*[dim];
   for(int i=0; i<dim; i++)
   {
      inputrehi[i] = new double[deg+1]; inputremi[i] = new double[deg+1];
      inputrelo[i] = new double[deg+1];
      inputimhi[i] = new double[deg+1]; inputimmi[i] = new double[deg+1];
      inputimlo[i] = new double[deg+1];
   }
   double **outputrehi_h = new double*[dim+1];
   double **outputremi_h = new double*[dim+1];
   double **outputrelo_h = new double*[dim+1];
   double **outputimhi_h = new double*[dim+1];
   double **outputimmi_h = new double*[dim+1];
   double **outputimlo_h = new double*[dim+1];
   for(int i=0; i<=dim; i++)
   {
      outputrehi_h[i] = new double[deg+1];
      outputremi_h[i] = new double[deg+1];
      outputrelo_h[i] = new double[deg+1];
      outputimhi_h[i] = new double[deg+1];
      outputimmi_h[i] = new double[deg+1];
      outputimlo_h[i] = new double[deg+1];
   }
   double **outputrehi_d = new double*[dim+1];
   double **outputremi_d = new double*[dim+1];
   double **outputrelo_d = new double*[dim+1];
   double **outputimhi_d = new double*[dim+1];
   double **outputimmi_d = new double*[dim+1];
   double **outputimlo_d = new double*[dim+1];
   for(int i=0; i<=dim; i++)
   {
      outputrehi_d[i] = new double[deg+1];
      outputremi_d[i] = new double[deg+1];
      outputrelo_d[i] = new double[deg+1];
      outputimhi_d[i] = new double[deg+1];
      outputimmi_d[i] = new double[deg+1];
      outputimlo_d[i] = new double[deg+1];
   }
   bool fail = make_complex3_monomial
      (dim,nvr,pwr,deg,idx,exp,
       cffrehi,cffremi,cffrelo,cffimhi,cffimmi,cffimlo);

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
      {
         cout << cffrehi[i] << "  " << cffremi[i] << "  " << cffrelo[i]
              << endl;
         cout << cffimhi[i] << "  " << cffimmi[i] << "  " << cffimlo[i]
              << endl;
      }
   }
   make_complex3_input
      (dim,deg,inputrehi,inputremi,inputrelo,inputimhi,inputimmi,inputimlo);

   if(verbose > 0)
   {
      cout << "Random input series :" << endl;
      for(int i=0; i<dim; i++)
      {
         cout << "-> coefficients of series " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
         {
            cout << inputrehi[i][j] << "  " << inputremi[i][j]
                                    << "  " << inputrelo[i][j] << endl;
            cout << inputimhi[i][j] << "  " << inputimmi[i][j]
                                    << "  " << inputimlo[i][j] << endl;
         }
      }
   }
   double errsum = 0.0;
   double errtot = 0.0;

   CPU_cmplx3_evaldiff(dim,nvr,deg,idx,
      cffrehi,cffremi,cffrelo,cffimhi,cffimmi,cffimlo,
      inputrehi,inputremi,inputrelo,inputimhi,inputimmi,inputimlo,
      outputrehi_h,outputremi_h,outputrelo_h,
      outputimhi_h,outputimmi_h,outputimlo_h);

   if(verbose > 0)
   {
      cout << "The value of the product :" << endl;
      for(int i=0; i<=deg; i++)
      {
         cout << outputrehi_h[dim][i] << "  "
              << outputremi_h[dim][i] << "  "
              << outputrelo_h[dim][i] << endl;
         cout << outputimhi_h[dim][i] << "  "
              << outputimmi_h[dim][i] << "  "
              << outputimlo_h[dim][i] << endl;
      }
   }
   if(nvr == dim) // the product of all input series equals one
   {
      for(int i=0; i<=deg; i++)
         errsum = errsum
                + abs(outputrehi_h[dim][i] - cffrehi[i])
                + abs(outputremi_h[dim][i] - cffremi[i])
                + abs(outputrelo_h[dim][i] - cffrelo[i])
                + abs(outputimhi_h[dim][i] - cffimhi[i])
                + abs(outputimmi_h[dim][i] - cffimmi[i])
                + abs(outputimlo_h[dim][i] - cffimlo[i]);

      if(verbose > 0) cout << "Coefficient error : " << errsum << endl;

      errtot += errsum;
   }
   if(nvr > 2)
   {
      GPU_cmplx3_evaldiff(deg+1,dim,nvr,deg,idx,
         cffrehi,cffremi,cffrelo,cffimhi,cffimmi,cffimlo,
         inputrehi,inputremi,inputrelo,inputimhi,inputimmi,inputimlo,
         outputrehi_d,outputremi_d,outputrelo_d,
         outputimhi_d,outputimmi_d,outputimlo_d);

      errsum = 0.0;

      if(verbose > 0)
         cout << "The value of the product computed on the GPU :" << endl;

      for(int i=0; i<=deg; i++) 
      {
         if(verbose > 0)
         {
            cout << outputrehi_d[dim][i] << "  "
                 << outputremi_d[dim][i] << "  "
                 << outputrelo_d[dim][i] << endl;
            cout << outputimhi_d[dim][i] << "  "
                 << outputimmi_d[dim][i] << "  "
                 << outputimlo_d[dim][i] << endl;
         }
         errsum = errsum
                + abs(outputrehi_h[dim][i] - outputrehi_d[dim][i])
                + abs(outputremi_h[dim][i] - outputremi_d[dim][i])
                + abs(outputrelo_h[dim][i] - outputrelo_d[dim][i])
                + abs(outputimhi_h[dim][i] - outputimhi_d[dim][i])
                + abs(outputimmi_h[dim][i] - outputimmi_d[dim][i])
                + abs(outputimlo_h[dim][i] - outputimlo_d[dim][i]);
      }
      if(verbose > 0) cout << "The sum of errors : " << errsum << endl;

      errtot += errsum;

      if(nvr == dim) // the product of all input series equals one
      {
         errsum = 0.0;
         for(int i=0; i<=deg; i++)
            errsum = errsum
                   + abs(outputrehi_d[dim][i] - cffrehi[i])
                   + abs(outputremi_d[dim][i] - cffremi[i])
                   + abs(outputrelo_d[dim][i] - cffrelo[i])
                   + abs(outputimhi_d[dim][i] - cffimhi[i])
                   + abs(outputimmi_d[dim][i] - cffimmi[i])
                   + abs(outputimlo_d[dim][i] - cffimlo[i]);

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
         {
            cout << outputrehi_h[idx[k]][i] << "  "
                 << outputremi_h[idx[k]][i] << "  "
                 << outputrelo_h[idx[k]][i] << endl;
            cout << outputimhi_h[idx[k]][i] << "  "
                 << outputimmi_h[idx[k]][i] << "  "
                 << outputimlo_h[idx[k]][i] << endl;
         }
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
               cout << outputrehi_d[idx[k]][i] << "  "
                    << outputremi_d[idx[k]][i] << "  "
                    << outputrelo_d[idx[k]][i] << endl;
               cout << outputimhi_d[idx[k]][i] << "  "
                    << outputimmi_d[idx[k]][i] << "  "
                    << outputimlo_d[idx[k]][i] << endl;
            }
            errsum = errsum
                   + abs(outputrehi_h[idx[k]][i] - outputrehi_d[idx[k]][i])
                   + abs(outputremi_h[idx[k]][i] - outputremi_d[idx[k]][i])
                   + abs(outputrelo_h[idx[k]][i] - outputrelo_d[idx[k]][i])
                   + abs(outputimhi_h[idx[k]][i] - outputimhi_d[idx[k]][i])
                   + abs(outputimmi_h[idx[k]][i] - outputimmi_d[idx[k]][i])
                   + abs(outputimlo_h[idx[k]][i] - outputimlo_d[idx[k]][i]);
         }
         if(verbose > 0) cout << "The sum of errors : " << errsum << endl;

         errtot += errsum;
      }
   }
   if(verbose > 0) cout << "Total sum of all errors : " << errtot << endl;

   return errtot;
}
