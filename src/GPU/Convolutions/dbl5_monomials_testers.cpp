/* The file dbl5_monomials_testers.cpp contains the definitions of
 * functions with prototypes in dbl5_monomials_testers.h. */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <vector_types.h>
#include "random_monomials.h"
#include "random5_monomials.h"
#include "dbl5_monomials_host.h"
#include "dbl5_monomials_kernels.h"
#include "dbl5_monomials_testers.h"

using namespace std;

int main_dbl5_test
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
   double realsum = test_dbl5_real(dim,nvr,pwr,deg,vrblvl-1);
   double complexsum = test_dbl5_complex(dim,nvr,pwr,deg,vrblvl-1);

   const double tol = 1.0e-76;

   int fail = int(realsum > tol) + int(complexsum > tol);

   if(vrblvl > 0)
   {
      cout << "Sum of all errors in penta double precision :" << endl;
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

      cout << "  Seed used : " << seedused << endl;
   }
   return fail;
}

double test_dbl5_real ( int dim, int nvr, int pwr, int deg, int verbose )
{
   int *idx = new int[nvr];           // indices of variables in the monomial
   int *exp = new int[nvr];           // exponents of the variables
   int *expfac = new int[nvr];        // exponents of common factor
   int nbrfac;                        // number of common factors
   double *cfftb = new double[deg+1]; // highest doubles of series coefficient
   double *cffix = new double[deg+1]; // second highest doubles
   double *cffmi = new double[deg+1]; // middle doubles of series coefficient
   double *cffrg = new double[deg+1]; // second lowest doubles
   double *cffpk = new double[deg+1]; // lowest doubles of series coefficient

   // The input are dim power series of degree deg,
   // the output are nvr+1 power series of degree deg,
   // for the evaluated and differentiated monomial.

   double **inputtb = new double*[dim];
   for(int i=0; i<dim; i++) inputtb[i] = new double[deg+1];
   double **inputix = new double*[dim];
   for(int i=0; i<dim; i++) inputix[i] = new double[deg+1];
   double **inputmi = new double*[dim];
   for(int i=0; i<dim; i++) inputmi[i] = new double[deg+1];
   double **inputrg = new double*[dim];
   for(int i=0; i<dim; i++) inputrg[i] = new double[deg+1];
   double **inputpk = new double*[dim];
   for(int i=0; i<dim; i++) inputpk[i] = new double[deg+1];
   double **outputtb_h = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputtb_h[i] = new double[deg+1];
   double **outputtb_d = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputtb_d[i] = new double[deg+1];
   double **outputix_h = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputix_h[i] = new double[deg+1];
   double **outputix_d = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputix_d[i] = new double[deg+1];
   double **outputmi_h = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputmi_h[i] = new double[deg+1];
   double **outputmi_d = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputmi_d[i] = new double[deg+1];
   double **outputrg_h = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputrg_h[i] = new double[deg+1];
   double **outputrg_d = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputrg_d[i] = new double[deg+1];
   double **outputpk_h = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputpk_h[i] = new double[deg+1];
   double **outputpk_d = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputpk_d[i] = new double[deg+1];

   srand(time(NULL));

   bool fail = make_real5_monomial
      (dim,nvr,pwr,deg,idx,exp,cfftb,cffix,cffmi,cffrg,cffpk);

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
         cout << cfftb[i] << "  " << cffix[i] << "  " << cffmi[i] << endl
              << cffrg[i] << "  " << cffpk[i] << endl;
   }
   make_real5_input(dim,deg,inputtb,inputix,inputmi,inputrg,inputpk);

   if(verbose > 0)
   {
      cout << "Random input series :" << endl;
      for(int i=0; i<dim; i++)
      {
         cout << "-> coefficients of series " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
            cout << inputtb[i][j] << "  " << inputix[i][j]
                                  << "  " << inputmi[i][j] << endl
                 << inputrg[i][j] << "  " << inputpk[i][j] << endl;
      }
   }
   double errsum = 0.0;
   double errtot = 0.0;

   CPU_dbl5_evaldiff(dim,nvr,deg,idx,cfftb,cffix,cffmi,cffrg,cffpk,
                     inputtb,inputix,inputmi,inputrg,inputpk,
                     outputtb_h,outputix_h,outputmi_h,outputrg_h,outputpk_h);
   if(verbose > 0)
   {
      cout << "The value of the product :" << endl;
      for(int i=0; i<=deg; i++)
         cout << outputtb_h[dim][i] << "  " << outputix_h[dim][i] 
                                    << "  " << outputmi_h[dim][i] << endl
              << outputrg_h[dim][i] << "  " << outputpk_h[dim][i] << endl;
   }
   if(nvr == dim) // the product of all input series equals one
   {
      for(int i=0; i<=deg; i++)
         errsum = errsum
                + abs(outputtb_h[dim][i] - cfftb[i])
                + abs(outputix_h[dim][i] - cffix[i])
                + abs(outputmi_h[dim][i] - cffmi[i])
                + abs(outputrg_h[dim][i] - cffrg[i])
                + abs(outputpk_h[dim][i] - cffpk[i]);

      if(verbose > 0) cout << "Coefficient error : " << errsum << endl;

      errtot += errsum;
   }
   if(nvr > 2)
   {
      GPU_dbl5_evaldiff(deg+1,dim,nvr,deg,idx,cfftb,cffix,cffmi,cffrg,cffpk,
         inputtb,inputix,inputmi,inputrg,inputpk,
         outputtb_d,outputix_d,outputmi_d,outputrg_d,outputpk_d);

      errsum = 0.0;

      if(verbose > 0)
         cout << "The value of the product computed on the GPU :" << endl;

      for(int i=0; i<=deg; i++)
      {
         if(verbose > 0)
         {
            cout << outputtb_d[dim][i] << "  "
                 << outputix_d[dim][i] << "  "
                 << outputmi_d[dim][i] << endl
                 << outputrg_d[dim][i] << "  "
                 << outputpk_d[dim][i] << endl;
         }
         errsum = errsum
                + abs(outputtb_h[dim][i] - outputtb_d[dim][i])
                + abs(outputix_h[dim][i] - outputix_d[dim][i])
                + abs(outputmi_h[dim][i] - outputmi_d[dim][i])
                + abs(outputrg_h[dim][i] - outputrg_d[dim][i])
                + abs(outputpk_h[dim][i] - outputpk_d[dim][i]);
      }
      if(verbose > 0) cout << "Sum of errors : " << errsum << endl;

      errtot += errsum;

      if(nvr == dim) // the product of all input series equals one
      {
         errsum = 0.0;
         for(int i=0; i<=deg; i++)
            errsum = errsum
                   + abs(outputtb_d[dim][i] - cfftb[i])
                   + abs(outputix_d[dim][i] - cffix[i])
                   + abs(outputmi_d[dim][i] - cffmi[i])
                   + abs(outputrg_d[dim][i] - cffrg[i])
                   + abs(outputpk_d[dim][i] - cffpk[i]);

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
            cout << outputtb_h[idx[k]][i] << "  "
                 << outputix_h[idx[k]][i] << "  "
                 << outputmi_h[idx[k]][i] << endl
                 << outputrg_h[idx[k]][i] << "  "
                 << outputpk_h[idx[k]][i] << endl;
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
               cout << outputtb_d[idx[k]][i] << "  "
                    << outputix_d[idx[k]][i] << "  "
                    << outputmi_d[idx[k]][i] << endl
                    << outputrg_d[idx[k]][i] << "  "
                    << outputpk_d[idx[k]][i] << endl;
            }
            errsum = errsum
                   + abs(outputtb_h[idx[k]][i] - outputtb_d[idx[k]][i])
                   + abs(outputix_h[idx[k]][i] - outputix_d[idx[k]][i])
                   + abs(outputmi_h[idx[k]][i] - outputmi_d[idx[k]][i])
                   + abs(outputrg_h[idx[k]][i] - outputrg_d[idx[k]][i])
                   + abs(outputpk_h[idx[k]][i] - outputpk_d[idx[k]][i]);
         }
         if(verbose > 0) cout << "Sum of errors : " << errsum << endl;

         errtot += errsum;
      }
      if(verbose > 0) cout << "Total sum of all errors : " << errtot << endl;
   }
   return errtot;
}

double test_dbl5_complex ( int dim, int nvr, int pwr, int deg, int verbose )
{
   int *idx = new int[nvr];       // indices of variables in the monomial
   int *exp = new int[nvr];       // exponents of the variables
   int *expfac = new int[nvr];    // exponents of common factor
   int nbrfac;                    // number of common factors
   double *cffretb = new double[deg+1]; // highest real parts of coefficients
   double *cffreix = new double[deg+1]; // second highest real parts
   double *cffremi = new double[deg+1]; // middle real parts
   double *cffrerg = new double[deg+1]; // second lowest real parts
   double *cffrepk = new double[deg+1]; // lowest real parts
   double *cffimtb = new double[deg+1]; // highest imaginary parts
   double *cffimix = new double[deg+1]; // second highest imaginary parts
   double *cffimmi = new double[deg+1]; // middle imaginary parts
   double *cffimrg = new double[deg+1]; // second lowest imaginary parts
   double *cffimpk = new double[deg+1]; // lowest imaginary parts

   // The input are dim power series of degree deg,
   // the output are nvr+1 power series of degree deg,
   // for the evaluated and differentiated monomial.

   double **inputretb = new double*[dim];
   double **inputreix = new double*[dim];
   double **inputremi = new double*[dim];
   double **inputrerg = new double*[dim];
   double **inputrepk = new double*[dim];
   double **inputimtb = new double*[dim];
   double **inputimix = new double*[dim];
   double **inputimmi = new double*[dim];
   double **inputimrg = new double*[dim];
   double **inputimpk = new double*[dim];
   for(int i=0; i<dim; i++)
   {
      inputretb[i] = new double[deg+1]; inputreix[i] = new double[deg+1];
      inputremi[i] = new double[deg+1]; inputrerg[i] = new double[deg+1];
      inputrepk[i] = new double[deg+1];
      inputimtb[i] = new double[deg+1]; inputimix[i] = new double[deg+1];
      inputimmi[i] = new double[deg+1]; inputimrg[i] = new double[deg+1];
      inputimpk[i] = new double[deg+1];
   }
   double **outputretb_h = new double*[dim+1];
   double **outputreix_h = new double*[dim+1];
   double **outputremi_h = new double*[dim+1];
   double **outputrerg_h = new double*[dim+1];
   double **outputrepk_h = new double*[dim+1];
   double **outputimtb_h = new double*[dim+1];
   double **outputimix_h = new double*[dim+1];
   double **outputimmi_h = new double*[dim+1];
   double **outputimrg_h = new double*[dim+1];
   double **outputimpk_h = new double*[dim+1];
   for(int i=0; i<=dim; i++)
   {
      outputretb_h[i] = new double[deg+1];
      outputreix_h[i] = new double[deg+1];
      outputremi_h[i] = new double[deg+1];
      outputrerg_h[i] = new double[deg+1];
      outputrepk_h[i] = new double[deg+1];
      outputimtb_h[i] = new double[deg+1];
      outputimix_h[i] = new double[deg+1];
      outputimmi_h[i] = new double[deg+1];
      outputimrg_h[i] = new double[deg+1];
      outputimpk_h[i] = new double[deg+1];
   }
   double **outputretb_d = new double*[dim+1];
   double **outputreix_d = new double*[dim+1];
   double **outputremi_d = new double*[dim+1];
   double **outputrerg_d = new double*[dim+1];
   double **outputrepk_d = new double*[dim+1];
   double **outputimtb_d = new double*[dim+1];
   double **outputimix_d = new double*[dim+1];
   double **outputimmi_d = new double*[dim+1];
   double **outputimrg_d = new double*[dim+1];
   double **outputimpk_d = new double*[dim+1];
   for(int i=0; i<=dim; i++)
   {
      outputretb_d[i] = new double[deg+1];
      outputreix_d[i] = new double[deg+1];
      outputremi_d[i] = new double[deg+1];
      outputrerg_d[i] = new double[deg+1];
      outputrepk_d[i] = new double[deg+1];
      outputimtb_d[i] = new double[deg+1];
      outputimix_d[i] = new double[deg+1];
      outputimmi_d[i] = new double[deg+1];
      outputimrg_d[i] = new double[deg+1];
      outputimpk_d[i] = new double[deg+1];
   }

   bool fail = make_complex5_monomial
      (dim,nvr,pwr,deg,idx,exp,
       cffretb,cffreix,cffremi,cffrerg,cffrepk,
       cffimtb,cffimix,cffimmi,cffimrg,cffimpk);

   if(fail) return 1;

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
         cout << cffretb[i] << "  " << cffreix[i] << "  " << cffremi[i]
              << cffrerg[i] << "  " << cffrepk[i] << endl;
         cout << cffimtb[i] << "  " << cffimix[i] << "  " << cffimmi[i]
              << cffimrg[i] << "  " << cffimpk[i] << endl;
      }
   }
   make_complex5_input
      (dim,deg,inputretb,inputreix,inputremi,inputrerg,inputrepk,
               inputimtb,inputimix,inputimmi,inputimrg,inputimpk);

   if(verbose > 0)
   {
      cout << "Random input series :" << endl;
      for(int i=0; i<dim; i++)
      {
         cout << "-> coefficients of series " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
         {
            cout << inputretb[i][j] << "  " << inputreix[i][j]
                                    << "  " << inputremi[i][j] << endl
                 << inputrerg[i][j] << "  " << inputrepk[i][j] << endl;
            cout << inputimtb[i][j] << "  " << inputimix[i][j]
                                    << "  " << inputimmi[i][j] << endl
                 << inputimrg[i][j] << "  " << inputimpk[i][j] << endl;
         }
      }
   }
   double errsum = 0.0;
   double errtot = 0.0;

   CPU_cmplx5_evaldiff(dim,nvr,deg,idx,
      cffretb,cffreix,cffremi,cffrerg,cffrepk,
      cffimtb,cffimix,cffimmi,cffimrg,cffimpk,
      inputretb,inputreix,inputremi,inputrerg,inputrepk,
      inputimtb,inputimix,inputimmi,inputimrg,inputimpk,
      outputretb_h,outputreix_h,outputremi_h,outputrerg_h,outputrepk_h,
      outputimtb_h,outputimix_h,outputimmi_h,outputimrg_h,outputimpk_h);

   if(verbose > 0)
   {
      cout << "The value of the product :" << endl;
      for(int i=0; i<=deg; i++)
      {
         cout << outputretb_h[dim][i] << "  " << outputreix_h[dim][i]
                                      << "  " << outputremi_h[dim][i] << endl
              << outputrerg_h[dim][i] << "  " << outputrepk_h[dim][i] << endl;
         cout << outputimtb_h[dim][i] << "  " << outputimix_h[dim][i]
                                      << "  " << outputimmi_h[dim][i] << endl
              << outputimrg_h[dim][i] << "  " << outputimpk_h[dim][i] << endl;
      }
   }
   if(nvr == dim) // the product of all input series equals one
   {
      for(int i=0; i<=deg; i++)
         errsum = errsum
                + abs(outputretb_h[dim][i] - cffretb[i])
                + abs(outputreix_h[dim][i] - cffreix[i])
                + abs(outputremi_h[dim][i] - cffremi[i])
                + abs(outputrerg_h[dim][i] - cffrerg[i])
                + abs(outputrepk_h[dim][i] - cffrepk[i])
                + abs(outputimtb_h[dim][i] - cffimtb[i])
                + abs(outputimix_h[dim][i] - cffimix[i])
                + abs(outputimmi_h[dim][i] - cffimmi[i])
                + abs(outputimrg_h[dim][i] - cffimrg[i])
                + abs(outputimpk_h[dim][i] - cffimpk[i]);

      if(verbose > 0) cout << "Coefficient error : " << errsum << endl;

      errtot += errsum;
   }
   if(nvr > 2)
   {
      GPU_cmplx5_evaldiff(deg+1,dim,nvr,deg,idx,
         cffretb,cffreix,cffremi,cffrerg,cffrepk,
         cffimtb,cffimix,cffimmi,cffimrg,cffimpk,
         inputretb,inputreix,inputremi,inputrerg,inputrepk,
         inputimtb,inputimix,inputimmi,inputimrg,inputimpk,
         outputretb_d,outputreix_d,outputremi_d,outputrerg_d,outputrepk_d,
         outputimtb_d,outputimix_d,outputimmi_d,outputimrg_d,outputimpk_d);

      errsum = 0.0;

      if(verbose > 0)
         cout << "The value of the product computed on the GPU :" << endl;

      for(int i=0; i<=deg; i++) 
      {
         if(verbose > 0)
         {
            cout << outputretb_d[dim][i] << "  "
                 << outputreix_d[dim][i] << "  "
                 << outputremi_d[dim][i] << endl
                 << outputrerg_d[dim][i] << "  "
                 << outputrepk_d[dim][i] << endl;
            cout << outputimtb_d[dim][i] << "  "
                 << outputimix_d[dim][i] << "  "
                 << outputimmi_d[dim][i] << endl
                 << outputimrg_d[dim][i] << "  "
                 << outputimpk_d[dim][i] << endl;
         }
         errsum = errsum
                + abs(outputretb_h[dim][i] - outputretb_d[dim][i])
                + abs(outputreix_h[dim][i] - outputreix_d[dim][i])
                + abs(outputremi_h[dim][i] - outputremi_d[dim][i])
                + abs(outputrerg_h[dim][i] - outputrerg_d[dim][i])
                + abs(outputrepk_h[dim][i] - outputrepk_d[dim][i])
                + abs(outputimtb_h[dim][i] - outputimtb_d[dim][i])
                + abs(outputimix_h[dim][i] - outputimix_d[dim][i])
                + abs(outputimmi_h[dim][i] - outputimmi_d[dim][i])
                + abs(outputimrg_h[dim][i] - outputimrg_d[dim][i])
                + abs(outputimpk_h[dim][i] - outputimpk_d[dim][i]);
      }
      if(verbose > 0) cout << "The sum of errors : " << errsum << endl;

      errtot += errsum;

      if(nvr == dim) // the product of all input series equals one
      {
         errsum = 0.0;
         for(int i=0; i<=deg; i++)
            errsum = errsum
                   + abs(outputretb_d[dim][i] - cffretb[i])
                   + abs(outputreix_d[dim][i] - cffreix[i])
                   + abs(outputremi_d[dim][i] - cffremi[i])
                   + abs(outputrerg_d[dim][i] - cffrerg[i])
                   + abs(outputrepk_d[dim][i] - cffrepk[i])
                   + abs(outputimtb_d[dim][i] - cffimtb[i])
                   + abs(outputimix_d[dim][i] - cffimix[i])
                   + abs(outputimmi_d[dim][i] - cffimmi[i])
                   + abs(outputimrg_d[dim][i] - cffimrg[i])
                   + abs(outputimpk_d[dim][i] - cffimpk[i]);

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
            cout << outputretb_h[idx[k]][i] << "  "
                 << outputreix_h[idx[k]][i] << "  "
                 << outputremi_h[idx[k]][i] << "  "
                 << outputrerg_h[idx[k]][i] << "  "
                 << outputrepk_h[idx[k]][i] << endl;
            cout << outputimtb_h[idx[k]][i] << "  "
                 << outputimix_h[idx[k]][i] << "  "
                 << outputimmi_h[idx[k]][i] << "  "
                 << outputimrg_h[idx[k]][i] << "  "
                 << outputimpk_h[idx[k]][i] << endl;
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
               cout << outputretb_d[idx[k]][i] << "  "
                    << outputreix_d[idx[k]][i] << "  "
                    << outputremi_d[idx[k]][i] << "  "
                    << outputrerg_d[idx[k]][i] << "  "
                    << outputrepk_d[idx[k]][i] << endl;
               cout << outputimtb_d[idx[k]][i] << "  "
                    << outputimix_d[idx[k]][i] << "  "
                    << outputimmi_d[idx[k]][i] << "  "
                    << outputimrg_d[idx[k]][i] << "  "
                    << outputimpk_d[idx[k]][i] << endl;
            }
            errsum = errsum
                   + abs(outputretb_h[idx[k]][i] - outputretb_d[idx[k]][i])
                   + abs(outputreix_h[idx[k]][i] - outputreix_d[idx[k]][i])
                   + abs(outputremi_h[idx[k]][i] - outputremi_d[idx[k]][i])
                   + abs(outputrerg_h[idx[k]][i] - outputrerg_d[idx[k]][i])
                   + abs(outputrepk_h[idx[k]][i] - outputrepk_d[idx[k]][i])
                   + abs(outputimtb_h[idx[k]][i] - outputimtb_d[idx[k]][i])
                   + abs(outputimix_h[idx[k]][i] - outputimix_d[idx[k]][i])
                   + abs(outputimmi_h[idx[k]][i] - outputimmi_d[idx[k]][i])
                   + abs(outputimrg_h[idx[k]][i] - outputimrg_d[idx[k]][i])
                   + abs(outputimpk_h[idx[k]][i] - outputimpk_d[idx[k]][i]);
         }
         if(verbose > 0) cout << "The sum of errors : " << errsum << endl;

         errtot += errsum;
      }
   }
   if(verbose > 0) cout << "Total sum of all errors : " << errtot << endl;

   return errtot;
}
