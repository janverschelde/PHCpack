// Tests the evaluation and differentiation of a monomial
// in deca double precision.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <vector_types.h>
#include "random_monomials.h"
#include "random10_monomials.h"
#include "dbl10_monomials_host.h"
#include "dbl10_monomials_kernels.h"

using namespace std;

double test_real ( int dim, int nvr, int pwr, int deg );
/*
 * DESCRIPTION :
 *   Tests the evaluation and differentiation for random real data.
 *   Returns the sum of all errors.
 * 
 * ON ENTRY :
 *   dim      dimension, total number of variables;
 *   nvr      number of variables in the product;
 *   pwr      highest power of each variable;
 *   deg      truncation degree of the series. */

// double test_complex ( int dim, int nvr, int pwr, int deg );
/*
 * DESCRIPTION :
 *   Tests the evaluation and differentiation for random complex data.
 *   Returns the sum of all errors.
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
   double realsum = test_real(dim,nvr,pwr,deg);
   // cout << endl << "Testing for complex input data ..." << endl;
   // double complexsum = test_complex(dim,nvr,pwr,deg);

   const double tol = 1.0e-152;

   cout << endl << "Sum of all errors :" << endl;
   cout << "  on real data : " << realsum;
   if(realsum < tol)
      cout << "  pass," << endl;
   else
      cout << "  fail!" << endl;
/*
   cout << "  on complex data : " << complexsum;
   if(complexsum < tol)
      cout << "  pass." << endl;
   else
      cout << "  fail!" << endl;
 */
   return 0;
}

double test_real ( int dim, int nvr, int pwr, int deg )
{
   int *idx = new int[nvr];           // indices of variables in the monomial
   int *exp = new int[nvr];           // exponents of the variables
   int *expfac = new int[nvr];        // exponents of common factor
   int nbrfac;                        // number of common factors
   double *cffrtb = new double[deg+1]; // highest coefficient doubles
   double *cffrix = new double[deg+1]; // second highest doubles
   double *cffrmi = new double[deg+1]; // third highest doubles
   double *cffrrg = new double[deg+1]; // fourth highest doubles
   double *cffrpk = new double[deg+1]; // fifth highest doubles
   double *cffltb = new double[deg+1]; // fifth lowest doubles
   double *cfflix = new double[deg+1]; // fourth lowest doubles
   double *cfflmi = new double[deg+1]; // third lowest doubles
   double *cfflrg = new double[deg+1]; // second lowest doubles
   double *cfflpk = new double[deg+1]; // lowest doubles

   // The input are dim power series of degree deg,
   // the output are nvr+1 power series of degree deg,
   // for the evaluated and differentiated monomial.

   double **inputrtb = new double*[dim];
   for(int i=0; i<dim; i++) inputrtb[i] = new double[deg+1];
   double **inputrix = new double*[dim];
   for(int i=0; i<dim; i++) inputrix[i] = new double[deg+1];
   double **inputrmi = new double*[dim];
   for(int i=0; i<dim; i++) inputrmi[i] = new double[deg+1];
   double **inputrrg = new double*[dim];
   for(int i=0; i<dim; i++) inputrrg[i] = new double[deg+1];
   double **inputrpk = new double*[dim];
   for(int i=0; i<dim; i++) inputrpk[i] = new double[deg+1];
   double **inputltb = new double*[dim];
   for(int i=0; i<dim; i++) inputltb[i] = new double[deg+1];
   double **inputlix = new double*[dim];
   for(int i=0; i<dim; i++) inputlix[i] = new double[deg+1];
   double **inputlmi = new double*[dim];
   for(int i=0; i<dim; i++) inputlmi[i] = new double[deg+1];
   double **inputlrg = new double*[dim];
   for(int i=0; i<dim; i++) inputlrg[i] = new double[deg+1];
   double **inputlpk = new double*[dim];
   for(int i=0; i<dim; i++) inputlpk[i] = new double[deg+1];
   double **outputrtb_h = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputrtb_h[i] = new double[deg+1];
   double **outputrtb_d = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputrtb_d[i] = new double[deg+1];
   double **outputrix_h = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputrix_h[i] = new double[deg+1];
   double **outputrix_d = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputrix_d[i] = new double[deg+1];
   double **outputrmi_h = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputrmi_h[i] = new double[deg+1];
   double **outputrmi_d = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputrmi_d[i] = new double[deg+1];
   double **outputrrg_h = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputrrg_h[i] = new double[deg+1];
   double **outputrrg_d = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputrrg_d[i] = new double[deg+1];
   double **outputrpk_h = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputrpk_h[i] = new double[deg+1];
   double **outputrpk_d = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputrpk_d[i] = new double[deg+1];
   double **outputltb_h = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputltb_h[i] = new double[deg+1];
   double **outputltb_d = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputltb_d[i] = new double[deg+1];
   double **outputlix_h = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputlix_h[i] = new double[deg+1];
   double **outputlix_d = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputlix_d[i] = new double[deg+1];
   double **outputlmi_h = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputlmi_h[i] = new double[deg+1];
   double **outputlmi_d = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputlmi_d[i] = new double[deg+1];
   double **outputlrg_h = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputlrg_h[i] = new double[deg+1];
   double **outputlrg_d = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputlrg_d[i] = new double[deg+1];
   double **outputlpk_h = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputlpk_h[i] = new double[deg+1];
   double **outputlpk_d = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputlpk_d[i] = new double[deg+1];

   srand(time(NULL));

   bool fail = make_real10_monomial
      (dim,nvr,pwr,deg,idx,exp,
       cffrtb,cffrix,cffrmi,cffrrg,cffrpk,
       cffltb,cfflix,cfflmi,cfflrg,cfflpk);

   if(fail) return 1.0;

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
      cout << cffrtb[i] << "  " << cffrix[i] << "  " << cffrmi[i] << endl
           << cffrrg[i] << "  " << cffrpk[i] << endl;
      cout << cffltb[i] << "  " << cfflix[i] << "  " << cfflmi[i] << endl
           << cfflrg[i] << "  " << cfflpk[i] << endl;
   }

   make_real10_input(dim,deg,
      inputrtb,inputrix,inputrmi,inputrrg,inputrpk,
      inputltb,inputlix,inputlmi,inputlrg,inputlpk);

   cout << "Random input series :" << endl;
   for(int i=0; i<dim; i++)
   {
      cout << "-> coefficients of series " << i << " :" << endl;
      for(int j=0; j<=deg; j++)
      {
         cout << inputrtb[i][j] << "  " << inputrix[i][j]
                                << "  " << inputrmi[i][j] << endl
              << inputrrg[i][j] << "  " << inputrpk[i][j] << endl;
         cout << inputltb[i][j] << "  " << inputlix[i][j]
                                << "  " << inputlmi[i][j] << endl
              << inputlrg[i][j] << "  " << inputlpk[i][j] << endl;
      }
   }
   double errsum = 0.0;
   double errtot = 0.0;

   CPU_dbl10_evaldiff(dim,nvr,deg,idx,
      cffrtb,cffrix,cffrmi,cffrrg,cffrpk,
      cffltb,cfflix,cfflmi,cfflrg,cfflpk,
      inputrtb,inputrix,inputrmi,inputrrg,inputrpk,
      inputltb,inputlix,inputlmi,inputlrg,inputlpk,
      outputrtb_h,outputrix_h,outputrmi_h,outputrrg_h,outputrpk_h,
      outputltb_h,outputlix_h,outputlmi_h,outputlrg_h,outputlpk_h);

   cout << "The value of the product :" << endl;
   for(int i=0; i<=deg; i++)
   {
      cout << outputrtb_h[dim][i] << "  " << outputrix_h[dim][i] 
                                  << "  " << outputrmi_h[dim][i] << endl
           << outputrrg_h[dim][i] << "  " << outputrpk_h[dim][i] << endl;
      cout << outputltb_h[dim][i] << "  " << outputlix_h[dim][i] 
                                  << "  " << outputlmi_h[dim][i] << endl
           << outputlrg_h[dim][i] << "  " << outputlpk_h[dim][i] << endl;
   }
   if(nvr == dim) // the product of all input series equals one
   {
      for(int i=0; i<=deg; i++)
         errsum = errsum
                + abs(outputrtb_h[dim][i] - cffrtb[i])
                + abs(outputrix_h[dim][i] - cffrix[i])
                + abs(outputrmi_h[dim][i] - cffrmi[i])
                + abs(outputrrg_h[dim][i] - cffrrg[i])
                + abs(outputrpk_h[dim][i] - cffrpk[i])
                + abs(outputltb_h[dim][i] - cffltb[i])
                + abs(outputlix_h[dim][i] - cfflix[i])
                + abs(outputlmi_h[dim][i] - cfflmi[i])
                + abs(outputlrg_h[dim][i] - cfflrg[i])
                + abs(outputlpk_h[dim][i] - cfflpk[i]);
      cout << "Coefficient error : " << errsum << endl; errtot += errsum;
   }

   if(nvr > 2)
   {
      GPU_dbl10_evaldiff(deg+1,dim,nvr,deg,idx,
         cffrtb,cffrix,cffrmi,cffrrg,cffrpk,
         cffltb,cfflix,cfflmi,cfflrg,cfflpk,
         inputrtb,inputrix,inputrmi,inputrrg,inputrpk,
         inputltb,inputlix,inputlmi,inputlrg,inputlpk,
         outputrtb_d,outputrix_d,outputrmi_d,outputrrg_d,outputrpk_d,
         outputltb_d,outputlix_d,outputlmi_d,outputlrg_d,outputlpk_d);

      cout << "The value of the product computed on the GPU :" << endl;
      errsum = 0.0;
      for(int i=0; i<=deg; i++)
      {
         cout << outputrtb_d[dim][i] << "  "
              << outputrix_d[dim][i] << "  "
              << outputrmi_d[dim][i] << endl
              << outputrrg_d[dim][i] << "  "
              << outputrpk_d[dim][i] << endl;
         cout << outputltb_d[dim][i] << "  "
              << outputlix_d[dim][i] << "  "
              << outputlmi_d[dim][i] << endl
              << outputlrg_d[dim][i] << "  "
              << outputlpk_d[dim][i] << endl;
         errsum = errsum
                + abs(outputrtb_h[dim][i] - outputrtb_d[dim][i])
                + abs(outputrix_h[dim][i] - outputrix_d[dim][i])
                + abs(outputrmi_h[dim][i] - outputrmi_d[dim][i])
                + abs(outputrrg_h[dim][i] - outputrrg_d[dim][i])
                + abs(outputrpk_h[dim][i] - outputrpk_d[dim][i])
                + abs(outputltb_h[dim][i] - outputltb_d[dim][i])
                + abs(outputlix_h[dim][i] - outputlix_d[dim][i])
                + abs(outputlmi_h[dim][i] - outputlmi_d[dim][i])
                + abs(outputlrg_h[dim][i] - outputlrg_d[dim][i])
                + abs(outputlpk_h[dim][i] - outputlpk_d[dim][i]);
      }
      cout << "Sum of errors : " << errsum << endl; errtot += errsum;

      if(nvr == dim) // the product of all input series equals one
      {
         errsum = 0.0;
         for(int i=0; i<=deg; i++)
            errsum = errsum
                   + abs(outputrtb_d[dim][i] - cffrtb[i])
                   + abs(outputrix_d[dim][i] - cffrix[i])
                   + abs(outputrmi_d[dim][i] - cffrmi[i])
                   + abs(outputrrg_d[dim][i] - cffrrg[i])
                   + abs(outputrpk_d[dim][i] - cffrpk[i])
                   + abs(outputltb_d[dim][i] - cffltb[i])
                   + abs(outputlix_d[dim][i] - cfflix[i])
                   + abs(outputlmi_d[dim][i] - cfflmi[i])
                   + abs(outputlrg_d[dim][i] - cfflrg[i])
                   + abs(outputlpk_d[dim][i] - cfflpk[i]);
         cout << "Coefficient error : " << errsum << endl; errtot += errsum;
      }
   }
   for(int k=0; k<nvr; k++)
   {
      cout << "-> derivative for index " << idx[k] << " :" << endl;
      for(int i=0; i<=deg; i++)
      {
         cout << outputrtb_h[idx[k]][i] << "  "
              << outputrix_h[idx[k]][i] << "  "
              << outputrmi_h[idx[k]][i] << endl
              << outputrrg_h[idx[k]][i] << "  "
              << outputrpk_h[idx[k]][i] << endl;
         cout << outputltb_h[idx[k]][i] << "  "
              << outputlix_h[idx[k]][i] << "  "
              << outputlmi_h[idx[k]][i] << endl
              << outputlrg_h[idx[k]][i] << "  "
              << outputlpk_h[idx[k]][i] << endl;
      }
      if(nvr > 2)
      {
         cout << "-> derivative for index " << idx[k]
              << " computed on GPU :" << endl;
         errsum = 0.0;
         for(int i=0; i<=deg; i++)
         {
            cout << outputrtb_d[idx[k]][i] << "  "
                 << outputrix_d[idx[k]][i] << "  "
                 << outputrmi_d[idx[k]][i] << endl
                 << outputrrg_d[idx[k]][i] << "  "
                 << outputrpk_d[idx[k]][i] << endl;
            cout << outputltb_d[idx[k]][i] << "  "
                 << outputlix_d[idx[k]][i] << "  "
                 << outputlmi_d[idx[k]][i] << endl
                 << outputlrg_d[idx[k]][i] << "  "
                 << outputlpk_d[idx[k]][i] << endl;
            errsum = errsum
                   + abs(outputrtb_h[idx[k]][i] - outputrtb_d[idx[k]][i])
                   + abs(outputrix_h[idx[k]][i] - outputrix_d[idx[k]][i])
                   + abs(outputrmi_h[idx[k]][i] - outputrmi_d[idx[k]][i])
                   + abs(outputrrg_h[idx[k]][i] - outputrrg_d[idx[k]][i])
                   + abs(outputrpk_h[idx[k]][i] - outputrpk_d[idx[k]][i])
                   + abs(outputltb_h[idx[k]][i] - outputltb_d[idx[k]][i])
                   + abs(outputlix_h[idx[k]][i] - outputlix_d[idx[k]][i])
                   + abs(outputlmi_h[idx[k]][i] - outputlmi_d[idx[k]][i])
                   + abs(outputlrg_h[idx[k]][i] - outputlrg_d[idx[k]][i])
                   + abs(outputlpk_h[idx[k]][i] - outputlpk_d[idx[k]][i]);
         }
         cout << "Sum of errors : " << errsum << endl; errtot += errsum;
      }
      cout << "Total sum of all errors : " << errtot << endl;
   }
   return errtot;
}

/*
double test_complex ( int dim, int nvr, int pwr, int deg )
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

   srand(time(NULL));

   bool fail = make_complex5_monomial
      (dim,nvr,pwr,deg,idx,exp,
       cffretb,cffreix,cffremi,cffrerg,cffrepk,
       cffimtb,cffimix,cffimmi,cffimrg,cffimpk);

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
      cout << cffretb[i] << "  " << cffreix[i] << "  " << cffremi[i]
           << cffrerg[i] << "  " << cffrepk[i] << endl;
      cout << cffimtb[i] << "  " << cffimix[i] << "  " << cffimmi[i]
           << cffimrg[i] << "  " << cffimpk[i] << endl;
   }
   make_complex5_input
      (dim,deg,inputretb,inputreix,inputremi,inputrerg,inputrepk,
               inputimtb,inputimix,inputimmi,inputimrg,inputimpk);

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
   double errsum = 0.0;
   double errtot = 0.0;

   CPU_cmplx5_evaldiff(dim,nvr,deg,idx,
      cffretb,cffreix,cffremi,cffrerg,cffrepk,
      cffimtb,cffimix,cffimmi,cffimrg,cffimpk,
      inputretb,inputreix,inputremi,inputrerg,inputrepk,
      inputimtb,inputimix,inputimmi,inputimrg,inputimpk,
      outputretb_h,outputreix_h,outputremi_h,outputrerg_h,outputrepk_h,
      outputimtb_h,outputimix_h,outputimmi_h,outputimrg_h,outputimpk_h);

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
      cout << "Coefficient error : " << errsum << endl; errtot += errsum;
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

      cout << "The value of the product computed on the GPU :" << endl;
      errsum = 0.0;
      for(int i=0; i<=deg; i++) 
      {
         cout << outputretb_d[dim][i] << "  " << outputreix_d[dim][i]
                                      << "  " << outputremi_d[dim][i] << endl
              << outputrerg_d[dim][i] << "  " << outputrepk_d[dim][i] << endl;
         cout << outputimtb_d[dim][i] << "  " << outputimix_d[dim][i] 
                                      << "  " << outputimmi_d[dim][i] << endl
              << outputimrg_d[dim][i] << "  " << outputimpk_d[dim][i] << endl;
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
      cout << "The sum of errors : " << errsum << endl; errtot += errsum;

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
         cout << "Coefficient error : " << errsum << endl; errtot += errsum;
      }
   }
   for(int k=0; k<nvr; k++)
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
      if(nvr > 2)
      {
         cout << "-> derivative for index " << idx[k]
              << " computed on GPU :" << endl;
         errsum = 0.0;
         for(int i=0; i<=deg; i++)
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
         cout << "The sum of errors : " << errsum << endl; errtot += errsum;
      }
   }
   cout << "Total sum of all errors : " << errtot << endl;

   return errtot;
}
*/
