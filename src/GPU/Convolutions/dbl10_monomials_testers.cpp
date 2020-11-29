/* The file dbl10_monomials_testers.cpp contains the definitions of
 * functions with prototypes in dbl10_monomials_testers.h. */

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
#include "dbl10_monomials_testers.h"

using namespace std;

int main_dbl10_test
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
   double realsum = test_dbl10_real(dim,nvr,pwr,deg,vrblvl-1);
   double complexsum = test_dbl10_complex(dim,nvr,pwr,deg,vrblvl-1);

   const double tol = 1.0e-152;

   int fail = int(realsum > tol) + int(complexsum > tol);

   if(vrblvl > 0)
   {
      cout << "Sum of all errors in deca double precision :" << endl;
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

double test_dbl10_real ( int dim, int nvr, int pwr, int deg, int verbose )
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

   bool fail = make_real10_monomial
      (dim,nvr,pwr,deg,idx,exp,
       cffrtb,cffrix,cffrmi,cffrrg,cffrpk,
       cffltb,cfflix,cfflmi,cfflrg,cfflpk);

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
         cout << cffrtb[i] << "  " << cffrix[i] << "  " << cffrmi[i] << endl
              << cffrrg[i] << "  " << cffrpk[i] << endl;
         cout << cffltb[i] << "  " << cfflix[i] << "  " << cfflmi[i] << endl
              << cfflrg[i] << "  " << cfflpk[i] << endl;
      }
   }
   make_real10_input(dim,deg,
      inputrtb,inputrix,inputrmi,inputrrg,inputrpk,
      inputltb,inputlix,inputlmi,inputlrg,inputlpk);

   if(verbose > 0)
   {
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

   if(verbose > 0)
   {
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

      if(verbose > 0) cout << "Coefficient error : " << errsum << endl;

      errtot += errsum;
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

      errsum = 0.0;

      if(verbose > 0)
         cout << "The value of the product computed on the GPU :" << endl;

      for(int i=0; i<=deg; i++)
      {
         if(verbose > 0)
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
         }
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
      if(verbose > 0) cout << "Sum of errors : " << errsum << endl;

      errtot += errsum;

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
            }
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
         if(verbose > 0) cout << "Sum of errors : " << errsum << endl;

         errtot += errsum;
      }
      if(verbose > 0) cout << "Total sum of all errors : " << errtot << endl;
   }
   return errtot;
}


double test_dbl10_complex ( int dim, int nvr, int pwr, int deg, int verbose )
{
   int *idx = new int[nvr];       // indices of variables in the monomial
   int *exp = new int[nvr];       // exponents of the variables
   int *expfac = new int[nvr];    // exponents of common factor
   int nbrfac;                    // number of common factors
   double *cffrertb = new double[deg+1]; // highest real parts of coefficients
   double *cffrerix = new double[deg+1]; // second highest real parts
   double *cffrermi = new double[deg+1]; // third highest real parts
   double *cffrerrg = new double[deg+1]; // fourth highest real parts
   double *cffrerpk = new double[deg+1]; // fifth highest real parts
   double *cffreltb = new double[deg+1]; // fifth lowest real parts
   double *cffrelix = new double[deg+1]; // fourth lowest real parts
   double *cffrelmi = new double[deg+1]; // third lowest real parts
   double *cffrelrg = new double[deg+1]; // second lowest real parts
   double *cffrelpk = new double[deg+1]; // lowest real parts
   double *cffimrtb = new double[deg+1]; // highest imaginary parts
   double *cffimrix = new double[deg+1]; // second highest imaginary parts
   double *cffimrmi = new double[deg+1]; // third highest imaginary parts
   double *cffimrrg = new double[deg+1]; // fourth highest imaginary parts
   double *cffimrpk = new double[deg+1]; // fifth highest imaginary parts
   double *cffimltb = new double[deg+1]; // fifth lowest imaginary parts
   double *cffimlix = new double[deg+1]; // fourth lowest imaginary parts
   double *cffimlmi = new double[deg+1]; // third lowest imaginary parts
   double *cffimlrg = new double[deg+1]; // second lowest imaginary parts
   double *cffimlpk = new double[deg+1]; // lowest imaginary parts

   // The input are dim power series of degree deg,
   // the output are nvr+1 power series of degree deg,
   // for the evaluated and differentiated monomial.

   double **inputrertb = new double*[dim];
   double **inputrerix = new double*[dim];
   double **inputrermi = new double*[dim];
   double **inputrerrg = new double*[dim];
   double **inputrerpk = new double*[dim];
   double **inputreltb = new double*[dim];
   double **inputrelix = new double*[dim];
   double **inputrelmi = new double*[dim];
   double **inputrelrg = new double*[dim];
   double **inputrelpk = new double*[dim];
   double **inputimrtb = new double*[dim];
   double **inputimrix = new double*[dim];
   double **inputimrmi = new double*[dim];
   double **inputimrrg = new double*[dim];
   double **inputimrpk = new double*[dim];
   double **inputimltb = new double*[dim];
   double **inputimlix = new double*[dim];
   double **inputimlmi = new double*[dim];
   double **inputimlrg = new double*[dim];
   double **inputimlpk = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      inputrertb[i] = new double[deg+1]; inputrerix[i] = new double[deg+1];
      inputrermi[i] = new double[deg+1]; inputrerrg[i] = new double[deg+1];
      inputrerpk[i] = new double[deg+1];
      inputreltb[i] = new double[deg+1]; inputrelix[i] = new double[deg+1];
      inputrelmi[i] = new double[deg+1]; inputrelrg[i] = new double[deg+1];
      inputrelpk[i] = new double[deg+1];
      inputimrtb[i] = new double[deg+1]; inputimrix[i] = new double[deg+1];
      inputimrmi[i] = new double[deg+1]; inputimrrg[i] = new double[deg+1];
      inputimrpk[i] = new double[deg+1];
      inputimltb[i] = new double[deg+1]; inputimlix[i] = new double[deg+1];
      inputimlmi[i] = new double[deg+1]; inputimlrg[i] = new double[deg+1];
      inputimlpk[i] = new double[deg+1];
   }
   double **outputrertb_h = new double*[dim+1];
   double **outputrerix_h = new double*[dim+1];
   double **outputrermi_h = new double*[dim+1];
   double **outputrerrg_h = new double*[dim+1];
   double **outputrerpk_h = new double*[dim+1];
   double **outputreltb_h = new double*[dim+1];
   double **outputrelix_h = new double*[dim+1];
   double **outputrelmi_h = new double*[dim+1];
   double **outputrelrg_h = new double*[dim+1];
   double **outputrelpk_h = new double*[dim+1];
   double **outputimrtb_h = new double*[dim+1];
   double **outputimrix_h = new double*[dim+1];
   double **outputimrmi_h = new double*[dim+1];
   double **outputimrrg_h = new double*[dim+1];
   double **outputimrpk_h = new double*[dim+1];
   double **outputimltb_h = new double*[dim+1];
   double **outputimlix_h = new double*[dim+1];
   double **outputimlmi_h = new double*[dim+1];
   double **outputimlrg_h = new double*[dim+1];
   double **outputimlpk_h = new double*[dim+1];

   for(int i=0; i<=dim; i++)
   {
      outputrertb_h[i] = new double[deg+1];
      outputrerix_h[i] = new double[deg+1];
      outputrermi_h[i] = new double[deg+1];
      outputrerrg_h[i] = new double[deg+1];
      outputrerpk_h[i] = new double[deg+1];
      outputreltb_h[i] = new double[deg+1];
      outputrelix_h[i] = new double[deg+1];
      outputrelmi_h[i] = new double[deg+1];
      outputrelrg_h[i] = new double[deg+1];
      outputrelpk_h[i] = new double[deg+1];
      outputimrtb_h[i] = new double[deg+1];
      outputimrix_h[i] = new double[deg+1];
      outputimrmi_h[i] = new double[deg+1];
      outputimrrg_h[i] = new double[deg+1];
      outputimrpk_h[i] = new double[deg+1];
      outputimltb_h[i] = new double[deg+1];
      outputimlix_h[i] = new double[deg+1];
      outputimlmi_h[i] = new double[deg+1];
      outputimlrg_h[i] = new double[deg+1];
      outputimlpk_h[i] = new double[deg+1];
   }
   double **outputrertb_d = new double*[dim+1];
   double **outputrerix_d = new double*[dim+1];
   double **outputrermi_d = new double*[dim+1];
   double **outputrerrg_d = new double*[dim+1];
   double **outputrerpk_d = new double*[dim+1];
   double **outputreltb_d = new double*[dim+1];
   double **outputrelix_d = new double*[dim+1];
   double **outputrelmi_d = new double*[dim+1];
   double **outputrelrg_d = new double*[dim+1];
   double **outputrelpk_d = new double*[dim+1];
   double **outputimrtb_d = new double*[dim+1];
   double **outputimrix_d = new double*[dim+1];
   double **outputimrmi_d = new double*[dim+1];
   double **outputimrrg_d = new double*[dim+1];
   double **outputimrpk_d = new double*[dim+1];
   double **outputimltb_d = new double*[dim+1];
   double **outputimlix_d = new double*[dim+1];
   double **outputimlmi_d = new double*[dim+1];
   double **outputimlrg_d = new double*[dim+1];
   double **outputimlpk_d = new double*[dim+1];

   for(int i=0; i<=dim; i++)
   {
      outputrertb_d[i] = new double[deg+1];
      outputrerix_d[i] = new double[deg+1];
      outputrermi_d[i] = new double[deg+1];
      outputrerrg_d[i] = new double[deg+1];
      outputrerpk_d[i] = new double[deg+1];
      outputreltb_d[i] = new double[deg+1];
      outputrelix_d[i] = new double[deg+1];
      outputrelmi_d[i] = new double[deg+1];
      outputrelrg_d[i] = new double[deg+1];
      outputrelpk_d[i] = new double[deg+1];
      outputimrtb_d[i] = new double[deg+1];
      outputimrix_d[i] = new double[deg+1];
      outputimrmi_d[i] = new double[deg+1];
      outputimrrg_d[i] = new double[deg+1];
      outputimrpk_d[i] = new double[deg+1];
      outputimltb_d[i] = new double[deg+1];
      outputimlix_d[i] = new double[deg+1];
      outputimlmi_d[i] = new double[deg+1];
      outputimlrg_d[i] = new double[deg+1];
      outputimlpk_d[i] = new double[deg+1];
   }
   bool fail = make_complex10_monomial
      (dim,nvr,pwr,deg,idx,exp,
       cffrertb,cffrerix,cffrermi,cffrerrg,cffrerpk,
       cffreltb,cffrelix,cffrelmi,cffrelrg,cffrelpk,
       cffimrtb,cffimrix,cffimrmi,cffimrrg,cffimrpk,
       cffimltb,cffimlix,cffimlmi,cffimlrg,cffimlpk);

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
         cout << cffrertb[i] << "  " << cffrerix[i] << "  " << cffrermi[i]
              << cffrerrg[i] << "  " << cffrerpk[i] << endl;
         cout << cffreltb[i] << "  " << cffrelix[i] << "  " << cffrelmi[i]
              << cffrelrg[i] << "  " << cffrelpk[i] << endl;
         cout << cffimrtb[i] << "  " << cffimrix[i] << "  " << cffimrmi[i]
              << cffimrrg[i] << "  " << cffimrpk[i] << endl;
         cout << cffimltb[i] << "  " << cffimlix[i] << "  " << cffimlmi[i]
              << cffimlrg[i] << "  " << cffimlpk[i] << endl;
      }
   }
   make_complex10_input
      (dim,deg,
       inputrertb,inputrerix,inputrermi,inputrerrg,inputrerpk,
       inputreltb,inputrelix,inputrelmi,inputrelrg,inputrelpk,
       inputimrtb,inputimrix,inputimrmi,inputimrrg,inputimrpk,
       inputimltb,inputimlix,inputimlmi,inputimlrg,inputimlpk);

   if(verbose > 0)
   {
      cout << "Random input series :" << endl;
      for(int i=0; i<dim; i++)
      {
         cout << "-> coefficients of series " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
         {
            cout << inputrertb[i][j] << "  " << inputrerix[i][j]
                                     << "  " << inputrermi[i][j] << endl
                 << inputrerrg[i][j] << "  " << inputrerpk[i][j] << endl;
            cout << inputreltb[i][j] << "  " << inputrelix[i][j]
                                     << "  " << inputrelmi[i][j] << endl
                 << inputrelrg[i][j] << "  " << inputrelpk[i][j] << endl;
            cout << inputimrtb[i][j] << "  " << inputimrix[i][j]
                                     << "  " << inputimrmi[i][j] << endl
                 << inputimrrg[i][j] << "  " << inputimrpk[i][j] << endl;
            cout << inputimltb[i][j] << "  " << inputimlix[i][j]
                                     << "  " << inputimlmi[i][j] << endl
                 << inputimlrg[i][j] << "  " << inputimlpk[i][j] << endl;
         }
      }
   }
   double errsum = 0.0;
   double errtot = 0.0;

   CPU_cmplx10_evaldiff(dim,nvr,deg,idx,
      cffrertb,cffrerix,cffrermi,cffrerrg,cffrerpk,
      cffreltb,cffrelix,cffrelmi,cffrelrg,cffrelpk,
      cffimrtb,cffimrix,cffimrmi,cffimrrg,cffimrpk,
      cffimltb,cffimlix,cffimlmi,cffimlrg,cffimlpk,
      inputrertb,inputrerix,inputrermi,inputrerrg,inputrerpk,
      inputreltb,inputrelix,inputrelmi,inputrelrg,inputrelpk,
      inputimrtb,inputimrix,inputimrmi,inputimrrg,inputimrpk,
      inputimltb,inputimlix,inputimlmi,inputimlrg,inputimlpk,
      outputrertb_h,outputrerix_h,outputrermi_h,outputrerrg_h,outputrerpk_h,
      outputreltb_h,outputrelix_h,outputrelmi_h,outputrelrg_h,outputrelpk_h,
      outputimrtb_h,outputimrix_h,outputimrmi_h,outputimrrg_h,outputimrpk_h,
      outputimltb_h,outputimlix_h,outputimlmi_h,outputimlrg_h,outputimlpk_h);

   if(verbose > 0)
   {
      cout << "The value of the product :" << endl;
      for(int i=0; i<=deg; i++)
      {
         cout << outputrertb_h[dim][i] << "  "
              << outputrerix_h[dim][i] << "  "
              << outputrermi_h[dim][i] << endl
              << outputrerrg_h[dim][i] << "  "
              << outputrerpk_h[dim][i] << endl;
         cout << outputreltb_h[dim][i] << "  "
              << outputrelix_h[dim][i] << "  "
              << outputrelmi_h[dim][i] << endl
              << outputrelrg_h[dim][i] << "  "
              << outputrelpk_h[dim][i] << endl;
         cout << outputimrtb_h[dim][i] << "  "
              << outputimrix_h[dim][i] << "  "
              << outputimrmi_h[dim][i] << endl
              << outputimrrg_h[dim][i] << "  "
              << outputimrpk_h[dim][i] << endl;
         cout << outputimltb_h[dim][i] << "  "
              << outputimlix_h[dim][i] << "  "
              << outputimlmi_h[dim][i] << endl
              << outputimlrg_h[dim][i] << "  "
              << outputimlpk_h[dim][i] << endl;
      }
   }
   if(nvr == dim) // the product of all input series equals one
   {
      for(int i=0; i<=deg; i++)
         errsum = errsum
                + abs(outputrertb_h[dim][i] - cffrertb[i])
                + abs(outputrerix_h[dim][i] - cffrerix[i])
                + abs(outputrermi_h[dim][i] - cffrermi[i])
                + abs(outputrerrg_h[dim][i] - cffrerrg[i])
                + abs(outputrerpk_h[dim][i] - cffrerpk[i])
                + abs(outputreltb_h[dim][i] - cffreltb[i])
                + abs(outputrelix_h[dim][i] - cffrelix[i])
                + abs(outputrelmi_h[dim][i] - cffrelmi[i])
                + abs(outputrelrg_h[dim][i] - cffrelrg[i])
                + abs(outputrelpk_h[dim][i] - cffrelpk[i])
                + abs(outputimrtb_h[dim][i] - cffimrtb[i])
                + abs(outputimrix_h[dim][i] - cffimrix[i])
                + abs(outputimrmi_h[dim][i] - cffimrmi[i])
                + abs(outputimrrg_h[dim][i] - cffimrrg[i])
                + abs(outputimrpk_h[dim][i] - cffimrpk[i])
                + abs(outputimltb_h[dim][i] - cffimltb[i])
                + abs(outputimlix_h[dim][i] - cffimlix[i])
                + abs(outputimlmi_h[dim][i] - cffimlmi[i])
                + abs(outputimlrg_h[dim][i] - cffimlrg[i])
                + abs(outputimlpk_h[dim][i] - cffimlpk[i]);

      if(verbose > 0) cout << "Coefficient error : " << errsum << endl;

      errtot += errsum;
   }

   if(nvr > 2)
   {
      GPU_cmplx10_evaldiff(deg+1,dim,nvr,deg,idx,
         cffrertb,cffrerix,cffrermi,cffrerrg,cffrerpk,
         cffreltb,cffrelix,cffrelmi,cffrelrg,cffrelpk,
         cffimrtb,cffimrix,cffimrmi,cffimrrg,cffimrpk,
         cffimltb,cffimlix,cffimlmi,cffimlrg,cffimlpk,
         inputrertb,inputrerix,inputrermi,inputrerrg,inputrerpk,
         inputreltb,inputrelix,inputrelmi,inputrelrg,inputrelpk,
         inputimrtb,inputimrix,inputimrmi,inputimrrg,inputimrpk,
         inputimltb,inputimlix,inputimlmi,inputimlrg,inputimlpk,
         outputrertb_d,outputrerix_d,outputrermi_d,
         outputrerrg_d,outputrerpk_d,
         outputreltb_d,outputrelix_d,outputrelmi_d,
         outputrelrg_d,outputrelpk_d,
         outputimrtb_d,outputimrix_d,outputimrmi_d,
         outputimrrg_d,outputimrpk_d,
         outputimltb_d,outputimlix_d,outputimlmi_d,
         outputimlrg_d,outputimlpk_d);

      errsum = 0.0;

      if(verbose > 0)
         cout << "The value of the product computed on the GPU :" << endl;

      for(int i=0; i<=deg; i++) 
      {
         if(verbose > 0)
         {
            cout << outputrertb_d[dim][i] << "  " << outputrerix_d[dim][i]
                                          << "  " << outputrermi_d[dim][i]
                 << endl
                 << outputrerrg_d[dim][i] << "  " << outputrerpk_d[dim][i]
                 << endl;
            cout << outputreltb_d[dim][i] << "  " << outputrelix_d[dim][i]
                                          << "  " << outputrelmi_d[dim][i]
                 << endl
                 << outputrelrg_d[dim][i] << "  " << outputrelpk_d[dim][i]
                 << endl;
            cout << outputimrtb_d[dim][i] << "  " << outputimrix_d[dim][i] 
                                          << "  " << outputimrmi_d[dim][i]
                 << endl
                 << outputimrrg_d[dim][i] << "  " << outputimrpk_d[dim][i]
                 << endl;
            cout << outputimltb_d[dim][i] << "  " << outputimlix_d[dim][i] 
                                          << "  " << outputimlmi_d[dim][i]
                 << endl
                 << outputimlrg_d[dim][i] << "  " << outputimlpk_d[dim][i]
                 << endl;
         }
         errsum = errsum
                + abs(outputrertb_h[dim][i] - outputrertb_d[dim][i])
                + abs(outputrerix_h[dim][i] - outputrerix_d[dim][i])
                + abs(outputrermi_h[dim][i] - outputrermi_d[dim][i])
                + abs(outputrerrg_h[dim][i] - outputrerrg_d[dim][i])
                + abs(outputrerpk_h[dim][i] - outputrerpk_d[dim][i])
                + abs(outputreltb_h[dim][i] - outputreltb_d[dim][i])
                + abs(outputrelix_h[dim][i] - outputrelix_d[dim][i])
                + abs(outputrelmi_h[dim][i] - outputrelmi_d[dim][i])
                + abs(outputrelrg_h[dim][i] - outputrelrg_d[dim][i])
                + abs(outputrelpk_h[dim][i] - outputrelpk_d[dim][i])
                + abs(outputimrtb_h[dim][i] - outputimrtb_d[dim][i])
                + abs(outputimrix_h[dim][i] - outputimrix_d[dim][i])
                + abs(outputimrmi_h[dim][i] - outputimrmi_d[dim][i])
                + abs(outputimrrg_h[dim][i] - outputimrrg_d[dim][i])
                + abs(outputimrpk_h[dim][i] - outputimrpk_d[dim][i])
                + abs(outputimltb_h[dim][i] - outputimltb_d[dim][i])
                + abs(outputimlix_h[dim][i] - outputimlix_d[dim][i])
                + abs(outputimlmi_h[dim][i] - outputimlmi_d[dim][i])
                + abs(outputimlrg_h[dim][i] - outputimlrg_d[dim][i])
                + abs(outputimlpk_h[dim][i] - outputimlpk_d[dim][i]);
      }
      if(verbose > 0) cout << "The sum of errors : " << errsum << endl;

      errtot += errsum;

      if(nvr == dim) // the product of all input series equals one
      {
         errsum = 0.0;
         for(int i=0; i<=deg; i++)
            errsum = errsum
                   + abs(outputrertb_d[dim][i] - cffrertb[i])
                   + abs(outputrerix_d[dim][i] - cffrerix[i])
                   + abs(outputrermi_d[dim][i] - cffrermi[i])
                   + abs(outputrerrg_d[dim][i] - cffrerrg[i])
                   + abs(outputrerpk_d[dim][i] - cffrerpk[i])
                   + abs(outputreltb_d[dim][i] - cffreltb[i])
                   + abs(outputrelix_d[dim][i] - cffrelix[i])
                   + abs(outputrelmi_d[dim][i] - cffrelmi[i])
                   + abs(outputrelrg_d[dim][i] - cffrelrg[i])
                   + abs(outputrelpk_d[dim][i] - cffrelpk[i])
                   + abs(outputimrtb_d[dim][i] - cffimrtb[i])
                   + abs(outputimrix_d[dim][i] - cffimrix[i])
                   + abs(outputimrmi_d[dim][i] - cffimrmi[i])
                   + abs(outputimrrg_d[dim][i] - cffimrrg[i])
                   + abs(outputimrpk_d[dim][i] - cffimrpk[i])
                   + abs(outputimltb_d[dim][i] - cffimltb[i])
                   + abs(outputimlix_d[dim][i] - cffimlix[i])
                   + abs(outputimlmi_d[dim][i] - cffimlmi[i])
                   + abs(outputimlrg_d[dim][i] - cffimlrg[i])
                   + abs(outputimlpk_d[dim][i] - cffimlpk[i]);

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
            cout << outputrertb_h[idx[k]][i] << "  "
                 << outputrerix_h[idx[k]][i] << "  "
                 << outputrermi_h[idx[k]][i] << "  "
                 << outputrerrg_h[idx[k]][i] << "  "
                 << outputrerpk_h[idx[k]][i] << endl;
            cout << outputreltb_h[idx[k]][i] << "  "
                 << outputrelix_h[idx[k]][i] << "  "
                 << outputrelmi_h[idx[k]][i] << "  "
                 << outputrelrg_h[idx[k]][i] << "  "
                 << outputrelpk_h[idx[k]][i] << endl;
            cout << outputimrtb_h[idx[k]][i] << "  "
                 << outputimrix_h[idx[k]][i] << "  "
                 << outputimrmi_h[idx[k]][i] << "  "
                 << outputimrrg_h[idx[k]][i] << "  "
                 << outputimrpk_h[idx[k]][i] << endl;
            cout << outputimltb_h[idx[k]][i] << "  "
                 << outputimlix_h[idx[k]][i] << "  "
                 << outputimlmi_h[idx[k]][i] << "  "
                 << outputimlrg_h[idx[k]][i] << "  "
                 << outputimlpk_h[idx[k]][i] << endl;
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
               cout << outputrertb_d[idx[k]][i] << "  "
                    << outputrerix_d[idx[k]][i] << "  "
                    << outputrermi_d[idx[k]][i] << "  "
                    << outputrerrg_d[idx[k]][i] << "  "
                    << outputrerpk_d[idx[k]][i] << endl;
               cout << outputreltb_d[idx[k]][i] << "  "
                    << outputrelix_d[idx[k]][i] << "  "
                    << outputrelmi_d[idx[k]][i] << "  "
                    << outputrelrg_d[idx[k]][i] << "  "
                    << outputrelpk_d[idx[k]][i] << endl;
               cout << outputimrtb_d[idx[k]][i] << "  "
                    << outputimrix_d[idx[k]][i] << "  "
                    << outputimrmi_d[idx[k]][i] << "  "
                    << outputimrrg_d[idx[k]][i] << "  "
                    << outputimrpk_d[idx[k]][i] << endl;
               cout << outputimltb_d[idx[k]][i] << "  "
                    << outputimlix_d[idx[k]][i] << "  "
                    << outputimlmi_d[idx[k]][i] << "  "
                    << outputimlrg_d[idx[k]][i] << "  "
                    << outputimlpk_d[idx[k]][i] << endl;
            }
            errsum = errsum
                   + abs(outputrertb_h[idx[k]][i] - outputrertb_d[idx[k]][i])
                   + abs(outputrerix_h[idx[k]][i] - outputrerix_d[idx[k]][i])
                   + abs(outputrermi_h[idx[k]][i] - outputrermi_d[idx[k]][i])
                   + abs(outputrerrg_h[idx[k]][i] - outputrerrg_d[idx[k]][i])
                   + abs(outputrerpk_h[idx[k]][i] - outputrerpk_d[idx[k]][i])
                   + abs(outputreltb_h[idx[k]][i] - outputreltb_d[idx[k]][i])
                   + abs(outputrelix_h[idx[k]][i] - outputrelix_d[idx[k]][i])
                   + abs(outputrelmi_h[idx[k]][i] - outputrelmi_d[idx[k]][i])
                   + abs(outputrelrg_h[idx[k]][i] - outputrelrg_d[idx[k]][i])
                   + abs(outputrelpk_h[idx[k]][i] - outputrelpk_d[idx[k]][i])
                   + abs(outputimrtb_h[idx[k]][i] - outputimrtb_d[idx[k]][i])
                   + abs(outputimrix_h[idx[k]][i] - outputimrix_d[idx[k]][i])
                   + abs(outputimrmi_h[idx[k]][i] - outputimrmi_d[idx[k]][i])
                   + abs(outputimrrg_h[idx[k]][i] - outputimrrg_d[idx[k]][i])
                   + abs(outputimrpk_h[idx[k]][i] - outputimrpk_d[idx[k]][i])
                   + abs(outputimltb_h[idx[k]][i] - outputimltb_d[idx[k]][i])
                   + abs(outputimlix_h[idx[k]][i] - outputimlix_d[idx[k]][i])
                   + abs(outputimlmi_h[idx[k]][i] - outputimlmi_d[idx[k]][i])
                   + abs(outputimlrg_h[idx[k]][i] - outputimlrg_d[idx[k]][i])
                   + abs(outputimlpk_h[idx[k]][i] - outputimlpk_d[idx[k]][i]);
         }
         if(verbose > 0) cout << "The sum of errors : " << errsum << endl;

         errtot += errsum;
      }
   }
   if(verbose > 0) cout << "Total sum of all errors : " << errtot << endl;

   return errtot;
}
