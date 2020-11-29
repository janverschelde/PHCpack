/* The file dbl4_monomials_testers.cpp contains the definitions of
 * functions with prototypes in dbl4_monomials_testers.h. */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <vector_types.h>
#include "random_monomials.h"
#include "random4_monomials.h"
#include "dbl4_monomials_host.h"
#include "dbl4_monomials_kernels.h"
#include "dbl4_monomials_testers.h"

using namespace std;

int main_dbl4_test
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
   double realsum = test_dbl4_real(dim,nvr,pwr,deg,vrblvl-1);
   double complexsum = test_dbl4_complex(dim,nvr,pwr,deg,vrblvl-1);

   const double tol = 1.0e-60;

   int fail = int(realsum > tol) + int(complexsum > tol);

   if(vrblvl > 0)
   {
      cout << "Sum of all errors in quad double precision :" << endl;
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

double test_dbl4_real ( int dim, int nvr, int pwr, int deg, int verbose )
{
   int *idx = new int[nvr];           // indices of variables in the monomial
   int *exp = new int[nvr];           // exponents of the variables
   int *expfac = new int[nvr];        // exponents of common factor
   int nbrfac;                        // number of common factors
   double *cffhihi = new double[deg+1]; // highest coefficient doubles
   double *cfflohi = new double[deg+1]; // second highest coefficient doubles
   double *cffhilo = new double[deg+1]; // second lowest coefficient doubles
   double *cfflolo = new double[deg+1]; // lowest coefficient doubles

   // The input are dim power series of degree deg,
   // the output are nvr+1 power series of degree deg,
   // for the evaluated and differentiated monomial.

   double **inputhihi = new double*[dim];
   for(int i=0; i<dim; i++) inputhihi[i] = new double[deg+1];
   double **inputlohi = new double*[dim];
   for(int i=0; i<dim; i++) inputlohi[i] = new double[deg+1];
   double **inputhilo = new double*[dim];
   for(int i=0; i<dim; i++) inputhilo[i] = new double[deg+1];
   double **inputlolo = new double*[dim];
   for(int i=0; i<dim; i++) inputlolo[i] = new double[deg+1];
   double **outputhihi_h = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputhihi_h[i] = new double[deg+1];
   double **outputlohi_h = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputlohi_h[i] = new double[deg+1];
   double **outputhilo_h = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputhilo_h[i] = new double[deg+1];
   double **outputlolo_h = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputlolo_h[i] = new double[deg+1];
   double **outputhihi_d = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputhihi_d[i] = new double[deg+1];
   double **outputlohi_d = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputlohi_d[i] = new double[deg+1];
   double **outputhilo_d = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputhilo_d[i] = new double[deg+1];
   double **outputlolo_d = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputlolo_d[i] = new double[deg+1];

   bool fail = make_real4_monomial(dim,nvr,pwr,deg,idx,exp,
                                   cffhihi,cfflohi,cffhilo,cfflolo);

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
         cout << " " << cffhihi[i] << " " << cfflohi[i] << endl;
         cout << " " << cffhilo[i] << " " << cfflolo[i] << endl;
      }
   }
   make_real4_input(dim,deg,inputhihi,inputlohi,inputhilo,inputlolo);

   if(verbose > 0)
   {
      cout << "Random input series :" << endl;
      for(int i=0; i<dim; i++)
      {
         cout << "-> coefficients of series " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
         {
            cout << inputhihi[i][j] << "  " << inputlohi[i][j] << endl;
            cout << inputhilo[i][j] << "  " << inputlolo[i][j] << endl;
         }
      }
   }
   double errsum = 0.0;
   double errtot = 0.0;
 
   CPU_dbl4_evaldiff(dim,nvr,deg,idx,cffhihi,cfflohi,cffhilo,cfflolo,
                     inputhihi,inputlohi,inputhilo,inputlolo,
                     outputhihi_h,outputlohi_h,outputhilo_h,outputlolo_h);
   if(verbose > 0)
   {
      cout << "The value of the product :" << endl;
      for(int i=0; i<=deg; i++)
      {
         cout << outputhihi_h[dim][i] << "  " << outputlohi_h[dim][i] << endl;
         cout << outputhilo_h[dim][i] << "  " << outputlolo_h[dim][i] << endl;
      }
   }
   if(nvr == dim) // the product of all input series equals one
   {
      for(int i=0; i<=deg; i++)
         errsum = errsum
                + abs(outputhihi_h[dim][i] - cffhihi[i])
                + abs(outputlohi_h[dim][i] - cfflohi[i])
                + abs(outputhilo_h[dim][i] - cffhilo[i])
                + abs(outputlolo_h[dim][i] - cfflolo[i]);

      if(verbose > 0) cout << "Coefficient error : " << errsum << endl;

      errtot += errsum;
   }
   if(nvr > 2)
   {
      GPU_dbl4_evaldiff(deg+1,dim,nvr,deg,idx,
         cffhihi,cfflohi,cffhilo,cfflolo,
         inputhihi,inputlohi,inputhilo,inputlolo,
         outputhihi_d,outputlohi_d,outputhilo_d,outputlolo_d);

      if(verbose > 0)
         cout << "The value of the product computed on the GPU :" << endl;

      for(int i=0; i<=deg; i++)
      {
         if(verbose > 0)
         {
            cout << outputhihi_d[dim][i] << "  "
                 << outputlohi_d[dim][i] << endl;
            cout << outputhilo_d[dim][i] << "  "
                 << outputlolo_d[dim][i] << endl;
         }
         errsum = errsum 
                + abs(outputhihi_h[dim][i] - outputhihi_d[dim][i])
                + abs(outputlohi_h[dim][i] - outputlohi_d[dim][i])
                + abs(outputhilo_h[dim][i] - outputhilo_d[dim][i])
                + abs(outputlolo_h[dim][i] - outputlolo_d[dim][i]);
      }
      if(verbose > 0) cout << "Sum of errors : " << errsum << endl;

      errtot += errsum;

      if(nvr == dim) // the product of all input series equals one
      {
         errsum = 0.0;
         for(int i=0; i<=deg; i++)
            errsum = errsum
                   + abs(outputhihi_d[dim][i] - cffhihi[i])
                   + abs(outputlohi_d[dim][i] - cfflohi[i])
                   + abs(outputhilo_d[dim][i] - cffhilo[i])
                   + abs(outputlolo_d[dim][i] - cfflolo[i]);

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
            cout << outputhihi_h[idx[k]][i] << "  "
                 << outputlohi_h[idx[k]][i] << endl;
            cout << outputhilo_h[idx[k]][i] << "  "
                 << outputlolo_h[idx[k]][i] << endl;
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
               cout << outputhihi_d[idx[k]][i] << "  "
                    << outputlohi_d[idx[k]][i] << endl;
               cout << outputhilo_d[idx[k]][i] << "  "
                    << outputlolo_d[idx[k]][i] << endl;
            }
            errsum = errsum
                   + abs(outputhihi_h[idx[k]][i] - outputhihi_d[idx[k]][i])
                   + abs(outputlohi_h[idx[k]][i] - outputlohi_d[idx[k]][i])
                   + abs(outputhilo_h[idx[k]][i] - outputhilo_d[idx[k]][i])
                   + abs(outputlolo_h[idx[k]][i] - outputlolo_d[idx[k]][i]);
         }
         if(verbose > 0) cout << "Sum of errors : " << errsum << endl;

         errtot += errsum;
      }
   }
   if(verbose > 0) cout << "Total sum of all errors : " << errtot << endl;

   return errtot;
}

double test_dbl4_complex ( int dim, int nvr, int pwr, int deg, int verbose )
{
   int *idx = new int[nvr];             // indices of variables in the monomial
   int *exp = new int[nvr];             // exponents of the variables
   int *expfac = new int[nvr];          // exponents of common factor
   int nbrfac;                          // number of common factors
   double *cffrehihi = new double[deg+1]; // highest real coefficients parts
   double *cffrelohi = new double[deg+1]; // second highest real parts
   double *cffrehilo = new double[deg+1]; // second lowest real parts 
   double *cffrelolo = new double[deg+1]; // lowest real parts
   double *cffimhihi = new double[deg+1]; // highest imaginary parts
   double *cffimlohi = new double[deg+1]; // second highest imaginary parts
   double *cffimhilo = new double[deg+1]; // second lowest imaginary parts 
   double *cffimlolo = new double[deg+1]; // lowest imaginary parts

   // The input are dim power series of degree deg,
   // the output are nvr+1 power series of degree deg,
   // for the evaluated and differentiated monomial.

   double **inputrehihi = new double*[dim];
   double **inputrelohi = new double*[dim];
   double **inputrehilo = new double*[dim];
   double **inputrelolo = new double*[dim];
   double **inputimhihi = new double*[dim];
   double **inputimlohi = new double*[dim];
   double **inputimhilo = new double*[dim];
   double **inputimlolo = new double*[dim];
   for(int i=0; i<dim; i++)
   {
      inputrehihi[i] = new double[deg+1];
      inputrelohi[i] = new double[deg+1];
      inputrehilo[i] = new double[deg+1];
      inputrelolo[i] = new double[deg+1];
      inputimhihi[i] = new double[deg+1];
      inputimlohi[i] = new double[deg+1];
      inputimhilo[i] = new double[deg+1];
      inputimlolo[i] = new double[deg+1];
   }
   double **outputrehihi_h = new double*[dim+1];
   double **outputrelohi_h = new double*[dim+1];
   double **outputrehilo_h = new double*[dim+1];
   double **outputrelolo_h = new double*[dim+1];
   double **outputimhihi_h = new double*[dim+1];
   double **outputimlohi_h = new double*[dim+1];
   double **outputimhilo_h = new double*[dim+1];
   double **outputimlolo_h = new double*[dim+1];
   for(int i=0; i<=dim; i++)
   {
      outputrehihi_h[i] = new double[deg+1];
      outputrelohi_h[i] = new double[deg+1];
      outputrehilo_h[i] = new double[deg+1];
      outputrelolo_h[i] = new double[deg+1];
      outputimhihi_h[i] = new double[deg+1];
      outputimlohi_h[i] = new double[deg+1];
      outputimhilo_h[i] = new double[deg+1];
      outputimlolo_h[i] = new double[deg+1];
   }
   double **outputrehihi_d = new double*[dim+1];
   double **outputrelohi_d = new double*[dim+1];
   double **outputrehilo_d = new double*[dim+1];
   double **outputrelolo_d = new double*[dim+1];
   double **outputimhihi_d = new double*[dim+1];
   double **outputimlohi_d = new double*[dim+1];
   double **outputimhilo_d = new double*[dim+1];
   double **outputimlolo_d = new double*[dim+1];
   for(int i=0; i<=dim; i++)
   {
      outputrehihi_d[i] = new double[deg+1];
      outputrelohi_d[i] = new double[deg+1];
      outputrehilo_d[i] = new double[deg+1];
      outputrelolo_d[i] = new double[deg+1];
      outputimhihi_d[i] = new double[deg+1];
      outputimlohi_d[i] = new double[deg+1];
      outputimhilo_d[i] = new double[deg+1];
      outputimlolo_d[i] = new double[deg+1];
   }
   bool fail = make_complex4_monomial
     (dim,nvr,pwr,deg,idx,exp,
      cffrehihi,cffrelohi,cffrehilo,cffrelolo,
      cffimhihi,cffimlohi,cffimhilo,cffimlolo);

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
         cout << cffrehihi[i] << "  " << cffrelohi[i] << endl;
         cout << cffrehilo[i] << "  " << cffrelolo[i] << endl;
         cout << cffimhihi[i] << "  " << cffimlohi[i] << endl;
         cout << cffimhilo[i] << "  " << cffimlolo[i] << endl;
      }
   }
   make_complex4_input(dim,deg,
      inputrehihi,inputrelohi,inputrehilo,inputrelolo,
      inputimhihi,inputimlohi,inputimhilo,inputimlolo);

   if(verbose > 0)
   {
      cout << "Random input series :" << endl;
      for(int i=0; i<dim; i++)
      {
         cout << "-> coefficients of series " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
         {
            cout << inputrehihi[i][j] << "  " << inputrelohi[i][j] << endl;
            cout << inputrehilo[i][j] << "  " << inputrelolo[i][j] << endl;
            cout << inputimhihi[i][j] << "  " << inputimlohi[i][j] << endl;
            cout << inputimhilo[i][j] << "  " << inputimlolo[i][j] << endl;
         }
      }
   }
   double errsum = 0.0;
   double errtot = 0.0;

   CPU_cmplx4_evaldiff
      (dim,nvr,deg,idx,
       cffrehihi,cffrelohi,cffrehilo,cffrelolo,
       cffimhihi,cffimlohi,cffimhilo,cffimlolo,
       inputrehihi,inputrelohi,inputrehilo,inputrelolo,
       inputimhihi,inputimlohi,inputimhilo,inputimlolo,
       outputrehihi_h,outputrelohi_h,outputrehilo_h,outputrelolo_h,
       outputimhihi_h,outputimlohi_h,outputimhilo_h,outputimlolo_h);

   if(verbose > 0)
   {
      cout << "The value of the product :" << endl;
      for(int i=0; i<=deg; i++)
      {
         cout << outputrehihi_h[dim][i] << "  "
              << outputrelohi_h[dim][i] << endl;
         cout << outputrehilo_h[dim][i] << "  "
              << outputrelolo_h[dim][i] << endl;
         cout << outputimhihi_h[dim][i] << "  "
              << outputimlohi_h[dim][i] << endl;
         cout << outputimhilo_h[dim][i] << "  "
              << outputimlolo_h[dim][i] << endl;
      }
   }
   if(nvr == dim) // the product of all input series equals one
   {
      for(int i=0; i<=deg; i++)
         errsum = errsum
                + abs(outputrehihi_h[dim][i] - cffrehihi[i])
                + abs(outputrelohi_h[dim][i] - cffrelohi[i])
                + abs(outputrehilo_h[dim][i] - cffrehilo[i])
                + abs(outputrelolo_h[dim][i] - cffrelolo[i])
                + abs(outputimhihi_h[dim][i] - cffimhihi[i])
                + abs(outputimlohi_h[dim][i] - cffimlohi[i])
                + abs(outputimhilo_h[dim][i] - cffimhilo[i])
                + abs(outputimlolo_h[dim][i] - cffimlolo[i]);

      if(verbose > 0) cout << "Coefficient error : " << errsum << endl;

      errtot += errsum;
   }
   if(nvr > 2)
   {
      GPU_cmplx4_evaldiff(deg+1,dim,nvr,deg,idx,
         cffrehihi,cffrelohi,cffrehilo,cffrelolo,
         cffimhihi,cffimlohi,cffimhilo,cffimlolo,
         inputrehihi,inputrelohi,inputrehilo,inputrelolo,
         inputimhihi,inputimlohi,inputimhilo,inputimlolo,
         outputrehihi_d,outputrelohi_d,outputrehilo_d,outputrelolo_d,
         outputimhihi_d,outputimlohi_d,outputimhilo_d,outputimlolo_d);

      errsum = 0.0;

      if(verbose > 0)
         cout << "The value of the product computed on the GPU :" << endl;

      for(int i=0; i<=deg; i++) 
      {
         if(verbose > 0)
         {
            cout << outputrehihi_d[dim][i] << "  "
                 << outputrelohi_d[dim][i] << endl;
            cout << outputrehilo_d[dim][i] << "  "
                 << outputrelolo_d[dim][i] << endl;
            cout << outputimhihi_d[dim][i] << "  "
                 << outputimlohi_d[dim][i] << endl;
            cout << outputimhilo_d[dim][i] << "  "
                 << outputimlolo_d[dim][i] << endl;
         }
         errsum = errsum
                + abs(outputrehihi_h[dim][i] - outputrehihi_d[dim][i])
                + abs(outputrelohi_h[dim][i] - outputrelohi_d[dim][i])
                + abs(outputrehilo_h[dim][i] - outputrehilo_d[dim][i])
                + abs(outputrelolo_h[dim][i] - outputrelolo_d[dim][i])
                + abs(outputimhihi_h[dim][i] - outputimhihi_d[dim][i])
                + abs(outputimlohi_h[dim][i] - outputimlohi_d[dim][i])
                + abs(outputimhilo_h[dim][i] - outputimhilo_d[dim][i])
                + abs(outputimlolo_h[dim][i] - outputimlolo_d[dim][i]);
      }
      if(verbose > 0) cout << "The sum of errors : " << errsum << endl;

      errtot += errsum;

      if(nvr == dim) // the product of all input series equals one
      {
         errsum = 0.0;
         for(int i=0; i<=deg; i++)
            errsum = errsum
                   + abs(outputrehihi_d[dim][i] - cffrehihi[i])
                   + abs(outputrelohi_d[dim][i] - cffrelohi[i])
                   + abs(outputrehilo_d[dim][i] - cffrehilo[i])
                   + abs(outputrelolo_d[dim][i] - cffrelolo[i])
                   + abs(outputimhihi_d[dim][i] - cffimhihi[i])
                   + abs(outputimlohi_d[dim][i] - cffimlohi[i])
                   + abs(outputimhilo_d[dim][i] - cffimhilo[i])
                   + abs(outputimlolo_d[dim][i] - cffimlolo[i]);

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
            cout << outputrehihi_h[idx[k]][i] << "  "
                 << outputrelohi_h[idx[k]][i] << endl;
            cout << outputrehilo_h[idx[k]][i] << "  "
                 << outputrelolo_h[idx[k]][i] << endl;
            cout << outputimhihi_h[idx[k]][i] << "  "
                 << outputimlohi_h[idx[k]][i] << endl;
            cout << outputimhilo_h[idx[k]][i] << "  "
                 << outputimlolo_h[idx[k]][i] << endl;
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
               cout << outputrehihi_d[idx[k]][i] << "  "
                    << outputrelohi_d[idx[k]][i] << endl;
               cout << outputrehilo_d[idx[k]][i] << "  "
                    << outputrelolo_d[idx[k]][i] << endl;
               cout << outputimhihi_d[idx[k]][i] << "  "
                    << outputimlohi_d[idx[k]][i] << endl;
               cout << outputimhilo_d[idx[k]][i] << "  "
                    << outputimlolo_d[idx[k]][i] << endl;
            }
            errsum = errsum
               + abs(outputrehihi_h[idx[k]][i] - outputrehihi_d[idx[k]][i])
               + abs(outputrelohi_h[idx[k]][i] - outputrelohi_d[idx[k]][i])
               + abs(outputrehilo_h[idx[k]][i] - outputrehilo_d[idx[k]][i])
               + abs(outputrelolo_h[idx[k]][i] - outputrelolo_d[idx[k]][i])
               + abs(outputimhihi_h[idx[k]][i] - outputimhihi_d[idx[k]][i])
               + abs(outputimlohi_h[idx[k]][i] - outputimlohi_d[idx[k]][i])
               + abs(outputimhilo_h[idx[k]][i] - outputimhilo_d[idx[k]][i])
               + abs(outputimlolo_h[idx[k]][i] - outputimlolo_d[idx[k]][i]);
         }
         if(verbose > 0) cout << "The sum of errors : " << errsum << endl;

         errtot += errsum;
      }
   }
   if(verbose > 0) cout << "Total sum of all errors : " << errtot << endl;

   return errtot;
}
