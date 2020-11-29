/* The file dbl8_monomials_testers.cpp contains the definitions of
 * functions with prototypes in dbl8_monomials_testers.h. */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <vector_types.h>
#include "random_monomials.h"
#include "random8_monomials.h"
#include "dbl8_monomials_host.h"
#include "dbl8_monomials_kernels.h"
#include "dbl8_monomials_testers.h"

using namespace std;

int main_dbl8_test
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
   double realsum = test_dbl8_real(dim,nvr,pwr,deg,vrblvl-1);
   double complexsum = test_dbl8_complex(dim,nvr,pwr,deg,vrblvl-1);

   const double tol = 1.0e-120;

   int fail = int(realsum > tol) + int(complexsum > tol);

   if(vrblvl > 0)
   {
      cout << "Sum of all errors in octo double precision :" << endl;
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

double test_dbl8_real ( int dim, int nvr, int pwr, int deg, int verbose )
{
   int *idx = new int[nvr];           // indices of variables in the monomial
   int *exp = new int[nvr];           // exponents of the variables
   int *expfac = new int[nvr];        // exponents of common factor
   int nbrfac;                        // number of common factors
   double *cffhihihi = new double[deg+1]; // highest coefficient doubles
   double *cfflohihi = new double[deg+1]; // second highest coefficients
   double *cffhilohi = new double[deg+1]; // third highest coefficients
   double *cfflolohi = new double[deg+1]; // fourth highest coefficients
   double *cffhihilo = new double[deg+1]; // fourth lowest coefficients
   double *cfflohilo = new double[deg+1]; // third lowest coefficients
   double *cffhilolo = new double[deg+1]; // second lowest coefficients
   double *cfflololo = new double[deg+1]; // lowest coefficient doubles

   // The input are dim power series of degree deg,
   // the output are nvr+1 power series of degree deg,
   // for the evaluated and differentiated monomial.

   double **inputhihihi = new double*[dim];
   for(int i=0; i<dim; i++) inputhihihi[i] = new double[deg+1];
   double **inputlohihi = new double*[dim];
   for(int i=0; i<dim; i++) inputlohihi[i] = new double[deg+1];
   double **inputhilohi = new double*[dim];
   for(int i=0; i<dim; i++) inputhilohi[i] = new double[deg+1];
   double **inputlolohi = new double*[dim];
   for(int i=0; i<dim; i++) inputlolohi[i] = new double[deg+1];
   double **inputhihilo = new double*[dim];
   for(int i=0; i<dim; i++) inputhihilo[i] = new double[deg+1];
   double **inputlohilo = new double*[dim];
   for(int i=0; i<dim; i++) inputlohilo[i] = new double[deg+1];
   double **inputhilolo = new double*[dim];
   for(int i=0; i<dim; i++) inputhilolo[i] = new double[deg+1];
   double **inputlololo = new double*[dim];
   for(int i=0; i<dim; i++) inputlololo[i] = new double[deg+1];
   double **outputhihihi_h = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputhihihi_h[i] = new double[deg+1];
   double **outputlohihi_h = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputlohihi_h[i] = new double[deg+1];
   double **outputhilohi_h = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputhilohi_h[i] = new double[deg+1];
   double **outputlolohi_h = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputlolohi_h[i] = new double[deg+1];
   double **outputhihilo_h = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputhihilo_h[i] = new double[deg+1];
   double **outputlohilo_h = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputlohilo_h[i] = new double[deg+1];
   double **outputhilolo_h = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputhilolo_h[i] = new double[deg+1];
   double **outputlololo_h = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputlololo_h[i] = new double[deg+1];
   double **outputhihihi_d = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputhihihi_d[i] = new double[deg+1];
   double **outputlohihi_d = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputlohihi_d[i] = new double[deg+1];
   double **outputhilohi_d = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputhilohi_d[i] = new double[deg+1];
   double **outputlolohi_d = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputlolohi_d[i] = new double[deg+1];
   double **outputhihilo_d = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputhihilo_d[i] = new double[deg+1];
   double **outputlohilo_d = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputlohilo_d[i] = new double[deg+1];
   double **outputhilolo_d = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputhilolo_d[i] = new double[deg+1];
   double **outputlololo_d = new double*[dim+1];
   for(int i=0; i<=dim; i++) outputlololo_d[i] = new double[deg+1];

   bool fail = make_real8_monomial
     (dim,nvr,pwr,deg,idx,exp,
      cffhihihi,cfflohihi,cffhilohi,cfflolohi,
      cffhihilo,cfflohilo,cffhilolo,cfflololo);

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
         cout << " " << cffhihihi[i] << " " << cfflohihi[i] << endl;
         cout << " " << cffhilohi[i] << " " << cfflolohi[i] << endl;
         cout << " " << cffhihilo[i] << " " << cfflohilo[i] << endl;
         cout << " " << cffhilolo[i] << " " << cfflololo[i] << endl;
      }
   }
   make_real8_input(dim,deg,
      inputhihihi,inputlohihi,inputhilohi,inputlolohi,
      inputhihilo,inputlohilo,inputhilolo,inputlololo);

   if(verbose > 0)
   {
      cout << "Random input series :" << endl;
      for(int i=0; i<dim; i++)
      {
         cout << "-> coefficients of series " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
         {
            cout << inputhihihi[i][j] << "  " << inputlohihi[i][j] << endl;
            cout << inputhilohi[i][j] << "  " << inputlolohi[i][j] << endl;
            cout << inputhihilo[i][j] << "  " << inputlohilo[i][j] << endl;
            cout << inputhilolo[i][j] << "  " << inputlololo[i][j] << endl;
         }
      }
   }
   double errsum = 0.0;
   double errtot = 0.0;
 
   CPU_dbl8_evaldiff(dim,nvr,deg,idx,
      cffhihihi,cfflohihi,cffhilohi,cfflolohi,
      cffhihilo,cfflohilo,cffhilolo,cfflololo,
      inputhihihi,inputlohihi,inputhilohi,inputlolohi,
      inputhihilo,inputlohilo,inputhilolo,inputlololo,
      outputhihihi_h,outputlohihi_h,outputhilohi_h,outputlolohi_h,
      outputhihilo_h,outputlohilo_h,outputhilolo_h,outputlololo_h);

   if(verbose > 0)
   {
      cout << "The value of the product :" << endl;
      for(int i=0; i<=deg; i++)
      {
         cout << outputhihihi_h[dim][i] << "  "
              << outputlohihi_h[dim][i] << endl;
         cout << outputlohihi_h[dim][i] << "  "
              << outputlolohi_h[dim][i] << endl;
         cout << outputhihilo_h[dim][i] << "  "
              << outputlohilo_h[dim][i] << endl;
         cout << outputhilolo_h[dim][i] << "  "
              << outputlololo_h[dim][i] << endl;
      }
   }
   if(nvr == dim) // the product of all input series equals one
   {
      for(int i=0; i<=deg; i++)
         errsum = errsum
                + abs(outputhihihi_h[dim][i] - cffhihihi[i])
                + abs(outputlohihi_h[dim][i] - cfflohihi[i])
                + abs(outputhilohi_h[dim][i] - cffhilohi[i])
                + abs(outputlolohi_h[dim][i] - cfflolohi[i])
                + abs(outputhihilo_h[dim][i] - cffhihilo[i])
                + abs(outputlohilo_h[dim][i] - cfflohilo[i])
                + abs(outputhilolo_h[dim][i] - cffhilolo[i])
                + abs(outputlololo_h[dim][i] - cfflololo[i]);

      if(verbose > 0) cout << "Coefficient error : " << errsum << endl;

      errtot += errsum;
   }

   if(nvr > 2)
   {
      GPU_dbl8_evaldiff(deg+1,dim,nvr,deg,idx,
         cffhihihi,cfflohihi,cffhilohi,cfflolohi,
         cffhihilo,cfflohilo,cffhilolo,cfflololo,
         inputhihihi,inputlohihi,inputhilohi,inputlolohi,
         inputhihilo,inputlohilo,inputhilolo,inputlololo,
         outputhihihi_d,outputlohihi_d,outputhilohi_d,outputlolohi_d,
         outputhihilo_d,outputlohilo_d,outputhilolo_d,outputlololo_d);

      if(verbose > 0)
         cout << "The value of the product computed on the GPU :" << endl;
       
      for(int i=0; i<=deg; i++)
      {
         if(verbose > 0)
         {
            cout << outputhihihi_d[dim][i] << "  "
                 << outputlohihi_d[dim][i] << endl;
            cout << outputhilohi_d[dim][i] << "  "
                 << outputlolohi_d[dim][i] << endl;
            cout << outputhihilo_d[dim][i] << "  "
                 << outputlohilo_d[dim][i] << endl;
            cout << outputhilolo_d[dim][i] << "  "
                 << outputlololo_d[dim][i] << endl;
         }
         errsum = errsum 
             + abs(outputhihihi_h[dim][i] - outputhihihi_d[dim][i])
             + abs(outputlohihi_h[dim][i] - outputlohihi_d[dim][i])
             + abs(outputhilohi_h[dim][i] - outputhilohi_d[dim][i])
             + abs(outputlolohi_h[dim][i] - outputlolohi_d[dim][i])
             + abs(outputhihilo_h[dim][i] - outputhihilo_d[dim][i])
             + abs(outputlohilo_h[dim][i] - outputlohilo_d[dim][i])
             + abs(outputhilolo_h[dim][i] - outputhilolo_d[dim][i])
             + abs(outputlololo_h[dim][i] - outputlololo_d[dim][i]);
      }
      if(verbose > 0) cout << "Sum of errors : " << errsum << endl;

      errtot += errsum;

      if(nvr == dim) // the product of all input series equals one
      {
         errsum = 0.0;
         for(int i=0; i<=deg; i++)
            errsum = errsum
                   + abs(outputhihihi_d[dim][i] - cffhihihi[i])
                   + abs(outputlohihi_d[dim][i] - cfflohihi[i])
                   + abs(outputhilohi_d[dim][i] - cffhilohi[i])
                   + abs(outputlolohi_d[dim][i] - cfflolohi[i])
                   + abs(outputhihilo_d[dim][i] - cffhihilo[i])
                   + abs(outputlohilo_d[dim][i] - cfflohilo[i])
                   + abs(outputhilolo_d[dim][i] - cffhilolo[i])
                   + abs(outputlololo_d[dim][i] - cfflololo[i]);

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
            cout << outputhihihi_h[idx[k]][i] << "  "
                 << outputlohihi_h[idx[k]][i] << endl;
            cout << outputhilohi_h[idx[k]][i] << "  "
                 << outputlolohi_h[idx[k]][i] << endl;
            cout << outputhihilo_h[idx[k]][i] << "  "
                 << outputlohilo_h[idx[k]][i] << endl;
            cout << outputhilolo_h[idx[k]][i] << "  "
                 << outputlololo_h[idx[k]][i] << endl;
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
               cout << outputhihihi_d[idx[k]][i] << "  "
                    << outputlohihi_d[idx[k]][i] << endl;
               cout << outputhilohi_d[idx[k]][i] << "  "
                    << outputlolohi_d[idx[k]][i] << endl;
               cout << outputhihilo_d[idx[k]][i] << "  "
                    << outputlohilo_d[idx[k]][i] << endl;
               cout << outputhilolo_d[idx[k]][i] << "  "
                    << outputlololo_d[idx[k]][i] << endl;
            }
            errsum = errsum
               + abs(outputhihihi_h[idx[k]][i] - outputhihihi_d[idx[k]][i])
               + abs(outputlohihi_h[idx[k]][i] - outputlohihi_d[idx[k]][i])
               + abs(outputhilohi_h[idx[k]][i] - outputhilohi_d[idx[k]][i])
               + abs(outputlolohi_h[idx[k]][i] - outputlolohi_d[idx[k]][i])
               + abs(outputhihilo_h[idx[k]][i] - outputhihilo_d[idx[k]][i])
               + abs(outputlohilo_h[idx[k]][i] - outputlohilo_d[idx[k]][i])
               + abs(outputhilolo_h[idx[k]][i] - outputhilolo_d[idx[k]][i])
               + abs(outputlololo_h[idx[k]][i] - outputlololo_d[idx[k]][i]);
         }
         if(verbose > 0) cout << "Sum of errors : " << errsum << endl;

         errtot += errsum;
      }
   }
   if(verbose > 0) cout << "Total sum of all errors : " << errtot << endl;

   return errtot;
}

double test_dbl8_complex ( int dim, int nvr, int pwr, int deg, int verbose )
{
   int *idx = new int[nvr];             // indices of variables in the monomial
   int *exp = new int[nvr];             // exponents of the variables
   int *expfac = new int[nvr];          // exponents of common factor
   int nbrfac;                          // number of common factors
   double *cffrehihihi = new double[deg+1]; // highest real coefficients parts
   double *cffrelohihi = new double[deg+1]; // second highest real parts
   double *cffrehilohi = new double[deg+1]; // third highest real parts 
   double *cffrelolohi = new double[deg+1]; // fourth highest real parts
   double *cffrehihilo = new double[deg+1]; // fourth lowest real parts
   double *cffrelohilo = new double[deg+1]; // third lowest real parts
   double *cffrehilolo = new double[deg+1]; // second lowest real parts 
   double *cffrelololo = new double[deg+1]; // lowest real parts
   double *cffimhihihi = new double[deg+1]; // highest imaginary parts
   double *cffimlohihi = new double[deg+1]; // second highest imaginary parts
   double *cffimhilohi = new double[deg+1]; // third highest imaginary parts 
   double *cffimlolohi = new double[deg+1]; // fourth highest imaginary parts
   double *cffimhihilo = new double[deg+1]; // fourth lowest imaginary parts
   double *cffimlohilo = new double[deg+1]; // third lowest imaginary parts
   double *cffimhilolo = new double[deg+1]; // second lowest imaginary parts 
   double *cffimlololo = new double[deg+1]; // lowest imaginary parts

   // The input are dim power series of degree deg,
   // the output are nvr+1 power series of degree deg,
   // for the evaluated and differentiated monomial.

   double **inputrehihihi = new double*[dim];
   double **inputrelohihi = new double*[dim];
   double **inputrehilohi = new double*[dim];
   double **inputrelolohi = new double*[dim];
   double **inputrehihilo = new double*[dim];
   double **inputrelohilo = new double*[dim];
   double **inputrehilolo = new double*[dim];
   double **inputrelololo = new double*[dim];
   double **inputimhihihi = new double*[dim];
   double **inputimlohihi = new double*[dim];
   double **inputimhilohi = new double*[dim];
   double **inputimlolohi = new double*[dim];
   double **inputimhihilo = new double*[dim];
   double **inputimlohilo = new double*[dim];
   double **inputimhilolo = new double*[dim];
   double **inputimlololo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      inputrehihihi[i] = new double[deg+1];
      inputrelohihi[i] = new double[deg+1];
      inputrehilohi[i] = new double[deg+1];
      inputrelolohi[i] = new double[deg+1];
      inputrehihilo[i] = new double[deg+1];
      inputrelohilo[i] = new double[deg+1];
      inputrehilolo[i] = new double[deg+1];
      inputrelololo[i] = new double[deg+1];
      inputimhihihi[i] = new double[deg+1];
      inputimlohihi[i] = new double[deg+1];
      inputimhilohi[i] = new double[deg+1];
      inputimlolohi[i] = new double[deg+1];
      inputimhihilo[i] = new double[deg+1];
      inputimlohilo[i] = new double[deg+1];
      inputimhilolo[i] = new double[deg+1];
      inputimlololo[i] = new double[deg+1];
   }
   double **outputrehihihi_h = new double*[dim+1];
   double **outputrelohihi_h = new double*[dim+1];
   double **outputrehilohi_h = new double*[dim+1];
   double **outputrelolohi_h = new double*[dim+1];
   double **outputrehihilo_h = new double*[dim+1];
   double **outputrelohilo_h = new double*[dim+1];
   double **outputrehilolo_h = new double*[dim+1];
   double **outputrelololo_h = new double*[dim+1];
   double **outputimhihihi_h = new double*[dim+1];
   double **outputimlohihi_h = new double*[dim+1];
   double **outputimhilohi_h = new double*[dim+1];
   double **outputimlolohi_h = new double*[dim+1];
   double **outputimhihilo_h = new double*[dim+1];
   double **outputimlohilo_h = new double*[dim+1];
   double **outputimhilolo_h = new double*[dim+1];
   double **outputimlololo_h = new double*[dim+1];

   for(int i=0; i<=dim; i++)
   {
      outputrehihihi_h[i] = new double[deg+1];
      outputrelohihi_h[i] = new double[deg+1];
      outputrehilohi_h[i] = new double[deg+1];
      outputrelolohi_h[i] = new double[deg+1];
      outputrehihilo_h[i] = new double[deg+1];
      outputrelohilo_h[i] = new double[deg+1];
      outputrehilolo_h[i] = new double[deg+1];
      outputrelololo_h[i] = new double[deg+1];
      outputimhihihi_h[i] = new double[deg+1];
      outputimlohihi_h[i] = new double[deg+1];
      outputimhilohi_h[i] = new double[deg+1];
      outputimlolohi_h[i] = new double[deg+1];
      outputimhihilo_h[i] = new double[deg+1];
      outputimlohilo_h[i] = new double[deg+1];
      outputimhilolo_h[i] = new double[deg+1];
      outputimlololo_h[i] = new double[deg+1];
   }
   double **outputrehihihi_d = new double*[dim+1];
   double **outputrelohihi_d = new double*[dim+1];
   double **outputrehilohi_d = new double*[dim+1];
   double **outputrelolohi_d = new double*[dim+1];
   double **outputrehihilo_d = new double*[dim+1];
   double **outputrelohilo_d = new double*[dim+1];
   double **outputrehilolo_d = new double*[dim+1];
   double **outputrelololo_d = new double*[dim+1];
   double **outputimhihihi_d = new double*[dim+1];
   double **outputimlohihi_d = new double*[dim+1];
   double **outputimhilohi_d = new double*[dim+1];
   double **outputimlolohi_d = new double*[dim+1];
   double **outputimhihilo_d = new double*[dim+1];
   double **outputimlohilo_d = new double*[dim+1];
   double **outputimhilolo_d = new double*[dim+1];
   double **outputimlololo_d = new double*[dim+1];

   for(int i=0; i<=dim; i++)
   {
      outputrehihihi_d[i] = new double[deg+1];
      outputrelohihi_d[i] = new double[deg+1];
      outputrehilohi_d[i] = new double[deg+1];
      outputrelolohi_d[i] = new double[deg+1];
      outputrehihilo_d[i] = new double[deg+1];
      outputrelohilo_d[i] = new double[deg+1];
      outputrehilolo_d[i] = new double[deg+1];
      outputrelololo_d[i] = new double[deg+1];
      outputimhihihi_d[i] = new double[deg+1];
      outputimlohihi_d[i] = new double[deg+1];
      outputimhilohi_d[i] = new double[deg+1];
      outputimlolohi_d[i] = new double[deg+1];
      outputimhihilo_d[i] = new double[deg+1];
      outputimlohilo_d[i] = new double[deg+1];
      outputimhilolo_d[i] = new double[deg+1];
      outputimlololo_d[i] = new double[deg+1];
   }
   bool fail = make_complex8_monomial
     (dim,nvr,pwr,deg,idx,exp,
      cffrehihihi,cffrelohihi,cffrehilohi,cffrelolohi,
      cffrehihilo,cffrelohilo,cffrehilolo,cffrelololo,
      cffimhihihi,cffimlohihi,cffimhilohi,cffimlolohi,
      cffimhihilo,cffimlohilo,cffimhilolo,cffimlololo);

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
         cout << cffrehihihi[i] << "  " << cffrelohihi[i] << endl;
         cout << cffrehilohi[i] << "  " << cffrelolohi[i] << endl;
         cout << cffrehihilo[i] << "  " << cffrelohilo[i] << endl;
         cout << cffrehilolo[i] << "  " << cffrelololo[i] << endl;
         cout << cffimhihihi[i] << "  " << cffimlohihi[i] << endl;
         cout << cffimhilohi[i] << "  " << cffimlolohi[i] << endl;
         cout << cffimhihilo[i] << "  " << cffimlohilo[i] << endl;
         cout << cffimhilolo[i] << "  " << cffimlololo[i] << endl;
      }
   }
   make_complex8_input(dim,deg,
      inputrehihihi,inputrelohihi,inputrehilohi,inputrelolohi,
      inputrehihilo,inputrelohilo,inputrehilolo,inputrelololo,
      inputimhihihi,inputimlohihi,inputimhilohi,inputimlolohi,
      inputimhihilo,inputimlohilo,inputimhilolo,inputimlololo);

   if(verbose > 0)
   {
      cout << "Random input series :" << endl;
      for(int i=0; i<dim; i++)
      {
         cout << "-> coefficients of series " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
         {
            cout << inputrehihihi[i][j] << "  "
                 << inputrelohihi[i][j] << endl;
            cout << inputrehilohi[i][j] << "  "
                 << inputrelolohi[i][j] << endl;
            cout << inputrehihilo[i][j] << "  "
                 << inputrelohilo[i][j] << endl;
            cout << inputrehilolo[i][j] << "  "
                 << inputrelololo[i][j] << endl;
            cout << inputimhihihi[i][j] << "  "
                 << inputimlohihi[i][j] << endl;
            cout << inputimhilohi[i][j] << "  "
                 << inputimlolohi[i][j] << endl;
            cout << inputimhihilo[i][j] << "  "
                 << inputimlohilo[i][j] << endl;
            cout << inputimhilolo[i][j] << "  "
                 << inputimlololo[i][j] << endl;
         }
      }
   }
   double errsum = 0.0;
   double errtot = 0.0;

   CPU_cmplx8_evaldiff
      (dim,nvr,deg,idx,
       cffrehihihi,cffrelohihi,cffrehilohi,cffrelolohi,
       cffrehihilo,cffrelohilo,cffrehilolo,cffrelololo,
       cffimhihihi,cffimlohihi,cffimhilohi,cffimlolohi,
       cffimhihilo,cffimlohilo,cffimhilolo,cffimlololo,
       inputrehihihi,inputrelohihi,inputrehilohi,inputrelolohi,
       inputrehihilo,inputrelohilo,inputrehilolo,inputrelololo,
       inputimhihihi,inputimlohihi,inputimhilohi,inputimlolohi,
       inputimhihilo,inputimlohilo,inputimhilolo,inputimlololo,
       outputrehihihi_h,outputrelohihi_h,outputrehilohi_h,outputrelolohi_h,
       outputrehihilo_h,outputrelohilo_h,outputrehilolo_h,outputrelololo_h,
       outputimhihihi_h,outputimlohihi_h,outputimhilohi_h,outputimlolohi_h,
       outputimhihilo_h,outputimlohilo_h,outputimhilolo_h,outputimlololo_h);

   if(verbose > 0)
   {
      cout << "The value of the product :" << endl;
      for(int i=0; i<=deg; i++)
      {
         cout << outputrehihihi_h[dim][i] << "  "
              << outputrelohihi_h[dim][i] << endl;
         cout << outputrehilohi_h[dim][i] << "  "
              << outputrelolohi_h[dim][i] << endl;
         cout << outputrehihilo_h[dim][i] << "  "
              << outputrelohilo_h[dim][i] << endl;
         cout << outputrehilolo_h[dim][i] << "  "
              << outputrelololo_h[dim][i] << endl;
         cout << outputimhihihi_h[dim][i] << "  "
              << outputimlohihi_h[dim][i] << endl;
         cout << outputimhilohi_h[dim][i] << "  "
              << outputimlolohi_h[dim][i] << endl;
         cout << outputimhihilo_h[dim][i] << "  "
              << outputimlohilo_h[dim][i] << endl;
         cout << outputimhilolo_h[dim][i] << "  "
              << outputimlololo_h[dim][i] << endl;
      }
   }
   if(nvr == dim) // the product of all input series equals one
   {
      for(int i=0; i<=deg; i++)
         errsum = errsum
                + abs(outputrehihihi_h[dim][i] - cffrehihihi[i])
                + abs(outputrelohihi_h[dim][i] - cffrelohihi[i])
                + abs(outputrehilohi_h[dim][i] - cffrehilohi[i])
                + abs(outputrelolohi_h[dim][i] - cffrelolohi[i])
                + abs(outputrehihilo_h[dim][i] - cffrehihilo[i])
                + abs(outputrelohilo_h[dim][i] - cffrelohilo[i])
                + abs(outputrehilolo_h[dim][i] - cffrehilolo[i])
                + abs(outputrelololo_h[dim][i] - cffrelololo[i])
                + abs(outputimhihihi_h[dim][i] - cffimhihihi[i])
                + abs(outputimlohihi_h[dim][i] - cffimlohihi[i])
                + abs(outputimhilohi_h[dim][i] - cffimhilohi[i])
                + abs(outputimlolohi_h[dim][i] - cffimlolohi[i])
                + abs(outputimhihilo_h[dim][i] - cffimhihilo[i])
                + abs(outputimlohilo_h[dim][i] - cffimlohilo[i])
                + abs(outputimhilolo_h[dim][i] - cffimhilolo[i])
                + abs(outputimlololo_h[dim][i] - cffimlololo[i]);

      if(verbose > 0) cout << "Coefficient error : " << errsum << endl;

      errtot += errsum;
   }
   if(nvr > 2)
   {
      GPU_cmplx8_evaldiff(deg+1,dim,nvr,deg,idx,
         cffrehihihi,cffrelohihi,cffrehilohi,cffrelolohi,
         cffrehihilo,cffrelohilo,cffrehilolo,cffrelololo,
         cffimhihihi,cffimlohihi,cffimhilohi,cffimlolohi,
         cffimhihilo,cffimlohilo,cffimhilolo,cffimlololo,
         inputrehihihi,inputrelohihi,inputrehilohi,inputrelolohi,
         inputrehihilo,inputrelohilo,inputrehilolo,inputrelololo,
         inputimhihihi,inputimlohihi,inputimhilohi,inputimlolohi,
         inputimhihilo,inputimlohilo,inputimhilolo,inputimlololo,
         outputrehihihi_d,outputrelohihi_d,outputrehilohi_d,outputrelolohi_d,
         outputrehihilo_d,outputrelohilo_d,outputrehilolo_d,outputrelololo_d,
         outputimhihihi_d,outputimlohihi_d,outputimhilohi_d,outputimlolohi_d,
         outputimhihilo_d,outputimlohilo_d,outputimhilolo_d,outputimlololo_d);

      errsum = 0.0;

      if(verbose > 0)
         cout << "The value of the product computed on the GPU :" << endl;

      for(int i=0; i<=deg; i++) 
      {
         if(verbose > 0)
         {
            cout << outputrehihihi_d[dim][i] << "  "
                 << outputrelohihi_d[dim][i] << endl;
            cout << outputrehilohi_d[dim][i] << "  "
                 << outputrelolohi_d[dim][i] << endl;
            cout << outputrehihilo_d[dim][i] << "  "
                 << outputrelohilo_d[dim][i] << endl;
            cout << outputrehilolo_d[dim][i] << "  "
                 << outputrelololo_d[dim][i] << endl;
            cout << outputimhihihi_d[dim][i] << "  "
                 << outputimlohihi_d[dim][i] << endl;
            cout << outputimhilohi_d[dim][i] << "  "
                 << outputimlolohi_d[dim][i] << endl;
            cout << outputimhihilo_d[dim][i] << "  "
                 << outputimlohilo_d[dim][i] << endl;
            cout << outputimhilolo_d[dim][i] << "  "
                 << outputimlololo_d[dim][i] << endl;
         }
         errsum = errsum
                + abs(outputrehihihi_h[dim][i] - outputrehihihi_d[dim][i])
                + abs(outputrelohihi_h[dim][i] - outputrelohihi_d[dim][i])
                + abs(outputrehilohi_h[dim][i] - outputrehilohi_d[dim][i])
                + abs(outputrelolohi_h[dim][i] - outputrelolohi_d[dim][i])
                + abs(outputrehihilo_h[dim][i] - outputrehihilo_d[dim][i])
                + abs(outputrelohilo_h[dim][i] - outputrelohilo_d[dim][i])
                + abs(outputrehilolo_h[dim][i] - outputrehilolo_d[dim][i])
                + abs(outputrelololo_h[dim][i] - outputrelololo_d[dim][i])
                + abs(outputimhihihi_h[dim][i] - outputimhihihi_d[dim][i])
                + abs(outputimlohihi_h[dim][i] - outputimlohihi_d[dim][i])
                + abs(outputimhilohi_h[dim][i] - outputimhilohi_d[dim][i])
                + abs(outputimlolohi_h[dim][i] - outputimlolohi_d[dim][i])
                + abs(outputimhihilo_h[dim][i] - outputimhihilo_d[dim][i])
                + abs(outputimlohilo_h[dim][i] - outputimlohilo_d[dim][i])
                + abs(outputimhilolo_h[dim][i] - outputimhilolo_d[dim][i])
                + abs(outputimlololo_h[dim][i] - outputimlololo_d[dim][i]);
      }
      if(verbose > 0) cout << "The sum of errors : " << errsum << endl;

      errtot += errsum;

      if(nvr == dim) // the product of all input series equals one
      {
         errsum = 0.0;
         for(int i=0; i<=deg; i++)
            errsum = errsum
                   + abs(outputrehihihi_d[dim][i] - cffrehihihi[i])
                   + abs(outputrelohihi_d[dim][i] - cffrelohihi[i])
                   + abs(outputrehilohi_d[dim][i] - cffrehilohi[i])
                   + abs(outputrelolohi_d[dim][i] - cffrelolohi[i])
                   + abs(outputrehihilo_d[dim][i] - cffrehihilo[i])
                   + abs(outputrelohilo_d[dim][i] - cffrelohilo[i])
                   + abs(outputrehilolo_d[dim][i] - cffrehilolo[i])
                   + abs(outputrelololo_d[dim][i] - cffrelololo[i])
                   + abs(outputimhihihi_d[dim][i] - cffimhihihi[i])
                   + abs(outputimlohihi_d[dim][i] - cffimlohihi[i])
                   + abs(outputimhilohi_d[dim][i] - cffimhilohi[i])
                   + abs(outputimlolohi_d[dim][i] - cffimlolohi[i])
                   + abs(outputimhihilo_d[dim][i] - cffimhihilo[i])
                   + abs(outputimlohilo_d[dim][i] - cffimlohilo[i])
                   + abs(outputimhilolo_d[dim][i] - cffimhilolo[i])
                   + abs(outputimlololo_d[dim][i] - cffimlololo[i]);

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
            cout << outputrehihihi_h[idx[k]][i] << "  "
                 << outputrelohihi_h[idx[k]][i] << endl;
            cout << outputrehilohi_h[idx[k]][i] << "  "
                 << outputrelolohi_h[idx[k]][i] << endl;
            cout << outputrehihilo_h[idx[k]][i] << "  "
                 << outputrelohilo_h[idx[k]][i] << endl;
            cout << outputrehilolo_h[idx[k]][i] << "  "
                 << outputrelololo_h[idx[k]][i] << endl;
            cout << outputimhihihi_h[idx[k]][i] << "  "
                 << outputimlohihi_h[idx[k]][i] << endl;
            cout << outputimhilohi_h[idx[k]][i] << "  "
                 << outputimlolohi_h[idx[k]][i] << endl;
            cout << outputimhihilo_h[idx[k]][i] << "  "
                 << outputimlohilo_h[idx[k]][i] << endl;
            cout << outputimhilolo_h[idx[k]][i] << "  "
                 << outputimlololo_h[idx[k]][i] << endl;
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
               cout << outputrehihihi_d[idx[k]][i] << "  "
                    << outputrelohihi_d[idx[k]][i] << endl;
               cout << outputrehilohi_d[idx[k]][i] << "  "
                    << outputrelolohi_d[idx[k]][i] << endl;
               cout << outputrehihilo_d[idx[k]][i] << "  "
                    << outputrelohilo_d[idx[k]][i] << endl;
               cout << outputrehilolo_d[idx[k]][i] << "  "
                    << outputrelololo_d[idx[k]][i] << endl;
               cout << outputimhihihi_d[idx[k]][i] << "  "
                    << outputimlohihi_d[idx[k]][i] << endl;
               cout << outputimhilohi_d[idx[k]][i] << "  "
                    << outputimlolohi_d[idx[k]][i] << endl;
               cout << outputimhihilo_d[idx[k]][i] << "  "
                    << outputimlohilo_d[idx[k]][i] << endl;
               cout << outputimhilolo_d[idx[k]][i] << "  "
                    << outputimlololo_d[idx[k]][i] << endl;
            }
            errsum = errsum
               + abs(outputrehihihi_h[idx[k]][i]
                   - outputrehihihi_d[idx[k]][i])
               + abs(outputrelohihi_h[idx[k]][i]
                   - outputrelohihi_d[idx[k]][i])
               + abs(outputrehilolo_h[idx[k]][i]
                   - outputrehilolo_d[idx[k]][i])
               + abs(outputrelololo_h[idx[k]][i]
                   - outputrelololo_d[idx[k]][i])
               + abs(outputrehihihi_h[idx[k]][i]
                   - outputrehihihi_d[idx[k]][i])
               + abs(outputrelohihi_h[idx[k]][i]
                   - outputrelohihi_d[idx[k]][i])
               + abs(outputrehilolo_h[idx[k]][i]
                   - outputrehilolo_d[idx[k]][i])
               + abs(outputrelololo_h[idx[k]][i]
                   - outputrelololo_d[idx[k]][i])
               + abs(outputimhihihi_h[idx[k]][i]
                   - outputimhihihi_d[idx[k]][i])
               + abs(outputimlohihi_h[idx[k]][i]
                   - outputimlohihi_d[idx[k]][i])
               + abs(outputimhilolo_h[idx[k]][i]
                   - outputimhilolo_d[idx[k]][i])
               + abs(outputimlololo_h[idx[k]][i]
                   - outputimlololo_d[idx[k]][i])
               + abs(outputimhihihi_h[idx[k]][i]
                   - outputimhihihi_d[idx[k]][i])
               + abs(outputimlohihi_h[idx[k]][i]
                   - outputimlohihi_d[idx[k]][i])
               + abs(outputimhilolo_h[idx[k]][i]
                   - outputimhilolo_d[idx[k]][i])
               + abs(outputimlololo_h[idx[k]][i]
                   - outputimlololo_d[idx[k]][i]);
         }
         if(verbose > 0) cout << "The sum of errors : " << errsum << endl;

         errtot += errsum;
      }
   }
   if(verbose > 0) cout << "Total sum of all errors : " << errtot << endl;

   return errtot;
}
