/* The file dbl4_polynomials_testers.cpp contains the definitions of
 * functions with prototypes in dbl4_polynomials_testers.h. */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <vector_types.h>
#include "random_polynomials.h"
#include "random4_monomials.h"
#include "random4_polynomials.h"
#include "convolution_jobs.h"
#include "addition_jobs.h"
#include "write_job_counts.h"
#include "write_gpu_timings.h"
#include "dbl4_polynomials_host.h"
#include "dbl4_polynomials_kernels.h"
#include "job_makers.h"
#include "dbl4_polynomials_testers.h"

using namespace std;

int main_dbl4_test_polynomial
 ( int seed, int dim, int nbr, int nva, int pwr, int deg, int vrblvl,
   double tol, bool jobrep, int mode )
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
   if(vrblvl > 0) cout << "  Seed used : " << seedused << endl;

   double realsum = test_dbl4_polynomial
                       (dim,nbr,nva,pwr,deg,vrblvl-1,jobrep,mode);
   double compsum = test_cmplx4_polynomial
                       (dim,nbr,nva,pwr,deg,vrblvl-1,jobrep,mode);

   int fail = int(realsum > tol) + int(compsum > tol);

   if(vrblvl > 0)
   {
      if(mode == 2)
      {
         cout << scientific << setprecision(2);
         cout << "Sum of all errors in quad double precision :" << endl;
         cout << "  on real data : " << realsum;
         if(realsum < tol)
            cout << "  pass." << endl;
         else
         {
            cout << " > " << tol;
            cout << "  fail!" << endl;
         }
         cout << "  on complex data : " << compsum;
         if(compsum < tol)
            cout << "  pass." << endl;
         else
         {
            cout << " > " << tol;
            cout << "  fail!" << endl;
         }
      }
      cout << "  Seed used : " <<  seedused << endl;
   }
   return fail;
}

int dbl4_make_input
 ( int dim, int nbr, int nva, int pwr, int deg,
   int *nvr, int **idx, int **exp,
   double **inputhihi, double **inputlohi,
   double **inputhilo, double **inputlolo,
   double *csthihi, double *cstlohi, double *csthilo, double *cstlolo,
   double **cffhihi, double **cfflohi, double **cffhilo, double **cfflolo,
   bool verbose )
{
   make_real4_input(dim,deg,inputhihi,inputlohi,inputhilo,inputlolo);

   if(verbose)
   {
      cout << scientific << setprecision(16);
      cout << "Random input series :" << endl;
      for(int i=0; i<dim; i++)
      {
         cout << "-> coefficients of series " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
            cout << inputhihi[i][j] << "  " << inputlohi[i][j] << endl
                 << inputhilo[i][j] << "  " << inputlolo[i][j] << endl;
      }
   }
   if(nva < 0) // read supports
   {
      read_supports(dim,nbr,nvr);
      for(int i=0; i<nbr; i++) idx[i] = new int[nvr[i]];
   }
   else if(nva == 0) // random supports
   {
      make_supports(dim,nbr,nvr); // random supports
      for(int i=0; i<nbr; i++) idx[i] = new int[nvr[i]];
   }
   else
   {
      for(int i=0; i<nbr; i++)
      {
         idx[i] = new int[nva];
         nvr[i] = nva;
      }
   }
   if(nva > 0)
   {
      if(nbr == dim)
         make_real4_cyclic
            (dim,nva,deg,idx,csthihi,cstlohi,csthilo,cstlolo,
                             cffhihi,cfflohi,cffhilo,cfflolo);
      else
         make_real4_products
            (dim,nbr,nva,deg,idx,csthihi,cstlohi,csthilo,cstlolo,
                                 cffhihi,cfflohi,cffhilo,cfflolo);
   }
   else
   {
      for(int i=0; i<nbr; i++) exp[i] = new int[nvr[i]];

      bool fail = make_real4_polynomial
                     (dim,nbr,pwr,deg,nvr,idx,exp,
                      csthihi,cstlohi,csthilo,cstlolo,
                      cffhihi,cfflohi,cffhilo,cfflolo);
   }
   if(verbose)
   {
      cout << "Coefficient series of the constant term :" << endl;
      for(int j=0; j<=deg; j++)
         cout << csthihi[j] << "  " << cstlohi[j] << endl
              << csthilo[j] << "  " << cstlolo[j] << endl;

      for(int i=0; i<nbr; i++)
      {
         cout << "Generated random monomial " << i << " :" << endl;
         cout << "   the indices :";
         for(int j=0; j<nvr[i]; j++) cout << " " << idx[i][j];
         cout << endl;
         if(nva == 0)
         {
            cout << " the exponents :";
            for(int j=0; j<nvr[i]; j++) cout << " " << exp[i][j];
            cout << endl;
         }
         cout << " coefficient series :" << endl;
         for(int j=0; j<=deg; j++)
            cout << cffhihi[i][j] << "  " << cfflohi[i][j] << endl
                 << cffhilo[i][j] << "  " << cfflolo[i][j] << endl;
      }
   }
   if(nva <= 0)
   {
      bool dup = duplicate_supports(dim,nbr,nvr,idx,verbose);
      if(dup)
      {
         cout << "Duplicates in support found." << endl;
         return 1;
      }
      else if(verbose)
         cout << "No duplicates in support found." << endl;
   }
   return 0;
}

int cmplx4_make_input
 ( int dim, int nbr, int nva, int pwr, int deg,
   int *nvr, int **idx, int **exp,
   double **inputrehihi, double **inputrelohi,
   double **inputrehilo, double **inputrelolo,
   double **inputimhihi, double **inputimlohi,
   double **inputimhilo, double **inputimlolo,
   double *cstrehihi, double *cstrelohi,
   double *cstrehilo, double *cstrelolo,
   double *cstimhihi, double *cstimlohi,
   double *cstimhilo, double *cstimlolo,
   double **cffrehihi, double **cffrelohi,
   double **cffrehilo, double **cffrelolo,
   double **cffimhihi, double **cffimlohi,
   double **cffimhilo, double **cffimlolo, bool verbose )
{
   make_complex4_input(dim,deg,
      inputrehihi,inputrelohi,inputrehilo,inputrelolo,
      inputimhihi,inputimlohi,inputimhilo,inputimlolo);

   if(verbose)
   {
      cout << scientific << setprecision(16);
      cout << "Random input series :" << endl;
      for(int i=0; i<dim; i++)
      {
         cout << "-> coefficients of series " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
         {
            cout << inputrehihi[i][j] << "  " << inputrelohi[i][j] << endl
                 << inputrehilo[i][j] << "  " << inputrelolo[i][j] << endl;
            cout << inputimhihi[i][j] << "  " << inputimlohi[i][j] << endl
                 << inputimhilo[i][j] << "  " << inputimlolo[i][j] << endl;
         }
      }
   }
   if(nva < 0) // read supports
   {
      read_supports(dim,nbr,nvr);
      for(int i=0; i<nbr; i++) idx[i] = new int[nvr[i]];
   }
   else if(nva == 0) // random supports
   {
      make_supports(dim,nbr,nvr);
      for(int i=0; i<nbr; i++) idx[i] = new int[nvr[i]];
   }
   else
   {
      for(int i=0; i<nbr; i++)
      {
         idx[i] = new int[nva];
         nvr[i] = nva;
      }
   }
   if(nva > 0)
   {
      if(nbr == dim)
         make_complex4_cyclic
            (dim,nva,deg,idx,
             cstrehihi,cstrelohi,cstrehilo,cstrelolo,
             cstimhihi,cstimlohi,cstimhilo,cstimlolo,
             cffrehihi,cffrelohi,cffrehilo,cffrelolo,
             cffimhihi,cffimlohi,cffimhilo,cffimlolo);
      else
         make_complex4_products
            (dim,nbr,nva,deg,idx,
             cstrehihi,cstrelohi,cstrehilo,cstrelolo,
             cstimhihi,cstimlohi,cstimhilo,cstimlolo,
             cffrehihi,cffrelohi,cffrehilo,cffrelolo,
             cffimhihi,cffimlohi,cffimhilo,cffimlolo);
   }
   else
   {
      for(int i=0; i<nbr; i++) exp[i] = new int[nvr[i]];

      bool fail = make_complex4_polynomial
                     (dim,nbr,pwr,deg,nvr,idx,exp,
                      cstrehihi,cstrelohi,cstrehilo,cstrelolo,
                      cstimhihi,cstimlohi,cstimhilo,cstimlolo,
                      cffrehihi,cffrelohi,cffrehilo,cffrelolo,
                      cffimhihi,cffimlohi,cffimhilo,cffimlolo);
   }
/*
   cout << "setting the constant term to zero ..." << endl;
   for(int j=0; j<=deg; j++)
   {
      cstrehihi[j] = 0.0; cstrelohi[j] = 0.0;
      cstrehilo[j] = 0.0; cstrelolo[j] = 0.0;
      cstimhihi[j] = 0.0; cstimlohi[j] = 0.0;
      cstimhilo[j] = 0.0; cstimlolo[j] = 0.0;
   }
 */
   if(verbose)
   {
      cout << "Coefficient series of the constant term :" << endl;
      for(int j=0; j<=deg; j++)
      {
         cout << cstrehihi[j] << "  " << cstrelohi[j] << endl
              << cstrehilo[j] << "  " << cstrelolo[j] << endl;
         cout << cstimhihi[j] << "  " << cstimlohi[j] << endl
              << cstimhilo[j] << "  " << cstimlolo[j] << endl;
      }
      for(int i=0; i<nbr; i++)
      {
         cout << "Generated random monomial " << i << " :" << endl;
         cout << "   the indices :";
         for(int j=0; j<nvr[i]; j++) cout << " " << idx[i][j];
         cout << endl;
         if(nva == 0)
         {
            cout << " the exponents :";
            for(int j=0; j<nvr[i]; j++) cout << " " << exp[i][j];
            cout << endl;
         }
         cout << " coefficient series :" << endl;
         for(int j=0; j<=deg; j++)
         {
            cout << cffrehihi[i][j] << "  " << cffrelohi[i][j] << endl
                 << cffrehilo[i][j] << "  " << cffrelolo[i][j] << endl;
            cout << cffimhihi[i][j] << "  " << cffimlohi[i][j] << endl
                 << cffimhilo[i][j] << "  " << cffimlolo[i][j] << endl;
         }
      }
    }
    if(nva <= 0)
    {
       bool dup = duplicate_supports(dim,nbr,nvr,idx,verbose);
       if(dup)
       {
          cout << "Duplicates in support found." << endl;
          return 1;
       }
       else if(verbose)
          cout << "No duplicates in support found." << endl;
    }
    return 0;
}

double dbl4_error_sum1
 ( int dim, int deg,
   double **resultshihi_h, double **resultslohi_h, 
   double **resultshilo_h, double **resultslolo_h,
   double **resultshihi_d, double **resultslohi_d,
   double **resultshilo_d, double **resultslolo_d, bool verbose )
{
   double err = 0.0;

   if(verbose) cout << "The value of the polynomial :" << endl;
   for(int i=0; i<=deg; i++)
   {
      if(verbose)
      {
         cout << resultshihi_h[dim][i] << "  "
              << resultslohi_h[dim][i] << endl
              << resultshilo_h[dim][i] << "  "
              << resultslolo_h[dim][i] << endl;
         cout << resultshihi_d[dim][i] << "  "
              << resultslohi_d[dim][i] << endl
              << resultshilo_d[dim][i] << "  "
              << resultslolo_d[dim][i] << endl;
      }
      err = err + abs(resultshihi_h[dim][i] - resultshihi_d[dim][i])
                + abs(resultslohi_h[dim][i] - resultslohi_d[dim][i])
                + abs(resultshilo_h[dim][i] - resultshilo_d[dim][i])
                + abs(resultslolo_h[dim][i] - resultslolo_d[dim][i]);
   }
   if(verbose) cout << "error : " << err << endl;

   double sumerr = err;

   for(int k=0; k<dim; k++)
   {
      if(verbose) cout << "Derivative " << k << " :" << endl;
      err = 0.0;
      for(int i=0; i<=deg; i++)
      {
         if(verbose)
         {
            cout << resultshihi_h[k][i] << "  "
                 << resultslohi_h[k][i] << endl
                 << resultshilo_h[k][i] << "  "
                 << resultslolo_h[k][i] << endl;
            cout << resultshihi_d[k][i] << "  "
                 << resultslohi_d[k][i] << endl
                 << resultshilo_d[k][i] << "  "
                 << resultslolo_d[k][i] << endl;
         }
         err = err 
                   + abs(resultshihi_h[k][i] - resultshihi_d[k][i])
                   + abs(resultslohi_h[k][i] - resultslohi_d[k][i])
                   + abs(resultshilo_h[k][i] - resultshilo_d[k][i])
                   + abs(resultslolo_h[k][i] - resultslolo_d[k][i]);
      }
      if(verbose) cout << "error : " << err << endl;
      sumerr = sumerr + err;
   }
   return sumerr;
}

double dbl4_error_sum
 ( int dim, int deg,
   double **results1hihi_h, double **results1lohi_h, 
   double **results1hilo_h, double **results1lolo_h,
   double **results2hihi_h, double **results2lohi_h,
   double **results2hilo_h, double **results2lolo_h,
   double **resultshihi_d, double **resultslohi_d,
   double **resultshilo_d, double **resultslolo_d, bool verbose )
{
   double err = 0.0;

   if(verbose) cout << "The value of the polynomial :" << endl;
   for(int i=0; i<=deg; i++)
   {
      if(verbose)
      {
         cout << results1hihi_h[dim][i] << "  "
              << results1lohi_h[dim][i] << endl
              << results1hilo_h[dim][i] << "  "
              << results1lolo_h[dim][i] << endl;
         cout << results2hihi_h[dim][i] << "  "
              << results2lohi_h[dim][i] << endl
              << results2hilo_h[dim][i] << "  "
              << results2lolo_h[dim][i] << endl;
         cout << resultshihi_d[dim][i] << "  "
              << resultslohi_d[dim][i] << endl
              << resultshilo_d[dim][i] << "  "
              << resultslolo_d[dim][i] << endl;
      }
      err = err + abs(results1hihi_h[dim][i] - results2hihi_h[dim][i])
                + abs(results1lohi_h[dim][i] - results2lohi_h[dim][i])
                + abs(results1hilo_h[dim][i] - results2hilo_h[dim][i])
                + abs(results1lolo_h[dim][i] - results2lolo_h[dim][i])
                + abs(results1hihi_h[dim][i] - resultshihi_d[dim][i])
                + abs(results1lohi_h[dim][i] - resultslohi_d[dim][i])
                + abs(results1hilo_h[dim][i] - resultshilo_d[dim][i])
                + abs(results1lolo_h[dim][i] - resultslolo_d[dim][i]);
   }
   if(verbose) cout << "error : " << err << endl;

   double sumerr = err;

   for(int k=0; k<dim; k++)
   {
      if(verbose) cout << "Derivative " << k << " :" << endl;
      err = 0.0;
      for(int i=0; i<=deg; i++)
      {
         if(verbose)
         {
            cout << results1hihi_h[k][i] << "  "
                 << results1lohi_h[k][i] << endl
                 << results1hilo_h[k][i] << "  "
                 << results1lolo_h[k][i] << endl;
            cout << results2hihi_h[k][i] << "  "
                 << results2lohi_h[k][i] << endl
                 << results2hilo_h[k][i] << "  "
                 << results2lolo_h[k][i] << endl;
            cout << resultshihi_d[k][i] << "  "
                 << resultslohi_d[k][i] << endl
                 << resultshilo_d[k][i] << "  "
                 << resultslolo_d[k][i] << endl;
         }
         err = err + abs(results1hihi_h[k][i] - results2hihi_h[k][i])
                   + abs(results1lohi_h[k][i] - results2lohi_h[k][i])
                   + abs(results1hilo_h[k][i] - results2hilo_h[k][i])
                   + abs(results1lolo_h[k][i] - results2lolo_h[k][i])
                   + abs(results1hihi_h[k][i] - resultshihi_d[k][i])
                   + abs(results1lohi_h[k][i] - resultslohi_d[k][i])
                   + abs(results1hilo_h[k][i] - resultshilo_d[k][i])
                   + abs(results1lolo_h[k][i] - resultslolo_d[k][i]);
      }
      if(verbose) cout << "error : " << err << endl;
      sumerr = sumerr + err;
   }
   return sumerr;
}

double cmplx4_error_sum1
 ( int dim, int deg,
   double **resultsrehihi_h, double **resultsrelohi_h,
   double **resultsrehilo_h, double **resultsrelolo_h,
   double **resultsimhihi_h, double **resultsimlohi_h,
   double **resultsimhilo_h, double **resultsimlolo_h,
   double **resultsrehihi_d, double **resultsrelohi_d,
   double **resultsrehilo_d, double **resultsrelolo_d,
   double **resultsimhihi_d, double **resultsimlohi_d,
   double **resultsimhilo_d, double **resultsimlolo_d, bool verbose )
{
   double err = 0.0;

   if(verbose) cout << "The value of the polynomial :" << endl;
   for(int i=0; i<=deg; i++)
   {
      if(verbose)
      {
         cout << resultsrehihi_h[dim][i] << "  "
              << resultsrelohi_h[dim][i] << endl
              << resultsrehilo_h[dim][i] << "  "
              << resultsrelolo_h[dim][i] << endl;
         cout << resultsimhihi_h[dim][i] << "  "
              << resultsimlohi_h[dim][i] << endl
              << resultsimhilo_h[dim][i] << "  "
              << resultsimlolo_h[dim][i] << endl;
         cout << resultsrehihi_d[dim][i] << "  "
              << resultsrelohi_d[dim][i] << endl
              << resultsrehilo_d[dim][i] << "  "
              << resultsrelolo_d[dim][i] << endl;
         cout << resultsimhihi_d[dim][i] << "  "
              << resultsimlohi_d[dim][i] << endl
              << resultsimhilo_d[dim][i] << "  "
              << resultsimlolo_d[dim][i] << endl;
      }
      err = err + abs(resultsrehihi_h[dim][i] - resultsrehihi_d[dim][i])
                + abs(resultsrelohi_h[dim][i] - resultsrelohi_d[dim][i])
                + abs(resultsrehilo_h[dim][i] - resultsrehilo_d[dim][i])
                + abs(resultsrelolo_h[dim][i] - resultsrelolo_d[dim][i])
                + abs(resultsimhihi_h[dim][i] - resultsimhihi_d[dim][i])
                + abs(resultsimlohi_h[dim][i] - resultsimlohi_d[dim][i])
                + abs(resultsimhilo_h[dim][i] - resultsimhilo_d[dim][i])
                + abs(resultsimlolo_h[dim][i] - resultsimlolo_d[dim][i]);
   }
   if(verbose) cout << "error : " << err << endl;

   double sumerr = err;

   for(int k=0; k<dim; k++)
   {
      if(verbose) cout << "Derivative " << k << " :" << endl;
      err = 0.0;
      for(int i=0; i<=deg; i++)
      {
         if(verbose)
         {
            cout << resultsrehihi_h[k][i] << "  "
                 << resultsrelohi_h[k][i] << endl
                 << resultsrehilo_h[k][i] << "  "
                 << resultsrelolo_h[k][i] << endl;
            cout << resultsimhihi_h[k][i] << "  "
                 << resultsimlohi_h[k][i] << endl
                 << resultsimhilo_h[k][i] << "  "
                 << resultsimlolo_h[k][i] << endl;
            cout << resultsrehihi_d[k][i] << "  "
                 << resultsrelohi_d[k][i] << endl
                 << resultsrehilo_d[k][i] << "  "
                 << resultsrelolo_d[k][i] << endl;
            cout << resultsimhihi_d[k][i] << "  "
                 << resultsimlohi_d[k][i] << endl
                 << resultsimhilo_d[k][i] << "  "
                 << resultsimlolo_d[k][i] << endl;
         }
         err = err + abs(resultsrehihi_h[k][i] - resultsrehihi_d[k][i])
                   + abs(resultsrelohi_h[k][i] - resultsrelohi_d[k][i])
                   + abs(resultsrehilo_h[k][i] - resultsrehilo_d[k][i])
                   + abs(resultsrelolo_h[k][i] - resultsrelolo_d[k][i])
                   + abs(resultsimhihi_h[k][i] - resultsimhihi_d[k][i])
                   + abs(resultsimlohi_h[k][i] - resultsimlohi_d[k][i])
                   + abs(resultsimhilo_h[k][i] - resultsimhilo_d[k][i])
                   + abs(resultsimlolo_h[k][i] - resultsimlolo_d[k][i]);
      }
      if(verbose) cout << "error : " << err << endl;
      sumerr = sumerr + err;
   }
   return sumerr;
}

double cmplx4_error_sum
 ( int dim, int deg,
   double **results1rehihi_h, double **results1relohi_h,
   double **results1rehilo_h, double **results1relolo_h,
   double **results1imhihi_h, double **results1imlohi_h,
   double **results1imhilo_h, double **results1imlolo_h,
   double **results2rehihi_h, double **results2relohi_h,
   double **results2rehilo_h, double **results2relolo_h,
   double **results2imhihi_h, double **results2imlohi_h,
   double **results2imhilo_h, double **results2imlolo_h,
   double **resultsrehihi_d, double **resultsrelohi_d,
   double **resultsrehilo_d, double **resultsrelolo_d,
   double **resultsimhihi_d, double **resultsimlohi_d,
   double **resultsimhilo_d, double **resultsimlolo_d, bool verbose )
{
   double err = 0.0;

   if(verbose) cout << "The value of the polynomial :" << endl;
   for(int i=0; i<=deg; i++)
   {
      if(verbose)
      {
         cout << results1rehihi_h[dim][i] << "  "
              << results1relohi_h[dim][i] << endl
              << results1rehilo_h[dim][i] << "  "
              << results1relolo_h[dim][i] << endl;
         cout << results1imhihi_h[dim][i] << "  "
              << results1imlohi_h[dim][i] << endl
              << results1imhilo_h[dim][i] << "  "
              << results1imlolo_h[dim][i] << endl;
         cout << results2rehihi_h[dim][i] << "  "
              << results2relohi_h[dim][i] << endl
              << results2rehilo_h[dim][i] << "  "
              << results2relolo_h[dim][i] << endl;
         cout << results2imhihi_h[dim][i] << "  "
              << results2imlohi_h[dim][i] << endl
              << results2imhilo_h[dim][i] << "  "
              << results2imlolo_h[dim][i] << endl;
         cout << resultsrehihi_d[dim][i] << "  "
              << resultsrelohi_d[dim][i] << endl
              << resultsrehilo_d[dim][i] << "  "
              << resultsrelolo_d[dim][i] << endl;
         cout << resultsimhihi_d[dim][i] << "  "
              << resultsimlohi_d[dim][i] << endl
              << resultsimhilo_d[dim][i] << "  "
              << resultsimlolo_d[dim][i] << endl;
      }
      err = err + abs(results1rehihi_h[dim][i] - results2rehihi_h[dim][i])
                + abs(results1relohi_h[dim][i] - results2relohi_h[dim][i])
                + abs(results1rehilo_h[dim][i] - results2rehilo_h[dim][i])
                + abs(results1relolo_h[dim][i] - results2relolo_h[dim][i])
                + abs(results1imhihi_h[dim][i] - results2imhihi_h[dim][i])
                + abs(results1imlohi_h[dim][i] - results2imlohi_h[dim][i])
                + abs(results1imhilo_h[dim][i] - results2imhilo_h[dim][i])
                + abs(results1imlolo_h[dim][i] - results2imlolo_h[dim][i])
                + abs(results1rehihi_h[dim][i] - resultsrehihi_d[dim][i])
                + abs(results1relohi_h[dim][i] - resultsrelohi_d[dim][i])
                + abs(results1rehilo_h[dim][i] - resultsrehilo_d[dim][i])
                + abs(results1relolo_h[dim][i] - resultsrelolo_d[dim][i])
                + abs(results1imhihi_h[dim][i] - resultsimhihi_d[dim][i])
                + abs(results1imlohi_h[dim][i] - resultsimlohi_d[dim][i])
                + abs(results1imhilo_h[dim][i] - resultsimhilo_d[dim][i])
                + abs(results1imlolo_h[dim][i] - resultsimlolo_d[dim][i]);
   }
   if(verbose) cout << "error : " << err << endl;

   double sumerr = err;

   for(int k=0; k<dim; k++)
   {
      if(verbose) cout << "Derivative " << k << " :" << endl;
      err = 0.0;
      for(int i=0; i<=deg; i++)
      {
         if(verbose)
         {
            cout << results1rehihi_h[k][i] << "  "
                 << results1relohi_h[k][i] << endl
                 << results1rehilo_h[k][i] << "  "
                 << results1relolo_h[k][i] << endl;
            cout << results1imhihi_h[k][i] << "  "
                 << results1imlohi_h[k][i] << endl
                 << results1imhilo_h[k][i] << "  "
                 << results1imlolo_h[k][i] << endl;
            cout << results2rehihi_h[k][i] << "  "
                 << results2relohi_h[k][i] << endl
                 << results2rehilo_h[k][i] << "  "
                 << results2relolo_h[k][i] << endl;
            cout << results2imhihi_h[k][i] << "  "
                 << results2imlohi_h[k][i] << endl
                 << results2imhilo_h[k][i] << "  "
                 << results2imlolo_h[k][i] << endl;
            cout << resultsrehihi_d[k][i] << "  "
                 << resultsrelohi_d[k][i] << endl
                 << resultsrehilo_d[k][i] << "  "
                 << resultsrelolo_d[k][i] << endl;
            cout << resultsimhihi_d[k][i] << "  "
                 << resultsimlohi_d[k][i] << endl
                 << resultsimhilo_d[k][i] << "  "
                 << resultsimlolo_d[k][i] << endl;
         }
         err = err + abs(results1rehihi_h[k][i] - results2rehihi_h[k][i])
                   + abs(results1relohi_h[k][i] - results2relohi_h[k][i])
                   + abs(results1rehilo_h[k][i] - results2rehilo_h[k][i])
                   + abs(results1relolo_h[k][i] - results2relolo_h[k][i])
                   + abs(results1imhihi_h[k][i] - results2imhihi_h[k][i])
                   + abs(results1imlohi_h[k][i] - results2imlohi_h[k][i])
                   + abs(results1imhilo_h[k][i] - results2imhilo_h[k][i])
                   + abs(results1imlolo_h[k][i] - results2imlolo_h[k][i])
                   + abs(results1rehihi_h[k][i] - resultsrehihi_d[k][i])
                   + abs(results1relohi_h[k][i] - resultsrelohi_d[k][i])
                   + abs(results1rehilo_h[k][i] - resultsrehilo_d[k][i])
                   + abs(results1relolo_h[k][i] - resultsrelolo_d[k][i])
                   + abs(results1imhihi_h[k][i] - resultsimhihi_d[k][i])
                   + abs(results1imlohi_h[k][i] - resultsimlohi_d[k][i])
                   + abs(results1imhilo_h[k][i] - resultsimhilo_d[k][i])
                   + abs(results1imlolo_h[k][i] - resultsimlolo_d[k][i]);
      }
      if(verbose) cout << "error : " << err << endl;
      sumerr = sumerr + err;
   }
   return sumerr;
}

double test_dbl4_polynomial
 ( int dim, int nbr, int nva, int pwr, int deg, int verbose, bool jobrep,
   int mode )
{
   cout << "*** testing the real arithmetic on real data ***" << endl;

   if(nbr < 1)
      return 0.0;
   else
   {
      double **inputhihi = new double*[dim]; // dim series of degree deg
      for(int i=0; i<dim; i++) inputhihi[i] = new double[deg+1];
      double **inputlohi = new double*[dim];
      for(int i=0; i<dim; i++) inputlohi[i] = new double[deg+1];
      double **inputhilo = new double*[dim];
      for(int i=0; i<dim; i++) inputhilo[i] = new double[deg+1];
      double **inputlolo = new double*[dim];
      for(int i=0; i<dim; i++) inputlolo[i] = new double[deg+1];
      // The output are dim+1 power series of degree deg
      // for the evaluated and differentiated polynomial.
      double **output1hihi_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output1hihi_h[i] = new double[deg+1];
      double **output1lohi_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output1lohi_h[i] = new double[deg+1];
      double **output1hilo_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output1hilo_h[i] = new double[deg+1];
      double **output1lolo_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output1lolo_h[i] = new double[deg+1];
      double **output2hihi_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output2hihi_h[i] = new double[deg+1];
      double **output2lohi_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output2lohi_h[i] = new double[deg+1];
      double **output2hilo_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output2hilo_h[i] = new double[deg+1];
      double **output2lolo_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output2lolo_h[i] = new double[deg+1];
      double **outputhihi_d = new double*[dim+1];
      for(int i=0; i<=dim; i++) outputhihi_d[i] = new double[deg+1];
      double **outputlohi_d = new double*[dim+1];
      for(int i=0; i<=dim; i++) outputlohi_d[i] = new double[deg+1];
      double **outputhilo_d = new double*[dim+1];
      for(int i=0; i<=dim; i++) outputhilo_d[i] = new double[deg+1];
      double **outputlolo_d = new double*[dim+1];
      for(int i=0; i<=dim; i++) outputlolo_d[i] = new double[deg+1];

      double *csthihi = new double[deg+1]; // constant coefficient series
      double *cstlohi = new double[deg+1];
      double *csthilo = new double[deg+1];
      double *cstlolo = new double[deg+1];
      double **cffhihi = new double*[nbr]; // coefficient series of terms
      for(int i=0; i<nbr; i++) cffhihi[i] = new double[deg+1];
      double **cfflohi = new double*[nbr];
      for(int i=0; i<nbr; i++) cfflohi[i] = new double[deg+1];
      double **cffhilo = new double*[nbr];
      for(int i=0; i<nbr; i++) cffhilo[i] = new double[deg+1];
      double **cfflolo = new double*[nbr];
      for(int i=0; i<nbr; i++) cfflolo[i] = new double[deg+1];
      int *nvr = new int[nbr]; // number of variables in each monomial
      int **idx = new int*[nbr];  // indices of variables in monomials
      int **exp = new int*[nbr];  // exponents of the variables

      bool vrb = (verbose > 1);

      int fail = dbl4_make_input(dim,nbr,nva,pwr,deg,nvr,idx,exp,
                                 inputhihi,inputlohi,inputhilo,inputlolo,
                                 csthihi,cstlohi,csthilo,cstlolo,
                                 cffhihi,cfflohi,cffhilo,cfflolo,vrb);
      if(fail == 1)
      {
          cout << "Duplicates in support, returning 0 as error." << endl;
          return 0.0;
      }
      ConvolutionJobs cnvjobs(dim);

      cnvjobs.make(nbr,nvr,idx,vrb);

      if(vrb)
      {
         write_convolution_counts(cnvjobs);

         for(int k=0; k<cnvjobs.get_depth(); k++)
         {
            cout << "jobs at layer " << k << " :" << endl;
            for(int i=0; i<cnvjobs.get_layer_count(k); i++)
               cout << cnvjobs.get_job(k,i) << endl;
         }
         cout << endl;
      }
      AdditionJobs addjobs(dim,nbr);

      addjobs.make(nbr,nvr,idx,vrb);

      if(vrb)
      {
         cout << "The differential indices :" << endl;
         for(int i=0; i<dim; i++)
         {
            cout << "variable " << i << " :";
            for(int j=0; j<=addjobs.get_differential_count(i); j++)
               cout << " " << addjobs.get_differential_index(i,j);
            cout << endl;
         }
         write_addition_counts(addjobs);

         for(int k=0; k<addjobs.get_depth(); k++)
         {
            cout << "jobs at layer " << k << " :" << endl;
            for(int i=0; i<addjobs.get_layer_count(k); i++)
               cout << addjobs.get_job(k,i) << endl;
         }
      }
      double timelapsec1_h,timelapsec2_h;
      double cnvlapms,addlapms,timelapms_d,walltimes_d;

      if((mode == 1) || (mode == 2))
      {
         if(vrb) cout << "Computing without convolution jobs ..." << endl;
         CPU_dbl4_poly_evaldiff
            (dim,nbr,deg,nvr,idx,csthihi,cstlohi,csthilo,cstlolo,
             cffhihi,cfflohi,cffhilo,cfflolo,
             inputhihi,inputlohi,inputhilo,inputlolo,
             output1hihi_h,output1lohi_h,output1hilo_h,output1lolo_h,
             &timelapsec1_h,vrb);
         if(vrb) cout << "Computing with convolution jobs ..." << endl;
         CPU_dbl4_poly_evaldiffjobs
            (dim,nbr,deg,nvr,idx,csthihi,cstlohi,csthilo,cstlolo,
             cffhihi,cfflohi,cffhilo,cfflolo,
             inputhihi,inputlohi,inputhilo,inputlolo,
             output2hihi_h,output2lohi_h,output2hilo_h,output2lolo_h,
             cnvjobs,addjobs,&timelapsec2_h,vrb);
      }
      if((mode == 0) || (mode == 2))
      {
         if(vrb) cout << "Computing on the device ..." << endl;
         GPU_dbl4_poly_evaldiff
            (deg+1,dim,nbr,deg,nvr,idx,csthihi,cstlohi,csthilo,cstlolo,
             cffhihi,cfflohi,cffhilo,cfflolo,
             inputhihi,inputlohi,inputhilo,inputlolo,
             outputhihi_d,outputlohi_d,outputhilo_d,outputlolo_d,
             cnvjobs,addjobs,&cnvlapms,&addlapms,&timelapms_d,
             &walltimes_d,vrb);
      }
      double sumerr = 0.0;
      if(mode == 2)
      {
         sumerr = dbl4_error_sum(dim,deg,
                     output1hihi_h,output1lohi_h,
                     output1hilo_h,output1lolo_h,
                     output2hihi_h,output2lohi_h,
                     output2hilo_h,output2lolo_h,
                     outputhihi_d,outputlohi_d,
                     outputhilo_d,outputlolo_d,vrb);
         cout << "sum of all errors " << sumerr << endl;
      }
      if(verbose > 0)
      {
         if(jobrep)
         {
            cout << "dimension : " << dim << endl;
            if(nva > 0)
            {
               cout << "number of variables per monomial : " << nva << endl;
            }
            cout << "number of monomials : " << nbr << endl;
            write_convolution_counts(cnvjobs);
            write_addition_counts(addjobs);
            write_operation_counts(deg,cnvjobs,addjobs);
         }
         if((mode == 1) || (mode == 2))
         {
            cout << fixed << setprecision(3);
            cout << "Elapsed CPU time (Linux), Wall time (Windows) : "
                 << endl;
            cout << "  (1) without jobs : " << timelapsec1_h << " seconds,"
                 << endl;
            cout << "  (2) cnv/add jobs : " << timelapsec2_h << " seconds."
                 << endl;
         }
         if((mode == 0) || (mode == 2))
         {
            cout << fixed << setprecision(2);
            cout << "Time spent by convolution kernels : "
                 << cnvlapms << " milliseconds." << endl;
            cout << "Time spent by addition kernels    : "
                 << addlapms << " milliseconds." << endl;
            cout << "Time spent by all kernels         : "
                 << timelapms_d << " milliseconds." << endl;
            cout << "Total wall clock computation time : ";
            cout << fixed << setprecision(3) << walltimes_d
                 << " seconds." << endl;
            cout << scientific << setprecision(16);
         }
      }
      return sumerr;
   }
}

double test_cmplx4_polynomial
 ( int dim, int nbr, int nva, int pwr, int deg, int verbose, bool jobrep,
   int mode )
{
   cout << "*** testing the complex arithmetic on complex data ***" << endl;

   if(nbr < 1)
      return 0.0;
   else
   {
      double **inputrehihi = new double*[dim]; // dim series of degree deg
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
      // The output are dim+1 power series of degree deg
      // for the evaluated and differentiated polynomial.
      double **output1rehihi_h = new double*[dim+1];
      double **output1relohi_h = new double*[dim+1];
      double **output1rehilo_h = new double*[dim+1];
      double **output1relolo_h = new double*[dim+1];
      double **output1imhihi_h = new double*[dim+1];
      double **output1imlohi_h = new double*[dim+1];
      double **output1imhilo_h = new double*[dim+1];
      double **output1imlolo_h = new double*[dim+1];
      for(int i=0; i<=dim; i++)
      {
         output1rehihi_h[i] = new double[deg+1];
         output1relohi_h[i] = new double[deg+1];
         output1rehilo_h[i] = new double[deg+1];
         output1relolo_h[i] = new double[deg+1];
         output1imhihi_h[i] = new double[deg+1];
         output1imlohi_h[i] = new double[deg+1];
         output1imhilo_h[i] = new double[deg+1];
         output1imlolo_h[i] = new double[deg+1];
      }
      double **output2rehihi_h = new double*[dim+1];
      double **output2relohi_h = new double*[dim+1];
      double **output2rehilo_h = new double*[dim+1];
      double **output2relolo_h = new double*[dim+1];
      double **output2imhihi_h = new double*[dim+1];
      double **output2imlohi_h = new double*[dim+1];
      double **output2imhilo_h = new double*[dim+1];
      double **output2imlolo_h = new double*[dim+1];
      for(int i=0; i<=dim; i++)
      {
         output2rehihi_h[i] = new double[deg+1];
         output2relohi_h[i] = new double[deg+1];
         output2rehilo_h[i] = new double[deg+1];
         output2relolo_h[i] = new double[deg+1];
         output2imhihi_h[i] = new double[deg+1];
         output2imlohi_h[i] = new double[deg+1];
         output2imhilo_h[i] = new double[deg+1];
         output2imlolo_h[i] = new double[deg+1];
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
      double *cstrehihi = new double[deg+1]; // constant coefficient series
      double *cstrelohi = new double[deg+1];
      double *cstrehilo = new double[deg+1];
      double *cstrelolo = new double[deg+1];
      double *cstimhihi = new double[deg+1];
      double *cstimlohi = new double[deg+1];
      double *cstimhilo = new double[deg+1];
      double *cstimlolo = new double[deg+1];
      double **cffrehihi = new double*[nbr]; // coefficient series of terms
      double **cffrelohi = new double*[nbr];
      double **cffrehilo = new double*[nbr];
      double **cffrelolo = new double*[nbr];
      double **cffimhihi = new double*[nbr];
      double **cffimlohi = new double*[nbr];
      double **cffimhilo = new double*[nbr];
      double **cffimlolo = new double*[nbr];
      for(int i=0; i<nbr; i++)
      {
         cffrehihi[i] = new double[deg+1];
         cffrelohi[i] = new double[deg+1];
         cffrehilo[i] = new double[deg+1];
         cffrelolo[i] = new double[deg+1];
         cffimhihi[i] = new double[deg+1];
         cffimlohi[i] = new double[deg+1];
         cffimhilo[i] = new double[deg+1];
         cffimlolo[i] = new double[deg+1];
      }
      int *nvr = new int[nbr]; // number of variables in each monomial
      int **idx = new int*[nbr];  // indices of variables in monomials
      int **exp = new int*[nbr];  // exponents of the variables

      bool vrb = (verbose > 1);

      int fail = cmplx4_make_input(dim,nbr,nva,pwr,deg,nvr,idx,exp,
                    inputrehihi,inputrelohi,inputrehilo,inputrelolo,
                    inputimhihi,inputimlohi,inputimhilo,inputimlolo,
                    cstrehihi,cstrelohi,cstrehilo,cstrelolo,
                    cstimhihi,cstimlohi,cstimhilo,cstimlolo,
                    cffrehihi,cffrelohi,cffrehilo,cffrelolo,
                    cffimhihi,cffimlohi,cffimhilo,cffimlolo,vrb);
      if(fail == 1)
      {
          cout << "Duplicates in support, returning 0 as error." << endl;
          return 0.0;
      }
      ConvolutionJobs cnvjobs(dim);
      AdditionJobs addjobs(dim,nbr);

      make_all_jobs(dim,nbr,nvr,idx,&cnvjobs,&addjobs,vrb);

      ComplexConvolutionJobs cmplxcnvjobs(dim);
      ComplexIncrementJobs cmplxincjobs(cmplxcnvjobs,vrb);
      ComplexAdditionJobs cmplxaddjobs(dim,nbr);

      make_all_complex_jobs
         (dim,nbr,nvr,idx,&cmplxcnvjobs,&cmplxincjobs,&cmplxaddjobs,vrb);

      double timelapsec1_h,timelapsec2_h;
      double cnvlapms,addlapms,timelapms_d,walltimes_d;

      if((mode == 1) || (mode == 2))
      {
         if(vrb) cout << "Computing without convolution jobs ..." << endl;
         CPU_cmplx4_poly_evaldiff
            (dim,nbr,deg,nvr,idx,
             cstrehihi,cstrelohi,cstrehilo,cstrelolo,
             cstimhihi,cstimlohi,cstimhilo,cstimlolo,
             cffrehihi,cffrelohi,cffrehilo,cffrelolo,
             cffimhihi,cffimlohi,cffimhilo,cffimlolo,
             inputrehihi,inputrelohi,inputrehilo,inputrelolo,
             inputimhihi,inputimlohi,inputimhilo,inputimlolo,
             output1rehihi_h,output1relohi_h,output1rehilo_h,output1relolo_h,
             output1imhihi_h,output1imlohi_h,output1imhilo_h,output1imlolo_h,
             &timelapsec1_h,vrb);
         if(vrb) cout << "Computing with convolution jobs ..." << endl;
         CPU_cmplx4_poly_evaldiffjobs
            (dim,nbr,deg,nvr,idx,
             cstrehihi,cstrelohi,cstrehilo,cstrelolo,
             cstimhihi,cstimlohi,cstimhilo,cstimlolo,
             cffrehihi,cffrelohi,cffrehilo,cffrelolo,
             cffimhihi,cffimlohi,cffimhilo,cffimlolo,
             inputrehihi,inputrelohi,inputrehilo,inputrelolo,
             inputimhihi,inputimlohi,inputimhilo,inputimlolo,
             output2rehihi_h,output2relohi_h,output2rehilo_h,output2relolo_h,
             output2imhihi_h,output2imlohi_h,output2imhilo_h,output2imlolo_h,
             cnvjobs,addjobs,&timelapsec2_h,vrb);
      }
      if((mode == 0) || (mode == 2))
      {
         if(vrb) cout << "Computing on the device ..." << endl;
/*
         GPU_cmplx4_poly_evaldiff
            (deg+1,dim,nbr,deg,nvr,idx,
             cstrehihi,cstrelohi,cstrehilo,cstrelolo,
             cstimhihi,cstimlohi,cstimhilo,cstimlolo,
             cffrehihi,cffrelohi,cffrehilo,cffrelolo,
             cffimhihi,cffimlohi,cffimhilo,cffimlolo,
             inputrehihi,inputrelohi,inputrehilo,inputrelolo,
             inputimhihi,inputimlohi,inputimhilo,inputimlolo,
             outputrehihi_d,outputrelohi_d,outputrehilo_d,outputrelolo_d,
             outputimhihi_d,outputimlohi_d,outputimhilo_d,outputimlolo_d,
             cnvjobs,addjobs,&cnvlapms,&addlapms,&timelapms_d,
             &walltimes_d,vrb);
 */
         GPU_cmplx4vectorized_poly_evaldiff
            (deg+1,dim,nbr,deg,nvr,idx,
             cstrehihi,cstrelohi,cstrehilo,cstrelolo,
             cstimhihi,cstimlohi,cstimhilo,cstimlolo,
             cffrehihi,cffrelohi,cffrehilo,cffrelolo,
             cffimhihi,cffimlohi,cffimhilo,cffimlolo,
             inputrehihi,inputrelohi,inputrehilo,inputrelolo,
             inputimhihi,inputimlohi,inputimhilo,inputimlolo,
             outputrehihi_d,outputrelohi_d,outputrehilo_d,outputrelolo_d,
             outputimhihi_d,outputimlohi_d,outputimhilo_d,outputimlolo_d,
             cmplxcnvjobs,cmplxincjobs,cmplxaddjobs,
             &cnvlapms,&addlapms,&timelapms_d,&walltimes_d,vrb);
      }
      double sumerr = 0.0;
      if(mode == 2)
      {
         sumerr = cmplx4_error_sum(dim,deg,
                     output1rehihi_h,output1relohi_h,
                     output1rehilo_h,output1relolo_h,
                     output1imhihi_h,output1imlohi_h,
                     output1imhilo_h,output1imlolo_h,
                     output1rehihi_h,output1relohi_h,
                     output1rehilo_h,output1relolo_h,
                     output1imhihi_h,output1imlohi_h,
                     output1imhilo_h,output1imlolo_h,
           /*        output2rehihi_h,output2relohi_h,
                     output2rehilo_h,output2relolo_h,
                     output2imhihi_h,output2imlohi_h,
                     output2imhilo_h,output2imlolo_h, */
                     outputrehihi_d,outputrelohi_d,
                     outputrehilo_d,outputrelolo_d,
                     outputimhihi_d,outputimlohi_d,
                     outputimhilo_d,outputimlolo_d,vrb);
         cout << "sum of all errors " << sumerr << endl;
      }
      if(verbose > 0)
      {
         if(jobrep) write_jobs_report(dim,nva,nbr,deg,cnvjobs,addjobs);
         if((mode == 1) || (mode == 2))
            write_CPU_timings(timelapsec1_h,timelapsec2_h);
         if((mode == 0) || (mode == 2))
            write_GPU_timings(cnvlapms,addlapms,timelapms_d,walltimes_d);
      } 
      return sumerr;
   }
}

int test_dbl4_sequence
 ( int seed, int dim, int nva, int nbr, int pwr, int vrblvl,
   bool jobrep, int mode )
{
   double const tol = 1.0e-56;

   int deg = 0;
   cout << "---> running in quad double precision for degree 0 ..." << endl;
   int fail = main_dbl4_test_polynomial
                 (seed,dim,nbr,nva,pwr,deg,vrblvl,tol,jobrep,mode);
   deg = 8;
   cout << "---> running for degree 8 ..." << endl;
   fail += main_dbl4_test_polynomial
              (seed,dim,nbr,nva,pwr,deg,vrblvl,tol,false,mode);
   deg = 15;
   cout << "---> running for degree 15 ..." << endl;
   fail += main_dbl4_test_polynomial
              (seed,dim,nbr,nva,pwr,deg,vrblvl,tol,false,mode);
   deg = 31;
   cout << "---> running for degree 31 ..." << endl;
   fail += main_dbl4_test_polynomial
              (seed,dim,nbr,nva,pwr,deg,vrblvl,tol,false,mode);
   deg = 63;
   cout << "---> running for degree 63 ..." << endl;
   fail += main_dbl4_test_polynomial
              (seed,dim,nbr,nva,pwr,deg,vrblvl,tol,false,mode);
   deg = 95;
   cout << "---> running for degree 95 ..." << endl;
   fail += main_dbl4_test_polynomial
              (seed,dim,nbr,nva,pwr,deg,vrblvl,tol,false,mode);
   deg = 127;
   cout << "---> running for degree 127 ..." << endl;
   fail += main_dbl4_test_polynomial
              (seed,dim,nbr,nva,pwr,deg,vrblvl,tol,false,mode);
   deg = 152;
   cout << "---> running for degree 152 ..." << endl;
   fail += main_dbl4_test_polynomial
              (seed,dim,nbr,nva,pwr,deg,vrblvl,tol,false,mode);
   deg = 159;
   cout << "---> running for degree 159 ..." << endl;
   fail += main_dbl4_test_polynomial
              (seed,dim,nbr,nva,pwr,deg,vrblvl,tol,false,mode);
   deg = 191;
   cout << "---> running for degree 191 ..." << endl;
   fail += main_dbl4_test_polynomial
              (seed,dim,nbr,nva,pwr,deg,vrblvl,tol,false,mode);

   if(mode == 2)
   {
      if(fail == 0)
         cout << "All tests passed." << endl;
      else
         cout << "Number of failed tests : " << fail << endl;
   }
   return 0;
}
