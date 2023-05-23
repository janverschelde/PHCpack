/* The file dbl8_polynomials_testers.cpp contains the definitions of
 * functions with prototypes in dbl8_polynomials_testers.h. */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <vector_types.h>
#include "random_polynomials.h"
#include "random8_monomials.h"
#include "random8_polynomials.h"
#include "convolution_jobs.h"
#include "addition_jobs.h"
#include "write_job_counts.h"
#include "dbl8_polynomials_host.h"
#include "dbl8_polynomials_kernels.h"
#include "dbl8_polynomials_testers.h"

using namespace std;

int main_dbl8_test_polynomial
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

   double realsum = test_dbl8_real_polynomial
                       (dim,nbr,nva,pwr,deg,vrblvl-1,jobrep,mode);

   int fail = int(realsum > tol);

   if(vrblvl > 0)
   {
      if(mode == 2)
      {
         cout << scientific << setprecision(2);
         cout << "Sum of all errors in octo double precision :" << endl;
         cout << "  on real data : " << realsum;
         if(realsum < tol)
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

void dbl8_make_input
 ( int dim, int nbr, int nva, int pwr, int deg,
   int *nvr, int **idx, int **exp,
   double **inputhihihi, double **inputlohihi,
   double **inputhilohi, double **inputlolohi,
   double **inputhihilo, double **inputlohilo,
   double **inputhilolo, double **inputlololo,
   double *csthihihi, double *cstlohihi,
   double *csthilohi, double *cstlolohi,
   double *csthihilo, double *cstlohilo,
   double *csthilolo, double *cstlololo,
   double **cffhihihi, double **cfflohihi,
   double **cffhilohi, double **cfflolohi,
   double **cffhihilo, double **cfflohilo,
   double **cffhilolo, double **cfflololo, bool verbose )
{
   make_real8_input(dim,deg,
      inputhihihi,inputlohihi,inputhilohi,inputlolohi,
      inputhihilo,inputlohilo,inputhilolo,inputlololo);

   if(verbose)
   {
      cout << scientific << setprecision(16);
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
   if(nva == 0) // random supports
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
         make_real8_cyclic
            (dim,nva,deg,idx,
             csthihihi,cstlohihi,csthilohi,cstlolohi,
             csthihilo,cstlohilo,csthilolo,cstlololo,
             cffhihihi,cfflohihi,cffhilohi,cfflolohi,
             cffhihilo,cfflohilo,cffhilolo,cfflololo);
      else
         make_real8_products
            (dim,nbr,nva,deg,idx,
             csthihihi,cstlohihi,csthilohi,cstlolohi,
             csthihilo,cstlohilo,csthilolo,cstlololo,
             cffhihihi,cfflohihi,cffhilohi,cfflolohi,
             cffhihilo,cfflohilo,cffhilolo,cfflololo);
   }
   else
   {
      for(int i=0; i<nbr; i++) exp[i] = new int[nvr[i]];

      bool fail = make_real8_polynomial
                     (dim,nbr,pwr,deg,nvr,idx,exp,
                      csthihihi,cstlohihi,csthilohi,cstlolohi,
                      csthihilo,cstlohilo,csthilolo,cstlololo,
                      cffhihihi,cfflohihi,cffhilohi,cffhilolo,
                      cffhihilo,cfflohilo,cffhilolo,cfflololo);
   }
   if(verbose)
   {
      cout << "Coefficient series of the constant term :" << endl;
      for(int j=0; j<=deg; j++)
      {
         cout << csthihihi[j] << " " << cstlohihi[j] << endl;
         cout << csthilohi[j] << " " << cstlolohi[j] << endl;
         cout << csthihilo[j] << " " << cstlohilo[j] << endl;
         cout << csthilolo[j] << " " << cstlololo[j] << endl;
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
            cout << cffhihihi[i][j] << " " << cfflohihi[i][j] << endl;
            cout << cffhilohi[i][j] << " " << cfflolohi[i][j] << endl;
            cout << cffhihilo[i][j] << " " << cfflohilo[i][j] << endl;
            cout << cffhilolo[i][j] << " " << cfflololo[i][j] << endl;
         }
      }
   }
   if(nva == 0)
   {
      bool dup = duplicate_supports(dim,nbr,nvr,idx,verbose);
      if(dup)
         cout << "Duplicate supports found." << endl;
      else if(verbose)
         cout << "No duplicate supports found." << endl;
   }
}

void cmplx8_make_input
 ( int dim, int nbr, int nva, int pwr, int deg,
   int *nvr, int **idx, int **exp,
   double **inputrehihihi, double **inputrelohihi,
   double **inputrehilohi, double **inputrelolohi,
   double **inputrehihilo, double **inputrelohilo,
   double **inputrehilolo, double **inputrelololo,
   double **inputimhihihi, double **inputimlohihi,
   double **inputimhilohi, double **inputimlolohi,
   double **inputimhihilo, double **inputimlohilo,
   double **inputimhilolo, double **inputimlololo,
   double *cstrehihihi, double *cstrelohihi,
   double *cstrehilohi, double *cstrelolohi,
   double *cstrehihilo, double *cstrelohilo,
   double *cstrehilolo, double *cstrelololo,
   double *cstimhihihi, double *cstimlohihi,
   double *cstimhilohi, double *cstimlolohi,
   double *cstimhihilo, double *cstimlohilo,
   double *cstimhilolo, double *cstimlololo,
   double **cffrehihihi, double **cffrelohihi,
   double **cffrehilohi, double **cffrelolohi,
   double **cffrehihilo, double **cffrelohilo,
   double **cffrehilolo, double **cffrelololo,
   double **cffimhihihi, double **cffimlohihi,
   double **cffimhilohi, double **cffimlolohi,
   double **cffimhihilo, double **cffimlohilo,
   double **cffimhilolo, double **cffimlololo, bool verbose )
{
   make_complex8_input(dim,deg,
      inputrehihihi,inputrelohihi,inputrehilohi,inputrelolohi,
      inputrehihilo,inputrelohilo,inputrehilolo,inputrelololo,
      inputimhihihi,inputimlohihi,inputimhilohi,inputimlolohi,
      inputimhihilo,inputimlohilo,inputimhilolo,inputimlololo);

   if(verbose)
   {
      cout << scientific << setprecision(16);
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
   if(nva == 0) // random supports
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
         make_cmplx8_cyclic
            (dim,nva,deg,idx,
             cstrehihihi,cstrelohihi,cstrehilohi,cstrelolohi,
             cstrehihilo,cstrelohilo,cstrehilolo,cstrelololo,
             cstimhihihi,cstimlohihi,cstimhilohi,cstimlolohi,
             cstimhihilo,cstimlohilo,cstimhilolo,cstimlololo,
             cffrehihihi,cffrelohihi,cffrehilohi,cffrelolohi,
             cffrehihilo,cffrelohilo,cffrehilolo,cffrelololo,
             cffimhihihi,cffimlohihi,cffimhilohi,cffimlolohi,
             cffimhihilo,cffimlohilo,cffimhilolo,cffimlololo);
      else
         make_cmplx8_products
            (dim,nbr,nva,deg,idx,
             cstrehihihi,cstrelohihi,cstrehilohi,cstrelolohi,
             cstrehihilo,cstrelohilo,cstrehilolo,cstrelololo,
             cstimhihihi,cstimlohihi,cstimhilohi,cstimlolohi,
             cstimhihilo,cstimlohilo,cstimhilolo,cstimlololo,
             cffrehihihi,cffrelohihi,cffrehilohi,cffrelolohi,
             cffrehihilo,cffrelohilo,cffrehilolo,cffrelololo,
             cffimhihihi,cffimlohihi,cffimhilohi,cffimlolohi,
             cffimhihilo,cffimlohilo,cffimhilolo,cffimlololo);
   }
   else
   {
      for(int i=0; i<nbr; i++) exp[i] = new int[nvr[i]];

      bool fail = make_cmplx8_polynomial
                     (dim,nbr,pwr,deg,nvr,idx,exp,
                      cstrehihihi,cstrelohihi,cstrehilohi,cstrelolohi,
                      cstrehihilo,cstrelohilo,cstrehilolo,cstrelololo,
                      cstimhihihi,cstimlohihi,cstimhilohi,cstimlolohi,
                      cstimhihilo,cstimlohilo,cstimhilolo,cstimlololo,
                      cffrehihihi,cffrelohihi,cffrehilohi,cffrehilolo,
                      cffrehihilo,cffrelohilo,cffrehilolo,cffrelololo,
                      cffimhihihi,cffimlohihi,cffimhilohi,cffimhilolo,
                      cffimhihilo,cffimlohilo,cffimhilolo,cffimlololo);
   }
   if(verbose)
   {
      cout << "Coefficient series of the constant term :" << endl;
      for(int j=0; j<=deg; j++)
      {
         cout << cstrehihihi[j] << " " << cstrelohihi[j] << endl;
         cout << cstrehilohi[j] << " " << cstrelolohi[j] << endl;
         cout << cstrehihilo[j] << " " << cstrelohilo[j] << endl;
         cout << cstrehilolo[j] << " " << cstrelololo[j] << endl;
         cout << cstimhihihi[j] << " " << cstimlohihi[j] << endl;
         cout << cstimhilohi[j] << " " << cstimlolohi[j] << endl;
         cout << cstimhihilo[j] << " " << cstimlohilo[j] << endl;
         cout << cstimhilolo[j] << " " << cstimlololo[j] << endl;
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
            cout << cffrehihihi[i][j] << " " << cffrelohihi[i][j] << endl;
            cout << cffrehilohi[i][j] << " " << cffrelolohi[i][j] << endl;
            cout << cffrehihilo[i][j] << " " << cffrelohilo[i][j] << endl;
            cout << cffrehilolo[i][j] << " " << cffrelololo[i][j] << endl;
            cout << cffimhihihi[i][j] << " " << cffimlohihi[i][j] << endl;
            cout << cffimhilohi[i][j] << " " << cffimlolohi[i][j] << endl;
            cout << cffimhihilo[i][j] << " " << cffimlohilo[i][j] << endl;
            cout << cffimhilolo[i][j] << " " << cffimlololo[i][j] << endl;
         }
      }
   }
   if(nva == 0)
   {
      bool dup = duplicate_supports(dim,nbr,nvr,idx,verbose);
      if(dup)
         cout << "Duplicate supports found." << endl;
      else if(verbose)
         cout << "No duplicate supports found." << endl;
   }
}

double dbl8_error_sum
 ( int dim, int deg,
   double **results1hihihi_h, double **results1hilohi_h, 
   double **results1hihilo_h, double **results1hilolo_h,
   double **results1lohihi_h, double **results1lolohi_h, 
   double **results1lohilo_h, double **results1lololo_h,
   double **results2hihihi_h, double **results2hilohi_h,
   double **results2hihilo_h, double **results2hilolo_h,
   double **results2lohihi_h, double **results2lolohi_h,
   double **results2lohilo_h, double **results2lololo_h,
   double **resultshihihi_d, double **resultshilohi_d,
   double **resultshihilo_d, double **resultshilolo_d,
   double **resultslohihi_d, double **resultslolohi_d,
   double **resultslohilo_d, double **resultslololo_d, bool verbose )
{
   double err = 0.0;

   if(verbose) cout << "The value of the polynomial :" << endl;
   for(int i=0; i<=deg; i++)
   {
      if(verbose)
      {
         cout << results1hihihi_h[dim][i] << "  "
              << results1hilohi_h[dim][i] << endl
              << results1hihilo_h[dim][i] << "  "
              << results1hilolo_h[dim][i] << endl;
         cout << results1lohihi_h[dim][i] << "  "
              << results1lolohi_h[dim][i] << endl
              << results1lohilo_h[dim][i] << "  "
              << results1lololo_h[dim][i] << endl;
         cout << results2hihihi_h[dim][i] << "  "
              << results2hilohi_h[dim][i] << endl
              << results2hihilo_h[dim][i] << "  "
              << results2hilolo_h[dim][i] << endl;
         cout << results2lohihi_h[dim][i] << "  "
              << results2lolohi_h[dim][i] << endl
              << results2lohilo_h[dim][i] << "  "
              << results2lololo_h[dim][i] << endl;
         cout << resultshihihi_d[dim][i] << "  "
              << resultshilohi_d[dim][i] << endl
              << resultshihilo_d[dim][i] << "  "
              << resultshilolo_d[dim][i] << endl;
         cout << resultslohihi_d[dim][i] << "  "
              << resultslolohi_d[dim][i] << endl
              << resultslohilo_d[dim][i] << "  "
              << resultslololo_d[dim][i] << endl;
      }
      err = err + abs(results1hihihi_h[dim][i] - results2hihihi_h[dim][i])
                + abs(results1hilohi_h[dim][i] - results2hilohi_h[dim][i])
                + abs(results1hihilo_h[dim][i] - results2hihilo_h[dim][i])
                + abs(results1hilolo_h[dim][i] - results2hilolo_h[dim][i])
                + abs(results1lohihi_h[dim][i] - results2lohihi_h[dim][i])
                + abs(results1lolohi_h[dim][i] - results2lolohi_h[dim][i])
                + abs(results1lohilo_h[dim][i] - results2lohilo_h[dim][i])
                + abs(results1lololo_h[dim][i] - results2lololo_h[dim][i])
                + abs(results1hihihi_h[dim][i] - resultshihihi_d[dim][i])
                + abs(results1hilohi_h[dim][i] - resultshilohi_d[dim][i])
                + abs(results1hihilo_h[dim][i] - resultshihilo_d[dim][i])
                + abs(results1hilolo_h[dim][i] - resultshilolo_d[dim][i])
                + abs(results1lohihi_h[dim][i] - resultslohihi_d[dim][i])
                + abs(results1lolohi_h[dim][i] - resultslolohi_d[dim][i])
                + abs(results1lohilo_h[dim][i] - resultslohilo_d[dim][i])
                + abs(results1lololo_h[dim][i] - resultslololo_d[dim][i]);
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
            cout << results1hihihi_h[k][i] << "  "
                 << results1hilohi_h[k][i] << endl
                 << results1hihilo_h[k][i] << "  "
                 << results1hilolo_h[k][i] << endl;
            cout << results1lohihi_h[k][i] << "  "
                 << results1lolohi_h[k][i] << endl
                 << results1lohilo_h[k][i] << "  "
                 << results1lololo_h[k][i] << endl;
            cout << results2hihihi_h[k][i] << "  "
                 << results2hilohi_h[k][i] << endl
                 << results2hihilo_h[k][i] << "  "
                 << results2hilolo_h[k][i] << endl;
            cout << results2lohihi_h[k][i] << "  "
                 << results2lolohi_h[k][i] << endl
                 << results2lohilo_h[k][i] << "  "
                 << results2lololo_h[k][i] << endl;
            cout << resultshihihi_d[k][i] << "  "
                 << resultshilohi_d[k][i] << endl
                 << resultshihilo_d[k][i] << "  "
                 << resultshilolo_d[k][i] << endl;
            cout << resultslohihi_d[k][i] << "  "
                 << resultslolohi_d[k][i] << endl
                 << resultslohilo_d[k][i] << "  "
                 << resultslololo_d[k][i] << endl;
         }
         err = err + abs(results1hihihi_h[k][i] - results2hihihi_h[k][i])
                   + abs(results1hilohi_h[k][i] - results2hilohi_h[k][i])
                   + abs(results1hihilo_h[k][i] - results2hihilo_h[k][i])
                   + abs(results1hilolo_h[k][i] - results2hilolo_h[k][i])
                   + abs(results1lohihi_h[k][i] - results2lohihi_h[k][i])
                   + abs(results1lolohi_h[k][i] - results2lolohi_h[k][i])
                   + abs(results1lohilo_h[k][i] - results2lohilo_h[k][i])
                   + abs(results1lololo_h[k][i] - results2lololo_h[k][i])
                   + abs(results1hihihi_h[k][i] - resultshihihi_d[k][i])
                   + abs(results1hilohi_h[k][i] - resultshilohi_d[k][i])
                   + abs(results1hihilo_h[k][i] - resultshihilo_d[k][i])
                   + abs(results1hilolo_h[k][i] - resultshilolo_d[k][i])
                   + abs(results1lohihi_h[k][i] - resultslohihi_d[k][i])
                   + abs(results1lolohi_h[k][i] - resultslolohi_d[k][i])
                   + abs(results1lohilo_h[k][i] - resultslohilo_d[k][i])
                   + abs(results1lololo_h[k][i] - resultslololo_d[k][i]);
      }
      if(verbose) cout << "error : " << err << endl;
      sumerr = sumerr + err;
   }
   return sumerr;
}

double cmplx8_error_sum
 ( int dim, int deg,
   double **resultsrehihihi_h, double **resultsrelohihi_h, 
   double **resultsrehilohi_h, double **resultsrelolohi_h,
   double **resultsrehihilo_h, double **resultsrelohilo_h, 
   double **resultsrehilolo_h, double **resultsrelololo_h,
   double **resultsimhihihi_h, double **resultsimlohihi_h, 
   double **resultsimhilohi_h, double **resultsimlolohi_h,
   double **resultsimhihilo_h, double **resultsimlohilo_h, 
   double **resultsimhilolo_h, double **resultsimlololo_h,
   double **resultsrehihihi_d, double **resultsrelohihi_d,
   double **resultsrehilohi_d, double **resultsrelolohi_d,
   double **resultsrehihilo_d, double **resultsrelohilo_d,
   double **resultsrehilolo_d, double **resultsrelololo_d,
   double **resultsimhihihi_d, double **resultsimlohihi_d,
   double **resultsimhilohi_d, double **resultsimlolohi_d,
   double **resultsimhihilo_d, double **resultsimlohilo_d,
   double **resultsimhilolo_d, double **resultsimlololo_d, bool verbose )
{
   double err = 0.0;

   if(verbose) cout << "The value of the polynomial :" << endl;
   for(int i=0; i<=deg; i++)
   {
      if(verbose)
      {
         cout << resultsrehihihi_h[dim][i] << "  "
              << resultsrelohihi_h[dim][i] << endl
              << resultsrehilohi_h[dim][i] << "  "
              << resultsrehihilo_h[dim][i] << endl;
         cout << resultsrehihilo_h[dim][i] << "  "
              << resultsrelohilo_h[dim][i] << endl
              << resultsrehilolo_h[dim][i] << "  "
              << resultsrelololo_h[dim][i] << endl;
         cout << resultsimhihihi_h[dim][i] << "  "
              << resultsimlohihi_h[dim][i] << endl
              << resultsimhilohi_h[dim][i] << "  "
              << resultsimlolohi_h[dim][i] << endl;
         cout << resultsimhihilo_h[dim][i] << "  "
              << resultsimlohilo_h[dim][i] << endl
              << resultsimhilolo_h[dim][i] << "  "
              << resultsimlololo_h[dim][i] << endl;
         cout << resultsrehihihi_d[dim][i] << "  "
              << resultsrelohihi_d[dim][i] << endl
              << resultsrehilohi_d[dim][i] << "  "
              << resultsrelolohi_d[dim][i] << endl;
         cout << resultsrehihilo_d[dim][i] << "  "
              << resultsrelohilo_d[dim][i] << endl
              << resultsrehilolo_d[dim][i] << "  "
              << resultsrelololo_d[dim][i] << endl;
         cout << resultsimhihihi_d[dim][i] << "  "
              << resultsimlohihi_d[dim][i] << endl
              << resultsimhilohi_d[dim][i] << "  "
              << resultsimlolohi_d[dim][i] << endl;
         cout << resultsimhihilo_d[dim][i] << "  "
              << resultsimlohilo_d[dim][i] << endl
              << resultsimhilolo_d[dim][i] << "  "
              << resultsimlololo_d[dim][i] << endl;
      }
      err = err + abs(resultsrehihihi_h[dim][i] - resultsrehihihi_d[dim][i])
                + abs(resultsrelohihi_h[dim][i] - resultsrelohihi_d[dim][i])
                + abs(resultsrehilohi_h[dim][i] - resultsrehilohi_d[dim][i])
                + abs(resultsrelolohi_h[dim][i] - resultsrelolohi_d[dim][i])
                + abs(resultsrehihilo_h[dim][i] - resultsrehihilo_d[dim][i])
                + abs(resultsrelohilo_h[dim][i] - resultsrelohilo_d[dim][i])
                + abs(resultsrehilolo_h[dim][i] - resultsrehilolo_d[dim][i])
                + abs(resultsrelololo_h[dim][i] - resultsrelololo_d[dim][i])
                + abs(resultsimhihihi_h[dim][i] - resultsimhihihi_d[dim][i])
                + abs(resultsimlohihi_h[dim][i] - resultsimlohihi_d[dim][i])
                + abs(resultsimhilohi_h[dim][i] - resultsimhilohi_d[dim][i])
                + abs(resultsimlolohi_h[dim][i] - resultsimlolohi_d[dim][i])
                + abs(resultsimhihilo_h[dim][i] - resultsimhihilo_d[dim][i])
                + abs(resultsimlohilo_h[dim][i] - resultsimlohilo_d[dim][i])
                + abs(resultsimhilolo_h[dim][i] - resultsimhilolo_d[dim][i])
                + abs(resultsimlololo_h[dim][i] - resultsimlololo_d[dim][i]);
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
            cout << resultsrehihihi_h[k][i] << "  "
                 << resultsrelohihi_h[k][i] << endl
                 << resultsrehilohi_h[k][i] << "  "
                 << resultsrelolohi_h[k][i] << endl;
            cout << resultsrehihilo_h[k][i] << "  "
                 << resultsrelohilo_h[k][i] << endl
                 << resultsrehilolo_h[k][i] << "  "
                 << resultsrelololo_h[k][i] << endl;
            cout << resultsimhihihi_h[k][i] << "  "
                 << resultsimlohihi_h[k][i] << endl
                 << resultsimhilohi_h[k][i] << "  "
                 << resultsimlolohi_h[k][i] << endl;
            cout << resultsimhihilo_h[k][i] << "  "
                 << resultsimlohilo_h[k][i] << endl
                 << resultsimhilolo_h[k][i] << "  "
                 << resultsimlololo_h[k][i] << endl;
            cout << resultsrehihihi_d[k][i] << "  "
                 << resultsrelohihi_d[k][i] << endl
                 << resultsrehilohi_d[k][i] << "  "
                 << resultsrelolohi_d[k][i] << endl;
            cout << resultsrehihilo_d[k][i] << "  "
                 << resultsrelohilo_d[k][i] << endl
                 << resultsrehilolo_d[k][i] << "  "
                 << resultsrelololo_d[k][i] << endl;
            cout << resultsimhihihi_d[k][i] << "  "
                 << resultsimlohihi_d[k][i] << endl
                 << resultsimhilohi_d[k][i] << "  "
                 << resultsimlolohi_d[k][i] << endl;
            cout << resultsimhihilo_d[k][i] << "  "
                 << resultsimlohilo_d[k][i] << endl
                 << resultsimhilolo_d[k][i] << "  "
                 << resultsimlololo_d[k][i] << endl;
         }
         err = err + abs(resultsrehihihi_h[k][i] - resultsrehihihi_d[k][i])
                   + abs(resultsrehilohi_h[k][i] - resultsrehilohi_d[k][i])
                   + abs(resultsrehihilo_h[k][i] - resultsrehihilo_d[k][i])
                   + abs(resultsrehilolo_h[k][i] - resultsrehilolo_d[k][i])
                   + abs(resultsrelohihi_h[k][i] - resultsrelohihi_d[k][i])
                   + abs(resultsrelolohi_h[k][i] - resultsrelolohi_d[k][i])
                   + abs(resultsrelohilo_h[k][i] - resultsrelohilo_d[k][i])
                   + abs(resultsrelololo_h[k][i] - resultsrelololo_d[k][i])
                   + abs(resultsimhihihi_h[k][i] - resultsimhihihi_d[k][i])
                   + abs(resultsimhilohi_h[k][i] - resultsimhilohi_d[k][i])
                   + abs(resultsimhihilo_h[k][i] - resultsimhihilo_d[k][i])
                   + abs(resultsimhilolo_h[k][i] - resultsimhilolo_d[k][i])
                   + abs(resultsimlohihi_h[k][i] - resultsimlohihi_d[k][i])
                   + abs(resultsimlolohi_h[k][i] - resultsimlolohi_d[k][i])
                   + abs(resultsimlohilo_h[k][i] - resultsimlohilo_d[k][i])
                   + abs(resultsimlololo_h[k][i] - resultsimlololo_d[k][i]);
      }
      if(verbose) cout << "error : " << err << endl;
      sumerr = sumerr + err;
   }
   return sumerr;
}

double test_dbl8_real_polynomial
 ( int dim, int nbr, int nva, int pwr, int deg, int verbose, bool jobrep,
   int mode )
{
   if(nbr < 1)
      return 0.0;
   else
   {
      double **inputhihihi = new double*[dim]; // dim series of degree deg
      for(int i=0; i<dim; i++) inputhihihi[i] = new double[deg+1];
      double **inputhilohi = new double*[dim];
      for(int i=0; i<dim; i++) inputhilohi[i] = new double[deg+1];
      double **inputhihilo = new double*[dim];
      for(int i=0; i<dim; i++) inputhihilo[i] = new double[deg+1];
      double **inputhilolo = new double*[dim];
      for(int i=0; i<dim; i++) inputhilolo[i] = new double[deg+1];
      double **inputlohihi = new double*[dim];
      for(int i=0; i<dim; i++) inputlohihi[i] = new double[deg+1];
      double **inputlolohi = new double*[dim];
      for(int i=0; i<dim; i++) inputlolohi[i] = new double[deg+1];
      double **inputlohilo = new double*[dim];
      for(int i=0; i<dim; i++) inputlohilo[i] = new double[deg+1];
      double **inputlololo = new double*[dim];
      for(int i=0; i<dim; i++) inputlololo[i] = new double[deg+1];
      // The output are dim+1 power series of degree deg
      // for the evaluated and differentiated polynomial.
      double **output1hihihi_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output1hihihi_h[i] = new double[deg+1];
      double **output1hilohi_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output1hilohi_h[i] = new double[deg+1];
      double **output1hihilo_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output1hihilo_h[i] = new double[deg+1];
      double **output1hilolo_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output1hilolo_h[i] = new double[deg+1];
      double **output1lohihi_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output1lohihi_h[i] = new double[deg+1];
      double **output1lolohi_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output1lolohi_h[i] = new double[deg+1];
      double **output1lohilo_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output1lohilo_h[i] = new double[deg+1];
      double **output1lololo_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output1lololo_h[i] = new double[deg+1];
      double **output2hihihi_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output2hihihi_h[i] = new double[deg+1];
      double **output2hilohi_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output2hilohi_h[i] = new double[deg+1];
      double **output2hihilo_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output2hihilo_h[i] = new double[deg+1];
      double **output2hilolo_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output2hilolo_h[i] = new double[deg+1];
      double **output2lohihi_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output2lohihi_h[i] = new double[deg+1];
      double **output2lolohi_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output2lolohi_h[i] = new double[deg+1];
      double **output2lohilo_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output2lohilo_h[i] = new double[deg+1];
      double **output2lololo_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output2lololo_h[i] = new double[deg+1];
      double **outputhihihi_d = new double*[dim+1];
      for(int i=0; i<=dim; i++) outputhihihi_d[i] = new double[deg+1];
      double **outputhilohi_d = new double*[dim+1];
      for(int i=0; i<=dim; i++) outputhilohi_d[i] = new double[deg+1];
      double **outputhihilo_d = new double*[dim+1];
      for(int i=0; i<=dim; i++) outputhihilo_d[i] = new double[deg+1];
      double **outputhilolo_d = new double*[dim+1];
      for(int i=0; i<=dim; i++) outputhilolo_d[i] = new double[deg+1];
      double **outputlohihi_d = new double*[dim+1];
      for(int i=0; i<=dim; i++) outputlohihi_d[i] = new double[deg+1];
      double **outputlolohi_d = new double*[dim+1];
      for(int i=0; i<=dim; i++) outputlolohi_d[i] = new double[deg+1];
      double **outputlohilo_d = new double*[dim+1];
      for(int i=0; i<=dim; i++) outputlohilo_d[i] = new double[deg+1];
      double **outputlololo_d = new double*[dim+1];
      for(int i=0; i<=dim; i++) outputlololo_d[i] = new double[deg+1];

      double *csthihihi = new double[deg+1]; // constant coefficient series
      double *csthilohi = new double[deg+1];
      double *csthihilo = new double[deg+1];
      double *csthilolo = new double[deg+1];
      double *cstlohihi = new double[deg+1]; 
      double *cstlolohi = new double[deg+1];
      double *cstlohilo = new double[deg+1];
      double *cstlololo = new double[deg+1];
      double **cffhihihi = new double*[nbr]; // coefficient series of terms
      for(int i=0; i<nbr; i++) cffhihihi[i] = new double[deg+1];
      double **cffhilohi = new double*[nbr];
      for(int i=0; i<nbr; i++) cffhilohi[i] = new double[deg+1];
      double **cffhihilo = new double*[nbr];
      for(int i=0; i<nbr; i++) cffhihilo[i] = new double[deg+1];
      double **cffhilolo = new double*[nbr];
      for(int i=0; i<nbr; i++) cffhilolo[i] = new double[deg+1];
      double **cfflohihi = new double*[nbr];
      for(int i=0; i<nbr; i++) cfflohihi[i] = new double[deg+1];
      double **cfflolohi = new double*[nbr];
      for(int i=0; i<nbr; i++) cfflolohi[i] = new double[deg+1];
      double **cfflohilo = new double*[nbr];
      for(int i=0; i<nbr; i++) cfflohilo[i] = new double[deg+1];
      double **cfflololo = new double*[nbr];
      for(int i=0; i<nbr; i++) cfflololo[i] = new double[deg+1];
      int *nvr = new int[nbr]; // number of variables in each monomial
      int **idx = new int*[nbr];  // indices of variables in monomials
      int **exp = new int*[nbr];  // exponents of the variables

      bool vrb = (verbose > 1);

      dbl8_make_input(dim,nbr,nva,pwr,deg,nvr,idx,exp,
                      inputhihihi,inputhilohi,inputhihilo,inputhilolo,
                      inputlohihi,inputlolohi,inputlohilo,inputlololo,
                      csthihihi,csthilohi,csthihilo,csthilolo,
                      cstlohihi,cstlolohi,cstlohilo,cstlololo,
                      cffhihihi,cffhilohi,cffhihilo,cffhilolo,
                      cfflohihi,cfflolohi,cfflohilo,cfflololo,vrb);

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
         CPU_dbl8_poly_evaldiff
            (dim,nbr,deg,nvr,idx,
             csthihihi,csthilohi,csthihilo,csthilolo,
             cstlohihi,cstlolohi,cstlohilo,cstlololo,
             cffhihihi,cffhilohi,cffhihilo,cffhilolo,
             cfflohihi,cfflolohi,cfflohilo,cfflololo,
             inputhihihi,inputhilohi,inputhihilo,inputhilolo,
             inputlohihi,inputlolohi,inputlohilo,inputlololo,
             output1hihihi_h,output1hilohi_h,output1hihilo_h,output1hilolo_h,
             output1lohihi_h,output1lolohi_h,output1lohilo_h,output1lololo_h,
             &timelapsec1_h,vrb);
         if(vrb) cout << "Computing with convolution jobs ..." << endl;
         CPU_dbl8_poly_evaldiffjobs
            (dim,nbr,deg,nvr,idx,
             csthihihi,csthilohi,csthihilo,csthilolo,
             cstlohihi,cstlolohi,cstlohilo,cstlololo,
             cffhihihi,cffhilohi,cffhihilo,cffhilolo,
             cfflohihi,cfflolohi,cfflohilo,cfflololo,
             inputhihihi,inputhilohi,inputhihilo,inputhilolo,
             inputlohihi,inputlolohi,inputlohilo,inputlololo,
             output2hihihi_h,output2hilohi_h,output2hihilo_h,output2hilolo_h,
             output2lohihi_h,output2lolohi_h,output2lohilo_h,output2lololo_h,
             cnvjobs,addjobs,&timelapsec2_h,vrb);
      }
      if((mode == 0) || (mode == 2))
      {
         if(vrb) cout << "Computing on the device ..." << endl;
         GPU_dbl8_poly_evaldiff
            (deg+1,dim,nbr,deg,nvr,idx,
             csthihihi,csthilohi,csthihilo,csthilolo,
             cstlohihi,cstlolohi,cstlohilo,cstlololo,
             cffhihihi,cffhilohi,cffhihilo,cffhilolo,
             cfflohihi,cfflolohi,cfflohilo,cfflololo,
             inputhihihi,inputhilohi,inputhihilo,inputhilolo,
             inputlohihi,inputlolohi,inputlohilo,inputlololo,
             outputhihihi_d,outputhilohi_d,outputhihilo_d,outputhilolo_d,
             outputlohihi_d,outputlolohi_d,outputlohilo_d,outputlololo_d,
             cnvjobs,addjobs,&cnvlapms,&addlapms,&timelapms_d,
             &walltimes_d,vrb);
      }
      double sumerr = 0.0;
      if(mode == 2)
         sumerr = dbl8_error_sum(dim,deg,
                     output1hihihi_h,output1hilohi_h,
                     output1hihilo_h,output1hilolo_h,
                     output1lohihi_h,output1lolohi_h,
                     output1lohilo_h,output1lololo_h,
                     output2hihihi_h,output2hilohi_h,
                     output2hihilo_h,output2hilolo_h,
                     output2lohihi_h,output2lolohi_h,
                     output2lohilo_h,output2lololo_h,
                     outputhihihi_d,outputhilohi_d,
                     outputhihilo_d,outputhilolo_d,
                     outputlohihi_d,outputlolohi_d,
                     outputlohilo_d,outputlololo_d,vrb);
 
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

int test_dbl8_sequence
 ( int seed, int dim, int nva, int nbr, int pwr, int vrblvl,
   bool jobrep, int mode )
{
   const double tol = 1.0e-120;

   int deg = 0;
   cout << "---> running in octo double precision for degree 0 ..." << endl;
   int fail = main_dbl8_test_polynomial
                 (seed,dim,nbr,nva,pwr,deg,vrblvl,tol,jobrep,mode);
   deg = 8;
   cout << "---> running for degree 8 ..." << endl;
   fail += main_dbl8_test_polynomial
              (seed,dim,nbr,nva,pwr,deg,vrblvl,tol,false,mode);
   deg = 15;
   cout << "---> running for degree 15 ..." << endl;
   fail += main_dbl8_test_polynomial
              (seed,dim,nbr,nva,pwr,deg,vrblvl,tol,false,mode);
   deg = 31;
   cout << "---> running for degree 31 ..." << endl;
   fail += main_dbl8_test_polynomial
              (seed,dim,nbr,nva,pwr,deg,vrblvl,tol,false,mode);
   deg = 63;
   cout << "---> running for degree 63 ..." << endl;
   fail += main_dbl8_test_polynomial
              (seed,dim,nbr,nva,pwr,deg,vrblvl,tol,false,mode);
   deg = 95;
   cout << "---> running for degree 95 ..." << endl;
   fail += main_dbl8_test_polynomial
              (seed,dim,nbr,nva,pwr,deg,vrblvl,tol,false,mode);
   deg = 127;
   cout << "---> running for degree 127 ..." << endl;
   fail += main_dbl8_test_polynomial
              (seed,dim,nbr,nva,pwr,deg,vrblvl,tol,false,mode);
   deg = 152;
   cout << "---> running for degree 152 ..." << endl;
   fail += main_dbl8_test_polynomial
              (seed,dim,nbr,nva,pwr,deg,vrblvl,tol,false,mode);
   deg = 159;
   cout << "---> running for degree 159 ..." << endl;
   fail += main_dbl8_test_polynomial
              (seed,dim,nbr,nva,pwr,deg,vrblvl,tol,false,mode);
   deg = 191;
   cout << "---> running for degree 191 ..." << endl;
   fail += main_dbl8_test_polynomial
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
