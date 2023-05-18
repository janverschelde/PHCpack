/* The file dbl3_polynomials_testers.cpp contains the definitions of
 * functions with prototypes in dbl3_polynomials_testers.h. */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <vector_types.h>
#include "random_polynomials.h"
#include "random3_monomials.h"
#include "random3_polynomials.h"
#include "convolution_jobs.h"
#include "addition_jobs.h"
#include "write_job_counts.h"
#include "write_gpu_timings.h"
#include "dbl3_polynomials_host.h"
#include "dbl3_polynomials_kernels.h"
#include "job_makers.h"
#include "dbl3_polynomials_testers.h"

using namespace std;

int main_dbl3_test_polynomial
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

   double realsum = test_dbl3_real_polynomial
                       (dim,nbr,nva,pwr,deg,vrblvl-1,jobrep,mode);
   double compsum = test_cmplx3_real_polynomial
                       (dim,nbr,nva,pwr,deg,vrblvl-1,jobrep,mode);

   int fail = int(realsum > tol) + int(compsum > tol);

   if(vrblvl > 0)
   {
      if(mode == 2)
      {
         cout << scientific << setprecision(2);
         cout << "Sum of all errors in triple double precision :" << endl;
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

void dbl3_make_input
 ( int dim, int nbr, int nva, int pwr, int deg,
   int *nvr, int **idx, int **exp,
   double **inputhi, double **inputmi, double **inputlo,
   double *csthi, double *cstmi, double *cstlo,
   double **cffhi, double **cffmi, double **cfflo, bool verbose )
{
   make_real3_input(dim,deg,inputhi,inputmi,inputlo);

   if(verbose)
   {
      cout << scientific << setprecision(16);
      cout << "Random input series :" << endl;
      for(int i=0; i<dim; i++)
      {
         cout << "-> coefficients of series " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
            cout << inputhi[i][j]
                 << "  " << inputmi[i][j]
                 << "  " << inputlo[i][j] << endl;
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
         make_real3_cyclic
            (dim,nva,deg,idx,csthi,cstmi,cstlo,cffhi,cffmi,cfflo);
      else
         make_real3_products
            (dim,nbr,nva,deg,idx,csthi,cstmi,cstlo,cffhi,cffmi,cfflo);
   }
   else
   {
      for(int i=0; i<nbr; i++) exp[i] = new int[nvr[i]];

      bool fail = make_real3_polynomial
                     (dim,nbr,pwr,deg,nvr,idx,exp,csthi,cstmi,cstlo,
                      cffhi,cffmi,cfflo);
   }
   if(verbose)
   {
      cout << "Coefficient series of the constant term :" << endl;
      for(int j=0; j<=deg; j++)
         cout << csthi[j] << "  " << cstmi[j]
                          << "  " << cstlo[j] << endl;

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
            cout << cffhi[i][j] << "  " << cffmi[i][j]
                                << "  " << cfflo[i][j] << endl;
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

void cmplx3_make_input
 ( int dim, int nbr, int nva, int pwr, int deg,
   int *nvr, int **idx, int **exp,
   double **inputrehi, double **inputremi, double **inputrelo,
   double **inputimhi, double **inputimmi, double **inputimlo,
   double *cstrehi, double *cstremi, double *cstrelo,
   double *cstimhi, double *cstimmi, double *cstimlo,
   double **cffrehi, double **cffremi, double **cffrelo,
   double **cffimhi, double **cffimmi, double **cffimlo, bool verbose )
{
   make_complex3_input(dim,deg,
      inputrehi,inputremi,inputrelo,inputimhi,inputimmi,inputimlo);

   if(verbose)
   {
      cout << scientific << setprecision(16);
      cout << "Random input series :" << endl;
      for(int i=0; i<dim; i++)
      {
         cout << "-> coefficients of series " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
         {
            cout << inputrehi[i][j]
                 << "  " << inputremi[i][j]
                 << "  " << inputrelo[i][j] << endl;
            cout << inputimhi[i][j]
                 << "  " << inputimmi[i][j]
                 << "  " << inputimlo[i][j] << endl;
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
         make_complex3_cyclic
            (dim,nva,deg,idx,
             cstrehi,cstremi,cstrelo,cstimhi,cstimmi,cstimlo,
             cffrehi,cffremi,cffrelo,cffimhi,cffimmi,cffimlo);
      else
         make_complex3_products
            (dim,nbr,nva,deg,idx,
             cstrehi,cstremi,cstrelo,cstimhi,cstimmi,cstimlo,
             cffrehi,cffremi,cffrelo,cffimhi,cffimmi,cffimlo);
   }
   else
   {
      for(int i=0; i<nbr; i++) exp[i] = new int[nvr[i]];

      bool fail = make_complex3_polynomial
                     (dim,nbr,pwr,deg,nvr,idx,exp,
                      cstrehi,cstremi,cstrelo,cstimhi,cstimmi,cstimlo,
                      cffrehi,cffremi,cffrelo,cffimhi,cffimmi,cffimlo);
   }
   if(verbose)
   {
      cout << "Coefficient series of the constant term :" << endl;
      for(int j=0; j<=deg; j++)
      {
         cout << cstrehi[j] << "  " << cstremi[j]
                            << "  " << cstrelo[j] << endl;
         cout << cstimhi[j] << "  " << cstimmi[j]
                            << "  " << cstimlo[j] << endl;
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
            cout << cffrehi[i][j] << "  " << cffremi[i][j]
                                  << "  " << cffrelo[i][j] << endl;
            cout << cffimhi[i][j] << "  " << cffimmi[i][j]
                                  << "  " << cffimlo[i][j] << endl;
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

double dbl3_error_sum
 ( int dim, int deg,
   double **results1hi_h, double **results1mi_h, double **results1lo_h,
   double **results2hi_h, double **results2mi_h, double **results2lo_h,
   double **resultshi_d, double **resultsmi_d, double **resultslo_d,
   bool verbose )
{
   double err = 0.0;

   if(verbose) cout << "The value of the polynomial :" << endl;
   for(int i=0; i<=deg; i++)
   {
      if(verbose)
      {
         cout << results1hi_h[dim][i] << "  "
              << results1mi_h[dim][i] << "  "
              << results1lo_h[dim][i] << endl;
         cout << results2hi_h[dim][i] << "  "
              << results2mi_h[dim][i] << "  "
              << results2lo_h[dim][i] << endl;
         cout << resultshi_d[dim][i] << "  "
              << resultsmi_d[dim][i] << "  "
              << resultslo_d[dim][i] << endl;
      }
      err = err + abs(results1hi_h[dim][i] - results2hi_h[dim][i])
                + abs(results1mi_h[dim][i] - results2mi_h[dim][i])
                + abs(results1lo_h[dim][i] - results2lo_h[dim][i])
                + abs(results1hi_h[dim][i] - resultshi_d[dim][i])
                + abs(results1mi_h[dim][i] - resultsmi_d[dim][i])
                + abs(results1lo_h[dim][i] - resultslo_d[dim][i]);
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
            cout << results1hi_h[k][i] << "  "
                 << results1mi_h[k][i] << "  "
                 << results1lo_h[k][i] << endl;
            cout << results2hi_h[k][i] << "  "
                 << results2mi_h[k][i] << "  "
                 << results2lo_h[k][i] << endl;
            cout << resultshi_d[k][i] << "  "
                 << resultsmi_d[k][i] << "  "
                 << resultslo_d[k][i] << endl;
         }
         err = err + abs(results1hi_h[k][i] - results2hi_h[k][i])
                   + abs(results1mi_h[k][i] - results2mi_h[k][i])
                   + abs(results1lo_h[k][i] - results2lo_h[k][i])
                   + abs(results1hi_h[k][i] - resultshi_d[k][i])
                   + abs(results1mi_h[k][i] - resultsmi_d[k][i])
                   + abs(results1lo_h[k][i] - resultslo_d[k][i]);
      }
      if(verbose) cout << "error : " << err << endl;
      sumerr = sumerr + err;
   }
   return sumerr;
}

double cmplx3_error_sum
 ( int dim, int deg,
   double **results1rehi_h, double **results1remi_h, double **results1relo_h,
   double **results1imhi_h, double **results1immi_h, double **results1imlo_h,
   double **results2rehi_h, double **results2remi_h, double **results2relo_h,
   double **results2imhi_h, double **results2immi_h, double **results2imlo_h,
   double **resultsrehi_d, double **resultsremi_d, double **resultsrelo_d,
   double **resultsimhi_d, double **resultsimmi_d, double **resultsimlo_d,
   bool verbose )
{
   double err = 0.0;

   if(verbose) cout << "The value of the polynomial :" << endl;
   for(int i=0; i<=deg; i++)
   {
      if(verbose)
      {
         cout << results1rehi_h[dim][i] << "  "
              << results1remi_h[dim][i] << "  "
              << results1relo_h[dim][i] << endl;
         cout << results1imhi_h[dim][i] << "  "
              << results1immi_h[dim][i] << "  "
              << results1imlo_h[dim][i] << endl;
         cout << results2rehi_h[dim][i] << "  "
              << results2remi_h[dim][i] << "  "
              << results2relo_h[dim][i] << endl;
         cout << results2imhi_h[dim][i] << "  "
              << results2immi_h[dim][i] << "  "
              << results2imlo_h[dim][i] << endl;
         cout << resultsrehi_d[dim][i] << "  "
              << resultsremi_d[dim][i] << "  "
              << resultsrelo_d[dim][i] << endl;
         cout << resultsimhi_d[dim][i] << "  "
              << resultsimmi_d[dim][i] << "  "
              << resultsimlo_d[dim][i] << endl;
      }
      err = err + abs(results1rehi_h[dim][i] - results2rehi_h[dim][i])
                + abs(results1remi_h[dim][i] - results2remi_h[dim][i])
                + abs(results1relo_h[dim][i] - results2relo_h[dim][i])
                + abs(results1imhi_h[dim][i] - results2imhi_h[dim][i])
                + abs(results1immi_h[dim][i] - results2immi_h[dim][i])
                + abs(results1imlo_h[dim][i] - results2imlo_h[dim][i])
                + abs(results1rehi_h[dim][i] - resultsrehi_d[dim][i])
                + abs(results1remi_h[dim][i] - resultsremi_d[dim][i])
                + abs(results1relo_h[dim][i] - resultsrelo_d[dim][i])
                + abs(results1imhi_h[dim][i] - resultsimhi_d[dim][i])
                + abs(results1immi_h[dim][i] - resultsimmi_d[dim][i])
                + abs(results1imlo_h[dim][i] - resultsimlo_d[dim][i]);
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
            cout << results1rehi_h[k][i] << "  "
                 << results1remi_h[k][i] << "  "
                 << results1relo_h[k][i] << endl;
            cout << results1imhi_h[k][i] << "  "
                 << results1immi_h[k][i] << "  "
                 << results1imlo_h[k][i] << endl;
            cout << results2rehi_h[k][i] << "  "
                 << results2remi_h[k][i] << "  "
                 << results2relo_h[k][i] << endl;
            cout << results2imhi_h[k][i] << "  "
                 << results2immi_h[k][i] << "  "
                 << results2imlo_h[k][i] << endl;
            cout << resultsrehi_d[k][i] << "  "
                 << resultsremi_d[k][i] << "  "
                 << resultsrelo_d[k][i] << endl;
            cout << resultsimhi_d[k][i] << "  "
                 << resultsimmi_d[k][i] << "  "
                 << resultsimlo_d[k][i] << endl;
         }
         err = err + abs(results1rehi_h[k][i] - results2rehi_h[k][i])
                   + abs(results1remi_h[k][i] - results2remi_h[k][i])
                   + abs(results1relo_h[k][i] - results2relo_h[k][i])
                   + abs(results1imhi_h[k][i] - results2imhi_h[k][i])
                   + abs(results1immi_h[k][i] - results2immi_h[k][i])
                   + abs(results1imlo_h[k][i] - results2imlo_h[k][i])
                   + abs(results1rehi_h[k][i] - resultsrehi_d[k][i])
                   + abs(results1remi_h[k][i] - resultsremi_d[k][i])
                   + abs(results1relo_h[k][i] - resultsrelo_d[k][i])
                   + abs(results1imhi_h[k][i] - resultsimhi_d[k][i])
                   + abs(results1immi_h[k][i] - resultsimmi_d[k][i])
                   + abs(results1imlo_h[k][i] - resultsimlo_d[k][i]);
      }
      if(verbose) cout << "error : " << err << endl;
      sumerr = sumerr + err;
   }
   return sumerr;
}

double test_dbl3_real_polynomial
 ( int dim, int nbr, int nva, int pwr, int deg, int verbose, bool jobrep,
   int mode )
{
   if(nbr < 1)
      return 0.0;
   else
   {
      double **inputhi = new double*[dim]; // dim series of degree deg
      for(int i=0; i<dim; i++) inputhi[i] = new double[deg+1];
      double **inputmi = new double*[dim];
      for(int i=0; i<dim; i++) inputmi[i] = new double[deg+1];
      double **inputlo = new double*[dim];
      for(int i=0; i<dim; i++) inputlo[i] = new double[deg+1];
      // The output are dim+1 power series of degree deg
      // for the evaluated and differentiated polynomial.
      double **output1hi_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output1hi_h[i] = new double[deg+1];
      double **output1mi_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output1mi_h[i] = new double[deg+1];
      double **output1lo_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output1lo_h[i] = new double[deg+1];
      double **output2hi_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output2hi_h[i] = new double[deg+1];
      double **output2mi_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output2mi_h[i] = new double[deg+1];
      double **output2lo_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output2lo_h[i] = new double[deg+1];
      double **outputhi_d = new double*[dim+1];
      for(int i=0; i<=dim; i++) outputhi_d[i] = new double[deg+1];
      double **outputmi_d = new double*[dim+1];
      for(int i=0; i<=dim; i++) outputmi_d[i] = new double[deg+1];
      double **outputlo_d = new double*[dim+1];
      for(int i=0; i<=dim; i++) outputlo_d[i] = new double[deg+1];

      double *csthi = new double[deg+1]; // constant coefficient series
      double *cstmi = new double[deg+1];
      double *cstlo = new double[deg+1];
      double **cffhi = new double*[nbr]; // coefficient series of terms
      for(int i=0; i<nbr; i++) cffhi[i] = new double[deg+1];
      double **cffmi = new double*[nbr];
      for(int i=0; i<nbr; i++) cffmi[i] = new double[deg+1];
      double **cfflo = new double*[nbr];
      for(int i=0; i<nbr; i++) cfflo[i] = new double[deg+1];
      int *nvr = new int[nbr]; // number of variables in each monomial
      int **idx = new int*[nbr];  // indices of variables in monomials
      int **exp = new int*[nbr];  // exponents of the variables

      bool vrb = (verbose > 1);

      dbl3_make_input(dim,nbr,nva,pwr,deg,nvr,idx,exp,
                      inputhi,inputmi,inputlo,csthi,cstmi,cstlo,
                      cffhi,cffmi,cfflo,vrb);

      ConvolutionJobs cnvjobs(dim);
      AdditionJobs addjobs(dim,nbr);

      make_all_jobs(dim,nbr,nvr,idx,&cnvjobs,&addjobs,vrb);

      double timelapsec1_h,timelapsec2_h;
      double cnvlapms,addlapms,timelapms_d,walltimes_d;

      if((mode == 1) || (mode == 2))
      {
         if(vrb) cout << "Computing without convolution jobs ..." << endl;
         CPU_dbl3_poly_evaldiff
            (dim,nbr,deg,nvr,idx,csthi,cstmi,cstlo,cffhi,cffmi,cfflo,
             inputhi,inputmi,inputlo,output1hi_h,output1mi_h,output1lo_h,
             &timelapsec1_h,vrb);
         if(vrb) cout << "Computing with convolution jobs ..." << endl;
         CPU_dbl3_poly_evaldiffjobs
            (dim,nbr,deg,nvr,idx,csthi,cstmi,cstlo,cffhi,cffmi,cfflo,
             inputhi,inputmi,inputlo,output2hi_h,output2mi_h,output2lo_h,
             cnvjobs,addjobs,&timelapsec2_h,vrb);
      }
      if((mode == 0) || (mode == 2))
      {
         if(vrb) cout << "Computing on the device ..." << endl;
         GPU_dbl3_poly_evaldiff
            (deg+1,dim,nbr,deg,nvr,idx,csthi,cstmi,cstlo,cffhi,cffmi,cfflo,
             inputhi,inputmi,inputlo,outputhi_d,outputmi_d,outputlo_d,
             cnvjobs,addjobs,&cnvlapms,&addlapms,&timelapms_d,
             &walltimes_d,vrb);
      }
      double sumerr = 0.0;
      if(mode == 2)
         sumerr = dbl3_error_sum(dim,deg,
                     output1hi_h,output1mi_h,output1lo_h,
                     output2hi_h,output2mi_h,output2lo_h,
                     outputhi_d,outputmi_d,outputlo_d,vrb);

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

double test_cmplx3_real_polynomial
 ( int dim, int nbr, int nva, int pwr, int deg, int verbose, bool jobrep,
   int mode )
{
   if(nbr < 1)
      return 0.0;
   else
   {
      double **inputrehi = new double*[dim]; // dim series of degree deg
      double **inputremi = new double*[dim];
      double **inputrelo = new double*[dim];
      double **inputimhi = new double*[dim];
      double **inputimmi = new double*[dim];
      double **inputimlo = new double*[dim];
      for(int i=0; i<dim; i++)
      {
         inputrehi[i] = new double[deg+1];
         inputremi[i] = new double[deg+1];
         inputrelo[i] = new double[deg+1];
         inputimhi[i] = new double[deg+1];
         inputimmi[i] = new double[deg+1];
         inputimlo[i] = new double[deg+1];
      }
      // The output are dim+1 power series of degree deg
      // for the evaluated and differentiated polynomial.
      double **output1rehi_h = new double*[dim+1];
      double **output1remi_h = new double*[dim+1];
      double **output1relo_h = new double*[dim+1];
      double **output1imhi_h = new double*[dim+1];
      double **output1immi_h = new double*[dim+1];
      double **output1imlo_h = new double*[dim+1];
      for(int i=0; i<=dim; i++)
      {
         output1rehi_h[i] = new double[deg+1];
         output1remi_h[i] = new double[deg+1];
         output1relo_h[i] = new double[deg+1];
         output1imhi_h[i] = new double[deg+1];
         output1immi_h[i] = new double[deg+1];
         output1imlo_h[i] = new double[deg+1];
      }
      double **output2rehi_h = new double*[dim+1];
      double **output2remi_h = new double*[dim+1];
      double **output2relo_h = new double*[dim+1];
      double **output2imhi_h = new double*[dim+1];
      double **output2immi_h = new double*[dim+1];
      double **output2imlo_h = new double*[dim+1];
      for(int i=0; i<=dim; i++)
      {
         output2rehi_h[i] = new double[deg+1];
         output2remi_h[i] = new double[deg+1];
         output2relo_h[i] = new double[deg+1];
         output2imhi_h[i] = new double[deg+1];
         output2immi_h[i] = new double[deg+1];
         output2imlo_h[i] = new double[deg+1];
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
      double *cstrehi = new double[deg+1]; // constant coefficient series
      double *cstremi = new double[deg+1];
      double *cstrelo = new double[deg+1];
      double *cstimhi = new double[deg+1];
      double *cstimmi = new double[deg+1];
      double *cstimlo = new double[deg+1];
      double **cffrehi = new double*[nbr]; // coefficient series of terms
      double **cffremi = new double*[nbr];
      double **cffrelo = new double*[nbr];
      double **cffimhi = new double*[nbr];
      double **cffimmi = new double*[nbr];
      double **cffimlo = new double*[nbr];
      for(int i=0; i<nbr; i++)
      {
         cffrehi[i] = new double[deg+1];
         cffremi[i] = new double[deg+1];
         cffrelo[i] = new double[deg+1];
         cffimhi[i] = new double[deg+1];
         cffimmi[i] = new double[deg+1];
         cffimlo[i] = new double[deg+1];
      }
      int *nvr = new int[nbr]; // number of variables in each monomial
      int **idx = new int*[nbr];  // indices of variables in monomials
      int **exp = new int*[nbr];  // exponents of the variables

      bool vrb = (verbose > 1);

      cmplx3_make_input(dim,nbr,nva,pwr,deg,nvr,idx,exp,
         inputrehi,inputremi,inputrelo,inputimhi,inputimmi,inputimlo,
         cstrehi,cstremi,cstrelo,cstimhi,cstimmi,cstimlo,
         cffrehi,cffremi,cffrelo,cffimhi,cffimmi,cffimlo,vrb);

      ConvolutionJobs cnvjobs(dim);
      AdditionJobs addjobs(dim,nbr);

      make_all_jobs(dim,nbr,nvr,idx,&cnvjobs,&addjobs,vrb);

      double timelapsec1_h,timelapsec2_h;
      double cnvlapms,addlapms,timelapms_d,walltimes_d;

      if((mode == 1) || (mode == 2))
      {
         if(vrb) cout << "Computing without convolution jobs ..." << endl;
         CPU_cmplx3_poly_evaldiff
            (dim,nbr,deg,nvr,idx,
             cstrehi,cstremi,cstrelo,cstimhi,cstimmi,cstimlo,
             cffrehi,cffremi,cffrelo,cffimhi,cffimmi,cffimlo,
             inputrehi,inputremi,inputrelo,inputimhi,inputimmi,inputimlo,
             output1rehi_h,output1remi_h,output1relo_h,
             output1imhi_h,output1immi_h,output1imlo_h,&timelapsec1_h,vrb);
         if(vrb) cout << "Computing with convolution jobs ..." << endl;
         CPU_cmplx3_poly_evaldiffjobs
            (dim,nbr,deg,nvr,idx,
             cstrehi,cstremi,cstrelo,cstimhi,cstimmi,cstimlo,
             cffrehi,cffremi,cffrelo,cffimhi,cffimmi,cffimlo,
             inputrehi,inputremi,inputrelo,inputimhi,inputimmi,inputimlo,
             output2rehi_h,output2remi_h,output2relo_h,
             output2imhi_h,output2immi_h,output2imlo_h,
             cnvjobs,addjobs,&timelapsec2_h,vrb);
      }
      if((mode == 0) || (mode == 2))
      {
         if(vrb) cout << "Computing on the device ..." << endl;
         GPU_cmplx3_poly_evaldiff
            (deg+1,dim,nbr,deg,nvr,idx,
             cstrehi,cstremi,cstrelo,cstimhi,cstimmi,cstimlo,
             cffrehi,cffremi,cffrelo,cffimhi,cffimmi,cffimlo,
             inputrehi,inputremi,inputrelo,inputimhi,inputimmi,inputimlo,
             outputrehi_d,outputremi_d,outputrelo_d,
             outputimhi_d,outputimmi_d,outputimlo_d,
             cnvjobs,addjobs,&cnvlapms,&addlapms,&timelapms_d,
             &walltimes_d,vrb);
      }
      double sumerr = 0.0;
      if(mode == 2)
         sumerr = cmplx3_error_sum(dim,deg,
                     output1rehi_h,output1remi_h,output1relo_h,
                     output1imhi_h,output1immi_h,output1imlo_h,
                     output2rehi_h,output2remi_h,output2relo_h,
                     output2imhi_h,output2immi_h,output2imlo_h,
                     outputrehi_d,outputremi_d,outputrelo_d,
                     outputimhi_d,outputimmi_d,outputimlo_d,vrb);

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

int test_dbl3_sequence
 ( int seed, int dim, int nva, int nbr, int pwr, int vrblvl,
   bool jobrep, int mode )
{
   const double tol = 1.0e-40;

   int deg = 0;
   cout << "---> running in triple double precision for degree 0 ..."
        << endl;
   int fail = main_dbl3_test_polynomial
                 (seed,dim,nbr,nva,pwr,deg,vrblvl,tol,jobrep,mode);
   deg = 8;
   cout << "---> running for degree 8 ..." << endl;
   fail += main_dbl3_test_polynomial
              (seed,dim,nbr,nva,pwr,deg,vrblvl,tol,false,mode);
   deg = 15;
   cout << "---> running for degree 15 ..." << endl;
   fail += main_dbl3_test_polynomial
              (seed,dim,nbr,nva,pwr,deg,vrblvl,tol,false,mode);
   deg = 31;
   cout << "---> running for degree 31 ..." << endl;
   fail += main_dbl3_test_polynomial
              (seed,dim,nbr,nva,pwr,deg,vrblvl,tol,false,mode);
   deg = 63;
   cout << "---> running for degree 63 ..." << endl;
   fail += main_dbl3_test_polynomial
              (seed,dim,nbr,nva,pwr,deg,vrblvl,tol,false,mode);
   deg = 95;
   cout << "---> running for degree 95 ..." << endl;
   fail += main_dbl3_test_polynomial
              (seed,dim,nbr,nva,pwr,deg,vrblvl,tol,false,mode);
   deg = 127;
   cout << "---> running for degree 127 ..." << endl;
   fail += main_dbl3_test_polynomial
              (seed,dim,nbr,nva,pwr,deg,vrblvl,tol,false,mode);
   deg = 152;
   cout << "---> running for degree 152 ..." << endl;
   fail += main_dbl3_test_polynomial
              (seed,dim,nbr,nva,pwr,deg,vrblvl,tol,false,mode);
   deg = 159;
   cout << "---> running for degree 159 ..." << endl;
   fail += main_dbl3_test_polynomial
              (seed,dim,nbr,nva,pwr,deg,vrblvl,tol,false,mode);
   deg = 191;
   cout << "---> running for degree 191 ..." << endl;
   fail += main_dbl3_test_polynomial
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
