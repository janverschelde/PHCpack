/* The file dbl2_polynomials_testers.cpp contains the definitions of
 * functions with prototypes in dbl2_polynomials_testers.h. */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <vector_types.h>
#include "random_polynomials.h"
#include "random2_monomials.h"
#include "random2_polynomials.h"
#include "convolution_jobs.h"
#include "addition_jobs.h"
#include "write_job_counts.h"
#include "write_gpu_timings.h"
#include "dbl2_polynomials_host.h"
#include "dbl2_polynomials_kernels.h"
#include "job_makers.h"
#include "dbl2_polynomials_testers.h"

using namespace std;

int main_dbl2_test_polynomial
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

   double realsum = test_dbl2_polynomial
                       (dim,nbr,nva,pwr,deg,vrblvl-1,jobrep,mode);
   double compsum = test_cmplx2_polynomial
                       (dim,nbr,nva,pwr,deg,vrblvl-1,jobrep,mode);

   int fail = int(realsum > tol) + int(compsum > tol);

   if(vrblvl > 0)
   {
      if(mode == 2)
      {
         cout << scientific << setprecision(2);
         cout << "Sum of all errors in double double precision :" << endl;
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

int dbl2_make_input
 ( int dim, int nbr, int nva, int pwr, int deg,
   int *nvr, int **idx, int **exp,
   double **inputhi, double **inputlo,
   double *csthi, double *cstlo, double **cffhi, double **cfflo,
   bool verbose )
{
   make_real2_input(dim,deg,inputhi,inputlo);

   if(verbose)
   {
      cout << scientific << setprecision(16);
      cout << "Random input series :" << endl;
      for(int i=0; i<dim; i++)
      {
         cout << "-> coefficients of series " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
            cout << inputhi[i][j] << "  " << inputlo[i][j] << endl;
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
         make_real2_cyclic(dim,nva,deg,idx,csthi,cstlo,cffhi,cfflo);
      else
         make_real2_products(dim,nbr,nva,deg,idx,csthi,cstlo,cffhi,cfflo);
   }
   else
   {
      for(int i=0; i<nbr; i++) exp[i] = new int[nvr[i]];

      bool fail = make_real2_polynomial(dim,nbr,pwr,deg,nvr,idx,exp,
                                        csthi,cstlo,cffhi,cfflo);
   }
   if(verbose)
   {
      cout << "Coefficient series of the constant term :" << endl;
      for(int j=0; j<=deg; j++)
         cout << csthi[j] << "  " << cstlo[j] << endl;

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
            cout << cffhi[i][j] << "  " << cfflo[i][j] << endl;
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

int cmplx2_make_input
 ( int dim, int nbr, int nva, int pwr, int deg,
   int *nvr, int **idx, int **exp,
   double **inputrehi, double **inputrelo,
   double **inputimhi, double **inputimlo,
   double *cstrehi, double *cstrelo, double *cstimhi, double *cstimlo,
   double **cffrehi, double **cffrelo, double **cffimhi, double **cffimlo,
   bool verbose )
{
   make_complex2_input
      (dim,deg,inputrehi,inputrelo,inputimhi,inputimlo);

   if(verbose)
   {
      cout << scientific << setprecision(16);
      cout << "Random input series :" << endl;
      for(int i=0; i<dim; i++)
      {
         cout << "-> coefficients of series " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
            cout << inputrehi[i][j] << "  " << inputrelo[i][j] << endl
                 << inputimhi[i][j] << "  " << inputimlo[i][j] << endl;
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
         make_complex2_cyclic
            (dim,nva,deg,idx,cstrehi,cstrelo,cstimhi,cstimlo,
             cffrehi,cffrelo,cffimhi,cffimlo);
      else
         make_complex2_products
            (dim,nbr,nva,deg,idx,cstrehi,cstrelo,cstimhi,cstimlo,
             cffrehi,cffrelo,cffimhi,cffimlo);
   }
   else
   {
      for(int i=0; i<nbr; i++) exp[i] = new int[nvr[i]];

      bool fail = make_complex2_polynomial
                     (dim,nbr,pwr,deg,nvr,idx,exp,
                      cstrehi,cstrelo,cstimhi,cstimlo,
                      cffrehi,cffrelo,cffimhi,cffimlo);
   }
   if(verbose)
   {
      cout << "Coefficient series of the constant term :" << endl;
      for(int j=0; j<=deg; j++)
         cout << cstrehi[j] << "  " << cstrelo[j] << endl
              << cstimhi[j] << "  " << cstimlo[j] << endl;

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
            cout << cffrehi[i][j] << "  " << cffrelo[i][j] << endl
                 << cffimhi[i][j] << "  " << cffimlo[i][j] << endl;
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

double dbl2_error_sum1
 ( int dim, int deg,
   double **resultshi_h, double **resultslo_h,
   double **resultshi_d, double **resultslo_d, bool verbose )

{
   double err = 0.0;

   if(verbose) cout << "The value of the polynomial :" << endl;
   for(int i=0; i<=deg; i++)
   {
      if(verbose)
      {
         cout << resultshi_h[dim][i] << "  "
              << resultslo_h[dim][i] << endl;
         cout << resultshi_d[dim][i] << "  "
              << resultslo_d[dim][i] << endl;
      }
      err = err + abs(resultshi_h[dim][i] - resultshi_d[dim][i])
                + abs(resultslo_h[dim][i] - resultslo_d[dim][i]);
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
            cout << resultshi_h[k][i] << "  "
                 << resultslo_h[k][i] << endl;
            cout << resultshi_d[k][i] << "  "
                 << resultslo_d[k][i] << endl;
         }
         err = err + abs(resultshi_h[k][i] - resultshi_d[k][i])
                   + abs(resultslo_h[k][i] - resultslo_d[k][i]);
      }
      if(verbose) cout << "error : " << err << endl;
      sumerr = sumerr + err;
   }
   return sumerr;
}

double dbl2_error_sum
 ( int dim, int deg,
   double **results1hi_h, double **results1lo_h,
   double **results2hi_h, double **results2lo_h,
   double **resultshi_d, double **resultslo_d, bool verbose )
{
   double err = 0.0;

   if(verbose) cout << "The value of the polynomial :" << endl;
   for(int i=0; i<=deg; i++)
   {
      if(verbose)
      {
         cout << results1hi_h[dim][i] << "  "
              << results1lo_h[dim][i] << endl;
         cout << results2hi_h[dim][i] << "  "
              << results2lo_h[dim][i] << endl;
         cout << resultshi_d[dim][i] << "  "
              << resultslo_d[dim][i] << endl;
      }
      err = err + abs(results1hi_h[dim][i] - results2hi_h[dim][i])
                + abs(results1lo_h[dim][i] - results2lo_h[dim][i])
                + abs(results1hi_h[dim][i] - resultshi_d[dim][i])
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
                 << results1lo_h[k][i] << endl;
            cout << results2hi_h[k][i] << "  "
                 << results2lo_h[k][i] << endl;
            cout << resultshi_d[k][i] << "  "
                 << resultslo_d[k][i] << endl;
         }
         err = err + abs(results1hi_h[k][i] - results2hi_h[k][i])
                   + abs(results1lo_h[k][i] - results2lo_h[k][i])
                   + abs(results1hi_h[k][i] - resultshi_d[k][i])
                   + abs(results1lo_h[k][i] - resultslo_d[k][i]);
      }
      if(verbose) cout << "error : " << err << endl;
      sumerr = sumerr + err;
   }
   return sumerr;
}

double cmplx2_error_sum1
 ( int dim, int deg,
   double **resultsrehi_h, double **resultsrelo_h,
   double **resultsimhi_h, double **resultsimlo_h,
   double **resultsrehi_d, double **resultsrelo_d,
   double **resultsimhi_d, double **resultsimlo_d, bool verbose )
{
   double err = 0.0;

   if(verbose) cout << "The value of the polynomial :" << endl;
   for(int i=0; i<=deg; i++)
   {
      if(verbose)
      {
         cout << resultsrehi_h[dim][i] << "  "
              << resultsrelo_h[dim][i] << endl;
         cout << resultsimhi_h[dim][i] << "  "
              << resultsimlo_h[dim][i] << endl;
         cout << resultsrehi_d[dim][i] << "  "
              << resultsrelo_d[dim][i] << endl;
         cout << resultsimhi_d[dim][i] << "  "
              << resultsimlo_d[dim][i] << endl;
      }
      err = err + abs(resultsrehi_h[dim][i] - resultsrehi_d[dim][i])
                + abs(resultsrelo_h[dim][i] - resultsrelo_d[dim][i])
                + abs(resultsimhi_h[dim][i] - resultsimhi_d[dim][i])
                + abs(resultsimlo_h[dim][i] - resultsimlo_d[dim][i]);
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
            cout << resultsrehi_h[k][i] << "  "
                 << resultsrelo_h[k][i] << endl;
            cout << resultsimhi_h[k][i] << "  "
                 << resultsimlo_h[k][i] << endl;
            cout << resultsrehi_d[k][i] << "  "
                 << resultsrelo_d[k][i] << endl;
            cout << resultsimhi_d[k][i] << "  "
                 << resultsimlo_d[k][i] << endl;
         }
         err = err + abs(resultsrehi_h[k][i] - resultsrehi_d[k][i])
                   + abs(resultsrelo_h[k][i] - resultsrelo_d[k][i])
                   + abs(resultsimhi_h[k][i] - resultsimhi_d[k][i])
                   + abs(resultsimlo_h[k][i] - resultsimlo_d[k][i]);
      }
      if(verbose) cout << "error : " << err << endl;
      sumerr = sumerr + err;
   }
   return sumerr;
}

double cmplx2_error_sum
 ( int dim, int deg,
   double **results1rehi_h, double **results1relo_h,
   double **results1imhi_h, double **results1imlo_h,
   double **results2rehi_h, double **results2relo_h,
   double **results2imhi_h, double **results2imlo_h,
   double **resultsrehi_d, double **resultsrelo_d,
   double **resultsimhi_d, double **resultsimlo_d, bool verbose )
{
   double err = 0.0;

   if(verbose) cout << "The value of the polynomial :" << endl;
   for(int i=0; i<=deg; i++)
   {
      if(verbose)
      {
         cout << results1rehi_h[dim][i] << "  "
              << results1relo_h[dim][i] << endl;
         cout << results1imhi_h[dim][i] << "  "
              << results1imlo_h[dim][i] << endl;
         cout << results2rehi_h[dim][i] << "  "
              << results2relo_h[dim][i] << endl;
         cout << results2imhi_h[dim][i] << "  "
              << results2imlo_h[dim][i] << endl;
         cout << resultsrehi_d[dim][i] << "  "
              << resultsrelo_d[dim][i] << endl;
         cout << resultsimhi_d[dim][i] << "  "
              << resultsimlo_d[dim][i] << endl;
      }
      err = err + abs(results1rehi_h[dim][i] - results2rehi_h[dim][i])
                + abs(results1relo_h[dim][i] - results2relo_h[dim][i])
                + abs(results1imhi_h[dim][i] - results2imhi_h[dim][i])
                + abs(results1imlo_h[dim][i] - results2imlo_h[dim][i])
                + abs(results1rehi_h[dim][i] - resultsrehi_d[dim][i])
                + abs(results1relo_h[dim][i] - resultsrelo_d[dim][i])
                + abs(results1imhi_h[dim][i] - resultsimhi_d[dim][i])
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
                 << results1relo_h[k][i] << endl;
            cout << results1imhi_h[k][i] << "  "
                 << results1imlo_h[k][i] << endl;
            cout << results2rehi_h[k][i] << "  "
                 << results2relo_h[k][i] << endl;
            cout << results2imhi_h[k][i] << "  "
                 << results2imlo_h[k][i] << endl;
            cout << resultsrehi_d[k][i] << "  "
                 << resultsrelo_d[k][i] << endl;
            cout << resultsimhi_d[k][i] << "  "
                 << resultsimlo_d[k][i] << endl;
         }
         err = err + abs(results1rehi_h[k][i] - results2rehi_h[k][i])
                   + abs(results1relo_h[k][i] - results2relo_h[k][i])
                   + abs(results1imhi_h[k][i] - results2imhi_h[k][i])
                   + abs(results1imlo_h[k][i] - results2imlo_h[k][i])
                   + abs(results1rehi_h[k][i] - resultsrehi_d[k][i])
                   + abs(results1relo_h[k][i] - resultsrelo_d[k][i])
                   + abs(results1imhi_h[k][i] - resultsimhi_d[k][i])
                   + abs(results1imlo_h[k][i] - resultsimlo_d[k][i]);
      }
      if(verbose) cout << "error : " << err << endl;
      sumerr = sumerr + err;
   }
   return sumerr;
}

double test_dbl2_polynomial
 ( int dim, int nbr, int nva, int pwr, int deg, int verbose, bool jobrep,
   int mode )
{
   cout << "*** testing the real arithmetic on real data ***" << endl;

   if(nbr < 1)
      return 0.0;
   else
   {
      double **inputhi = new double*[dim]; // dim series of degree deg
      for(int i=0; i<dim; i++) inputhi[i] = new double[deg+1];
      double **inputlo = new double*[dim];
      for(int i=0; i<dim; i++) inputlo[i] = new double[deg+1];
      // The output are dim+1 power series of degree deg
      // for the evaluated and differentiated polynomial.
      double **output1hi_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output1hi_h[i] = new double[deg+1];
      double **output1lo_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output1lo_h[i] = new double[deg+1];
      double **output2hi_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output2hi_h[i] = new double[deg+1];
      double **output2lo_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output2lo_h[i] = new double[deg+1];
      double **outputhi_d = new double*[dim+1];
      for(int i=0; i<=dim; i++) outputhi_d[i] = new double[deg+1];
      double **outputlo_d = new double*[dim+1];
      for(int i=0; i<=dim; i++) outputlo_d[i] = new double[deg+1];

      double *csthi = new double[deg+1]; // constant coefficient series
      double *cstlo = new double[deg+1];
      double **cffhi = new double*[nbr]; // coefficient series of terms
      for(int i=0; i<nbr; i++) cffhi[i] = new double[deg+1];
      double **cfflo = new double*[nbr]; // coefficient series of terms
      for(int i=0; i<nbr; i++) cfflo[i] = new double[deg+1];
      int *nvr = new int[nbr]; // number of variables in each monomial
      int **idx = new int*[nbr];  // indices of variables in monomials
      int **exp = new int*[nbr];  // exponents of the variables

      bool vrb = (verbose > 1);

      int fail = dbl2_make_input(dim,nbr,nva,pwr,deg,nvr,idx,exp,inputhi,
                                 inputlo,csthi,cstlo,cffhi,cfflo,vrb);
      if(fail == 1)
      {
          cout << "Duplicates in support, returning 0 as error." << endl;
          return 0.0;
      }
      ConvolutionJobs cnvjobs(dim);
      AdditionJobs addjobs(dim,nbr);

      make_all_jobs(dim,nbr,nvr,idx,&cnvjobs,&addjobs,vrb);

      double timelapsec1_h,timelapsec2_h;
      double cnvlapms,addlapms,timelapms_d,walltimes_d;

      if((mode == 1) || (mode == 2))
      {
         if(vrb) cout << "Computing without convolution jobs ..." << endl;
         CPU_dbl2_poly_evaldiff
            (dim,nbr,deg,nvr,idx,csthi,cstlo,cffhi,cfflo,inputhi,inputlo,
             output1hi_h,output1lo_h,&timelapsec1_h,vrb);
         if(vrb) cout << "Computing with convolution jobs ..." << endl;
         CPU_dbl2_poly_evaldiffjobs
            (dim,nbr,deg,nvr,idx,csthi,cstlo,cffhi,cfflo,inputhi,inputlo,
             output2hi_h,output2lo_h,cnvjobs,addjobs,&timelapsec2_h,vrb);
      }
      if((mode == 0) || (mode == 2))
      {
         if(vrb) cout << "Computing on the device ..." << endl;
         GPU_dbl2_poly_evaldiff
            (deg+1,dim,nbr,deg,nvr,idx,csthi,cstlo,cffhi,cfflo,
             inputhi,inputlo,outputhi_d,outputlo_d,cnvjobs,addjobs,
             &cnvlapms,&addlapms,&timelapms_d,&walltimes_d,vrb);
      }
      double sumerr = 0.0;
      if(mode == 2)
      {
         sumerr = dbl2_error_sum(dim,deg,output1hi_h,output1lo_h,
                     output2hi_h,output2lo_h,outputhi_d,outputlo_d,vrb);
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

double test_cmplx2_polynomial
 ( int dim, int nbr, int nva, int pwr, int deg, int verbose, bool jobrep,
   int mode )
{
   cout << "*** testing the complex arithmetic on complex data ***" << endl;

   if(nbr < 1)
      return 0.0;
   else
   {
      double **inputrehi = new double*[dim]; // dim series of degree deg
      double **inputrelo = new double*[dim]; 
      double **inputimhi = new double*[dim];
      double **inputimlo = new double*[dim];
      for(int i=0; i<dim; i++) 
      {
         inputrehi[i] = new double[deg+1]; inputrelo[i] = new double[deg+1];
         inputimhi[i] = new double[deg+1]; inputimlo[i] = new double[deg+1];
      }
      // The output are dim+1 power series of degree deg
      // for the evaluated and differentiated polynomial.
      double **output1rehi_h = new double*[dim+1];
      double **output1relo_h = new double*[dim+1];
      double **output1imhi_h = new double*[dim+1];
      double **output1imlo_h = new double*[dim+1];
      for(int i=0; i<=dim; i++)
      {
         output1rehi_h[i] = new double[deg+1];
         output1relo_h[i] = new double[deg+1];
         output1imhi_h[i] = new double[deg+1];
         output1imlo_h[i] = new double[deg+1];
      }
      double **output2rehi_h = new double*[dim+1];
      double **output2relo_h = new double*[dim+1];
      double **output2imhi_h = new double*[dim+1];
      double **output2imlo_h = new double*[dim+1];
      for(int i=0; i<=dim; i++)
      {
         output2rehi_h[i] = new double[deg+1];
         output2relo_h[i] = new double[deg+1];
         output2imhi_h[i] = new double[deg+1];
         output2imlo_h[i] = new double[deg+1];
      }
      double **outputrehi_d = new double*[dim+1];
      double **outputrelo_d = new double*[dim+1];
      double **outputimhi_d = new double*[dim+1];
      double **outputimlo_d = new double*[dim+1];
      for(int i=0; i<=dim; i++)
      {
         outputrehi_d[i] = new double[deg+1];
         outputrelo_d[i] = new double[deg+1];
         outputimhi_d[i] = new double[deg+1];
         outputimlo_d[i] = new double[deg+1];
      }
      double *cstrehi = new double[deg+1]; // constant coefficient series
      double *cstrelo = new double[deg+1];
      double *cstimhi = new double[deg+1];
      double *cstimlo = new double[deg+1];
      double **cffrehi = new double*[nbr]; // coefficient series of terms
      double **cffrelo = new double*[nbr];
      double **cffimhi = new double*[nbr];
      double **cffimlo = new double*[nbr];
      for(int i=0; i<nbr; i++)
      {
         cffrehi[i] = new double[deg+1]; cffrelo[i] = new double[deg+1];
         cffimhi[i] = new double[deg+1]; cffimlo[i] = new double[deg+1];
      }
      int *nvr = new int[nbr]; // number of variables in each monomial
      int **idx = new int*[nbr];  // indices of variables in monomials
      int **exp = new int*[nbr];  // exponents of the variables

      bool vrb = (verbose > 1);

      int fail = cmplx2_make_input(dim,nbr,nva,pwr,deg,nvr,idx,exp,
                                   inputrehi,inputrelo,inputimhi,inputimlo,
                                   cstrehi,cstrelo,cstimhi,cstimlo,
                                   cffrehi,cffrelo,cffimhi,cffimlo,vrb);
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
         CPU_cmplx2_poly_evaldiff
            (dim,nbr,deg,nvr,idx,cstrehi,cstrelo,cstimhi,cstimlo,
             cffrehi,cffrelo,cffimhi,cffimlo,
             inputrehi,inputrelo,inputimhi,inputimlo,
             output1rehi_h,output1relo_h,output1imhi_h,output1imlo_h,
             &timelapsec1_h,vrb);
         if(vrb) cout << "Computing with convolution jobs ..." << endl;
         CPU_cmplx2_poly_evaldiffjobs
            (dim,nbr,deg,nvr,idx,cstrehi,cstrelo,cstimhi,cstimlo,
             cffrehi,cffrelo,cffimhi,cffimlo,
             inputrehi,inputrelo,inputimhi,inputimlo,
             output2rehi_h,output2relo_h,output2imhi_h,output2imlo_h,
             cnvjobs,addjobs,&timelapsec2_h,vrb);
      }
      if((mode == 0) || (mode == 2))
      {
         if(vrb) cout << "Computing on the device ..." << endl;
/*
         GPU_cmplx2_poly_evaldiff
            (deg+1,dim,nbr,deg,nvr,idx,
             cstrehi,cstrelo,cstimhi,cstimlo,
             cffrehi,cffrelo,cffimhi,cffimlo,
             inputrehi,inputrelo,inputimhi,inputimlo,
             outputrehi_d,outputrelo_d,outputimhi_d,outputimlo_d,
             cnvjobs,addjobs,&cnvlapms,&addlapms,&timelapms_d,
             &walltimes_d,vrb);
 */
         GPU_cmplx2vectorized_poly_evaldiff
            (deg+1,dim,nbr,deg,nvr,idx,
             cstrehi,cstrelo,cstimhi,cstimlo,cffrehi,cffrelo,cffimhi,cffimlo,
             inputrehi,inputrelo,inputimhi,inputimlo,
             outputrehi_d,outputrelo_d,outputimhi_d,outputimlo_d,
             cmplxcnvjobs,cmplxincjobs,cmplxaddjobs,
             &cnvlapms,&addlapms,&timelapms_d,&walltimes_d,vrb);
      }
      double sumerr = 0.0;
      if(mode == 2)
      {
         sumerr = cmplx2_error_sum(dim,deg,
                     output1rehi_h,output1relo_h,output1imhi_h,output1imlo_h,
                     output2rehi_h,output2relo_h,output2imhi_h,output2imlo_h,
                     outputrehi_d,outputrelo_d,outputimhi_d,outputimlo_d,vrb);
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

int test_dbl2_sequence
 ( int seed, int dim, int nva, int nbr, int pwr, int vrblvl,
   bool jobrep, int mode )
{
   const double tol = 1.0e-24;

   int deg = 0;
   cout << "---> running in double double precision for degree 0 ..."
        << endl;
   int fail = main_dbl2_test_polynomial
                 (seed,dim,nbr,nva,pwr,deg,vrblvl,tol,jobrep,mode);
   deg = 8;
   cout << "---> running for degree 8 ..." << endl;
   fail += main_dbl2_test_polynomial
              (seed,dim,nbr,nva,pwr,deg,vrblvl,tol,false,mode);
   deg = 15;
   cout << "---> running for degree 15 ..." << endl;
   fail += main_dbl2_test_polynomial
              (seed,dim,nbr,nva,pwr,deg,vrblvl,tol,false,mode);
   deg = 31;
   cout << "---> running for degree 31 ..." << endl;
   fail += main_dbl2_test_polynomial
              (seed,dim,nbr,nva,pwr,deg,vrblvl,tol,false,mode);
   deg = 63;
   cout << "---> running for degree 63 ..." << endl;
   fail += main_dbl2_test_polynomial
              (seed,dim,nbr,nva,pwr,deg,vrblvl,tol,false,mode);
   deg = 95;
   cout << "---> running for degree 95 ..." << endl;
   fail += main_dbl2_test_polynomial
              (seed,dim,nbr,nva,pwr,deg,vrblvl,tol,false,mode);
   deg = 127;
   cout << "---> running for degree 127 ..." << endl;
   fail += main_dbl2_test_polynomial
              (seed,dim,nbr,nva,pwr,deg,vrblvl,tol,false,mode);
   deg = 152;
   cout << "---> running for degree 152 ..." << endl;
   fail += main_dbl2_test_polynomial
              (seed,dim,nbr,nva,pwr,deg,vrblvl,tol,false,mode);
   deg = 159;
   cout << "---> running for degree 159 ..." << endl;
   fail += main_dbl2_test_polynomial
              (seed,dim,nbr,nva,pwr,deg,vrblvl,tol,false,mode);
   deg = 191;
   cout << "---> running for degree 191 ..." << endl;
   fail += main_dbl2_test_polynomial
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
