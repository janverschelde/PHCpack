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
#include "dbl2_polynomials_host.h"
#include "dbl2_polynomials_kernels.h"
#include "dbl2_polynomials_testers.h"

using namespace std;

int main_dbl2_test_polynomial
 ( int seed, int dim, int nbr, int nva, int pwr, int deg, int vrblvl,
   double tol )
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

   double realsum = test_dbl2_real_polynomial(dim,nbr,nva,pwr,deg,vrblvl-1);

   int fail = int(realsum > tol);

   if(vrblvl > 0)
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
      cout << "  Seed used : " <<  seedused << endl;
   }
   return fail;
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

double test_dbl2_real_polynomial
 ( int dim, int nbr, int nva, int pwr, int deg, int verbose )
{
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

      make_real2_input(dim,deg,inputhi,inputlo);

      if(verbose > 1)
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
      double *csthi = new double[deg+1]; // constant coefficient series
      double *cstlo = new double[deg+1];
      double **cffhi = new double*[nbr]; // coefficient series of terms
      for(int i=0; i<nbr; i++) cffhi[i] = new double[deg+1];
      double **cfflo = new double*[nbr]; // coefficient series of terms
      for(int i=0; i<nbr; i++) cfflo[i] = new double[deg+1];
      int *nvr = new int[nbr]; // number of variables in each monomial

      if(nva == 0) make_supports(dim,nbr,nvr); // random supports

      int **idx = new int*[nbr];  // indices of variables in monomials

      if(nva == 0)
         for(int i=0; i<nbr; i++) idx[i] = new int[nvr[i]];
      else
      {
         for(int i=0; i<nbr; i++)
         {
            idx[i] = new int[nva];
            nvr[i] = nva;
         }
      }
      int **exp = new int*[nbr];  // exponents of the variables
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
      if(verbose > 1)
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
      bool vrb = (verbose > 1);
      if(nva == 0)
      {
         bool dup = duplicate_supports(dim,nbr,nvr,idx,vrb);
         if(dup)
            cout << "Duplicate supports found." << endl;
         else if(vrb)
            cout << "No duplicate supports found." << endl;
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
      double timelapsec1_h,timelapsec2_h,timelapms_d;

      if(vrb) cout << "Computing without convolution jobs ..." << endl;
      CPU_dbl2_poly_evaldiff
         (dim,nbr,deg,nvr,idx,csthi,cstlo,cffhi,cfflo,inputhi,inputlo,
          output1hi_h,output1lo_h,&timelapsec1_h,vrb);
      if(vrb) cout << "Computing with convolution jobs ..." << endl;
      CPU_dbl2_poly_evaldiffjobs
         (dim,nbr,deg,nvr,idx,csthi,cstlo,cffhi,cfflo,inputhi,inputlo,
          output2hi_h,output2lo_h,cnvjobs,addjobs,&timelapsec2_h,vrb);
      if(vrb) cout << "Computing on the device ..." << endl;
      GPU_dbl2_poly_evaldiff
         (deg+1,dim,nbr,deg,nvr,idx,csthi,cstlo,cffhi,cfflo,
          inputhi,inputlo,outputhi_d,outputlo_d,cnvjobs,addjobs,
          &timelapms_d,vrb);

      double sumerr = dbl2_error_sum(dim,deg,output1hi_h,output1lo_h,
                         output2hi_h,output2lo_h,outputhi_d,outputlo_d,vrb);

      if(verbose > 0)
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

         cout << fixed << setprecision(3);
         cout << "Elapsed CPU time (Linux), Wall time (Windows) : " << endl;
         cout << "  (1) without jobs : " << timelapsec1_h << " seconds,"
              << endl;
         cout << "  (2) cnv/add jobs : " << timelapsec2_h << " seconds."
              << endl;
         cout << "Time spent by all kernels in milliseconds : ";
         cout << fixed << setprecision(2) << timelapms_d << endl;
         cout << scientific << setprecision(16);
      }
      return sumerr;
   }
}
