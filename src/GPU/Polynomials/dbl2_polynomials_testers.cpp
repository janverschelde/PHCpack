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

   double realsum = test_dbl2_real_polynomial
                       (dim,nbr,nva,pwr,deg,vrblvl-1,jobrep,mode);

   int fail = int(realsum > tol);

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
      }
      cout << "  Seed used : " <<  seedused << endl;
   }
   return fail;
}

void dbl2_make_input
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
   if(nva == 0)
   {
      bool dup = duplicate_supports(dim,nbr,nvr,idx,verbose);
      if(dup)
         cout << "Duplicate supports found." << endl;
      else if(verbose)
         cout << "No duplicate supports found." << endl;
   }
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
 ( int dim, int nbr, int nva, int pwr, int deg, int verbose, bool jobrep,
   int mode )
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

      dbl2_make_input(dim,nbr,nva,pwr,deg,nvr,idx,exp,
                      inputhi,inputlo,csthi,cstlo,cffhi,cfflo,vrb);

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
         sumerr = dbl2_error_sum(dim,deg,output1hi_h,output1lo_h,
                     output2hi_h,output2lo_h,outputhi_d,outputlo_d,vrb);

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

int test_dbl2_sequence
 ( int seed, int dim, int nva, int nbr, int pwr, int vrblvl,
   bool jobrep, int mode )
{
   int deg = 15;
   int fail = main_dbl2_test_polynomial
                (seed,dim,nbr,nva,pwr,deg,vrblvl,1.0e-24,jobrep,mode);
   deg = 31;
   cout << "---> running for degree 31 ..." << endl;
   fail += main_dbl2_test_polynomial
              (seed,dim,nbr,nva,pwr,deg,vrblvl,1.0e-24,false,mode);
   deg = 63;
   cout << "---> running for degree 63 ..." << endl;
   fail += main_dbl2_test_polynomial
              (seed,dim,nbr,nva,pwr,deg,vrblvl,1.0e-24,false,mode);
   deg = 95;
   cout << "---> running for degree 95 ..." << endl;
   fail += main_dbl2_test_polynomial
              (seed,dim,nbr,nva,pwr,deg,vrblvl,1.0e-24,false,mode);
   deg = 127;
   cout << "---> running for degree 127 ..." << endl;
   fail += main_dbl2_test_polynomial
              (seed,dim,nbr,nva,pwr,deg,vrblvl,1.0e-24,false,mode);
   deg = 152;
   cout << "---> running for degree 152 ..." << endl;
   fail += main_dbl2_test_polynomial
              (seed,dim,nbr,nva,pwr,deg,vrblvl,1.0e-24,false,mode);
   deg = 159;
   cout << "---> running for degree 159 ..." << endl;
   fail += main_dbl2_test_polynomial
              (seed,dim,nbr,nva,pwr,deg,vrblvl,1.0e-24,false,mode);
   deg = 191;
   cout << "---> running for degree 191 ..." << endl;
   fail += main_dbl2_test_polynomial
              (seed,dim,nbr,nva,pwr,deg,vrblvl,1.0e-24,false,mode);

   if(mode == 2)
   {
      if(fail == 0)
         cout << "All tests passed." << endl;
      else
         cout << "Number of failed tests : " << fail << endl;
   }
   return 0;
}
