/* The file dbl_polynomials_testers.cpp contains the definitions of
 * functions with prototypes in dbl_polynomials_testers.h. */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <vector_types.h>
#include "random_numbers.h"
#include "random_monomials.h"
#include "random_polynomials.h"
#include "convolution_jobs.h"
#include "addition_jobs.h"
#include "write_job_counts.h"
#include "write_gpu_timings.h"
#include "job_makers.h"
#include "dbl_polynomials_host.h"
#include "dbl_polynomials_kernels.h"
#include "dbl_polynomials_testers.h"

using namespace std;

int main_dbl_test_polynomial
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

   double realsum = test_dbl_polynomial
                       (dim,nbr,nva,pwr,deg,vrblvl-1,jobrep,mode);
   double compsum = test_cmplx_polynomial
                       (dim,nbr,nva,pwr,deg,vrblvl-1,jobrep,mode);

   int fail = int(realsum > tol) + int(compsum > tol);

   if(vrblvl > 0)
   {
      if(mode == 2)
      {
         cout << scientific << setprecision(2);
         cout << "Sum of all errors in double precision :" << endl;
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

int dbl_make_input
 ( int dim, int nbr, int nva, int pwr, int deg,
   int *nvr, int **idx, int **exp,
   double **input, double *cst, double **cff, bool verbose )
 { 
   make_real_input(dim,deg,input);

   if(verbose)
   {
      cout << scientific << setprecision(16);
      cout << "Random input series :" << endl;
      for(int i=0; i<dim; i++)
      {
         cout << "-> coefficients of series " << i << " :" << endl;
         for(int j=0; j<=deg; j++) cout << input[i][j] << endl;
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
         make_real_cyclic(dim,nva,deg,idx,cst,cff);
      else
         make_real_products(dim,nbr,nva,deg,idx,cst,cff);
   }
   else
   {
      for(int i=0; i<nbr; i++) exp[i] = new int[nvr[i]];

      bool fail = make_real_polynomial(dim,nbr,pwr,deg,nvr,idx,exp,cst,cff);
   }
   if(verbose)
   {
      cout << "Coefficient series of the constant term :" << endl;
      for(int j=0; j<=deg; j++) cout << cst[j] << endl;

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
         for(int j=0; j<=deg; j++) cout << " " << cff[i][j] << endl;
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

int cmplx_make_input
 ( int dim, int nbr, int nva, int pwr, int deg,
   int *nvr, int **idx, int **exp,
   double **inputre, double **inputim, double *cstre, double *cstim,
   double **cffre, double **cffim, bool verbose )
 { 
   make_complex_input(dim,deg,inputre,inputim);

   if(verbose)
   {
      cout << scientific << setprecision(16);
      cout << "Random input series :" << endl;
      for(int i=0; i<dim; i++)
      {
         cout << "-> coefficients of series " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
            cout << inputre[i][j] << "  " << inputim[i][j] << endl;
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
         make_complex_cyclic(dim,nva,deg,idx,cstre,cstim,cffre,cffim);
      else
         make_complex_products(dim,nbr,nva,deg,idx,cstre,cstim,cffre,cffim);
   }
   else
   {
      for(int i=0; i<nbr; i++) exp[i] = new int[nvr[i]];

      bool fail = make_complex_polynomial(dim,nbr,pwr,deg,nvr,idx,exp,
                                          cstre,cstim,cffre,cffim);
   }
   if(verbose)
   {
      cout << "Coefficient series of the constant term :" << endl;
      for(int j=0; j<=deg; j++)
         cout << cstre[j] << "  " << cstim[j] << endl;

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
            cout << cffre[i][j] << "  " << cffim[i][j] << endl;
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

double dbl_error_sum1
 ( int dim, int deg, double **results_h, double **results_d, bool verbose )
{
   double err = 0.0;

   if(verbose) cout << "The value of the polynomial :" << endl;
   for(int i=0; i<=deg; i++)
   {
      if(verbose)
      {
         cout << results_h[dim][i] << endl;
         cout << results_d[dim][i] << endl;
      }
      err = err + abs(results_h[dim][i] - results_d[dim][i]);
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
            cout << results_h[k][i] << endl;
            cout << results_d[k][i] << endl;
         }
         err = err + abs(results_h[k][i] - results_d[k][i]);
      }
      if(verbose) cout << "error : " << err << endl;
      sumerr = sumerr + err;
   }
   return sumerr;
}

double dbl_error_sum
 ( int dim, int deg, double **results1_h, double **results2_h,
   double **results_d, bool verbose )
{
   double err = 0.0;

   if(verbose) cout << "The value of the polynomial :" << endl;
   for(int i=0; i<=deg; i++)
   {
      if(verbose)
      {
         cout << results1_h[dim][i] << endl;
         cout << results2_h[dim][i] << endl;
         cout << results_d[dim][i] << endl;
      }
      err = err + abs(results1_h[dim][i] - results2_h[dim][i])
                + abs(results1_h[dim][i] - results_d[dim][i]);
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
            cout << results1_h[k][i] << endl;
            cout << results2_h[k][i] << endl;
            cout << results_d[k][i] << endl;
         }
         err = err + abs(results1_h[k][i] - results2_h[k][i])
                   + abs(results1_h[k][i] - results_d[k][i]);
      }
      if(verbose) cout << "error : " << err << endl;
      sumerr = sumerr + err;
   }
   return sumerr;
}

double cmplx_error_sum1
 ( int dim, int deg,
   double **resultsre_h, double **resultsim_h,
   double **resultsre_d, double **resultsim_d, bool verbose )
{
   double err = 0.0;

   if(verbose) cout << "The value of the polynomial :" << endl;
   for(int i=0; i<=deg; i++)
   {
      if(verbose)
      {
         cout << resultsre_h[dim][i] << "  "
              << resultsim_h[dim][i] << endl;
         cout << resultsre_d[dim][i] <<  "  "
              << resultsim_d[dim][i] << endl;
      }
      err = err + abs(resultsre_h[dim][i] - resultsre_d[dim][i])
                + abs(resultsim_h[dim][i] - resultsim_d[dim][i]);
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
            cout << resultsre_h[k][i] << "  "
                 << resultsim_h[k][i] << endl;
            cout << resultsre_d[k][i] << "  "
                 << resultsim_d[k][i] << endl;
         }
         err = err + abs(resultsre_h[k][i] - resultsre_d[k][i])
                   + abs(resultsim_h[k][i] - resultsim_d[k][i]);
      }
      if(verbose) cout << "error : " << err << endl;
      sumerr = sumerr + err;
   }
   return sumerr;
}

double cmplx_error_sum
 ( int dim, int deg,
   double **results1re_h, double **results1im_h,
   double **results2re_h, double **results2im_h,
   double **resultsre_d, double **resultsim_d, bool verbose )
{
   double err = 0.0;

   if(verbose) cout << "The value of the polynomial :" << endl;
   for(int i=0; i<=deg; i++)
   {
      if(verbose)
      {
         cout << results1re_h[dim][i] << "  "
              << results1im_h[dim][i] << endl;
         cout << results2re_h[dim][i] << "  "
              << results2im_h[dim][i] << endl;
         cout << resultsre_d[dim][i] <<  "  "
              << resultsim_d[dim][i] << endl;
      }
      err = err + abs(results1re_h[dim][i] - results2re_h[dim][i])
                + abs(results1im_h[dim][i] - results2im_h[dim][i])
                + abs(results1re_h[dim][i] - resultsre_d[dim][i])
                + abs(results1im_h[dim][i] - resultsim_d[dim][i]);
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
            cout << results1re_h[k][i] << "  "
                 << results1im_h[k][i] << endl;
            cout << results2re_h[k][i] << "  "
                 << results2im_h[k][i] << endl;
            cout << resultsre_d[k][i] << "  "
                 << resultsim_d[k][i] << endl;
         }
         err = err + abs(results1re_h[k][i] - results2re_h[k][i])
                   + abs(results1im_h[k][i] - results2im_h[k][i])
                   + abs(results1re_h[k][i] - resultsre_d[k][i])
                   + abs(results1im_h[k][i] - resultsim_d[k][i]);
      }
      if(verbose) cout << "error : " << err << endl;
      sumerr = sumerr + err;
   }
   return sumerr;
}

double test_dbl_polynomial
 ( int dim, int nbr, int nva, int pwr, int deg, int verbose, bool jobrep,
   int mode )
{
   cout << "*** testing the real arithmetic on real data ***" << endl;

   if(nbr < 1)
      return 0.0;
   else
   {
      double **input = new double*[dim]; // dim series of degree deg
      for(int i=0; i<dim; i++) input[i] = new double[deg+1];
      // The output are dim+1 power series of degree deg
      // for the evaluated and differentiated polynomial.
      double **output1_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output1_h[i] = new double[deg+1];
      double **output2_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output2_h[i] = new double[deg+1];
      double **output_d = new double*[dim+1];
      for(int i=0; i<=dim; i++) output_d[i] = new double[deg+1];

      double *cst = new double[deg+1]; // constant coefficient series
      double **cff = new double*[nbr]; // coefficient series of terms
      for(int i=0; i<nbr; i++) cff[i] = new double[deg+1];
      int *nvr = new int[nbr]; // number of variables in each monomial
      int **idx = new int*[nbr];  // indices of variables in monomials
      int **exp = new int*[nbr];  // exponents of the variables

      const bool vrb = (verbose > 1);

      int fail = dbl_make_input(dim,nbr,nva,pwr,deg,nvr,idx,exp,
                                input,cst,cff,vrb);
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
         CPU_dbl_poly_evaldiff
            (dim,nbr,deg,nvr,idx,cst,cff,input,output1_h,&timelapsec1_h,vrb);
         if(vrb) cout << "Computing with convolution jobs ..." << endl;
         CPU_dbl_poly_evaldiffjobs
            (dim,nbr,deg,nvr,idx,cst,cff,input,output2_h,cnvjobs,addjobs,
             &timelapsec2_h,vrb);
      }
      if((mode == 0) || (mode == 2))
      {
         if(vrb) cout << "Computing on the device ..." << endl;
         GPU_dbl_poly_evaldiff
            (deg+1,dim,nbr,deg,nvr,idx,cst,cff,input,output_d,
             cnvjobs,addjobs,&cnvlapms,&addlapms,&timelapms_d,
             &walltimes_d,vrb);
      }
      double sumerr = 0.0;
      if(mode == 2)
      {
         sumerr = dbl_error_sum(dim,deg,output1_h,output2_h,output_d,vrb);
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

double test_cmplx_polynomial
 ( int dim, int nbr, int nva, int pwr, int deg, int verbose, bool jobrep,
   int mode )
{
   cout << "*** testing the complex arithmetic on complex data ***" << endl;

   if(nbr < 1)
      return 0.0;
   else
   {
      double **inputre = new double*[dim]; // dim series of degree deg
      double **inputim = new double*[dim];
      for(int i=0; i<dim; i++)
      {
         inputre[i] = new double[deg+1];
         inputim[i] = new double[deg+1];
      }
      // The output are dim+1 power series of degree deg
      // for the evaluated and differentiated polynomial.
      double **output1re_h = new double*[dim+1];
      double **output1im_h = new double*[dim+1];
      for(int i=0; i<=dim; i++)
      {
         output1re_h[i] = new double[deg+1];
         output1im_h[i] = new double[deg+1];
      }
      double **output2re_h = new double*[dim+1];
      double **output2im_h = new double*[dim+1];
      for(int i=0; i<=dim; i++)
      {
         output2re_h[i] = new double[deg+1];
         output2im_h[i] = new double[deg+1];
      }
      double **outputre_d = new double*[dim+1];
      double **outputim_d = new double*[dim+1];
      for(int i=0; i<=dim; i++)
      {
         outputre_d[i] = new double[deg+1];
         outputim_d[i] = new double[deg+1];
      }
      double *cstre = new double[deg+1]; // constant coefficient series
      double *cstim = new double[deg+1];
      double **cffre = new double*[nbr]; // coefficient series of terms
      double **cffim = new double*[nbr];
      for(int i=0; i<nbr; i++)
      {
         cffre[i] = new double[deg+1];
         cffim[i] = new double[deg+1];
      }
      int *nvr = new int[nbr]; // number of variables in each monomial
      int **idx = new int*[nbr];  // indices of variables in monomials
      int **exp = new int*[nbr];  // exponents of the variables

      const bool vrb = (verbose > 1);

      int fail = cmplx_make_input(dim,nbr,nva,pwr,deg,nvr,idx,exp,inputre,
                                  inputim,cstre,cstim,cffre,cffim,vrb);
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
         CPU_cmplx_poly_evaldiff
            (dim,nbr,deg,nvr,idx,cstre,cstim,cffre,cffim,inputre,inputim,
             output1re_h,output1im_h,&timelapsec1_h,vrb);
         if(vrb) cout << "Computing with convolution jobs ..." << endl;
         CPU_cmplx_poly_evaldiffjobs
            (dim,nbr,deg,nvr,idx,cstre,cstim,cffre,cffim,inputre,inputim,
             output2re_h,output2im_h,cnvjobs,addjobs,&timelapsec2_h,vrb);
      }
      if((mode == 0) || (mode == 2))
      {
         if(vrb) cout << "Computing on the device ..." << endl;
/*
         GPU_cmplx_poly_evaldiff
            (deg+1,dim,nbr,deg,nvr,idx,cstre,cstim,cffre,cffim,
             inputre,inputim,outputre_d,outputim_d,cnvjobs,addjobs,
             &cnvlapms,&addlapms,&timelapms_d,&walltimes_d,vrb);
 */
         GPU_cmplxvectorized_poly_evaldiff
            (deg+1,dim,nbr,deg,nvr,idx,cstre,cstim,cffre,cffim,
             inputre,inputim,outputre_d,outputim_d,
             cmplxcnvjobs,cmplxincjobs,cmplxaddjobs,
             &cnvlapms,&addlapms,&timelapms_d,&walltimes_d,vrb);
      }
      double sumerr = 0.0;
      if(mode == 2)
      {
         sumerr = cmplx_error_sum(dim,deg,output1re_h,output1im_h,
                     output2re_h,output2im_h,outputre_d,outputim_d,vrb);
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

int test_dbl_sequence
 ( int seed, int dim, int nva, int nbr, int pwr, int vrblvl,
   bool jobrep, int mode )
{
   const double tol = 1.0e-8;

   int deg = 0;
   cout << "---> running in double precision for degree 0 ..." << endl;
   int fail = main_dbl_test_polynomial
                 (seed,dim,nbr,nva,pwr,deg,vrblvl,tol,jobrep,mode);
   deg = 8;
   cout << "---> running for degree 8 ..." << endl;
   fail += main_dbl_test_polynomial
              (seed,dim,nbr,nva,pwr,deg,vrblvl,tol,false,mode);
   deg = 15;
   cout << "---> running for degree 15 ..." << endl;
   fail += main_dbl_test_polynomial
              (seed,dim,nbr,nva,pwr,deg,vrblvl,tol,false,mode);
   deg = 31;
   cout << "---> running for degree 31 ..." << endl;
   fail += main_dbl_test_polynomial
              (seed,dim,nbr,nva,pwr,deg,vrblvl,tol,false,mode);
   deg = 63;
   cout << "---> running for degree 63 ..." << endl;
   fail += main_dbl_test_polynomial
              (seed,dim,nbr,nva,pwr,deg,vrblvl,tol,false,mode);
   deg = 95;
   cout << "---> running for degree 95 ..." << endl;
   fail += main_dbl_test_polynomial
              (seed,dim,nbr,nva,pwr,deg,vrblvl,tol,false,mode);
   deg = 127;
   cout << "---> running for degree 127 ..." << endl;
   fail += main_dbl_test_polynomial
              (seed,dim,nbr,nva,pwr,deg,vrblvl,tol,false,mode);
   deg = 152;
   cout << "---> running for degree 152 ..." << endl;
   fail += main_dbl_test_polynomial
              (seed,dim,nbr,nva,pwr,deg,vrblvl,tol,false,mode);
   deg = 159;
   cout << "---> running for degree 159 ..." << endl;
   fail += main_dbl_test_polynomial
              (seed,dim,nbr,nva,pwr,deg,vrblvl,tol,false,mode);
   deg = 191;
   cout << "---> running for degree 191 ..." << endl;
   fail += main_dbl_test_polynomial
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
