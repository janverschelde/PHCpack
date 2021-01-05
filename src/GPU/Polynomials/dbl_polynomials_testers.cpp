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
#include "dbl_polynomials_host.h"
#include "dbl_polynomials_kernels.h"
#include "dbl_polynomials_testers.h"

using namespace std;

int main_dbl_test_polynomial
 ( int seed, int dim, int nbr, int nva, int pwr, int deg, int vrblvl )
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

   double realsum = test_dbl_real_polynomial(dim,nbr,nva,pwr,deg,vrblvl-1);

   const double tol = 1.0e-12;

   int fail = int(realsum > tol);

   if(vrblvl > 0)
   {
      cout << scientific << setprecision(2);
      cout << "Sum of all errors in double precision :" << endl;
      cout << "  on real data : " << realsum;
      if(realsum < tol)
         cout << "  pass." << endl;
      else
         cout << "  fail!" << endl;

      cout << "  Seed used : " <<  seedused << endl;
   }
   return fail;
}

double test_dbl_real_polynomial
 ( int dim, int nbr, int nva, int pwr, int deg, int verbose )
{
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

      make_real_input(dim,deg,input);

      if(verbose > 1)
      {
         cout << scientific << setprecision(16);
         cout << "Random input series :" << endl;
         for(int i=0; i<dim; i++)
         {
            cout << "-> coefficients of series " << i << " :" << endl;
            for(int j=0; j<=deg; j++) cout << input[i][j] << endl;
         }
      }
      double *cst = new double[deg+1]; // constant coefficient series
      double **cff = new double*[nbr]; // coefficient series of terms
      for(int i=0; i<nbr; i++) cff[i] = new double[deg+1];
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
            make_real_cyclic(dim,nva,deg,idx,cst,cff);
         else
            make_real_products(dim,nbr,nva,deg,idx,cst,cff);
      }
      else
      {
         for(int i=0; i<nbr; i++) exp[i] = new int[nvr[i]];

         bool fail = make_real_polynomial(dim,nbr,pwr,deg,nvr,idx,exp,cst,cff);
      }
      if(verbose > 1)
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
      CPU_dbl_poly_evaldiff
         (dim,nbr,deg,nvr,idx,cst,cff,input,output1_h,&timelapsec1_h,vrb);
      if(vrb) cout << "Computing with convolution jobs ..." << endl;
      CPU_dbl_poly_evaldiffjobs
         (dim,nbr,deg,nvr,idx,cst,cff,input,output2_h,cnvjobs,addjobs,
          &timelapsec2_h,vrb);
      if(vrb) cout << "Computing on the device ..." << endl;
      GPU_dbl_poly_evaldiff
         (deg+1,dim,nbr,deg,nvr,idx,cst,cff,input,output_d,
          cnvjobs,addjobs,&timelapms_d,vrb);

      double err = 0.0;

      if(vrb) cout << "The value of the polynomial :" << endl;
      for(int i=0; i<=deg; i++)
      {
         if(vrb)
         {
            cout << output1_h[dim][i] << endl;
            cout << output2_h[dim][i] << endl;
            cout << output_d[dim][i] << endl;
         }
         err = err + abs(output1_h[dim][i] - output2_h[dim][i])
                   + abs(output1_h[dim][i] - output_d[dim][i]);
      }
      if(vrb) cout << "error : " << err << endl;

      double sumerr = err;

      for(int k=0; k<dim; k++)
      {
         if(vrb) cout << "Derivative " << k << " :" << endl;
         err = 0.0;
         for(int i=0; i<=deg; i++)
         {
            if(vrb)
            {
               cout << output1_h[k][i] << endl;
               cout << output2_h[k][i] << endl;
               cout << output_d[k][i] << endl;
            }
            err = err + abs(output1_h[k][i] - output2_h[k][i])
                      + abs(output1_h[k][i] - output_d[k][i]);
         }
         if(vrb) cout << "error : " << err << endl;
         sumerr = sumerr + err;
      }
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
         cout << "Time spent by all kernels : ";
         cout << fixed << setprecision(2) << timelapms_d
              << " milliseconds." << endl;
         cout << scientific << setprecision(16);
      } 
      return sumerr;
   }
}
