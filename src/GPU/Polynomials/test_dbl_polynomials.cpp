/* Tests polynomial evaluation and differentiation in double precision. */

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
#include "dbl_polynomials_host.h"
#include "dbl_polynomials_kernels.h"

using namespace std;

double test_dbl_real_polynomial
 ( int dim, int nbr, int pwr, int deg, int verbose );
/*
 * DESCRIPTION :
 *   Tests the evaluation and differentiation for random real data.
 *   Returns the sum of all errors.
 * 
 * ON ENTRY :
 *   dim      dimension, total number of variables;
 *   nbr      number of terms in the polynomial;
 *   pwr      highest power of each variable;
 *   deg      truncation degree of the series;
 *   verbose  if zero, then no output is written. */

int main_dbl_test_polynomial
 ( int seed, int dim, int nbr, int pwr, int deg, int vrblvl );
/*
 * DESCRIPTION :
 *   Runs tests on a random polynomial in double precision.
 *   Returns 0 if all tests passed,
 *   otherwise, returns the number of failed tests.
 *
 * ON ENTRY :
 *   seed     seed for the random number generator;
 *   dim      dimension, total number of variables;
 *   nbr      number of terms in the polynomial;
 *   pwr      highest power of each variable;
 *   deg      truncation degree of the series;
 *   vrblvl   is the verbose level, if 0 then no output. */

int main ( void )
{
   cout << "Give the seed (0 for time) : ";
   int seed; cin >> seed;

   cout << "Give the dimension : ";
   int dim;  cin >> dim;

   cout << "Give the number of terms : ";
   int nbr; cin >> nbr;

   // cout << "Give the largest power of each variable : "; cin >> pwr;
   const int pwr=1;

   cout << "Give the degree of the series : ";
   int deg; cin >> deg;

   cout << "Give the verbose level : ";
   int vrb; cin >> vrb;

   int fail = main_dbl_test_polynomial(seed,dim,nbr,pwr,deg,vrb);

   if(fail == 0)
      cout << "All tests passed." << endl;
   else
      cout << "Number of failed tests : " << fail << endl;

   return 0;
}

int main_dbl_test_polynomial
 ( int seed, int dim, int nbr, int pwr, int deg, int vrblvl )
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

   double realsum = test_dbl_real_polynomial(dim,nbr,pwr,deg,vrblvl-1);

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
 ( int dim, int nbr, int pwr, int deg, int verbose )
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

      if(verbose > 0)
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

      make_supports(dim,nbr,nvr); // define supports of polynomial

      int **idx = new int*[nbr];  // indices of variables in monomials
      for(int i=0; i<nbr; i++) idx[i] = new int[nvr[i]];
      int **exp = new int*[nbr];  // exponents of the variables
      for(int i=0; i<nbr; i++) exp[i] = new int[nvr[i]];

      bool fail = make_real_polynomial(dim,nbr,pwr,deg,nvr,idx,exp,cst,cff);

      if(verbose > 0)
      {
         cout << "Coefficient series of the constant term :" << endl;
         for(int j=0; j<=deg; j++) cout << cst[j] << endl;

         for(int i=0; i<nbr; i++)
         {
            cout << "Generated random monomial " << i << " :" << endl;
            cout << "   the indices :";
            for(int j=0; j<nvr[i]; j++) cout << " " << idx[i][j];
            cout << endl;
            cout << " the exponents :";
            for(int j=0; j<nvr[i]; j++) cout << " " << exp[i][j];
            cout << endl;
            cout << " coefficient series :" << endl;
            for(int j=0; j<=deg; j++) cout << " " << cff[i][j] << endl;
         }
      }
      bool vrb = (verbose > 0);
      bool dup = duplicate_supports(dim,nbr,nvr,idx,vrb);
      if(dup)
         cout << "Duplicate supports found." << endl;
      else
         cout << "No duplicate supports found." << endl;

      ConvolutionJobs cnvjobs(dim);

      cnvjobs.make(nbr,nvr,idx,vrb);

      cout << "number of convolution jobs : " << cnvjobs.get_count() << endl;
      cout << "number of layers : " << cnvjobs.get_depth() << endl;
      cout << "frequency of layer counts :" << endl;
      int checksum = 0;
      for(int i=0; i<cnvjobs.get_depth(); i++)
      {
         cout << i << " : " << cnvjobs.get_layer_count(i) << endl;
         checksum = checksum + cnvjobs.get_layer_count(i); 
      }
      cout << "layer count sum : " << checksum << endl;

      for(int k=0; k<cnvjobs.get_depth(); k++)
      {
         cout << "jobs at layer " << k << " :" << endl;
         for(int i=0; i<cnvjobs.get_layer_count(k); i++)
            cout << cnvjobs.get_job(k,i) << endl;
      }
      cout << endl;

      AdditionJobs addjobs(dim,nbr);

      addjobs.make(nbr,nvr,idx,true);

      cout << "The differential indices :" << endl;
      for(int i=0; i<dim; i++)
      {
         cout << "variable " << i << " :";
         for(int j=0; j<=addjobs.get_differential_count(i); j++)
            cout << " " << addjobs.get_differential_index(i,j);
         cout << endl;
      }
      cout << "number of addition jobs : " << addjobs.get_count() << endl;
      cout << "number of layers : " << addjobs.get_depth() << endl;
      cout << "frequency of layer counts :" << endl;
      checksum = 0;
      for(int i=0; i<addjobs.get_depth(); i++)
      {
         cout << i << " : " << addjobs.get_layer_count(i) << endl;
         checksum = checksum + addjobs.get_layer_count(i); 
      }
      cout << "layer count sum : " << checksum << endl;

      for(int k=0; k<addjobs.get_depth(); k++)
      {
         cout << "jobs at layer " << k << " :" << endl;
         for(int i=0; i<addjobs.get_layer_count(k); i++)
            cout << addjobs.get_job(k,i) << endl;
      }
      if(vrb) cout << "Computing without convolution jobs ..." << endl;
      CPU_dbl_poly_evaldiff(dim,nbr,deg,nvr,idx,cst,cff,input,output1_h,vrb);
      if(vrb) cout << "Computing with convolution jobs ..." << endl;
      CPU_dbl_poly_evaldiffjobs
         (dim,nbr,deg,nvr,idx,cst,cff,input,output2_h,cnvjobs,addjobs,vrb);
      if(vrb) cout << "Computing on the device ..." << endl;
      GPU_dbl_poly_evaldiff
         (deg+1,dim,nbr,deg,nvr,idx,cst,cff,input,output_d,cnvjobs,vrb);

      double err = 0.0;

      if(verbose > 0) cout << "The value of the polynomial :" << endl;
      for(int i=0; i<=deg; i++)
      {
         if(verbose > 0)
         {
            cout << output1_h[dim][i] << endl;
            cout << output2_h[dim][i] << endl;
            cout << output_d[dim][i] << endl;
         }
         err = err + abs(output1_h[dim][i] - output2_h[dim][i])
                   + abs(output1_h[dim][i] - output_d[dim][i]);
      }
      if(verbose > 0) cout << "error : " << err << endl;

      double sumerr = err;

      for(int k=0; k<dim; k++)
      {
         if(verbose > 0) cout << "Derivative " << k << " :" << endl;
         err = 0.0;
         for(int i=0; i<=deg; i++)
         {
            if(verbose > 0)
            {
               cout << output1_h[k][i] << endl;
               cout << output2_h[k][i] << endl;
               cout << output_d[k][i] << endl;
            }
            err = err + abs(output1_h[k][i] - output2_h[k][i])
                      + abs(output1_h[k][i] - output_d[k][i]);
         }
         if(verbose > 0) cout << "error : " << err << endl;
         sumerr = sumerr + err;
      }
      return sumerr;
   }
}
