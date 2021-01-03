/* Tests polynomial evaluation and differentiation
 * in deca double precision. */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <cmath>
// #include <vector_types.h>
#include "random_polynomials.h"
#include "random10_monomials.h"
#include "random10_polynomials.h"
#include "convolution_jobs.h"
#include "addition_jobs.h"
#include "write_job_counts.h"
// #include "dbl5_polynomials_host.h"
// #include "dbl5_polynomials_kernels.h"

using namespace std;

double test_dbl10_real_polynomial
 ( int dim, int nbr, int nva, int pwr, int deg, int verbose );
/*
 * DESCRIPTION :
 *   Tests the evaluation and differentiation for random real data.
 *   Returns the sum of all errors.
 * 
 * ON ENTRY :
 *   dim      dimension, total number of variables;
 *   nbr      number of terms in the polynomial;
 *   nva      number of variables per monomial (for products and cyclic);
 *   pwr      highest power of each variable;
 *   deg      truncation degree of the series;
 *   verbose  if zero, then no output is written. */

int main_dbl10_test_polynomial
 ( int seed, int dim, int nbr, int nva, int pwr, int deg, int vrblvl );
/*
 * DESCRIPTION :
 *   Runs tests on a random polynomial in triple double precision.
 *   Returns 0 if all tests passed,
 *   otherwise, returns the number of failed tests.
 *
 * ON ENTRY :
 *   seed     seed for the random number generator;
 *   dim      dimension, total number of variables;
 *   nbr      number of terms in the polynomial;
 *   nva      number of variables in each monomial (for products, cyclic);
 *   pwr      highest power of each variable;
 *   deg      truncation degree of the series;
 *   vrblvl   is the verbose level, if 0 then no output. */

int main ( void )
{
   cout << "Give the seed (0 for time) : ";
   int seed; cin >> seed;

   cout << "Give the dimension : ";
   int dim;  cin >> dim;

   cout << "Give the variables per monomial (0 for random polynomial) : ";
   int nva; cin >> nva;

   int nbr; // number of monomials, not counting the constant

   if(nva > 0)
   {
      cout << "Enter 0 for products, other number of cyclic : ";
      cin >> nbr;

      if(nbr == 0)
         nbr = products_count(dim,nva);
      else
         nbr = dim;

      cout << "-> number of monomials : " << nbr << endl;
   }
   else
   {
      cout << "Give the number of terms : ";
      cin >> nbr;
   }
   // cout << "Give the largest power of each variable : "; cin >> pwr;
   const int pwr=1;

   cout << "Give the degree of the series : ";
   int deg; cin >> deg;

   cout << "Give the verbose level : ";
   int vrb; cin >> vrb;

   int fail = main_dbl10_test_polynomial(seed,dim,nbr,nva,pwr,deg,vrb);

   if(fail == 0)
      cout << "All tests passed." << endl;
   else
      cout << "Number of failed tests : " << fail << endl;

   return 0;
}

int main_dbl10_test_polynomial
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

   double realsum = test_dbl10_real_polynomial(dim,nbr,nva,pwr,deg,vrblvl-1);

   const double tol = 1.0e-152;

   int fail = int(realsum > tol);

   if(vrblvl > 0)
   {
      cout << scientific << setprecision(2);
      cout << "Sum of all errors in deca double precision :" << endl;
      cout << "  on real data : " << realsum;
      if(realsum < tol)
         cout << "  pass." << endl;
      else
         cout << "  fail!" << endl;

      cout << "  Seed used : " <<  seedused << endl;
   }
   return fail;
}

double test_dbl10_real_polynomial
 ( int dim, int nbr, int nva, int pwr, int deg, int verbose )
{
   if(nbr < 1)
      return 0.0;
   else
   {
      double **inputrtb = new double*[dim]; // dim series of degree deg
      for(int i=0; i<dim; i++) inputrtb[i] = new double[deg+1];
      double **inputrix = new double*[dim];
      for(int i=0; i<dim; i++) inputrix[i] = new double[deg+1];
      double **inputrmi = new double*[dim];
      for(int i=0; i<dim; i++) inputrmi[i] = new double[deg+1];
      double **inputrrg = new double*[dim];
      for(int i=0; i<dim; i++) inputrrg[i] = new double[deg+1];
      double **inputrpk = new double*[dim];
      for(int i=0; i<dim; i++) inputrpk[i] = new double[deg+1];
      double **inputltb = new double*[dim];
      for(int i=0; i<dim; i++) inputltb[i] = new double[deg+1];
      double **inputlix = new double*[dim];
      for(int i=0; i<dim; i++) inputlix[i] = new double[deg+1];
      double **inputlmi = new double*[dim];
      for(int i=0; i<dim; i++) inputlmi[i] = new double[deg+1];
      double **inputlrg = new double*[dim];
      for(int i=0; i<dim; i++) inputlrg[i] = new double[deg+1];
      double **inputlpk = new double*[dim];
      for(int i=0; i<dim; i++) inputlpk[i] = new double[deg+1];
      // The output are dim+1 power series of degree deg
      // for the evaluated and differentiated polynomial.
      double **output1rtb_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output1rtb_h[i] = new double[deg+1];
      double **output1rix_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output1rix_h[i] = new double[deg+1];
      double **output1rmi_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output1rmi_h[i] = new double[deg+1];
      double **output1rrg_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output1rrg_h[i] = new double[deg+1];
      double **output1rpk_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output1rpk_h[i] = new double[deg+1];
      double **output1ltb_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output1ltb_h[i] = new double[deg+1];
      double **output1lix_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output1lix_h[i] = new double[deg+1];
      double **output1lmi_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output1lmi_h[i] = new double[deg+1];
      double **output1lrg_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output1lrg_h[i] = new double[deg+1];
      double **output1lpk_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output1lpk_h[i] = new double[deg+1];
      double **output2rtb_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output2rtb_h[i] = new double[deg+1];
      double **output2rix_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output2rix_h[i] = new double[deg+1];
      double **output2rmi_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output2rmi_h[i] = new double[deg+1];
      double **output2rrg_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output2rrg_h[i] = new double[deg+1];
      double **output2rpk_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output2rpk_h[i] = new double[deg+1];
      double **output2ltb_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output2ltb_h[i] = new double[deg+1];
      double **output2lix_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output2lix_h[i] = new double[deg+1];
      double **output2lmi_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output2lmi_h[i] = new double[deg+1];
      double **output2lrg_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output2lrg_h[i] = new double[deg+1];
      double **output2lpk_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output2lpk_h[i] = new double[deg+1];
      double **outputrtb_d = new double*[dim+1];
      for(int i=0; i<=dim; i++) outputrtb_d[i] = new double[deg+1];
      double **outputrix_d = new double*[dim+1];
      for(int i=0; i<=dim; i++) outputrix_d[i] = new double[deg+1];
      double **outputrmi_d = new double*[dim+1];
      for(int i=0; i<=dim; i++) outputrmi_d[i] = new double[deg+1];
      double **outputrrg_d = new double*[dim+1];
      for(int i=0; i<=dim; i++) outputrrg_d[i] = new double[deg+1];
      double **outputrpk_d = new double*[dim+1];
      for(int i=0; i<=dim; i++) outputrpk_d[i] = new double[deg+1];
      double **outputltb_d = new double*[dim+1];
      for(int i=0; i<=dim; i++) outputltb_d[i] = new double[deg+1];
      double **outputlix_d = new double*[dim+1];
      for(int i=0; i<=dim; i++) outputlix_d[i] = new double[deg+1];
      double **outputlmi_d = new double*[dim+1];
      for(int i=0; i<=dim; i++) outputlmi_d[i] = new double[deg+1];
      double **outputlrg_d = new double*[dim+1];
      for(int i=0; i<=dim; i++) outputlrg_d[i] = new double[deg+1];
      double **outputlpk_d = new double*[dim+1];
      for(int i=0; i<=dim; i++) outputlpk_d[i] = new double[deg+1];

      make_real10_input(dim,deg,
         inputrtb,inputrix,inputrmi,inputrrg,inputrpk,
         inputltb,inputlix,inputlmi,inputlrg,inputlpk);

      if(verbose > 0)
      {
         cout << scientific << setprecision(16);
         cout << "Random input series :" << endl;
         for(int i=0; i<dim; i++)
         {
            cout << "-> coefficients of series " << i << " :" << endl;
            for(int j=0; j<=deg; j++)
            {
               cout << inputrtb[i][j] << "  " << inputrix[i][j]
                                      << "  " << inputrmi[i][j] << endl
                    << inputrrg[i][j] << "  " << inputrpk[i][j] << endl;
               cout << inputltb[i][j] << "  " << inputlix[i][j]
                                      << "  " << inputlmi[i][j] << endl
                    << inputlrg[i][j] << "  " << inputlpk[i][j] << endl;
            }
         }
      }
      double *cstrtb = new double[deg+1]; // constant coefficient series
      double *cstrix = new double[deg+1];
      double *cstrmi = new double[deg+1];
      double *cstrrg = new double[deg+1];
      double *cstrpk = new double[deg+1];
      double *cstltb = new double[deg+1];
      double *cstlix = new double[deg+1];
      double *cstlmi = new double[deg+1];
      double *cstlrg = new double[deg+1];
      double *cstlpk = new double[deg+1];
      double **cffrtb = new double*[nbr]; // coefficient series of terms
      for(int i=0; i<nbr; i++) cffrtb[i] = new double[deg+1];
      double **cffrix = new double*[nbr];
      for(int i=0; i<nbr; i++) cffrix[i] = new double[deg+1];
      double **cffrmi = new double*[nbr];
      for(int i=0; i<nbr; i++) cffrmi[i] = new double[deg+1];
      double **cffrrg = new double*[nbr];
      for(int i=0; i<nbr; i++) cffrrg[i] = new double[deg+1];
      double **cffrpk = new double*[nbr];
      for(int i=0; i<nbr; i++) cffrpk[i] = new double[deg+1];
      double **cffltb = new double*[nbr];
      for(int i=0; i<nbr; i++) cffltb[i] = new double[deg+1];
      double **cfflix = new double*[nbr];
      for(int i=0; i<nbr; i++) cfflix[i] = new double[deg+1];
      double **cfflmi = new double*[nbr];
      for(int i=0; i<nbr; i++) cfflmi[i] = new double[deg+1];
      double **cfflrg = new double*[nbr];
      for(int i=0; i<nbr; i++) cfflrg[i] = new double[deg+1];
      double **cfflpk = new double*[nbr];
      for(int i=0; i<nbr; i++) cfflpk[i] = new double[deg+1];
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
            make_real10_cyclic
               (dim,nva,deg,idx,
                cstrtb,cstrix,cstrmi,cstrrg,cstrpk,
                cstltb,cstlix,cstlmi,cstlrg,cstlpk,
                cffrtb,cffrix,cffrmi,cffrrg,cffrpk,
                cffltb,cfflix,cfflmi,cfflrg,cfflpk);
         else
            make_real10_products
               (dim,nbr,nva,deg,idx,
                cstrtb,cstrix,cstrmi,cstrrg,cstrpk,
                cstltb,cstlix,cstlmi,cstlrg,cstlpk,
                cffrtb,cffrix,cffrmi,cffrrg,cffrpk,
                cffltb,cfflix,cfflmi,cfflrg,cfflpk);
      }
      else
      {
         for(int i=0; i<nbr; i++) exp[i] = new int[nvr[i]];

         bool fail = make_real10_polynomial
                        (dim,nbr,pwr,deg,nvr,idx,exp,
                         cstrtb,cstrix,cstrmi,cstrrg,cstrpk,
                         cstltb,cstlix,cstlmi,cstlrg,cstlpk,
                         cffrtb,cffrix,cffrmi,cffrrg,cffrpk,
                         cffltb,cfflix,cfflmi,cfflrg,cfflpk);
      }
      if(verbose > 0)
      {
         cout << "Coefficient series of the constant term :" << endl;
         for(int j=0; j<=deg; j++)
         {
            cout << cstrtb[j] << "  " << cstrix[j]
                              << "  " << cstrmi[j] << endl
                 << cstrrg[j] << "  " << cstrpk[j] << endl;
            cout << cstltb[j] << "  " << cstlix[j]
                              << "  " << cstlmi[j] << endl
                 << cstlrg[j] << "  " << cstlpk[j] << endl;
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
               cout << cffrtb[i][j] << "  " << cffrix[i][j]
                                    << "  " << cffrmi[i][j] << endl
                    << cffrrg[i][j] << "  " << cffrpk[i][j] << endl;
               cout << cffltb[i][j] << "  " << cfflix[i][j]
                                    << "  " << cfflmi[i][j] << endl
                    << cfflrg[i][j] << "  " << cfflpk[i][j] << endl;
            }
         }
      }
      bool vrb = (verbose > 0);
      if(nva == 0)
      {
         bool dup = duplicate_supports(dim,nbr,nvr,idx,vrb);
         if(dup)
            cout << "Duplicate supports found." << endl;
         else
            cout << "No duplicate supports found." << endl;
      }
 /*
      ConvolutionJobs cnvjobs(dim);

      cnvjobs.make(nbr,nvr,idx,vrb);

      write_convolution_counts(cnvjobs);

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
      write_addition_counts(addjobs);

      for(int k=0; k<addjobs.get_depth(); k++)
      {
         cout << "jobs at layer " << k << " :" << endl;
         for(int i=0; i<addjobs.get_layer_count(k); i++)
            cout << addjobs.get_job(k,i) << endl;
      }
      if(vrb) cout << "Computing without convolution jobs ..." << endl;
      CPU_dbl10_poly_evaldiff
         (dim,nbr,deg,nvr,idx,csttb,cstix,cstrg,cstpk,
          cfftb,cffix,cffrg,cffpk,
          inputtb,inputix,inputrg,inputpk,
          output1tb_h,output1ix_h,output1rg_h,output1pk_h,vrb);
      if(vrb) cout << "Computing with convolution jobs ..." << endl;
      CPU_dbl10_poly_evaldiffjobs
         (dim,nbr,deg,nvr,idx,csttb,cstix,cstrg,cstpk,
          cfftb,cffix,cffrg,cffpk,
          inputtb,inputix,inputrg,inputpk,
          output2tb_h,output2ix_h,output2rg_h,output2pk_h,
          cnvjobs,addjobs,vrb);
      if(vrb) cout << "Computing on the device ..." << endl;
      GPU_dbl10_poly_evaldiff
         (deg+1,dim,nbr,deg,nvr,idx,csttb,cstix,cstrg,cstpk,
          cfftb,cffix,cffrg,cffpk,
          inputtb,inputix,inputrg,inputpk,
          outputtb_d,outputix_d,outputrg_d,outputpk_d,
          cnvjobs,addjobs,vrb);

      double err = 0.0;

      if(verbose > 0) cout << "The value of the polynomial :" << endl;
      for(int i=0; i<=deg; i++)
      {
         if(verbose > 0)
         {
            cout << output1tb_h[dim][i] << "  "
                 << output1ix_h[dim][i] << endl
                 << output1rg_h[dim][i] << "  "
                 << output1pk_h[dim][i] << endl;
            cout << output2tb_h[dim][i] << "  "
                 << output2ix_h[dim][i] << endl
                 << output2rg_h[dim][i] << "  "
                 << output2pk_h[dim][i] << endl;
            cout << outputtb_d[dim][i] << "  "
                 << outputix_d[dim][i] << endl
                 << outputrg_d[dim][i] << "  "
                 << outputpk_d[dim][i] << endl;
         }
         err = err + abs(output1tb_h[dim][i] - output2tb_h[dim][i])
                   + abs(output1ix_h[dim][i] - output2ix_h[dim][i])
                   + abs(output1rg_h[dim][i] - output2rg_h[dim][i])
                   + abs(output1pk_h[dim][i] - output2pk_h[dim][i])
                   + abs(output1tb_h[dim][i] - outputtb_d[dim][i])
                   + abs(output1ix_h[dim][i] - outputix_d[dim][i])
                   + abs(output1rg_h[dim][i] - outputrg_d[dim][i])
                   + abs(output1pk_h[dim][i] - outputpk_d[dim][i]);
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
               cout << output1tb_h[k][i] << "  "
                    << output1ix_h[k][i] << endl
                    << output1rg_h[k][i] << "  "
                    << output1pk_h[k][i] << endl;
               cout << output2tb_h[k][i] << "  "
                    << output2ix_h[k][i] << endl
                    << output2rg_h[k][i] << "  "
                    << output2pk_h[k][i] << endl;
               cout << outputtb_d[k][i] << "  "
                    << outputix_d[k][i] << endl
                    << outputrg_d[k][i] << "  "
                    << outputpk_d[k][i] << endl;
            }
            err = err + abs(output1tb_h[k][i] - output2tb_h[k][i])
                      + abs(output1ix_h[k][i] - output2ix_h[k][i])
                      + abs(output1rg_h[k][i] - output2rg_h[k][i])
                      + abs(output1pk_h[k][i] - output2pk_h[k][i])
                      + abs(output1tb_h[k][i] - outputtb_d[k][i])
                      + abs(output1ix_h[k][i] - outputix_d[k][i])
                      + abs(output1rg_h[k][i] - outputrg_d[k][i])
                      + abs(output1pk_h[k][i] - outputpk_d[k][i]);
         }
         if(verbose > 0) cout << "error : " << err << endl;
         sumerr = sumerr + err;
      }
      cout << "dimension : " << dim << endl;
      if(nva > 0)
      {
         cout << "number of variables per monomial : " << nva << endl;
      }
      cout << "number of monomials : " << nbr << endl;
      write_convolution_counts(cnvjobs);
      write_addition_counts(addjobs);
      write_operation_counts(deg,cnvjobs,addjobs);

      return sumerr;
  */
      return 0.0;
   }
}
