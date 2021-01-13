/* The file dbl10_polynomials_testers.cpp contains the definitions of
 * functions with prototypes in dbl10_polynomials_testers.h. */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <vector_types.h>
#include "random_polynomials.h"
#include "random10_monomials.h"
#include "random10_polynomials.h"
#include "convolution_jobs.h"
#include "addition_jobs.h"
#include "write_job_counts.h"
#include "dbl10_polynomials_host.h"
#include "dbl10_polynomials_kernels.h"
#include "dbl10_polynomials_testers.h"

using namespace std;

int main_dbl10_test_polynomial
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

   double realsum = test_dbl10_real_polynomial
                       (dim,nbr,nva,pwr,deg,vrblvl-1,jobrep,mode);

   int fail = int(realsum > tol);

   if(vrblvl > 0)
   {
      if(mode == 2)
      {
         cout << scientific << setprecision(2);
         cout << "Sum of all errors in deca double precision :" << endl;
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

void dbl10_make_input
 ( int dim, int nbr, int nva, int pwr, int deg,
   int *nvr, int **idx, int **exp,
   double **inputrtb, double **inputrix, double **inputrmi,
   double **inputrrg, double **inputrpk,
   double **inputltb, double **inputlix, double **inputlmi,
   double **inputlrg, double **inputlpk,
   double *cstrtb, double *cstrix, double *cstrmi, 
   double *cstrrg, double *cstrpk,
   double *cstltb, double *cstlix, double *cstlmi, 
   double *cstlrg, double *cstlpk,
   double **cffrtb, double **cffrix, double **cffrmi, 
   double **cffrrg, double **cffrpk,
   double **cffltb, double **cfflix, double **cfflmi, 
   double **cfflrg, double **cfflpk, bool verbose )
{
   make_real10_input(dim,deg,
      inputrtb,inputrix,inputrmi,inputrrg,inputrpk,
      inputltb,inputlix,inputlmi,inputlrg,inputlpk);

   if(verbose)
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
   if(nva == 0) // random supports
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
   if(verbose)
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
   if(nva == 0)
   {
      bool dup = duplicate_supports(dim,nbr,nvr,idx,verbose);
      if(dup)
         cout << "Duplicate supports found." << endl;
      else if(verbose)
         cout << "No duplicate supports found." << endl;
   }
}

double dbl10_error_sum
 ( int dim, int deg,
   double **results1rtb_h, double **results1rix_h, double **results1rmi_h, 
   double **results1rrg_h, double **results1rpk_h,
   double **results1ltb_h, double **results1lix_h, double **results1lmi_h, 
   double **results1lrg_h, double **results1lpk_h,
   double **results2rtb_h, double **results2rix_h, double **results2rmi_h, 
   double **results2rrg_h, double **results2rpk_h,
   double **results2ltb_h, double **results2lix_h, double **results2lmi_h, 
   double **results2lrg_h, double **results2lpk_h,
   double **resultsrtb_d, double **resultsrix_d, double **resultsrmi_d, 
   double **resultsrrg_d, double **resultsrpk_d,
   double **resultsltb_d, double **resultslix_d, double **resultslmi_d, 
   double **resultslrg_d, double **resultslpk_d, bool verbose )
{
   double err = 0.0;

   if(verbose) cout << "The value of the polynomial :" << endl;
   for(int i=0; i<=deg; i++)
   {
      if(verbose)
      {
         cout << results1rtb_h[dim][i] << "  "
              << results1rix_h[dim][i] << "  "
              << results1rmi_h[dim][i] << endl
              << results1rrg_h[dim][i] << "  "
              << results1rpk_h[dim][i] << endl;
         cout << results1ltb_h[dim][i] << "  "
              << results1lix_h[dim][i] << "  "
              << results1lmi_h[dim][i] << endl
              << results1lrg_h[dim][i] << "  "
              << results1lpk_h[dim][i] << endl;
         cout << results2rtb_h[dim][i] << "  "
              << results2rix_h[dim][i] << "  "
              << results2rmi_h[dim][i] << endl
              << results2rrg_h[dim][i] << "  "
              << results2rpk_h[dim][i] << endl;
         cout << results2ltb_h[dim][i] << "  "
              << results2lix_h[dim][i] << "  "
              << results2lmi_h[dim][i] << endl
              << results2lrg_h[dim][i] << "  "
              << results2lpk_h[dim][i] << endl;
         cout << resultsrtb_d[dim][i] << "  "
              << resultsrix_d[dim][i] << "  "
              << resultsrmi_d[dim][i] << endl
              << resultsrrg_d[dim][i] << "  "
              << resultsrpk_d[dim][i] << endl;
         cout << resultsltb_d[dim][i] << "  "
              << resultslix_d[dim][i] << "  "
              << resultslmi_d[dim][i] << endl
              << resultslrg_d[dim][i] << "  "
              << resultslpk_d[dim][i] << endl;
      }
      err = err + abs(results1rtb_h[dim][i] - results2rtb_h[dim][i])
                + abs(results1rix_h[dim][i] - results2rix_h[dim][i])
                + abs(results1rmi_h[dim][i] - results2rmi_h[dim][i])
                + abs(results1rrg_h[dim][i] - results2rrg_h[dim][i])
                + abs(results1rpk_h[dim][i] - results2rpk_h[dim][i])
                + abs(results1ltb_h[dim][i] - results2ltb_h[dim][i])
                + abs(results1lix_h[dim][i] - results2lix_h[dim][i])
                + abs(results1lmi_h[dim][i] - results2lmi_h[dim][i])
                + abs(results1lrg_h[dim][i] - results2lrg_h[dim][i])
                + abs(results1lpk_h[dim][i] - results2lpk_h[dim][i])
                + abs(results1rtb_h[dim][i] - resultsrtb_d[dim][i])
                + abs(results1rix_h[dim][i] - resultsrix_d[dim][i])
                + abs(results1rmi_h[dim][i] - resultsrmi_d[dim][i])
                + abs(results1rrg_h[dim][i] - resultsrrg_d[dim][i])
                + abs(results1rpk_h[dim][i] - resultsrpk_d[dim][i])
                + abs(results1ltb_h[dim][i] - resultsltb_d[dim][i])
                + abs(results1lix_h[dim][i] - resultslix_d[dim][i])
                + abs(results1lmi_h[dim][i] - resultslmi_d[dim][i])
                + abs(results1lrg_h[dim][i] - resultslrg_d[dim][i])
                + abs(results1lpk_h[dim][i] - resultslpk_d[dim][i]);
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
            cout << results1rtb_h[k][i] << "  "
                 << results1rix_h[k][i] << "  "
                 << results1rmi_h[k][i] << endl
                 << results1rrg_h[k][i] << "  "
                 << results1rpk_h[k][i] << endl;
            cout << results1ltb_h[k][i] << "  "
                 << results1lix_h[k][i] << "  "
                 << results1lmi_h[k][i] << endl
                 << results1lrg_h[k][i] << "  "
                 << results1lpk_h[k][i] << endl;
            cout << results2rtb_h[k][i] << "  "
                 << results2rix_h[k][i] << "  "
                 << results2rmi_h[k][i] << endl
                 << results2rrg_h[k][i] << "  "
                 << results2rpk_h[k][i] << endl;
            cout << results2ltb_h[k][i] << "  "
                 << results2lix_h[k][i] << "  "
                 << results2lmi_h[k][i] << endl
                 << results2lrg_h[k][i] << "  "
                 << results2lpk_h[k][i] << endl;
            cout << resultsrtb_d[k][i] << "  "
                 << resultsrix_d[k][i] << "  "
                 << resultsrmi_d[k][i] << endl
                 << resultsrrg_d[k][i] << "  "
                 << resultsrpk_d[k][i] << endl;
            cout << resultsltb_d[k][i] << "  "
                 << resultslix_d[k][i] << "  "
                 << resultslmi_d[k][i] << endl
                 << resultslrg_d[k][i] << "  "
                 << resultslpk_d[k][i] << endl;
         }
         err = err + abs(results1rtb_h[k][i] - results2rtb_h[k][i])
                   + abs(results1rix_h[k][i] - results2rix_h[k][i])
                   + abs(results1rmi_h[k][i] - results2rmi_h[k][i])
                   + abs(results1rrg_h[k][i] - results2rrg_h[k][i])
                   + abs(results1rpk_h[k][i] - results2rpk_h[k][i])
                   + abs(results1ltb_h[k][i] - results2ltb_h[k][i])
                   + abs(results1lix_h[k][i] - results2lix_h[k][i])
                   + abs(results1lmi_h[k][i] - results2lmi_h[k][i])
                   + abs(results1lrg_h[k][i] - results2lrg_h[k][i])
                   + abs(results1lpk_h[k][i] - results2lpk_h[k][i])
                   + abs(results1rtb_h[k][i] - resultsrtb_d[k][i])
                   + abs(results1rix_h[k][i] - resultsrix_d[k][i])
                   + abs(results1rmi_h[k][i] - resultsrmi_d[k][i])
                   + abs(results1rrg_h[k][i] - resultsrrg_d[k][i])
                   + abs(results1rpk_h[k][i] - resultsrpk_d[k][i])
                   + abs(results1ltb_h[k][i] - resultsltb_d[k][i])
                   + abs(results1lix_h[k][i] - resultslix_d[k][i])
                   + abs(results1lmi_h[k][i] - resultslmi_d[k][i])
                   + abs(results1lrg_h[k][i] - resultslrg_d[k][i])
                   + abs(results1lpk_h[k][i] - resultslpk_d[k][i]);
      }
      if(verbose) cout << "error : " << err << endl;
      sumerr = sumerr + err;
   }
   return sumerr;
}

double test_dbl10_real_polynomial
 ( int dim, int nbr, int nva, int pwr, int deg, int verbose, bool jobrep,
   int mode )
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
      int **idx = new int*[nbr];  // indices of variables in monomials
      int **exp = new int*[nbr];  // exponents of the variables

      bool vrb = (verbose > 1);

      dbl10_make_input(dim,nbr,nva,pwr,deg,nvr,idx,exp,
                       inputrtb,inputrix,inputrmi,inputrrg,inputrpk,
                       inputltb,inputlix,inputlmi,inputlrg,inputlpk,
                       cstrtb,cstrix,cstrmi,cstrrg,cstrpk,
                       cstltb,cstlix,cstlmi,cstlrg,cstlpk,
                       cffrtb,cffrix,cffrmi,cffrrg,cffrpk,
                       cffltb,cfflix,cfflmi,cfflrg,cfflpk,vrb);

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
         CPU_dbl10_poly_evaldiff
            (dim,nbr,deg,nvr,idx,
             cstrtb,cstrix,cstrmi,cstrrg,cstrpk,
             cstltb,cstlix,cstlmi,cstlrg,cstlpk,
             cffrtb,cffrix,cffrmi,cffrrg,cffrpk,
             cffltb,cfflix,cfflmi,cfflrg,cfflpk,
             inputrtb,inputrix,inputrmi,inputrrg,inputrpk,
             inputltb,inputlix,inputlmi,inputlrg,inputlpk,
             output1rtb_h,output1rix_h,output1rmi_h,output1rrg_h,output1rpk_h,
             output1ltb_h,output1lix_h,output1lmi_h,output1lrg_h,output1lpk_h,
             &timelapsec1_h,vrb);
         if(vrb) cout << "Computing with convolution jobs ..." << endl;
         CPU_dbl10_poly_evaldiffjobs
            (dim,nbr,deg,nvr,idx,
             cstrtb,cstrix,cstrmi,cstrrg,cstrpk,
             cstltb,cstlix,cstlmi,cstlrg,cstlpk,
             cffrtb,cffrix,cffrmi,cffrrg,cffrpk,
             cffltb,cfflix,cfflmi,cfflrg,cfflpk,
             inputrtb,inputrix,inputrmi,inputrrg,inputrpk,
             inputltb,inputlix,inputlmi,inputlrg,inputlpk,
             output2rtb_h,output2rix_h,output2rmi_h,output2rrg_h,output2rpk_h,
             output2ltb_h,output2lix_h,output2lmi_h,output2lrg_h,output2lpk_h,
             cnvjobs,addjobs,&timelapsec2_h,vrb);
      }
      if((mode == 0) || (mode == 2))
      {
         if(vrb) cout << "Computing on the device ..." << endl;
         GPU_dbl10_poly_evaldiff
            (deg+1,dim,nbr,deg,nvr,idx,
             cstrtb,cstrix,cstrmi,cstrrg,cstrpk,
             cstltb,cstlix,cstlmi,cstlrg,cstlpk,
             cffrtb,cffrix,cffrmi,cffrrg,cffrpk,
             cffltb,cfflix,cfflmi,cfflrg,cfflpk,
             inputrtb,inputrix,inputrmi,inputrrg,inputrpk,
             inputltb,inputlix,inputlmi,inputlrg,inputlpk,
             outputrtb_d,outputrix_d,outputrmi_d,outputrrg_d,outputrpk_d,
             outputltb_d,outputlix_d,outputlmi_d,outputlrg_d,outputlpk_d,
             cnvjobs,addjobs,&cnvlapms,&addlapms,&timelapms_d,
             &walltimes_d,vrb);
      }
      double sumerr = 0.0;
      if(mode == 2)
         sumerr = dbl10_error_sum(dim,deg,
                     output1rtb_h,output1rix_h,output1rmi_h,
                     output1rrg_h,output1rpk_h,
                     output1ltb_h,output1lix_h,output1lmi_h,
                     output1lrg_h,output1lpk_h,
                     output2rtb_h,output2rix_h,output2rmi_h,
                     output2rrg_h,output2rpk_h,
                     output2ltb_h,output2lix_h,output2lmi_h,
                     output2lrg_h,output2lpk_h,
                     outputrtb_d,outputrix_d,outputrmi_d,
                     outputrrg_d,outputrpk_d,
                     outputltb_d,outputlix_d,outputlmi_d,
                     outputlrg_d,outputlpk_d,vrb);

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

int test_dbl10_sequence
 ( int seed, int dim, int nva, int nbr, int pwr, int vrblvl,
   bool jobrep, int mode )
{
   const double tol = 1.0e-152;

   int deg = 0;
   cout << "---> running in deca double precision for degree 0 ..." << endl;
   int fail = main_dbl10_test_polynomial
                 (seed,dim,nbr,nva,pwr,deg,vrblvl,tol,jobrep,mode);
   deg = 8;
   cout << "---> running for degree 8 ..." << endl;
   fail += main_dbl10_test_polynomial
              (seed,dim,nbr,nva,pwr,deg,vrblvl,tol,false,mode);
   deg = 15;
   cout << "---> running for degree 15 ..." << endl;
   fail += main_dbl10_test_polynomial
              (seed,dim,nbr,nva,pwr,deg,vrblvl,tol,false,mode);
   deg = 31;
   cout << "---> running for degree 31 ..." << endl;
   fail += main_dbl10_test_polynomial
              (seed,dim,nbr,nva,pwr,deg,vrblvl,tol,false,mode);
   deg = 63;
   cout << "---> running for degree 63 ..." << endl;
   fail += main_dbl10_test_polynomial
              (seed,dim,nbr,nva,pwr,deg,vrblvl,tol,false,mode);
   deg = 95;
   cout << "---> running for degree 95 ..." << endl;
   fail += main_dbl10_test_polynomial
              (seed,dim,nbr,nva,pwr,deg,vrblvl,tol,false,mode);
   deg = 127;
   cout << "---> running for degree 127 ..." << endl;
   fail += main_dbl10_test_polynomial
              (seed,dim,nbr,nva,pwr,deg,vrblvl,tol,false,mode);
   deg = 152;
   cout << "---> running for degree 152 ..." << endl;
   fail += main_dbl10_test_polynomial
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
