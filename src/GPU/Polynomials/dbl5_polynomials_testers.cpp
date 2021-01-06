/* The file dbl5_polynomials_testers.cpp contains the definitions of
 * functions with prototypes in dbl5_polynomials_testers.h. */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <vector_types.h>
#include "random_polynomials.h"
#include "random5_monomials.h"
#include "random5_polynomials.h"
#include "convolution_jobs.h"
#include "addition_jobs.h"
#include "write_job_counts.h"
#include "dbl5_polynomials_host.h"
#include "dbl5_polynomials_kernels.h"
#include "dbl5_polynomials_testers.h"

using namespace std;

int main_dbl5_test_polynomial
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

   double realsum = test_dbl5_real_polynomial(dim,nbr,nva,pwr,deg,vrblvl-1);

   int fail = int(realsum > tol);

   if(vrblvl > 0)
   {
      cout << scientific << setprecision(2);
      cout << "Sum of all errors in penta double precision :" << endl;
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

void dbl5_make_input
 ( int dim, int nbr, int nva, int pwr, int deg,
   int *nvr, int **idx, int **exp,
   double **inputtb, double **inputix, double **inputmi,
   double **inputrg, double **inputpk,
   double *csttb, double *cstix, double *cstmi, 
   double *cstrg, double *cstpk,
   double **cfftb, double **cffix, double **cffmi, 
   double **cffrg, double **cffpk, bool verbose )
{
   make_real5_input(dim,deg,inputtb,inputix,inputmi,inputrg,inputpk);

   if(verbose)
   {
      cout << scientific << setprecision(16);
      cout << "Random input series :" << endl;
      for(int i=0; i<dim; i++)
      {
         cout << "-> coefficients of series " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
            cout << inputtb[i][j] << "  " << inputix[i][j]
                                  << "  " << inputmi[i][j] << endl
                 << inputrg[i][j] << "  " << inputpk[i][j] << endl;
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
         make_real5_cyclic
            (dim,nva,deg,idx,csttb,cstix,cstmi,cstrg,cstpk,
                             cfftb,cffix,cffmi,cffrg,cffpk);
      else
         make_real5_products
            (dim,nbr,nva,deg,idx,csttb,cstix,cstmi,cstrg,cstpk,
                                 cfftb,cffix,cffmi,cffrg,cffpk);
   }
   else
   {
      for(int i=0; i<nbr; i++) exp[i] = new int[nvr[i]];

      bool fail = make_real5_polynomial
                     (dim,nbr,pwr,deg,nvr,idx,exp,
                      csttb,cstix,cstmi,cstrg,cstpk,
                      cfftb,cffix,cffmi,cffrg,cffpk);
   }
   if(verbose)
   {
      cout << "Coefficient series of the constant term :" << endl;
      for(int j=0; j<=deg; j++)
         cout << csttb[j] << "  " << cstix[j] << "  " << cstmi[j] << endl
              << cstrg[j] << "  " << cstpk[j] << endl;

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
            cout << cfftb[i][j] << "  " << cffix[i][j]
                                << "  " << cffmi[i][j] << endl
                 << cffrg[i][j] << "  " << cffpk[i][j] << endl;
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

double dbl5_error_sum
 ( int dim, int deg,
   double **results1tb_h, double **results1ix_h, double **results1mi_h, 
   double **results1rg_h, double **results1pk_h,
   double **results2tb_h, double **results2ix_h, double **results2mi_h, 
   double **results2rg_h, double **results2pk_h,
   double **resultstb_d, double **resultsix_d, double **resultsmi_d, 
   double **resultsrg_d, double **resultspk_d, bool verbose )
{
   double err = 0.0;

   if(verbose) cout << "The value of the polynomial :" << endl;
   for(int i=0; i<=deg; i++)
   {
      if(verbose)
      {
         cout << results1tb_h[dim][i] << "  "
              << results1ix_h[dim][i] << "  "
              << results1mi_h[dim][i] << endl
              << results1rg_h[dim][i] << "  "
              << results1pk_h[dim][i] << endl;
         cout << results2tb_h[dim][i] << "  "
              << results2ix_h[dim][i] << "  "
              << results2mi_h[dim][i] << endl
              << results2rg_h[dim][i] << "  "
              << results2pk_h[dim][i] << endl;
         cout << resultstb_d[dim][i] << "  "
              << resultsix_d[dim][i] << "  "
              << resultsmi_d[dim][i] << endl
              << resultsrg_d[dim][i] << "  "
              << resultspk_d[dim][i] << endl;
      }
      err = err + abs(results1tb_h[dim][i] - results2tb_h[dim][i])
                + abs(results1ix_h[dim][i] - results2ix_h[dim][i])
                + abs(results1mi_h[dim][i] - results2mi_h[dim][i])
                + abs(results1rg_h[dim][i] - results2rg_h[dim][i])
                + abs(results1pk_h[dim][i] - results2pk_h[dim][i])
                + abs(results1tb_h[dim][i] - resultstb_d[dim][i])
                + abs(results1ix_h[dim][i] - resultsix_d[dim][i])
                + abs(results1mi_h[dim][i] - resultsmi_d[dim][i])
                + abs(results1rg_h[dim][i] - resultsrg_d[dim][i])
                + abs(results1pk_h[dim][i] - resultspk_d[dim][i]);
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
            cout << results1tb_h[k][i] << "  "
                 << results1ix_h[k][i] << "  "
                 << results1mi_h[k][i] << endl
                 << results1rg_h[k][i] << "  "
                 << results1pk_h[k][i] << endl;
            cout << results2tb_h[k][i] << "  "
                 << results2ix_h[k][i] << "  "
                 << results2mi_h[k][i] << endl
                 << results2rg_h[k][i] << "  "
                 << results2pk_h[k][i] << endl;
            cout << resultstb_d[k][i] << "  "
                 << resultsix_d[k][i] << "  "
                 << resultsmi_d[k][i] << endl
                 << resultsrg_d[k][i] << "  "
                 << resultspk_d[k][i] << endl;
         }
         err = err + abs(results1tb_h[k][i] - results2tb_h[k][i])
                   + abs(results1ix_h[k][i] - results2ix_h[k][i])
                   + abs(results1mi_h[k][i] - results2mi_h[k][i])
                   + abs(results1rg_h[k][i] - results2rg_h[k][i])
                   + abs(results1pk_h[k][i] - results2pk_h[k][i])
                   + abs(results1tb_h[k][i] - resultstb_d[k][i])
                   + abs(results1ix_h[k][i] - resultsix_d[k][i])
                   + abs(results1mi_h[k][i] - resultsmi_d[k][i])
                   + abs(results1rg_h[k][i] - resultsrg_d[k][i])
                   + abs(results1pk_h[k][i] - resultspk_d[k][i]);
      }
      if(verbose) cout << "error : " << err << endl;
      sumerr = sumerr + err;
   }
   return sumerr;
}

double test_dbl5_real_polynomial
 ( int dim, int nbr, int nva, int pwr, int deg, int verbose )
{
   if(nbr < 1)
      return 0.0;
   else
   {
      double **inputtb = new double*[dim]; // dim series of degree deg
      for(int i=0; i<dim; i++) inputtb[i] = new double[deg+1];
      double **inputix = new double*[dim];
      for(int i=0; i<dim; i++) inputix[i] = new double[deg+1];
      double **inputmi = new double*[dim];
      for(int i=0; i<dim; i++) inputmi[i] = new double[deg+1];
      double **inputrg = new double*[dim];
      for(int i=0; i<dim; i++) inputrg[i] = new double[deg+1];
      double **inputpk = new double*[dim];
      for(int i=0; i<dim; i++) inputpk[i] = new double[deg+1];
      // The output are dim+1 power series of degree deg
      // for the evaluated and differentiated polynomial.
      double **output1tb_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output1tb_h[i] = new double[deg+1];
      double **output1ix_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output1ix_h[i] = new double[deg+1];
      double **output1mi_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output1mi_h[i] = new double[deg+1];
      double **output1rg_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output1rg_h[i] = new double[deg+1];
      double **output1pk_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output1pk_h[i] = new double[deg+1];
      double **output2tb_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output2tb_h[i] = new double[deg+1];
      double **output2ix_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output2ix_h[i] = new double[deg+1];
      double **output2mi_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output2mi_h[i] = new double[deg+1];
      double **output2rg_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output2rg_h[i] = new double[deg+1];
      double **output2pk_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output2pk_h[i] = new double[deg+1];
      double **outputtb_d = new double*[dim+1];
      for(int i=0; i<=dim; i++) outputtb_d[i] = new double[deg+1];
      double **outputix_d = new double*[dim+1];
      for(int i=0; i<=dim; i++) outputix_d[i] = new double[deg+1];
      double **outputmi_d = new double*[dim+1];
      for(int i=0; i<=dim; i++) outputmi_d[i] = new double[deg+1];
      double **outputrg_d = new double*[dim+1];
      for(int i=0; i<=dim; i++) outputrg_d[i] = new double[deg+1];
      double **outputpk_d = new double*[dim+1];
      for(int i=0; i<=dim; i++) outputpk_d[i] = new double[deg+1];

      double *csttb = new double[deg+1]; // constant coefficient series
      double *cstix = new double[deg+1];
      double *cstmi = new double[deg+1];
      double *cstrg = new double[deg+1];
      double *cstpk = new double[deg+1];
      double **cfftb = new double*[nbr]; // coefficient series of terms
      for(int i=0; i<nbr; i++) cfftb[i] = new double[deg+1];
      double **cffix = new double*[nbr];
      for(int i=0; i<nbr; i++) cffix[i] = new double[deg+1];
      double **cffmi = new double*[nbr];
      for(int i=0; i<nbr; i++) cffmi[i] = new double[deg+1];
      double **cffrg = new double*[nbr];
      for(int i=0; i<nbr; i++) cffrg[i] = new double[deg+1];
      double **cffpk = new double*[nbr];
      for(int i=0; i<nbr; i++) cffpk[i] = new double[deg+1];
      int *nvr = new int[nbr]; // number of variables in each monomial
      int **idx = new int*[nbr];  // indices of variables in monomials
      int **exp = new int*[nbr];  // exponents of the variables

      bool vrb = (verbose > 1);

      dbl5_make_input(dim,nbr,nva,pwr,deg,nvr,idx,exp,
                      inputtb,inputix,inputmi,inputrg,inputpk,
                      csttb,cstix,cstmi,cstrg,cstpk,
                      cfftb,cffix,cffmi,cffrg,cffpk,vrb);

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
      CPU_dbl5_poly_evaldiff
         (dim,nbr,deg,nvr,idx,csttb,cstix,cstmi,cstrg,cstpk,
          cfftb,cffix,cffmi,cffrg,cffpk,
          inputtb,inputix,inputmi,inputrg,inputpk,
          output1tb_h,output1ix_h,output1mi_h,output1rg_h,output1pk_h,
          &timelapsec1_h,vrb);
      if(vrb) cout << "Computing with convolution jobs ..." << endl;
      CPU_dbl5_poly_evaldiffjobs
         (dim,nbr,deg,nvr,idx,csttb,cstix,cstmi,cstrg,cstpk,
          cfftb,cffix,cffmi,cffrg,cffpk,
          inputtb,inputix,inputmi,inputrg,inputpk,
          output2tb_h,output2ix_h,output2mi_h,output2rg_h,output2pk_h,
          cnvjobs,addjobs,&timelapsec2_h,vrb);
      if(vrb) cout << "Computing on the device ..." << endl;
      GPU_dbl5_poly_evaldiff
         (deg+1,dim,nbr,deg,nvr,idx,csttb,cstix,cstmi,cstrg,cstpk,
          cfftb,cffix,cffmi,cffrg,cffpk,
          inputtb,inputix,inputmi,inputrg,inputpk,
          outputtb_d,outputix_d,outputmi_d,outputrg_d,outputpk_d,
          cnvjobs,addjobs,&timelapms_d,vrb);

      double sumerr = dbl5_error_sum(dim,deg,
                         output1tb_h,output1ix_h,output1mi_h,
                         output1rg_h,output1pk_h,
                         output2tb_h,output2ix_h,output2mi_h,
                         output2rg_h,output2pk_h,
                         outputtb_d,outputix_d,outputmi_d,
                         outputrg_d,outputpk_d,vrb);

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
