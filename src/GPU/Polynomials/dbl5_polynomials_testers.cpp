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
#include "job_makers.h"
#include "dbl5_polynomials_testers.h"

using namespace std;

int main_dbl5_test_polynomial
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

   double realsum = test_dbl5_real_polynomial
                       (dim,nbr,nva,pwr,deg,vrblvl-1,jobrep,mode);

   int fail = int(realsum > tol);

   if(vrblvl > 0)
   {
      if(mode == 2)
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

void cmplx5_make_input
 ( int dim, int nbr, int nva, int pwr, int deg,
   int *nvr, int **idx, int **exp,
   double **inputretb, double **inputreix, double **inputremi,
   double **inputrerg, double **inputrepk,
   double **inputimtb, double **inputimix, double **inputimmi,
   double **inputimrg, double **inputimpk,
   double *cstretb, double *cstreix, double *cstremi,
   double *cstrerg, double *cstrepk,
   double *cstimtb, double *cstimix, double *cstimmi,
   double *cstimrg, double *cstimpk,
   double **cffretb, double **cffreix, double **cffremi,
   double **cffrerg, double **cffrepk,
   double **cffimtb, double **cffimix, double **cffimmi,
   double **cffimrg, double **cffimpk, bool verbose )
{
   make_complex5_input(dim,deg,
      inputretb,inputreix,inputremi,inputrerg,inputrepk,
      inputimtb,inputimix,inputimmi,inputimrg,inputimpk);

   if(verbose)
   {
      cout << scientific << setprecision(16);
      cout << "Random input series :" << endl;
      for(int i=0; i<dim; i++)
      {
         cout << "-> coefficients of series " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
         {
            cout << inputretb[i][j] << "  " << inputreix[i][j]
                                    << "  " << inputremi[i][j] << endl
                 << inputrerg[i][j] << "  " << inputrepk[i][j] << endl;
            cout << inputimtb[i][j] << "  " << inputimix[i][j]
                                    << "  " << inputimmi[i][j] << endl
                 << inputimrg[i][j] << "  " << inputimpk[i][j] << endl;
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
         make_complex5_cyclic
            (dim,nva,deg,idx,
             cstretb,cstreix,cstremi,cstrerg,cstrepk,
             cstimtb,cstimix,cstimmi,cstimrg,cstimpk,
             cffretb,cffreix,cffremi,cffrerg,cffrepk,
             cffimtb,cffimix,cffimmi,cffimrg,cffimpk);
      else
         make_complex5_products
            (dim,nbr,nva,deg,idx,
             cstretb,cstreix,cstremi,cstrerg,cstrepk,
             cstimtb,cstimix,cstimmi,cstimrg,cstimpk,
             cffretb,cffreix,cffremi,cffrerg,cffrepk,
             cffimtb,cffimix,cffimmi,cffimrg,cffimpk);
   }
   else
   {
      for(int i=0; i<nbr; i++) exp[i] = new int[nvr[i]];

      bool fail = make_complex5_polynomial
                     (dim,nbr,pwr,deg,nvr,idx,exp,
                      cstretb,cstreix,cstremi,cstrerg,cstrepk,
                      cstimtb,cstimix,cstimmi,cstimrg,cstimpk,
                      cffretb,cffreix,cffremi,cffrerg,cffrepk,
                      cffimtb,cffimix,cffimmi,cffimrg,cffimpk);
   }
   if(verbose)
   {
      cout << "Coefficient series of the constant term :" << endl;
      for(int j=0; j<=deg; j++)
      {
         cout << cstretb[j] << "  " << cstreix[j]
                            << "  " << cstremi[j] << endl
              << cstrerg[j] << "  " << cstrepk[j] << endl;
         cout << cstimtb[j] << "  " << cstimix[j]
                            << "  " << cstimix[j] << endl
              << cstimrg[j] << "  " << cstimpk[j] << endl;
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
            cout << cffretb[i][j] << "  " << cffreix[i][j]
                                  << "  " << cffremi[i][j] << endl
                 << cffrerg[i][j] << "  " << cffrepk[i][j] << endl;
            cout << cffimtb[i][j] << "  " << cffimix[i][j]
                                  << "  " << cffimmi[i][j] << endl
                 << cffimrg[i][j] << "  " << cffimpk[i][j] << endl;
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
 ( int dim, int nbr, int nva, int pwr, int deg, int verbose, bool jobrep,
   int mode )
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
      double timelapsec1_h,timelapsec2_h;
      double cnvlapms,addlapms,timelapms_d,walltimes_d;

      if((mode == 1) || (mode == 2))
      {
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
      }
      if((mode == 0) || (mode == 2))
      {
         if(vrb) cout << "Computing on the device ..." << endl;
         GPU_dbl5_poly_evaldiff
            (deg+1,dim,nbr,deg,nvr,idx,csttb,cstix,cstmi,cstrg,cstpk,
             cfftb,cffix,cffmi,cffrg,cffpk,
             inputtb,inputix,inputmi,inputrg,inputpk,
             outputtb_d,outputix_d,outputmi_d,outputrg_d,outputpk_d,
             cnvjobs,addjobs,&cnvlapms,&addlapms,&timelapms_d,
             &walltimes_d,vrb);
      }
      double sumerr = 0.0;
      if(mode == 2)
         sumerr = dbl5_error_sum(dim,deg,
                     output1tb_h,output1ix_h,output1mi_h,
                     output1rg_h,output1pk_h,
                     output2tb_h,output2ix_h,output2mi_h,
                     output2rg_h,output2pk_h,
                     outputtb_d,outputix_d,outputmi_d,
                     outputrg_d,outputpk_d,vrb);

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

double test_cmplx5_real_polynomial
 ( int dim, int nbr, int nva, int pwr, int deg, int verbose, bool jobrep,
   int mode )
{
   if(nbr < 1)
      return 0.0;
   else
   {
      double **inputretb = new double*[dim]; // dim series of degree deg
      double **inputreix = new double*[dim];
      double **inputremi = new double*[dim];
      double **inputrerg = new double*[dim];
      double **inputrepk = new double*[dim];
      double **inputimtb = new double*[dim];
      double **inputimix = new double*[dim];
      double **inputimmi = new double*[dim];
      double **inputimrg = new double*[dim];
      double **inputimpk = new double*[dim];
      for(int i=0; i<dim; i++)
      {
         inputretb[i] = new double[deg+1];
         inputreix[i] = new double[deg+1];
         inputremi[i] = new double[deg+1];
         inputrerg[i] = new double[deg+1];
         inputrepk[i] = new double[deg+1];
         inputimtb[i] = new double[deg+1];
         inputimix[i] = new double[deg+1];
         inputimmi[i] = new double[deg+1];
         inputimrg[i] = new double[deg+1];
         inputimpk[i] = new double[deg+1];
      }
      // The output are dim+1 power series of degree deg
      // for the evaluated and differentiated polynomial.
      double **output1retb_h = new double*[dim+1];
      double **output1reix_h = new double*[dim+1];
      double **output1remi_h = new double*[dim+1];
      double **output1rerg_h = new double*[dim+1];
      double **output1repk_h = new double*[dim+1];
      double **output1imtb_h = new double*[dim+1];
      double **output1imix_h = new double*[dim+1];
      double **output1immi_h = new double*[dim+1];
      double **output1imrg_h = new double*[dim+1];
      double **output1impk_h = new double*[dim+1];
      for(int i=0; i<=dim; i++)
      {
         output1retb_h[i] = new double[deg+1];
         output1reix_h[i] = new double[deg+1];
         output1remi_h[i] = new double[deg+1];
         output1rerg_h[i] = new double[deg+1];
         output1repk_h[i] = new double[deg+1];
         output1imtb_h[i] = new double[deg+1];
         output1imix_h[i] = new double[deg+1];
         output1immi_h[i] = new double[deg+1];
         output1imrg_h[i] = new double[deg+1];
         output1impk_h[i] = new double[deg+1];
      }
      double **output2retb_h = new double*[dim+1];
      double **output2reix_h = new double*[dim+1];
      double **output2remi_h = new double*[dim+1];
      double **output2rerg_h = new double*[dim+1];
      double **output2repk_h = new double*[dim+1];
      double **output2imtb_h = new double*[dim+1];
      double **output2imix_h = new double*[dim+1];
      double **output2immi_h = new double*[dim+1];
      double **output2imrg_h = new double*[dim+1];
      double **output2impk_h = new double*[dim+1];
      for(int i=0; i<=dim; i++)
      {
         output2retb_h[i] = new double[deg+1];
         output2reix_h[i] = new double[deg+1];
         output2remi_h[i] = new double[deg+1];
         output2rerg_h[i] = new double[deg+1];
         output2repk_h[i] = new double[deg+1];
         output2imtb_h[i] = new double[deg+1];
         output2imix_h[i] = new double[deg+1];
         output2immi_h[i] = new double[deg+1];
         output2imrg_h[i] = new double[deg+1];
         output2impk_h[i] = new double[deg+1];
      }
      double **outputretb_d = new double*[dim+1];
      double **outputreix_d = new double*[dim+1];
      double **outputremi_d = new double*[dim+1];
      double **outputrerg_d = new double*[dim+1];
      double **outputrepk_d = new double*[dim+1];
      double **outputimtb_d = new double*[dim+1];
      double **outputimix_d = new double*[dim+1];
      double **outputimmi_d = new double*[dim+1];
      double **outputimrg_d = new double*[dim+1];
      double **outputimpk_d = new double*[dim+1];
      for(int i=0; i<=dim; i++)
      {
         outputretb_d[i] = new double[deg+1];
         outputreix_d[i] = new double[deg+1];
         outputremi_d[i] = new double[deg+1];
         outputrerg_d[i] = new double[deg+1];
         outputrepk_d[i] = new double[deg+1];
         outputimtb_d[i] = new double[deg+1];
         outputimix_d[i] = new double[deg+1];
         outputimmi_d[i] = new double[deg+1];
         outputimrg_d[i] = new double[deg+1];
         outputimpk_d[i] = new double[deg+1];
      }
      double *cstretb = new double[deg+1]; // constant coefficient series
      double *cstreix = new double[deg+1];
      double *cstremi = new double[deg+1];
      double *cstrerg = new double[deg+1];
      double *cstrepk = new double[deg+1];
      double *cstimtb = new double[deg+1];
      double *cstimix = new double[deg+1];
      double *cstimmi = new double[deg+1];
      double *cstimrg = new double[deg+1];
      double *cstimpk = new double[deg+1];
      double **cffretb = new double*[nbr]; // coefficient series of terms
      double **cffreix = new double*[nbr];
      double **cffremi = new double*[nbr];
      double **cffrerg = new double*[nbr];
      double **cffrepk = new double*[nbr];
      double **cffimtb = new double*[nbr];
      double **cffimix = new double*[nbr];
      double **cffimmi = new double*[nbr];
      double **cffimrg = new double*[nbr];
      double **cffimpk = new double*[nbr];
      for(int i=0; i<nbr; i++)
      {
         cffretb[i] = new double[deg+1];
         cffreix[i] = new double[deg+1];
         cffremi[i] = new double[deg+1];
         cffrerg[i] = new double[deg+1];
         cffrepk[i] = new double[deg+1];
         cffimtb[i] = new double[deg+1];
         cffimix[i] = new double[deg+1];
         cffimmi[i] = new double[deg+1];
         cffimrg[i] = new double[deg+1];
         cffimpk[i] = new double[deg+1];
      }
      int *nvr = new int[nbr]; // number of variables in each monomial
      int **idx = new int*[nbr];  // indices of variables in monomials
      int **exp = new int*[nbr];  // exponents of the variables

      bool vrb = (verbose > 1);

      cmplx5_make_input(dim,nbr,nva,pwr,deg,nvr,idx,exp,
         inputretb,inputreix,inputremi,inputrerg,inputrepk,
         inputimtb,inputimix,inputimmi,inputimrg,inputimpk,
         cstretb,cstreix,cstremi,cstrerg,cstrepk,
         cstimtb,cstimix,cstimmi,cstimrg,cstimpk,
         cffretb,cffreix,cffremi,cffrerg,cffrepk,
         cffimtb,cffimix,cffimmi,cffimrg,cffimpk,vrb);

      ConvolutionJobs cnvjobs(dim);
      AdditionJobs addjobs(dim,nbr);

      make_all_jobs(dim,nbr,nvr,idx,&cnvjobs,&addjobs,vrb);

      double timelapsec1_h,timelapsec2_h;
      double cnvlapms,addlapms,timelapms_d,walltimes_d;

      if((mode == 1) || (mode == 2))
      {
         if(vrb) cout << "Computing without convolution jobs ..." << endl;
         CPU_cmplx5_poly_evaldiff
            (dim,nbr,deg,nvr,idx,
             cstretb,cstreix,cstremi,cstrerg,cstrepk,
             cstimtb,cstimix,cstimmi,cstimrg,cstimpk,
             cffretb,cffreix,cffremi,cffrerg,cffrepk,
             cffimtb,cffimix,cffimmi,cffimrg,cffimpk,
             inputretb,inputreix,inputremi,inputrerg,inputrepk,
             inputimtb,inputimix,inputimmi,inputimrg,inputimpk,
             output1retb_h,output1reix_h,output1remi_h,
             output1rerg_h,output1repk_h,
             output1imtb_h,output1imix_h,output1imrg_h,
             output1imrg_h,output1impk_h,&timelapsec1_h,vrb);
   /*
         if(vrb) cout << "Computing with convolution jobs ..." << endl;
         CPU_cmplx5_poly_evaldiffjobs
            (dim,nbr,deg,nvr,idx,
             cstretb,cstreix,cstremi,cstrerg,cstrepk,
             cstimtb,cstimix,cstimmi,cstimrg,cstimpk,
             cffretb,cffreix,cffremi,cffrerg,cffrepk,
             cffimtb,cffimix,cffimmi,cffimrg,cffimpk,
             inputretb,inputreix,inputremi,inputrerg,inputrepk,
             inputimtb,inputimix,inputimmi,inputimrg,inputimpk,
             output2retb_h,output2reix_h,output2remi_h,
             output2rerg_h,output2repk_h,
             output2imtb_h,output2imix_h,output2immi_h,
             output2imrg_h,output2impk_h,
             cnvjobs,addjobs,&timelapsec2_h,vrb);
    */
      }
    /*
      if((mode == 0) || (mode == 2))
      {
         if(vrb) cout << "Computing on the device ..." << endl;
         GPU_cmplx5_poly_evaldiff
            (deg+1,dim,nbr,deg,nvr,idx,
             cstretb,cstreix,cstremi,cstrerg,cstrepk,
             cstimtb,cstimix,cstimmi,cstimrg,cstimpk,
             cffretb,cffreix,cffremi,cffrerg,cffrepk,
             cffimtb,cffimix,cffimmi,cffimrg,cffimpk,
             inputretb,inputreix,inputremi,inputrerg,inputrepk,
             inputimtb,inputimix,inputimmi,inputimrg,inputimpk,
             outputretb_d,outputreix_d,outputremi_d,outputrerg_d,outputrepk_d,
             outputimtb_d,outputimix_d,outputimmi_d,outputimrg_d,outputimpk_d,
             cnvjobs,addjobs,&cnvlapms,&addlapms,&timelapms_d,
             &walltimes_d,vrb);
      }
     */
      double sumerr = 0.0;
    /*
      if(mode == 2)
         sumerr = cmplx5_error_sum(dim,deg,
                     output1retb_h,output1reix_h,output1remi_h,
                     output1rerg_h,output1repk_h,
                     output1imtb_h,output1imix_h,output1immi_h,
                     output1imrg_h,output1impk_h,
                     output2retb_h,output2reix_h,output2remi_h,
                     output2rerg_h,output2repk_h,
                     output2imtb_h,output2imix_h,output2immi_h,
                     output2imrg_h,output2impk_h,
                     outputretb_d,outputreix_d,outputremi_d,
                     outputrerg_d,outputrepk_d,
                     outputimtb_d,outputimix_d,outputimmi_d,
                     outputimrg_d,outputimpk_d,vrb);

      if(verbose > 0)
      {
         if(jobrep) write_jobs_report(dim,nva,nbr,deg,cnvjobs,addjobs);
         if((mode == 1) || (mode == 2))
            write_CPU_timings(timelapsec1_h,timelapsec2_h);
         if((mode == 0) || (mode == 2))
            write_GPU_timings(cnvlapms,addlapms,timelapms_d,walltimes_d);
      } 
      */
      return sumerr;
   }
}

int test_dbl5_sequence
 ( int seed, int dim, int nva, int nbr, int pwr, int vrblvl,
   bool jobrep, int mode )
{
   const double tol = 1.0e-72;

   int deg = 0;
   cout << "---> running in penta double precision for degree 0 ..." << endl;
   int fail = main_dbl5_test_polynomial
                (seed,dim,nbr,nva,pwr,deg,vrblvl,tol,jobrep,mode);
   deg = 8;
   cout << "---> running for degree 8 ..." << endl;
   fail += main_dbl5_test_polynomial
              (seed,dim,nbr,nva,pwr,deg,vrblvl,tol,false,mode);
   deg = 15;
   cout << "---> running for degree 15 ..." << endl;
   fail += main_dbl5_test_polynomial
              (seed,dim,nbr,nva,pwr,deg,vrblvl,tol,false,mode);
   deg = 31;
   cout << "---> running for degree 31 ..." << endl;
   fail += main_dbl5_test_polynomial
              (seed,dim,nbr,nva,pwr,deg,vrblvl,tol,false,mode);
   deg = 63;
   cout << "---> running for degree 63 ..." << endl;
   fail += main_dbl5_test_polynomial
              (seed,dim,nbr,nva,pwr,deg,vrblvl,tol,false,mode);
   deg = 95;
   cout << "---> running for degree 95 ..." << endl;
   fail += main_dbl5_test_polynomial
              (seed,dim,nbr,nva,pwr,deg,vrblvl,tol,false,mode);
   deg = 127;
   cout << "---> running for degree 127 ..." << endl;
   fail += main_dbl5_test_polynomial
              (seed,dim,nbr,nva,pwr,deg,vrblvl,tol,false,mode);
   deg = 152;
   cout << "---> running for degree 152 ..." << endl;
   fail += main_dbl5_test_polynomial
              (seed,dim,nbr,nva,pwr,deg,vrblvl,tol,false,mode);
   deg = 159;
   cout << "---> running for degree 159 ..." << endl;
   fail += main_dbl5_test_polynomial
              (seed,dim,nbr,nva,pwr,deg,vrblvl,tol,false,mode);
   deg = 191;
   cout << "---> running for degree 191 ..." << endl;
   fail += main_dbl5_test_polynomial
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
