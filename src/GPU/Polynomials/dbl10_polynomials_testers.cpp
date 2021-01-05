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

      if(verbose > 1)
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
      if(verbose > 1)
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
          cnvjobs,addjobs,&timelapms_d,vrb);

      double err = 0.0;

      if(vrb) cout << "The value of the polynomial :" << endl;
      for(int i=0; i<=deg; i++)
      {
         if(vrb)
         {
            cout << output1rtb_h[dim][i] << "  "
                 << output1rix_h[dim][i] << "  "
                 << output1rmi_h[dim][i] << endl
                 << output1rrg_h[dim][i] << "  "
                 << output1rpk_h[dim][i] << endl;
            cout << output1ltb_h[dim][i] << "  "
                 << output1lix_h[dim][i] << "  "
                 << output1lmi_h[dim][i] << endl
                 << output1lrg_h[dim][i] << "  "
                 << output1lpk_h[dim][i] << endl;
            cout << output2rtb_h[dim][i] << "  "
                 << output2rix_h[dim][i] << "  "
                 << output2rmi_h[dim][i] << endl
                 << output2rrg_h[dim][i] << "  "
                 << output2rpk_h[dim][i] << endl;
            cout << output2ltb_h[dim][i] << "  "
                 << output2lix_h[dim][i] << "  "
                 << output2lmi_h[dim][i] << endl
                 << output2lrg_h[dim][i] << "  "
                 << output2lpk_h[dim][i] << endl;
            cout << outputrtb_d[dim][i] << "  "
                 << outputrix_d[dim][i] << "  "
                 << outputrmi_d[dim][i] << endl
                 << outputrrg_d[dim][i] << "  "
                 << outputrpk_d[dim][i] << endl;
            cout << outputltb_d[dim][i] << "  "
                 << outputlix_d[dim][i] << "  "
                 << outputlmi_d[dim][i] << endl
                 << outputlrg_d[dim][i] << "  "
                 << outputlpk_d[dim][i] << endl;
         }
         err = err + abs(output1rtb_h[dim][i] - output2rtb_h[dim][i])
                   + abs(output1rix_h[dim][i] - output2rix_h[dim][i])
                   + abs(output1rmi_h[dim][i] - output2rmi_h[dim][i])
                   + abs(output1rrg_h[dim][i] - output2rrg_h[dim][i])
                   + abs(output1rpk_h[dim][i] - output2rpk_h[dim][i])
                   + abs(output1ltb_h[dim][i] - output2ltb_h[dim][i])
                   + abs(output1lix_h[dim][i] - output2lix_h[dim][i])
                   + abs(output1lmi_h[dim][i] - output2lmi_h[dim][i])
                   + abs(output1lrg_h[dim][i] - output2lrg_h[dim][i])
                   + abs(output1lpk_h[dim][i] - output2lpk_h[dim][i])
                   + abs(output1rtb_h[dim][i] - outputrtb_d[dim][i])
                   + abs(output1rix_h[dim][i] - outputrix_d[dim][i])
                   + abs(output1rmi_h[dim][i] - outputrmi_d[dim][i])
                   + abs(output1rrg_h[dim][i] - outputrrg_d[dim][i])
                   + abs(output1rpk_h[dim][i] - outputrpk_d[dim][i])
                   + abs(output1ltb_h[dim][i] - outputltb_d[dim][i])
                   + abs(output1lix_h[dim][i] - outputlix_d[dim][i])
                   + abs(output1lmi_h[dim][i] - outputlmi_d[dim][i])
                   + abs(output1lrg_h[dim][i] - outputlrg_d[dim][i])
                   + abs(output1lpk_h[dim][i] - outputlpk_d[dim][i]);
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
               cout << output1rtb_h[k][i] << "  "
                    << output1rix_h[k][i] << "  "
                    << output1rmi_h[k][i] << endl
                    << output1rrg_h[k][i] << "  "
                    << output1rpk_h[k][i] << endl;
               cout << output1ltb_h[k][i] << "  "
                    << output1lix_h[k][i] << "  "
                    << output1lmi_h[k][i] << endl
                    << output1lrg_h[k][i] << "  "
                    << output1lpk_h[k][i] << endl;
               cout << output2rtb_h[k][i] << "  "
                    << output2rix_h[k][i] << "  "
                    << output2rmi_h[k][i] << endl
                    << output2rrg_h[k][i] << "  "
                    << output2rpk_h[k][i] << endl;
               cout << output2ltb_h[k][i] << "  "
                    << output2lix_h[k][i] << "  "
                    << output2lmi_h[k][i] << endl
                    << output2lrg_h[k][i] << "  "
                    << output2lpk_h[k][i] << endl;
               cout << outputrtb_d[k][i] << "  "
                    << outputrix_d[k][i] << "  "
                    << outputrmi_d[k][i] << endl
                    << outputrrg_d[k][i] << "  "
                    << outputrpk_d[k][i] << endl;
               cout << outputltb_d[k][i] << "  "
                    << outputlix_d[k][i] << "  "
                    << outputlmi_d[k][i] << endl
                    << outputlrg_d[k][i] << "  "
                    << outputlpk_d[k][i] << endl;
            }
            err = err + abs(output1rtb_h[k][i] - output2rtb_h[k][i])
                      + abs(output1rix_h[k][i] - output2rix_h[k][i])
                      + abs(output1rmi_h[k][i] - output2rmi_h[k][i])
                      + abs(output1rrg_h[k][i] - output2rrg_h[k][i])
                      + abs(output1rpk_h[k][i] - output2rpk_h[k][i])
                      + abs(output1ltb_h[k][i] - output2ltb_h[k][i])
                      + abs(output1lix_h[k][i] - output2lix_h[k][i])
                      + abs(output1lmi_h[k][i] - output2lmi_h[k][i])
                      + abs(output1lrg_h[k][i] - output2lrg_h[k][i])
                      + abs(output1lpk_h[k][i] - output2lpk_h[k][i])
                      + abs(output1rtb_h[k][i] - outputrtb_d[k][i])
                      + abs(output1rix_h[k][i] - outputrix_d[k][i])
                      + abs(output1rmi_h[k][i] - outputrmi_d[k][i])
                      + abs(output1rrg_h[k][i] - outputrrg_d[k][i])
                      + abs(output1rpk_h[k][i] - outputrpk_d[k][i])
                      + abs(output1ltb_h[k][i] - outputltb_d[k][i])
                      + abs(output1lix_h[k][i] - outputlix_d[k][i])
                      + abs(output1lmi_h[k][i] - outputlmi_d[k][i])
                      + abs(output1lrg_h[k][i] - outputlrg_d[k][i])
                      + abs(output1lpk_h[k][i] - outputlpk_d[k][i]);
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
         cout << "Time spent by all kernels in milliseconds : ";
         cout << fixed << setprecision(2) << timelapms_d << endl;
         cout << scientific << setprecision(16);
      }
      return sumerr;
   }
}
