/* The file dbl4_polynomials_testers.cpp contains the definitions of
 * functions with prototypes in dbl4_polynomials_testers.h. */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <vector_types.h>
#include "random_polynomials.h"
#include "random4_monomials.h"
#include "random4_polynomials.h"
#include "convolution_jobs.h"
#include "addition_jobs.h"
#include "write_job_counts.h"
#include "dbl4_polynomials_host.h"
#include "dbl4_polynomials_kernels.h"
#include "dbl4_polynomials_testers.h"

using namespace std;

int main_dbl4_test_polynomial
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

   double realsum = test_dbl4_real_polynomial(dim,nbr,nva,pwr,deg,vrblvl-1);

   const double tol = 1.0e-60;

   int fail = int(realsum > tol);

   if(vrblvl > 0)
   {
      cout << scientific << setprecision(2);
      cout << "Sum of all errors in quad double precision :" << endl;
      cout << "  on real data : " << realsum;
      if(realsum < tol)
         cout << "  pass." << endl;
      else
         cout << "  fail!" << endl;

      cout << "  Seed used : " <<  seedused << endl;
   }
   return fail;
}

double test_dbl4_real_polynomial
 ( int dim, int nbr, int nva, int pwr, int deg, int verbose )
{
   if(nbr < 1)
      return 0.0;
   else
   {
      double **inputhihi = new double*[dim]; // dim series of degree deg
      for(int i=0; i<dim; i++) inputhihi[i] = new double[deg+1];
      double **inputlohi = new double*[dim];
      for(int i=0; i<dim; i++) inputlohi[i] = new double[deg+1];
      double **inputhilo = new double*[dim];
      for(int i=0; i<dim; i++) inputhilo[i] = new double[deg+1];
      double **inputlolo = new double*[dim];
      for(int i=0; i<dim; i++) inputlolo[i] = new double[deg+1];
      // The output are dim+1 power series of degree deg
      // for the evaluated and differentiated polynomial.
      double **output1hihi_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output1hihi_h[i] = new double[deg+1];
      double **output1lohi_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output1lohi_h[i] = new double[deg+1];
      double **output1hilo_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output1hilo_h[i] = new double[deg+1];
      double **output1lolo_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output1lolo_h[i] = new double[deg+1];
      double **output2hihi_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output2hihi_h[i] = new double[deg+1];
      double **output2lohi_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output2lohi_h[i] = new double[deg+1];
      double **output2hilo_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output2hilo_h[i] = new double[deg+1];
      double **output2lolo_h = new double*[dim+1];
      for(int i=0; i<=dim; i++) output2lolo_h[i] = new double[deg+1];
      double **outputhihi_d = new double*[dim+1];
      for(int i=0; i<=dim; i++) outputhihi_d[i] = new double[deg+1];
      double **outputlohi_d = new double*[dim+1];
      for(int i=0; i<=dim; i++) outputlohi_d[i] = new double[deg+1];
      double **outputhilo_d = new double*[dim+1];
      for(int i=0; i<=dim; i++) outputhilo_d[i] = new double[deg+1];
      double **outputlolo_d = new double*[dim+1];
      for(int i=0; i<=dim; i++) outputlolo_d[i] = new double[deg+1];

      make_real4_input(dim,deg,inputhihi,inputlohi,inputhilo,inputlolo);

      if(verbose > 0)
      {
         cout << scientific << setprecision(16);
         cout << "Random input series :" << endl;
         for(int i=0; i<dim; i++)
         {
            cout << "-> coefficients of series " << i << " :" << endl;
            for(int j=0; j<=deg; j++)
               cout << inputhihi[i][j] << "  " << inputlohi[i][j] << endl
                    << inputhilo[i][j] << "  " << inputlolo[i][j] << endl;
         }
      }
      double *csthihi = new double[deg+1]; // constant coefficient series
      double *cstlohi = new double[deg+1];
      double *csthilo = new double[deg+1];
      double *cstlolo = new double[deg+1];
      double **cffhihi = new double*[nbr]; // coefficient series of terms
      for(int i=0; i<nbr; i++) cffhihi[i] = new double[deg+1];
      double **cfflohi = new double*[nbr];
      for(int i=0; i<nbr; i++) cfflohi[i] = new double[deg+1];
      double **cffhilo = new double*[nbr];
      for(int i=0; i<nbr; i++) cffhilo[i] = new double[deg+1];
      double **cfflolo = new double*[nbr];
      for(int i=0; i<nbr; i++) cfflolo[i] = new double[deg+1];
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
            make_real4_cyclic
               (dim,nva,deg,idx,csthihi,cstlohi,csthilo,cstlolo,
                                cffhihi,cfflohi,cffhilo,cfflolo);
         else
            make_real4_products
               (dim,nbr,nva,deg,idx,csthihi,cstlohi,csthilo,cstlolo,
                                    cffhihi,cfflohi,cffhilo,cfflolo);
      }
      else
      {
         for(int i=0; i<nbr; i++) exp[i] = new int[nvr[i]];

         bool fail = make_real4_polynomial
                        (dim,nbr,pwr,deg,nvr,idx,exp,
                         csthihi,cstlohi,csthilo,cstlolo,
                         cffhihi,cfflohi,cffhilo,cfflolo);
      }
      if(verbose > 0)
      {
         cout << "Coefficient series of the constant term :" << endl;
         for(int j=0; j<=deg; j++)
            cout << csthihi[j] << "  " << cstlohi[j] << endl
                 << csthilo[j] << "  " << cstlolo[j] << endl;

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
               cout << cffhihi[i][j] << "  " << cfflohi[i][j] << endl
                    << cffhilo[i][j] << "  " << cfflolo[i][j] << endl;
         }
      }
      bool vrb = (verbose > 0);
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
      if(vrb) cout << "Computing without convolution jobs ..." << endl;
      CPU_dbl4_poly_evaldiff
         (dim,nbr,deg,nvr,idx,csthihi,cstlohi,csthilo,cstlolo,
          cffhihi,cfflohi,cffhilo,cfflolo,
          inputhihi,inputlohi,inputhilo,inputlolo,
          output1hihi_h,output1lohi_h,output1hilo_h,output1lolo_h,vrb);
      if(vrb) cout << "Computing with convolution jobs ..." << endl;
      CPU_dbl4_poly_evaldiffjobs
         (dim,nbr,deg,nvr,idx,csthihi,cstlohi,csthilo,cstlolo,
          cffhihi,cfflohi,cffhilo,cfflolo,
          inputhihi,inputlohi,inputhilo,inputlolo,
          output2hihi_h,output2lohi_h,output2hilo_h,output2lolo_h,
          cnvjobs,addjobs,vrb);
      if(vrb) cout << "Computing on the device ..." << endl;
      GPU_dbl4_poly_evaldiff
         (deg+1,dim,nbr,deg,nvr,idx,csthihi,cstlohi,csthilo,cstlolo,
          cffhihi,cfflohi,cffhilo,cfflolo,
          inputhihi,inputlohi,inputhilo,inputlolo,
          outputhihi_d,outputlohi_d,outputhilo_d,outputlolo_d,
          cnvjobs,addjobs,vrb);

      double err = 0.0;

      if(vrb) cout << "The value of the polynomial :" << endl;
      for(int i=0; i<=deg; i++)
      {
         if(vrb)
         {
            cout << output1hihi_h[dim][i] << "  "
                 << output1lohi_h[dim][i] << endl
                 << output1hilo_h[dim][i] << "  "
                 << output1lolo_h[dim][i] << endl;
            cout << output2hihi_h[dim][i] << "  "
                 << output2lohi_h[dim][i] << endl
                 << output2hilo_h[dim][i] << "  "
                 << output2lolo_h[dim][i] << endl;
            cout << outputhihi_d[dim][i] << "  "
                 << outputlohi_d[dim][i] << endl
                 << outputhilo_d[dim][i] << "  "
                 << outputlolo_d[dim][i] << endl;
         }
         err = err + abs(output1hihi_h[dim][i] - output2hihi_h[dim][i])
                   + abs(output1lohi_h[dim][i] - output2lohi_h[dim][i])
                   + abs(output1hilo_h[dim][i] - output2hilo_h[dim][i])
                   + abs(output1lolo_h[dim][i] - output2lolo_h[dim][i])
                   + abs(output1hihi_h[dim][i] - outputhihi_d[dim][i])
                   + abs(output1lohi_h[dim][i] - outputlohi_d[dim][i])
                   + abs(output1hilo_h[dim][i] - outputhilo_d[dim][i])
                   + abs(output1lolo_h[dim][i] - outputlolo_d[dim][i]);
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
               cout << output1hihi_h[k][i] << "  "
                    << output1lohi_h[k][i] << endl
                    << output1hilo_h[k][i] << "  "
                    << output1lolo_h[k][i] << endl;
               cout << output2hihi_h[k][i] << "  "
                    << output2lohi_h[k][i] << endl
                    << output2hilo_h[k][i] << "  "
                    << output2lolo_h[k][i] << endl;
               cout << outputhihi_d[k][i] << "  "
                    << outputlohi_d[k][i] << endl
                    << outputhilo_d[k][i] << "  "
                    << outputlolo_d[k][i] << endl;
            }
            err = err + abs(output1hihi_h[k][i] - output2hihi_h[k][i])
                      + abs(output1lohi_h[k][i] - output2lohi_h[k][i])
                      + abs(output1hilo_h[k][i] - output2hilo_h[k][i])
                      + abs(output1lolo_h[k][i] - output2lolo_h[k][i])
                      + abs(output1hihi_h[k][i] - outputhihi_d[k][i])
                      + abs(output1lohi_h[k][i] - outputlohi_d[k][i])
                      + abs(output1hilo_h[k][i] - outputhilo_d[k][i])
                      + abs(output1lolo_h[k][i] - outputlolo_d[k][i]);
         }
         if(vrb) cout << "error : " << err << endl;
         sumerr = sumerr + err;
      }
      if(vrb)
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
      return sumerr;
   }
}
