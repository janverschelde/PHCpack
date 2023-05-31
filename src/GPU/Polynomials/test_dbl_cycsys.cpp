/* Tests evaluation and differentiation of the cyclic n-roots system 
 * in double precision with the polynomial system data structure. */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <vector_types.h>
#include "random_monomials.h"
#include "random_series.h"
#include "random_polynomials.h"
#include "job_makers.h"
#include "dbl_polynomials_host.h"
#include "dbl_polynomials_kernels.h"
#include "dbl_polynomials_testers.h"

using namespace std;

void write_polynomial_indices ( int dim, int *nbr, int **nvr, int ***idx );
/*
 * DESCRIPTION :
 *   Writes the indices of a polynomial system.
 *
 * ON ENTRY :
 *   dim      number of polynomials in the system;
 *   nbr      nbr[i] is the number of monomials in the i-th polynomial;
 *   nvr      nvr[i][j] is the number of variables of the j-th monomial
 *            of the i-th polynomial;
 *   idx      idx[i][j][k] is the index to the k-th variable
 *            in the j-th monomial of the i-th polynomial.    */

double test_dbl_sysevaldiff
 ( int dim, int deg, int *nbr, int **nvr, int ***idx, int vrblvl );
/*
 * DESCRIPTION :
 *   Tests the evaluation and differentiation for random real data.
 *   Returns the sum of all errors.
 * 
 * ON ENTRY :
 *   dim      dimension, total number of variables and polynomials;
 *   deg      truncation degree of the series;
 *   nbr      nbr[i] is the number of monomials in the i-th polynomial;
 *   nvr      nvr[i][j] is the number of variables of the j-th monomial
 *            of the i-th polynomial;
 *   idx      idx[i][j][k] is the index to the k-th variable
 *            in the j-th monomial of the i-th polynomial;
 *   vrblvl   is the  verbose level. */

double test_cmplx_sysevaldiff
 ( int dim, int deg, int *nbr, int **nvr, int ***idx, int vrblvl );
/*
 * DESCRIPTION :
 *   Tests the evaluation and differentiation for random complex data.
 *   Returns the sum of all errors.
 * 
 * ON ENTRY :
 *   dim      dimension, total number of variables and polynomials;
 *   deg      truncation degree of the series;
 *   nbr      nbr[i] is the number of monomials in the i-th polynomial;
 *   nvr      nvr[i][j] is the number of variables of the j-th monomial
 *            of the i-th polynomial;
 *   idx      idx[i][j][k] is the index to the k-th variable
 *            in the j-th monomial of the i-th polynomial;
 *   vrblvl   is the  verbose level. */

int main ( void )
{
   cout << "evaluation and differentiation of cyclic n-roots ..." << endl;

   cout << "-> give the seed (0 for time) : ";
   int seed; cin >> seed;
   cout << "-> give the dimension : ";
   int dim; cin >> dim;
   cout << "-> give the degree : ";
   int deg; cin >> deg;
   cout << "-> give the verbose level : ";
   int vrblvl; cin >> vrblvl;

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

   int *nbr = new int[dim];
   for(int i=0; i<dim-1; i++) nbr[i] = dim;
   nbr[dim-1] = 1;

   int **nvr = new int*[dim]; // number of variables in dim polynomials
   for(int i=0; i<dim-1; i++)
   {
      nvr[i] = new int[dim];
      for(int j=0; j<dim; j++) nvr[i][j] = i+1;
   }
   nvr[dim-1] = new int[1]; // 2 monomials in the last polynomial
   nvr[dim-1][0] = dim;     // but constant is stored separately

   int ***idx = new int**[dim]; // we have dim polynomials and
   for(int i=0; i<dim-1; i++)   // dim monomials in each polynomial
   {
      idx[i] = new int*[dim];
      for(int j=0; j<dim; j++)
      {
         idx[i][j] = new int[nvr[i][j]];
         for(int k=0; k<nvr[i][j]; k++) idx[i][j][k] = (j + k) % dim;
         insert_sort(nvr[i][j],idx[i][j]);
      }
   }
   idx[dim-1] = new int*[1];
   idx[dim-1][0] = new int[dim]; // except for the last monomial
   for(int k=0; k<dim; k++) idx[dim-1][0][k] = k;

   write_polynomial_indices(dim,nbr,nvr,idx);

   const double tol = 1.0e-10;
   double realsum = test_dbl_sysevaldiff(dim,deg,nbr,nvr,idx,vrblvl);
   double compsum = test_cmplx_sysevaldiff(dim,deg,nbr,nvr,idx,vrblvl);

   int fail = int(realsum > tol) + int(compsum > tol);

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
   cout << "  Seed used : " <<  seedused << endl;

   return fail;
}

void write_polynomial_indices ( int dim, int *nbr, int **nvr, int ***idx )
{
   for(int i=0; i<dim; i++)
   {
      cout << "polynomial " << i
           << " has " << nbr[i] << " monomials :" << endl;
      for(int j=0; j<nbr[i]; j++)
      {
         for(int k=0; k<nvr[i][j]; k++) cout << " " << idx[i][j][k];
         cout << endl;
      }
   }
}

double test_dbl_sysevaldiff
 ( int dim, int deg, int *nbr, int **nvr, int ***idx, int vrblvl )
{
   const int degp1 = deg+1;
/* generate constant and coefficients */
   double rnd;
   double **cst = new double*[dim];
   for(int i=0; i<dim; i++)
   {
      cst[i] = new double[deg+1];
      random_dbl_exponential(deg,&rnd,cst[i]);
   }
   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
      cout << "Constant coefficient series :" << endl;
      for(int i=0; i<dim; i++)
      {
         cout << "-> coefficients of series " << i << " :" << endl;
         for(int j=0; j<=deg; j++) cout << cst[i][j] << endl;
      }
   }
   double ***cff = new double**[dim];
   for(int i=0; i<dim; i++)
   {
      cff[i] = new double*[nbr[i]];
      for(int j=0; j<nbr[i]; j++)
      {
         cff[i][j] = new double[deg+1];
         random_dbl_exponential(deg,&rnd,cff[i][j]);
      }
   }
   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
      for(int i=0; i<dim; i++)
      {
         for(int j=0; j<nbr[i]; j++)
         {
            cout << "-> coefficients of monomial " << j
                 << " of polynomial " << i << " :" << endl;
            for(int k=0; k<=deg; k++) cout << cff[i][j][k] << endl;
         }
      }
   }
/* define the input and allocate the output */
   double **input = new double*[dim]; // dim series of degree deg
   for(int i=0; i<dim; i++) input[i] = new double[deg+1];
   make_real_input(dim,deg,input);

   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
      cout << "Random input series :" << endl;
      for(int i=0; i<dim; i++)
      {
         cout << "-> coefficients of series " << i << " :" << endl;
         for(int j=0; j<=deg; j++) cout << input[i][j] << endl;
      }
   }
   double ***output_h = new double**[dim];
   double ***output_d = new double**[dim];

   for(int i=0; i<dim; i++)
   {
      output_h[i] = new double*[dim+1];
      output_d[i] = new double*[dim+1];

      for(int j=0; j<=dim; j++)
      {
         output_h[i][j] = new double[degp1];
         output_d[i][j] = new double[degp1];
      }
   }
   if(vrblvl > 0) cout << "computing on the host ..." << endl;

   double timelapsed_h = 0.0;

   for(int i=0; i<dim; i++)
   {
      double lapsed;

      CPU_dbl_poly_evaldiff
        (dim,nbr[i],deg,nvr[i],idx[i],cst[i],cff[i],input,
         output_h[i],&lapsed,vrblvl);

      timelapsed_h += lapsed;
   }
   if(vrblvl > 0) cout << "computing on the device ..." << endl;

   double timelapsed_d = 0.0;

   for(int i=0; i<dim; i++)
   {
      ConvolutionJobs cnvjobs(dim);
      AdditionJobs addjobs(dim,nbr[i]);

      make_all_jobs(dim,nbr[i],nvr[i],idx[i],&cnvjobs,&addjobs,vrblvl);

      double cnvlapms,addlapms,timelapms_d,walltimes_d;

      GPU_dbl_poly_evaldiff
         (degp1,dim,nbr[i],deg,nvr[i],idx[i],cst[i],cff[i],input,
          output_d[i],cnvjobs,addjobs,&cnvlapms,&addlapms,&timelapms_d,
          &walltimes_d,vrblvl);

      timelapsed_d += walltimes_d;
   }
   if(vrblvl > 0) cout << "computing the errors ..." << endl;

   double sumerr = 0.0;

   for(int i=0; i<dim; i++)
   {
      sumerr += dbl_error_sum1(dim,deg,output_h[i],output_d[i],vrblvl);
   }
   cout << "sum of all errors " << sumerr << endl;

   return sumerr;
}

double test_cmplx_sysevaldiff
 ( int dim, int deg, int *nbr, int **nvr, int ***idx, int vrblvl )
{
   const int degp1 = deg+1;
/* generate constant and coefficients */
   double rndre,rndim;
   double **cstre = new double*[dim];
   double **cstim = new double*[dim];
   for(int i=0; i<dim; i++)
   {
      cstre[i] = new double[degp1];
      cstim[i] = new double[degp1];
      random_cmplx_exponential(deg,&rndre,&rndim,cstre[i],cstim[i]);
   }
   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
      cout << "Constant coefficient series :" << endl;
      for(int i=0; i<dim; i++)
      {
         cout << "-> coefficients of series " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
            cout << cstre[i][j] << "  " << cstim[i][j] << endl;
      }
   }
   double ***cffre = new double**[dim];
   double ***cffim = new double**[dim];
   for(int i=0; i<dim; i++)
   {
      cffre[i] = new double*[nbr[i]];
      cffim[i] = new double*[nbr[i]];
      for(int j=0; j<nbr[i]; j++)
      {
         cffre[i][j] = new double[degp1];
         cffim[i][j] = new double[degp1];
         random_cmplx_exponential(deg,&rndre,&rndim,cffre[i][j],cffim[i][j]);
      }
   }
   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
      for(int i=0; i<dim; i++)
      {
         for(int j=0; j<nbr[i]; j++)
         {
            cout << "-> coefficients of monomial " << j
                 << " of polynomial " << i << " :" << endl;
            for(int k=0; k<=deg; k++)
               cout << cffre[i][j][k] << "  " << cffim[i][j][k] << endl;
         }
      }
   }
/* generate input series and allocate the output */
   double **inputre = new double*[dim]; // dim series of degree deg
   double **inputim = new double*[dim]; // dim series of degree deg
   for(int i=0; i<dim; i++)
   {
      inputre[i] = new double[degp1];
      inputim[i] = new double[degp1];
   }
   make_complex_input(dim,deg,inputre,inputim);

   if(vrblvl > 1)
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
   double ***outputre_h = new double**[dim];
   double ***outputim_h = new double**[dim];
   double ***outputre_d = new double**[dim];
   double ***outputim_d = new double**[dim];

   for(int i=0; i<dim; i++)
   {
      outputre_h[i] = new double*[dim+1];
      outputim_h[i] = new double*[dim+1];
      outputre_d[i] = new double*[dim+1];
      outputim_d[i] = new double*[dim+1];

      for(int j=0; j<=dim; j++)
      {
         outputre_h[i][j] = new double[degp1];
         outputim_h[i][j] = new double[degp1];
         outputre_d[i][j] = new double[degp1];
         outputim_d[i][j] = new double[degp1];
      }
   }
   if(vrblvl > 0) cout << "evaluating on the host ..." << endl;

   double timelapsed_h = 0.0;

   for(int i=0; i<dim; i++)
   {
      double lapsed;

      CPU_cmplx_poly_evaldiff
        (dim,nbr[i],deg,nvr[i],idx[i],cstre[i],cstim[i],cffre[i],cffim[i],
         inputre,inputim,outputre_h[i],outputim_h[i],&lapsed,vrblvl);

      timelapsed_h += lapsed;
   }
   if(vrblvl > 0) cout << "computing on the device ..." << endl;

   double timelapsed_d = 0.0;

   for(int i=0; i<dim; i++)
   {
      ComplexConvolutionJobs cnvjobs(dim);
      ComplexIncrementJobs incjobs(cnvjobs,vrblvl);
      ComplexAdditionJobs addjobs(dim,nbr[i]);

      make_all_complex_jobs
         (dim,nbr[i],nvr[i],idx[i],&cnvjobs,&incjobs,&addjobs,vrblvl);

      double cnvlapms,addlapms,timelapms_d,walltimes_d;

      GPU_cmplxvectorized_poly_evaldiff
         (degp1,dim,nbr[i],deg,nvr[i],idx[i],cstre[i],cstim[i],
          cffre[i],cffim[i],inputre,inputim,outputre_d[i],outputim_d[i],
          cnvjobs,incjobs,addjobs,&cnvlapms,&addlapms,&timelapms_d,
          &walltimes_d,vrblvl);

      timelapsed_d += walltimes_d;
   }
   if(vrblvl > 0) cout << "computing the errors ..." << endl;

   double sumerr = 0.0;

   for(int i=0; i<dim; i++)
   {
      sumerr += cmplx_error_sum1
                   (dim,deg,outputre_h[i],outputim_h[i],
                            outputre_d[i],outputim_d[i],vrblvl);
   }
   cout << "sum of all errors " << sumerr << endl;

   return sumerr;
}
