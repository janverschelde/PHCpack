/* Tests evaluation and differentiation of the cyclic n-roots system 
 * in double double precision with the polynomial system data structure. */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <vector_types.h>
#include "random2_monomials.h"
#include "random2_series.h"
#include "random2_polynomials.h"
#include "job_makers.h"
#include "dbl2_polynomials_host.h"
#include "dbl2_polynomials_kernels.h"
#include "dbl2_polynomials_testers.h"
#include "cyclic_indices.h"
#include "dbl2_indexed_coefficients.h"

using namespace std;

double test_dbl2_sysevaldiff
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

double test_cmplx2_sysevaldiff
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

   int *nbr = new int[dim];     // number of monomials in each polynomial
   int **nvr = new int*[dim];   // number of variables in dim polynomials
   int ***idx = new int**[dim]; // we have dim polynomials 

   make_polynomial_indices(dim,nbr,nvr,idx);
   write_polynomial_indices(dim,nbr,nvr,idx);

   const double tol = 1.0e-24;
   double realsum = test_dbl2_sysevaldiff(dim,deg,nbr,nvr,idx,vrblvl);
   double compsum = test_cmplx2_sysevaldiff(dim,deg,nbr,nvr,idx,vrblvl);

   int fail = int(realsum > tol) + int(compsum > tol);

   cout << scientific << setprecision(2);
   cout << "Sum of all errors in double double precision :" << endl;
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

double test_dbl2_sysevaldiff
 ( int dim, int deg, int *nbr, int **nvr, int ***idx, int vrblvl )
{
   const int degp1 = deg+1;
   double **csthi = new double*[dim];
   double **cstlo = new double*[dim];
   double ***cffhi = new double**[dim];
   double ***cfflo = new double**[dim];

   dbl2_make_coefficients
      (dim,deg,nbr,nvr,idx,csthi,cstlo,cffhi,cfflo,vrblvl);

/* define the input and allocate the output */
   double **inputhi = new double*[dim]; // dim series of degree deg
   double **inputlo = new double*[dim]; // dim series of degree deg
   for(int i=0; i<dim; i++)
   {
      inputhi[i] = new double[degp1];
      inputlo[i] = new double[degp1];
   }
   make_real2_input(dim,deg,inputhi,inputlo);

   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
      cout << "Random input series :" << endl;
      for(int i=0; i<dim; i++)
      {
         cout << "-> coefficients of series " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
            cout << inputhi[i][j] << "  " << inputlo[i][j] << endl;
      }
   }
   double ***outputhi_h = new double**[dim];
   double ***outputlo_h = new double**[dim];
   double ***outputhi_d = new double**[dim];
   double ***outputlo_d = new double**[dim];

   for(int i=0; i<dim; i++)
   {
      outputhi_h[i] = new double*[dim+1];
      outputlo_h[i] = new double*[dim+1];
      outputhi_d[i] = new double*[dim+1];
      outputlo_d[i] = new double*[dim+1];

      for(int j=0; j<=dim; j++)
      {
         outputhi_h[i][j] = new double[degp1];
         outputlo_h[i][j] = new double[degp1];
         outputhi_d[i][j] = new double[degp1];
         outputlo_d[i][j] = new double[degp1];
      }
   }
   if(vrblvl > 0) cout << "computing on the host ..." << endl;

   double timelapsed_h = 0.0;

   for(int i=0; i<dim; i++)
   {
      double lapsed;

      CPU_dbl2_poly_evaldiff
        (dim,nbr[i],deg,nvr[i],idx[i],csthi[i],cstlo[i],
         cffhi[i],cfflo[i],inputhi,inputlo,outputhi_h[i],outputlo_h[i],
         &lapsed,vrblvl);

      timelapsed_h += lapsed;
   }
   if(vrblvl > 0) cout << "computing on the device ..." << endl;

   double timelapsed_d = 0.0;
   bool vrb = (vrblvl > 1);

   for(int i=0; i<dim; i++)
   {
      ConvolutionJobs cnvjobs(dim);
      AdditionJobs addjobs(dim,nbr[i]);

      make_all_jobs(dim,nbr[i],nvr[i],idx[i],&cnvjobs,&addjobs,vrb);

      double cnvlapms,addlapms,timelapms_d,walltimes_d;

      GPU_dbl2_poly_evaldiff
         (degp1,dim,nbr[i],deg,nvr[i],idx[i],csthi[i],cstlo[i],
          cffhi[i],cfflo[i],inputhi,inputlo,outputhi_d[i],outputlo_d[i],
          cnvjobs,addjobs,&cnvlapms,&addlapms,&timelapms_d,
          &walltimes_d,vrblvl);

      timelapsed_d += walltimes_d;
   }
   if(vrblvl > 0) cout << "computing the errors ..." << endl;

   double sumerr = 0.0;

   for(int i=0; i<dim; i++)
   {
      sumerr += dbl2_error_sum1(dim,deg,outputhi_h[i],outputlo_h[i],
                                        outputhi_d[i],outputlo_d[i],vrb);
   }
   cout << scientific << setprecision(2)
        << "sum of all errors " << sumerr << endl;

   return sumerr;
}

double test_cmplx2_sysevaldiff
 ( int dim, int deg, int *nbr, int **nvr, int ***idx, int vrblvl )
{
   const int degp1 = deg+1;
   double **cstrehi = new double*[dim];
   double **cstrelo = new double*[dim];
   double **cstimhi = new double*[dim];
   double **cstimlo = new double*[dim];
   double ***cffrehi = new double**[dim];
   double ***cffrelo = new double**[dim];
   double ***cffimhi = new double**[dim];
   double ***cffimlo = new double**[dim];

   cmplx2_make_coefficients
      (dim,deg,nbr,nvr,idx,cstrehi,cstrelo,cstimhi,cstimlo,
       cffrehi,cffrelo,cffimhi,cffimlo,vrblvl);

/* generate input series and allocate the output */
   double **inputrehi = new double*[dim]; // dim series of degree deg
   double **inputrelo = new double*[dim];
   double **inputimhi = new double*[dim];
   double **inputimlo = new double*[dim];
   for(int i=0; i<dim; i++)
   {
      inputrehi[i] = new double[degp1];
      inputrelo[i] = new double[degp1];
      inputimhi[i] = new double[degp1];
      inputimlo[i] = new double[degp1];
   }
   make_complex2_input(dim,deg,inputrehi,inputrelo,inputimhi,inputimlo);

   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
      cout << "Random input series :" << endl;
      for(int i=0; i<dim; i++)
      {
         cout << "-> coefficients of series " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
            cout << inputrehi[i][j] << "  " << inputrelo[i][j] << endl
                 << inputimhi[i][j] << "  " << inputimlo[i][j] << endl;
      }
   }
   double ***outputrehi_h = new double**[dim];
   double ***outputrelo_h = new double**[dim];
   double ***outputimhi_h = new double**[dim];
   double ***outputimlo_h = new double**[dim];
   double ***outputrehi_d = new double**[dim];
   double ***outputrelo_d = new double**[dim];
   double ***outputimhi_d = new double**[dim];
   double ***outputimlo_d = new double**[dim];

   for(int i=0; i<dim; i++)
   {
      outputrehi_h[i] = new double*[dim+1];
      outputrelo_h[i] = new double*[dim+1];
      outputimhi_h[i] = new double*[dim+1];
      outputimlo_h[i] = new double*[dim+1];
      outputrehi_d[i] = new double*[dim+1];
      outputrelo_d[i] = new double*[dim+1];
      outputimhi_d[i] = new double*[dim+1];
      outputimlo_d[i] = new double*[dim+1];

      for(int j=0; j<=dim; j++)
      {
         outputrehi_h[i][j] = new double[degp1];
         outputrelo_h[i][j] = new double[degp1];
         outputimhi_h[i][j] = new double[degp1];
         outputimlo_h[i][j] = new double[degp1];
         outputrehi_d[i][j] = new double[degp1];
         outputrelo_d[i][j] = new double[degp1];
         outputimhi_d[i][j] = new double[degp1];
         outputimlo_d[i][j] = new double[degp1];
      }
   }
   if(vrblvl > 0) cout << "evaluating on the host ..." << endl;

   double timelapsed_h = 0.0;

   for(int i=0; i<dim; i++)
   {
      double lapsed;

      CPU_cmplx2_poly_evaldiff
        (dim,nbr[i],deg,nvr[i],idx[i],
         cstrehi[i],cstrelo[i],cstimhi[i],cstimlo[i],
         cffrehi[i],cffrelo[i],cffimhi[i],cffimlo[i],
         inputrehi,inputrelo,inputimhi,inputimlo,
         outputrehi_h[i],outputrelo_h[i],outputimhi_h[i],outputimlo_h[i],
         &lapsed,vrblvl);

      timelapsed_h += lapsed;
   }
   if(vrblvl > 0) cout << "computing on the device ..." << endl;

   double timelapsed_d = 0.0;
   bool vrb = (vrblvl > 1);

   for(int i=0; i<dim; i++)
   {
      ComplexConvolutionJobs cnvjobs(dim);
      ComplexIncrementJobs incjobs(cnvjobs,vrb);
      ComplexAdditionJobs addjobs(dim,nbr[i]);

      make_all_complex_jobs
         (dim,nbr[i],nvr[i],idx[i],&cnvjobs,&incjobs,&addjobs,vrb);

      double cnvlapms,addlapms,timelapms_d,walltimes_d;

      GPU_cmplx2vectorized_poly_evaldiff
         (degp1,dim,nbr[i],deg,nvr[i],idx[i],
          cstrehi[i],cstrelo[i],cstimhi[i],cstimlo[i],
          cffrehi[i],cffrelo[i],cffimhi[i],cffimlo[i],
          inputrehi,inputrelo,inputimhi,inputimlo,
          outputrehi_d[i],outputrelo_d[i],outputimhi_d[i],outputimlo_d[i],
          cnvjobs,incjobs,addjobs,&cnvlapms,&addlapms,&timelapms_d,
          &walltimes_d,vrblvl);

      timelapsed_d += walltimes_d;
   }
   if(vrblvl > 0) cout << "computing the errors ..." << endl;

   double sumerr = 0.0;

   for(int i=0; i<dim; i++)
   {
      sumerr += cmplx2_error_sum1(dim,deg,
                   outputrehi_h[i],outputrelo_h[i],
                   outputimhi_h[i],outputimlo_h[i],
                   outputrehi_d[i],outputrelo_d[i],
                   outputimhi_d[i],outputimlo_d[i],vrb);
   }
   cout << scientific << setprecision(2)
        << "sum of all errors " << sumerr << endl;

   return sumerr;
}
