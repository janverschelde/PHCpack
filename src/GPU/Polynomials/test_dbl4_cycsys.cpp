/* Tests evaluation and differentiation of the cyclic n-roots system 
 * in quad double precision with the polynomial system data structure. */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <vector_types.h>
#include "random4_vectors.h"
#include "random4_monomials.h"
#include "random4_series.h"
#include "random4_polynomials.h"
#include "job_makers.h"
#include "dbl4_polynomials_host.h"
#include "dbl4_polynomials_kernels.h"
#include "dbl4_polynomials_testers.h"
#include "cyclic_indices.h"
#include "dbl4_indexed_coefficients.h"

using namespace std;

double test_dbl4_sysevaldiff
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

double test_cmplx4_sysevaldiff
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

   const double tol = 1.0e-48;
   double realsum = test_dbl4_sysevaldiff(dim,deg,nbr,nvr,idx,vrblvl);
   double compsum = test_cmplx4_sysevaldiff(dim,deg,nbr,nvr,idx,vrblvl);

   int fail = int(realsum > tol) + int(compsum > tol);

   cout << scientific << setprecision(2);
   cout << "Sum of all errors in quad double precision :" << endl;
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

double test_dbl4_sysevaldiff
 ( int dim, int deg, int *nbr, int **nvr, int ***idx, int vrblvl )
{
   const int degp1 = deg+1;
   double **csthihi = new double*[dim];
   double **cstlohi = new double*[dim];
   double **csthilo = new double*[dim];
   double **cstlolo = new double*[dim];
   double ***cffhihi = new double**[dim];
   double ***cfflohi = new double**[dim];
   double ***cffhilo = new double**[dim];
   double ***cfflolo = new double**[dim];

   dbl4_make_coefficients
      (dim,deg,nbr,nvr,idx,csthihi,cstlohi,csthilo,cstlolo,
       cffhihi,cfflohi,cffhilo,cfflolo,vrblvl);

/* define the input and allocate the output */

   double **inputhihi = new double*[dim]; // dim series of degree deg
   double **inputlohi = new double*[dim];
   double **inputhilo = new double*[dim];
   double **inputlolo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      inputhihi[i] = new double[deg+1];
      inputlohi[i] = new double[deg+1];
      inputhilo[i] = new double[deg+1];
      inputlolo[i] = new double[deg+1];
   }
   make_real4_input(dim,deg,inputhihi,inputlohi,inputhilo,inputlolo);

   if(vrblvl > 1)
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
   double ***outputhihi_h = new double**[dim];
   double ***outputlohi_h = new double**[dim];
   double ***outputhilo_h = new double**[dim];
   double ***outputlolo_h = new double**[dim];
   double ***outputhihi_d = new double**[dim];
   double ***outputlohi_d = new double**[dim];
   double ***outputhilo_d = new double**[dim];
   double ***outputlolo_d = new double**[dim];

   for(int i=0; i<dim; i++)
   {
      outputhihi_h[i] = new double*[dim+1];
      outputlohi_h[i] = new double*[dim+1];
      outputhilo_h[i] = new double*[dim+1];
      outputlolo_h[i] = new double*[dim+1];
      outputhihi_d[i] = new double*[dim+1];
      outputlohi_d[i] = new double*[dim+1];
      outputhilo_d[i] = new double*[dim+1];
      outputlolo_d[i] = new double*[dim+1];

      for(int j=0; j<=dim; j++)
      {
         outputhihi_h[i][j] = new double[degp1];
         outputlohi_h[i][j] = new double[degp1];
         outputhilo_h[i][j] = new double[degp1];
         outputlolo_h[i][j] = new double[degp1];
         outputhihi_d[i][j] = new double[degp1];
         outputlohi_d[i][j] = new double[degp1];
         outputhilo_d[i][j] = new double[degp1];
         outputlolo_d[i][j] = new double[degp1];
      }
   }
   if(vrblvl > 0) cout << "computing on the host ..." << endl;

   double timelapsed_h = 0.0;

   for(int i=0; i<dim; i++)
   {
      double lapsed;

      CPU_dbl4_poly_evaldiff
        (dim,nbr[i],deg,nvr[i],idx[i],
         csthihi[i],cstlohi[i],csthilo[i],cstlolo[i],
         cffhihi[i],cfflohi[i],cffhilo[i],cfflolo[i],
         inputhihi,inputlohi,inputhilo,inputlolo,
         outputhihi_h[i],outputlohi_h[i],outputhilo_h[i],outputlolo_h[i],
         &lapsed,vrblvl);

      timelapsed_h += lapsed;
   }
   if(vrblvl > 0) cout << "computing on the device ..." << endl;

   double timelapsed_d = 0.0;
   const bool vrb = (vrblvl > 1);

   for(int i=0; i<dim; i++)
   {
      ConvolutionJobs cnvjobs(dim);
      AdditionJobs addjobs(dim,nbr[i]);

      make_all_jobs(dim,nbr[i],nvr[i],idx[i],&cnvjobs,&addjobs,vrb);

      double cnvlapms,addlapms,timelapms_d,walltimes_d;

      GPU_dbl4_poly_evaldiff
         (degp1,dim,nbr[i],deg,nvr[i],idx[i],
          csthihi[i],cstlohi[i],csthilo[i],cstlolo[i],
          cffhihi[i],cfflohi[i],cffhilo[i],cfflolo[i],
          inputhihi,inputlohi,inputhilo,inputlolo,
          outputhihi_d[i],outputlohi_d[i],outputhilo_d[i],outputlolo_d[i],
          cnvjobs,addjobs,&cnvlapms,&addlapms,&timelapms_d,
          &walltimes_d,vrblvl);

      timelapsed_d += walltimes_d;
   }
   if(vrblvl > 0) cout << "computing the errors ..." << endl;

   double sumerr = 0.0;

   for(int i=0; i<dim; i++)
   {
      sumerr += dbl4_error_sum1(dim,deg,
                   outputhihi_h[i],outputlohi_h[i],
                   outputhilo_h[i],outputlolo_h[i],
                   outputhihi_d[i],outputlohi_d[i],
                   outputhilo_d[i],outputlolo_d[i],vrb);
   }
   cout << scientific << setprecision(2)
        << "sum of all errors " << sumerr << endl;

   return sumerr;
}

double test_cmplx4_sysevaldiff
 ( int dim, int deg, int *nbr, int **nvr, int ***idx, int vrblvl )
{
   const int degp1 = deg+1;
   double **cstrehihi = new double*[dim];
   double **cstrelohi = new double*[dim];
   double **cstrehilo = new double*[dim];
   double **cstrelolo = new double*[dim];
   double **cstimhihi = new double*[dim];
   double **cstimlohi = new double*[dim];
   double **cstimhilo = new double*[dim];
   double **cstimlolo = new double*[dim];
   double ***cffrehihi = new double**[dim];
   double ***cffrelohi = new double**[dim];
   double ***cffrehilo = new double**[dim];
   double ***cffrelolo = new double**[dim];
   double ***cffimhihi = new double**[dim];
   double ***cffimlohi = new double**[dim];
   double ***cffimhilo = new double**[dim];
   double ***cffimlolo = new double**[dim];

   cmplx4_make_coefficients
      (dim,deg,nbr,nvr,idx,
       cstrehihi,cstrelohi,cstrehilo,cstrelolo,
       cstimhihi,cstimlohi,cstimhilo,cstimlolo,
       cffrehihi,cffrelohi,cffrehilo,cffrelolo,
       cffimhihi,cffimlohi,cffimhilo,cffimlolo,vrblvl);

/* generate input series and allocate the output */

   double **inputrehihi = new double*[dim]; // dim series of degree deg
   double **inputrelohi = new double*[dim];
   double **inputrehilo = new double*[dim];
   double **inputrelolo = new double*[dim];
   double **inputimhihi = new double*[dim];
   double **inputimlohi = new double*[dim];
   double **inputimhilo = new double*[dim];
   double **inputimlolo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      inputrehihi[i] = new double[degp1];
      inputrelohi[i] = new double[degp1];
      inputrehilo[i] = new double[degp1];
      inputrelolo[i] = new double[degp1];
      inputimhihi[i] = new double[degp1];
      inputimlohi[i] = new double[degp1];
      inputimhilo[i] = new double[degp1];
      inputimlolo[i] = new double[degp1];
   }
   make_complex4_input
      (dim,deg,inputrehihi,inputrelohi,inputrehilo,inputrelolo,
               inputimhihi,inputimlohi,inputimhilo,inputimlolo);

   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
      cout << "Random input series :" << endl;
      for(int i=0; i<dim; i++)
      {
         cout << "-> coefficients of series " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
            cout << inputrehihi[i][j] << "  " << inputrelohi[i][j] << endl
                 << inputrehilo[i][j] << "  " << inputrelolo[i][j] << endl
                 << inputimhihi[i][j] << "  " << inputimlohi[i][j] << endl
                 << inputimhilo[i][j] << "  " << inputimlolo[i][j] << endl;
      }
   }
   double ***outputrehihi_h = new double**[dim];
   double ***outputrelohi_h = new double**[dim];
   double ***outputrehilo_h = new double**[dim];
   double ***outputrelolo_h = new double**[dim];
   double ***outputimhihi_h = new double**[dim];
   double ***outputimlohi_h = new double**[dim];
   double ***outputimhilo_h = new double**[dim];
   double ***outputimlolo_h = new double**[dim];
   double ***outputrehihi_d = new double**[dim];
   double ***outputrelohi_d = new double**[dim];
   double ***outputrehilo_d = new double**[dim];
   double ***outputrelolo_d = new double**[dim];
   double ***outputimhihi_d = new double**[dim];
   double ***outputimlohi_d = new double**[dim];
   double ***outputimhilo_d = new double**[dim];
   double ***outputimlolo_d = new double**[dim];

   for(int i=0; i<dim; i++)
   {
      outputrehihi_h[i] = new double*[dim+1];
      outputrelohi_h[i] = new double*[dim+1];
      outputrehilo_h[i] = new double*[dim+1];
      outputrelolo_h[i] = new double*[dim+1];
      outputimhihi_h[i] = new double*[dim+1];
      outputimlohi_h[i] = new double*[dim+1];
      outputimhilo_h[i] = new double*[dim+1];
      outputimlolo_h[i] = new double*[dim+1];
      outputrehihi_d[i] = new double*[dim+1];
      outputrelohi_d[i] = new double*[dim+1];
      outputrehilo_d[i] = new double*[dim+1];
      outputrelolo_d[i] = new double*[dim+1];
      outputimhihi_d[i] = new double*[dim+1];
      outputimlohi_d[i] = new double*[dim+1];
      outputimhilo_d[i] = new double*[dim+1];
      outputimlolo_d[i] = new double*[dim+1];

      for(int j=0; j<=dim; j++)
      {
         outputrehihi_h[i][j] = new double[degp1];
         outputrelohi_h[i][j] = new double[degp1];
         outputrehilo_h[i][j] = new double[degp1];
         outputrelolo_h[i][j] = new double[degp1];
         outputimhihi_h[i][j] = new double[degp1];
         outputimlohi_h[i][j] = new double[degp1];
         outputimhilo_h[i][j] = new double[degp1];
         outputimlolo_h[i][j] = new double[degp1];
         outputrehihi_d[i][j] = new double[degp1];
         outputrelohi_d[i][j] = new double[degp1];
         outputrehilo_d[i][j] = new double[degp1];
         outputrelolo_d[i][j] = new double[degp1];
         outputimhihi_d[i][j] = new double[degp1];
         outputimlohi_d[i][j] = new double[degp1];
         outputimhilo_d[i][j] = new double[degp1];
         outputimlolo_d[i][j] = new double[degp1];
      }
   }
   if(vrblvl > 0) cout << "evaluating on the host ..." << endl;

   double timelapsed_h = 0.0;

   for(int i=0; i<dim; i++)
   {
      double lapsed;

      CPU_cmplx4_poly_evaldiff
        (dim,nbr[i],deg,nvr[i],idx[i],
         cstrehihi[i],cstrelohi[i],cstrehilo[i],cstrelolo[i],
         cstimhihi[i],cstimlohi[i],cstimhilo[i],cstimlolo[i],
         cffrehihi[i],cffrelohi[i],cffrehilo[i],cffrelolo[i],
         cffimhihi[i],cffimlohi[i],cffimhilo[i],cffimlolo[i],
         inputrehihi,inputrelohi,inputrehilo,inputrelolo,
         inputimhihi,inputimlohi,inputimhilo,inputimlolo,
         outputrehihi_h[i],outputrelohi_h[i],
         outputrehilo_h[i],outputrelolo_h[i],
         outputimhihi_h[i],outputimlohi_h[i],
         outputimhilo_h[i],outputimlolo_h[i],&lapsed,vrblvl);

      timelapsed_h += lapsed;
   }
   if(vrblvl > 0) cout << "computing on the device ..." << endl;

   double timelapsed_d = 0.0;
   const bool vrb = (vrblvl > 1);

   for(int i=0; i<dim; i++)
   {
      ComplexConvolutionJobs cnvjobs(dim);
      ComplexIncrementJobs incjobs(cnvjobs,vrb);
      ComplexAdditionJobs addjobs(dim,nbr[i]);

      make_all_complex_jobs
         (dim,nbr[i],nvr[i],idx[i],&cnvjobs,&incjobs,&addjobs,vrb);

      double cnvlapms,addlapms,timelapms_d,walltimes_d;

      GPU_cmplx4vectorized_poly_evaldiff
         (degp1,dim,nbr[i],deg,nvr[i],idx[i],
          cstrehihi[i],cstrelohi[i],cstrehilo[i],cstrelolo[i],
          cstimhihi[i],cstimlohi[i],cstimhilo[i],cstimlolo[i],
          cffrehihi[i],cffrelohi[i],cffrehilo[i],cffrelolo[i],
          cffimhihi[i],cffimlohi[i],cffimhilo[i],cffimlolo[i],
          inputrehihi,inputrelohi,inputrehilo,inputrelolo,
          inputimhihi,inputimlohi,inputimhilo,inputimlolo,
          outputrehihi_d[i],outputrelohi_d[i],
          outputrehilo_d[i],outputrelolo_d[i],
          outputimhihi_d[i],outputimlohi_d[i],
          outputimhilo_d[i],outputimlolo_d[i],
          cnvjobs,incjobs,addjobs,&cnvlapms,&addlapms,&timelapms_d,
          &walltimes_d,vrblvl);

      timelapsed_d += walltimes_d;
   }
   if(vrblvl > 0) cout << "computing the errors ..." << endl;

   double sumerr = 0.0;

   for(int i=0; i<dim; i++)
   {
      double err = cmplx4_error_sum1(dim,deg,
                   outputrehihi_h[i],outputrelohi_h[i],
                   outputrehilo_h[i],outputrelolo_h[i],
                   outputimhihi_h[i],outputimlohi_h[i],
                   outputimhilo_h[i],outputimlolo_h[i],
                   outputrehihi_d[i],outputrelohi_d[i],
                   outputrehilo_d[i],outputrelolo_d[i],
                   outputimhihi_d[i],outputimlohi_d[i],
                   outputimhilo_d[i],outputimlolo_d[i],vrb);
      sumerr += err;
      if(err > 1.0e-48) cout << "large error at i = " << i << endl;
   }

   cout << scientific << setprecision(2)
        << "sum of all errors " << sumerr << endl;

   return sumerr;
}
