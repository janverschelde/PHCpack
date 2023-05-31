/* Tests evaluation and differentiation of the cyclic n-roots system 
 * in octo double precision with the polynomial system data structure. */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <vector_types.h>
#include "random8_vectors.h"
#include "random8_monomials.h"
#include "random8_series.h"
#include "random8_polynomials.h"
#include "job_makers.h"
#include "dbl8_polynomials_host.h"
#include "dbl8_polynomials_kernels.h"
#include "dbl8_polynomials_testers.h"
#include "cyclic_indices.h"

using namespace std;

double test_dbl8_sysevaldiff
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

double test_cmplx8_sysevaldiff
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

   const double tol = 1.0e-100;
   double realsum = test_dbl8_sysevaldiff(dim,deg,nbr,nvr,idx,vrblvl);
   double compsum = test_cmplx8_sysevaldiff(dim,deg,nbr,nvr,idx,vrblvl);

   int fail = int(realsum > tol) + int(compsum > tol);

   cout << scientific << setprecision(2);
   cout << "Sum of all errors in octo double precision :" << endl;
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

double test_dbl8_sysevaldiff
 ( int dim, int deg, int *nbr, int **nvr, int ***idx, int vrblvl )
{
   const int degp1 = deg+1;

/* generate constant and coefficients */

   double rndhihihi,rndlohihi,rndhilohi,rndlolohi;
   double rndhihilo,rndlohilo,rndhilolo,rndlololo;
   double **csthihihi = new double*[dim];
   double **cstlohihi = new double*[dim];
   double **csthilohi = new double*[dim];
   double **cstlolohi = new double*[dim];
   double **csthihilo = new double*[dim];
   double **cstlohilo = new double*[dim];
   double **csthilolo = new double*[dim];
   double **cstlololo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      csthihihi[i] = new double[degp1];
      cstlohihi[i] = new double[degp1];
      csthilohi[i] = new double[degp1];
      cstlolohi[i] = new double[degp1];
      csthihilo[i] = new double[degp1];
      cstlohilo[i] = new double[degp1];
      csthilolo[i] = new double[degp1];
      cstlololo[i] = new double[degp1];

      random_dbl8_exponential
         (deg,&rndhihihi,&rndlohihi,&rndhilohi,&rndlolohi,
              &rndhihilo,&rndlohilo,&rndhilolo,&rndlololo,
          csthihihi[i],cstlohihi[i],csthilohi[i],cstlolohi[i],
          csthihilo[i],cstlohilo[i],csthilolo[i],cstlololo[i]);
   }
   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
      cout << "Constant coefficient series :" << endl;
      for(int i=0; i<dim; i++)
      {
         cout << "-> coefficients of series " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
            cout << csthihihi[i][j] << "  " << cstlohihi[i][j] << endl
                 << csthilohi[i][j] << "  " << cstlolohi[i][j] << endl
                 << csthihilo[i][j] << "  " << cstlohilo[i][j] << endl
                 << csthilolo[i][j] << "  " << cstlololo[i][j] << endl;
      }
   }
   double ***cffhihihi = new double**[dim];
   double ***cfflohihi = new double**[dim];
   double ***cffhilohi = new double**[dim];
   double ***cfflolohi = new double**[dim];
   double ***cffhihilo = new double**[dim];
   double ***cfflohilo = new double**[dim];
   double ***cffhilolo = new double**[dim];
   double ***cfflololo = new double**[dim];

   for(int i=0; i<dim; i++)
   {
      cffhihihi[i] = new double*[nbr[i]];
      cfflohihi[i] = new double*[nbr[i]];
      cffhilohi[i] = new double*[nbr[i]];
      cfflolohi[i] = new double*[nbr[i]];
      cffhihilo[i] = new double*[nbr[i]];
      cfflohilo[i] = new double*[nbr[i]];
      cffhilolo[i] = new double*[nbr[i]];
      cfflololo[i] = new double*[nbr[i]];

      for(int j=0; j<nbr[i]; j++)
      {
         cffhihihi[i][j] = new double[degp1];
         cfflohihi[i][j] = new double[degp1];
         cffhilohi[i][j] = new double[degp1];
         cfflolohi[i][j] = new double[degp1];
         cffhihilo[i][j] = new double[degp1];
         cfflohilo[i][j] = new double[degp1];
         cffhilolo[i][j] = new double[degp1];
         cfflololo[i][j] = new double[degp1];

         random_dbl8_exponential
            (deg,&rndhihihi,&rndlohihi,&rndhilohi,&rndlolohi,
                 &rndhihilo,&rndlohilo,&rndhilolo,&rndlololo,
             cffhihihi[i][j],cfflohihi[i][j],cffhilohi[i][j],cfflolohi[i][j],
             cffhihilo[i][j],cfflohilo[i][j],cffhilolo[i][j],cfflololo[i][j]);
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
               cout << cffhihihi[i][j][k] << "  "
                    << cfflohihi[i][j][k] << endl
                    << cffhilohi[i][j][k] << "  "
                    << cfflolohi[i][j][k] << endl
                    << cffhihilo[i][j][k] << "  "
                    << cfflohilo[i][j][k] << endl
                    << cffhilolo[i][j][k] << "  "
                    << cfflololo[i][j][k] << endl;
         }
      }
   }

/* define the input and allocate the output */

   double **inputhihihi = new double*[dim]; // dim series of degree deg
   double **inputlohihi = new double*[dim];
   double **inputhilohi = new double*[dim];
   double **inputlolohi = new double*[dim];
   double **inputhihilo = new double*[dim];
   double **inputlohilo = new double*[dim];
   double **inputhilolo = new double*[dim];
   double **inputlololo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      inputhihihi[i] = new double[degp1];
      inputlohihi[i] = new double[degp1];
      inputhilohi[i] = new double[degp1];
      inputlolohi[i] = new double[degp1];
      inputhihilo[i] = new double[degp1];
      inputlohilo[i] = new double[degp1];
      inputhilolo[i] = new double[degp1];
      inputlololo[i] = new double[degp1];
   }
   make_real8_input(dim,deg,inputhihihi,inputlohihi,inputhilohi,inputlolohi,
                            inputhihilo,inputlohilo,inputhilolo,inputlololo);

   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
      cout << "Random input series :" << endl;
      for(int i=0; i<dim; i++)
      {
         cout << "-> coefficients of series " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
            cout << inputhihihi[i][j] << "  " << inputlohihi[i][j] << endl
                 << inputhilohi[i][j] << "  " << inputlolohi[i][j] << endl
                 << inputhihilo[i][j] << "  " << inputlohilo[i][j] << endl
                 << inputhilolo[i][j] << "  " << inputlololo[i][j] << endl;
      }
   }
   double ***outputhihihi_h = new double**[dim];
   double ***outputlohihi_h = new double**[dim];
   double ***outputhilohi_h = new double**[dim];
   double ***outputlolohi_h = new double**[dim];
   double ***outputhihilo_h = new double**[dim];
   double ***outputlohilo_h = new double**[dim];
   double ***outputhilolo_h = new double**[dim];
   double ***outputlololo_h = new double**[dim];
   double ***outputhihihi_d = new double**[dim];
   double ***outputlohihi_d = new double**[dim];
   double ***outputhilohi_d = new double**[dim];
   double ***outputlolohi_d = new double**[dim];
   double ***outputhihilo_d = new double**[dim];
   double ***outputlohilo_d = new double**[dim];
   double ***outputhilolo_d = new double**[dim];
   double ***outputlololo_d = new double**[dim];

   for(int i=0; i<dim; i++)
   {
      outputhihihi_h[i] = new double*[dim+1];
      outputlohihi_h[i] = new double*[dim+1];
      outputhilohi_h[i] = new double*[dim+1];
      outputlolohi_h[i] = new double*[dim+1];
      outputhihilo_h[i] = new double*[dim+1];
      outputlohilo_h[i] = new double*[dim+1];
      outputhilolo_h[i] = new double*[dim+1];
      outputlololo_h[i] = new double*[dim+1];
      outputhihihi_d[i] = new double*[dim+1];
      outputlohihi_d[i] = new double*[dim+1];
      outputhilohi_d[i] = new double*[dim+1];
      outputlolohi_d[i] = new double*[dim+1];
      outputhihilo_d[i] = new double*[dim+1];
      outputlohilo_d[i] = new double*[dim+1];
      outputhilolo_d[i] = new double*[dim+1];
      outputlololo_d[i] = new double*[dim+1];

      for(int j=0; j<=dim; j++)
      {
         outputhihihi_h[i][j] = new double[degp1];
         outputlohihi_h[i][j] = new double[degp1];
         outputhilohi_h[i][j] = new double[degp1];
         outputlolohi_h[i][j] = new double[degp1];
         outputhihilo_h[i][j] = new double[degp1];
         outputlohilo_h[i][j] = new double[degp1];
         outputhilolo_h[i][j] = new double[degp1];
         outputlololo_h[i][j] = new double[degp1];
         outputhihihi_d[i][j] = new double[degp1];
         outputlohihi_d[i][j] = new double[degp1];
         outputhilohi_d[i][j] = new double[degp1];
         outputlolohi_d[i][j] = new double[degp1];
         outputhihilo_d[i][j] = new double[degp1];
         outputlohilo_d[i][j] = new double[degp1];
         outputhilolo_d[i][j] = new double[degp1];
         outputlololo_d[i][j] = new double[degp1];
      }
   }
   if(vrblvl > 0) cout << "computing on the host ..." << endl;

   double timelapsed_h = 0.0;

   for(int i=0; i<dim; i++)
   {
      double lapsed;

      CPU_dbl8_poly_evaldiff
        (dim,nbr[i],deg,nvr[i],idx[i],
         csthihihi[i],cstlohihi[i],csthilohi[i],cstlolohi[i],
         csthihilo[i],cstlohilo[i],csthilolo[i],cstlololo[i],
         cffhihihi[i],cfflohihi[i],cffhilohi[i],cfflolohi[i],
         cffhihilo[i],cfflohilo[i],cffhilolo[i],cfflololo[i],
         inputhihihi,inputlohihi,inputhilohi,inputlolohi,
         inputhihilo,inputlohilo,inputhilolo,inputlololo,
         outputhihihi_h[i],outputlohihi_h[i],
         outputhilohi_h[i],outputlolohi_h[i],
         outputhihilo_h[i],outputlohilo_h[i],
         outputhilolo_h[i],outputlololo_h[i],&lapsed,vrblvl);

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

      GPU_dbl8_poly_evaldiff
         (degp1,dim,nbr[i],deg,nvr[i],idx[i],
          csthihihi[i],cstlohihi[i],csthilohi[i],cstlolohi[i],
          csthihilo[i],cstlohilo[i],csthilolo[i],cstlololo[i],
          cffhihihi[i],cfflohihi[i],cffhilohi[i],cfflolohi[i],
          cffhihilo[i],cfflohilo[i],cffhilolo[i],cfflololo[i],
          inputhihihi,inputlohihi,inputhilohi,inputlolohi,
          inputhihilo,inputlohilo,inputhilolo,inputlololo,
          outputhihihi_d[i],outputlohihi_d[i],
          outputhilohi_d[i],outputlolohi_d[i],
          outputhihilo_d[i],outputlohilo_d[i],
          outputhilolo_d[i],outputlololo_d[i],
          cnvjobs,addjobs,&cnvlapms,&addlapms,&timelapms_d,
          &walltimes_d,vrblvl);

      timelapsed_d += walltimes_d;
   }
   if(vrblvl > 0) cout << "computing the errors ..." << endl;

   double sumerr = 0.0;

   for(int i=0; i<dim; i++)
   {
      sumerr += dbl8_error_sum1(dim,deg,
                   outputhihihi_h[i],outputlohihi_h[i],
                   outputhilohi_h[i],outputlolohi_h[i],
                   outputhihilo_h[i],outputlohilo_h[i],
                   outputhilolo_h[i],outputlololo_h[i],
                   outputhihihi_d[i],outputlohihi_d[i],
                   outputhilohi_d[i],outputlolohi_d[i],
                   outputhihilo_d[i],outputlohilo_d[i],
                   outputhilolo_d[i],outputlololo_d[i],vrb);
   }
   cout << scientific << setprecision(2)
        << "sum of all errors " << sumerr << endl;

   return sumerr;
}

double test_cmplx8_sysevaldiff
 ( int dim, int deg, int *nbr, int **nvr, int ***idx, int vrblvl )
{
   const int degp1 = deg+1;

/* generate constant and coefficients */

   double rndrehihihi,rndrelohihi,rndrehilohi,rndrelolohi;
   double rndrehihilo,rndrelohilo,rndrehilolo,rndrelololo;
   double rndimhihihi,rndimlohihi,rndimhilohi,rndimlolohi;
   double rndimhihilo,rndimlohilo,rndimhilolo,rndimlololo;
   double **cstrehihihi = new double*[dim];
   double **cstrelohihi = new double*[dim];
   double **cstrehilohi = new double*[dim];
   double **cstrelolohi = new double*[dim];
   double **cstrehihilo = new double*[dim];
   double **cstrelohilo = new double*[dim];
   double **cstrehilolo = new double*[dim];
   double **cstrelololo = new double*[dim];
   double **cstimhihihi = new double*[dim];
   double **cstimlohihi = new double*[dim];
   double **cstimhilohi = new double*[dim];
   double **cstimlolohi = new double*[dim];
   double **cstimhihilo = new double*[dim];
   double **cstimlohilo = new double*[dim];
   double **cstimhilolo = new double*[dim];
   double **cstimlololo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      cstrehihihi[i] = new double[degp1];
      cstrelohihi[i] = new double[degp1];
      cstrehilohi[i] = new double[degp1];
      cstrelolohi[i] = new double[degp1];
      cstrehihilo[i] = new double[degp1];
      cstrelohilo[i] = new double[degp1];
      cstrehilolo[i] = new double[degp1];
      cstrelololo[i] = new double[degp1];
      cstimhihihi[i] = new double[degp1];
      cstimlohihi[i] = new double[degp1];
      cstimhilohi[i] = new double[degp1];
      cstimlolohi[i] = new double[degp1];
      cstimhihilo[i] = new double[degp1];
      cstimlohilo[i] = new double[degp1];
      cstimhilolo[i] = new double[degp1];
      cstimlololo[i] = new double[degp1];
/*
      random_cmplx8_exponential
         (deg,&rndrehihihi,&rndrelohihi,&rndrehilohi,&rndrelolohi,
              &rndrehihilo,&rndrelohilo,&rndrehilolo,&rndrelololo,
              &rndimhihihi,&rndimlohihi,&rndimhilohi,&rndimlolohi,
              &rndimhihilo,&rndimlohilo,&rndimhilolo,&rndimlololo,
          cstrehihihi[i],cstrelohihi[i],cstrehilohi[i],cstrelolohi[i],
          cstrehihilo[i],cstrelohilo[i],cstrehilolo[i],cstrelololo[i],
          cstimhihihi[i],cstimlohihi[i],cstimhilohi[i],cstimlolohi[i],
          cstimhihilo[i],cstimlohilo[i],cstimhilolo[i],cstimlololo[i]);
 */
      for(int k=0; k<=deg; k++)
         random_octo_complex
            (&cstrehihihi[i][k],&cstrelohihi[i][k],
             &cstrehilohi[i][k],&cstrelolohi[i][k],
             &cstrehihilo[i][k],&cstrelohilo[i][k],
             &cstrehilolo[i][k],&cstrelololo[i][k],
             &cstimhihihi[i][k],&cstimlohihi[i][k],
             &cstimhilohi[i][k],&cstimlolohi[i][k],
             &cstimhihilo[i][k],&cstimlohilo[i][k],
             &cstimhilolo[i][k],&cstimlololo[i][k]);
   }
   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
      cout << "Constant coefficient series :" << endl;
      for(int i=0; i<dim; i++)
      {
         cout << "-> coefficients of series " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
            cout << cstrehihihi[i][j] << "  " << cstrelohihi[i][j] << endl
                 << cstrehilohi[i][j] << "  " << cstrelolohi[i][j] << endl
                 << cstrehihilo[i][j] << "  " << cstrelohilo[i][j] << endl
                 << cstrehilolo[i][j] << "  " << cstrelololo[i][j] << endl
                 << cstimhihihi[i][j] << "  " << cstimlohihi[i][j] << endl
                 << cstimhilohi[i][j] << "  " << cstimlolohi[i][j] << endl
                 << cstimhihilo[i][j] << "  " << cstimlohilo[i][j] << endl
                 << cstimhilolo[i][j] << "  " << cstimlololo[i][j] << endl;
      }
   }
   double ***cffrehihihi = new double**[dim];
   double ***cffrelohihi = new double**[dim];
   double ***cffrehilohi = new double**[dim];
   double ***cffrelolohi = new double**[dim];
   double ***cffrehihilo = new double**[dim];
   double ***cffrelohilo = new double**[dim];
   double ***cffrehilolo = new double**[dim];
   double ***cffrelololo = new double**[dim];
   double ***cffimhihihi = new double**[dim];
   double ***cffimlohihi = new double**[dim];
   double ***cffimhilohi = new double**[dim];
   double ***cffimlolohi = new double**[dim];
   double ***cffimhihilo = new double**[dim];
   double ***cffimlohilo = new double**[dim];
   double ***cffimhilolo = new double**[dim];
   double ***cffimlololo = new double**[dim];

   for(int i=0; i<dim; i++)
   {
      cffrehihihi[i] = new double*[nbr[i]];
      cffrelohihi[i] = new double*[nbr[i]];
      cffrehilohi[i] = new double*[nbr[i]];
      cffrelolohi[i] = new double*[nbr[i]];
      cffrehihilo[i] = new double*[nbr[i]];
      cffrelohilo[i] = new double*[nbr[i]];
      cffrehilolo[i] = new double*[nbr[i]];
      cffrelololo[i] = new double*[nbr[i]];
      cffimhihihi[i] = new double*[nbr[i]];
      cffimlohihi[i] = new double*[nbr[i]];
      cffimhilohi[i] = new double*[nbr[i]];
      cffimlolohi[i] = new double*[nbr[i]];
      cffimhihilo[i] = new double*[nbr[i]];
      cffimlohilo[i] = new double*[nbr[i]];
      cffimhilolo[i] = new double*[nbr[i]];
      cffimlololo[i] = new double*[nbr[i]];

      for(int j=0; j<nbr[i]; j++)
      {
         cffrehihihi[i][j] = new double[degp1];
         cffrelohihi[i][j] = new double[degp1];
         cffrehilohi[i][j] = new double[degp1];
         cffrelolohi[i][j] = new double[degp1];
         cffrehihilo[i][j] = new double[degp1];
         cffrelohilo[i][j] = new double[degp1];
         cffrehilolo[i][j] = new double[degp1];
         cffrelololo[i][j] = new double[degp1];
         cffimhihihi[i][j] = new double[degp1];
         cffimlohihi[i][j] = new double[degp1];
         cffimhilohi[i][j] = new double[degp1];
         cffimlolohi[i][j] = new double[degp1];
         cffimhihilo[i][j] = new double[degp1];
         cffimlohilo[i][j] = new double[degp1];
         cffimhilolo[i][j] = new double[degp1];
         cffimlololo[i][j] = new double[degp1];
/*
         random_cmplx8_exponential
            (deg,&rndrehihihi,&rndrelohihi,&rndrehilohi,&rndrelolohi,
                 &rndrehihilo,&rndrelohilo,&rndrehilolo,&rndrelololo,
                 &rndimhihihi,&rndimlohihi,&rndimhilohi,&rndimlolohi,
                 &rndimhihilo,&rndimlohilo,&rndimhilolo,&rndimlololo,
             cffrehihihi[i][j],cffrelohihi[i][j],
             cffrehilohi[i][j],cffrelolohi[i][j],
             cffrehihilo[i][j],cffrelohilo[i][j],
             cffrehilolo[i][j],cffrelololo[i][j],
             cffimhihihi[i][j],cffimlohihi[i][j],
             cffimhilohi[i][j],cffimlolohi[i][j],
             cffimhihilo[i][j],cffimlohilo[i][j],
             cffimhilolo[i][j],cffimlololo[i][j]);
 */
         for(int k=0; k<=deg; k++)
            random_octo_complex
               (&cffrehihihi[i][j][k],&cffrelohihi[i][j][k],
                &cffrehilohi[i][j][k],&cffrelolohi[i][j][k],
                &cffrehihilo[i][j][k],&cffrelohilo[i][j][k],
                &cffrehilolo[i][j][k],&cffrelololo[i][j][k],
                &cffimhihihi[i][j][k],&cffimlohihi[i][j][k],
                &cffimhilohi[i][j][k],&cffimlolohi[i][j][k],
                &cffimhihilo[i][j][k],&cffimlohilo[i][j][k],
                &cffimhilolo[i][j][k],&cffimlololo[i][j][k]);
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
               cout << cffrehihihi[i][j][k] << "  "
                    << cffrelohihi[i][j][k] << endl
                    << cffrehilohi[i][j][k] << "  "
                    << cffrelolohi[i][j][k] << endl
                    << cffrehihilo[i][j][k] << "  "
                    << cffrelohilo[i][j][k] << endl
                    << cffrehilolo[i][j][k] << "  "
                    << cffrelololo[i][j][k] << endl
                    << cffimhihihi[i][j][k] << "  "
                    << cffimlohihi[i][j][k] << endl
                    << cffimhilohi[i][j][k] << "  "
                    << cffimlolohi[i][j][k] << endl
                    << cffimhihilo[i][j][k] << "  "
                    << cffimlohilo[i][j][k] << endl
                    << cffimhilolo[i][j][k] << "  "
                    << cffimlololo[i][j][k] << endl;
         }
      }
   }

/* generate input series and allocate the output */

   double **inputrehihihi = new double*[dim]; // dim series of degree deg
   double **inputrelohihi = new double*[dim];
   double **inputrehilohi = new double*[dim];
   double **inputrelolohi = new double*[dim];
   double **inputrehihilo = new double*[dim];
   double **inputrelohilo = new double*[dim];
   double **inputrehilolo = new double*[dim];
   double **inputrelololo = new double*[dim];
   double **inputimhihihi = new double*[dim];
   double **inputimlohihi = new double*[dim];
   double **inputimhilohi = new double*[dim];
   double **inputimlolohi = new double*[dim];
   double **inputimhihilo = new double*[dim];
   double **inputimlohilo = new double*[dim];
   double **inputimhilolo = new double*[dim];
   double **inputimlololo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      inputrehihihi[i] = new double[degp1];
      inputrelohihi[i] = new double[degp1];
      inputrehilohi[i] = new double[degp1];
      inputrelolohi[i] = new double[degp1];
      inputrehihilo[i] = new double[degp1];
      inputrelohilo[i] = new double[degp1];
      inputrehilolo[i] = new double[degp1];
      inputrelololo[i] = new double[degp1];
      inputimhihihi[i] = new double[degp1];
      inputimlohihi[i] = new double[degp1];
      inputimhilohi[i] = new double[degp1];
      inputimlolohi[i] = new double[degp1];
      inputimhihilo[i] = new double[degp1];
      inputimlohilo[i] = new double[degp1];
      inputimhilolo[i] = new double[degp1];
      inputimlololo[i] = new double[degp1];
   }
   make_complex8_input
      (dim,deg,inputrehihihi,inputrelohihi,inputrehilohi,inputrelolohi,
               inputrehihilo,inputrelohilo,inputrehilolo,inputrelololo,
               inputimhihihi,inputimlohihi,inputimhilohi,inputimlolohi,
               inputimhihilo,inputimlohilo,inputimhilolo,inputimlololo);

   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
      cout << "Random input series :" << endl;
      for(int i=0; i<dim; i++)
      {
         cout << "-> coefficients of series " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
            cout << inputrehihihi[i][j] << "  " << inputrelohihi[i][j] << endl
                 << inputrehilohi[i][j] << "  " << inputrelolohi[i][j] << endl
                 << inputrehihilo[i][j] << "  " << inputrelohilo[i][j] << endl
                 << inputrehilolo[i][j] << "  " << inputrelololo[i][j] << endl
                 << inputimhihihi[i][j] << "  " << inputimlohihi[i][j] << endl
                 << inputimhilohi[i][j] << "  " << inputimlolohi[i][j] << endl
                 << inputimhihilo[i][j] << "  " << inputimlohilo[i][j] << endl
                 << inputimhilolo[i][j] << "  " << inputimlololo[i][j] << endl;
      }
   }
   double ***outputrehihihi_h = new double**[dim];
   double ***outputrelohihi_h = new double**[dim];
   double ***outputrehilohi_h = new double**[dim];
   double ***outputrelolohi_h = new double**[dim];
   double ***outputrehihilo_h = new double**[dim];
   double ***outputrelohilo_h = new double**[dim];
   double ***outputrehilolo_h = new double**[dim];
   double ***outputrelololo_h = new double**[dim];
   double ***outputimhihihi_h = new double**[dim];
   double ***outputimlohihi_h = new double**[dim];
   double ***outputimhilohi_h = new double**[dim];
   double ***outputimlolohi_h = new double**[dim];
   double ***outputimhihilo_h = new double**[dim];
   double ***outputimlohilo_h = new double**[dim];
   double ***outputimhilolo_h = new double**[dim];
   double ***outputimlololo_h = new double**[dim];
   double ***outputrehihihi_d = new double**[dim];
   double ***outputrelohihi_d = new double**[dim];
   double ***outputrehilohi_d = new double**[dim];
   double ***outputrelolohi_d = new double**[dim];
   double ***outputrehihilo_d = new double**[dim];
   double ***outputrelohilo_d = new double**[dim];
   double ***outputrehilolo_d = new double**[dim];
   double ***outputrelololo_d = new double**[dim];
   double ***outputimhihihi_d = new double**[dim];
   double ***outputimlohihi_d = new double**[dim];
   double ***outputimhilohi_d = new double**[dim];
   double ***outputimlolohi_d = new double**[dim];
   double ***outputimhihilo_d = new double**[dim];
   double ***outputimlohilo_d = new double**[dim];
   double ***outputimhilolo_d = new double**[dim];
   double ***outputimlololo_d = new double**[dim];

   for(int i=0; i<dim; i++)
   {
      outputrehihihi_h[i] = new double*[dim+1];
      outputrelohihi_h[i] = new double*[dim+1];
      outputrehilohi_h[i] = new double*[dim+1];
      outputrelolohi_h[i] = new double*[dim+1];
      outputrehihilo_h[i] = new double*[dim+1];
      outputrelohilo_h[i] = new double*[dim+1];
      outputrehilolo_h[i] = new double*[dim+1];
      outputrelololo_h[i] = new double*[dim+1];
      outputimhihihi_h[i] = new double*[dim+1];
      outputimlohihi_h[i] = new double*[dim+1];
      outputimhilohi_h[i] = new double*[dim+1];
      outputimlolohi_h[i] = new double*[dim+1];
      outputimhihilo_h[i] = new double*[dim+1];
      outputimlohilo_h[i] = new double*[dim+1];
      outputimhilolo_h[i] = new double*[dim+1];
      outputimlololo_h[i] = new double*[dim+1];
      outputrehihihi_d[i] = new double*[dim+1];
      outputrelohihi_d[i] = new double*[dim+1];
      outputrehilohi_d[i] = new double*[dim+1];
      outputrelolohi_d[i] = new double*[dim+1];
      outputrehihilo_d[i] = new double*[dim+1];
      outputrelohilo_d[i] = new double*[dim+1];
      outputrehilolo_d[i] = new double*[dim+1];
      outputrelololo_d[i] = new double*[dim+1];
      outputimhihihi_d[i] = new double*[dim+1];
      outputimlohihi_d[i] = new double*[dim+1];
      outputimhilohi_d[i] = new double*[dim+1];
      outputimlolohi_d[i] = new double*[dim+1];
      outputimhihilo_d[i] = new double*[dim+1];
      outputimlohilo_d[i] = new double*[dim+1];
      outputimhilolo_d[i] = new double*[dim+1];
      outputimlololo_d[i] = new double*[dim+1];

      for(int j=0; j<=dim; j++)
      {
         outputrehihihi_h[i][j] = new double[degp1];
         outputrelohihi_h[i][j] = new double[degp1];
         outputrehilohi_h[i][j] = new double[degp1];
         outputrelolohi_h[i][j] = new double[degp1];
         outputrehihilo_h[i][j] = new double[degp1];
         outputrelohilo_h[i][j] = new double[degp1];
         outputrehilolo_h[i][j] = new double[degp1];
         outputrelololo_h[i][j] = new double[degp1];
         outputimhihihi_h[i][j] = new double[degp1];
         outputimlohihi_h[i][j] = new double[degp1];
         outputimhilohi_h[i][j] = new double[degp1];
         outputimlolohi_h[i][j] = new double[degp1];
         outputimhihilo_h[i][j] = new double[degp1];
         outputimlohilo_h[i][j] = new double[degp1];
         outputimhilolo_h[i][j] = new double[degp1];
         outputimlololo_h[i][j] = new double[degp1];
         outputrehihihi_d[i][j] = new double[degp1];
         outputrelohihi_d[i][j] = new double[degp1];
         outputrehilohi_d[i][j] = new double[degp1];
         outputrelolohi_d[i][j] = new double[degp1];
         outputrehihilo_d[i][j] = new double[degp1];
         outputrelohilo_d[i][j] = new double[degp1];
         outputrehilolo_d[i][j] = new double[degp1];
         outputrelololo_d[i][j] = new double[degp1];
         outputimhihihi_d[i][j] = new double[degp1];
         outputimlohihi_d[i][j] = new double[degp1];
         outputimhilohi_d[i][j] = new double[degp1];
         outputimlolohi_d[i][j] = new double[degp1];
         outputimhihilo_d[i][j] = new double[degp1];
         outputimlohilo_d[i][j] = new double[degp1];
         outputimhilolo_d[i][j] = new double[degp1];
         outputimlololo_d[i][j] = new double[degp1];
      }
   }
   if(vrblvl > 0) cout << "evaluating on the host ..." << endl;

   double timelapsed_h = 0.0;

   for(int i=0; i<dim; i++)
   {
      double lapsed;

      CPU_cmplx8_poly_evaldiff
        (dim,nbr[i],deg,nvr[i],idx[i],
         cstrehihihi[i],cstrelohihi[i],cstrehilohi[i],cstrelolohi[i],
         cstrehihilo[i],cstrelohilo[i],cstrehilolo[i],cstrelololo[i],
         cstimhihihi[i],cstimlohihi[i],cstimhilohi[i],cstimlolohi[i],
         cstimhihilo[i],cstimlohilo[i],cstimhilolo[i],cstimlololo[i],
         cffrehihihi[i],cffrelohihi[i],cffrehilohi[i],cffrelolohi[i],
         cffrehihilo[i],cffrelohilo[i],cffrehilolo[i],cffrelololo[i],
         cffimhihihi[i],cffimlohihi[i],cffimhilohi[i],cffimlolohi[i],
         cffimhihilo[i],cffimlohilo[i],cffimhilolo[i],cffimlololo[i],
         inputrehihihi,inputrelohihi,inputrehilohi,inputrelolohi,
         inputrehihilo,inputrelohilo,inputrehilolo,inputrelololo,
         inputimhihihi,inputimlohihi,inputimhilohi,inputimlolohi,
         inputimhihilo,inputimlohilo,inputimhilolo,inputimlololo,
         outputrehihihi_h[i],outputrelohihi_h[i],
         outputrehilohi_h[i],outputrelolohi_h[i],
         outputrehihilo_h[i],outputrelohilo_h[i],
         outputrehilolo_h[i],outputrelololo_h[i],
         outputimhihihi_h[i],outputimlohihi_h[i],
         outputimhilohi_h[i],outputimlolohi_h[i],
         outputimhihilo_h[i],outputimlohilo_h[i],
         outputimhilolo_h[i],outputimlololo_h[i],&lapsed,vrblvl);

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

      GPU_cmplx8vectorized_poly_evaldiff
         (degp1,dim,nbr[i],deg,nvr[i],idx[i],
          cstrehihihi[i],cstrelohihi[i],cstrehilohi[i],cstrelolohi[i],
          cstrehihilo[i],cstrelohilo[i],cstrehilolo[i],cstrelololo[i],
          cstimhihihi[i],cstimlohihi[i],cstimhilohi[i],cstimlolohi[i],
          cstimhihilo[i],cstimlohilo[i],cstimhilolo[i],cstimlololo[i],
          cffrehihihi[i],cffrelohihi[i],cffrehilohi[i],cffrelolohi[i],
          cffrehihilo[i],cffrelohilo[i],cffrehilolo[i],cffrelololo[i],
          cffimhihihi[i],cffimlohihi[i],cffimhilohi[i],cffimlolohi[i],
          cffimhihilo[i],cffimlohilo[i],cffimhilolo[i],cffimlololo[i],
          inputrehihihi,inputrelohihi,inputrehilohi,inputrelolohi,
          inputrehihilo,inputrelohilo,inputrehilolo,inputrelololo,
          inputimhihihi,inputimlohihi,inputimhilohi,inputimlolohi,
          inputimhihilo,inputimlohilo,inputimhilolo,inputimlololo,
          outputrehihihi_d[i],outputrelohihi_d[i],
          outputrehilohi_d[i],outputrelolohi_d[i],
          outputrehihilo_d[i],outputrelohilo_d[i],
          outputrehilolo_d[i],outputrelololo_d[i],
          outputimhihihi_d[i],outputimlohihi_d[i],
          outputimhilohi_d[i],outputimlolohi_d[i],
          outputimhihilo_d[i],outputimlohilo_d[i],
          outputimhilolo_d[i],outputimlololo_d[i],
          cnvjobs,incjobs,addjobs,&cnvlapms,&addlapms,&timelapms_d,
          &walltimes_d,vrblvl);

      timelapsed_d += walltimes_d;
   }
   if(vrblvl > 0) cout << "computing the errors ..." << endl;

   double sumerr = 0.0;

   for(int i=0; i<dim; i++)
   {
      sumerr += cmplx8_error_sum(dim,deg,
                   outputrehihihi_h[i],outputrelohihi_h[i],
                   outputrehilohi_h[i],outputrelolohi_h[i],
                   outputrehihilo_h[i],outputrelohilo_h[i],
                   outputrehilolo_h[i],outputrelololo_h[i],
                   outputimhihihi_h[i],outputimlohihi_h[i],
                   outputimhilohi_h[i],outputimlolohi_h[i],
                   outputimhihilo_h[i],outputimlohilo_h[i],
                   outputimhilolo_h[i],outputimlololo_h[i],
                   outputrehihihi_d[i],outputrelohihi_d[i],
                   outputrehilohi_d[i],outputrelolohi_d[i],
                   outputrehihilo_d[i],outputrelohilo_d[i],
                   outputrehilolo_d[i],outputrelololo_d[i],
                   outputimhihihi_d[i],outputimlohihi_d[i],
                   outputimhilohi_d[i],outputimlolohi_d[i],
                   outputimhihilo_d[i],outputimlohilo_d[i],
                   outputimhilolo_d[i],outputimlololo_d[i],vrb);
   }
   cout << scientific << setprecision(2)
        << "sum of all errors " << sumerr << endl;

   return sumerr;
}
