/* Tests operations on Newton for series in octo double precision */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <time.h>
#include "unimodular_matrices.h"
#include "dbl8_newton_testers.h"

using namespace std;

int main ( void )
{
   srand(time(NULL));

   cout << "testing Newton ..." << endl;

   int dim,deg,size,posvals,vrblvl,nbritr,nbsteps;
   prompt_dimensions(&dim,&deg,&size,&posvals,&vrblvl,&nbritr,&nbsteps);

   cout << "generating an integer matrix of dimension " << dim
        << " ..." << endl;
/*
 * 1. defining the unimodular matrix of monomials
 */
   int **rowsA = new int*[dim];  // exponents in the rows
   int *nvr = new int[dim];      // number of variables in each monomial
   int **idx = new int*[dim];    // indexes of variables in each monomial
   int **exp = new int*[dim];    // exponents of the variables
   int *nbrfac = new int[dim];   // number of exponents > 1 in each monomial
   int **expfac = new int*[dim]; // exponents of the common factors

   make_monomial_system
      (dim,size,posvals,nbritr,nvr,idx,exp,nbrfac,expfac,rowsA,vrblvl);
/*
 * 2. allocating input and output space for evaluation and differentiation
 */

   const int degp1 = deg+1;
   double **inputhihihi = new double*[dim];
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
   // allocate memory for coefficients and the output
   double *acchihihi = new double[degp1]; // accumulated power series
   double *acclohihi = new double[degp1];
   double *acchilohi = new double[degp1];
   double *acclolohi = new double[degp1];
   double *acchihilo = new double[degp1];
   double *acclohilo = new double[degp1];
   double *acchilolo = new double[degp1];
   double *acclololo = new double[degp1];
   double **cffhihihi = new double*[dim]; // the coefficients of monomials
   double **cfflohihi = new double*[dim];
   double **cffhilohi = new double*[dim];
   double **cfflolohi = new double*[dim];
   double **cffhihilo = new double*[dim];
   double **cfflohilo = new double*[dim];
   double **cffhilolo = new double*[dim];
   double **cfflololo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      cffhihihi[i] = new double[degp1];
      cfflohihi[i] = new double[degp1];
      cffhilohi[i] = new double[degp1];
      cfflolohi[i] = new double[degp1];
      cffhihilo[i] = new double[degp1];
      cfflohilo[i] = new double[degp1];
      cffhilolo[i] = new double[degp1];
      cfflololo[i] = new double[degp1];
   }
   double ***outputhihihi = new double**[dim];
   double ***outputlohihi = new double**[dim];
   double ***outputhilohi = new double**[dim];
   double ***outputlolohi = new double**[dim];
   double ***outputhihilo = new double**[dim];
   double ***outputlohilo = new double**[dim];
   double ***outputhilolo = new double**[dim];
   double ***outputlololo = new double**[dim];

   for(int i=0; i<dim; i++)
   {
      outputhihihi[i] = new double*[dim+1];
      outputlohihi[i] = new double*[dim+1];
      outputhilohi[i] = new double*[dim+1];
      outputlolohi[i] = new double*[dim+1];
      outputhihilo[i] = new double*[dim+1];
      outputlohilo[i] = new double*[dim+1];
      outputhilolo[i] = new double*[dim+1];
      outputlololo[i] = new double*[dim+1];

      for(int j=0; j<=dim; j++)
      {
         outputhihihi[i][j] = new double[degp1];
         outputlohihi[i][j] = new double[degp1];
         outputhilohi[i][j] = new double[degp1];
         outputlolohi[i][j] = new double[degp1];
         outputhihilo[i][j] = new double[degp1];
         outputlohilo[i][j] = new double[degp1];
         outputhilolo[i][j] = new double[degp1];
         outputlololo[i][j] = new double[degp1];
      }
   }
   // The function values are power series truncated at degree deg.
   double **funvalhihihi = new double*[dim];
   double **funvallohihi = new double*[dim];
   double **funvalhilohi = new double*[dim];
   double **funvallolohi = new double*[dim];
   double **funvalhihilo = new double*[dim];
   double **funvallohilo = new double*[dim];
   double **funvalhilolo = new double*[dim];
   double **funvallololo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      funvalhihihi[i] = new double[degp1];
      funvallohihi[i] = new double[degp1];
      funvalhilohi[i] = new double[degp1];
      funvallolohi[i] = new double[degp1];
      funvalhihilo[i] = new double[degp1];
      funvallohilo[i] = new double[degp1];
      funvalhilolo[i] = new double[degp1];
      funvallololo[i] = new double[degp1];
   }
   // The derivatives in the output are a series truncated at degree deg.
   // The coefficients of the series are matrices of dimension dim.
   double ***jacvalhihihi = new double**[degp1];
   double ***jacvallohihi = new double**[degp1];
   double ***jacvalhilohi = new double**[degp1];
   double ***jacvallolohi = new double**[degp1];
   double ***jacvalhihilo = new double**[degp1];
   double ***jacvallohilo = new double**[degp1];
   double ***jacvalhilolo = new double**[degp1];
   double ***jacvallololo = new double**[degp1];

   for(int i=0; i<degp1; i++) // jacval[i] is matrix of dimension dim
   {
      jacvalhihihi[i] = new double*[dim];
      jacvallohihi[i] = new double*[dim];
      jacvalhilohi[i] = new double*[dim];
      jacvallolohi[i] = new double*[dim];
      jacvalhihilo[i] = new double*[dim];
      jacvallohilo[i] = new double*[dim];
      jacvalhilolo[i] = new double*[dim];
      jacvallololo[i] = new double*[dim];

      for(int j=0; j<dim; j++)
      {
         jacvalhihihi[i][j] = new double[dim];
         jacvallohihi[i][j] = new double[dim];
         jacvalhilohi[i][j] = new double[dim];
         jacvallolohi[i][j] = new double[dim];
         jacvalhihilo[i][j] = new double[dim];
         jacvallohilo[i][j] = new double[dim];
         jacvalhilolo[i][j] = new double[dim];
         jacvallololo[i][j] = new double[dim];
      }
   }
/*
 * 3. allocate space to solve the linearized power series system
 */
   // The solution x(t) to jacval(t)*x(t) = -funval(t) in linearized
   // format is a series truncated at degree deg, with as coefficients
   // arrays of dimension dim.
   double **solhihihi = new double*[degp1];
   double **sollohihi = new double*[degp1];
   double **solhilohi = new double*[degp1];
   double **sollolohi = new double*[degp1];
   double **solhihilo = new double*[degp1];
   double **sollohilo = new double*[degp1];
   double **solhilolo = new double*[degp1];
   double **sollololo = new double*[degp1];

   for(int i=0; i<degp1; i++)
   {
      solhihihi[i] = new double[dim];
      sollohihi[i] = new double[dim];
      solhilohi[i] = new double[dim];
      sollolohi[i] = new double[dim];
      solhihilo[i] = new double[dim];
      sollohilo[i] = new double[dim];
      solhilolo[i] = new double[dim];
      sollololo[i] = new double[dim];
   }
   // The right hand side -funval(t) in linearized format is a series
   // truncated at degree deg, with arrays of dimension dim as coefficients.
   double **rhshihihi = new double*[degp1];
   double **rhslohihi = new double*[degp1];
   double **rhshilohi = new double*[degp1];
   double **rhslolohi = new double*[degp1];
   double **rhshihilo = new double*[degp1];
   double **rhslohilo = new double*[degp1];
   double **rhshilolo = new double*[degp1];
   double **rhslololo = new double*[degp1];

   for(int i=0; i<degp1; i++)
   {
      rhshihihi[i] = new double[dim];
      rhslohihi[i] = new double[dim];
      rhshilohi[i] = new double[dim];
      rhslolohi[i] = new double[dim];
      rhshihilo[i] = new double[dim];
      rhslohilo[i] = new double[dim];
      rhshilolo[i] = new double[dim];
      rhslololo[i] = new double[dim];
   }
   // Allocate work space for the inplace LU solver.
   double **workmathihihi = new double*[dim];
   double **workmatlohihi = new double*[dim];
   double **workmathilohi = new double*[dim];
   double **workmatlolohi = new double*[dim];
   double **workmathihilo = new double*[dim];
   double **workmatlohilo = new double*[dim];
   double **workmathilolo = new double*[dim];
   double **workmatlololo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      workmathihihi[i] = new double[dim];
      workmatlohihi[i] = new double[dim];
      workmathilohi[i] = new double[dim];
      workmatlolohi[i] = new double[dim];
      workmathihilo[i] = new double[dim];
      workmatlohilo[i] = new double[dim];
      workmathilolo[i] = new double[dim];
      workmatlololo[i] = new double[dim];
   }
   int *ipvt = new int[dim];
   double *workvechihihi = new double[dim];
   double *workveclohihi = new double[dim];
   double *workvechilohi = new double[dim];
   double *workveclolohi = new double[dim];
   double *workvechihilo = new double[dim];
   double *workveclohilo = new double[dim];
   double *workvechilolo = new double[dim];
   double *workveclololo = new double[dim];
   // Copy the rhs vector into work space for inplace solver.
   double **workrhshihihi = new double*[degp1];
   double **workrhslohihi = new double*[degp1];
   double **workrhshilohi = new double*[degp1];
   double **workrhslolohi = new double*[degp1];
   double **workrhshihilo = new double*[degp1];
   double **workrhslohilo = new double*[degp1];
   double **workrhshilolo = new double*[degp1];
   double **workrhslololo = new double*[degp1];

   for(int i=0; i<degp1; i++)
   {
      workrhshihihi[i] = new double[dim];
      workrhslohihi[i] = new double[dim];
      workrhshilohi[i] = new double[dim];
      workrhslolohi[i] = new double[dim];
      workrhshihilo[i] = new double[dim];
      workrhslohilo[i] = new double[dim];
      workrhshilolo[i] = new double[dim];
      workrhslololo[i] = new double[dim];
   }
   double **resvechihihi = new double*[degp1];
   double **resveclohihi = new double*[degp1];
   double **resvechilohi = new double*[degp1];
   double **resveclolohi = new double*[degp1];
   double **resvechihilo = new double*[degp1];
   double **resveclohilo = new double*[degp1];
   double **resvechilolo = new double*[degp1];
   double **resveclololo = new double*[degp1];

   for(int i=0; i<degp1; i++)
   {
      resvechihihi[i] = new double[dim];
      resveclohihi[i] = new double[dim];
      resvechilohi[i] = new double[dim];
      resveclolohi[i] = new double[dim];
      resvechihilo[i] = new double[dim];
      resveclohilo[i] = new double[dim];
      resvechilolo[i] = new double[dim];
      resveclololo[i] = new double[dim];
   }
   double resmaxhihihi,resmaxlohihi,resmaxhilohi,resmaxlolohi;
   double resmaxhihilo,resmaxlohilo,resmaxhilolo,resmaxlololo;

/*
 * 4. initialize input, coefficient, evaluate, differentiate, and solve
 */
   // Define the initial input, a vector of ones.
   dbl8_unit_series_vector
      (dim,deg,inputhihihi,inputlohihi,inputhilohi,inputlolohi,
               inputhihilo,inputlohilo,inputhilolo,inputlololo);

   if(vrblvl > 0)
   {
      cout << scientific << setprecision(16);
      cout << "The leading coefficients of the input series :" << endl;
      for(int i=0; i<dim; i++)
      {
         cout << i << " : " << inputhihihi[i][0] << "  "
                            << inputlohihi[i][0] << endl;
         cout << "     " << inputhilohi[i][0] << "  "
                         << inputlolohi[i][0] << endl;
         cout << "     " << inputhihilo[i][0] << "  "
                         << inputlohilo[i][0] << endl;
         cout << "     " << inputhilolo[i][0] << "  "
                         << inputlololo[i][0] << endl;
      }
   }
   for(int step=0; step<nbsteps; step++)
   {
      cout << "step " << step << " ..." << endl;

      dbl8_newton_step
         (dim,deg,nvr,idx,exp,nbrfac,expfac,
          cffhihihi,cfflohihi,cffhilohi,cfflolohi,
          cffhihilo,cfflohilo,cffhilolo,cfflololo,
          acchihihi,acclohihi,acchilohi,acclolohi,
          acchihilo,acclohilo,acchilolo,acclololo,
           inputhihihi, inputlohihi, inputhilohi, inputlolohi,
           inputhihilo, inputlohilo, inputhilolo, inputlololo,
          outputhihihi,outputlohihi,outputhilohi,outputlolohi,
          outputhihilo,outputlohilo,outputhilolo,outputlololo,
          funvalhihihi,funvallohihi,funvalhilohi,funvallolohi,
          funvalhihilo,funvallohilo,funvalhilolo,funvallololo,
          jacvalhihihi,jacvallohihi,jacvalhilohi,jacvallolohi,
          jacvalhihilo,jacvallohilo,jacvalhilolo,jacvallololo,
          rhshihihi,rhslohihi,rhshilohi,rhslolohi,
          rhshihilo,rhslohilo,rhshilolo,rhslololo,
          solhihihi,sollohihi,solhilohi,sollolohi,
          solhihilo,sollohilo,solhilolo,sollololo,
          workmathihihi,workmatlohihi,workmathilohi,workmatlolohi,
          workmathihilo,workmatlohilo,workmathilolo,workmatlololo,
          workvechihihi,workveclohihi,workvechilohi,workveclolohi,
          workvechihilo,workveclohilo,workvechilolo,workveclololo,
          workrhshihihi,workrhslohihi,workrhshilohi,workrhslolohi,
          workrhshihilo,workrhslohilo,workrhshilolo,workrhslololo,
          resvechihihi,resveclohihi,resvechilohi,resveclolohi,
          resvechihilo,resveclohilo,resvechilolo,resveclololo,
          &resmaxhihihi,&resmaxlohihi,&resmaxhilohi,&resmaxlolohi,
          &resmaxhihilo,&resmaxlohilo,&resmaxhilolo,&resmaxlololo,
          ipvt,vrblvl);
   }
   return 0;
}
