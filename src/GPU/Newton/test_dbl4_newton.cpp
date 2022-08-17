/* Tests operations on Newton for series in quad double precision */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <time.h>
#include "unimodular_matrices.h"
#include "dbl4_newton_testers.h"

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
   double **inputhihi = new double*[dim];
   double **inputlohi = new double*[dim];
   double **inputhilo = new double*[dim];
   double **inputlolo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      inputhihi[i] = new double[degp1];
      inputlohi[i] = new double[degp1];
      inputhilo[i] = new double[degp1];
      inputlolo[i] = new double[degp1];
   }
   // allocate memory for coefficients and the output
   double *acchihi = new double[degp1]; // accumulated power series
   double *acclohi = new double[degp1];
   double *acchilo = new double[degp1];
   double *acclolo = new double[degp1];
   double **cffhihi = new double*[dim]; // the coefficients of monomials
   double **cfflohi = new double*[dim];
   double **cffhilo = new double*[dim];
   double **cfflolo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      cffhihi[i] = new double[degp1];
      cfflohi[i] = new double[degp1];
      cffhilo[i] = new double[degp1];
      cfflolo[i] = new double[degp1];
   }
   double ***outputhihi = new double**[dim];
   double ***outputlohi = new double**[dim];
   double ***outputhilo = new double**[dim];
   double ***outputlolo = new double**[dim];

   for(int i=0; i<dim; i++)
   {
      outputhihi[i] = new double*[dim+1];
      outputlohi[i] = new double*[dim+1];
      outputhilo[i] = new double*[dim+1];
      outputlolo[i] = new double*[dim+1];

      for(int j=0; j<=dim; j++)
      {
         outputhihi[i][j] = new double[degp1];
         outputlohi[i][j] = new double[degp1];
         outputhilo[i][j] = new double[degp1];
         outputlolo[i][j] = new double[degp1];
      }
   }
   // The function values are power series truncated at degree deg.
   double **funvalhihi = new double*[dim];
   double **funvallohi = new double*[dim];
   double **funvalhilo = new double*[dim];
   double **funvallolo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      funvalhihi[i] = new double[degp1];
      funvallohi[i] = new double[degp1];
      funvalhilo[i] = new double[degp1];
      funvallolo[i] = new double[degp1];
   }
   // The derivatives in the output are a series truncated at degree deg.
   // The coefficients of the series are matrices of dimension dim.
   double ***jacvalhihi = new double**[degp1];
   double ***jacvallohi = new double**[degp1];
   double ***jacvalhilo = new double**[degp1];
   double ***jacvallolo = new double**[degp1];

   for(int i=0; i<degp1; i++) // jacval[i] is matrix of dimension dim
   {
      jacvalhihi[i] = new double*[dim];
      jacvallohi[i] = new double*[dim];
      jacvalhilo[i] = new double*[dim];
      jacvallolo[i] = new double*[dim];

      for(int j=0; j<dim; j++)
      {
         jacvalhihi[i][j] = new double[dim];
         jacvallohi[i][j] = new double[dim];
         jacvalhilo[i][j] = new double[dim];
         jacvallolo[i][j] = new double[dim];
      }
   }
/*
 * 3. allocate space to solve the linearized power series system
 */
   // The solution x(t) to jacval(t)*x(t) = -funval(t) in linearized
   // format is a series truncated at degree deg, with as coefficients
   // arrays of dimension dim.
   double **solhihi = new double*[degp1];
   double **sollohi = new double*[degp1];
   double **solhilo = new double*[degp1];
   double **sollolo = new double*[degp1];

   for(int i=0; i<degp1; i++)
   {
      solhihi[i] = new double[dim];
      sollohi[i] = new double[dim];
      solhilo[i] = new double[dim];
      sollolo[i] = new double[dim];
   }
   // The right hand side -funval(t) in linearized format is a series
   // truncated at degree deg, with arrays of dimension dim as coefficients.
   double **rhshihi = new double*[degp1];
   double **rhslohi = new double*[degp1];
   double **rhshilo = new double*[degp1];
   double **rhslolo = new double*[degp1];

   for(int i=0; i<degp1; i++)
   {
      rhshihi[i] = new double[dim];
      rhslohi[i] = new double[dim];
      rhshilo[i] = new double[dim];
      rhslolo[i] = new double[dim];
   }
   // Allocate work space for the inplace LU solver.
   double **workmathihi = new double*[dim];
   double **workmatlohi = new double*[dim];
   double **workmathilo = new double*[dim];
   double **workmatlolo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      workmathihi[i] = new double[dim];
      workmatlohi[i] = new double[dim];
      workmathilo[i] = new double[dim];
      workmatlolo[i] = new double[dim];
   }
   int *ipvt = new int[dim];
   double *workvechihi = new double[dim];
   double *workveclohi = new double[dim];
   double *workvechilo = new double[dim];
   double *workveclolo = new double[dim];
   // Copy the rhs vector into work space for inplace solver.
   double **workrhshihi = new double*[degp1];
   double **workrhslohi = new double*[degp1];
   double **workrhshilo = new double*[degp1];
   double **workrhslolo = new double*[degp1];

   for(int i=0; i<degp1; i++)
   {
      workrhshihi[i] = new double[dim];
      workrhslohi[i] = new double[dim];
      workrhshilo[i] = new double[dim];
      workrhslolo[i] = new double[dim];
   }
   double **resvechihi = new double*[degp1];
   double **resveclohi = new double*[degp1];
   double **resvechilo = new double*[degp1];
   double **resveclolo = new double*[degp1];

   for(int i=0; i<degp1; i++)
   {
      resvechihi[i] = new double[dim];
      resveclohi[i] = new double[dim];
      resvechilo[i] = new double[dim];
      resveclolo[i] = new double[dim];
   }
   double resmaxhihi,resmaxlohi,resmaxhilo,resmaxlolo;

/*
 * 4. initialize input, coefficient, evaluate, differentiate, and solve
 */
   // Define the initial input, a vector of ones.
   dbl4_unit_series_vector
      (dim,deg,inputhihi,inputlohi,inputhilo,inputlolo);

   if(vrblvl > 0)
   {
      cout << scientific << setprecision(16);
      cout << "The leading coefficients of the input series :" << endl;
      for(int i=0; i<dim; i++)
      {
         cout << i << " : " << inputhihi[i][0] << "  "
                            << inputlohi[i][0] << endl;
         cout << "     " << inputhilo[i][0] << "  "
                         << inputlolo[i][0] << endl;
      }
   }
   for(int step=0; step<nbsteps; step++)
   {
      cout << "step " << step << " ..." << endl;

      dbl4_newton_step
         (dim,deg,nvr,idx,exp,nbrfac,expfac,
          cffhihi,cfflohi,cffhilo,cfflolo,acchihi,acclohi,acchilo,acclolo,
           inputhihi, inputlohi, inputhilo, inputlolo,
          outputhihi,outputlohi,outputhilo,outputlolo,
          funvalhihi,funvallohi,funvalhilo,funvallolo,
          jacvalhihi,jacvallohi,jacvalhilo,jacvallolo,
          rhshihi,rhslohi,rhshilo,rhslolo,solhihi,sollohi,solhilo,sollolo,
          workmathihi,workmatlohi,workmathilo,workmatlolo,
          workvechihi,workveclohi,workvechilo,workveclolo,
          workrhshihi,workrhslohi,workrhshilo,workrhslolo,
          resvechihi,resveclohi,resvechilo,resveclolo,
          &resmaxhihi,&resmaxlohi,&resmaxhilo,&resmaxlolo,ipvt,vrblvl);
   }
   return 0;
}
