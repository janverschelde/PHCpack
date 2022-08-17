/* Tests operations on Newton for series in double double precision */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <time.h>
#include "unimodular_matrices.h"
#include "dbl2_newton_testers.h"

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
   double **inputhi = new double*[dim];
   double **inputlo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      inputhi[i] = new double[degp1];
      inputlo[i] = new double[degp1];
   }
   // allocate memory for coefficients and the output
   double *acchi = new double[degp1]; // accumulated power series
   double *acclo = new double[degp1];
   double **cffhi = new double*[dim]; // the coefficients of monomials
   double **cfflo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      cffhi[i] = new double[degp1];
      cfflo[i] = new double[degp1];
   }
   double ***outputhi = new double**[dim];
   double ***outputlo = new double**[dim];

   for(int i=0; i<dim; i++)
   {
      outputhi[i] = new double*[dim+1];
      outputlo[i] = new double*[dim+1];

      for(int j=0; j<=dim; j++)
      {
         outputhi[i][j] = new double[degp1];
         outputlo[i][j] = new double[degp1];
      }
   }
   // The function values are power series truncated at degree deg.
   double **funvalhi = new double*[dim];
   double **funvallo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      funvalhi[i] = new double[degp1];
      funvallo[i] = new double[degp1];
   }
   // The derivatives in the output are a series truncated at degree deg.
   // The coefficients of the series are matrices of dimension dim.
   double ***jacvalhi = new double**[degp1];
   double ***jacvallo = new double**[degp1];

   for(int i=0; i<degp1; i++) // jacval[i] is matrix of dimension dim
   {
      jacvalhi[i] = new double*[dim];
      jacvallo[i] = new double*[dim];

      for(int j=0; j<dim; j++)
      {
         jacvalhi[i][j] = new double[dim];
         jacvallo[i][j] = new double[dim];
      }
   }
/*
 * 3. allocate space to solve the linearized power series system
 */
   // The solution x(t) to jacval(t)*x(t) = -funval(t) in linearized
   // format is a series truncated at degree deg, with as coefficients
   // arrays of dimension dim.

   double **solhi = new double*[degp1];
   double **sollo = new double*[degp1];

   for(int i=0; i<degp1; i++)
   {
      solhi[i] = new double[dim];
      sollo[i] = new double[dim];
   }
   // The right hand side -funval(t) in linearized format is a series
   // truncated at degree deg, with arrays of dimension dim as coefficients.
   double **rhshi = new double*[degp1];
   double **rhslo = new double*[degp1];

   for(int i=0; i<degp1; i++)
   {
      rhshi[i] = new double[dim];
      rhslo[i] = new double[dim];
   }
   // Allocate work space for the inplace LU solver.
   double **workmathi = new double*[dim];
   double **workmatlo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      workmathi[i] = new double[dim];
      workmatlo[i] = new double[dim];
   }
   int *ipvt = new int[dim];
   double *workvechi = new double[dim];
   double *workveclo = new double[dim];
   // Copy the rhs vector into work space for inplace solver.
   double **workrhshi = new double*[degp1];
   double **workrhslo = new double*[degp1];

   for(int i=0; i<degp1; i++)
   {
      workrhshi[i] = new double[dim];
      workrhslo[i] = new double[dim];
   }
   double **resvechi = new double*[degp1];
   double **resveclo = new double*[degp1];

   for(int i=0; i<degp1; i++)
   {
      resvechi[i] = new double[dim];
      resveclo[i] = new double[dim];
   }
   double resmaxhi,resmaxlo;

/*
 * 4. initialize input, coefficient, evaluate, differentiate, and solve
 */
   // Define the initial input, a vector of ones.
   dbl2_unit_series_vector(dim,deg,inputhi,inputlo);

   if(vrblvl > 0)
   {
      cout << scientific << setprecision(16);
      cout << "The leading coefficients of the input series :" << endl;
      for(int i=0; i<dim; i++)
         cout << i << " : " << inputhi[i][0] << "  "
                            << inputlo[i][0] << endl;
   }
   for(int step=0; step<nbsteps; step++)
   {
      cout << "step " << step << " ..." << endl;

      dbl2_newton_step
         (dim,deg,nvr,idx,exp,nbrfac,expfac,cffhi,cfflo,acchi,acclo,
          inputhi,inputlo,outputhi,outputlo,funvalhi,funvallo,
          jacvalhi,jacvallo,rhshi,rhslo,solhi,sollo,workmathi,workmatlo,
          workvechi,workveclo,workrhshi,workrhslo,resvechi,resveclo,
          &resmaxhi,&resmaxlo,ipvt,vrblvl);
   }
   return 0;
}
