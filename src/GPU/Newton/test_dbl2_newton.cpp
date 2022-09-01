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

   int dim,deg,size,posvals,vrblvl,nbritr,nbsteps,szt,nbt,mode;
   prompt_dimensions(&dim,&deg,&size,&posvals,&vrblvl,&nbritr,&nbsteps);

   szt = dim; nbt = 1; mode = 1; // stub values

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
   double **inputhi_h = new double*[dim];
   double **inputlo_h = new double*[dim];
   double **inputhi_d = new double*[dim];
   double **inputlo_d = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      inputhi_h[i] = new double[degp1];
      inputlo_h[i] = new double[degp1];
      inputhi_d[i] = new double[degp1];
      inputlo_d[i] = new double[degp1];
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
   // The function values are power series truncated at degree deg.
   double **funvalhi_h = new double*[dim];
   double **funvallo_h = new double*[dim];
   double **funvalhi_d = new double*[dim];
   double **funvallo_d = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      funvalhi_h[i] = new double[degp1];
      funvallo_h[i] = new double[degp1];
      funvalhi_d[i] = new double[degp1];
      funvallo_d[i] = new double[degp1];
   }
   // The derivatives in the output are a series truncated at degree deg.
   // The coefficients of the series are matrices of dimension dim.
   double ***jacvalhi_h = new double**[degp1];
   double ***jacvallo_h = new double**[degp1];
   double ***jacvalhi_d = new double**[degp1];
   double ***jacvallo_d = new double**[degp1];

   for(int i=0; i<degp1; i++) // jacval[i] is matrix of dimension dim
   {
      jacvalhi_h[i] = new double*[dim];
      jacvallo_h[i] = new double*[dim];
      jacvalhi_d[i] = new double*[dim];
      jacvallo_d[i] = new double*[dim];

      for(int j=0; j<dim; j++)
      {
         jacvalhi_h[i][j] = new double[dim];
         jacvallo_h[i][j] = new double[dim];
         jacvalhi_d[i][j] = new double[dim];
         jacvallo_d[i][j] = new double[dim];
      }
   }
/*
 * 3. allocate space to solve the linearized power series system
 */
   // The solution x(t) to jacval(t)*x(t) = -funval(t) in linearized
   // format is a series truncated at degree deg, with as coefficients
   // arrays of dimension dim.

   double **solhi_h = new double*[degp1];
   double **sollo_h = new double*[degp1];
   double **solhi_d = new double*[degp1];
   double **sollo_d = new double*[degp1];

   for(int i=0; i<degp1; i++)
   {
      solhi_h[i] = new double[dim];
      sollo_h[i] = new double[dim];
      solhi_d[i] = new double[dim];
      sollo_d[i] = new double[dim];
   }
   // The right hand side -funval(t) in linearized format is a series
   // truncated at degree deg, with arrays of dimension dim as coefficients.
   double **rhshi_h = new double*[degp1];
   double **rhslo_h = new double*[degp1];
   double **rhshi_d = new double*[degp1];
   double **rhslo_d = new double*[degp1];

   for(int i=0; i<degp1; i++)
   {
      rhshi_h[i] = new double[dim];
      rhslo_h[i] = new double[dim];
      rhshi_d[i] = new double[dim];
      rhslo_d[i] = new double[dim];
   }
   // Copy the rhs vector into work space for inplace solver.
   double **urhshi_h = new double*[degp1];
   double **urhslo_h = new double*[degp1];
   double **urhshi_d = new double*[degp1];
   double **urhslo_d = new double*[degp1];
   for(int i=0; i<degp1; i++)
   {
      urhshi_h[i] = new double[dim];
      urhslo_h[i] = new double[dim];
      urhshi_d[i] = new double[dim];
      urhslo_d[i] = new double[dim];
   }
   // Allocate work space for the inplace LU solver.
   double **workmathi = new double*[dim];
   double **workmatlo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      workmathi[i] = new double[dim];
      workmatlo[i] = new double[dim];
   }
   // int *ipvt = new int[dim];
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
   double **Qhi_h = new double*[dim];
   double **Qlo_h = new double*[dim];
   double **Qhi_d = new double*[dim];
   double **Qlo_d = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      Qhi_h[i] = new double[dim];
      Qlo_h[i] = new double[dim];
      Qhi_d[i] = new double[dim];
      Qlo_d[i] = new double[dim];
   }
   double **Rhi_h = new double*[dim];
   double **Rlo_h = new double*[dim];
   double **Rhi_d = new double*[dim];
   double **Rlo_d = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      Rhi_h[i] = new double[dim];
      Rlo_h[i] = new double[dim];
      Rhi_d[i] = new double[dim];
      Rlo_d[i] = new double[dim];
   }
/*
 * 4. initialize input, coefficient, evaluate, differentiate, and solve
 */
   // Define the initial input, a vector of ones.
   dbl2_unit_series_vector(dim,deg,inputhi_h,inputlo_h);

   if(vrblvl > 0)
   {
      cout << scientific << setprecision(16);
      cout << "The leading coefficients of the input series :" << endl;
      for(int i=0; i<dim; i++)
         cout << i << " : " << inputhi_h[i][0] << "  "
                            << inputlo_h[i][0] << endl;
   }
   for(int step=0; step<nbsteps; step++)
   {
      cout << "step " << step << " ..." << endl;
/*
      dbl2_newton_lustep
         (dim,deg,nvr,idx,exp,nbrfac,expfac,cffhi,cfflo,acchi,acclo,
          inputhi,inputlo,outputhi,outputlo,funvalhi,funvallo,
          jacvalhi,jacvallo,rhshi,rhslo,solhi,sollo,workmathi,workmatlo,
          workvechi,workveclo,workrhshi,workrhslo,resvechi,resveclo,
          &resmaxhi,&resmaxlo,ipvt,vrblvl);
 */
      dbl2_newton_qrstep
         (szt,nbt,dim,deg,nvr,idx,exp,nbrfac,expfac,cffhi,cfflo,acchi,acclo,
          inputhi_h,inputlo_h,inputhi_d,inputlo_d,
          outputhi_h,outputlo_h,outputhi_d,outputlo_d,
          funvalhi_h,funvallo_h,funvalhi_d,funvallo_d,
          jacvalhi_h,jacvallo_h,jacvalhi_d,jacvallo_d,
          rhshi_h,rhslo_h,rhshi_d,rhslo_d,urhshi_h,urhslo_h,urhshi_d,urhslo_d,
          solhi_h,sollo_h,solhi_d,sollo_d,
          Qhi_h,Qlo_h,Qhi_d,Qlo_d,Rhi_h,Rlo_h,Rhi_d,Rlo_d,
          workmathi,workmatlo,workvechi,workveclo,
          resvechi,resveclo,&resmaxhi,&resmaxlo,vrblvl,mode);
   }
   return 0;
}
