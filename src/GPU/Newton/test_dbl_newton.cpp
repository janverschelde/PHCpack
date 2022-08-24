/* Tests operations on Newton for series in double precision */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <time.h>
#include "unimodular_matrices.h"
#include "dbl_newton_testers.h"

using namespace std;

void exponents_check ( int dim, int **rowsA );
/*
 * DESCRIPTION :
 *   Writes the exponents of the solution for checking ... */

int main ( void )
{
   srand(time(NULL));

   cout << "testing Newton ..." << endl;

   int dim,deg,size,posvals,vrblvl,nbritr,nbsteps;
   int seed,szt,nbt,mode;

   prompt_newton_setup
      (&seed,&szt,&nbt,&dim,&deg,&size,&posvals,&vrblvl,&mode,
       &nbritr,&nbsteps);

   if(seed == 0)
      srand(time(NULL));
   else
      srand(seed);

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
   double **input_h = new double*[dim];
   double **input_d = new double*[dim];
   for(int i=0; i<dim; i++)
   {
       input_h[i] = new double[degp1];
       input_d[i] = new double[degp1];
   }
   // allocate memory for coefficients and the output
   double *acc = new double[degp1]; // accumulated power series
   double **cff = new double*[dim]; // the coefficients of monomials
   for(int i=0; i<dim; i++) cff[i] = new double[degp1];
   double ***output = new double**[dim];
   for(int i=0; i<dim; i++)
   {
      output[i] = new double*[dim+1];
      for(int j=0; j<=dim; j++) output[i][j] = new double[degp1];
   }
   // The function values are power series truncated at degree deg.
   double **funval = new double*[dim];
   for(int i=0; i<dim; i++) funval[i] = new double[degp1];
   // The derivatives in the output are a series truncated at degree deg.
   // The coefficients of the series are matrices of dimension dim.
   double ***jacval = new double**[degp1];
   for(int i=0; i<degp1; i++) // jacval[i] is matrix of dimension dim
   {
      jacval[i] = new double*[dim];
      for(int j=0; j<dim; j++) jacval[i][j] = new double[dim];
   }
/*
 * 3. allocate space to solve the linearized power series system
 */
   // The solution x(t) to jacval(t)*x(t) = -funval(t) in linearized
   // format is a series truncated at degree deg, with as coefficients
   // arrays of dimension dim.
   double **sol_h = new double*[degp1];
   double **sol_d = new double*[degp1];
   for(int i=0; i<degp1; i++) 
   {
      sol_h[i] = new double[dim];
      sol_d[i] = new double[dim];
   }
   // The right hand side -funval(t) in linearized format is a series
   // truncated at degree deg, with arrays of dimension dim as coefficients.
   double **rhs = new double*[degp1];
   for(int i=0; i<degp1; i++) rhs[i] = new double[dim];
   // Allocate work space for the inplace LU solver.
   double **workmat = new double*[dim];
   for(int i=0; i<dim; i++) workmat[i] = new double[dim];
   int *ipvt = new int[dim];
   double *workvec = new double[dim];
   // Copy the rhs vector into work space for inplace solver.
   double **urhs_h = new double*[degp1];
   double **urhs_d = new double*[degp1];
   for(int i=0; i<degp1; i++)
   {
      urhs_h[i] = new double[dim];
      urhs_d[i] = new double[dim];
   }
   double **resvec = new double*[degp1];
   for(int i=0; i<degp1; i++) resvec[i] = new double[dim];
   double resmax;
   double **Q_h = new double*[dim];
   double **Q_d = new double*[dim];
   for(int i=0; i<dim; i++)
   {
      Q_h[i] = new double[dim];
      Q_d[i] = new double[dim];
   }
   double **R_h = new double*[dim];
   double **R_d = new double*[dim];
   for(int i=0; i<dim; i++)
   {
      R_h[i] = new double[dim];
      R_d[i] = new double[dim];
   }
/*
 * 4. initialize input, coefficient, evaluate, differentiate, and solve
 */
   // Define the initial input, a vector of ones.
   dbl_unit_series_vector(dim,deg,input_h);

   if(vrblvl > 0)
   {
      cout << scientific << setprecision(16);
      cout << "The leading coefficients of the input series :" << endl;
      for(int i=0; i<dim; i++)
         cout << i << " : " << input_h[i][0] << endl;
   }
   for(int step=0; step<nbsteps; step++)
   {
      cout << "*** running Newton step " << step << " ***" << endl;

      dbl_newton_qrstep
         (szt,nbt,dim,deg,nvr,idx,exp,nbrfac,expfac,cff,acc,
          input_h,input_d,output,funval,jacval,
          rhs,urhs_h,urhs_d,sol_h,sol_d,Q_h,Q_d,R_h,R_d,
          workmat,workvec,resvec,&resmax,vrblvl,mode);
   }
   exponents_check(dim, rowsA);

   return 0;
}

void exponents_check ( int dim, int **rowsA )
{
   cout << "The matrix :" << endl;
   write_exponent_matrix(dim, rowsA);

   int *expsol = new int[dim];
   int **copyA = new int*[dim];  // copy of A
   int **unimd = new int*[dim];  // unimodular transformation

   for(int i=0; i<dim; i++)      // initialize the data
   {
      unimd[i] = new int[dim];
      copyA[i] = new int[dim];
   }
   copy_integer_matrix(dim, rowsA, copyA);
   int sing = lower_triangulate(dim, copyA, unimd, 0);
   exponent_forward_substitution(dim, copyA, expsol);
   cout << "exponents after forward substitution :" << endl;
   for(int i=0; i<dim; i++) cout << " " << expsol[i];
   cout << endl;
   exponent_unimodular_transformation(dim, unimd, expsol);
   cout << "exponents after unimodular transformation :" << endl;
   for(int i=0; i<dim; i++) cout << " " << expsol[i];
   cout << endl;

   for(int i=0; i<dim; i++)
   {
      free(unimd[i]);
      free(copyA[i]);
   }
   free(unimd); free(copyA); free(expsol);
}
