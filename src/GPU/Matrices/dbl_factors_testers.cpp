/* The file dbl_factors_testers.cpp define the functions specified in
   the file dbl_factors_testers.h. */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include "random_matrices.h"
#include "dbl_factorizations.h"

using namespace std;

void test_factors_real_lufac ( void )
{
   cout << "Give the dimension : ";
   int dim; cin >> dim;

   cout << "Generating a random matrix of dimension " << dim
        << " ..." << endl;

   double **A = new double*[dim];
   for(int i=0; i<dim; i++) A[i] = new double[dim];

   random_dbl_matrix(dim,dim,A);

   cout << scientific << setprecision(16);

   cout << "A random matrix :" << endl;
   for(int i=0; i<dim; i++)
      for(int j=0; j<dim; j++)
         cout << "A[" << i << "][" << j << "] : " << A[i][j] << endl;

   double *sol = new double[dim];
   for(int i=0; i<dim; i++) sol[i] = 1.0;

   double *rhs = new double[dim];
   for(int i=0; i<dim; i++)
   {
      rhs[i] = 0.0;
      for(int j=0; j<dim; j++) rhs[i] = rhs[i] + A[i][j]*sol[j];
   }
   cout << "The sums of the columns :" << endl;
   for(int i=0; i<dim; i++)
      cout << "b[" << i << "] : " << rhs[i] << endl;

   double *x = new double[dim];
   int *pivots = new int[dim];

   CPU_dbl_factors_lusolve(dim,A,pivots,rhs,x);

   cout << "The computed solution :" << endl;
   for(int i=0; i<dim; i++)
      cout << "x[" << i << "] : " << x[i] << endl;
}

void test_factors_cmplx_lufac ( void )
{
   cout << "Give the dimension : ";
   int dim; cin >> dim;

   cout << "Generating a random matrix of dimension " << dim
        << " ..." << endl;

   double **Are = new double*[dim];
   double **Aim = new double*[dim];
   for(int i=0; i<dim; i++)
   {
      Are[i] = new double[dim];
      Aim[i] = new double[dim];
   }
   random_cmplx_matrix(dim,dim,Are,Aim);

   cout << scientific << setprecision(16);

   cout << "A random matrix :" << endl;
   for(int i=0; i<dim; i++)
      for(int j=0; j<dim; j++)
         cout << "A[" << i << "][" << j << "] : "
              << Are[i][j] << "  " << Aim[i][j] << endl;

   double *solre = new double[dim];
   double *solim = new double[dim];
   for(int i=0; i<dim; i++)
   {
      solre[i] = 1.0;
      solim[i] = 0.0;
   }
   double *rhsre = new double[dim];
   double *rhsim = new double[dim];
   double accre,accim;

   for(int i=0; i<dim; i++)
   {
      rhsre[i] = 0.0;
      rhsim[i] = 0.0;
      for(int j=0; j<dim; j++) // rhs[i] = rhs[i] + A[i][j]*sol[j];
      {
         accre = Are[i][j]*solre[j] - Aim[i][j]*solim[j];
         accim = Aim[i][j]*solre[j] + Are[i][j]*solim[j];
         rhsre[i] = rhsre[i] + accre;
         rhsim[i] = rhsim[i] + accim;
      }
   }
   cout << "The sums of the columns :" << endl;
   for(int i=0; i<dim; i++)
      cout << "b[" << i << "] : " << rhsre[i] << "  " << rhsim[i] << endl;

   double *xre = new double[dim];
   double *xim = new double[dim];
   int *pivots = new int[dim];

   CPU_cmplx_factors_lusolve(dim,Are,Aim,pivots,rhsre,rhsim,xre,xim);

   cout << "The computed solution :" << endl;
   for(int i=0; i<dim; i++)
      cout << "x[" << i << "] : " << xre[i] << "  " << xim[i] << endl;
}
