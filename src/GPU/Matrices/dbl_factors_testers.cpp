/* The file dbl_factors_testers.cpp define the functions specified in
   the file dbl_factors_testers.h. */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
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

void test_factors_real_houseqr ( void )
{
   cout << "Give the number of rows : ";
   int nrows; cin >> nrows;

   cout << "Give the number of columns : ";
   int ncols; cin >> ncols;

   cout << "Generating a random " << nrows
        << "-by-" << ncols << " matrix ..." << endl;

   double **A = new double*[nrows];
   double **Q = new double*[nrows];
   double **QT = new double*[nrows];
   double **QTQ = new double*[nrows];
   double **R = new double*[nrows];
   double **QTA = new double*[nrows];

   for(int i=0; i<nrows; i++)
   {
      A[i] = new double[ncols];
      Q[i] = new double[nrows];
      QT[i] = new double[nrows];
      QTQ[i] = new double[nrows];
      R[i] = new double[ncols];
      QTA[i] = new double[ncols];
   }

   random_dbl_matrix(nrows,ncols,A);

   cout << scientific << setprecision(16);

   cout << "A random matrix :" << endl;
   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
         cout << "A[" << i << "][" << j << "] : " << A[i][j] << endl;

   CPU_dbl_factors_houseqr(nrows,ncols,A,Q,R);

   cout << "The matrix Q :" << endl;
   for(int i=0; i<nrows; i++)
      for(int j=0; j<nrows; j++)
      {
         cout << "Q[" << i << "][" << j << "] : " << Q[i][j] << endl;
         QT[j][i] = Q[i][j];
      }

   cout << "The matrix R :" << endl;
   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
         cout << "R[" << i << "][" << j << "] : " << R[i][j] << endl;

   CPU_dbl_factors_matmatmul(nrows,nrows,nrows,QT,Q,QTQ);

   double error = 0.0;

   cout << "The matrix transpose(Q)*Q :" << endl;
   for(int i=0; i<nrows; i++)
      for(int j=0; j<nrows; j++)
      {
         cout << "Q'*Q[" << i << "][" << j << "] : " << QTQ[i][j] << endl;
         if(i == j)
            error = error + abs(QTQ[i][j] - 1.0);
         else
            error = error + abs(QTQ[i][j]);
      }

   cout << scientific << setprecision(2);
   cout << "Sum of errors : " << error << endl;

   CPU_dbl_factors_matmatmul(nrows,nrows,ncols,QT,A,QTA);

   error = 0.0;

   cout << scientific << setprecision(16);

   cout << "The matrix transpose(Q)*A :" << endl;
   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
      {
         cout << "Q'*A[" << i << "][" << j << "] : " << QTA[i][j] << endl;
         error = error + abs(R[i][j] - QTA[i][j]);
      }

   cout << scientific << setprecision(2);
   cout << "Sum of errors : " << error << endl;
}
