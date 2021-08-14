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

   cout << "Give the verbose level (1 to see all numbers) : ";
   int verbose; cin >> verbose;

   cout << "Generating a random matrix of dimension " << dim
        << " ..." << endl;

   double **A = new double*[dim];
   for(int i=0; i<dim; i++) A[i] = new double[dim];

   random_dbl_matrix(dim,dim,A);

   if(verbose > 0)
   {
      cout << scientific << setprecision(16);

      cout << "A random matrix :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
            cout << "A[" << i << "][" << j << "] : " << A[i][j] << endl;
   }
   double *sol = new double[dim];
   for(int i=0; i<dim; i++) sol[i] = 1.0;

   double *rhs = new double[dim];
   for(int i=0; i<dim; i++)
   {
      rhs[i] = 0.0;
      for(int j=0; j<dim; j++) rhs[i] = rhs[i] + A[i][j]*sol[j];
   }
   if(verbose > 0)
   {
      cout << "The sums of the columns :" << endl;
      for(int i=0; i<dim; i++)
         cout << "b[" << i << "] : " << rhs[i] << endl;
   }
   double *x = new double[dim];
   int *pivots = new int[dim];

   CPU_dbl_factors_lusolve(dim,A,pivots,rhs,x);

   if(verbose > 0)
   {
      cout << "The computed solution :" << endl;
      for(int i=0; i<dim; i++)
         cout << "x[" << i << "] : " << x[i] << endl;
   }
   double error = 0.0;
   for(int i=0; i<dim; i++)
      error = error + abs(x[i] - 1.0);

   cout << scientific << setprecision(2);
   cout << "Sum of errors : " << error << endl;
}

void test_factors_cmplx_lufac ( void )
{
   cout << "Give the dimension : ";
   int dim; cin >> dim;

   cout << "Give the verbose level (1 to see all numbers) : ";
   int verbose; cin >> verbose;

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

   if(verbose > 0)
   {
      cout << scientific << setprecision(16);

      cout << "A random matrix :" << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
            cout << "A[" << i << "][" << j << "] : "
                 << Are[i][j] << "  " << Aim[i][j] << endl;
   }
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
   if(verbose > 0)
   {
      cout << "The sums of the columns :" << endl;
      for(int i=0; i<dim; i++)
         cout << "b[" << i << "] : " << rhsre[i] << "  " << rhsim[i] << endl;
   }
   double *xre = new double[dim];
   double *xim = new double[dim];
   int *pivots = new int[dim];

   CPU_cmplx_factors_lusolve(dim,Are,Aim,pivots,rhsre,rhsim,xre,xim);

   if(verbose > 0)
   {
      cout << "The computed solution :" << endl;
      for(int i=0; i<dim; i++)
         cout << "x[" << i << "] : " << xre[i] << "  " << xim[i] << endl;
   }
   double error = 0.0;
   for(int i=0; i<dim; i++)
      error = error + abs(xre[i] - 1.0) + abs(xim[i]);

   cout << scientific << setprecision(2);
   cout << "Sum of errors : " << error << endl;
}

int test_real_qr_factors
 ( int nrows, int ncols, double **A, double **Q, double **R,
   double tol, int verbose )
{
   double **QT = new double*[nrows];
   double **QTQ = new double*[nrows];
   double **QTA = new double*[nrows];

   for(int i=0; i<nrows; i++)
   {
      QT[i] = new double[nrows];
      QTQ[i] = new double[nrows];
      QTA[i] = new double[ncols];
   }

   if(verbose > 0) 
   {
      cout << scientific << setprecision(16);
      cout << "The matrix Q :" << endl;
   }
   for(int i=0; i<nrows; i++)
      for(int j=0; j<nrows; j++)
      {
         if(verbose > 0)
            cout << "Q[" << i << "][" << j << "] : " << Q[i][j] << endl;
         QT[j][i] = Q[i][j];
      }
   CPU_dbl_factors_matmatmul(nrows,nrows,nrows,QT,Q,QTQ);

   double errorQ = 0.0;

   if(verbose > 0) cout << "The matrix transpose(Q)*Q :" << endl;
   for(int i=0; i<nrows; i++)
      for(int j=0; j<nrows; j++)
      {
         if(verbose > 0)
            cout << "Q'*Q[" << i << "][" << j << "] : "
                 << QTQ[i][j] << endl;
         if(i == j)
            errorQ = errorQ + abs(QTQ[i][j] - 1.0);
         else
            errorQ = errorQ + abs(QTQ[i][j]);
      }

   cout << scientific << setprecision(2);
   cout << "Sum of errors on |Q^T*Q - I| : " << errorQ << endl;

   if(verbose > 0)
   {
      cout << scientific << setprecision(16);
      cout << "The matrix R :" << endl;
      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++)
            cout << "R[" << i << "][" << j << "] : " << R[i][j] << endl;
   }
   CPU_dbl_factors_matmatmul(nrows,nrows,ncols,QT,A,QTA);

   double errorR = 0.0;

   if(verbose > 0) cout << "The matrix transpose(Q)*A :" << endl;
   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
      {
         if(verbose > 0)
            cout << "Q'*A[" << i << "][" << j << "] : "
                 << QTA[i][j] << endl;
         errorR = errorR + abs(R[i][j] - QTA[i][j]);
      }

   cout << scientific << setprecision(2);
   cout << "Sum of errors on |Q^T*A - R| : " << errorR << endl;

   for(int i=0; i<nrows; i++)
   {
      free(QT[i]); free(QTQ[i]); free(QTA[i]);
   }
   free(QT); free(QTQ); free(QTA);

   return int(errorQ + errorR > tol);
}

int test_real_qr_factors_probe
 ( int nrows, int ncols, double **A, double **Q, double **R,
   double tol, int nbprobes, int verbose )
{
   int rowidx,colidx;
   double Qsum,Rsum;
   double errorQ = 0.0;
   double errorR = 0.0;

   for(int p=0; p<nbprobes; p++)
   {
      rowidx = rand() % nrows;
      colidx = rand() % ncols;

      if(verbose > 0)
      {
         cout << "Probing row index : " << rowidx
              << ", column index : " << colidx << "." << endl;
      }
      Qsum = 0.0;
      Rsum = 0.0;

      for(int i=0; i<nrows; i++)
      {
         Qsum = Qsum + Q[i][rowidx]*Q[i][colidx];
         Rsum = Rsum + Q[i][rowidx]*A[i][colidx];
      }
      if(verbose > 0)
      {
         cout << scientific << setprecision(16);
         cout << "Q^T*Q[" << rowidx << "][" << rowidx << "] : "
              << Qsum << endl;
         cout << "Q^T*A[" << rowidx << "][" << colidx << "] : "
              << Rsum << endl;
         cout << "    R[" << rowidx << "][" << colidx << "] : "
              << R[rowidx][colidx] << endl;
      }
      if(rowidx == colidx)
         errorQ = errorQ + fabs(Qsum - 1.0);
      else
         errorQ = errorQ + fabs(Qsum);

      errorR = errorR + fabs(Rsum - R[rowidx][colidx]);
   }
   cout << scientific << setprecision(2);
   cout << "Sum of errors on |Q^T*Q - I| : " << errorQ << endl;
   cout << "Sum of errors on |Q^T*A - R| : " << errorR << endl;

   return int(errorQ + errorR > tol);
}

int test_cmplx_qr_factors
 ( int nrows, int ncols, double **Are, double **Aim,
   double **Qre, double **Qim, double **Rre, double **Rim,
   double tol, int verbose )
{
   double **QHre = new double*[nrows];
   double **QHim = new double*[nrows];
   double **QHQre = new double*[nrows];
   double **QHQim = new double*[nrows];
   double **QHAre = new double*[nrows];
   double **QHAim = new double*[nrows];

   for(int i=0; i<nrows; i++)
   {
      QHre[i] = new double[nrows];
      QHim[i] = new double[nrows];
      QHQre[i] = new double[nrows];
      QHQim[i] = new double[nrows];
      QHAre[i] = new double[ncols];
      QHAim[i] = new double[ncols];
   }
   if(verbose > 0)
   {
      cout << scientific << setprecision(16);
      cout << "The matrix Q :" << endl;
   }
   for(int i=0; i<nrows; i++)
      for(int j=0; j<nrows; j++)
      {
         if(verbose > 0)
            cout << "Q[" << i << "][" << j << "] : "
                 << Qre[i][j] << "  " << Qim[i][j] << endl;
         QHre[j][i] = Qre[i][j];
         QHim[j][i] = -Qim[i][j]; // Hermitian transpose
      }

   CPU_cmplx_factors_matmatmul
      (nrows,nrows,nrows,QHre,QHim,Qre,Qim,QHQre,QHQim);

   double errorQ = 0.0;

   if(verbose > 0) cout << "The matrix transpose(Q)*Q :" << endl;
   for(int i=0; i<nrows; i++)
      for(int j=0; j<nrows; j++)
      {
         if(verbose > 0)
            cout << "Q^H*Q[" << i << "][" << j << "] : "
                 << QHQre[i][j] << "  " << QHQim[i][j] << endl;
         if(i == j)
            errorQ = errorQ + abs(QHQre[i][j] - 1.0) + abs(QHQim[i][j]);
         else
            errorQ = errorQ + abs(QHQre[i][j]) + abs(QHQim[i][j]);
      }

   cout << scientific << setprecision(2);
   cout << "Sum of errors on |Q^H*Q - I| : " << errorQ << endl;

   if(verbose > 0)
   {
      cout << scientific << setprecision(16);
      cout << "The matrix R :" << endl;
      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++)
            cout << "R[" << i << "][" << j << "] : "
                 << Rre[i][j] << "  " << Rim[i][j] << endl;
   }
   CPU_cmplx_factors_matmatmul
      (nrows,nrows,ncols,QHre,QHim,Are,Aim,QHAre,QHAim);

   double errorR = 0.0;

   if(verbose > 0) cout << "The matrix transpose(Q)*A :" << endl;
   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
      {
         if(verbose > 0)
            cout << "Q^H*A[" << i << "][" << j << "] : "
                 << QHAre[i][j] << "  " << QHAim[i][j] << endl;
         errorR = errorR + abs(Rre[i][j] - QHAre[i][j])
                         + abs(Rim[i][j] - QHAim[i][j]);
      }

   cout << scientific << setprecision(2);
   cout << "Sum of errors on |Q^H*A - R| : " << errorR << endl;

   for(int i=0; i<nrows; i++)
   {
      free(QHre[i]); free(QHQre[i]); free(QHAre[i]);
      free(QHim[i]); free(QHQim[i]); free(QHAim[i]);
   }
   free(QHre); free(QHQre); free(QHAre);
   free(QHim); free(QHQim); free(QHAim);

   return int(errorQ + errorR > tol);
}

int test_cmplx_qr_factors_probe
 ( int nrows, int ncols, double **Are, double **Aim,
   double **Qre, double **Qim, double **Rre, double **Rim,
   double tol, int nbprobes, int verbose )
{
   int rowidx,colidx;
   double Qsumre,Rsumre;
   double Qsumim,Rsumim;
   double acc;
   double errorQ = 0.0;
   double errorR = 0.0;

   for(int p=0; p<nbprobes; p++)
   {
      rowidx = rand() % nrows;
      colidx = rand() % ncols;

      if(verbose > 0)
      {
         cout << "Probing row index : " << rowidx
              << ", column index : " << colidx << "." << endl;
      }
      Qsumre = 0.0;
      Qsumim = 0.0;
      Rsumre = 0.0;
      Rsumim = 0.0;

      for(int i=0; i<nrows; i++)
      {
         // multiply Q^H with Q
         Qsumre = Qsumre + Qre[i][rowidx]*Qre[i][colidx];
         Qsumre = Qsumre + Qim[i][rowidx]*Qim[i][colidx];
         Qsumim = Qsumim - Qim[i][rowidx]*Qre[i][colidx];
         Qsumim = Qsumim + Qre[i][rowidx]*Qim[i][colidx];
         // multiply Q^H with A
         Rsumre = Rsumre + Qre[i][rowidx]*Are[i][colidx];
         Rsumre = Rsumre + Qim[i][rowidx]*Aim[i][colidx];
         Rsumim = Rsumim - Qim[i][rowidx]*Are[i][colidx];
         Rsumim = Rsumim + Qre[i][rowidx]*Aim[i][colidx];
      }
      if(verbose > 0)
      {
         cout << scientific << setprecision(16);
         cout << "Q^T*Q[" << rowidx << "][" << rowidx << "] : "
              << Qsumre << "  " << Qsumim << endl;
         cout << "Q^T*A[" << rowidx << "][" << colidx << "] : "
              << Rsumre << "  " << Rsumim << endl;
         cout << "    R[" << rowidx << "][" << colidx << "] : "
              << Rre[rowidx][colidx] << "  "
              << Rim[rowidx][colidx] << endl;
      }
      if(rowidx == colidx)
      {
         errorQ = errorQ + fabs(Qsumre - 1.0) + fabs(Qsumim);
      }
      else
      {
         errorQ = errorQ + fabs(Qsumre) + fabs(Qsumim);
      }
      errorR = errorR + fabs(Rsumre - Rre[rowidx][colidx])
                      + fabs(Rsumim - Rim[rowidx][colidx]);
   }
   cout << scientific << setprecision(2);
   cout << "Sum of errors on |Q^T*Q - I| : " << errorQ << endl;
   cout << "Sum of errors on |Q^T*A - R| : " << errorR << endl;

   return int(errorQ + errorR > tol);
}

void test_factors_real_houseqr ( void )
{
   cout << "Give the number of rows : ";
   int nrows; cin >> nrows;

   cout << "Give the number of columns : ";
   int ncols; cin >> ncols;

   cout << "Give the verbose level (1 to see all numbers) : ";
   int verbose; cin >> verbose;

   cout << "Generating a random " << nrows
        << "-by-" << ncols << " matrix ..." << endl;

   double **A = new double*[nrows];
   double **Q = new double*[nrows];
   double **R = new double*[nrows];

   for(int i=0; i<nrows; i++)
   {
      A[i] = new double[ncols];
      Q[i] = new double[nrows];
      R[i] = new double[ncols];
   }
   random_dbl_matrix(nrows,ncols,A);

   if(verbose > 0)
   {
      cout << scientific << setprecision(16);

      cout << "A random matrix :" << endl;
      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++)
            cout << "A[" << i << "][" << j << "] : " << A[i][j] << endl;
   }
   CPU_dbl_factors_houseqr(nrows,ncols,A,Q,R);

   const double tol = 1.0e-12;
   const int fail = test_real_qr_factors(nrows,ncols,A,Q,R,tol,verbose);
   if(fail == 0)
      cout << "The test succeeded." << endl;
   else
   {
      cout << scientific << setprecision(2);
      cout << "The test failed for tol = " << tol << "." << endl;
   }
   for(int i=0; i<nrows; i++)
   {
      free(A[i]); free(Q[i]); free(R[i]);
   }
   free(A); free(Q); free(R);
}

void test_factors_cmplx_houseqr ( void )
{
   cout << "Give the number of rows : ";
   int nrows; cin >> nrows;

   cout << "Give the number of columns : ";
   int ncols; cin >> ncols;

   cout << "Give the verbose level (1 to see all numbers) : ";
   int verbose; cin >> verbose;

   cout << "Generating a random " << nrows
        << "-by-" << ncols << " matrix ..." << endl;

   double **Are = new double*[nrows];
   double **Aim = new double*[nrows];
   double **Qre = new double*[nrows];
   double **Qim = new double*[nrows];
   double **Rre = new double*[nrows];
   double **Rim = new double*[nrows];

   for(int i=0; i<nrows; i++)
   {
      Are[i] = new double[ncols];
      Aim[i] = new double[ncols];
      Qre[i] = new double[nrows];
      Qim[i] = new double[nrows];
      Rre[i] = new double[ncols];
      Rim[i] = new double[ncols];
   }
   random_cmplx_matrix(nrows,ncols,Are,Aim);

   if(verbose > 0)
   {
      cout << scientific << setprecision(16);

      cout << "A random matrix :" << endl;
      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++)
            cout << "A[" << i << "][" << j << "] : "
                 << Are[i][j] << "  " << Aim[i][j] << endl;
   }
   CPU_cmplx_factors_houseqr(nrows,ncols,Are,Aim,Qre,Qim,Rre,Rim);

   const double tol = 1.0e-12;
   const int fail = test_cmplx_qr_factors
      (nrows,ncols,Are,Aim,Qre,Qim,Rre,Rim,tol,verbose);
   if(fail == 0)
      cout << "The test succeeded." << endl;
   else
   {
      cout << scientific << setprecision(2);
      cout << "The test failed for tol = " << tol << "." << endl;
   }
   for(int i=0; i<nrows; i++)
   {
      free(Are[i]); free(Qre[i]); free(Rre[i]);
      free(Aim[i]); free(Qim[i]); free(Rim[i]);
   }
   free(Are); free(Qre); free(Rre);
   free(Aim); free(Qim); free(Rim);
}
