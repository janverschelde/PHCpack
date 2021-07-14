// The file dbl_baqr_testers.cpp defines the function with prototypes in
// the file dbl_baqr_testers.h.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#include "random_matrices.h"
#include "dbl_factorizations.h"
#include "dbl_test_utilities.h"
#include "dbl_baqr_host.h"

using namespace std;

void test_real_blocked_qr ( void )
{
   cout << "Give the size of each tile : ";
   int sizetile; cin >> sizetile;

   cout << "Give the number of tiles : ";
   int numtiles; cin >> numtiles;

   const int ncols = sizetile*numtiles;

   cout << "Give the number of rows (" << " >= " << ncols << " ) : ";
   int nrows; cin >> nrows;

   cout << "Give the verbose level (1 to see all numbers) : ";
   int verbose; cin >> verbose;

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

   if(verbose > 0)
   {
      cout << scientific << setprecision(16);

      cout << "A random matrix :" << endl;
      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++)
            cout << "A[" << i << "][" << j << "] : " << A[i][j] << endl;
   }
   bool vrb = (verbose > 0);
   CPU_dbl_blocked_houseqr(nrows,ncols,sizetile,numtiles,A,Q,R,vrb);

   if(verbose > 0) cout << "The matrix Q :" << endl;
   for(int i=0; i<nrows; i++)
      for(int j=0; j<nrows; j++)
      {
         if(verbose > 0)
            cout << "Q[" << i << "][" << j << "] : " << Q[i][j] << endl;
         QT[j][i] = Q[i][j];
      }
  
   if(verbose > 0)
   {
      cout << "The matrix R :" << endl;
      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++)
            cout << "R[" << i << "][" << j << "] : " << R[i][j] << endl;
   }
   CPU_dbl_factors_matmatmul(nrows,nrows,nrows,QT,Q,QTQ);

   double error = 0.0;

   if(verbose > 0) cout << "The matrix transpose(Q)*Q :" << endl;
   for(int i=0; i<nrows; i++)
      for(int j=0; j<nrows; j++)
      {
         if(verbose > 0)
            cout << "Q'*Q[" << i << "][" << j << "] : "
                 << QTQ[i][j] << endl;
         if(i == j)
            error = error + abs(QTQ[i][j] - 1.0);
         else
            error = error + abs(QTQ[i][j]);
      }

   cout << scientific << setprecision(2);
   cout << "Sum of errors on |Q^T*Q - I| : " << error << endl;

   CPU_dbl_factors_matmatmul(nrows,nrows,ncols,QT,A,QTA);

   error = 0.0;

   cout << scientific << setprecision(16);

   if(verbose > 0) cout << "The matrix transpose(Q)*A :" << endl;
   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
      {
         if(verbose > 0)
            cout << "Q'*A[" << i << "][" << j << "] : "
                 << QTA[i][j] << endl;
         error = error + abs(R[i][j] - QTA[i][j]);
      }

   cout << scientific << setprecision(2);
   cout << "Sum of errors on |Q^T*A - R| : " << error << endl;
}

void test_cmplx_blocked_qr ( void )
{
   cout << "Give the size of each tile : ";
   int sizetile; cin >> sizetile;

   cout << "Give the number of tiles : ";
   int numtiles; cin >> numtiles;

   const int ncols = sizetile*numtiles;

   cout << "Give the number of rows (" << " >= " << ncols << " ) : ";
   int nrows; cin >> nrows;

   cout << "Give the verbose level (1 to see all numbers) : ";
   int verbose; cin >> verbose;

   cout << "Generating a random " << nrows
        << "-by-" << ncols << " matrix ..." << endl;

   double **Are = new double*[nrows];
   double **Aim = new double*[nrows];
   double **Qre = new double*[nrows];
   double **Qim = new double*[nrows];
   double **QHre = new double*[nrows];
   double **QHim = new double*[nrows];
   double **QHQre = new double*[nrows];
   double **QHQim = new double*[nrows];
   double **Rre = new double*[nrows];
   double **Rim = new double*[nrows];
   double **QHAre = new double*[nrows];
   double **QHAim = new double*[nrows];

   for(int i=0; i<nrows; i++)
   {
      Are[i] = new double[ncols];
      Aim[i] = new double[ncols];
      Qre[i] = new double[nrows];
      Qim[i] = new double[nrows];
      QHre[i] = new double[nrows];
      QHim[i] = new double[nrows];
      QHQre[i] = new double[nrows];
      QHQim[i] = new double[nrows];
      Rre[i] = new double[ncols];
      Rim[i] = new double[ncols];
      QHAre[i] = new double[ncols];
      QHAim[i] = new double[ncols];
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
   bool vrb = (verbose > 0);
   CPU_cmplx_blocked_houseqr
      (nrows,ncols,sizetile,numtiles,Are,Aim,Qre,Qim,Rre,Rim,vrb);

   if(verbose > 0) cout << "The matrix Q :" << endl;
   for(int i=0; i<nrows; i++)
      for(int j=0; j<nrows; j++)
      {
         if(verbose > 0)
            cout << "Q[" << i << "][" << j << "] : "
                 << Qre[i][j] << "  " << Qim[i][j] << endl;
         QHre[j][i] = Qre[i][j];
         QHim[j][i] = -Qim[i][j]; // Hermitian transpose
      }

   if(verbose > 0)
   {
      cout << "The matrix R :" << endl;
      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++)
            cout << "R[" << i << "][" << j << "] : "
                 << Rre[i][j] << "  " << Rim[i][j] << endl;
   }
   CPU_cmplx_factors_matmatmul
      (nrows,nrows,nrows,QHre,QHim,Qre,Qim,QHQre,QHQim);

   double error = 0.0;

   if(verbose > 0) cout << "The matrix transpose(Q)*Q :" << endl;
   for(int i=0; i<nrows; i++)
      for(int j=0; j<nrows; j++)
      {
         if(verbose > 0)
            cout << "Q^H*Q[" << i << "][" << j << "] : "
                 << QHQre[i][j] << "  " << QHQim[i][j] << endl;
         if(i == j)
            error = error + abs(QHQre[i][j] - 1.0) + abs(QHQim[i][j]);
         else
            error = error + abs(QHQre[i][j]) + abs(QHQim[i][j]);
      }

   cout << scientific << setprecision(2);
   cout << "Sum of errors on |Q^H*Q - I| : " << error << endl;

   CPU_cmplx_factors_matmatmul
      (nrows,nrows,ncols,QHre,QHim,Are,Aim,QHAre,QHAim);

   error = 0.0;

   cout << scientific << setprecision(16);

   if(verbose > 0) cout << "The matrix transpose(Q)*A :" << endl;
   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
      {
         if(verbose > 0)
            cout << "Q^H*A[" << i << "][" << j << "] : "
                 << QHAre[i][j] << "  " << QHAim[i][j] << endl;
         error = error + abs(Rre[i][j] - QHAre[i][j])
                       + abs(Rim[i][j] - QHAim[i][j]);
      }

   cout << scientific << setprecision(2);
   cout << "Sum of errors on |Q^H*A - R| : " << error << endl;
}
