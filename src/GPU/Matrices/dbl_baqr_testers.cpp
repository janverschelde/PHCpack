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
   cout << "Give the number of rows : ";
   int nrows; cin >> nrows;

   cout << "Give the size of each tile : ";
   int sizetile; cin >> sizetile;

   cout << "Give the number of tiles : ";
   int numtiles; cin >> numtiles;

   cout << "Give the verbose level (1 to see all numbers) : ";
   int verbose; cin >> verbose;

   const int ncols = sizetile*numtiles;

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
   CPU_dbl_blocked_houseqr(nrows,ncols,sizetile,numtiles,A,Q,R);

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
