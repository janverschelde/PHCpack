/* Tests if the single indexed matrix matrix multiplication
 * corresponds to the double indexed matrix matrix multiplication. */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include "random_matrices.h"
#include "double_matrix_multiplications.h"

using namespace std;

int main ( void )
{
   srand(time(NULL));

   cout << "Give #rows of the product A*B : ";
   int m; cin >> m;
   cout << "Give #columns of the product A*B : ";
   int n; cin >> n;
   cout << "Give #columns of A, #rows of B : ";
   int k; cin >> k;

   double **CC = new double*[m];
   for(int i=0; i<m; i++) CC[i] = new double[n];

   double **AA = new double*[m];
   for(int i=0; i<m; i++) AA[i] = new double[k];

   double **BB = new double*[k];
   for(int i=0; i<k; i++) BB[i] = new double[n];

   random_dbl_matrix(m, k, AA);

   double *A = new double[m*k];
   double2single_row_major(m, k, AA, A);

   cout << scientific << setprecision(16);

   double error = 0.0;

   cout << "A random " << n << "-by-" << k << " matrix A :" << endl;
   for(int i=0, idx=0; i<m; i++)
      for(int j=0; j<k; j++)
      {
         cout << "AA[" << i << "][" << j << "] : " << AA[i][j] << endl
              << " A[" << idx << "]    : " << A[idx] << endl;
         error = error + abs(AA[i][j] - A[idx++]);
      }

   cout << scientific << setprecision(3)
        << "sum of errors : " << error << endl;

   random_dbl_matrix(k, n, BB);

   cout << scientific << setprecision(16);

   cout << "A random " << k << "-by-" << n << " matrix B :" << endl;
   for(int i=0; i<k; i++)
      for(int j=0; j<n; j++)
         cout << "BB[" << i << "][" << j << "] : " << BB[i][j] << endl;

   double_indexed_matrix_multiplication(m, n, k, AA, BB, CC);

   double **TT = new double*[n];
   for(int i=0; i<n; i++) TT[i] = new double[k];
   transpose_rows_columns(k, n, BB, TT);

   cout << "Transpose of the matrix B :" << endl;
   for(int i=0; i<n; i++)
      for(int j=0; j<k; j++)
         cout << "TT[" << i << "][" << j << "] : " << TT[i][j] << endl;

   double *B = new double[k*n];
   double2single_column_major(k, n, TT, B);

   cout << scientific << setprecision(16);

   error = 0.0;

   cout << "A random " << n << "-by-" << k << " matrix TT :" << endl;
   for(int i=0, idx=0; i<n; i++)
      for(int j=0; j<k; j++)
      {
         cout << "TT[" << i << "][" << j << "] : " << TT[i][j] << endl
              << " B[" << idx << "]    : " << B[idx] << endl;
         error = error + abs(TT[i][j] - B[idx++]);
      }

   cout << scientific << setprecision(3)
        << "sum of errors : " << error << endl;

   double *C = new double[m*n];
   single_indexed_matrix_multiplication(m, n, k, A, B, C);

   cout << scientific << setprecision(16);

   error = 0.0;
   cout << "comparing matrix products ..." << endl;
   for(int i=0, idx=0; i<m; i++)
      for(int j=0; j<n; j++)
      {
         cout << "CC[" << i << "][" << j << "] : " << CC[i][j] << endl
              << " C[" << idx << "]    : " << C[idx] << endl;
         error = error + abs(CC[i][j] - C[idx++]);
      }

   cout << scientific << setprecision(3)
        << "sum of errors : " << error << endl;

   return 0;
}
