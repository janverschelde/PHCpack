// The file dbl_baqr_testers.cpp defines the function with prototypes in
// the file dbl_baqr_testers.h.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#include "random_matrices.h"
#include "dbl_factorizations.h"
#include "dbl_factors_testers.h"
#include "dbl_baqr_host.h"
#include "dbl_baqr_kernels.h"

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

   cout << "-> Generating a random " << nrows
        << "-by-" << ncols << " matrix ..." << endl;

   double **A = new double*[nrows];
   double **Q_h = new double*[nrows];
   double **Q_d = new double*[nrows];
   double **R_h = new double*[nrows];
   double **R_d = new double*[nrows];

   for(int i=0; i<nrows; i++)
   {
      A[i] = new double[ncols];
      Q_h[i] = new double[nrows];
      Q_d[i] = new double[nrows];
      R_h[i] = new double[ncols];
      R_d[i] = new double[ncols];
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
   double timelapsed_h,timelapsed_d;
   double houselapsedms,tileRlapsedms,vb2Wlapsedms;

   cout << "-> CPU computes the block Householder QR ..." << endl;

   CPU_dbl_blocked_houseqr
      (nrows,ncols,sizetile,numtiles,A,Q_h,R_h,&timelapsed_h,vrb);

   cout << "-> Testing the QR factorization ..." << endl;

   test_real_qr_factors(nrows,ncols,A,Q_h,R_h,verbose);

   cout << "-> GPU computes the block Householder QR ..." << endl;

   GPU_dbl_blocked_houseqr
      (nrows,ncols,sizetile,numtiles,A,Q_d,R_d,
       &houselapsedms,&tileRlapsedms,&vb2Wlapsedms,&timelapsed_d,vrb);

   cout << fixed << setprecision(3);
   cout << "Elapsed CPU time (Linux), Wall time (Windows) : "
        << timelapsed_h << " seconds." << endl;
   cout << "         Time spent by the Householder kernel : ";
   cout << houselapsedms << " milliseconds." << endl;
   cout << "  Time spent by the kernel to reduce one tile : ";
   cout << tileRlapsedms << " milliseconds." << endl;
   cout << "    Time spent by the kernel for the W matrix : ";
   cout << vb2Wlapsedms << " milliseconds." << endl;
   cout << "        Total GPU wall clock computation time : ";
   cout << fixed << setprecision(3) << timelapsed_d << " seconds." << endl;

   for(int i=0; i<nrows; i++)
   {
      free(A[i]); free(Q_h[i]); free(Q_d[i]); free(R_h[i]); free(R_d[i]);
   }
   free(A); free(Q_h); free(Q_d); free(R_h); free(R_d);
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

   cout << "-> Generating a random " << nrows
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
   double timelapsed_h;
   bool vrb = (verbose > 0);

   cout << "-> CPU computes the block Householder QR ..." << endl;

   CPU_cmplx_blocked_houseqr
      (nrows,ncols,sizetile,numtiles,Are,Aim,Qre,Qim,Rre,Rim,
       &timelapsed_h,vrb);

   cout << "-> Testing the QR factorization ..." << endl;

   test_cmplx_qr_factors(nrows,ncols,Are,Aim,Qre,Qim,Rre,Rim,verbose);

   cout << fixed << setprecision(3);
   cout << "Elapsed CPU time (Linux), Wall time (Windows) : "
        << timelapsed_h << " seconds." << endl;

   for(int i=0; i<nrows; i++)
   {
      free(Are[i]); free(Qre[i]); free(Rre[i]);
      free(Aim[i]); free(Qim[i]); free(Rim[i]);
   }
   free(Are); free(Qre); free(Rre);
   free(Aim); free(Qim); free(Rim);
}
