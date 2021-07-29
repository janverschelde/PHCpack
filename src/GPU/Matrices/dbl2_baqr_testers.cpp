// The file dbl2_baqr_testers.cpp defines the function with prototypes in
// the file dbl2_baqr_testers.h.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#include "random2_matrices.h"
#include "dbl2_factorizations.h"
#include "dbl2_factors_testers.h"
#include "dbl2_baqr_host.h"
#include "dbl2_baqr_kernels.h"

using namespace std;

void test_real2_blocked_qr ( void )
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

   double **Ahi = new double*[nrows];
   double **Alo = new double*[nrows];
   double **Qhi_h = new double*[nrows];
   double **Qlo_h = new double*[nrows];
   double **Qhi_d = new double*[nrows];
   double **Qlo_d = new double*[nrows];
   double **Rhi_h = new double*[nrows];
   double **Rlo_h = new double*[nrows];
   double **Rhi_d = new double*[nrows];
   double **Rlo_d = new double*[nrows];

   for(int i=0; i<nrows; i++)
   {
      Ahi[i] = new double[ncols];
      Alo[i] = new double[ncols];
      Qhi_h[i] = new double[nrows];
      Qlo_h[i] = new double[nrows];
      Qhi_d[i] = new double[nrows];
      Qlo_d[i] = new double[nrows];
      Rhi_h[i] = new double[ncols];
      Rlo_h[i] = new double[ncols];
      Rhi_d[i] = new double[ncols];
      Rlo_d[i] = new double[ncols];
   }
   random_dbl2_matrix(nrows,ncols,Ahi,Alo);

   if(verbose > 0)
   {
      cout << scientific << setprecision(16);

      cout << "A random matrix :" << endl;
      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++)
            cout << "A[" << i << "][" << j << "] : "
                 << Ahi[i][j] << "  " << Alo[i][j] << endl;
   }
   double timelapsed_h;
   bool vrb = (verbose > 0);

   cout << "-> CPU computes the block Householder QR ..." << endl;

   CPU_dbl2_blocked_houseqr
      (nrows,ncols,sizetile,numtiles,Ahi,Alo,Qhi_h,Qlo_h,Rhi_h,Rlo_h,
       &timelapsed_h,vrb);

   cout << "-> Testing the QR factorization ..." << endl;

   const double tol = 1.0e-26;
   int fail = test_real2_qr_factors
      (nrows,ncols,Ahi,Alo,Qhi_h,Qlo_h,Rhi_h,Rlo_h,tol,verbose);
   if(fail == 0)
      cout << "The test succeeded." << endl;
   else
   {
      cout << scientific << setprecision(2);
      cout << "The test failed for tol = " << tol << "." << endl;
   }
   double timelapsed_d;
   double houselapsedms,tileRlapsedms,vb2Wlapsedms;
   double WYTlapsedms,QWYTlapsedms,Qaddlapsedms;
   double YWTlapsedms,YWTClapsedms,Raddlapsedms;

   cout << "-> GPU computes the block Householder QR ..." << endl;

   GPU_dbl2_blocked_houseqr
      (nrows,ncols,sizetile,numtiles,Ahi,Alo,Qhi_d,Qlo_d,Rhi_d,Rlo_d,
       &houselapsedms,&tileRlapsedms,&vb2Wlapsedms,
       &WYTlapsedms,&QWYTlapsedms,&Qaddlapsedms,
       &YWTlapsedms,&YWTClapsedms,&Raddlapsedms,&timelapsed_d,vrb);

   fail = test_real2_qr_factors
             (nrows,ncols,Ahi,Alo,Qhi_d,Qlo_d,Rhi_d,Rlo_d,tol,verbose);
   if(fail == 0)
      cout << "The test succeeded." << endl;
   else
   {
      cout << scientific << setprecision(2);
      cout << "The test failed for tol = " << tol << "." << endl;
   }

   cout << fixed << setprecision(3);
   cout << "Elapsed CPU time (Linux), Wall time (Windows) : "
        << timelapsed_h << " seconds." << endl;
   cout << "         Time spent by the Householder kernel : ";
   cout << houselapsedms << " milliseconds." << endl;
   cout << "  Time spent by the kernel to reduce one tile : ";
   cout << tileRlapsedms << " milliseconds." << endl;
   cout << "    Time spent by the kernel for the W matrix : ";
   cout << vb2Wlapsedms << " milliseconds." << endl;
   cout << " Time spent by the kernel for computing W*Y^T : ";
   cout << WYTlapsedms << " milliseconds." << endl;
   cout << " Time spent by the kernel for computing Y*W^T : ";
   cout << YWTlapsedms << " milliseconds." << endl;
   cout << " Time spent by the kernel for computing Q*WYT : ";
   cout << QWYTlapsedms << " milliseconds." << endl;
   cout << " Time spent by the kernel for computing YWT*C : ";
   cout << YWTClapsedms << " milliseconds." << endl;
   cout << "Time spent by the kernel for adding QWYT to Q : ";
   cout << Qaddlapsedms << " milliseconds." << endl;
   cout << "Time spent by the kernel for adding R to YWTC : ";
   cout << Raddlapsedms << " milliseconds." << endl;
   cout << "        Total GPU wall clock computation time : ";
   cout << fixed << setprecision(3) << timelapsed_d << " seconds." << endl;

   for(int i=0; i<nrows; i++)
   {
      free(Ahi[i]); free(Alo[i]);
      free(Qhi_h[i]); free(Rhi_h[i]);
      free(Qhi_d[i]); free(Rhi_d[i]);
      free(Qlo_h[i]); free(Rlo_h[i]);
      free(Qlo_d[i]); free(Rlo_d[i]);
   }
   free(Ahi); free(Qhi_h); free(Rhi_h); free(Qhi_d); free(Rhi_d);
   free(Alo); free(Qlo_h); free(Rlo_h); free(Qlo_d); free(Rlo_d);
}

void test_cmplx2_blocked_qr ( void )
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

   double **Arehi = new double*[nrows];
   double **Arelo = new double*[nrows];
   double **Aimhi = new double*[nrows];
   double **Aimlo = new double*[nrows];
   double **Qrehi = new double*[nrows];
   double **Qrelo = new double*[nrows];
   double **Qimhi = new double*[nrows];
   double **Qimlo = new double*[nrows];
   double **Rrehi = new double*[nrows];
   double **Rrelo = new double*[nrows];
   double **Rimhi = new double*[nrows];
   double **Rimlo = new double*[nrows];

   for(int i=0; i<nrows; i++)
   {
      Arehi[i] = new double[ncols];
      Arelo[i] = new double[ncols];
      Aimhi[i] = new double[ncols];
      Aimlo[i] = new double[ncols];
      Qrehi[i] = new double[nrows];
      Qrelo[i] = new double[nrows];
      Qimhi[i] = new double[nrows];
      Qimlo[i] = new double[nrows];
      Rrehi[i] = new double[ncols];
      Rrelo[i] = new double[ncols];
      Rimhi[i] = new double[ncols];
      Rimlo[i] = new double[ncols];
   }
   random_cmplx2_matrix(nrows,ncols,Arehi,Arelo,Aimhi,Aimlo);

   if(verbose > 0)
   {
      cout << scientific << setprecision(16);

      cout << "A random matrix :" << endl;
      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++)
         {
            cout << "A[" << i << "][" << j << "]re : "
                 << Arehi[i][j] << "  " << Arelo[i][j] << endl;
            cout << "A[" << i << "][" << j << "]im : "
                 << Aimhi[i][j] << "  " << Aimlo[i][j] << endl;
         }
   }
   double timelapsed_h;
   bool vrb = (verbose > 0);

   cout << "-> Computed the block Householder QR ..." << endl;

   CPU_cmplx2_blocked_houseqr
      (nrows,ncols,sizetile,numtiles,
       Arehi,Arelo,Aimhi,Aimlo,
       Qrehi,Qrelo,Qimhi,Qimlo,
       Rrehi,Rrelo,Rimhi,Rimlo,&timelapsed_h,vrb);

   cout << "-> Testing the QR factorization ..." << endl;

   const double tol = 1.0e-26;
   const int fail = test_cmplx2_qr_factors
      (nrows,ncols,Arehi,Arelo,Aimhi,Aimlo,
                   Qrehi,Qrelo,Qimhi,Qimlo,
                   Rrehi,Rrelo,Rimhi,Rimlo,tol,verbose);
   if(fail == 0)
      cout << "The test succeeded." << endl;
   else
   {
      cout << scientific << setprecision(2);
      cout << "The test failed for tol = " << tol << "." << endl;
   }
   cout << fixed << setprecision(3);
   cout << "Elapsed CPU time (Linux), Wall time (Windows) : "
        << timelapsed_h << " seconds." << endl;

   for(int i=0; i<nrows; i++)
   {
      free(Arehi[i]); free(Qrehi[i]); free(Rrehi[i]);
      free(Arelo[i]); free(Qrelo[i]); free(Rrelo[i]);
      free(Aimhi[i]); free(Qimhi[i]); free(Rimhi[i]);
      free(Aimlo[i]); free(Qimlo[i]); free(Rimlo[i]);
   }
   free(Arehi); free(Qrehi); free(Rrehi);
   free(Arelo); free(Qrelo); free(Rrelo);
   free(Aimhi); free(Qimhi); free(Rimhi);
   free(Aimlo); free(Qimlo); free(Rimlo);
}
