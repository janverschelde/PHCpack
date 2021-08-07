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

void test_real_blocked_qr
 ( int seed, int szt, int nbt, int nrows, int vrb, int mode )
{
   const int sizetile = szt;
   const int numtiles = nbt;
   const int ncols = sizetile*numtiles;
   const int verbose = vrb;

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
   bool bvrb = (verbose > 0);
   double timelapsed_h;
   const double tol = 1.0e-8;
   int fail;

   if((mode == 1) || (mode == 2))
   {
      cout << "-> CPU computes the block Householder QR ..." << endl;

      CPU_dbl_blocked_houseqr
         (nrows,ncols,sizetile,numtiles,A,Q_h,R_h,&timelapsed_h,bvrb);

      cout << "-> Testing the QR factorization ..." << endl;

      fail = test_real_qr_factors(nrows,ncols,A,Q_h,R_h,tol,verbose);
      if(fail == 0)
         cout << "The test succeeded." << endl;
      else
      {
         cout << scientific << setprecision(2);
         cout << "The test failed for tol = " << tol << "." << endl;
      }
   }
   double timelapsed_d;
   double houselapsedms,tileRlapsedms,vb2Wlapsedms;
   double WYTlapsedms,QWYTlapsedms,Qaddlapsedms;
   double YWTlapsedms,YWTClapsedms,Raddlapsedms;
   long int addcnt = 0;
   long int mulcnt = 0;
   long int divcnt = 0;
   long int sqrtcnt = 0;

   if((mode == 0) || (mode == 2))
   {
      cout << "-> GPU computes the blocked Householder QR ..." << endl;

      GPU_dbl_blocked_houseqr
         (nrows,ncols,sizetile,numtiles,A,Q_d,R_d,
          &houselapsedms,&tileRlapsedms,&vb2Wlapsedms,
          &WYTlapsedms,&QWYTlapsedms,&Qaddlapsedms,
          &YWTlapsedms,&YWTClapsedms,&Raddlapsedms,&timelapsed_d,
          &addcnt,&mulcnt,&divcnt,&sqrtcnt,bvrb);
 
      fail = test_real_qr_factors(nrows,ncols,A,Q_d,R_d,tol,verbose);
      if(fail == 0)
         cout << "The test succeeded." << endl;
      else
      {
         cout << scientific << setprecision(2);
         cout << "The test failed for tol = " << tol << "." << endl;
      }
   }
   cout << endl;
   cout << fixed << setprecision(3);
   if(mode == 1)
   {
      cout << "Elapsed CPU time (Linux), Wall time (Windows) : "
           << timelapsed_h << " seconds." << endl;
   }
   if(mode == 0)
   {
      cout << "         Time spent by the Householder kernel : ";
      cout << houselapsedms << " milliseconds." << endl;
      cout << "  Time spent by the kernel to reduce one tile : ";
      cout << tileRlapsedms << " milliseconds." << endl;
      cout << "    Time spent by the kernel for the W matrix : ";
      cout << vb2Wlapsedms << " milliseconds." << endl;
      // cout << " Time spent by the kernel for computing W*Y^T : ";
      // cout << WYTlapsedms << " milliseconds." << endl;
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
      cout << endl;
      cout << "             Number of additions/subtractions : "
           << addcnt << endl;
      cout << "                    Number of multiplications : "
           << mulcnt << endl;
      cout << "                          Number of divisions : "
           << divcnt << endl;
      cout << "                    Number of calls to sqrt() : "
           << sqrtcnt << endl;
      long int flopcnt = addcnt + mulcnt + divcnt + sqrtcnt;
      cout << "    Total number of floating-point operations : "
           << flopcnt << endl;
   }
   for(int i=0; i<nrows; i++)
   {
      free(A[i]); free(Q_h[i]); free(Q_d[i]); free(R_h[i]); free(R_d[i]);
   }
   free(A); free(Q_h); free(Q_d); free(R_h); free(R_d);
}

void test_cmplx_blocked_qr
 ( int seed, int szt, int nbt, int nrows, int vrb, int mode )
{
   const int sizetile = szt;
   const int numtiles = nbt;
   const int ncols = sizetile*numtiles;
   const int verbose = vrb;

   cout << "-> Generating a random " << nrows
        << "-by-" << ncols << " matrix ..." << endl;

   double **Are = new double*[nrows];
   double **Aim = new double*[nrows];
   double **Qre_h = new double*[nrows];
   double **Qre_d = new double*[nrows];
   double **Qim_h = new double*[nrows];
   double **Qim_d = new double*[nrows];
   double **Rre_h = new double*[nrows];
   double **Rre_d = new double*[nrows];
   double **Rim_h = new double*[nrows];
   double **Rim_d = new double*[nrows];

   for(int i=0; i<nrows; i++)
   {
      Are[i] = new double[ncols];
      Aim[i] = new double[ncols];
      Qre_h[i] = new double[nrows];
      Qre_d[i] = new double[nrows];
      Qim_h[i] = new double[nrows];
      Qim_d[i] = new double[nrows];
      Rre_h[i] = new double[ncols];
      Rre_d[i] = new double[ncols];
      Rim_h[i] = new double[ncols];
      Rim_d[i] = new double[ncols];
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
   bool bvrb = (verbose > 0);

   cout << "-> CPU computes the block Householder QR ..." << endl;

   CPU_cmplx_blocked_houseqr
      (nrows,ncols,sizetile,numtiles,Are,Aim,Qre_h,Qim_h,Rre_h,Rim_h,
       &timelapsed_h,bvrb);

   cout << "-> Testing the QR factorization ..." << endl;

   const double tol = 1.0e-8;
   int fail = test_cmplx_qr_factors
      (nrows,ncols,Are,Aim,Qre_h,Qim_h,Rre_h,Rim_h,tol,verbose);
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

   cout << "-> GPU computes the blocked Householder QR ..." << endl;

   GPU_cmplx_blocked_houseqr
      (nrows,ncols,sizetile,numtiles,Are,Aim,Qre_d,Qim_d,Rre_d,Rim_d,
       &houselapsedms,&tileRlapsedms,&vb2Wlapsedms,
       &WYTlapsedms,&QWYTlapsedms,&Qaddlapsedms,
       &YWTlapsedms,&YWTClapsedms,&Raddlapsedms,&timelapsed_d,bvrb);

   fail = test_cmplx_qr_factors
             (nrows,ncols,Are,Aim,Qre_d,Qim_d,Rre_d,Rim_d,tol,verbose);
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
   // cout << " Time spent by the kernel for computing W*Y^T : ";
   // cout << WYTlapsedms << " milliseconds." << endl;
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
      free(Are[i]); free(Qre_h[i]); free(Rre_h[i]);
                    free(Qre_d[i]); free(Rre_d[i]);
      free(Aim[i]); free(Qim_h[i]); free(Rim_h[i]);
                    free(Qim_d[i]); free(Rim_d[i]);
   }
   free(Are); free(Qre_h); free(Qre_d); free(Rre_h); free(Rre_d);
   free(Aim); free(Qim_h); free(Qim_d); free(Rim_h); free(Rim_d);
}
