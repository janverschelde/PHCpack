// The file dbl_qrbs_testers.cpp defines the function with prototypes in
// the file dbl_qrbs_testers.h.

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
#include "dbl_tabs_host.h"
#include "dbl_tabs_kernels.h"
#include "dbl_test_utilities.h"

using namespace std;

void test_real_blocked_qrbs
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
   double *sol = new double[ncols];
   for(int i=0; i<ncols; i++) sol[i] = 1.0;

   double *rhs = new double[nrows];
   for(int i=0; i<nrows; i++)
   {
      rhs[i] = 0.0;
      for(int j=0; j<ncols; j++) rhs[i] = rhs[i] + A[i][j]*sol[j];
   }
   if(verbose > 0)
   {
      cout << "The sums of the columns :" << endl;
      for(int i=0; i<nrows; i++)
         cout << "b[" << i << "] : " << rhs[i] << endl;
   }
   bool bvrb = (verbose > 0);
   double qrtimelapsed_h,bstimelapsed_h;
   const double tol = 1.0e-8;
   int fail;

   if((mode == 1) || (mode == 2))
   {
      double *x = new double[ncols];
      double *qTrhs = new double[nrows];

      cout << "-> CPU computes the block Householder QR ..." << endl;

      CPU_dbl_blocked_houseqr
         (nrows,ncols,sizetile,numtiles,A,Q_h,R_h,&qrtimelapsed_h,bvrb);

      cout << "-> Testing the QR factorization ..." << endl;

      // fail = test_real_qr_factors(nrows,ncols,A,Q_h,R_h,tol,verbose);
      fail = test_real_qr_factors_probe(nrows,ncols,A,Q_h,R_h,tol,2,1);
      if(fail == 0)
         cout << "The test succeeded." << endl;
      else
      {
         cout << scientific << setprecision(2);
         cout << "The test failed for tol = " << tol << "." << endl;
      }
      cout << "-> CPU multiplies Q^T with b ..." << endl;

      for(int i=0; i<nrows; i++)
      {
         qTrhs[i] = 0.0;
         for(int j=0; j<nrows; j++)
            qTrhs[i] = qTrhs[i] + Q_h[j][i]*rhs[j];
      }
      cout << "-> CPU solves an upper triangular system ..." << endl;

      CPU_dbl_upper_tiled_solver
         (ncols,sizetile,numtiles,R_h,qTrhs,x,&bstimelapsed_h);

      if(verbose > 0)
      {
         cout << "CPU solution computed with tiling :" << endl;
         cout << scientific << setprecision(16);
         for(int i=0; i<ncols; i++)
            cout << "x[" << i << "] : " << x[i] << endl;
      }
      cout << scientific << setprecision(2);
      cout << "   Sum of CPU errors on solution : "
           << dbl_Difference_Sum(ncols,sol,x) << endl;

      free(x); free(qTrhs);
   }
   double qrtimelapsed_d,bstimelapsed_d;
   double houselapsedms,RTvlapsedms,tileRlapsedms,vb2Wlapsedms;
   double WYTlapsedms,QWYTlapsedms,Qaddlapsedms;
   double YWTlapsedms,YWTClapsedms,Raddlapsedms;
   double invlapsed,mullapsed,sublapsed,bselapsedms;
   long long int qraddcnt = 0;
   long long int qrmulcnt = 0;
   long long int qrdivcnt = 0;
   long long int sqrtcnt = 0;
   long long int bsaddcnt = 0;
   long long int bsmulcnt = 0;
   long long int bsdivcnt = 0;

   if((mode == 0) || (mode == 2))
   {
      double *x_d = new double[ncols];
      double *qTrhs_d = new double[nrows];

      cout << "-> GPU computes the blocked Householder QR ..." << endl;

      GPU_dbl_blocked_houseqr
         (nrows,ncols,sizetile,numtiles,A,Q_d,R_d,
          &houselapsedms,&RTvlapsedms,&tileRlapsedms,&vb2Wlapsedms,
          &WYTlapsedms,&QWYTlapsedms,&Qaddlapsedms,
          &YWTlapsedms,&YWTClapsedms,&Raddlapsedms,&qrtimelapsed_d,
          &qraddcnt,&qrmulcnt,&qrdivcnt,&sqrtcnt,bvrb);

      cout << "-> Testing the QR factorization ..." << endl;
 
      // fail = test_real_qr_factors(nrows,ncols,A,Q_d,R_d,tol,verbose);
      fail = test_real_qr_factors_probe(nrows,ncols,A,Q_d,R_d,tol,2,1);
      if(fail == 0)
         cout << "The test succeeded." << endl;
      else
      {
         cout << scientific << setprecision(2);
         cout << "The test failed for tol = " << tol << "." << endl;
      }
      // preliminary CPU computation, or for testing purposes only
      cout << "-> CPU multiplies Q^T with b ..." << endl;

      for(int i=0; i<nrows; i++)
      {
         qTrhs_d[i] = 0.0;
         for(int j=0; j<nrows; j++)
            qTrhs_d[i] = qTrhs_d[i] + Q_d[j][i]*rhs[j];
      }
      cout << "-> GPU solves an upper triangular system ..." << endl;

      GPU_dbl_upper_tiled_solver
         (ncols,sizetile,numtiles,R_d,qTrhs_d,x_d,
          &invlapsed,&mullapsed,&sublapsed,&bselapsedms,&bstimelapsed_d,
          &bsaddcnt,&bsmulcnt,&bsdivcnt);

      cout << scientific << setprecision(2);
      cout << "   Sum of GPU errors on solution : "
           << dbl_Difference_Sum(ncols,sol,x_d) << endl;

      free(x_d); free(qTrhs_d);
   }
   cout << endl;
   cout << fixed << setprecision(3);
   if((mode == 1) || (mode == 2))
   {
      cout << "QR Elapsed CPU time (Linux), Wall time (Windows) : "
           << qrtimelapsed_h << " seconds." << endl;
      cout << "BS Elapsed CPU time (Linux), Wall time (Windows) : "
           << bstimelapsed_h << " seconds." << endl;
   }
   if((mode == 0) || (mode == 2))
   {
      cout << "         Time spent by the Householder kernel : "
           << houselapsedms << " milliseconds." << endl;
      cout << "      Time spent by the kernel for beta*R^T*v : "
           << RTvlapsedms << " milliseconds." << endl;
      cout << "  Time spent by the kernel to reduce one tile : "
           << tileRlapsedms << " milliseconds." << endl;
      cout << "    Time spent by the kernel for the W matrix : "
           << vb2Wlapsedms << " milliseconds." << endl;
      // cout << " Time spent by the kernel for computing W*Y^T : ";
      // cout << WYTlapsedms << " milliseconds." << endl;
      cout << " Time spent by the kernel for computing Y*W^T : "
           << YWTlapsedms << " milliseconds." << endl;
      cout << " Time spent by the kernel for computing Q*WYT : "
           << QWYTlapsedms << " milliseconds." << endl;
      cout << " Time spent by the kernel for computing YWT*C : "
           << YWTClapsedms << " milliseconds." << endl;
      cout << "Time spent by the kernel for adding QWYT to Q : "
           << Qaddlapsedms << " milliseconds." << endl;
      cout << "Time spent by the kernel for adding R to YWTC : "
           << Raddlapsedms << " milliseconds." << endl;
      const double totlapsedms = houselapsedms + RTvlapsedms
         + tileRlapsedms + vb2Wlapsedms + YWTlapsedms + QWYTlapsedms
         + YWTClapsedms + Qaddlapsedms + Raddlapsedms;
      cout << "                    Time spent by all kernels : "
           << totlapsedms << " milliseconds." << endl;
      cout << "     Total QR GPU wall clock computation time : ";
      cout << fixed << setprecision(3)
           << qrtimelapsed_d << " seconds." << endl;
      cout << endl;
      cout << "             Number of additions/subtractions : "
           << qraddcnt << endl;
      cout << "                    Number of multiplications : "
           << qrmulcnt << endl;
      cout << "                          Number of divisions : "
           << qrdivcnt << endl;
      cout << "                    Number of calls to sqrt() : "
           << sqrtcnt << endl;
      long long int qrflopcnt = qraddcnt + qrmulcnt + qrdivcnt + sqrtcnt;
      cout << "    Total number of floating-point operations : "
           << qrflopcnt << endl;
      cout << endl;
      double qrkernflops = 1000.0*((double) qrflopcnt)/totlapsedms;
      double qrwallflops = ((double) qrflopcnt)/qrtimelapsed_d;
      const int gigacnt = pow(2.0,30);
      cout << "QR Kernel Time Flops : "
           << scientific << setprecision(3) << qrkernflops;
      cout << fixed << setprecision(3)
           << " = " << qrkernflops/gigacnt << " Gigaflops" << endl;
      cout << " QR Wall Clock Flops : "
           << scientific << setprecision(3) << qrwallflops;
      cout << fixed << setprecision(3)
           << " = " << qrwallflops/gigacnt << " Gigaflops" << endl;
      cout << endl;
      cout << "          Time spent to invert diagonal tiles : ";
      cout << invlapsed << " milliseconds." << endl;
      cout << "   Time spent to multiply with inverted tiles : ";
      cout << mullapsed << " milliseconds." << endl;
      cout << "             Time spent for back substitution : ";
      cout << sublapsed << " milliseconds." << endl;
      cout << "                    Time spent by all kernels : ";
      cout << bselapsedms << " milliseconds." << endl;
      cout << "     Total BS GPU wall clock computation time : ";
      cout << fixed << setprecision(3)
           << bstimelapsed_d << " seconds." << endl;
      cout << endl;
      cout << "             Number of additions/subtractions : "
           << bsaddcnt << endl;
      cout << "                    Number of multiplications : "
           << bsmulcnt << endl;
      cout << "                          Number of divisions : "
           << bsdivcnt << endl;
      long long int bsflopcnt = bsaddcnt + bsmulcnt + bsdivcnt;
      cout << "    Total number of floating-point operations : "
           << bsflopcnt << endl;
      cout << endl;
      double bskernflops = 1000.0*((double) bsflopcnt)/bselapsedms;
      double bswallflops = ((double) bsflopcnt)/bstimelapsed_d;
      // const int gigacnt = pow(2.0,30);
      cout << "BS Kernel Time Flops : "
           << scientific << setprecision(3) << bskernflops;
      cout << fixed << setprecision(3)
           << " = " << bskernflops/gigacnt << " Gigaflops" << endl;
      cout << " BS Wall Clock Flops : "
           << scientific << setprecision(3) << bswallflops;
      cout << fixed << setprecision(3)
           << " = " << bswallflops/gigacnt << " Gigaflops" << endl;
      cout << endl;
      long long int totalflopcnt = qrflopcnt + bsflopcnt;
      double totalelapsedms = totlapsedms + bselapsedms;
      double totaltimelapsed = qrtimelapsed_d + bstimelapsed_d;
      double totalkernflops = 1000.0*((double) totalflopcnt)/totalelapsedms;
      double totalwallflops = ((double) totalflopcnt)/totaltimelapsed;
      // const int gigacnt = pow(2.0,30);
      cout << "Total Kernel Time Flops : "
           << scientific << setprecision(3) << totalkernflops;
      cout << fixed << setprecision(3)
           << " = " << totalkernflops/gigacnt << " Gigaflops" << endl;
      cout << " Total Wall Clock Flops : "
           << scientific << setprecision(3) << totalwallflops;
      cout << fixed << setprecision(3)
           << " = " << totalwallflops/gigacnt << " Gigaflops" << endl;
   }
   for(int i=0; i<nrows; i++)
   {
      free(A[i]); free(Q_h[i]); free(Q_d[i]); free(R_h[i]); free(R_d[i]);
   }
   free(A); free(Q_h); free(Q_d); free(R_h); free(R_d);
}

void test_cmplx_blocked_qrbs
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
   double *solre = new double[ncols];
   double *solim = new double[ncols];
   for(int i=0; i<ncols; i++)
   {
      solre[i] = 1.0;
      solim[i] = 0.0;
   }
   double *rhsre = new double[nrows];
   double *rhsim = new double[nrows];
   double accre,accim;

   for(int i=0; i<nrows; i++)
   {
      rhsre[i] = 0.0;
      rhsim[i] = 0.0;
      for(int j=0; j<ncols; j++) // rhs[i] = rhs[i] + A[i][j]*sol[j];
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
      for(int i=0; i<nrows; i++)
         cout << "b[" << i << "] : "
              << rhsre[i] << "  " << rhsim[i] << endl;
   }
   double qrtimelapsed_h,bstimelapsed_h;
   bool bvrb = (verbose > 0);
   const double tol = 1.0e-8;
   int fail;

   if((mode == 1) || (mode == 2))
   {
      double *xre = new double[ncols];
      double *xim = new double[ncols];
      double *qHrhsre = new double[nrows];
      double *qHrhsim = new double[nrows];

      cout << "-> CPU computes the block Householder QR ..." << endl;

      CPU_cmplx_blocked_houseqr
         (nrows,ncols,sizetile,numtiles,Are,Aim,Qre_h,Qim_h,Rre_h,Rim_h,
          &qrtimelapsed_h,bvrb);

      cout << "-> Testing the QR factorization ..." << endl;

      // fail = test_cmplx_qr_factors
      //   (nrows,ncols,Are,Aim,Qre_h,Qim_h,Rre_h,Rim_h,tol,verbose);
      fail = test_cmplx_qr_factors_probe
         (nrows,ncols,Are,Aim,Qre_h,Qim_h,Rre_h,Rim_h,tol,2,1);
      if(fail == 0)
         cout << "The test succeeded." << endl;
      else
      {
         cout << scientific << setprecision(2);
         cout << "The test failed for tol = " << tol << "." << endl;
      }
      cout << "-> CPU multiplies Q^H with b ..." << endl;

      for(int i=0; i<nrows; i++)
      {
         qHrhsre[i] = 0.0;
         qHrhsim[i] = 0.0;
         for(int j=0; j<nrows; j++) // rhs[i] = rhs[i] + Q[j][i]*sol[j];
         {
            accre =   Qre_h[j][i]*rhsre[j] + Qim_h[j][i]*rhsim[j];
            accim = - Qim_h[j][i]*rhsre[j] + Qre_h[j][i]*rhsim[j];
            qHrhsre[i] = qHrhsre[i] + accre;
            qHrhsim[i] = qHrhsim[i] + accim;
         }
      }
      cout << "-> CPU solves an upper triangular system ..." << endl;

      CPU_cmplx_upper_tiled_solver
         (ncols,sizetile,numtiles,Rre_h,Rim_h,qHrhsre,qHrhsim,xre,xim,
          &bstimelapsed_h);

      if(verbose > 0)
      {
         cout << "CPU solution computed with tiling :" << endl;
         cout << scientific << setprecision(16);
         for(int i=0; i<ncols; i++)
            cout << "x[" << i << "] : " << xre[i] << "  " << xim[i] << endl;
      }
      cout << scientific << setprecision(2);
      cout << "   Sum of CPU errors on solution : "
           << cmplx_Difference_Sum(ncols,solre,solim,xre,xim) << endl;

      free(xre); free(xim); free(qHrhsre); free(qHrhsim);
   }
   double qrtimelapsed_d,bstimelapsed_d;
   double houselapsedms,RHvlapsedms,tileRlapsedms,vb2Wlapsedms;
   double WYTlapsedms,QWYTlapsedms,Qaddlapsedms;
   double YWTlapsedms,YWTClapsedms,Raddlapsedms;
   double invlapsed,mullapsed,sublapsed,bselapsedms;
   long long int qraddcnt = 0;
   long long int qrmulcnt = 0;
   long long int qrdivcnt = 0;
   long long int sqrtcnt = 0;
   long long int bsaddcnt = 0;
   long long int bsmulcnt = 0;
   long long int bsdivcnt = 0;

   if((mode == 0) || (mode == 2))
   {
      double *xre_d = new double[ncols];
      double *xim_d = new double[ncols];
      double *qHrhsre_d = new double[nrows];
      double *qHrhsim_d = new double[nrows];

      cout << "-> GPU computes the blocked Householder QR ..." << endl;

      GPU_cmplx_blocked_houseqr
         (nrows,ncols,sizetile,numtiles,Are,Aim,Qre_d,Qim_d,Rre_d,Rim_d,
          &houselapsedms,&RHvlapsedms,&tileRlapsedms,&vb2Wlapsedms,
          &WYTlapsedms,&QWYTlapsedms,&Qaddlapsedms,
          &YWTlapsedms,&YWTClapsedms,&Raddlapsedms,&qrtimelapsed_d,
          &qraddcnt,&qrmulcnt,&qrdivcnt,&sqrtcnt,bvrb);

      cout << "-> Testing the QR factorization ..." << endl;

      // fail = test_cmplx_qr_factors
      //           (nrows,ncols,Are,Aim,Qre_d,Qim_d,Rre_d,Rim_d,tol,verbose);
      fail = test_cmplx_qr_factors_probe
                (nrows,ncols,Are,Aim,Qre_d,Qim_d,Rre_d,Rim_d,tol,2,1);
      if(fail == 0)
         cout << "The test succeeded." << endl;
      else
      {
         cout << scientific << setprecision(2);
         cout << "The test failed for tol = " << tol << "." << endl;
      }
      // preliminary, for testing purposes only
      cout << "-> CPU multiplies Q^H with b ..." << endl;

      for(int i=0; i<nrows; i++)
      {
         qHrhsre_d[i] = 0.0;
         qHrhsim_d[i] = 0.0;
         for(int j=0; j<nrows; j++) // rhs[i] = rhs[i] + Q[j][i]*sol[j];
         {
            accre =   Qre_d[j][i]*rhsre[j] + Qim_d[j][i]*rhsim[j];
            accim = - Qim_d[j][i]*rhsre[j] + Qre_d[j][i]*rhsim[j];
            qHrhsre_d[i] = qHrhsre_d[i] + accre;
            qHrhsim_d[i] = qHrhsim_d[i] + accim;
         }
      }
      cout << "-> GPU solves an upper triangular system ..." << endl;

      GPU_cmplx_upper_tiled_solver
         (ncols,sizetile,numtiles,Rre_d,Rim_d,qHrhsre_d,qHrhsim_d,xre_d,xim_d,
          &invlapsed,&mullapsed,&sublapsed,&bselapsedms,&bstimelapsed_d,
          &bsaddcnt,&bsmulcnt,&bsdivcnt);

      if(verbose > 0)
      {
         cout << "GPU solution computed with tiling :" << endl;
         for(int i=0; i<ncols; i++)
            cout << "x[" << i << "] : "
                 << xre_d[i] << "  " << xim_d[i] << endl;
      }
      cout << scientific << setprecision(2);
      cout << "   Sum of GPU errors on solution : "
           << cmplx_Difference_Sum(ncols,solre,solim,xre_d,xim_d) << endl;
   }
   cout << endl;
   cout << fixed << setprecision(3);
   if((mode == 1) || (mode == 2))
   {
      cout << "QR Elapsed CPU time (Linux), Wall time (Windows) : "
           << qrtimelapsed_h << " seconds." << endl;
      cout << "BS Elapsed CPU time (Linux), Wall time (Windows) : "
           << bstimelapsed_h << " seconds." << endl;
   }
   if((mode == 0) || (mode == 2))
   {
      cout << "         Time spent by the Householder kernel : "
           << houselapsedms << " milliseconds." << endl;
      cout << "      Time spent by the kernel for beta*R^H*v : "
           << RHvlapsedms << " milliseconds." << endl;
      cout << "  Time spent by the kernel to reduce one tile : "
           << tileRlapsedms << " milliseconds." << endl;
      cout << "    Time spent by the kernel for the W matrix : "
           << vb2Wlapsedms << " milliseconds." << endl;
      // cout << " Time spent by the kernel for computing W*Y^H : ";
      // cout << WYTlapsedms << " milliseconds." << endl;
      cout << " Time spent by the kernel for computing Y*W^H : "
           << YWTlapsedms << " milliseconds." << endl;
      cout << " Time spent by the kernel for computing Q*WYH : "
           << QWYTlapsedms << " milliseconds." << endl;
      cout << " Time spent by the kernel for computing YWH*C : "
           << YWTClapsedms << " milliseconds." << endl;
      cout << "Time spent by the kernel for adding QWYH to Q : "
           << Qaddlapsedms << " milliseconds." << endl;
      cout << "Time spent by the kernel for adding R to YWHC : "
           << Raddlapsedms << " milliseconds." << endl;
      const double totlapsedms = houselapsedms + RHvlapsedms
         + tileRlapsedms + vb2Wlapsedms + YWTlapsedms + QWYTlapsedms
         + YWTClapsedms + Qaddlapsedms + Raddlapsedms;
      cout << "                    Time spent by all kernels : "
           << totlapsedms << " milliseconds." << endl;
      cout << "     Total QR GPU wall clock computation time : ";
      cout << fixed << setprecision(3)
           << qrtimelapsed_d << " seconds." << endl;
      cout << endl;
      cout << "             Number of additions/subtractions : "
           << qraddcnt << endl;
      cout << "                    Number of multiplications : "
           << qrmulcnt << endl;
      cout << "                          Number of divisions : "
           << qrdivcnt << endl;
      cout << "                    Number of calls to sqrt() : "
           << sqrtcnt << endl;
      long long int qrflopcnt = qraddcnt + qrmulcnt + qrdivcnt + sqrtcnt;
      cout << "    Total number of floating-point operations : "
           << qrflopcnt << endl;
      cout << endl;
      double qrkernflops = 1000.0*((double) qrflopcnt)/totlapsedms;
      double qrwallflops = ((double) qrflopcnt)/qrtimelapsed_d;
      const int gigacnt = pow(2.0,30);
      cout << "QR Kernel Time Flops : "
           << scientific << setprecision(3) << qrkernflops;
      cout << fixed << setprecision(3)
           << " = " << qrkernflops/gigacnt << " Gigaflops" << endl;
      cout << " QR Wall Clock Flops : "
           << scientific << setprecision(3) << qrwallflops;
      cout << fixed << setprecision(3)
           << " = " << qrwallflops/gigacnt << " Gigaflops" << endl;
      cout << endl;
      cout << "          Time spent to invert diagonal tiles : ";
      cout << invlapsed << " milliseconds." << endl;
      cout << "   Time spent to multiply with inverted tiles : ";
      cout << mullapsed << " milliseconds." << endl;
      cout << "             Time spent for back substitution : ";
      cout << sublapsed << " milliseconds." << endl;
      cout << "                    Time spent by all kernels : ";
      cout << bselapsedms << " milliseconds." << endl;
      cout << "     Total BS GPU wall clock computation time : ";
      cout << fixed << setprecision(3) << bstimelapsed_d
           << " seconds." << endl;
      cout << endl;
      cout << "             Number of additions/subtractions : "
           << bsaddcnt << endl;
      cout << "                    Number of multiplications : "
           << bsmulcnt << endl;
      cout << "                          Number of divisions : "
           << bsdivcnt << endl;
      long long int bsflopcnt = bsaddcnt + bsmulcnt + bsdivcnt;
      cout << "    Total number of floating-point operations : "
           << bsflopcnt << endl;
      cout << endl;
      double bskernflops = 1000.0*((double) bsflopcnt)/bselapsedms;
      double bswallflops = ((double) bsflopcnt)/bstimelapsed_d;
      // const int gigacnt = pow(2.0,30);
      cout << "BS Kernel Time Flops : "
           << scientific << setprecision(3) << bskernflops;
      cout << fixed << setprecision(3)
           << " = " << bskernflops/gigacnt << " Gigaflops" << endl;
      cout << " BS Wall Clock Flops : "
           << scientific << setprecision(3) << bswallflops;
      cout << fixed << setprecision(3)
           << " = " << bswallflops/gigacnt << " Gigaflops" << endl;
      cout << endl;
      long long int totalflopcnt = qrflopcnt + bsflopcnt;
      double totalelapsedms = totlapsedms + bselapsedms;
      double totaltimelapsed = qrtimelapsed_d + bstimelapsed_d;
      double totalkernflops = 1000.0*((double) totalflopcnt)/totalelapsedms;
      double totalwallflops = ((double) totalflopcnt)/totaltimelapsed;
      // const int gigacnt = pow(2.0,30);
      cout << "Total Kernel Time Flops : "
           << scientific << setprecision(3) << totalkernflops;
      cout << fixed << setprecision(3)
           << " = " << totalkernflops/gigacnt << " Gigaflops" << endl;
      cout << " Total Wall Clock Flops : "
           << scientific << setprecision(3) << totalwallflops;
      cout << fixed << setprecision(3)
           << " = " << totalwallflops/gigacnt << " Gigaflops" << endl;
   }
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
