// The file dbl2_baqr_testers.cpp defines the function with prototypes in
// the file dbl2_baqr_testers.h.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#include "random2_matrices.h"
#include "double_double_functions.h"
#include "dbl2_factorizations.h"
#include "dbl2_factors_testers.h"
#include "dbl2_baqr_host.h"
#include "dbl2_baqr_kernels.h"
#include "dbl2_tabs_host.h"
#include "dbl2_tabs_kernels.h"
#include "dbl2_tabs_testers.h"
#include "dbl2_test_utilities.h"

using namespace std;

void test_real2_blocked_qrbs
 ( int seed, int szt, int nbt, int nrows, int vrb, int mode )
{
   const int sizetile = szt;
   const int numtiles = nbt;
   const int ncols = sizetile*numtiles;
   const int verbose = vrb;

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
   double *solhi = new double[ncols];
   double *sollo = new double[ncols];
   for(int i=0; i<ncols; i++)
   {
      solhi[i] = 1.0;
      sollo[i] = 0.0;
   }
   double *rhshi = new double[nrows];
   double *rhslo = new double[nrows];
   double acchi,acclo;

   for(int i=0; i<nrows; i++)
   {
      rhshi[i] = 0.0;
      rhslo[i] = 0.0;
      for(int j=0; j<ncols; j++) // rhs[i] = rhs[i] + A[i][j]*sol[j];
      {
         ddf_mul(Ahi[i][j],Alo[i][j],solhi[j],sollo[j],&acchi,&acclo);
         ddf_inc(&rhshi[i],&rhslo[i],acchi,acclo);
      }
   }
   if(verbose > 0)
   {
      cout << "The sums of the columns :" << endl;
      for(int i=0; i<nrows; i++)
         cout << "b[" << i << "] : "
              << rhshi[i] << "  " << rhslo[i] << endl;
   }
   double qrtimelapsed_h,bstimelapsed_h;
   bool bvrb = (verbose > 0);
   const double tol = 1.0e-20;
   int fail;

   if((mode == 1) || (mode == 2))
   {
      double *xhi = new double[ncols];
      double *xlo = new double[ncols];
      double *qTrhshi = new double[nrows];
      double *qTrhslo = new double[nrows];

      cout << "-> CPU computes the blocked Householder QR ..." << endl;

      CPU_dbl2_blocked_houseqr
         (nrows,ncols,sizetile,numtiles,Ahi,Alo,Qhi_h,Qlo_h,Rhi_h,Rlo_h,
          &qrtimelapsed_h,bvrb);

      cout << "-> Testing the QR factorization ..." << endl;

      // fail = test_real2_qr_factors
      //   (nrows,ncols,Ahi,Alo,Qhi_h,Qlo_h,Rhi_h,Rlo_h,tol,verbose);
      fail = test_real2_qr_factors_probe
         (nrows,ncols,Ahi,Alo,Qhi_h,Qlo_h,Rhi_h,Rlo_h,tol,2,true);
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
         qTrhshi[i] = 0.0;
         qTrhslo[i] = 0.0;
         for(int j=0; j<nrows; j++) // qTrhs[i] = qTrhs[i] + Q[j][i]*rhs[j];
         {
            ddf_mul(Qhi_h[j][i],Qlo_h[j][i],rhshi[j],rhslo[j],&acchi,&acclo);
            ddf_inc(&qTrhshi[i],&qTrhslo[i],acchi,acclo);
         }
      }
      cout << "-> CPU solves an upper triangular system ..." << endl;

      CPU_dbl2_upper_tiled_solver
         (ncols,sizetile,numtiles,Rhi_h,Rlo_h,qTrhshi,qTrhslo,xhi,xlo,
          &bstimelapsed_h);

      if(verbose > 0)
      {
         cout << "CPU solution computed with tiling :" << endl;
         cout << scientific << setprecision(16);
         for(int i=0; i<ncols; i++)
            cout << "x[" << i << "] : "
                 << xhi[i] << "  " << xlo[i] << endl;
      }
      cout << scientific << setprecision(2);
      cout << "   Sum of CPU errors on solution : "
           << dbl2_Difference_Sum(ncols,solhi,sollo,xhi,xlo) << endl;

      free(xhi); free(xlo); free(qTrhshi); free(qTrhslo);
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
      double *xhi_d = new double[ncols];
      double *xlo_d = new double[ncols];
      double *qTrhshi_d = new double[nrows];
      double *qTrhslo_d = new double[nrows];

      cout << "-> GPU computes the blocked Householder QR ..." << endl;

      GPU_dbl2_blocked_houseqr
         (nrows,ncols,sizetile,numtiles,Ahi,Alo,Qhi_d,Qlo_d,Rhi_d,Rlo_d,
          &houselapsedms,&RTvlapsedms,&tileRlapsedms,&vb2Wlapsedms,
          &WYTlapsedms,&QWYTlapsedms,&Qaddlapsedms,
          &YWTlapsedms,&YWTClapsedms,&Raddlapsedms,&qrtimelapsed_d,
          &qraddcnt,&qrmulcnt,&qrdivcnt,&sqrtcnt,bvrb);

      cout << "-> Testing the QR factorization ..." << endl;

      // fail = test_real2_qr_factors
      //          (nrows,ncols,Ahi,Alo,Qhi_d,Qlo_d,Rhi_d,Rlo_d,tol,verbose);
      fail = test_real2_qr_factors_probe
                (nrows,ncols,Ahi,Alo,Qhi_d,Qlo_d,Rhi_d,Rlo_d,tol,2,true);
      if(fail == 0)
         cout << "The test succeeded." << endl;
      else
      {
         cout << scientific << setprecision(2);
         cout << "The test failed for tol = " << tol << "." << endl;
      }
      // preliminary CPU computation, for testing purposes only
      cout << "-> CPU multiplies Q^T with b ..." << endl;

      for(int i=0; i<nrows; i++)
      {
         qTrhshi_d[i] = 0.0;
         qTrhslo_d[i] = 0.0;
         for(int j=0; j<ncols; j++) // qTrhs[i] = qTrhs[i] + Q[j][i]*rhs[j];
         {
            ddf_mul(Qhi_d[j][i],Qlo_d[j][i],rhshi[j],rhslo[j],&acchi,&acclo);
            ddf_inc(&qTrhshi_d[i],&qTrhslo_d[i],acchi,acclo);
         }
      }
      cout << "-> GPU solves an upper triangular system ..." << endl;

      GPU_dbl2_upper_tiled_solver
         (ncols,sizetile,numtiles,Rhi_d,Rlo_d,qTrhshi_d,qTrhslo_d,xhi_d,xlo_d,
          &invlapsed,&mullapsed,&sublapsed,&bselapsedms,&bstimelapsed_d,
          &bsaddcnt,&bsmulcnt,&bsdivcnt);

      if(verbose > 0)
      {
         cout << scientific << setprecision(16);
         cout << "GPU solution computed with tiling :" << endl;
         for(int i=0; i<ncols; i++)
            cout << "x[" << i << "] : "
                 << xhi_d[i] << "  " << xlo_d[i] << endl;
      }
      cout << scientific << setprecision(2);
      cout << "   Sum of GPU errors on solution : "
           << dbl2_Difference_Sum(ncols,solhi,sollo,xhi_d,xlo_d) << endl;

      free(xhi_d); free(xlo_d); free(qTrhshi_d); free(qTrhslo_d);
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
           << qraddcnt << " x 20" << endl;
      cout << "                    Number of multiplications : "
           << qrmulcnt << " x 23" << endl;
      cout << "                          Number of divisions : "
           << qrdivcnt << " x 70" << endl;
      cout << "                    Number of calls to sqrt() : "
           << sqrtcnt << " x 50" << endl;
      long long int qrflopcnt
          = 20*qraddcnt + 23*qrmulcnt + 70*qrdivcnt + 50*sqrtcnt;
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
           << bsaddcnt << " x 20" << endl;
      cout << "                    Number of multiplications : "
           << bsmulcnt << " x 23" << endl;
      cout << "                          Number of divisions : "
           << bsdivcnt << " x 70" << endl;
      long long int bsflopcnt = 20*bsaddcnt + 23*bsmulcnt + 70*bsdivcnt;
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
      free(Ahi[i]); free(Alo[i]);
      free(Qhi_h[i]); free(Rhi_h[i]);
      free(Qhi_d[i]); free(Rhi_d[i]);
      free(Qlo_h[i]); free(Rlo_h[i]);
      free(Qlo_d[i]); free(Rlo_d[i]);
   }
   free(Ahi); free(Qhi_h); free(Rhi_h); free(Qhi_d); free(Rhi_d);
   free(Alo); free(Qlo_h); free(Rlo_h); free(Qlo_d); free(Rlo_d);
   free(solhi); free(sollo); free(rhshi); free(rhslo);
}

void test_cmplx2_blocked_qrbs
 ( int seed, int szt, int nbt, int nrows, int vrb, int mode )
{
   const int sizetile = szt;
   const int numtiles = nbt;
   const int ncols = sizetile*numtiles;
   const int verbose = vrb;

   cout << "-> Generating a random " << nrows
        << "-by-" << ncols << " matrix ..." << endl;

   double **Arehi = new double*[nrows];
   double **Arelo = new double*[nrows];
   double **Aimhi = new double*[nrows];
   double **Aimlo = new double*[nrows];
   double **Qrehi_h = new double*[nrows];
   double **Qrelo_h = new double*[nrows];
   double **Qimhi_h = new double*[nrows];
   double **Qimlo_h = new double*[nrows];
   double **Qrehi_d = new double*[nrows];
   double **Qrelo_d = new double*[nrows];
   double **Qimhi_d = new double*[nrows];
   double **Qimlo_d = new double*[nrows];
   double **Rrehi_h = new double*[nrows];
   double **Rrelo_h = new double*[nrows];
   double **Rimhi_h = new double*[nrows];
   double **Rimlo_h = new double*[nrows];
   double **Rrehi_d = new double*[nrows];
   double **Rrelo_d = new double*[nrows];
   double **Rimhi_d = new double*[nrows];
   double **Rimlo_d = new double*[nrows];

   for(int i=0; i<nrows; i++)
   {
      Arehi[i] = new double[ncols];
      Arelo[i] = new double[ncols];
      Aimhi[i] = new double[ncols];
      Aimlo[i] = new double[ncols];
      Qrehi_h[i] = new double[nrows];
      Qrelo_h[i] = new double[nrows];
      Qimhi_h[i] = new double[nrows];
      Qimlo_h[i] = new double[nrows];
      Qrehi_d[i] = new double[nrows];
      Qrelo_d[i] = new double[nrows];
      Qimhi_d[i] = new double[nrows];
      Qimlo_d[i] = new double[nrows];
      Rrehi_h[i] = new double[ncols];
      Rrelo_h[i] = new double[ncols];
      Rimhi_h[i] = new double[ncols];
      Rimlo_h[i] = new double[ncols];
      Rrehi_d[i] = new double[ncols];
      Rrelo_d[i] = new double[ncols];
      Rimhi_d[i] = new double[ncols];
      Rimlo_d[i] = new double[ncols];
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
   double *solrehi = new double[ncols];
   double *solrelo = new double[ncols];
   double *solimhi = new double[ncols];
   double *solimlo = new double[ncols];

   for(int i=0; i<ncols; i++)
   {
      solrehi[i] = 1.0; solrelo[i] = 0.0;
      solimhi[i] = 0.0; solimlo[i] = 0.0;
   }
   double *rhsrehi = new double[nrows];
   double *rhsrelo = new double[nrows];
   double *rhsimhi = new double[nrows];
   double *rhsimlo = new double[nrows];
   double acc1hi,acc1lo,acc2hi,acc2lo;
   double acc3hi,acc3lo,acc4hi,acc4lo;

   for(int i=0; i<nrows; i++)
   {
      rhsrehi[i] = 0.0; rhsrelo[i] = 0.0;
      rhsimhi[i] = 0.0; rhsimlo[i] = 0.0;

      for(int j=0; j<ncols; j++) // rhs[i] = rhs[i] + A[i][j]*sol[j];
      {
         ddf_mul(Arehi[i][j],Arelo[i][j],solrehi[j],solrelo[j],
                 &acc1hi,&acc1lo);
         ddf_mul(Aimhi[i][j],Aimlo[i][j],solimhi[j],solimlo[j],
                 &acc2hi,&acc2lo);
         ddf_mul(Aimhi[i][j],Aimlo[i][j],solrehi[j],solrelo[j],
                 &acc3hi,&acc3lo);
         ddf_mul(Arehi[i][j],Arelo[i][j],solimhi[j],solimlo[j],
                 &acc4hi,&acc4lo);
         ddf_dec(&acc1hi,&acc1lo,acc2hi,acc2lo);
         ddf_inc(&rhsrehi[i],&rhsrelo[i],acc1hi,acc1lo);
         ddf_inc(&acc3hi,&acc3lo,acc4hi,acc4lo);
         ddf_inc(&rhsimhi[i],&rhsimlo[i],acc3hi,acc3lo);
      }
   }
   if(verbose > 0)
   {
      cout << "The sums of the columns :" << endl;
      for(int i=0; i<nrows; i++)
      {
         cout << "b[" << i << "]re : "
              << rhsrehi[i] << "  " << rhsrelo[i] << endl;
         cout << "b[" << i << "]im : "
              << rhsimhi[i] << "  " << rhsimlo[i] << endl;
      }
   }
   double qrtimelapsed_h,bstimelapsed_h;
   bool bvrb = (verbose > 0);
   const double tol = 1.0e-20;
   int fail;

   if((mode == 1) || (mode == 2))
   {
      double *xrehi = new double[ncols];
      double *xrelo = new double[ncols];
      double *ximhi = new double[ncols];
      double *ximlo = new double[ncols];
      double *qHrhsrehi = new double[nrows];
      double *qHrhsrelo = new double[nrows];
      double *qHrhsimhi = new double[nrows];
      double *qHrhsimlo = new double[nrows];

      cout << "-> CPU computes the blocked Householder QR ..." << endl;

      CPU_cmplx2_blocked_houseqr
         (nrows,ncols,sizetile,numtiles,
          Arehi,  Arelo,  Aimhi,  Aimlo,
          Qrehi_h,Qrelo_h,Qimhi_h,Qimlo_h,
          Rrehi_h,Rrelo_h,Rimhi_h,Rimlo_h,&qrtimelapsed_h,bvrb);

      cout << "-> Testing the QR factorization ..." << endl;

      // fail = test_cmplx2_qr_factors
      //    (nrows,ncols,Arehi,  Arelo,  Aimhi,  Aimlo,
      //                 Qrehi_h,Qrelo_h,Qimhi_h,Qimlo_h,
      //                 Rrehi_h,Rrelo_h,Rimhi_h,Rimlo_h,tol,verbose);
      fail = test_cmplx2_qr_factors_probe
         (nrows,ncols,Arehi,  Arelo,  Aimhi,  Aimlo,
                      Qrehi_h,Qrelo_h,Qimhi_h,Qimlo_h,
                      Rrehi_h,Rrelo_h,Rimhi_h,Rimlo_h,tol,2,true);
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
         qHrhsrehi[i] = 0.0; qHrhsrelo[i] = 0.0;
         qHrhsimhi[i] = 0.0; qHrhsimlo[i] = 0.0;

         for(int j=0; j<nrows; j++) // qHrhs[i] = qHrhs[i] + Q[j][i]*rhs[j];
         {
            ddf_mul(Qrehi_h[j][i],Qrelo_h[j][i],rhsrehi[j],rhsrelo[j],
                    &acc1hi,&acc1lo);
            ddf_mul(Qimhi_h[j][i],Qimlo_h[j][i],rhsimhi[j],rhsimlo[j],
                    &acc2hi,&acc2lo);
            ddf_mul(Qimhi_h[j][i],Qimlo_h[j][i],rhsrehi[j],rhsrelo[j],
                    &acc3hi,&acc3lo);
            ddf_mul(Qrehi_h[j][i],Qrelo_h[j][i],rhsimhi[j],rhsimlo[j],
                    &acc4hi,&acc4lo);
            ddf_inc(&acc1hi,&acc1lo,acc2hi,acc2lo);
            ddf_inc(&qHrhsrehi[i],&qHrhsrelo[i],acc1hi,acc1lo);
            ddf_dec(&acc4hi,&acc4lo,acc3hi,acc3lo);
            ddf_inc(&qHrhsimhi[i],&qHrhsimlo[i],acc4hi,acc4lo);
         }
      }
      cout << "-> CPU solves an upper triangular system ..." << endl;

      CPU_cmplx2_upper_tiled_solver
         (ncols,sizetile,numtiles,Rrehi_h,Rrelo_h,Rimhi_h,Rimlo_h,
          qHrhsrehi,qHrhsrelo,qHrhsimhi,qHrhsimlo,xrehi,xrelo,ximhi,ximlo,
          &bstimelapsed_h);

      if(verbose > 0)
      {
         cout << "CPU solution computed with tiling :" << endl;
         cout << scientific << setprecision(16);
         for(int i=0; i<ncols; i++)
         {
            cout << "x[" << i << "]re : "
                 << xrehi[i] << "  " << xrelo[i] << endl;
            cout << "x[" << i << "]im : "
                 << ximhi[i] << "  " << ximlo[i] << endl;
         }
      }
      cout << scientific << setprecision(2);
      cout << "   Sum of CPU errors on solution : "
           << cmplx2_Difference_Sum(ncols,solrehi,solrelo,solimhi,solimlo,
                                            xrehi,  xrelo,  ximhi,  ximlo)
           << endl;

      free(xrehi); free(xrelo); free(ximhi); free(ximlo);
      free(qHrhsrehi); free(qHrhsrelo); free(qHrhsimhi); free(qHrhsimlo);
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
      double *xrehi_d = new double[ncols];
      double *xrelo_d = new double[ncols];
      double *ximhi_d = new double[ncols];
      double *ximlo_d = new double[ncols];
      double *qHrhsrehi_d = new double[nrows];
      double *qHrhsrelo_d = new double[nrows];
      double *qHrhsimhi_d = new double[nrows];
      double *qHrhsimlo_d = new double[nrows];

      cout << "-> GPU computes the blocked Householder QR ..." << endl;

      if(verbose > 0) // to verify that A has not changed ...
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
      GPU_cmplx2_blocked_houseqr
         (nrows,ncols,sizetile,numtiles,
          Arehi,  Arelo,  Aimhi,  Aimlo,
          Qrehi_d,Qrelo_d,Qimhi_d,Qimlo_d,
          Rrehi_d,Rrelo_d,Rimhi_d,Rimlo_d,
          &houselapsedms,&RHvlapsedms,&tileRlapsedms,&vb2Wlapsedms,
          &WYTlapsedms,&QWYTlapsedms,&Qaddlapsedms,
          &YWTlapsedms,&YWTClapsedms,&Raddlapsedms,&qrtimelapsed_d,
          &qraddcnt,&qrmulcnt,&qrdivcnt,&sqrtcnt,bvrb);

      cout << "-> Testing the QR factorization ..." << endl;

      // fail = test_cmplx2_qr_factors
      //           (nrows,ncols,Arehi,  Arelo,  Aimhi,  Aimlo,
      //                        Qrehi_d,Qrelo_d,Qimhi_d,Qimlo_d,
      //                        Rrehi_d,Rrelo_d,Rimhi_d,Rimlo_d,tol,verbose);
      fail = test_cmplx2_qr_factors_probe
                (nrows,ncols,Arehi,  Arelo,  Aimhi,  Aimlo,
                             Qrehi_d,Qrelo_d,Qimhi_d,Qimlo_d,
                             Rrehi_d,Rrelo_d,Rimhi_d,Rimlo_d,tol,2,true);
      if(fail == 0)
         cout << "The test succeeded." << endl;
      else
      {
         cout << scientific << setprecision(2);
         cout << "The test failed for tol = " << tol << "." << endl;
      }
      // preliminary CPU computation, for testing purposes
      cout << "-> CPU multiplies Q^H with b ..." << endl;

      for(int i=0; i<nrows; i++)
      {
         qHrhsrehi_d[i] = 0.0; qHrhsrelo_d[i] = 0.0;
         qHrhsimhi_d[i] = 0.0; qHrhsimlo_d[i] = 0.0;

         for(int j=0; j<nrows; j++) // qHrhs[i] = qHrhs[i] + Q[j][i]*rhs[j];
         {
            ddf_mul(Qrehi_d[j][i],Qrelo_d[j][i],rhsrehi[j],rhsrelo[j],
                    &acc1hi,&acc1lo);
            ddf_mul(Qimhi_d[j][i],Qimlo_d[j][i],rhsimhi[j],rhsimlo[j],
                    &acc2hi,&acc2lo);
            ddf_mul(Qimhi_d[j][i],Qimlo_d[j][i],rhsrehi[j],rhsrelo[j],
                    &acc3hi,&acc3lo);
            ddf_mul(Qrehi_d[j][i],Qrelo_d[j][i],rhsimhi[j],rhsimlo[j],
                    &acc4hi,&acc4lo);
            ddf_inc(&acc1hi,&acc1lo,acc2hi,acc2lo);
            ddf_inc(&qHrhsrehi_d[i],&qHrhsrelo_d[i],acc1hi,acc1lo);
            ddf_dec(&acc4hi,&acc4lo,acc3hi,acc3lo);
            ddf_inc(&qHrhsimhi_d[i],&qHrhsimlo_d[i],acc4hi,acc4lo);
         }
      }
      cout << "-> GPU solves an upper triangular system ..." << endl;

      GPU_cmplx2_upper_tiled_solver
         (ncols,sizetile,numtiles,Rrehi_d,Rrelo_d,Rimhi_d,Rimlo_d,
          qHrhsrehi_d,qHrhsrelo_d,qHrhsimhi_d,qHrhsimlo_d,
              xrehi_d,    xrelo_d,    ximhi_d,    ximlo_d,
          &invlapsed,&mullapsed,&sublapsed,&bselapsedms,&bstimelapsed_d,
          &bsaddcnt,&bsmulcnt,&bsdivcnt);

      if(verbose > 0)
      {
         cout << scientific << setprecision(16);
         cout << "GPU solution computed with tiling :" << endl;
         for(int i=0; i<ncols; i++)
         {
            cout << "x[" << i << "]re : "
                 << xrehi_d[i] << "  " << xrelo_d[i] << endl;
            cout << "x[" << i << "]im : "
                 << ximhi_d[i] << "  " << ximlo_d[i] << endl;
         }
      }
      cout << scientific << setprecision(2);
      cout << "   Sum of GPU errors on solution : "
           << cmplx2_Difference_Sum
                 (ncols,solrehi,solrelo,solimhi,solimlo,
                          xrehi_d,xrelo_d,ximhi_d,ximlo_d)
           << endl;
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
           << qraddcnt << " x 20 " << endl;
      cout << "                    Number of multiplications : "
           << qrmulcnt << " x 23 " << endl;
      cout << "                          Number of divisions : "
           << qrdivcnt << " x 70 " << endl;
      cout << "                    Number of calls to sqrt() : "
           << sqrtcnt << " x 50 " << endl;
      long long int qrflopcnt
          = 20*qraddcnt + 23*qrmulcnt + 70*qrdivcnt + 50*sqrtcnt;
      cout << "    Total number of floating-point operations : "
           << qrflopcnt << endl;
      cout << endl;
      double qrkernflops = 1000.0*((double) qrflopcnt)/totlapsedms;
      double qrwallflops = ((double) qrflopcnt)/qrtimelapsed_d;
      const int gigacnt = pow(2.0,30);
      cout << "Kernel Time Flops : "
           << scientific << setprecision(3) << qrkernflops;
      cout << fixed << setprecision(3)
           << " = " << qrkernflops/gigacnt << " Gigaflops" << endl;
      cout << " Wall Clock Flops : "
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
           << bsaddcnt << " x 20 " << endl;
      cout << "                    Number of multiplications : "
           << bsmulcnt << " x 23 " << endl;
      cout << "                          Number of divisions : "
           << bsdivcnt << " x 70 " << endl;
      long long int bsflopcnt = 20*bsaddcnt + 23*bsmulcnt + 70*bsdivcnt;
      cout << "    Total number of floating-point operations : "
           << bsflopcnt << endl;
      cout << endl;
      double bskernflops = 1000.0*((double) bsflopcnt)/bselapsedms;
      double bswallflops = ((double) bsflopcnt)/bstimelapsed_d;
      // const int gigacnt = pow(2.0,30);
      cout << "Kernel Time Flops : "
           << scientific << setprecision(3) << bskernflops;
      cout << fixed << setprecision(3)
           << " = " << bskernflops/gigacnt << " Gigaflops" << endl;
      cout << " Wall Clock Flops : "
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
      free(Arehi[i]);
      free(Qrehi_h[i]); free(Rrehi_h[i]);
      free(Qrehi_d[i]); free(Rrehi_d[i]);
      free(Arelo[i]);
      free(Qrelo_h[i]); free(Rrelo_h[i]);
      free(Qrelo_d[i]); free(Rrelo_d[i]);
      free(Aimhi[i]);
      free(Qimhi_h[i]); free(Rimhi_h[i]);
      free(Qimhi_d[i]); free(Rimhi_d[i]);
      free(Aimlo[i]);
      free(Qimlo_h[i]); free(Rimlo_h[i]);
      free(Qimlo_d[i]); free(Rimlo_d[i]);
   }
   free(Arehi); free(Qrehi_h); free(Rrehi_h); free(Qrehi_d); free(Rrehi_d);
   free(Arelo); free(Qrelo_h); free(Rrelo_h); free(Qrelo_d); free(Rrelo_d);
   free(Aimhi); free(Qimhi_h); free(Rimhi_h); free(Qimhi_d); free(Rimhi_d);
   free(Aimlo); free(Qimlo_h); free(Rimlo_h); free(Qimlo_d); free(Rimlo_d);
   free(solrehi); free(solrelo); free(solimhi); free(solimlo);
   free(rhsrehi); free(rhsrelo); free(rhsimhi); free(rhsimlo);
}
