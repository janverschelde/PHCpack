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
   double timelapsed_d;
   double houselapsedms,RTvlapsedms,tileRlapsedms,vb2Wlapsedms;
   double WYTlapsedms,QWYTlapsedms,Qaddlapsedms;
   double YWTlapsedms,YWTClapsedms,Raddlapsedms;
   long long int addcnt = 0;
   long long int mulcnt = 0;
   long long int divcnt = 0;
   long long int sqrtcnt = 0;

   if((mode == 0) || (mode == 2))
   {
      cout << "-> GPU computes the blocked Householder QR ..." << endl;

      GPU_dbl2_blocked_houseqr
         (nrows,ncols,sizetile,numtiles,Ahi,Alo,Qhi_d,Qlo_d,Rhi_d,Rlo_d,
          &houselapsedms,&RTvlapsedms,&tileRlapsedms,&vb2Wlapsedms,
          &WYTlapsedms,&QWYTlapsedms,&Qaddlapsedms,
          &YWTlapsedms,&YWTClapsedms,&Raddlapsedms,&timelapsed_d,
          &addcnt,&mulcnt,&divcnt,&sqrtcnt,bvrb);

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
      cout << "        Total GPU wall clock computation time : ";
      cout << fixed << setprecision(3) << timelapsed_d << " seconds." << endl;
      cout << endl;
      cout << "             Number of additions/subtractions : "
           << addcnt << " x 20 " << endl;
      cout << "                    Number of multiplications : "
           << mulcnt << " x 23 " << endl;
      cout << "                          Number of divisions : "
           << divcnt << " x 70 " << endl;
      cout << "                    Number of calls to sqrt() : "
           << sqrtcnt << " x 50 " << endl;
      long long int flopcnt = 20*addcnt + 23*mulcnt + 70*divcnt + 50*sqrtcnt;
      cout << "    Total number of floating-point operations : "
           << flopcnt << endl;
      cout << endl;
      double kernflops = 1000.0*((double) flopcnt)/totlapsedms;
      double wallflops = ((double) flopcnt)/timelapsed_d;
      const int gigacnt = pow(2.0,30);
      cout << "Kernel Time Flops : "
           << scientific << setprecision(3) << kernflops;
      cout << fixed << setprecision(3)
           << " = " << kernflops/gigacnt << " Gigaflops" << endl;
      cout << " Wall Clock Flops : "
           << scientific << setprecision(3) << wallflops;
      cout << fixed << setprecision(3)
           << " = " << wallflops/gigacnt << " Gigaflops" << endl;
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
   double timelapsed_d;
   double houselapsedms,RHvlapsedms,tileRlapsedms,vb2Wlapsedms;
   double WYTlapsedms,QWYTlapsedms,Qaddlapsedms;
   double YWTlapsedms,YWTClapsedms,Raddlapsedms;
   long long int addcnt = 0;
   long long int mulcnt = 0;
   long long int divcnt = 0;
   long long int sqrtcnt = 0;

   if((mode == 0) || (mode == 2))
   {
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
          &YWTlapsedms,&YWTClapsedms,&Raddlapsedms,&timelapsed_d,
          &addcnt,&mulcnt,&divcnt,&sqrtcnt,bvrb);

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
      cout << "        Total GPU wall clock computation time : ";
      cout << fixed << setprecision(3) << timelapsed_d << " seconds." << endl;
      cout << endl;
      cout << "             Number of additions/subtractions : "
           << addcnt << " x 20 " << endl;
      cout << "                    Number of multiplications : "
           << mulcnt << " x 23 " << endl;
      cout << "                          Number of divisions : "
           << divcnt << " x 70 " << endl;
      cout << "                    Number of calls to sqrt() : "
           << sqrtcnt << " x 50 " << endl;
      long long int flopcnt = 20*addcnt + 23*mulcnt + 70*divcnt + 50*sqrtcnt;
      cout << "    Total number of floating-point operations : "
           << flopcnt << endl;
      cout << endl;
      double kernflops = 1000.0*((double) flopcnt)/totlapsedms;
      double wallflops = ((double) flopcnt)/timelapsed_d;
      const int gigacnt = pow(2.0,30);
      cout << "Kernel Time Flops : "
           << scientific << setprecision(3) << kernflops;
      cout << fixed << setprecision(3)
           << " = " << kernflops/gigacnt << " Gigaflops" << endl;
      cout << " Wall Clock Flops : "
           << scientific << setprecision(3) << wallflops;
      cout << fixed << setprecision(3)
           << " = " << wallflops/gigacnt << " Gigaflops" << endl;
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
