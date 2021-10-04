// The file dbl4_baqr_testers.cpp defines the function with prototypes in
// the file dbl4_baqr_testers.h.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#include "random4_matrices.h"
#include "quad_double_functions.h"
#include "dbl4_factorizations.h"
#include "dbl4_factors_testers.h"
#include "dbl4_baqr_host.h"
#include "dbl4_baqr_kernels.h"
#include "dbl4_tabs_host.h"
#include "dbl4_tabs_kernels.h"
#include "dbl4_tabs_testers.h"
#include "dbl4_test_utilities.h"

using namespace std;

void test_real4_blocked_qrbs
 ( int seed, int szt, int nbt, int nrows, int vrb, int mode )
{
   const int sizetile = szt;
   const int numtiles = nbt;
   const int ncols = sizetile*numtiles;
   const int verbose = vrb;

   cout << "-> Generating a random " << nrows
        << "-by-" << ncols << " matrix ..." << endl;

   double **Ahihi = new double*[nrows];
   double **Alohi = new double*[nrows];
   double **Ahilo = new double*[nrows];
   double **Alolo = new double*[nrows];
   double **Qhihi_h = new double*[nrows];
   double **Qlohi_h = new double*[nrows];
   double **Qhilo_h = new double*[nrows];
   double **Qlolo_h = new double*[nrows];
   double **Qhihi_d = new double*[nrows];
   double **Qlohi_d = new double*[nrows];
   double **Qhilo_d = new double*[nrows];
   double **Qlolo_d = new double*[nrows];
   double **Rhihi_h = new double*[nrows];
   double **Rlohi_h = new double*[nrows];
   double **Rhilo_h = new double*[nrows];
   double **Rlolo_h = new double*[nrows];
   double **Rhihi_d = new double*[nrows];
   double **Rlohi_d = new double*[nrows];
   double **Rhilo_d = new double*[nrows];
   double **Rlolo_d = new double*[nrows];

   for(int i=0; i<nrows; i++)
   {
      Ahihi[i] = new double[ncols];
      Alohi[i] = new double[ncols];
      Ahilo[i] = new double[ncols];
      Alolo[i] = new double[ncols];
      Qhihi_h[i] = new double[nrows];
      Qlohi_h[i] = new double[nrows];
      Qhilo_h[i] = new double[nrows];
      Qlolo_h[i] = new double[nrows];
      Qhihi_d[i] = new double[nrows];
      Qlohi_d[i] = new double[nrows];
      Qhilo_d[i] = new double[nrows];
      Qlolo_d[i] = new double[nrows];
      Rhihi_h[i] = new double[ncols];
      Rlohi_h[i] = new double[ncols];
      Rhilo_h[i] = new double[ncols];
      Rlolo_h[i] = new double[ncols];
      Rhihi_d[i] = new double[ncols];
      Rlohi_d[i] = new double[ncols];
      Rhilo_d[i] = new double[ncols];
      Rlolo_d[i] = new double[ncols];
   }
   random_dbl4_matrix(nrows,ncols,Ahihi,Alohi,Ahilo,Alolo);

   if(verbose > 0)
   {
      cout << scientific << setprecision(16);

      cout << "A random matrix :" << endl;
      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++)
            cout << "A[" << i << "][" << j << "] : "
                 << Ahihi[i][j] << "  " << Alohi[i][j] << endl
                 << Ahilo[i][j] << "  " << Alolo[i][j] << endl;
   }
   double *solhihi = new double[ncols];
   double *sollohi = new double[ncols];
   double *solhilo = new double[ncols];
   double *sollolo = new double[ncols];

   for(int i=0; i<ncols; i++)
   {
      solhihi[i] = 1.0;
      sollohi[i] = 0.0;
      solhilo[i] = 0.0;
      sollolo[i] = 0.0;
   }
   double *rhshihi = new double[nrows];
   double *rhslohi = new double[nrows];
   double *rhshilo = new double[nrows];
   double *rhslolo = new double[nrows];
   double acchihi,acclohi,acchilo,acclolo;

   for(int i=0; i<nrows; i++)
   {
      rhshihi[i] = 0.0;
      rhslohi[i] = 0.0;
      rhshilo[i] = 0.0;
      rhslolo[i] = 0.0;

      for(int j=0; j<ncols; j++) // rhs[i] = rhs[i] + A[i][j]*sol[j];
      {
         qdf_mul(   Ahihi[i][j],Alohi[i][j],Ahilo[i][j],Alolo[i][j],
                  solhihi[j], sollohi[j], solhilo[j], sollolo[j],
                 &acchihi,   &acclohi,   &acchilo,   &acclolo);
         qdf_inc(&rhshihi[i],&rhslohi[i],&rhshilo[i],&rhslolo[i],
                  acchihi,    acclohi,    acchilo,    acclolo);
      }
   }
   if(verbose > 0)
   {
      cout << "The sums of the columns :" << endl;
      for(int i=0; i<nrows; i++)
         cout << "b[" << i << "] : "
              << rhshihi[i] << "  " << rhslohi[i] << endl
              << rhshilo[i] << "  " << rhslolo[i] << endl;
   }
   double qrtimelapsed_h,bstimelapsed_h;
   bool bvrb = (verbose > 0);
   const double tol = 1.0e-20;
   int fail;

   if((mode == 1) || (mode == 2))
   {
      double *xhihi = new double[ncols];
      double *xlohi = new double[ncols];
      double *xhilo = new double[ncols];
      double *xlolo = new double[ncols];
      double *qTrhshihi = new double[nrows];
      double *qTrhslohi = new double[nrows];
      double *qTrhshilo = new double[nrows];
      double *qTrhslolo = new double[nrows];

      cout << "-> CPU computes the blocked Householder QR ..." << endl;

      CPU_dbl4_blocked_houseqr
         (nrows,ncols,sizetile,numtiles,
          Ahihi,  Alohi,  Ahilo,  Alolo,
          Qhihi_h,Qlohi_h,Qhilo_h,Qlolo_h,
          Rhihi_h,Rlohi_h,Rhilo_h,Rlolo_h,&qrtimelapsed_h,bvrb);

      cout << "-> Testing the QR factorization ..." << endl;

      // fail = test_real2_qr_factors
      //   (nrows,ncols,Ahi,Alo,Qhi_h,Qlo_h,Rhi_h,Rlo_h,tol,verbose);
      fail = test_real4_qr_factors_probe
         (nrows,ncols,Ahihi,  Alohi,  Ahilo,  Alolo,
                      Qhihi_h,Qlohi_h,Qhilo_h,Qlolo_h,
                      Rhihi_h,Rlohi_h,Rhilo_h,Rlolo_h,tol,2,true);
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
         qTrhshihi[i] = 0.0;
         qTrhslohi[i] = 0.0;
         qTrhshilo[i] = 0.0;
         qTrhslolo[i] = 0.0;

         for(int j=0; j<nrows; j++) // qTrhs[i] = qTrhs[i] + Q[j][i]*rhs[j];
         {
            qdf_mul(Qhihi_h[j][i],Qlohi_h[j][i],Qhilo_h[j][i],Qlolo_h[j][i],
                  rhshihi[j],   rhslohi[j],   rhshilo[j],   rhslolo[j],
                 &acchihi,     &acclohi,     &acchilo,     &acclolo);
            qdf_inc(&qTrhshihi[i],&qTrhslohi[i],&qTrhshilo[i],&qTrhslolo[i],
                       acchihi,      acclohi,      acchilo,      acclolo);
         }
      }
      cout << "-> CPU solves an upper triangular system ..." << endl;

      CPU_dbl4_upper_tiled_solver
         (ncols,sizetile,numtiles,
              Rhihi_h,  Rlohi_h,  Rhilo_h,  Rlolo_h,
          qTrhshihi,qTrhslohi,qTrhshilo,qTrhslolo,
              xhihi,    xlohi,    xhilo,    xlolo,&bstimelapsed_h);

      if(verbose > 0)
      {
         cout << "CPU solution computed with tiling :" << endl;
         cout << scientific << setprecision(16);
         for(int i=0; i<ncols; i++)
            cout << "x[" << i << "] : "
                 << xhihi[i] << "  " << xlohi[i] << endl
                 << xhilo[i] << "  " << xlolo[i] << endl;
      }
      cout << scientific << setprecision(2);
      cout << "   Sum of CPU errors on solution : "
           << dbl4_Difference_Sum
                 (ncols,solhihi,sollohi,solhilo,sollolo,
                          xhihi,  xlohi,  xhilo,  xlolo) << endl;

      free(xhihi); free(xlohi); free(xhilo); free(xlolo);
      free(qTrhshihi); free(qTrhslohi); free(qTrhshilo); free(qTrhslolo);
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
   double addover = 0.0;
   double mulover = 0.0;
   double divover = 0.0;

   if((mode == 0) || (mode == 2))
   {
      double *xhihi_d = new double[ncols];
      double *xlohi_d = new double[ncols];
      double *xhilo_d = new double[ncols];
      double *xlolo_d = new double[ncols];
      double *qTrhshihi_d = new double[nrows];
      double *qTrhslohi_d = new double[nrows];
      double *qTrhshilo_d = new double[nrows];
      double *qTrhslolo_d = new double[nrows];

      cout << "-> GPU computes the blocked Householder QR ..." << endl;

      GPU_dbl4_blocked_houseqr
         (nrows,ncols,sizetile,numtiles,
          Ahihi,  Alohi,  Ahilo,  Alolo,
          Qhihi_d,Qlohi_d,Qhilo_d,Qlolo_d,
          Rhihi_d,Rlohi_d,Rhilo_d,Rlolo_d,
          &houselapsedms,&RTvlapsedms,&tileRlapsedms,&vb2Wlapsedms,
          &WYTlapsedms,&QWYTlapsedms,&Qaddlapsedms,
          &YWTlapsedms,&YWTClapsedms,&Raddlapsedms,&qrtimelapsed_d,
          &qraddcnt,&qrmulcnt,&qrdivcnt,&sqrtcnt,bvrb);

      cout << "-> Testing the QR factorization ..." << endl;

      // fail = test_real2_qr_factors
      //          (nrows,ncols,Ahi,Alo,Qhi_d,Qlo_d,Rhi_d,Rlo_d,tol,verbose);
      fail = test_real4_qr_factors_probe
                (nrows,ncols,Ahihi,  Alohi,  Ahilo,  Alolo,
                             Qhihi_d,Qlohi_d,Qhilo_d,Qlolo_d,
                             Rhihi_d,Rlohi_d,Rhilo_d,Rlolo_d,tol,2,true);
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
         qTrhshihi_d[i] = 0.0;
         qTrhslohi_d[i] = 0.0;
         qTrhshilo_d[i] = 0.0;
         qTrhslolo_d[i] = 0.0;

         for(int j=0; j<ncols; j++) // qTrhs[i] = qTrhs[i] + Q[j][i]*rhs[j];
         {
            qdf_mul(Qhihi_d[j][i],Qlohi_d[j][i],Qhilo_d[j][i],Qlolo_d[j][i],
                  rhshihi[j],   rhslohi[j],   rhshilo[j],   rhslolo[j],
                 &acchihi,     &acclohi,     &acchilo,     &acclolo);
            qdf_inc(&qTrhshihi_d[i],&qTrhslohi_d[i],
                    &qTrhshilo_d[i],&qTrhslolo_d[i],
                    acchihi,acclohi,acchilo,acclolo);
         }
      }
      cout << "-> GPU solves an upper triangular system ..." << endl;

      GPU_dbl4_upper_tiled_solver
         (ncols,sizetile,numtiles,
              Rhihi_d,    Rlohi_d,    Rhilo_d,    Rlolo_d,
          qTrhshihi_d,qTrhslohi_d,qTrhshilo_d,qTrhslolo_d,
              xhihi_d,    xlohi_d,    xhilo_d,    xlolo_d,
          &invlapsed,&mullapsed,&sublapsed,&bselapsedms,&bstimelapsed_d,
          &bsaddcnt,&addover,&bsmulcnt,&mulover,&bsdivcnt,&divover);

      if(verbose > 0)
      {
         cout << scientific << setprecision(16);
         cout << "GPU solution computed with tiling :" << endl;
         for(int i=0; i<ncols; i++)
            cout << "x[" << i << "] : "
                 << xhihi_d[i] << "  " << xlohi_d[i] << endl
                 << xhilo_d[i] << "  " << xlolo_d[i] << endl;
      }
      cout << scientific << setprecision(2);
      cout << "   Sum of GPU errors on solution : "
           << dbl4_Difference_Sum
                 (ncols,solhihi,sollohi,solhilo,sollolo,
                          xhihi_d,xlohi_d,xhilo_d,xlolo_d) << endl;

      free(xhihi_d); free(xlohi_d); free(xhilo_d); free(xlolo_d);
      free(qTrhshihi_d); free(qTrhslohi_d);
      free(qTrhshilo_d); free(qTrhslolo_d);
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
           << qraddcnt << " x 89" << endl;
      cout << "                    Number of multiplications : "
           << qrmulcnt << " x 336" << endl;
      cout << "                          Number of divisions : "
           << qrdivcnt << " x 893" << endl;
      cout << "                    Number of calls to sqrt() : "
           << sqrtcnt << " x 1345" << endl;
      long long int qrflopcnt
          = 89*qraddcnt + 336*qrmulcnt + 893*qrdivcnt + 1345*sqrtcnt;
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
           << bsaddcnt << " x 89" << endl;
      cout << "                    Number of multiplications : "
           << bsmulcnt << " x 336" << endl;
      cout << "                          Number of divisions : "
           << bsdivcnt << " x 893" << endl;
      long long int bsflopcnt = 89*bsaddcnt + 336*bsmulcnt + 893*bsdivcnt;
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
      free(Ahihi[i]); free(Alohi[i]); free(Ahilo[i]); free(Alolo[i]);
      free(Qhihi_h[i]); free(Rhihi_h[i]);
      free(Qhilo_d[i]); free(Rhilo_d[i]);
      free(Qlohi_h[i]); free(Rlohi_h[i]);
      free(Qlolo_d[i]); free(Rlolo_d[i]);
   }
   free(Ahihi); free(Alohi); free(Ahilo); free(Alolo);
   free(Qhihi_h); free(Qlohi_h); free(Qhilo_h); free(Qlolo_h);
   free(Rhihi_h); free(Rlohi_h); free(Rhilo_h); free(Rlolo_h);
   free(Qhihi_d); free(Qlohi_d); free(Qhilo_d); free(Qlolo_d);
   free(Rhihi_d); free(Rlohi_d); free(Rhilo_d); free(Rlolo_d);
   free(solhihi); free(sollohi); free(solhilo); free(sollolo);
   free(rhshihi); free(rhslohi); free(rhshilo); free(rhslolo);
}

void test_cmplx4_blocked_qrbs
 ( int seed, int szt, int nbt, int nrows, int vrb, int mode )
{
   const int sizetile = szt;
   const int numtiles = nbt;
   const int ncols = sizetile*numtiles;
   const int verbose = vrb;

   cout << "-> Generating a random " << nrows
        << "-by-" << ncols << " matrix ..." << endl;

   double **Arehihi = new double*[nrows];
   double **Arelohi = new double*[nrows];
   double **Arehilo = new double*[nrows];
   double **Arelolo = new double*[nrows];
   double **Aimhihi = new double*[nrows];
   double **Aimlohi = new double*[nrows];
   double **Aimhilo = new double*[nrows];
   double **Aimlolo = new double*[nrows];
   double **Qrehihi_h = new double*[nrows];
   double **Qrelohi_h = new double*[nrows];
   double **Qrehilo_h = new double*[nrows];
   double **Qrelolo_h = new double*[nrows];
   double **Qimhihi_h = new double*[nrows];
   double **Qimlohi_h = new double*[nrows];
   double **Qimhilo_h = new double*[nrows];
   double **Qimlolo_h = new double*[nrows];
   double **Qrehihi_d = new double*[nrows];
   double **Qrelohi_d = new double*[nrows];
   double **Qrehilo_d = new double*[nrows];
   double **Qrelolo_d = new double*[nrows];
   double **Qimhihi_d = new double*[nrows];
   double **Qimlohi_d = new double*[nrows];
   double **Qimhilo_d = new double*[nrows];
   double **Qimlolo_d = new double*[nrows];
   double **Rrehihi_h = new double*[nrows];
   double **Rrelohi_h = new double*[nrows];
   double **Rrehilo_h = new double*[nrows];
   double **Rrelolo_h = new double*[nrows];
   double **Rimhihi_h = new double*[nrows];
   double **Rimlohi_h = new double*[nrows];
   double **Rimhilo_h = new double*[nrows];
   double **Rimlolo_h = new double*[nrows];
   double **Rrehihi_d = new double*[nrows];
   double **Rrelohi_d = new double*[nrows];
   double **Rrehilo_d = new double*[nrows];
   double **Rrelolo_d = new double*[nrows];
   double **Rimhihi_d = new double*[nrows];
   double **Rimlohi_d = new double*[nrows];
   double **Rimhilo_d = new double*[nrows];
   double **Rimlolo_d = new double*[nrows];

   for(int i=0; i<nrows; i++)
   {
      Arehihi[i] = new double[ncols];
      Arelohi[i] = new double[ncols];
      Arehilo[i] = new double[ncols];
      Arelolo[i] = new double[ncols];
      Aimhihi[i] = new double[ncols];
      Aimlohi[i] = new double[ncols];
      Aimhilo[i] = new double[ncols];
      Aimlolo[i] = new double[ncols];
      Qrehihi_h[i] = new double[nrows];
      Qrelohi_h[i] = new double[nrows];
      Qrehilo_h[i] = new double[nrows];
      Qrelolo_h[i] = new double[nrows];
      Qimhihi_h[i] = new double[nrows];
      Qimlohi_h[i] = new double[nrows];
      Qimhilo_h[i] = new double[nrows];
      Qimlolo_h[i] = new double[nrows];
      Qrehihi_d[i] = new double[nrows];
      Qrelohi_d[i] = new double[nrows];
      Qrehilo_d[i] = new double[nrows];
      Qrelolo_d[i] = new double[nrows];
      Qimhihi_d[i] = new double[nrows];
      Qimlohi_d[i] = new double[nrows];
      Qimhilo_d[i] = new double[nrows];
      Qimlolo_d[i] = new double[nrows];
      Rrehihi_h[i] = new double[ncols];
      Rrelohi_h[i] = new double[ncols];
      Rrehilo_h[i] = new double[ncols];
      Rrelolo_h[i] = new double[ncols];
      Rimhihi_h[i] = new double[ncols];
      Rimlohi_h[i] = new double[ncols];
      Rimhilo_h[i] = new double[ncols];
      Rimlolo_h[i] = new double[ncols];
      Rrehihi_d[i] = new double[ncols];
      Rrelohi_d[i] = new double[ncols];
      Rrehilo_d[i] = new double[ncols];
      Rrelolo_d[i] = new double[ncols];
      Rimhihi_d[i] = new double[ncols];
      Rimlohi_d[i] = new double[ncols];
      Rimhilo_d[i] = new double[ncols];
      Rimlolo_d[i] = new double[ncols];
   }
   random_cmplx4_matrix
      (nrows,ncols,Arehihi,Arelohi,Arehilo,Arelolo,
                   Aimhihi,Aimlohi,Aimhilo,Aimlolo);

   if(verbose > 0)
   {
      cout << scientific << setprecision(16);

      cout << "A random matrix :" << endl;
      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++)
         {
            cout << "A[" << i << "][" << j << "]re : "
                 << Arehihi[i][j] << "  " << Arelohi[i][j] << endl
                 << Arehilo[i][j] << "  " << Arelolo[i][j] << endl;
            cout << "A[" << i << "][" << j << "]im : "
                 << Aimhihi[i][j] << "  " << Aimlohi[i][j] << endl
                 << Aimhilo[i][j] << "  " << Aimlolo[i][j] << endl;
         }
   }
   double *solrehihi = new double[ncols];
   double *solrelohi = new double[ncols];
   double *solrehilo = new double[ncols];
   double *solrelolo = new double[ncols];
   double *solimhihi = new double[ncols];
   double *solimlohi = new double[ncols];
   double *solimhilo = new double[ncols];
   double *solimlolo = new double[ncols];

   for(int i=0; i<ncols; i++)
   {
      solrehihi[i] = 1.0; solrelohi[i] = 0.0;
      solrehilo[i] = 0.0; solrelolo[i] = 0.0;
      solimhihi[i] = 0.0; solimlohi[i] = 0.0;
      solimhilo[i] = 0.0; solimlolo[i] = 0.0;
   }
   double *rhsrehihi = new double[nrows];
   double *rhsrelohi = new double[nrows];
   double *rhsrehilo = new double[nrows];
   double *rhsrelolo = new double[nrows];
   double *rhsimhihi = new double[nrows];
   double *rhsimlohi = new double[nrows];
   double *rhsimhilo = new double[nrows];
   double *rhsimlolo = new double[nrows];
   double acc1hihi,acc1lohi,acc1hilo,acc1lolo;
   double acc2hihi,acc2lohi,acc2hilo,acc2lolo;
   double acc3hihi,acc3lohi,acc3hilo,acc3lolo;
   double acc4hihi,acc4lohi,acc4hilo,acc4lolo;

   for(int i=0; i<nrows; i++)
   {
      rhsrehihi[i] = 0.0; rhsrelohi[i] = 0.0;
      rhsrehilo[i] = 0.0; rhsrelolo[i] = 0.0;
      rhsimhihi[i] = 0.0; rhsimlohi[i] = 0.0;
      rhsimhilo[i] = 0.0; rhsimlolo[i] = 0.0;

      for(int j=0; j<ncols; j++) // rhs[i] = rhs[i] + A[i][j]*sol[j];
      {
         qdf_mul(Arehihi[i][j],Arelohi[i][j],Arehilo[i][j],Arelolo[i][j],
               solrehihi[j], solrelohi[j], solrehilo[j], solrelolo[j],
               &acc1hihi,    &acc1lohi,    &acc1hilo,    &acc1lolo);
         qdf_mul(Aimhihi[i][j],Aimlohi[i][j],Aimhilo[i][j],Aimlolo[i][j],
               solimhihi[j], solimlohi[j], solimhilo[j], solimlolo[j],
               &acc2hihi,    &acc2lohi,    &acc2hilo,    &acc2lolo);
         qdf_mul(Aimhihi[i][j],Aimlohi[i][j],Aimhilo[i][j],Aimlolo[i][j],
               solrehihi[j], solrelohi[j], solrehilo[j], solrelolo[j],
               &acc3hihi,    &acc3lohi,    &acc3hilo,    &acc3lolo);
         qdf_mul(Arehihi[i][j],Arelohi[i][j],Arehilo[i][j],Arelolo[i][j],
               solimhihi[j], solimlohi[j], solimhilo[j], solimlolo[j],
               &acc4hihi,    &acc4lohi,    &acc4hilo,    &acc4lolo);
         qdf_dec(&acc1hihi,&acc1lohi,&acc1hilo,&acc1lolo,
                  acc2hihi, acc2lohi, acc2hilo, acc2lolo);
         qdf_inc(&rhsrehihi[i],&rhsrelohi[i],&rhsrehilo[i],&rhsrelolo[i],
                   acc1hihi,     acc1lohi,     acc1hilo,     acc1lolo);
         qdf_inc(&acc3hihi,&acc3lohi,&acc3hilo,&acc3lolo,
                  acc4hihi, acc4lohi, acc4hilo, acc4lolo);
         qdf_inc(&rhsimhihi[i],&rhsimlohi[i],&rhsimhilo[i],&rhsimlolo[i],
                   acc3hihi,     acc3lohi,     acc3hilo,     acc3lolo);
      }
   }
   if(verbose > 0)
   {
      cout << "The sums of the columns :" << endl;
      for(int i=0; i<nrows; i++)
      {
         cout << "b[" << i << "]re : "
              << rhsrehihi[i] << "  " << rhsrelohi[i] << endl
              << rhsrehilo[i] << "  " << rhsrelolo[i] << endl;
         cout << "b[" << i << "]im : "
              << rhsimhihi[i] << "  " << rhsimlohi[i] << endl
              << rhsimhilo[i] << "  " << rhsimlolo[i] << endl;
      }
   }
   double qrtimelapsed_h,bstimelapsed_h;
   bool bvrb = (verbose > 0);
   const double tol = 1.0e-20;
   int fail;

   if((mode == 1) || (mode == 2))
   {
      double *xrehihi = new double[ncols];
      double *xrelohi = new double[ncols];
      double *xrehilo = new double[ncols];
      double *xrelolo = new double[ncols];
      double *ximhihi = new double[ncols];
      double *ximlohi = new double[ncols];
      double *ximhilo = new double[ncols];
      double *ximlolo = new double[ncols];
      double *qHrhsrehihi = new double[nrows];
      double *qHrhsrelohi = new double[nrows];
      double *qHrhsrehilo = new double[nrows];
      double *qHrhsrelolo = new double[nrows];
      double *qHrhsimhihi = new double[nrows];
      double *qHrhsimlohi = new double[nrows];
      double *qHrhsimhilo = new double[nrows];
      double *qHrhsimlolo = new double[nrows];

      cout << "-> CPU computes the blocked Householder QR ..." << endl;

      CPU_cmplx4_blocked_houseqr
         (nrows,ncols,sizetile,numtiles,
          Arehihi,  Arelohi,  Arehilo,  Arelolo,
          Aimhihi,  Aimlohi,  Aimhilo,  Aimlolo,
          Qrehihi_h,Qrelohi_h,Qrehilo_h,Qrelolo_h,
          Qimhihi_h,Qimlohi_h,Qimhilo_h,Qimlolo_h,
          Rrehihi_h,Rrelohi_h,Rrehilo_h,Rrelolo_h,
          Rimhihi_h,Rimlohi_h,Rimhilo_h,Rimlolo_h,&qrtimelapsed_h,bvrb);

      cout << "-> Testing the QR factorization ..." << endl;

      // fail = test_cmplx2_qr_factors
      //    (nrows,ncols,Arehi,  Arelo,  Aimhi,  Aimlo,
      //                 Qrehi_h,Qrelo_h,Qimhi_h,Qimlo_h,
      //                 Rrehi_h,Rrelo_h,Rimhi_h,Rimlo_h,tol,verbose);
      fail = test_cmplx4_qr_factors_probe
         (nrows,ncols,Arehihi,  Arelohi,  Arehilo,  Arelolo,
                      Aimhihi,  Aimlohi,  Aimhilo,  Aimlolo,
                      Qrehihi_h,Qrelohi_h,Qrehilo_h,Qrelolo_h,
                      Qimhihi_h,Qimlohi_h,Qimhilo_h,Qimlolo_h,
                      Rrehihi_h,Rrelohi_h,Rrehilo_h,Rrelolo_h,
                      Rimhihi_h,Rimlohi_h,Rimhilo_h,Rimlolo_h,tol,2,true);
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
         qHrhsrehihi[i] = 0.0; qHrhsrelohi[i] = 0.0;
         qHrhsrehilo[i] = 0.0; qHrhsrelolo[i] = 0.0;
         qHrhsimhihi[i] = 0.0; qHrhsimlohi[i] = 0.0;
         qHrhsimhilo[i] = 0.0; qHrhsimlolo[i] = 0.0;

         for(int j=0; j<nrows; j++) // qHrhs[i] = qHrhs[i] + Q[j][i]*rhs[j];
         {
            qdf_mul(Qrehihi_h[j][i],Qrelohi_h[j][i],
                    Qrehilo_h[j][i],Qrelolo_h[j][i],
                    rhsrehihi[j],rhsrelohi[j],rhsrehilo[j],rhsrelolo[j],
                    &acc1hihi,&acc1lohi,&acc1hilo,&acc1lolo);
            qdf_mul(Qimhihi_h[j][i],Qimlohi_h[j][i],
                    Qimhilo_h[j][i],Qimlolo_h[j][i],
                    rhsimhihi[j],rhsimlohi[j],rhsimhilo[j],rhsimlolo[j],
                    &acc2hihi,&acc2lohi,&acc2hilo,&acc2lolo);
            qdf_mul(Qimhihi_h[j][i],Qimlohi_h[j][i],
                    Qimhilo_h[j][i],Qimlolo_h[j][i],
                    rhsrehihi[j],rhsrelohi[j],rhsrehilo[j],rhsrelolo[j],
                    &acc3hihi,&acc3lohi,&acc3hilo,&acc3lolo);
            qdf_mul(Qrehihi_h[j][i],Qrelohi_h[j][i],
                    Qrehilo_h[j][i],Qrelolo_h[j][i],
                    rhsimhihi[j],rhsimlohi[j],rhsimhilo[j],rhsimlolo[j],
                    &acc4hihi,&acc4lohi,&acc4hilo,&acc4lolo);
            qdf_inc(&acc1hihi,&acc1lohi,&acc1hilo,&acc1lolo,
                     acc2hihi, acc2lohi, acc2hilo, acc2lolo);
            qdf_inc(&qHrhsrehihi[i],&qHrhsrelohi[i],
                    &qHrhsrehilo[i],&qHrhsrelolo[i],
                    acc1hihi,acc1lohi,acc1hilo,acc1lolo);
            qdf_dec(&acc4hihi,&acc4lohi,&acc4hilo,&acc4lolo,
                     acc3hihi, acc3lohi, acc3hilo, acc3lolo);
            qdf_inc(&qHrhsimhihi[i],&qHrhsimlohi[i],
                    &qHrhsimhilo[i],&qHrhsimlolo[i],
                    acc4hihi,acc4lohi,acc4hilo,acc4lolo);
         }
      }
      cout << "-> CPU solves an upper triangular system ..." << endl;

      CPU_cmplx4_upper_tiled_solver
         (ncols,sizetile,numtiles,
              Rrehihi_h,  Rrelohi_h,  Rrehilo_h,  Rrelolo_h,
              Rimhihi_h,  Rimlohi_h,  Rimhilo_h,  Rimlolo_h,
          qHrhsrehihi,qHrhsrelohi,qHrhsrehilo,qHrhsrelolo,
          qHrhsimhihi,qHrhsimlohi,qHrhsimhilo,qHrhsimlolo,
              xrehihi,    xrelohi,    xrehilo,    xrelolo,
              ximhihi,    ximlohi,    ximhilo,    ximlolo,&bstimelapsed_h);

      if(verbose > 0)
      {
         cout << "CPU solution computed with tiling :" << endl;
         cout << scientific << setprecision(16);
         for(int i=0; i<ncols; i++)
         {
            cout << "x[" << i << "]re : "
                 << xrehihi[i] << "  " << xrelohi[i] << endl
                 << xrehilo[i] << "  " << xrelolo[i] << endl;
            cout << "x[" << i << "]im : "
                 << ximhihi[i] << "  " << ximlohi[i] << endl
                 << ximhilo[i] << "  " << ximlolo[i] << endl;
         }
      }
      cout << scientific << setprecision(2);
      cout << "   Sum of CPU errors on solution : "
           << cmplx4_Difference_Sum
                 (ncols,solrehihi,solrelohi,solrehilo,solrelolo,
                        solimhihi,solimlohi,solimhilo,solimlolo,
                          xrehihi,  xrelohi,  xrehilo,  xrelolo,
                          ximhihi,  ximlohi,  ximhilo,  ximlolo) << endl;

      free(xrehihi); free(xrelohi); free(xrehilo); free(xrelolo);
      free(ximhihi); free(ximlohi); free(ximhilo); free(ximlolo);
      free(qHrhsrehihi); free(qHrhsrelohi);
      free(qHrhsrehilo); free(qHrhsrelolo);
      free(qHrhsimhihi); free(qHrhsimlohi);
      free(qHrhsimhilo); free(qHrhsimlolo);
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
   double addover = 0.0;
   double mulover = 0.0;
   double divover = 0.0;

   if((mode == 0) || (mode == 2))
   {
      double *xrehihi_d = new double[ncols];
      double *xrelohi_d = new double[ncols];
      double *xrehilo_d = new double[ncols];
      double *xrelolo_d = new double[ncols];
      double *ximhihi_d = new double[ncols];
      double *ximlohi_d = new double[ncols];
      double *ximhilo_d = new double[ncols];
      double *ximlolo_d = new double[ncols];
      double *qHrhsrehihi_d = new double[nrows];
      double *qHrhsrelohi_d = new double[nrows];
      double *qHrhsrehilo_d = new double[nrows];
      double *qHrhsrelolo_d = new double[nrows];
      double *qHrhsimhihi_d = new double[nrows];
      double *qHrhsimlohi_d = new double[nrows];
      double *qHrhsimhilo_d = new double[nrows];
      double *qHrhsimlolo_d = new double[nrows];

      cout << "-> GPU computes the blocked Householder QR ..." << endl;

      if(verbose > 0) // to verify that A has not changed ...
      {
         cout << scientific << setprecision(16);
 
         cout << "A random matrix :" << endl;
         for(int i=0; i<nrows; i++)
            for(int j=0; j<ncols; j++)
            {
               cout << "A[" << i << "][" << j << "]re : "
                    << Arehihi[i][j] << "  " << Arelohi[i][j] << endl
                    << Arehilo[i][j] << "  " << Arelolo[i][j] << endl;
               cout << "A[" << i << "][" << j << "]im : "
                    << Aimhihi[i][j] << "  " << Aimlohi[i][j] << endl
                    << Aimhilo[i][j] << "  " << Aimlolo[i][j] << endl;
            }
      }
      GPU_cmplx4_blocked_houseqr
         (nrows,ncols,sizetile,numtiles,
          Arehihi,  Arelohi,  Arehilo,  Arelolo,
          Aimhihi,  Aimlohi,  Aimhilo,  Aimlolo,
          Qrehihi_d,Qrelohi_d,Qrehilo_d,Qrelolo_d,
          Qimhihi_d,Qimlohi_d,Qimhilo_d,Qimlolo_d,
          Rrehihi_d,Rrelohi_d,Rrehilo_d,Rrelolo_d,
          Rimhihi_d,Rimlohi_d,Rimhilo_d,Rimlolo_d,
          &houselapsedms,&RHvlapsedms,&tileRlapsedms,&vb2Wlapsedms,
          &WYTlapsedms,&QWYTlapsedms,&Qaddlapsedms,
          &YWTlapsedms,&YWTClapsedms,&Raddlapsedms,&qrtimelapsed_d,
          &qraddcnt,&qrmulcnt,&qrdivcnt,&sqrtcnt,bvrb);

      cout << "-> Testing the QR factorization ..." << endl;

      // fail = test_cmplx2_qr_factors
      //           (nrows,ncols,Arehi,  Arelo,  Aimhi,  Aimlo,
      //                        Qrehi_d,Qrelo_d,Qimhi_d,Qimlo_d,
      //                        Rrehi_d,Rrelo_d,Rimhi_d,Rimlo_d,tol,verbose);
      fail = test_cmplx4_qr_factors_probe
                (nrows,ncols,
                 Arehihi,  Arelohi,  Arehilo,  Arelolo,
                 Aimhihi,  Aimlohi,  Aimhilo,  Aimlolo,
                 Qrehihi_d,Qrelohi_d,Qrehilo_d,Qrelolo_d,
                 Qimhihi_d,Qimlohi_d,Qimhilo_d,Qimlolo_d,
                 Rrehihi_d,Rrelohi_d,Rrehilo_d,Rrelolo_d,
                 Rimhihi_d,Rimlohi_d,Rimhilo_d,Rimlolo_d,tol,2,true);
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
         qHrhsrehihi_d[i] = 0.0; qHrhsrelohi_d[i] = 0.0;
         qHrhsrehilo_d[i] = 0.0; qHrhsrelolo_d[i] = 0.0;
         qHrhsimhihi_d[i] = 0.0; qHrhsimlohi_d[i] = 0.0;
         qHrhsimhilo_d[i] = 0.0; qHrhsimlolo_d[i] = 0.0;

         for(int j=0; j<nrows; j++) // qHrhs[i] = qHrhs[i] + Q[j][i]*rhs[j];
         {
            qdf_mul(Qrehihi_d[j][i],Qrelohi_d[j][i],
                    Qrehilo_d[j][i],Qrelolo_d[j][i],
                    rhsrehihi[j],rhsrelohi[j],rhsrehilo[j],rhsrelolo[j],
                    &acc1hihi,   &acc1lohi,   &acc1hilo,   &acc1lolo);
            qdf_mul(Qimhihi_d[j][i],Qimlohi_d[j][i],
                    Qimhilo_d[j][i],Qimlolo_d[j][i],
                    rhsimhihi[j],rhsimlohi[j],rhsimhilo[j],rhsimlolo[j],
                    &acc2hihi,   &acc2lohi,   &acc2hilo,   &acc2lolo);
            qdf_mul(Qimhihi_d[j][i],Qimlohi_d[j][i],
                    Qimhilo_d[j][i],Qimlolo_d[j][i],
                    rhsrehihi[j],rhsrelohi[j],rhsrehilo[j],rhsrelolo[j],
                    &acc3hihi,   &acc3lohi,   &acc3hilo,   &acc3lolo);
            qdf_mul(Qrehihi_d[j][i],Qrelohi_d[j][i],
                    Qrehilo_d[j][i],Qrelolo_d[j][i],
                    rhsimhihi[j],rhsimlohi[j],rhsimhilo[j],rhsimlolo[j],
                    &acc4hihi,   &acc4lohi,   &acc4hilo,   &acc4lolo);
            qdf_inc(&acc1hihi,&acc1lohi,&acc1hilo,&acc1lolo,
                     acc2hihi, acc2lohi, acc2hilo, acc2lolo);
            qdf_inc(&qHrhsrehihi_d[i],&qHrhsrelohi_d[i],
                    &qHrhsrehilo_d[i],&qHrhsrelolo_d[i],
                    acc1hihi,acc1lohi,acc1hilo,acc1lolo);
            qdf_dec(&acc4hihi,&acc4lohi,&acc4hilo,&acc4lolo,
                     acc3hihi, acc3lohi, acc3hilo, acc3lolo);
            qdf_inc(&qHrhsimhihi_d[i],&qHrhsimlohi_d[i],
                    &qHrhsimhilo_d[i],&qHrhsimlolo_d[i],
                    acc4hihi,acc4lohi,acc4hilo,acc4lolo);
         }
      }
      cout << "-> GPU solves an upper triangular system ..." << endl;

      GPU_cmplx4_upper_tiled_solver
         (ncols,sizetile,numtiles,
              Rrehihi_d,    Rrelohi_d,    Rrehilo_d,    Rrelolo_d,
              Rimhihi_d,    Rimlohi_d,    Rimhilo_d,    Rimlolo_d,
          qHrhsrehihi_d,qHrhsrelohi_d,qHrhsrehilo_d,qHrhsrelolo_d,
          qHrhsimhihi_d,qHrhsimlohi_d,qHrhsimhilo_d,qHrhsimlolo_d,
              xrehihi_d,    xrelohi_d,    xrehilo_d,    xrelolo_d,
              ximhihi_d,    ximlohi_d,    ximhilo_d,    ximlolo_d,
          &invlapsed,&mullapsed,&sublapsed,&bselapsedms,&bstimelapsed_d,
          &bsaddcnt,&addover,&bsmulcnt,&mulover,&bsdivcnt,&divover);

      if(verbose > 0)
      {
         cout << scientific << setprecision(16);
         cout << "GPU solution computed with tiling :" << endl;
         for(int i=0; i<ncols; i++)
         {
            cout << "x[" << i << "]re : "
                 << xrehihi_d[i] << "  " << xrelohi_d[i] << endl
                 << xrehilo_d[i] << "  " << xrelolo_d[i] << endl;
            cout << "x[" << i << "]im : "
                 << ximhihi_d[i] << "  " << ximlohi_d[i] << endl
                 << ximhilo_d[i] << "  " << ximlolo_d[i] << endl;
         }
      }
      cout << scientific << setprecision(2);
      cout << "   Sum of GPU errors on solution : "
           << cmplx4_Difference_Sum
                 (ncols,solrehihi,solrelohi,solrehilo,solrelolo,
                        solimhihi,solimlohi,solimhilo,solimlolo,
                          xrehihi_d,xrelohi_d,xrehilo_d,xrelolo_d,
                          ximhihi_d,ximlohi_d,ximhilo_d,ximlolo_d) << endl;
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
           << qraddcnt << " x 89 " << endl;
      cout << "                    Number of multiplications : "
           << qrmulcnt << " x 336 " << endl;
      cout << "                          Number of divisions : "
           << qrdivcnt << " x 893 " << endl;
      cout << "                    Number of calls to sqrt() : "
           << sqrtcnt << " x 1345 " << endl;
      long long int qrflopcnt
          = 89*qraddcnt + 336*qrmulcnt + 893*qrdivcnt + 1345*sqrtcnt;
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
           << bsaddcnt << " x 89 " << endl;
      cout << "                    Number of multiplications : "
           << bsmulcnt << " x 336 " << endl;
      cout << "                          Number of divisions : "
           << bsdivcnt << " x 893 " << endl;
      long long int bsflopcnt = 89*bsaddcnt + 336*bsmulcnt + 893*bsdivcnt;
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
      free(Arehihi[i]); free(Arelohi[i]);
      free(Arehilo[i]); free(Arelolo[i]);
      free(Qrehihi_h[i]); free(Qrelohi_h[i]);
      free(Qrehilo_h[i]); free(Qrelolo_h[i]);
      free(Qrehihi_d[i]); free(Qrelohi_d[i]);
      free(Qrehilo_d[i]); free(Qrelolo_d[i]); 
      free(Rrehihi_h[i]); free(Rrelohi_h[i]);
      free(Rrehilo_h[i]); free(Rrelolo_h[i]);
      free(Rrehihi_d[i]); free(Rrelohi_d[i]);
      free(Rrehilo_d[i]); free(Rrelolo_d[i]);
      free(Aimhihi[i]); free(Aimlohi[i]);
      free(Aimhilo[i]); free(Aimlolo[i]);
      free(Qimhihi_h[i]); free(Qimlohi_h[i]);
      free(Qimhilo_h[i]); free(Qimlolo_h[i]);
      free(Qimhihi_d[i]); free(Qimlohi_d[i]);
      free(Qimhilo_d[i]); free(Qimlolo_d[i]);
      free(Rimhihi_h[i]); free(Rimlohi_h[i]);
      free(Rimhilo_h[i]); free(Rimlolo_h[i]);
      free(Rimhihi_d[i]); free(Rimlohi_d[i]);
      free(Rimhilo_d[i]); free(Rimlolo_d[i]);
   }
   free(Arehihi); free(Arelohi); free(Arehilo); free(Arelolo);
   free(Qrehihi_h); free(Qrelohi_h); free(Qrehilo_h); free(Qrelolo_h);
   free(Rrehihi_h); free(Rrelohi_h); free(Rrehilo_h); free(Rrelolo_h);
   free(Qrehihi_d); free(Qrelohi_d); free(Qrehilo_d); free(Qrelolo_d);
   free(Rrehihi_d); free(Rrelohi_d); free(Rrehilo_d); free(Rrelolo_d);
   free(Aimhihi); free(Aimlohi); free(Aimhilo); free(Aimlolo);
   free(Qimhihi_h); free(Qimlohi_h); free(Qimhilo_h); free(Qimlolo_h);
   free(Rimhihi_h); free(Rimlohi_h); free(Rimhilo_h); free(Rimlolo_h);
   free(Qimhihi_d); free(Qimlohi_d); free(Qimhilo_d); free(Qimlolo_d);
   free(Rimhihi_d); free(Rimlohi_d); free(Rimhilo_d); free(Rimlolo_d);
   free(solrehihi); free(solrelohi); free(solrehilo); free(solrelolo);
   free(solimhihi); free(solimlohi); free(solimhilo); free(solimlolo);
   free(rhsrehihi); free(rhsrelohi); free(rhsrehilo); free(rhsrelolo);
   free(rhsimhihi); free(rhsimlohi); free(rhsimhilo); free(rhsimlolo);
}
