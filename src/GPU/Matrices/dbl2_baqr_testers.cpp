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
#include "write_dbl2_qrtimeflops.h"
#include "dbl2_baqr_testers.h"

using namespace std;

void test_real2_blocked_qr
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
   double timelapsed_h;
   bool bvrb = (verbose > 0);
   const double tol = 1.0e-20;
   int fail;

   if((mode == 1) || (mode == 2))
   {

      cout << "-> CPU computes the blocked Householder QR ..." << endl;

      CPU_dbl2_blocked_houseqr
         (nrows,ncols,sizetile,numtiles,Ahi,Alo,Qhi_h,Qlo_h,Rhi_h,Rlo_h,
          &timelapsed_h,bvrb);

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
      cout << "Elapsed CPU time (Linux), Wall time (Windows) : "
           << timelapsed_h << " seconds." << endl;
   }
   if((mode == 0) || (mode == 2))
   {
      write_dbl2_qrtimeflops
         (0,nrows,ncols,houselapsedms,RTvlapsedms,tileRlapsedms,vb2Wlapsedms,
          WYTlapsedms,QWYTlapsedms,Qaddlapsedms,YWTlapsedms,YWTClapsedms,
          Raddlapsedms,timelapsed_d,addcnt,mulcnt,divcnt,sqrtcnt);
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
}

void test_cmplx2_blocked_qr
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
   double timelapsed_h;
   bool bvrb = (verbose > 0);
   const double tol = 1.0e-20;
   int fail;

   if((mode == 1) || (mode == 2))
   {
      cout << "-> CPU computes the blocked Householder QR ..." << endl;

      CPU_cmplx2_blocked_houseqr
         (nrows,ncols,sizetile,numtiles,
          Arehi,  Arelo,  Aimhi,  Aimlo,
          Qrehi_h,Qrelo_h,Qimhi_h,Qimlo_h,
          Rrehi_h,Rrelo_h,Rimhi_h,Rimlo_h,&timelapsed_h,bvrb);

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
      cout << "Elapsed CPU time (Linux), Wall time (Windows) : "
           << timelapsed_h << " seconds." << endl;
   }
   if((mode == 0) || (mode == 2))
   {
      write_dbl2_qrtimeflops
         (1,nrows,ncols,houselapsedms,RHvlapsedms,tileRlapsedms,vb2Wlapsedms,
          WYTlapsedms,QWYTlapsedms,Qaddlapsedms,YWTlapsedms,YWTClapsedms,
          Raddlapsedms,timelapsed_d,addcnt,mulcnt,divcnt,sqrtcnt);
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
}
