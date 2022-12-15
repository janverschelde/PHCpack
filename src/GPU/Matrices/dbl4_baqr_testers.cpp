// The file dbl4_baqr_testers.cpp defines the function with prototypes in
// the file dbl4_baqr_testers.h.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#include "random4_matrices.h"
#include "dbl4_factorizations.h"
#include "dbl4_factors_testers.h"
#include "dbl4_baqr_host.h"
#include "dbl4_baqr_kernels.h"
#include "write_dbl4_qrtimeflops.h"
#include "dbl4_baqr_testers.h"

using namespace std;

void test_real4_blocked_qr
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
                 << "          "
                 << Ahilo[i][j] << "  " << Alolo[i][j] << endl;
   }
   double timelapsed_h;
   bool bvrb = (verbose > 0);
   const double tol = 1.0e-40;
   int fail;

   if((mode == 1) || (mode == 2))
   {

      cout << "-> CPU computes the blocked Householder QR ..." << endl;

      CPU_dbl4_blocked_houseqr
         (nrows,ncols,sizetile,numtiles,Ahihi,Alohi,Ahilo,Alolo,
          Qhihi_h,Qlohi_h,Qhilo_h,Qlolo_h,Rhihi_h,Rlohi_h,Rhilo_h,Rlolo_h,
          &timelapsed_h,bvrb);

      cout << "-> Testing the QR factorization ..." << endl;

      fail = test_real4_qr_factors_probe
         (nrows,ncols,Ahihi,Alohi,Ahilo,Alolo,
          Qhihi_h,Qlohi_h,Qhilo_h,Qlolo_h,Rhihi_h,Rlohi_h,Rhilo_h,Rlolo_h,
          tol,2,true);

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

      GPU_dbl4_blocked_houseqr
         (nrows,ncols,sizetile,numtiles,
          Ahihi,  Alohi,  Ahilo,  Alolo,
          Qhihi_d,Qlohi_d,Qhilo_d,Qlolo_d,
          Rhihi_d,Rlohi_d,Rhilo_d,Rlolo_d,
          &houselapsedms,&RTvlapsedms,&tileRlapsedms,&vb2Wlapsedms,
          &WYTlapsedms,&QWYTlapsedms,&Qaddlapsedms,
          &YWTlapsedms,&YWTClapsedms,&Raddlapsedms,&timelapsed_d,
          &addcnt,&mulcnt,&divcnt,&sqrtcnt,bvrb);

      cout << "-> Testing the QR factorization ..." << endl;

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
      write_dbl4_qrtimeflops
         (0,nrows,ncols,houselapsedms,RTvlapsedms,tileRlapsedms,vb2Wlapsedms,
          WYTlapsedms,QWYTlapsedms,Qaddlapsedms,YWTlapsedms,YWTClapsedms,
          Raddlapsedms,timelapsed_d,addcnt,mulcnt,divcnt,sqrtcnt);
   }
   for(int i=0; i<nrows; i++)
   {
      free(Ahihi[i]);   free(Alohi[i]);   free(Ahilo[i]);   free(Alolo[i]);
      free(Qhihi_h[i]); free(Qlohi_h[i]); free(Qhilo_h[i]); free(Qlolo_h[i]);
      free(Qhihi_d[i]); free(Qlohi_d[i]); free(Qhilo_d[i]); free(Qlolo_d[i]);
      free(Rhihi_h[i]); free(Rlohi_h[i]); free(Rhilo_h[i]); free(Rlolo_h[i]);
      free(Rhihi_d[i]); free(Rlohi_d[i]); free(Rhilo_d[i]); free(Rlolo_d[i]);
   }
   free(Ahihi);   free(Alohi);   free(Ahilo);   free(Alolo);
   free(Qhihi_h); free(Qlohi_h); free(Qhilo_h); free(Qlolo_h);
   free(Qhihi_d); free(Qlohi_d); free(Qhilo_d); free(Qlolo_d);
   free(Rhihi_h); free(Rlohi_h); free(Rhilo_h); free(Rlolo_h);
   free(Rhihi_d); free(Rlohi_d); free(Rhilo_d); free(Rlolo_d);
}

void test_cmplx4_blocked_qr
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

/*
   for(int i=0; i<nrows; i++)      // simplify input to real data
      for(int j=0; j<ncols; j++)
      {
         Aimhihi[i][j] = 0.0;
         Aimlohi[i][j] = 0.0;
         Aimhilo[i][j] = 0.0;
         Aimlolo[i][j] = 0.0;
         // Arelohi[i][j] = 0.0;
         // Arehilo[i][j] = 0.0;
         // Arelolo[i][j] = 0.0;
      }
 */
   if(verbose > 0)
   {
      cout << scientific << setprecision(16);

      cout << "A random matrix :" << endl;
      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++)
         {
            cout << "A[" << i << "][" << j << "]re : "
                 << Arehihi[i][j] << "  " << Arelohi[i][j] << endl
                 << "            "
                 << Arehilo[i][j] << "  " << Arelolo[i][j] << endl;
            cout << "A[" << i << "][" << j << "]im : "
                 << Aimhihi[i][j] << "  " << Aimlohi[i][j] << endl
                 << "            "
                 << Aimhilo[i][j] << "  " << Aimlolo[i][j] << endl;
         }
   }
   double timelapsed_h;
   bool bvrb = (verbose > 0);
   const double tol = 1.0e-40;
   int fail;

   if((mode == 1) || (mode == 2))
   {
      cout << "-> CPU computes the blocked Householder QR ..." << endl;

      CPU_cmplx4_blocked_houseqr
         (nrows,ncols,sizetile,numtiles,
          Arehihi,  Arelohi,  Arehilo,  Arelolo,
          Aimhihi,  Aimlohi,  Aimhilo,  Aimlolo,
          Qrehihi_h,Qrelohi_h,Qrehilo_h,Qrelolo_h,
          Qimhihi_h,Qimlohi_h,Qimhilo_h,Qimlolo_h,
          Rrehihi_h,Rrelohi_h,Rrehilo_h,Rrelolo_h,
          Rimhihi_h,Rimlohi_h,Rimhilo_h,Rimlolo_h,&timelapsed_h,bvrb);

      cout << "-> Testing the QR factorization ..." << endl;
 /*
      fail = test_cmplx4_qr_factors
         (nrows,ncols,Arehihi,  Arelohi,  Arehilo,  Arelolo,
                      Aimhihi,  Aimlohi,  Aimhilo,  Aimlolo,
                      Qrehihi_h,Qrelohi_h,Qrehilo_h,Qrelolo_h,
                      Qimhihi_h,Qimlohi_h,Qimhilo_h,Qimlolo_h,
                      Rrehihi_h,Rrelohi_h,Rrehilo_h,Rrelolo_h,
                      Rimhihi_h,Rimlohi_h,Rimhilo_h,Rimlolo_h,tol,bvrb);
  */
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
                    << Arehihi[i][j] << "  " << Arelohi[i][j] << endl
                    << "            "
                    << Arehilo[i][j] << "  " << Arelolo[i][j] << endl;
               cout << "A[" << i << "][" << j << "]im : "
                    << Aimhihi[i][j] << "  " << Aimlohi[i][j] << endl
                    << "            "
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
          &YWTlapsedms,&YWTClapsedms,&Raddlapsedms,&timelapsed_d,
          &addcnt,&mulcnt,&divcnt,&sqrtcnt,bvrb);

      cout << "-> Testing the QR factorization ..." << endl;
  /*
      fail = test_cmplx4_qr_factors
                (nrows,ncols,Arehihi,  Arelohi,  Arehilo,  Arelolo,
                             Aimhihi,  Aimlohi,  Aimhilo,  Aimlolo,
                             Qrehihi_d,Qrelohi_d,Qrehilo_d,Qrelolo_d,
                             Qimhihi_d,Qimlohi_d,Qimhilo_d,Qimlolo_d,
                             Rrehihi_d,Rrelohi_d,Rrehilo_d,Rrelolo_d,
                             Rimhihi_d,Rimlohi_d,Rimhilo_d,Rimlolo_d,
                 tol,bvrb);
   */
      fail = test_cmplx4_qr_factors_probe
                (nrows,ncols,Arehihi,  Arelohi,  Arehilo,  Arelolo,
                             Aimhihi,  Aimlohi,  Aimhilo,  Aimlolo,
                             Qrehihi_d,Qrelohi_d,Qrehilo_d,Qrelolo_d,
                             Qimhihi_d,Qimlohi_d,Qimhilo_d,Qimlolo_d,
                             Rrehihi_d,Rrelohi_d,Rrehilo_d,Rrelolo_d,
                             Rimhihi_d,Rimlohi_d,Rimhilo_d,Rimlolo_d,
                 tol,2,true);
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
      write_dbl4_qrtimeflops
         (1,nrows,ncols,houselapsedms,RHvlapsedms,tileRlapsedms,vb2Wlapsedms,
          WYTlapsedms,QWYTlapsedms,Qaddlapsedms,YWTlapsedms,YWTClapsedms,
          Raddlapsedms,timelapsed_d,addcnt,mulcnt,divcnt,sqrtcnt);
   }
   for(int i=0; i<nrows; i++)
   {
      free(Arehihi[i]); free(Arelohi[i]); free(Arehilo[i]); free(Arelolo[i]);
      free(Aimhihi[i]); free(Aimlohi[i]); free(Aimhilo[i]); free(Aimlolo[i]);
      free(Qrehihi_h[i]); free(Qrelohi_h[i]);
      free(Qrehilo_h[i]); free(Qrelolo_h[i]);
      free(Qrehihi_d[i]); free(Qrelohi_d[i]);
      free(Qrehilo_d[i]); free(Qrelolo_d[i]);
      free(Qimhihi_h[i]); free(Qimlohi_h[i]);
      free(Qimhilo_h[i]); free(Qimlolo_h[i]);
      free(Qimhihi_d[i]); free(Qimlohi_d[i]);
      free(Qimhilo_d[i]); free(Qimlolo_d[i]);
      free(Rrehihi_h[i]); free(Rrelohi_h[i]);
      free(Rrehilo_h[i]); free(Rrelolo_h[i]);
      free(Rrehihi_d[i]); free(Rrelohi_d[i]);
      free(Rrehilo_d[i]); free(Rrelolo_d[i]);
      free(Rimhihi_h[i]); free(Rimlohi_h[i]);
      free(Rimhilo_h[i]); free(Rimlolo_h[i]);
      free(Rimhihi_d[i]); free(Rimlohi_d[i]);
      free(Rimhilo_d[i]); free(Rimlolo_d[i]);
   }
   free(Arehihi);   free(Arelohi);   free(Arehilo);   free(Arelolo);
   free(Aimhihi);   free(Aimlohi);   free(Aimhilo);   free(Aimlolo);
   free(Qrehihi_h); free(Qrelohi_h); free(Qrehilo_h); free(Qrelolo_h);
   free(Qrehihi_d); free(Qrelohi_d); free(Qrehilo_d); free(Qrelolo_d);
   free(Rrehihi_h); free(Rrelohi_h); free(Rrehilo_h); free(Rrelolo_h);
   free(Rrehihi_d); free(Rrelohi_d); free(Rrehilo_d); free(Rrelolo_d);
   free(Qimhihi_h); free(Qimlohi_h); free(Qimhilo_h); free(Qimlolo_h);
   free(Qimhihi_d); free(Qimlohi_d); free(Qimhilo_d); free(Qimlolo_d);
   free(Rimhihi_h); free(Rimlohi_h); free(Rimhilo_h); free(Rimlolo_h);
   free(Rimhihi_d); free(Rimlohi_d); free(Rimhilo_d); free(Rimlolo_d);
}
