// The file dbl8_baqr_testers.cpp defines the function with prototypes in
// the file dbl8_baqr_testers.h.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#include "random8_matrices.h"
#include "dbl8_factorizations.h"
#include "dbl8_factors_testers.h"
#include "dbl8_baqr_host.h"
#include "dbl8_baqr_kernels.h"
#include "write_dbl8_qrtimeflops.h"
#include "dbl8_baqr_testers.h"

using namespace std;

void test_real8_blocked_qr
 ( int seed, int szt, int nbt, int nrows, int vrb, int mode )
{
   const int sizetile = szt;
   const int numtiles = nbt;
   const int ncols = sizetile*numtiles;
   const int verbose = vrb;

   cout << "-> Generating a random " << nrows
        << "-by-" << ncols << " matrix ..." << endl;

   double **Ahihihi = new double*[nrows];
   double **Alohihi = new double*[nrows];
   double **Ahilohi = new double*[nrows];
   double **Alolohi = new double*[nrows];
   double **Ahihilo = new double*[nrows];
   double **Alohilo = new double*[nrows];
   double **Ahilolo = new double*[nrows];
   double **Alololo = new double*[nrows];
   double **Qhihihi_h = new double*[nrows];
   double **Qlohihi_h = new double*[nrows];
   double **Qhilohi_h = new double*[nrows];
   double **Qlolohi_h = new double*[nrows];
   double **Qhihilo_h = new double*[nrows];
   double **Qlohilo_h = new double*[nrows];
   double **Qhilolo_h = new double*[nrows];
   double **Qlololo_h = new double*[nrows];
   double **Qhihihi_d = new double*[nrows];
   double **Qlohihi_d = new double*[nrows];
   double **Qhilohi_d = new double*[nrows];
   double **Qlolohi_d = new double*[nrows];
   double **Qhihilo_d = new double*[nrows];
   double **Qlohilo_d = new double*[nrows];
   double **Qhilolo_d = new double*[nrows];
   double **Qlololo_d = new double*[nrows];
   double **Rhihihi_h = new double*[nrows];
   double **Rlohihi_h = new double*[nrows];
   double **Rhilohi_h = new double*[nrows];
   double **Rlolohi_h = new double*[nrows];
   double **Rhihilo_h = new double*[nrows];
   double **Rlohilo_h = new double*[nrows];
   double **Rhilolo_h = new double*[nrows];
   double **Rlololo_h = new double*[nrows];
   double **Rhihihi_d = new double*[nrows];
   double **Rlohihi_d = new double*[nrows];
   double **Rhilohi_d = new double*[nrows];
   double **Rlolohi_d = new double*[nrows];
   double **Rhihilo_d = new double*[nrows];
   double **Rlohilo_d = new double*[nrows];
   double **Rhilolo_d = new double*[nrows];
   double **Rlololo_d = new double*[nrows];

   for(int i=0; i<nrows; i++)
   {
      Ahihihi[i] = new double[ncols];
      Alohihi[i] = new double[ncols];
      Ahilohi[i] = new double[ncols];
      Alolohi[i] = new double[ncols];
      Ahihilo[i] = new double[ncols];
      Alohilo[i] = new double[ncols];
      Ahilolo[i] = new double[ncols];
      Alololo[i] = new double[ncols];
      Qhihihi_h[i] = new double[nrows];
      Qlohihi_h[i] = new double[nrows];
      Qhilohi_h[i] = new double[nrows];
      Qlolohi_h[i] = new double[nrows];
      Qhihilo_h[i] = new double[nrows];
      Qlohilo_h[i] = new double[nrows];
      Qhilolo_h[i] = new double[nrows];
      Qlololo_h[i] = new double[nrows];
      Qhihihi_d[i] = new double[nrows];
      Qlohihi_d[i] = new double[nrows];
      Qhilohi_d[i] = new double[nrows];
      Qlolohi_d[i] = new double[nrows];
      Qhihilo_d[i] = new double[nrows];
      Qlohilo_d[i] = new double[nrows];
      Qhilolo_d[i] = new double[nrows];
      Qlololo_d[i] = new double[nrows];
      Rhihihi_h[i] = new double[ncols];
      Rlohihi_h[i] = new double[ncols];
      Rhilohi_h[i] = new double[ncols];
      Rlolohi_h[i] = new double[ncols];
      Rhihilo_h[i] = new double[ncols];
      Rlohilo_h[i] = new double[ncols];
      Rhilolo_h[i] = new double[ncols];
      Rlololo_h[i] = new double[ncols];
      Rhihihi_d[i] = new double[ncols];
      Rlohihi_d[i] = new double[ncols];
      Rhilohi_d[i] = new double[ncols];
      Rlolohi_d[i] = new double[ncols];
      Rhihilo_d[i] = new double[ncols];
      Rlohilo_d[i] = new double[ncols];
      Rhilolo_d[i] = new double[ncols];
      Rlololo_d[i] = new double[ncols];
   }
   random_dbl8_matrix(nrows,ncols,Ahihihi,Alohihi,Ahilohi,Alolohi,
                                  Ahihilo,Alohilo,Ahilolo,Alololo);

   if(verbose > 0)
   {
      cout << scientific << setprecision(16);

      cout << "A random matrix :" << endl;
      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++)
            cout << "A[" << i << "][" << j << "] : "
                 << Ahihihi[i][j] << "  " << Alohihi[i][j] << endl
                 << "          "
                 << Ahilohi[i][j] << "  " << Alolohi[i][j] << endl
                 << "          "
                 << Ahihilo[i][j] << "  " << Alohilo[i][j] << endl
                 << "          "
                 << Ahilolo[i][j] << "  " << Alololo[i][j] << endl;
   }
   double timelapsed_h;
   bool bvrb = (verbose > 0);
   const double tol = 1.0e-40;
   int fail;

   if((mode == 1) || (mode == 2))
   {

      cout << "-> CPU computes the blocked Householder QR ..." << endl;

      CPU_dbl8_blocked_houseqr
         (nrows,ncols,sizetile,numtiles,
          Ahihihi,Alohihi,Ahilohi,Alolohi,
          Ahihilo,Alohilo,Ahilolo,Alololo,
          Qhihihi_h,Qlohihi_h,Qhilohi_h,Qlolohi_h,
          Qhihilo_h,Qlohilo_h,Qhilolo_h,Qlololo_h,
          Rhihihi_h,Rlohihi_h,Rhilohi_h,Rlolohi_h,
          Rhihilo_h,Rlohilo_h,Rhilolo_h,Rlololo_h,
          &timelapsed_h,bvrb);

      cout << "-> Testing the QR factorization ..." << endl;

      fail = test_real8_qr_factors_probe
         (nrows,ncols,
          Ahihihi,Alohihi,Ahilohi,Alolohi,
          Ahihilo,Alohilo,Ahilolo,Alololo,
          Qhihihi_h,Qlohihi_h,Qhilohi_h,Qlolohi_h,
          Qhihilo_h,Qlohilo_h,Qhilolo_h,Qlololo_h,
          Rhihihi_h,Rlohihi_h,Rhilohi_h,Rlolohi_h,
          Rhihilo_h,Rlohilo_h,Rhilolo_h,Rlololo_h,tol,2,true);

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

      GPU_dbl8_blocked_houseqr
         (nrows,ncols,sizetile,numtiles,
          Ahihihi,  Alohihi,  Ahilohi,  Alolohi,
          Ahihilo,  Alohilo,  Ahilolo,  Alololo,
          Qhihihi_d,Qlohihi_d,Qhilohi_d,Qlolohi_d,
          Qhihilo_d,Qlohilo_d,Qhilolo_d,Qlololo_d,
          Rhihihi_d,Rlohihi_d,Rhilohi_d,Rlolohi_d,
          Rhihilo_d,Rlohilo_d,Rhilolo_d,Rlololo_d,
          &houselapsedms,&RTvlapsedms,&tileRlapsedms,&vb2Wlapsedms,
          &WYTlapsedms,&QWYTlapsedms,&Qaddlapsedms,
          &YWTlapsedms,&YWTClapsedms,&Raddlapsedms,&timelapsed_d,
          &addcnt,&mulcnt,&divcnt,&sqrtcnt,bvrb);

      cout << "-> Testing the QR factorization ..." << endl;

      fail = test_real8_qr_factors_probe
                (nrows,ncols,
                 Ahihihi,  Alohihi,  Ahilohi,  Alolohi,
                 Ahihilo,  Alohilo,  Ahilolo,  Alololo,
                 Qhihihi_d,Qlohihi_d,Qhilohi_d,Qlolohi_d,
                 Qhihilo_d,Qlohilo_d,Qhilolo_d,Qlololo_d,
                 Rhihihi_d,Rlohihi_d,Rhilohi_d,Rlolohi_d,
                 Rhihilo_d,Rlohilo_d,Rhilolo_d,Rlololo_d,tol,2,true);
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
      write_dbl8_qrtimeflops
         (0,nrows,ncols,houselapsedms,RTvlapsedms,tileRlapsedms,vb2Wlapsedms,
          WYTlapsedms,QWYTlapsedms,Qaddlapsedms,YWTlapsedms,YWTClapsedms,
          Raddlapsedms,timelapsed_d,addcnt,mulcnt,divcnt,sqrtcnt);
   }
   for(int i=0; i<nrows; i++)
   {
      free(Ahihihi[i]);   free(Alohihi[i]);
      free(Ahilohi[i]);   free(Alolohi[i]);
      free(Ahihilo[i]);   free(Alohilo[i]);
      free(Ahilolo[i]);   free(Alololo[i]);
      free(Qhihihi_h[i]); free(Qlohihi_h[i]);
      free(Qhilohi_h[i]); free(Qlolohi_h[i]);
      free(Qhihilo_h[i]); free(Qlohilo_h[i]);
      free(Qhilolo_h[i]); free(Qlololo_h[i]);
      free(Qhihihi_d[i]); free(Qlohihi_d[i]);
      free(Qhilohi_d[i]); free(Qlolohi_d[i]);
      free(Qhihilo_d[i]); free(Qlohilo_d[i]);
      free(Qhilolo_d[i]); free(Qlololo_d[i]);
      free(Rhihihi_h[i]); free(Rlohihi_h[i]);
      free(Rhilohi_h[i]); free(Rlolohi_h[i]);
      free(Rhihilo_h[i]); free(Rlohilo_h[i]);
      free(Rhilolo_h[i]); free(Rlololo_h[i]);
      free(Rhihihi_d[i]); free(Rlohihi_d[i]);
      free(Rhilohi_d[i]); free(Rlolohi_d[i]);
      free(Rhihilo_d[i]); free(Rlohilo_d[i]);
      free(Rhilolo_d[i]); free(Rlololo_d[i]);
   }
   free(Ahihihi);   free(Alohihi);   free(Ahilohi);   free(Alolohi);
   free(Ahihilo);   free(Alohilo);   free(Ahilolo);   free(Alololo);
   free(Qhihihi_h); free(Qlohihi_h); free(Qhilohi_h); free(Qlolohi_h);
   free(Qhihilo_h); free(Qlohilo_h); free(Qhilolo_h); free(Qlololo_h);
   free(Qhihihi_d); free(Qlohihi_d); free(Qhilohi_d); free(Qlolohi_d);
   free(Qhihilo_d); free(Qlohilo_d); free(Qhilolo_d); free(Qlololo_d);
   free(Rhihihi_h); free(Rlohihi_h); free(Rhilohi_h); free(Rlolohi_h);
   free(Rhihilo_h); free(Rlohilo_h); free(Rhilolo_h); free(Rlololo_h);
   free(Rhihihi_d); free(Rlohihi_d); free(Rhilohi_d); free(Rlolohi_d);
   free(Rhihilo_d); free(Rlohilo_d); free(Rhilolo_d); free(Rlololo_d);
}

void test_cmplx8_blocked_qr
 ( int seed, int szt, int nbt, int nrows, int vrb, int mode )
{
   const int sizetile = szt;
   const int numtiles = nbt;
   const int ncols = sizetile*numtiles;
   const int verbose = vrb;

   cout << "-> Generating a random " << nrows
        << "-by-" << ncols << " matrix ..." << endl;

   double **Arehihihi = new double*[nrows];
   double **Arelohihi = new double*[nrows];
   double **Arehilohi = new double*[nrows];
   double **Arelolohi = new double*[nrows];
   double **Arehihilo = new double*[nrows];
   double **Arelohilo = new double*[nrows];
   double **Arehilolo = new double*[nrows];
   double **Arelololo = new double*[nrows];
   double **Aimhihihi = new double*[nrows];
   double **Aimlohihi = new double*[nrows];
   double **Aimhilohi = new double*[nrows];
   double **Aimlolohi = new double*[nrows];
   double **Aimhihilo = new double*[nrows];
   double **Aimlohilo = new double*[nrows];
   double **Aimhilolo = new double*[nrows];
   double **Aimlololo = new double*[nrows];
   double **Qrehihihi_h = new double*[nrows];
   double **Qrelohihi_h = new double*[nrows];
   double **Qrehilohi_h = new double*[nrows];
   double **Qrelolohi_h = new double*[nrows];
   double **Qrehihilo_h = new double*[nrows];
   double **Qrelohilo_h = new double*[nrows];
   double **Qrehilolo_h = new double*[nrows];
   double **Qrelololo_h = new double*[nrows];
   double **Qimhihihi_h = new double*[nrows];
   double **Qimlohihi_h = new double*[nrows];
   double **Qimhilohi_h = new double*[nrows];
   double **Qimlolohi_h = new double*[nrows];
   double **Qimhihilo_h = new double*[nrows];
   double **Qimlohilo_h = new double*[nrows];
   double **Qimhilolo_h = new double*[nrows];
   double **Qimlololo_h = new double*[nrows];
   double **Qrehihihi_d = new double*[nrows];
   double **Qrelohihi_d = new double*[nrows];
   double **Qrehilohi_d = new double*[nrows];
   double **Qrelolohi_d = new double*[nrows];
   double **Qrehihilo_d = new double*[nrows];
   double **Qrelohilo_d = new double*[nrows];
   double **Qrehilolo_d = new double*[nrows];
   double **Qrelololo_d = new double*[nrows];
   double **Qimhihihi_d = new double*[nrows];
   double **Qimlohihi_d = new double*[nrows];
   double **Qimhilohi_d = new double*[nrows];
   double **Qimlolohi_d = new double*[nrows];
   double **Qimhihilo_d = new double*[nrows];
   double **Qimlohilo_d = new double*[nrows];
   double **Qimhilolo_d = new double*[nrows];
   double **Qimlololo_d = new double*[nrows];
   double **Rrehihihi_h = new double*[nrows];
   double **Rrelohihi_h = new double*[nrows];
   double **Rrehilohi_h = new double*[nrows];
   double **Rrelolohi_h = new double*[nrows];
   double **Rrehihilo_h = new double*[nrows];
   double **Rrelohilo_h = new double*[nrows];
   double **Rrehilolo_h = new double*[nrows];
   double **Rrelololo_h = new double*[nrows];
   double **Rimhihihi_h = new double*[nrows];
   double **Rimlohihi_h = new double*[nrows];
   double **Rimhilohi_h = new double*[nrows];
   double **Rimlolohi_h = new double*[nrows];
   double **Rimhihilo_h = new double*[nrows];
   double **Rimlohilo_h = new double*[nrows];
   double **Rimhilolo_h = new double*[nrows];
   double **Rimlololo_h = new double*[nrows];
   double **Rrehihihi_d = new double*[nrows];
   double **Rrelohihi_d = new double*[nrows];
   double **Rrehilohi_d = new double*[nrows];
   double **Rrelolohi_d = new double*[nrows];
   double **Rrehihilo_d = new double*[nrows];
   double **Rrelohilo_d = new double*[nrows];
   double **Rrehilolo_d = new double*[nrows];
   double **Rrelololo_d = new double*[nrows];
   double **Rimhihihi_d = new double*[nrows];
   double **Rimlohihi_d = new double*[nrows];
   double **Rimhilohi_d = new double*[nrows];
   double **Rimlolohi_d = new double*[nrows];
   double **Rimhihilo_d = new double*[nrows];
   double **Rimlohilo_d = new double*[nrows];
   double **Rimhilolo_d = new double*[nrows];
   double **Rimlololo_d = new double*[nrows];

   for(int i=0; i<nrows; i++)
   {
      Arehihihi[i] = new double[ncols];
      Arelohihi[i] = new double[ncols];
      Arehilohi[i] = new double[ncols];
      Arelolohi[i] = new double[ncols];
      Arehihilo[i] = new double[ncols];
      Arelohilo[i] = new double[ncols];
      Arehilolo[i] = new double[ncols];
      Arelololo[i] = new double[ncols];
      Aimhihihi[i] = new double[ncols];
      Aimlohihi[i] = new double[ncols];
      Aimhilohi[i] = new double[ncols];
      Aimlolohi[i] = new double[ncols];
      Aimhihilo[i] = new double[ncols];
      Aimlohilo[i] = new double[ncols];
      Aimhilolo[i] = new double[ncols];
      Aimlololo[i] = new double[ncols];
      Qrehihihi_h[i] = new double[nrows];
      Qrelohihi_h[i] = new double[nrows];
      Qrehilohi_h[i] = new double[nrows];
      Qrelolohi_h[i] = new double[nrows];
      Qrehihilo_h[i] = new double[nrows];
      Qrelohilo_h[i] = new double[nrows];
      Qrehilolo_h[i] = new double[nrows];
      Qrelololo_h[i] = new double[nrows];
      Qimhihihi_h[i] = new double[nrows];
      Qimlohihi_h[i] = new double[nrows];
      Qimhilohi_h[i] = new double[nrows];
      Qimlolohi_h[i] = new double[nrows];
      Qimhihilo_h[i] = new double[nrows];
      Qimlohilo_h[i] = new double[nrows];
      Qimhilolo_h[i] = new double[nrows];
      Qimlololo_h[i] = new double[nrows];
      Qrehihihi_d[i] = new double[nrows];
      Qrelohihi_d[i] = new double[nrows];
      Qrehilohi_d[i] = new double[nrows];
      Qrelolohi_d[i] = new double[nrows];
      Qrehihilo_d[i] = new double[nrows];
      Qrelohilo_d[i] = new double[nrows];
      Qrehilolo_d[i] = new double[nrows];
      Qrelololo_d[i] = new double[nrows];
      Qimhihihi_d[i] = new double[nrows];
      Qimlohihi_d[i] = new double[nrows];
      Qimhilohi_d[i] = new double[nrows];
      Qimlolohi_d[i] = new double[nrows];
      Qimhihilo_d[i] = new double[nrows];
      Qimlohilo_d[i] = new double[nrows];
      Qimhilolo_d[i] = new double[nrows];
      Qimlololo_d[i] = new double[nrows];
      Rrehihihi_h[i] = new double[ncols];
      Rrelohihi_h[i] = new double[ncols];
      Rrehilohi_h[i] = new double[ncols];
      Rrelolohi_h[i] = new double[ncols];
      Rrehihilo_h[i] = new double[ncols];
      Rrelohilo_h[i] = new double[ncols];
      Rrehilolo_h[i] = new double[ncols];
      Rrelololo_h[i] = new double[ncols];
      Rimhihihi_h[i] = new double[ncols];
      Rimlohihi_h[i] = new double[ncols];
      Rimhilohi_h[i] = new double[ncols];
      Rimlolohi_h[i] = new double[ncols];
      Rimhihilo_h[i] = new double[ncols];
      Rimlohilo_h[i] = new double[ncols];
      Rimhilolo_h[i] = new double[ncols];
      Rimlololo_h[i] = new double[ncols];
      Rrehihihi_d[i] = new double[ncols];
      Rrelohihi_d[i] = new double[ncols];
      Rrehilohi_d[i] = new double[ncols];
      Rrelolohi_d[i] = new double[ncols];
      Rrehihilo_d[i] = new double[ncols];
      Rrelohilo_d[i] = new double[ncols];
      Rrehilolo_d[i] = new double[ncols];
      Rrelololo_d[i] = new double[ncols];
      Rimhihihi_d[i] = new double[ncols];
      Rimlohihi_d[i] = new double[ncols];
      Rimhilohi_d[i] = new double[ncols];
      Rimlolohi_d[i] = new double[ncols];
      Rimhihilo_d[i] = new double[ncols];
      Rimlohilo_d[i] = new double[ncols];
      Rimhilolo_d[i] = new double[ncols];
      Rimlololo_d[i] = new double[ncols];
   }
   random_cmplx8_matrix
      (nrows,ncols,Arehihihi,Arelohihi,Arehilohi,Arelolohi,
                   Arehihilo,Arelohilo,Arehilolo,Arelololo,
                   Aimhihihi,Aimlohihi,Aimhilohi,Aimlolohi,
                   Aimhihilo,Aimlohilo,Aimhilolo,Aimlololo);

   if(verbose > 0)
   {
      cout << scientific << setprecision(16);

      cout << "A random matrix :" << endl;
      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++)
         {
            cout << "A[" << i << "][" << j << "]re : "
                 << Arehihihi[i][j] << "  " << Arelohihi[i][j] << endl
                 << "            "
                 << Arehilohi[i][j] << "  " << Arelolohi[i][j] << endl
                 << "            "
                 << Arehihilo[i][j] << "  " << Arelohilo[i][j] << endl
                 << "            "
                 << Arehilolo[i][j] << "  " << Arelololo[i][j] << endl;
            cout << "A[" << i << "][" << j << "]im : "
                 << Aimhihihi[i][j] << "  " << Aimlohihi[i][j] << endl
                 << "            "
                 << Aimhilohi[i][j] << "  " << Aimlolohi[i][j] << endl
                 << "            "
                 << Aimhihilo[i][j] << "  " << Aimlohilo[i][j] << endl
                 << "            "
                 << Aimhilolo[i][j] << "  " << Aimlololo[i][j] << endl;
         }
   }
   double timelapsed_h;
   bool bvrb = (verbose > 0);
   const double tol = 1.0e-40;
   int fail;

   if((mode == 1) || (mode == 2))
   {
      cout << "-> CPU computes the blocked Householder QR ..." << endl;

      CPU_cmplx8_blocked_houseqr
         (nrows,ncols,sizetile,numtiles,
          Arehihihi,  Arelohihi,  Arehilohi,  Arelolohi,
          Arehihilo,  Arelohilo,  Arehilolo,  Arelololo,
          Aimhihihi,  Aimlohihi,  Aimhilohi,  Aimlolohi,
          Aimhihilo,  Aimlohilo,  Aimhilolo,  Aimlololo,
          Qrehihihi_h,Qrelohihi_h,Qrehilohi_h,Qrelolohi_h,
          Qrehihilo_h,Qrelohilo_h,Qrehilolo_h,Qrelololo_h,
          Qimhihihi_h,Qimlohihi_h,Qimhilohi_h,Qimlolohi_h,
          Qimhihilo_h,Qimlohilo_h,Qimhilolo_h,Qimlololo_h,
          Rrehihihi_h,Rrelohihi_h,Rrehilohi_h,Rrelolohi_h,
          Rrehihilo_h,Rrelohilo_h,Rrehilolo_h,Rrelololo_h,
          Rimhihihi_h,Rimlohihi_h,Rimhilohi_h,Rimlolohi_h,
          Rimhihilo_h,Rimlohilo_h,Rimhilolo_h,Rimlololo_h,&timelapsed_h,bvrb);

      cout << "-> Testing the QR factorization ..." << endl;

      fail = test_cmplx8_qr_factors_probe
         (nrows,ncols,
          Arehihihi,  Arelohihi,  Arehilohi,  Arelolohi,
          Arehihilo,  Arelohilo,  Arehilolo,  Arelololo,
          Aimhihihi,  Aimlohihi,  Aimhilohi,  Aimlolohi,
          Aimhihilo,  Aimlohilo,  Aimhilolo,  Aimlololo,
          Qrehihihi_h,Qrelohihi_h,Qrehilohi_h,Qrelolohi_h,
          Qrehihilo_h,Qrelohilo_h,Qrehilolo_h,Qrelololo_h,
          Qimhihihi_h,Qimlohihi_h,Qimhilohi_h,Qimlolohi_h,
          Qimhihilo_h,Qimlohilo_h,Qimhilolo_h,Qimlololo_h,
          Rrehihihi_h,Rrelohihi_h,Rrehilohi_h,Rrelolohi_h,
          Rrehihilo_h,Rrelohilo_h,Rrehilolo_h,Rrelololo_h,
          Rimhihihi_h,Rimlohihi_h,Rimhilohi_h,Rimlolohi_h,
          Rimhihilo_h,Rimlohilo_h,Rimhilolo_h,Rimlololo_h,tol,2,true);

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
                    << Arehihihi[i][j] << "  " << Arelohihi[i][j] << endl
                    << "            "
                    << Arehilohi[i][j] << "  " << Arelolohi[i][j] << endl
                    << "            "
                    << Arehihilo[i][j] << "  " << Arelohilo[i][j] << endl
                    << "            "
                    << Arehilolo[i][j] << "  " << Arelololo[i][j] << endl;
               cout << "A[" << i << "][" << j << "]im : "
                    << Aimhihihi[i][j] << "  " << Aimlohihi[i][j] << endl
                    << "            "
                    << Aimhilohi[i][j] << "  " << Aimlolohi[i][j] << endl
                    << "            "
                    << Aimhihilo[i][j] << "  " << Aimlohilo[i][j] << endl
                    << "            "
                    << Aimhilolo[i][j] << "  " << Aimlololo[i][j] << endl;
            }
      }
      GPU_cmplx8_blocked_houseqr
         (nrows,ncols,sizetile,numtiles,
          Arehihihi,  Arelohihi,  Arehilohi,  Arelolohi,
          Arehihilo,  Arelohilo,  Arehilolo,  Arelololo,
          Aimhihihi,  Aimlohihi,  Aimhilohi,  Aimlolohi,
          Aimhihilo,  Aimlohilo,  Aimhilolo,  Aimlololo,
          Qrehihihi_d,Qrelohihi_d,Qrehilohi_d,Qrelolohi_d,
          Qrehihilo_d,Qrelohilo_d,Qrehilolo_d,Qrelololo_d,
          Qimhihihi_d,Qimlohihi_d,Qimhilohi_d,Qimlolohi_d,
          Qimhihilo_d,Qimlohilo_d,Qimhilolo_d,Qimlololo_d,
          Rrehihihi_d,Rrelohihi_d,Rrehilohi_d,Rrelolohi_d,
          Rrehihilo_d,Rrelohilo_d,Rrehilolo_d,Rrelololo_d,
          Rimhihihi_d,Rimlohihi_d,Rimhilohi_d,Rimlolohi_d,
          Rimhihilo_d,Rimlohilo_d,Rimhilolo_d,Rimlololo_d,
          &houselapsedms,&RHvlapsedms,&tileRlapsedms,&vb2Wlapsedms,
          &WYTlapsedms,&QWYTlapsedms,&Qaddlapsedms,
          &YWTlapsedms,&YWTClapsedms,&Raddlapsedms,&timelapsed_d,
          &addcnt,&mulcnt,&divcnt,&sqrtcnt,bvrb);

      cout << "-> Testing the QR factorization ..." << endl;

      fail = test_cmplx8_qr_factors_probe
                (nrows,ncols,Arehihihi,  Arelohihi,  Arehilohi,  Arelolohi,
                             Arehihilo,  Arelohilo,  Arehilolo,  Arelololo,
                             Aimhihihi,  Aimlohihi,  Aimhilohi,  Aimlolohi,
                             Aimhihilo,  Aimlohilo,  Aimhilolo,  Aimlololo,
                             Qrehihihi_d,Qrelohihi_d,Qrehilohi_d,Qrelolohi_d,
                             Qrehihilo_d,Qrelohilo_d,Qrehilolo_d,Qrelololo_d,
                             Qimhihihi_d,Qimlohihi_d,Qimhilohi_d,Qimlolohi_d,
                             Qimhihilo_d,Qimlohilo_d,Qimhilolo_d,Qimlololo_d,
                             Rrehihihi_d,Rrelohihi_d,Rrehilohi_d,Rrelolohi_d,
                             Rrehihilo_d,Rrelohilo_d,Rrehilolo_d,Rrelololo_d,
                             Rimhihihi_d,Rimlohihi_d,Rimhilohi_d,Rimlolohi_d,
                             Rimhihilo_d,Rimlohilo_d,Rimhilolo_d,Rimlololo_d,
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
      write_dbl8_qrtimeflops
         (1,nrows,ncols,houselapsedms,RHvlapsedms,tileRlapsedms,vb2Wlapsedms,
          WYTlapsedms,QWYTlapsedms,Qaddlapsedms,YWTlapsedms,YWTClapsedms,
          Raddlapsedms,timelapsed_d,addcnt,mulcnt,divcnt,sqrtcnt);
   }
   for(int i=0; i<nrows; i++)
   {
      free(Arehihihi[i]); free(Arelohihi[i]);
      free(Arehilohi[i]); free(Arelolohi[i]);
      free(Arehihilo[i]); free(Arelohilo[i]);
      free(Arehilolo[i]); free(Arelololo[i]);
      free(Aimhihihi[i]); free(Aimlohihi[i]);
      free(Aimhilohi[i]); free(Aimlolohi[i]);
      free(Aimhihilo[i]); free(Aimlohilo[i]);
      free(Aimhilolo[i]); free(Aimlololo[i]);
      free(Qrehihihi_h[i]); free(Qrelohihi_h[i]);
      free(Qrehilohi_h[i]); free(Qrelolohi_h[i]);
      free(Qrehihilo_h[i]); free(Qrelohilo_h[i]);
      free(Qrehilolo_h[i]); free(Qrelololo_h[i]);
      free(Qrehihihi_d[i]); free(Qrelohihi_d[i]);
      free(Qrehilohi_d[i]); free(Qrelolohi_d[i]);
      free(Qrehihilo_d[i]); free(Qrelohilo_d[i]);
      free(Qrehilolo_d[i]); free(Qrelololo_d[i]);
      free(Qimhihihi_h[i]); free(Qimlohihi_h[i]);
      free(Qimhilohi_h[i]); free(Qimlolohi_h[i]);
      free(Qimhihilo_h[i]); free(Qimlohilo_h[i]);
      free(Qimhilolo_h[i]); free(Qimlololo_h[i]);
      free(Qimhihihi_d[i]); free(Qimlohihi_d[i]);
      free(Qimhilohi_d[i]); free(Qimlolohi_d[i]);
      free(Qimhihilo_d[i]); free(Qimlohilo_d[i]);
      free(Qimhilolo_d[i]); free(Qimlololo_d[i]);
      free(Rrehihihi_h[i]); free(Rrelohihi_h[i]);
      free(Rrehilohi_h[i]); free(Rrelolohi_h[i]);
      free(Rrehihilo_h[i]); free(Rrelohilo_h[i]);
      free(Rrehilolo_h[i]); free(Rrelololo_h[i]);
      free(Rrehihihi_d[i]); free(Rrelohihi_d[i]);
      free(Rrehilohi_d[i]); free(Rrelolohi_d[i]);
      free(Rrehihilo_d[i]); free(Rrelohilo_d[i]);
      free(Rrehilolo_d[i]); free(Rrelololo_d[i]);
      free(Rimhihihi_h[i]); free(Rimlohihi_h[i]);
      free(Rimhilohi_h[i]); free(Rimlolohi_h[i]);
      free(Rimhihilo_h[i]); free(Rimlohilo_h[i]);
      free(Rimhilolo_h[i]); free(Rimlololo_h[i]);
      free(Rimhihihi_d[i]); free(Rimlohihi_d[i]);
      free(Rimhilohi_d[i]); free(Rimlolohi_d[i]);
      free(Rimhihilo_d[i]); free(Rimlohilo_d[i]);
      free(Rimhilolo_d[i]); free(Rimlololo_d[i]);
   }
   free(Arehihihi);   free(Arelohihi);
   free(Arehilohi);   free(Arelolohi);
   free(Arehihilo);   free(Arelohilo);
   free(Arehilolo);   free(Arelololo);
   free(Aimhihihi);   free(Aimlohihi);
   free(Aimhilohi);   free(Aimlolohi);
   free(Aimhihilo);   free(Aimlohilo);
   free(Aimhilolo);   free(Aimlololo);
   free(Qrehihihi_h); free(Qrelohihi_h);
   free(Qrehilohi_h); free(Qrelolohi_h);
   free(Qrehihilo_h); free(Qrelohilo_h);
   free(Qrehilolo_h); free(Qrelololo_h);
   free(Qrehihihi_d); free(Qrelohihi_d);
   free(Qrehilohi_d); free(Qrelolohi_d);
   free(Qrehihilo_d); free(Qrelohilo_d);
   free(Qrehilolo_d); free(Qrelololo_d);
   free(Rrehihihi_h); free(Rrelohihi_h);
   free(Rrehilohi_h); free(Rrelolohi_h);
   free(Rrehihilo_h); free(Rrelohilo_h);
   free(Rrehilolo_h); free(Rrelololo_h);
   free(Rrehihihi_d); free(Rrelohihi_d);
   free(Rrehilohi_d); free(Rrelolohi_d);
   free(Rrehihilo_d); free(Rrelohilo_d);
   free(Rrehilolo_d); free(Rrelololo_d);
   free(Qimhihihi_h); free(Qimlohihi_h);
   free(Qimhilohi_h); free(Qimlolohi_h);
   free(Qimhihilo_h); free(Qimlohilo_h);
   free(Qimhilolo_h); free(Qimlololo_h);
   free(Qimhihihi_d); free(Qimlohihi_d);
   free(Qimhilohi_d); free(Qimlolohi_d);
   free(Qimhihilo_d); free(Qimlohilo_d);
   free(Qimhilolo_d); free(Qimlololo_d);
   free(Rimhihihi_h); free(Rimlohihi_h);
   free(Rimhilohi_h); free(Rimlolohi_h);
   free(Rimhihilo_h); free(Rimlohilo_h);
   free(Rimhilolo_h); free(Rimlololo_h);
   free(Rimhihihi_d); free(Rimlohihi_d);
   free(Rimhilohi_d); free(Rimlolohi_d);
   free(Rimhihilo_d); free(Rimlohilo_d);
   free(Rimhilolo_d); free(Rimlololo_d);
}
