// The file dbl8_baqr_testers.cpp defines the function with prototypes in
// the file dbl8_baqr_testers.h.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#include "random8_matrices.h"
#include "octo_double_functions.h"
#include "dbl8_factorizations.h"
#include "dbl8_factors_testers.h"
#include "dbl8_baqr_host.h"
#include "dbl8_baqr_kernels.h"
#include "dbl8_tabs_host.h"
#include "dbl8_tabs_kernels.h"
#include "dbl8_tabs_testers.h"
#include "dbl8_test_utilities.h"

using namespace std;

void test_real8_blocked_qrbs
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
   random_dbl8_matrix
      (nrows,ncols,Ahihihi,Alohihi,Ahilohi,Alolohi,
                   Ahihilo,Alohilo,Ahilolo,Alololo);

   if(verbose > 0)
   {
      cout << scientific << setprecision(16);

      cout << "A random matrix :" << endl;
      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++)
            cout << "A[" << i << "][" << j << "] : "
                 << Ahihihi[i][j] << "  " << Alohihi[i][j] << endl
                 << Ahilohi[i][j] << "  " << Alolohi[i][j] << endl
                 << Ahihilo[i][j] << "  " << Alohilo[i][j] << endl
                 << Ahilolo[i][j] << "  " << Alololo[i][j] << endl;
   }
   double *solhihihi = new double[ncols];
   double *sollohihi = new double[ncols];
   double *solhilohi = new double[ncols];
   double *sollolohi = new double[ncols];
   double *solhihilo = new double[ncols];
   double *sollohilo = new double[ncols];
   double *solhilolo = new double[ncols];
   double *sollololo = new double[ncols];

   for(int i=0; i<ncols; i++)
   {
      solhihihi[i] = 1.0;
      sollohihi[i] = 0.0;
      solhilohi[i] = 0.0;
      sollolohi[i] = 0.0;
      solhihilo[i] = 0.0;
      sollohilo[i] = 0.0;
      solhilolo[i] = 0.0;
      sollololo[i] = 0.0;
   }
   double *rhshihihi = new double[nrows];
   double *rhslohihi = new double[nrows];
   double *rhshilohi = new double[nrows];
   double *rhslolohi = new double[nrows];
   double *rhshihilo = new double[nrows];
   double *rhslohilo = new double[nrows];
   double *rhshilolo = new double[nrows];
   double *rhslololo = new double[nrows];
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;

   for(int i=0; i<nrows; i++)
   {
      rhshihihi[i] = 0.0;
      rhslohihi[i] = 0.0;
      rhshilohi[i] = 0.0;
      rhslolohi[i] = 0.0;
      rhshihilo[i] = 0.0;
      rhslohilo[i] = 0.0;
      rhshilolo[i] = 0.0;
      rhslololo[i] = 0.0;

      for(int j=0; j<ncols; j++) // rhs[i] = rhs[i] + A[i][j]*sol[j];
      {
         odf_mul(   Ahihihi[i][j],Alohihi[i][j],Ahilohi[i][j],Alolohi[i][j],
                    Ahihilo[i][j],Alohilo[i][j],Ahilolo[i][j],Alololo[i][j],
                  solhihihi[j], sollohihi[j], solhilohi[j], sollolohi[j],
                  solhihilo[j], sollohilo[j], solhilolo[j], sollololo[j],
                 &acchihihi,   &acclohihi,   &acchilohi,   &acclolohi,
                 &acchihilo,   &acclohilo,   &acchilolo,   &acclololo);
         odf_inc(&rhshihihi[i],&rhslohihi[i],&rhshilohi[i],&rhslolohi[i],
                 &rhshihilo[i],&rhslohilo[i],&rhshilolo[i],&rhslololo[i],
                  acchihihi,    acclohihi,    acchilohi,    acclolohi,
                  acchihilo,    acclohilo,    acchilolo,    acclololo);
      }
   }
   if(verbose > 0)
   {
      cout << "The sums of the columns :" << endl;
      for(int i=0; i<nrows; i++)
         cout << "b[" << i << "] : "
              << rhshihihi[i] << "  " << rhslohihi[i] << endl
              << rhshilolo[i] << "  " << rhslololo[i] << endl
              << rhshihihi[i] << "  " << rhslohihi[i] << endl
              << rhshilolo[i] << "  " << rhslololo[i] << endl;
   }
   double qrtimelapsed_h,bstimelapsed_h;
   bool bvrb = (verbose > 0);
   const double tol = 1.0e-20;
   int fail;

   if((mode == 1) || (mode == 2))
   {
      double *xhihihi = new double[ncols];
      double *xlohihi = new double[ncols];
      double *xhilohi = new double[ncols];
      double *xlolohi = new double[ncols];
      double *xhihilo = new double[ncols];
      double *xlohilo = new double[ncols];
      double *xhilolo = new double[ncols];
      double *xlololo = new double[ncols];
      double *qTrhshihihi = new double[nrows];
      double *qTrhslohihi = new double[nrows];
      double *qTrhshilohi = new double[nrows];
      double *qTrhslolohi = new double[nrows];
      double *qTrhshihilo = new double[nrows];
      double *qTrhslohilo = new double[nrows];
      double *qTrhshilolo = new double[nrows];
      double *qTrhslololo = new double[nrows];

      cout << "-> CPU computes the blocked Householder QR ..." << endl;

      CPU_dbl8_blocked_houseqr
         (nrows,ncols,sizetile,numtiles,
          Ahihihi,  Alohihi,  Ahilohi,  Alolohi,
          Ahihilo,  Alohilo,  Ahilolo,  Alololo,
          Qhihihi_h,Qlohihi_h,Qhilohi_h,Qlolohi_h,
          Qhihilo_h,Qlohilo_h,Qhilolo_h,Qlololo_h,
          Rhihihi_h,Rlohihi_h,Rhilohi_h,Rlolohi_h,
          Rhihilo_h,Rlohilo_h,Rhilolo_h,Rlololo_h,&qrtimelapsed_h,bvrb);

      cout << "-> Testing the QR factorization ..." << endl;

      // fail = test_real2_qr_factors
      //   (nrows,ncols,Ahi,Alo,Qhi_h,Qlo_h,Rhi_h,Rlo_h,tol,verbose);
      fail = test_real8_qr_factors_probe
         (nrows,ncols,Ahihihi,  Alohihi,  Ahilohi,  Alolohi,
                      Ahihilo,  Alohilo,  Ahilolo,  Alololo,
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
      cout << "-> CPU multiplies Q^T with b ..." << endl;

      for(int i=0; i<nrows; i++)
      {
         qTrhshihihi[i] = 0.0;
         qTrhslohihi[i] = 0.0;
         qTrhshilohi[i] = 0.0;
         qTrhslolohi[i] = 0.0;
         qTrhshihilo[i] = 0.0;
         qTrhslohilo[i] = 0.0;
         qTrhshilolo[i] = 0.0;
         qTrhslololo[i] = 0.0;

         for(int j=0; j<nrows; j++) // qTrhs[i] = qTrhs[i] + Q[j][i]*rhs[j];
         {
            odf_mul(Qhihihi_h[j][i],Qlohihi_h[j][i],
                    Qhilohi_h[j][i],Qlolohi_h[j][i],
                    Qhihilo_h[j][i],Qlohilo_h[j][i],
                    Qhilolo_h[j][i],Qlololo_h[j][i],
                  rhshihihi[j],   rhslohihi[j], rhshilohi[j], rhslolohi[j],
                  rhshihilo[j],   rhslohilo[j], rhshilolo[j], rhslololo[j],
                 &acchihihi,     &acclohihi,   &acchilohi,   &acclolohi,
                 &acchihilo,     &acclohilo,   &acchilolo,   &acclololo);
            odf_inc(&qTrhshihihi[i],&qTrhslohihi[i],
                    &qTrhshilohi[i],&qTrhslolohi[i],
                    &qTrhshihilo[i],&qTrhslohilo[i],
                    &qTrhshilolo[i],&qTrhslololo[i],
                       acchihihi,      acclohihi,    acchilohi,    acclolohi,
                       acchihilo,      acclohilo,    acchilolo,    acclololo);
         }
      }
      cout << "-> CPU solves an upper triangular system ..." << endl;

      CPU_dbl8_upper_tiled_solver
         (ncols,sizetile,numtiles,
              Rhihihi_h,  Rlohihi_h,  Rhilohi_h,  Rlolohi_h,
              Rhihilo_h,  Rlohilo_h,  Rhilolo_h,  Rlololo_h,
          qTrhshihihi,qTrhslohihi,qTrhshilohi,qTrhslolohi,
          qTrhshihilo,qTrhslohilo,qTrhshilolo,qTrhslololo,
              xhihihi,    xlohihi,    xhilohi,    xlolohi,
              xhihilo,    xlohilo,    xhilolo,    xlololo,&bstimelapsed_h);

      if(verbose > 0)
      {
         cout << "CPU solution computed with tiling :" << endl;
         cout << scientific << setprecision(16);
         for(int i=0; i<ncols; i++)
            cout << "x[" << i << "] : "
                 << xhihihi[i] << "  " << xlohihi[i] << endl
                 << xhilohi[i] << "  " << xlolohi[i] << endl
                 << xhihilo[i] << "  " << xlohilo[i] << endl
                 << xhilolo[i] << "  " << xlololo[i] << endl;
      }
      cout << scientific << setprecision(2);
      cout << "   Sum of CPU errors on solution : "
           << dbl8_Difference_Sum
                 (ncols,solhihihi,sollohihi,solhilohi,sollolohi,
                        solhihilo,sollohilo,solhilolo,sollololo,
                          xhihihi,  xlohihi,  xhilohi,  xlolohi,
                          xhihilo,  xlohilo,  xhilolo,  xlololo) << endl;

      free(xhihihi); free(xlohihi); free(xhilohi); free(xlolohi);
      free(xhihilo); free(xlohilo); free(xhilolo); free(xlololo);
      free(qTrhshihihi); free(qTrhslohihi);
      free(qTrhshilohi); free(qTrhslolohi);
      free(qTrhshihilo); free(qTrhslohilo);
      free(qTrhshilolo); free(qTrhslololo);
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
      double *xhihihi_d = new double[ncols];
      double *xlohihi_d = new double[ncols];
      double *xhilohi_d = new double[ncols];
      double *xlolohi_d = new double[ncols];
      double *xhihilo_d = new double[ncols];
      double *xlohilo_d = new double[ncols];
      double *xhilolo_d = new double[ncols];
      double *xlololo_d = new double[ncols];
      double *qTrhshihihi_d = new double[nrows];
      double *qTrhslohihi_d = new double[nrows];
      double *qTrhshilohi_d = new double[nrows];
      double *qTrhslolohi_d = new double[nrows];
      double *qTrhshihilo_d = new double[nrows];
      double *qTrhslohilo_d = new double[nrows];
      double *qTrhshilolo_d = new double[nrows];
      double *qTrhslololo_d = new double[nrows];

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
          &YWTlapsedms,&YWTClapsedms,&Raddlapsedms,&qrtimelapsed_d,
          &qraddcnt,&qrmulcnt,&qrdivcnt,&sqrtcnt,bvrb);

      cout << "-> Testing the QR factorization ..." << endl;

      // fail = test_real2_qr_factors
      //          (nrows,ncols,Ahi,Alo,Qhi_d,Qlo_d,Rhi_d,Rlo_d,tol,verbose);
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
      // preliminary CPU computation, for testing purposes only
      cout << "-> CPU multiplies Q^T with b ..." << endl;

      for(int i=0; i<nrows; i++)
      {
         qTrhshihihi_d[i] = 0.0;
         qTrhslohihi_d[i] = 0.0;
         qTrhshilohi_d[i] = 0.0;
         qTrhslolohi_d[i] = 0.0;
         qTrhshihilo_d[i] = 0.0;
         qTrhslohilo_d[i] = 0.0;
         qTrhshilolo_d[i] = 0.0;
         qTrhslololo_d[i] = 0.0;

         for(int j=0; j<ncols; j++) // qTrhs[i] = qTrhs[i] + Q[j][i]*rhs[j];
         {
            odf_mul(Qhihihi_d[j][i],Qlohihi_d[j][i],
                    Qhilohi_d[j][i],Qlolohi_d[j][i],
                    Qhihilo_d[j][i],Qlohilo_d[j][i],
                    Qhilolo_d[j][i],Qlololo_d[j][i],
                  rhshihihi[j],   rhslohihi[j], rhshilohi[j], rhslolohi[j],
                  rhshihilo[j],   rhslohilo[j], rhshilolo[j], rhslololo[j],
                 &acchihihi,     &acclohihi,   &acchilohi,   &acclolohi,
                 &acchihilo,     &acclohilo,   &acchilolo,   &acclololo);
            odf_inc(&qTrhshihihi_d[i],&qTrhslohihi_d[i],
                    &qTrhshilohi_d[i],&qTrhslolohi_d[i],
                    &qTrhshihilo_d[i],&qTrhslohilo_d[i],
                    &qTrhshilolo_d[i],&qTrhslololo_d[i],
                    acchihihi,acclohihi,acchilohi,acclolohi,
                    acchihilo,acclohilo,acchilolo,acclololo);
         }
      }
      cout << "-> GPU solves an upper triangular system ..." << endl;

      GPU_dbl8_upper_tiled_solver
         (ncols,sizetile,numtiles,
              Rhihihi_d,    Rlohihi_d,    Rhilohi_d,    Rlolohi_d,
              Rhihilo_d,    Rlohilo_d,    Rhilolo_d,    Rlololo_d,
          qTrhshihihi_d,qTrhslohihi_d,qTrhshilohi_d,qTrhslolohi_d,
          qTrhshihilo_d,qTrhslohilo_d,qTrhshilolo_d,qTrhslololo_d,
              xhihihi_d,    xlohihi_d,    xhilohi_d,    xlolohi_d,
              xhihilo_d,    xlohilo_d,    xhilolo_d,    xlololo_d,
          &invlapsed,&mullapsed,&sublapsed,&bselapsedms,&bstimelapsed_d,
          &bsaddcnt,&bsmulcnt,&bsdivcnt);

      if(verbose > 0)
      {
         cout << scientific << setprecision(16);
         cout << "GPU solution computed with tiling :" << endl;
         for(int i=0; i<ncols; i++)
            cout << "x[" << i << "] : "
                 << xhihihi_d[i] << "  " << xlohihi_d[i] << endl
                 << xhilohi_d[i] << "  " << xlolohi_d[i] << endl
                 << xhihilo_d[i] << "  " << xlohilo_d[i] << endl
                 << xhilolo_d[i] << "  " << xlololo_d[i] << endl;
      }
      cout << scientific << setprecision(2);
      cout << "   Sum of GPU errors on solution : "
           << dbl8_Difference_Sum
                 (ncols,solhihihi,sollohihi,solhilohi,sollolohi,
                        solhihilo,sollohilo,solhilolo,sollololo,
                          xhihihi_d,xlohihi_d,xhilohi_d,xlolohi_d,
                          xhihilo_d,xlohilo_d,xhilolo_d,xlololo_d) << endl;

      free(xhihihi_d); free(xlohihi_d); free(xhilohi_d); free(xlolohi_d);
      free(xhihilo_d); free(xlohilo_d); free(xhilolo_d); free(xlololo_d);
      free(qTrhshihihi_d); free(qTrhslohihi_d);
      free(qTrhshilohi_d); free(qTrhslolohi_d);
      free(qTrhshihilo_d); free(qTrhslohilo_d);
      free(qTrhshilolo_d); free(qTrhslololo_d);
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
           << qraddcnt << " x 270" << endl;
      cout << "                    Number of multiplications : "
           << qrmulcnt << " x 1742" << endl;
      cout << "                          Number of divisions : "
           << qrdivcnt << " x 5126" << endl;
      cout << "                    Number of calls to sqrt() : "
           << sqrtcnt << " x 8491" << endl;
      long long int qrflopcnt
          = 270*qraddcnt + 1742*qrmulcnt + 5126*qrdivcnt + 8491*sqrtcnt;
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
           << bsaddcnt << " x 270" << endl;
      cout << "                    Number of multiplications : "
           << bsmulcnt << " x 1742" << endl;
      cout << "                          Number of divisions : "
           << bsdivcnt << " x 5126" << endl;
      long long int bsflopcnt = 270*bsaddcnt + 1742*bsmulcnt + 5126*bsdivcnt;
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
      free(Ahihihi[i]); free(Alohihi[i]); free(Ahilohi[i]); free(Alolohi[i]);
      free(Ahihilo[i]); free(Alohilo[i]); free(Ahilolo[i]); free(Alololo[i]);
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
   free(Ahihihi); free(Alohihi); free(Ahilohi); free(Alolohi);
   free(Ahihilo); free(Alohilo); free(Ahilolo); free(Alololo);
   free(Qhihihi_h); free(Qlohihi_h); free(Qhilohi_h); free(Qlolohi_h);
   free(Qhihilo_h); free(Qlohilo_h); free(Qhilolo_h); free(Qlololo_h);
   free(Qhihihi_d); free(Qlohihi_d); free(Qhilohi_d); free(Qlolohi_d);
   free(Qhihilo_d); free(Qlohilo_d); free(Qhilolo_d); free(Qlololo_d);
   free(Rhihihi_h); free(Rlohihi_h); free(Rhilohi_h); free(Rlolohi_h);
   free(Rhihilo_h); free(Rlohilo_h); free(Rhilolo_h); free(Rlololo_h);
   free(Rhihihi_d); free(Rlohihi_d); free(Rhilohi_d); free(Rlolohi_d);
   free(Rhihilo_d); free(Rlohilo_d); free(Rhilolo_d); free(Rlololo_d);
   free(solhihihi); free(sollohihi); free(solhilohi); free(sollolohi);
   free(solhihilo); free(sollohilo); free(solhilolo); free(sollololo);
   free(rhshihihi); free(rhslohihi); free(rhshilohi); free(rhslolohi);
   free(rhshihilo); free(rhslohilo); free(rhshilolo); free(rhslololo);
}

void test_cmplx8_blocked_qrbs
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
                 << Arehilohi[i][j] << "  " << Arelolohi[i][j] << endl
                 << Arehihilo[i][j] << "  " << Arelohilo[i][j] << endl
                 << Arehilolo[i][j] << "  " << Arelololo[i][j] << endl;
            cout << "A[" << i << "][" << j << "]im : "
                 << Aimhihihi[i][j] << "  " << Aimlohihi[i][j] << endl
                 << Aimhilohi[i][j] << "  " << Aimlolohi[i][j] << endl
                 << Aimhihilo[i][j] << "  " << Aimlohilo[i][j] << endl
                 << Aimhilolo[i][j] << "  " << Aimlololo[i][j] << endl;
         }
   }
   double *solrehihihi = new double[ncols];
   double *solrelohihi = new double[ncols];
   double *solrehilohi = new double[ncols];
   double *solrelolohi = new double[ncols];
   double *solrehihilo = new double[ncols];
   double *solrelohilo = new double[ncols];
   double *solrehilolo = new double[ncols];
   double *solrelololo = new double[ncols];
   double *solimhihihi = new double[ncols];
   double *solimlohihi = new double[ncols];
   double *solimhilohi = new double[ncols];
   double *solimlolohi = new double[ncols];
   double *solimhihilo = new double[ncols];
   double *solimlohilo = new double[ncols];
   double *solimhilolo = new double[ncols];
   double *solimlololo = new double[ncols];

   for(int i=0; i<ncols; i++)
   {
      solrehihihi[i] = 1.0; solrelohihi[i] = 0.0;
      solrehilohi[i] = 0.0; solrelolohi[i] = 0.0;
      solrehihilo[i] = 0.0; solrelohilo[i] = 0.0;
      solrehilolo[i] = 0.0; solrelololo[i] = 0.0;
      solimhihihi[i] = 0.0; solimlohihi[i] = 0.0;
      solimhilohi[i] = 0.0; solimlolohi[i] = 0.0;
      solimhihilo[i] = 0.0; solimlohilo[i] = 0.0;
      solimhilolo[i] = 0.0; solimlololo[i] = 0.0;
   }
   double *rhsrehihihi = new double[nrows];
   double *rhsrelohihi = new double[nrows];
   double *rhsrehilohi = new double[nrows];
   double *rhsrelolohi = new double[nrows];
   double *rhsrehihilo = new double[nrows];
   double *rhsrelohilo = new double[nrows];
   double *rhsrehilolo = new double[nrows];
   double *rhsrelololo = new double[nrows];
   double *rhsimhihihi = new double[nrows];
   double *rhsimlohihi = new double[nrows];
   double *rhsimhilohi = new double[nrows];
   double *rhsimlolohi = new double[nrows];
   double *rhsimhihilo = new double[nrows];
   double *rhsimlohilo = new double[nrows];
   double *rhsimhilolo = new double[nrows];
   double *rhsimlololo = new double[nrows];
   double acc1hihihi,acc1lohihi,acc1hilohi,acc1lolohi;
   double acc1hihilo,acc1lohilo,acc1hilolo,acc1lololo;
   double acc2hihihi,acc2lohihi,acc2hilohi,acc2lolohi;
   double acc2hihilo,acc2lohilo,acc2hilolo,acc2lololo;
   double acc3hihihi,acc3lohihi,acc3hilohi,acc3lolohi;
   double acc3hihilo,acc3lohilo,acc3hilolo,acc3lololo;
   double acc4hihihi,acc4lohihi,acc4hilohi,acc4lolohi;
   double acc4hihilo,acc4lohilo,acc4hilolo,acc4lololo;

   for(int i=0; i<nrows; i++)
   {
      rhsrehihihi[i] = 0.0; rhsrelohihi[i] = 0.0;
      rhsrehilohi[i] = 0.0; rhsrelolohi[i] = 0.0;
      rhsrehihilo[i] = 0.0; rhsrelohilo[i] = 0.0;
      rhsrehilolo[i] = 0.0; rhsrelololo[i] = 0.0;
      rhsimhihihi[i] = 0.0; rhsimlohihi[i] = 0.0;
      rhsimhilohi[i] = 0.0; rhsimlolohi[i] = 0.0;
      rhsimhihilo[i] = 0.0; rhsimlohilo[i] = 0.0;
      rhsimhilolo[i] = 0.0; rhsimlololo[i] = 0.0;

      for(int j=0; j<ncols; j++) // rhs[i] = rhs[i] + A[i][j]*sol[j];
      {
         odf_mul(Arehihihi[i][j],Arelohihi[i][j],
                 Arehilohi[i][j],Arelolohi[i][j],
                 Arehihilo[i][j],Arelohilo[i][j],
                 Arehilolo[i][j],Arelololo[i][j],
               solrehihihi[j], solrelohihi[j], solrehilohi[j], solrelolohi[j],
               solrehihilo[j], solrelohilo[j], solrehilolo[j], solrelololo[j],
               &acc1hihihi,    &acc1lohihi,    &acc1hilohi,    &acc1lolohi,
               &acc1hihilo,    &acc1lohilo,    &acc1hilolo,    &acc1lololo);
         odf_mul(Aimhihihi[i][j],Aimlohihi[i][j],
                 Aimhilohi[i][j],Aimlolohi[i][j],
                 Aimhihilo[i][j],Aimlohilo[i][j],
                 Aimhilolo[i][j],Aimlololo[i][j],
               solimhihihi[j], solimlohihi[j], solimhilohi[j], solimlolohi[j],
               solimhihilo[j], solimlohilo[j], solimhilolo[j], solimlololo[j],
               &acc2hihihi,    &acc2lohihi,    &acc2hilohi,    &acc2lolohi,
               &acc2hihilo,    &acc2lohilo,    &acc2hilolo,    &acc2lololo);
         odf_mul(Aimhihihi[i][j],Aimlohihi[i][j],
                 Aimhilohi[i][j],Aimlolohi[i][j],
                 Aimhihilo[i][j],Aimlohilo[i][j],
                 Aimhilolo[i][j],Aimlololo[i][j],
               solrehihihi[j], solrelohihi[j], solrehilohi[j], solrelolohi[j],
               solrehihilo[j], solrelohilo[j], solrehilolo[j], solrelololo[j],
               &acc3hihihi,    &acc3lohihi,    &acc3hilohi,    &acc3lolohi,
               &acc3hihilo,    &acc3lohilo,    &acc3hilolo,    &acc3lololo);
         odf_mul(Arehihihi[i][j],Arelohihi[i][j],
                 Arehilohi[i][j],Arelolohi[i][j],
                 Arehihilo[i][j],Arelohilo[i][j],
                 Arehilolo[i][j],Arelololo[i][j],
               solimhihihi[j], solimlohihi[j], solimhilohi[j], solimlolohi[j],
               solimhihilo[j], solimlohilo[j], solimhilolo[j], solimlololo[j],
               &acc4hihihi,    &acc4lohihi,    &acc4hilohi,    &acc4lolohi,
               &acc4hihilo,    &acc4lohilo,    &acc4hilolo,    &acc4lololo);

         odf_dec(&acc1hihihi,&acc1lohihi,&acc1hilohi,&acc1lolohi,
                 &acc1hihilo,&acc1lohilo,&acc1hilolo,&acc1lololo,
                  acc2hihihi, acc2lohihi, acc2hilohi, acc2lolohi,
                  acc2hihilo, acc2lohilo, acc2hilolo, acc2lololo);
         odf_inc(&rhsrehihihi[i],&rhsrelohihi[i],
                 &rhsrehilohi[i],&rhsrelolohi[i],
                 &rhsrehihilo[i],&rhsrelohilo[i],
                 &rhsrehilolo[i],&rhsrelololo[i],
                   acc1hihihi,     acc1lohihi,   acc1hilohi,   acc1lolohi,
                   acc1hihilo,     acc1lohilo,   acc1hilolo,   acc1lololo);
         odf_inc(&acc3hihihi,&acc3lohihi,&acc3hilohi,&acc3lolohi,
                 &acc3hihilo,&acc3lohilo,&acc3hilolo,&acc3lololo,
                  acc4hihihi, acc4lohihi, acc4hilohi, acc4lolohi,
                  acc4hihilo, acc4lohilo, acc4hilolo, acc4lololo);
         odf_inc(&rhsimhihihi[i],&rhsimlohihi[i],
                 &rhsimhilohi[i],&rhsimlolohi[i],
                 &rhsimhihilo[i],&rhsimlohilo[i],
                 &rhsimhilolo[i],&rhsimlololo[i],
                   acc3hihihi,     acc3lohihi,   acc3hilohi,   acc3lolohi,
                   acc3hihilo,     acc3lohilo,   acc3hilolo,   acc3lololo);
      }
   }
   if(verbose > 0)
   {
      cout << "The sums of the columns :" << endl;
      for(int i=0; i<nrows; i++)
      {
         cout << "b[" << i << "]re : "
              << rhsrehihihi[i] << "  " << rhsrelohihi[i] << endl
              << rhsrehilohi[i] << "  " << rhsrelolohi[i] << endl
              << rhsrehihilo[i] << "  " << rhsrelohilo[i] << endl
              << rhsrehilolo[i] << "  " << rhsrelololo[i] << endl;
         cout << "b[" << i << "]im : "
              << rhsimhihihi[i] << "  " << rhsimlohihi[i] << endl
              << rhsimhilohi[i] << "  " << rhsimlolohi[i] << endl
              << rhsimhihilo[i] << "  " << rhsimlohilo[i] << endl
              << rhsimhilolo[i] << "  " << rhsimlololo[i] << endl;
      }
   }
   double qrtimelapsed_h,bstimelapsed_h;
   bool bvrb = (verbose > 0);
   const double tol = 1.0e-20;
   int fail;

   if((mode == 1) || (mode == 2))
   {
      double *xrehihihi = new double[ncols];
      double *xrelohihi = new double[ncols];
      double *xrehilohi = new double[ncols];
      double *xrelolohi = new double[ncols];
      double *xrehihilo = new double[ncols];
      double *xrelohilo = new double[ncols];
      double *xrehilolo = new double[ncols];
      double *xrelololo = new double[ncols];
      double *ximhihihi = new double[ncols];
      double *ximlohihi = new double[ncols];
      double *ximhilohi = new double[ncols];
      double *ximlolohi = new double[ncols];
      double *ximhihilo = new double[ncols];
      double *ximlohilo = new double[ncols];
      double *ximhilolo = new double[ncols];
      double *ximlololo = new double[ncols];
      double *qHrhsrehihihi = new double[nrows];
      double *qHrhsrelohihi = new double[nrows];
      double *qHrhsrehilohi = new double[nrows];
      double *qHrhsrelolohi = new double[nrows];
      double *qHrhsrehihilo = new double[nrows];
      double *qHrhsrelohilo = new double[nrows];
      double *qHrhsrehilolo = new double[nrows];
      double *qHrhsrelololo = new double[nrows];
      double *qHrhsimhihihi = new double[nrows];
      double *qHrhsimlohihi = new double[nrows];
      double *qHrhsimhilohi = new double[nrows];
      double *qHrhsimlolohi = new double[nrows];
      double *qHrhsimhihilo = new double[nrows];
      double *qHrhsimlohilo = new double[nrows];
      double *qHrhsimhilolo = new double[nrows];
      double *qHrhsimlololo = new double[nrows];

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
          Rimhihilo_h,Rimlohilo_h,Rimhilolo_h,Rimlololo_h,
          &qrtimelapsed_h,bvrb);

      cout << "-> Testing the QR factorization ..." << endl;

      // fail = test_cmplx2_qr_factors
      //    (nrows,ncols,Arehi,  Arelo,  Aimhi,  Aimlo,
      //                 Qrehi_h,Qrelo_h,Qimhi_h,Qimlo_h,
      //                 Rrehi_h,Rrelo_h,Rimhi_h,Rimlo_h,tol,verbose);
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
      cout << "-> CPU multiplies Q^H with b ..." << endl;

      for(int i=0; i<nrows; i++)
      {
         qHrhsrehihihi[i] = 0.0; qHrhsrelohihi[i] = 0.0;
         qHrhsrehilohi[i] = 0.0; qHrhsrelolohi[i] = 0.0;
         qHrhsrehihilo[i] = 0.0; qHrhsrelohilo[i] = 0.0;
         qHrhsrehilolo[i] = 0.0; qHrhsrelololo[i] = 0.0;
         qHrhsimhihihi[i] = 0.0; qHrhsimlohihi[i] = 0.0;
         qHrhsimhilohi[i] = 0.0; qHrhsimlolohi[i] = 0.0;
         qHrhsimhihilo[i] = 0.0; qHrhsimlohilo[i] = 0.0;
         qHrhsimhilolo[i] = 0.0; qHrhsimlololo[i] = 0.0;

         for(int j=0; j<nrows; j++) // qHrhs[i] = qHrhs[i] + Q[j][i]*rhs[j];
         {
            odf_mul(Qrehihihi_h[j][i],Qrelohihi_h[j][i],
                    Qrehilohi_h[j][i],Qrelolohi_h[j][i],
                    Qrehihilo_h[j][i],Qrelohilo_h[j][i],
                    Qrehilolo_h[j][i],Qrelololo_h[j][i],
                    rhsrehihihi[j],rhsrelohihi[j],
                    rhsrehilohi[j],rhsrelolohi[j],
                    rhsrehihilo[j],rhsrelohilo[j],
                    rhsrehilolo[j],rhsrelololo[j],
                    &acc1hihihi,&acc1lohihi,&acc1hilohi,&acc1lolohi,
                    &acc1hihilo,&acc1lohilo,&acc1hilolo,&acc1lololo);
            odf_mul(Qimhihihi_h[j][i],Qimlohihi_h[j][i],
                    Qimhilohi_h[j][i],Qimlolohi_h[j][i],
                    Qimhihilo_h[j][i],Qimlohilo_h[j][i],
                    Qimhilolo_h[j][i],Qimlololo_h[j][i],
                    rhsimhihihi[j],rhsimlohihi[j],
                    rhsimhilohi[j],rhsimlolohi[j],
                    rhsimhihilo[j],rhsimlohilo[j],
                    rhsimhilolo[j],rhsimlololo[j],
                    &acc2hihihi,&acc2lohihi,&acc2hilohi,&acc2lolohi,
                    &acc2hihilo,&acc2lohilo,&acc2hilolo,&acc2lololo);
            odf_mul(Qimhihihi_h[j][i],Qimlohihi_h[j][i],
                    Qimhilohi_h[j][i],Qimlolohi_h[j][i],
                    Qimhihilo_h[j][i],Qimlohilo_h[j][i],
                    Qimhilolo_h[j][i],Qimlololo_h[j][i],
                    rhsrehihihi[j],rhsrelohihi[j],
                    rhsrehilohi[j],rhsrelolohi[j],
                    rhsrehihilo[j],rhsrelohilo[j],
                    rhsrehilolo[j],rhsrelololo[j],
                    &acc3hihihi,&acc3lohihi,&acc3hilohi,&acc3lolohi,
                    &acc3hihilo,&acc3lohilo,&acc3hilolo,&acc3lololo);
            odf_mul(Qrehihihi_h[j][i],Qrelohihi_h[j][i],
                    Qrehilohi_h[j][i],Qrelolohi_h[j][i],
                    Qrehihilo_h[j][i],Qrelohilo_h[j][i],
                    Qrehilolo_h[j][i],Qrelololo_h[j][i],
                    rhsimhihihi[j],rhsimlohihi[j],
                    rhsimhilohi[j],rhsimlolohi[j],
                    rhsimhihilo[j],rhsimlohilo[j],
                    rhsimhilolo[j],rhsimlololo[j],
                    &acc4hihihi,&acc4lohihi,&acc4hilohi,&acc4lolohi,
                    &acc4hihilo,&acc4lohilo,&acc4hilolo,&acc4lololo);

            odf_inc(&acc1hihihi,&acc1lohihi,&acc1hilohi,&acc1lolohi,
                    &acc1hihilo,&acc1lohilo,&acc1hilolo,&acc1lololo,
                     acc2hihihi, acc2lohihi, acc2hilohi, acc2lolohi,
                     acc2hihilo, acc2lohilo, acc2hilolo, acc2lololo);
            odf_inc(&qHrhsrehihihi[i],&qHrhsrelohihi[i],
                    &qHrhsrehilohi[i],&qHrhsrelolohi[i],
                    &qHrhsrehihilo[i],&qHrhsrelohilo[i],
                    &qHrhsrehilolo[i],&qHrhsrelololo[i],
                    acc1hihihi,acc1lohihi,acc1hilohi,acc1lolohi,
                    acc1hihilo,acc1lohilo,acc1hilolo,acc1lololo);
            odf_dec(&acc4hihihi,&acc4lohihi,&acc4hilohi,&acc4lolohi,
                    &acc4hihilo,&acc4lohilo,&acc4hilolo,&acc4lololo,
                     acc3hihihi, acc3lohihi, acc3hilohi, acc3lolohi,
                     acc3hihilo, acc3lohilo, acc3hilolo, acc3lololo);
            odf_inc(&qHrhsimhihihi[i],&qHrhsimlohihi[i],
                    &qHrhsimhilohi[i],&qHrhsimlolohi[i],
                    &qHrhsimhihilo[i],&qHrhsimlohilo[i],
                    &qHrhsimhilolo[i],&qHrhsimlololo[i],
                    acc4hihihi,acc4lohihi,acc4hilohi,acc4lolohi,
                    acc4hihilo,acc4lohilo,acc4hilolo,acc4lololo);
         }
      }
      cout << "-> CPU solves an upper triangular system ..." << endl;

      CPU_cmplx8_upper_tiled_solver
         (ncols,sizetile,numtiles,
              Rrehihihi_h,  Rrelohihi_h,  Rrehilohi_h,  Rrelolohi_h,
              Rrehihilo_h,  Rrelohilo_h,  Rrehilolo_h,  Rrelololo_h,
              Rimhihihi_h,  Rimlohihi_h,  Rimhilohi_h,  Rimlolohi_h,
              Rimhihilo_h,  Rimlohilo_h,  Rimhilolo_h,  Rimlololo_h,
          qHrhsrehihihi,qHrhsrelohihi,qHrhsrehilohi,qHrhsrelolohi,
          qHrhsrehihilo,qHrhsrelohilo,qHrhsrehilolo,qHrhsrelololo,
          qHrhsimhihihi,qHrhsimlohihi,qHrhsimhilohi,qHrhsimlolohi,
          qHrhsimhihilo,qHrhsimlohilo,qHrhsimhilolo,qHrhsimlololo,
              xrehihihi,    xrelohihi,    xrehilohi,    xrelolohi,
              xrehihilo,    xrelohilo,    xrehilolo,    xrelololo,
              ximhihihi,    ximlohihi,    ximhilohi,    ximlolohi,
              ximhihilo,    ximlohilo,    ximhilolo,    ximlololo,
          &bstimelapsed_h);

      if(verbose > 0)
      {
         cout << "CPU solution computed with tiling :" << endl;
         cout << scientific << setprecision(16);

         for(int i=0; i<ncols; i++)
         {
            cout << "x[" << i << "]re : "
                 << xrehihihi[i] << "  " << xrelohihi[i] << endl
                 << xrehilohi[i] << "  " << xrelolohi[i] << endl
                 << xrehihilo[i] << "  " << xrelohilo[i] << endl
                 << xrehilolo[i] << "  " << xrelololo[i] << endl;
            cout << "x[" << i << "]im : "
                 << ximhihihi[i] << "  " << ximlohihi[i] << endl
                 << ximhilohi[i] << "  " << ximlolohi[i] << endl
                 << ximhihilo[i] << "  " << ximlohilo[i] << endl
                 << ximhilolo[i] << "  " << ximlololo[i] << endl;
         }
      }
      cout << scientific << setprecision(2);
      cout << "   Sum of CPU errors on solution : "
           << cmplx8_Difference_Sum
                 (ncols,solrehihihi,solrelohihi,solrehilohi,solrelolohi,
                        solrehihilo,solrelohilo,solrehilolo,solrelololo,
                        solimhihihi,solimlohihi,solimhilohi,solimlolohi,
                        solimhihilo,solimlohilo,solimhilolo,solimlololo,
                          xrehihihi,  xrelohihi,  xrehilohi,  xrelolohi,
                          xrehihilo,  xrelohilo,  xrehilolo,  xrelololo,
                          ximhihihi,  ximlohihi,  ximhilohi,  ximlolohi,
                          ximhihilo,  ximlohilo,  ximhilolo,  ximlololo)
           << endl;

      free(xrehihihi); free(xrelohihi); free(xrehilohi); free(xrelolohi);
      free(xrehihilo); free(xrelohilo); free(xrehilolo); free(xrelololo);
      free(ximhihihi); free(ximlohihi); free(ximhilohi); free(ximlolohi);
      free(ximhihilo); free(ximlohilo); free(ximhilolo); free(ximlololo);
      free(qHrhsrehihihi); free(qHrhsrelohihi);
      free(qHrhsrehilohi); free(qHrhsrelolohi);
      free(qHrhsrehihilo); free(qHrhsrelohilo);
      free(qHrhsrehilolo); free(qHrhsrelololo);
      free(qHrhsimhihihi); free(qHrhsimlohihi);
      free(qHrhsimhilohi); free(qHrhsimlolohi);
      free(qHrhsimhihilo); free(qHrhsimlohilo);
      free(qHrhsimhilolo); free(qHrhsimlololo);
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
      double *xrehihihi_d = new double[ncols];
      double *xrelohihi_d = new double[ncols];
      double *xrehilohi_d = new double[ncols];
      double *xrelolohi_d = new double[ncols];
      double *xrehihilo_d = new double[ncols];
      double *xrelohilo_d = new double[ncols];
      double *xrehilolo_d = new double[ncols];
      double *xrelololo_d = new double[ncols];
      double *ximhihihi_d = new double[ncols];
      double *ximlohihi_d = new double[ncols];
      double *ximhilohi_d = new double[ncols];
      double *ximlolohi_d = new double[ncols];
      double *ximhihilo_d = new double[ncols];
      double *ximlohilo_d = new double[ncols];
      double *ximhilolo_d = new double[ncols];
      double *ximlololo_d = new double[ncols];
      double *qHrhsrehihihi_d = new double[nrows];
      double *qHrhsrelohihi_d = new double[nrows];
      double *qHrhsrehilohi_d = new double[nrows];
      double *qHrhsrelolohi_d = new double[nrows];
      double *qHrhsrehihilo_d = new double[nrows];
      double *qHrhsrelohilo_d = new double[nrows];
      double *qHrhsrehilolo_d = new double[nrows];
      double *qHrhsrelololo_d = new double[nrows];
      double *qHrhsimhihihi_d = new double[nrows];
      double *qHrhsimlohihi_d = new double[nrows];
      double *qHrhsimhilohi_d = new double[nrows];
      double *qHrhsimlolohi_d = new double[nrows];
      double *qHrhsimhihilo_d = new double[nrows];
      double *qHrhsimlohilo_d = new double[nrows];
      double *qHrhsimhilolo_d = new double[nrows];
      double *qHrhsimlololo_d = new double[nrows];

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
                    << Arehilohi[i][j] << "  " << Arelolohi[i][j] << endl
                    << Arehihilo[i][j] << "  " << Arelohilo[i][j] << endl
                    << Arehilolo[i][j] << "  " << Arelololo[i][j] << endl;
               cout << "A[" << i << "][" << j << "]im : "
                    << Aimhihihi[i][j] << "  " << Aimlohihi[i][j] << endl
                    << Aimhilohi[i][j] << "  " << Aimlolohi[i][j] << endl
                    << Aimhihilo[i][j] << "  " << Aimlohilo[i][j] << endl
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
          &YWTlapsedms,&YWTClapsedms,&Raddlapsedms,&qrtimelapsed_d,
          &qraddcnt,&qrmulcnt,&qrdivcnt,&sqrtcnt,bvrb);

      cout << "-> Testing the QR factorization ..." << endl;

      // fail = test_cmplx2_qr_factors
      //           (nrows,ncols,Arehi,  Arelo,  Aimhi,  Aimlo,
      //                        Qrehi_d,Qrelo_d,Qimhi_d,Qimlo_d,
      //                        Rrehi_d,Rrelo_d,Rimhi_d,Rimlo_d,tol,verbose);
      fail = test_cmplx8_qr_factors_probe
                (nrows,ncols,
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
                 Rimhihilo_d,Rimlohilo_d,Rimhilolo_d,Rimlololo_d,tol,2,true);

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
         qHrhsrehihihi_d[i] = 0.0; qHrhsrelohihi_d[i] = 0.0;
         qHrhsrehilohi_d[i] = 0.0; qHrhsrelolohi_d[i] = 0.0;
         qHrhsrehihilo_d[i] = 0.0; qHrhsrelohilo_d[i] = 0.0;
         qHrhsrehilolo_d[i] = 0.0; qHrhsrelololo_d[i] = 0.0;
         qHrhsimhihihi_d[i] = 0.0; qHrhsimlohihi_d[i] = 0.0;
         qHrhsimhilohi_d[i] = 0.0; qHrhsimlolohi_d[i] = 0.0;
         qHrhsimhihilo_d[i] = 0.0; qHrhsimlohilo_d[i] = 0.0;
         qHrhsimhilolo_d[i] = 0.0; qHrhsimlololo_d[i] = 0.0;

         for(int j=0; j<nrows; j++) // qHrhs[i] = qHrhs[i] + Q[j][i]*rhs[j];
         {
            odf_mul(Qrehihihi_d[j][i],Qrelohihi_d[j][i],
                    Qrehilohi_d[j][i],Qrelolohi_d[j][i],
                    Qrehihilo_d[j][i],Qrelohilo_d[j][i],
                    Qrehilolo_d[j][i],Qrelololo_d[j][i],
                    rhsrehihihi[j],rhsrelohihi[j],
                    rhsrehilohi[j],rhsrelolohi[j],
                    rhsrehihilo[j],rhsrelohilo[j],
                    rhsrehilolo[j],rhsrelololo[j],
                    &acc1hihihi,   &acc1lohihi,   &acc1hilohi,   &acc1lolohi,
                    &acc1hihilo,   &acc1lohilo,   &acc1hilolo,   &acc1lololo);
            odf_mul(Qimhihihi_d[j][i],Qimlohihi_d[j][i],
                    Qimhilohi_d[j][i],Qimlolohi_d[j][i],
                    Qimhihilo_d[j][i],Qimlohilo_d[j][i],
                    Qimhilolo_d[j][i],Qimlololo_d[j][i],
                    rhsimhihihi[j],rhsimlohihi[j],
                    rhsimhilohi[j],rhsimlolohi[j],
                    rhsimhihilo[j],rhsimlohilo[j],
                    rhsimhilolo[j],rhsimlololo[j],
                    &acc2hihihi,   &acc2lohihi,   &acc2hilohi,   &acc2lolohi,
                    &acc2hihilo,   &acc2lohilo,   &acc2hilolo,   &acc2lololo);
            odf_mul(Qimhihihi_d[j][i],Qimlohihi_d[j][i],
                    Qimhilohi_d[j][i],Qimlolohi_d[j][i],
                    Qimhihilo_d[j][i],Qimlohilo_d[j][i],
                    Qimhilolo_d[j][i],Qimlololo_d[j][i],
                    rhsrehihihi[j],rhsrelohihi[j],
                    rhsrehilohi[j],rhsrelolohi[j],
                    rhsrehihilo[j],rhsrelohilo[j],
                    rhsrehilolo[j],rhsrelololo[j],
                    &acc3hihihi,   &acc3lohihi,   &acc3hilohi,   &acc3lolohi,
                    &acc3hihilo,   &acc3lohilo,   &acc3hilolo,   &acc3lololo);
            odf_mul(Qrehihihi_d[j][i],Qrelohihi_d[j][i],
                    Qrehilohi_d[j][i],Qrelolohi_d[j][i],
                    Qrehihilo_d[j][i],Qrelohilo_d[j][i],
                    Qrehilolo_d[j][i],Qrelololo_d[j][i],
                    rhsimhihihi[j],rhsimlohihi[j],
                    rhsimhilohi[j],rhsimlolohi[j],
                    rhsimhihilo[j],rhsimlohilo[j],
                    rhsimhilolo[j],rhsimlololo[j],
                    &acc4hihihi,   &acc4lohihi,   &acc4hilohi,   &acc4lolohi,
                    &acc4hihilo,   &acc4lohilo,   &acc4hilolo,   &acc4lololo);

            odf_inc(&acc1hihihi,&acc1lohihi,&acc1hilohi,&acc1lolohi,
                    &acc1hihilo,&acc1lohilo,&acc1hilolo,&acc1lololo,
                     acc2hihihi, acc2lohihi, acc2hilohi, acc2lolohi,
                     acc2hihilo, acc2lohilo, acc2hilolo, acc2lololo);
            odf_inc(&qHrhsrehihihi_d[i],&qHrhsrelohihi_d[i],
                    &qHrhsrehilohi_d[i],&qHrhsrelolohi_d[i],
                    &qHrhsrehihilo_d[i],&qHrhsrelohilo_d[i],
                    &qHrhsrehilolo_d[i],&qHrhsrelololo_d[i],
                    acc1hihihi,acc1lohihi,acc1hilohi,acc1lolohi,
                    acc1hihilo,acc1lohilo,acc1hilolo,acc1lololo);
            odf_dec(&acc4hihihi,&acc4lohihi,&acc4hilohi,&acc4lolohi,
                    &acc4hihilo,&acc4lohilo,&acc4hilolo,&acc4lololo,
                     acc3hihihi, acc3lohihi, acc3hilohi, acc3lolohi,
                     acc3hihilo, acc3lohilo, acc3hilolo, acc3lololo);
            odf_inc(&qHrhsimhihihi_d[i],&qHrhsimlohihi_d[i],
                    &qHrhsimhilohi_d[i],&qHrhsimlolohi_d[i],
                    &qHrhsimhihilo_d[i],&qHrhsimlohilo_d[i],
                    &qHrhsimhilolo_d[i],&qHrhsimlololo_d[i],
                    acc4hihihi,acc4lohihi,acc4hilohi,acc4lolohi,
                    acc4hihilo,acc4lohilo,acc4hilolo,acc4lololo);
         }
      }
      cout << "-> GPU solves an upper triangular system ..." << endl;

      GPU_cmplx8_upper_tiled_solver
         (ncols,sizetile,numtiles,
              Rrehihihi_d,    Rrelohihi_d,    Rrehilohi_d,    Rrelolohi_d,
              Rrehihilo_d,    Rrelohilo_d,    Rrehilolo_d,    Rrelololo_d,
              Rimhihihi_d,    Rimlohihi_d,    Rimhilohi_d,    Rimlolohi_d,
              Rimhihilo_d,    Rimlohilo_d,    Rimhilolo_d,    Rimlololo_d,
          qHrhsrehihihi_d,qHrhsrelohihi_d,qHrhsrehilohi_d,qHrhsrelolohi_d,
          qHrhsrehihilo_d,qHrhsrelohilo_d,qHrhsrehilolo_d,qHrhsrelololo_d,
          qHrhsimhihihi_d,qHrhsimlohihi_d,qHrhsimhilohi_d,qHrhsimlolohi_d,
          qHrhsimhihilo_d,qHrhsimlohilo_d,qHrhsimhilolo_d,qHrhsimlololo_d,
              xrehihihi_d,    xrelohihi_d,    xrehilohi_d,    xrelolohi_d,
              xrehihilo_d,    xrelohilo_d,    xrehilolo_d,    xrelololo_d,
              ximhihihi_d,    ximlohihi_d,    ximhilohi_d,    ximlolohi_d,
              ximhihilo_d,    ximlohilo_d,    ximhilolo_d,    ximlololo_d,
          &invlapsed,&mullapsed,&sublapsed,&bselapsedms,&bstimelapsed_d,
          &bsaddcnt,&bsmulcnt,&bsdivcnt);

      if(verbose > 0)
      {
         cout << scientific << setprecision(16);
         cout << "GPU solution computed with tiling :" << endl;
         for(int i=0; i<ncols; i++)
         {
            cout << "x[" << i << "]re : "
                 << xrehihihi_d[i] << "  " << xrelohihi_d[i] << endl
                 << xrehilohi_d[i] << "  " << xrelolohi_d[i] << endl
                 << xrehihilo_d[i] << "  " << xrelohilo_d[i] << endl
                 << xrehilolo_d[i] << "  " << xrelololo_d[i] << endl;
            cout << "x[" << i << "]im : "
                 << ximhihihi_d[i] << "  " << ximlohihi_d[i] << endl
                 << ximhilohi_d[i] << "  " << ximlolohi_d[i] << endl
                 << ximhihilo_d[i] << "  " << ximlohilo_d[i] << endl
                 << ximhilolo_d[i] << "  " << ximlololo_d[i] << endl;
         }
      }
      cout << scientific << setprecision(2);
      cout << "   Sum of GPU errors on solution : "
           << cmplx8_Difference_Sum
                 (ncols,solrehihihi,solrelohihi,solrehilohi,solrelolohi,
                        solrehihilo,solrelohilo,solrehilolo,solrelololo,
                        solimhihihi,solimlohihi,solimhilohi,solimlolohi,
                        solimhihilo,solimlohilo,solimhilolo,solimlololo,
                          xrehihihi_d,xrelohihi_d,xrehilohi_d,xrelolohi_d,
                          xrehihilo_d,xrelohilo_d,xrehilolo_d,xrelololo_d,
                          ximhihihi_d,ximlohihi_d,ximhilohi_d,ximlolohi_d,
                          ximhihilo_d,ximlohilo_d,ximhilolo_d,ximlololo_d)
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
           << qraddcnt << " x 270 " << endl;
      cout << "                    Number of multiplications : "
           << qrmulcnt << " x 1742 " << endl;
      cout << "                          Number of divisions : "
           << qrdivcnt << " x 5126 " << endl;
      cout << "                    Number of calls to sqrt() : "
           << sqrtcnt << " x 8491 " << endl;
      long long int qrflopcnt
          = 270*qraddcnt + 1742*qrmulcnt + 5126*qrdivcnt + 8491*sqrtcnt;
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
           << bsaddcnt << " x 270 " << endl;
      cout << "                    Number of multiplications : "
           << bsmulcnt << " x 1742 " << endl;
      cout << "                          Number of divisions : "
           << bsdivcnt << " x 5126 " << endl;
      long long int bsflopcnt = 270*bsaddcnt + 1742*bsmulcnt + 5126*bsdivcnt;
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
      free(Arehihihi[i]); free(Arelohihi[i]);
      free(Arehilohi[i]); free(Arelolohi[i]);
      free(Arehihilo[i]); free(Arelohilo[i]);
      free(Arehilolo[i]); free(Arelololo[i]);
      free(Qrehihihi_h[i]); free(Qrelohihi_h[i]);
      free(Qrehilohi_h[i]); free(Qrelolohi_h[i]);
      free(Qrehihilo_h[i]); free(Qrelohilo_h[i]);
      free(Qrehilolo_h[i]); free(Qrelololo_h[i]);
      free(Qrehihihi_d[i]); free(Qrelohihi_d[i]);
      free(Qrehilohi_d[i]); free(Qrelolohi_d[i]); 
      free(Qrehihilo_d[i]); free(Qrelohilo_d[i]);
      free(Qrehilolo_d[i]); free(Qrelololo_d[i]); 
      free(Rrehihihi_h[i]); free(Rrelohihi_h[i]);
      free(Rrehilohi_h[i]); free(Rrelolohi_h[i]);
      free(Rrehihilo_h[i]); free(Rrelohilo_h[i]);
      free(Rrehilolo_h[i]); free(Rrelololo_h[i]);
      free(Rrehihihi_d[i]); free(Rrelohihi_d[i]);
      free(Rrehilohi_d[i]); free(Rrelolohi_d[i]);
      free(Rrehihilo_d[i]); free(Rrelohilo_d[i]);
      free(Rrehilolo_d[i]); free(Rrelololo_d[i]);
      free(Aimhihihi[i]); free(Aimlohihi[i]);
      free(Aimhilohi[i]); free(Aimlolohi[i]);
      free(Aimhihilo[i]); free(Aimlohilo[i]);
      free(Aimhilolo[i]); free(Aimlololo[i]);
      free(Qimhihihi_h[i]); free(Qimlohihi_h[i]);
      free(Qimhilohi_h[i]); free(Qimlolohi_h[i]);
      free(Qimhihilo_h[i]); free(Qimlohilo_h[i]);
      free(Qimhilolo_h[i]); free(Qimlololo_h[i]);
      free(Qimhihihi_d[i]); free(Qimlohihi_d[i]);
      free(Qimhilohi_d[i]); free(Qimlolohi_d[i]);
      free(Qimhihilo_d[i]); free(Qimlohilo_d[i]);
      free(Qimhilolo_d[i]); free(Qimlololo_d[i]);
      free(Rimhihihi_h[i]); free(Rimlohihi_h[i]);
      free(Rimhilohi_h[i]); free(Rimlolohi_h[i]);
      free(Rimhihilo_h[i]); free(Rimlohilo_h[i]);
      free(Rimhilolo_h[i]); free(Rimlololo_h[i]);
      free(Rimhihihi_d[i]); free(Rimlohihi_d[i]);
      free(Rimhilohi_d[i]); free(Rimlolohi_d[i]);
      free(Rimhihilo_d[i]); free(Rimlohilo_d[i]);
      free(Rimhilolo_d[i]); free(Rimlololo_d[i]);
   }
   free(Arehihihi); free(Arelohihi); free(Arehilohi); free(Arelolohi);
   free(Arehihilo); free(Arelohilo); free(Arehilolo); free(Arelololo);
   free(Qrehihihi_h); free(Qrelohihi_h); free(Qrehilohi_h); free(Qrelolohi_h);
   free(Qrehihilo_h); free(Qrelohilo_h); free(Qrehilolo_h); free(Qrelololo_h);
   free(Rrehihihi_h); free(Rrelohihi_h); free(Rrehilohi_h); free(Rrelolohi_h);
   free(Rrehihilo_h); free(Rrelohilo_h); free(Rrehilolo_h); free(Rrelololo_h);
   free(Qrehihihi_d); free(Qrelohihi_d); free(Qrehilohi_d); free(Qrelolohi_d);
   free(Qrehihilo_d); free(Qrelohilo_d); free(Qrehilolo_d); free(Qrelololo_d);
   free(Rrehihihi_d); free(Rrelohihi_d); free(Rrehilohi_d); free(Rrelolohi_d);
   free(Rrehihilo_d); free(Rrelohilo_d); free(Rrehilolo_d); free(Rrelololo_d);
   free(Aimhihihi); free(Aimlohihi); free(Aimhilohi); free(Aimlolohi);
   free(Aimhihilo); free(Aimlohilo); free(Aimhilolo); free(Aimlololo);
   free(Qimhihihi_h); free(Qimlohihi_h); free(Qimhilohi_h); free(Qimlolohi_h);
   free(Qimhihilo_h); free(Qimlohilo_h); free(Qimhilolo_h); free(Qimlololo_h);
   free(Rimhihihi_h); free(Rimlohihi_h); free(Rimhilohi_h); free(Rimlolohi_h);
   free(Rimhihilo_h); free(Rimlohilo_h); free(Rimhilolo_h); free(Rimlololo_h);
   free(Qimhihihi_d); free(Qimlohihi_d); free(Qimhilohi_d); free(Qimlolohi_d);
   free(Qimhihilo_d); free(Qimlohilo_d); free(Qimhilolo_d); free(Qimlololo_d);
   free(Rimhihihi_d); free(Rimlohihi_d); free(Rimhilohi_d); free(Rimlolohi_d);
   free(Rimhihilo_d); free(Rimlohilo_d); free(Rimhilolo_d); free(Rimlololo_d);
   free(solrehihihi); free(solrelohihi); free(solrehilohi); free(solrelolohi);
   free(solrehihilo); free(solrelohilo); free(solrehilolo); free(solrelololo);
   free(solimhihihi); free(solimlohihi); free(solimhilohi); free(solimlolohi);
   free(solimhihilo); free(solimlohilo); free(solimhilolo); free(solimlololo);
   free(rhsrehihihi); free(rhsrelohihi); free(rhsrehilohi); free(rhsrelolohi);
   free(rhsrehihilo); free(rhsrelohilo); free(rhsrehilolo); free(rhsrelololo);
   free(rhsimhihihi); free(rhsimlohihi); free(rhsimhilohi); free(rhsimlolohi);
   free(rhsimhihilo); free(rhsimlohilo); free(rhsimhilolo); free(rhsimlololo);
}
