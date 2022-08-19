// The file dbl_bals_kernels.cu defines the functions with prototypes in
// the file dbl_bals_kernels.h.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#include "dbl_baqr_kernels.h"
#include "dbl_tabs_kernels.h"

using namespace std;

void GPU_dbl_qrbs_solve
 ( int nrows, int ncols, int szt, int nbt,
   double **A, double **Q, double **R, double *b, double *x, bool verbose )
{
   double qrtimelapsed_d;
   double houselapsedms,RTvlapsedms,tileRlapsedms,vb2Wlapsedms;
   double WYTlapsedms,QWYTlapsedms,Qaddlapsedms;
   double YWTlapsedms,YWTClapsedms,Raddlapsedms;
   long long int qraddcnt = 0;
   long long int qrmulcnt = 0;
   long long int qrdivcnt = 0;
   long long int sqrtcnt = 0;

   if(verbose > 0) 
      cout << "-> GPU computes the blocked Householder QR ..." << endl;

   GPU_dbl_blocked_houseqr
      (nrows,ncols,szt,nbt,A,Q,R,
       &houselapsedms,&RTvlapsedms,&tileRlapsedms,&vb2Wlapsedms,
       &WYTlapsedms,&QWYTlapsedms,&Qaddlapsedms,
       &YWTlapsedms,&YWTClapsedms,&Raddlapsedms,&qrtimelapsed_d,
       &qraddcnt,&qrmulcnt,&qrdivcnt,&sqrtcnt,verbose);

   double bstimelapsed_d;
   double elapsedms,invlapsed,mullapsed,sublapsed;
   long long int bsaddcnt = 0;
   long long int bsmulcnt = 0;
   long long int bsdivcnt = 0;

   if(verbose > 0)
      cout << "-> GPU solves an upper triangular system ..." << endl;

   if(verbose > 0)
   {
      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++)
            cout << "R[" << i << "][" << j << "] : " << R[i][j] << endl;

      for(int i=0; i<nrows; i++)
         cout << "b[" << i << "] : " << b[i] << endl;
   }

   GPU_dbl_upper_tiled_solver
      (ncols,szt,nbt,R,b,x,
       &invlapsed,&mullapsed,&sublapsed,&elapsedms,&bstimelapsed_d,
       &bsaddcnt,&bsmulcnt,&bsdivcnt);
}

void GPU_dbl_update_rhs
 ( int nrows, int ncols, int szt, int nbt,
   double ***mat, double **rhs, bool verbose )
{
   if(verbose) cout << "in side GPU_dbl_update_rhs ..." << endl;
}

void GPU_dbl_bals_solve
 ( int dim, int degp1, int szt, int nbt,
   double ***mat, double **rhs, double **sol, int vrblvl )
{
   const int nrows = dim;
   const int ncols = dim;
   const bool bvrb = (vrblvl > 0);

   double **A = new double*[nrows];
   double **Q = new double*[nrows];
   double **R = new double*[nrows];

   double *b = new double[nrows];
   double *x = new double[ncols];

   for(int i=0; i<nrows; i++)
   {
      A[i] = new double[ncols];
      for(int j=0; j<ncols; j++) A[i][j] = mat[0][i][j];
      b[i] = rhs[0][i];
      Q[i] = new double[nrows];
      R[i] = new double[ncols];
      for(int j=0; j<ncols; j++) R[i][j] = mat[0][i][j];
   }

   GPU_dbl_qrbs_solve(nrows,ncols,szt,nbt,A,Q,R,b,x,bvrb);

   for(int j=0; j<ncols; j++) sol[0][j] = x[j];

   for(int stage=1; stage<nbt; stage++)
   {
      if(vrblvl > 0)
         cout << "stage " << stage << " in solve tail ..." << endl;

      GPU_dbl_update_rhs(nrows,ncols,szt,nbt,mat,rhs,bvrb);

      for(int i=0; i<nrows; i++) b[i] = rhs[stage][i];

      double bstimelapsed_d;
      double elapsedms,invlapsed,mullapsed,sublapsed;
      long long int bsaddcnt = 0;
      long long int bsmulcnt = 0;
      long long int bsdivcnt = 0;

      if(bvrb > 0)
         cout << "-> GPU solves an upper triangular system ..." << endl;

      GPU_dbl_upper_tiled_solver
         (ncols,szt,nbt,R,b,x,
          &invlapsed,&mullapsed,&sublapsed,&elapsedms,&bstimelapsed_d,
          &bsaddcnt,&bsmulcnt,&bsdivcnt);

      for(int j=0; j<ncols; j++) sol[stage][j] = x[j];
   }

   for(int i=0; i<nrows; i++)
   {
      free(A[i]); free(Q[i]); free(R[i]);
   }
   free(A); free(Q); free(R); free(b); free(x);
}
