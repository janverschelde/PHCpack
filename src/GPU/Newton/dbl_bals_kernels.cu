// The file dbl_bals_kernels.cu defines the functions with prototypes in
// the file dbl_bals_kernels.h.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#include "dbl_baqr_kernels.h"
#include "dbl_tabs_kernels.h"
#include "dbl_bals_kernels.h"

using namespace std;

__global__ void dbl_bals_tail
 ( int ncols, int szt, double *A, double *x, double *b )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int idx = bdx*szt + tdx; // thread tdx updates b[idx]

   double Aj;           // register for A[idx][j]
   double xj;           // register for hold x[j]
   double bi = b[idx];  // register for b[idx]

   int offset = idx*ncols;

   for(int j=0; j<ncols; j++)
   {
      Aj = A[offset+j];
      xj = x[j];
      bi = bi - Aj*xj;
   }
   b[idx] = bi;
}

__global__ void dbl_bals_qtb
 ( int ncols, int szt, double *Qt, double *b, double *r )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int idx = bdx*szt + tdx; // thread tdx updates b[idx]

   double Qj;           // register for Q^T[idx][j]
   double bj;           // register for hold b[j]
   double ri = 0.0;     // register for result, b[idx]

   int offset = idx*ncols;

   for(int j=0; j<ncols; j++)
   {
      Qj = Qt[offset+j];
      bj = b[j];
      ri = ri + Qj*bj;
   }
   r[idx] = ri;
}

void GPU_dbl_bals_head
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

void GPU_dbl_bals_tail
 ( int nrows, int ncols, int szt, int nbt, int degp1, int stage,
   double ***mat, double **rhs, double **sol, bool verbose )
{
   if(verbose)
   {
      cout << "GPU_dbl_bals_tail input blocks of rhs :" << endl;
      for(int k=0; k<degp1; k++)
      {
         for(int i=0; i<nrows; i++)
            cout << "rhs[" << k << "][" << i << "] : " << rhs[k][i] << endl;
      }
   }

   double *b_d;
   const size_t szrhs = nrows*sizeof(double);
   cudaMalloc((void**)&b_d,szrhs);

   double *x_d;
   const size_t szsol = ncols*sizeof(double);
   cudaMalloc((void**)&x_d,szsol);
   cudaMemcpy(x_d,&sol[stage-1],szsol,cudaMemcpyHostToDevice);

   double *A_d;
   const size_t szmat = nrows*ncols*sizeof(double);
   cudaMalloc((void**)&A_d,szmat);

   double *A_h = new double[szmat];

   for(int k=stage; k<degp1; k++)
   {
      if(verbose)
         cout << "GPU_dbl_bals_tail launches " << nbt
              << " thread blocks in step " << k-stage << endl;

      int idx=0;
      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++) A_h[idx++] = mat[k][i][j];
      
      cudaMemcpy(b_d,&(rhs[k][0]),szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(A_d,A_h,szmat,cudaMemcpyHostToDevice);

      if(verbose)
         cout << "nbt = " << nbt << ", szt = " << szt
              << ", ncols = " << ncols << endl;

      dbl_bals_tail<<<nbt,szt>>>(ncols,szt,A_d,x_d,b_d);
      
      if(verbose)
         cout << "copying block " << k << " of right hand side ..." << endl;

      cudaMemcpy(&(rhs[k][0]),b_d,szrhs,cudaMemcpyDeviceToHost);
   }
   free(A_h);

   if(verbose)
   {
      cout << "GPU_dbl_bals_tail copied blocks of rhs :" << endl;
      for(int k=0; k<degp1; k++)
      {
         for(int i=0; i<nrows; i++)
            cout << "rhs[" << k << "][" << i << "] : " << rhs[k][i] << endl;
      }
   }
}

void GPU_dbl_bals_qtb
 ( int ncols, int szt, int nbt, double **Q, double *b, bool verbose )
{
   double *b_d;
   const size_t szrhs = ncols*sizeof(double);
   cudaMalloc((void**)&b_d,szrhs);

   double *r_d;
   const size_t szsol = ncols*sizeof(double);
   cudaMalloc((void**)&r_d,szsol);

   double *Qt_d;
   const size_t szmat = ncols*ncols*sizeof(double);
   cudaMalloc((void**)&Qt_d,szmat);

   double *Qt_h = new double[szmat];

   int idx=0;
   for(int i=0; i<ncols; i++)
      for(int j=0; j<ncols; j++) Qt_h[idx++] = Q[j][i];

   cudaMemcpy(b_d,b,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(Qt_d,Qt_h,szmat,cudaMemcpyHostToDevice);

   dbl_bals_qtb<<<nbt,szt>>>(ncols,szt,Qt_d,b_d,r_d);

   cudaMemcpy(b,r_d,szrhs,cudaMemcpyDeviceToHost);

   free(Qt_h);
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

   if(vrblvl)
   {
      cout << "GPU_dbl_bals_solve blocks of rhs :" << endl;
      for(int k=0; k<degp1; k++)
      {
         for(int i=0; i<nrows; i++)
            cout << "rhs[" << k << "][" << i << "] : " << rhs[k][i] << endl;
      }
   }

   for(int i=0; i<nrows; i++)
   {
      A[i] = new double[ncols];
      for(int j=0; j<ncols; j++) A[i][j] = mat[0][i][j];
      b[i] = rhs[0][i];
      Q[i] = new double[nrows];
      R[i] = new double[ncols];
      for(int j=0; j<ncols; j++) R[i][j] = mat[0][i][j];
   }

   GPU_dbl_bals_head(nrows,ncols,szt,nbt,A,Q,R,b,x,bvrb);

   for(int j=0; j<ncols; j++) sol[0][j] = x[j];

   for(int stage=1; stage<degp1; stage++)
   {
      if(vrblvl > 0)
         cout << "stage " << stage << " in solve tail ..." << endl;

      GPU_dbl_bals_tail(nrows,ncols,szt,nbt,degp1,stage,mat,rhs,sol,bvrb);

      if(vrblvl)
      {
         cout << "blocks of rhs before assignment :" << endl;
         for(int k=0; k<degp1; k++)
         {
            for(int i=0; i<nrows; i++)
               cout << "rhs[" << k << "][" << i << "] : " << rhs[k][i] << endl;
         }
      }

      for(int i=0; i<nrows; i++) 
      {
         cout << "assigning component " << i
              << ", stage = " << stage << endl;
         b[i] = rhs[stage][i];
      }
      double bstimelapsed_d;
      double elapsedms,invlapsed,mullapsed,sublapsed;
      long long int bsaddcnt = 0;
      long long int bsmulcnt = 0;
      long long int bsdivcnt = 0;

      if(bvrb > 0)
         cout << "-> GPU multiplies rhs with Q^T ..." << endl;

      GPU_dbl_bals_qtb(ncols,szt,nbt,Q,b,bvrb);

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
