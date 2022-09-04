// The file dbl2_bals_kernels.cu defines the functions with prototypes in
// the file dbl2_bals_kernels.h.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#ifdef gpufun
#include "double_double_gpufun.cu"
#endif
#include "dbl2_baqr_kernels.h"
#include "dbl2_tabs_kernels.h"
#include "dbl2_bals_kernels.h"

using namespace std;

__global__ void dbl2_bals_tail
 ( int ncols, int szt, double *Ahi, double *Alo,
   double *xhi, double *xlo, double *bhi, double *blo )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int idx = bdx*szt + tdx; // thread tdx updates b[idx]

   double Ajhi;             // register for Ahi[idx][j]
   double Ajlo;             // register for Alo[idx][j]
   double xjhi;             // register for xhi[j]
   double xjlo;             // register for xlo[j]
   double bihi = bhi[idx];  // register for bhi[idx]
   double bilo = blo[idx];  // register for blo[idx]
   double acchi,acclo;

   int offset = idx*ncols;

   for(int j=0; j<ncols; j++)
   {
      Ajhi = Ahi[offset+j];
      Ajlo = Alo[offset+j];
      xjhi = xhi[j];
      xjlo = xlo[j];
      // bi = bi - Aj*xj;
      ddg_mul(Ajhi,Ajlo,xjhi,xjlo,&acchi,&acclo);
      ddg_dec(&bihi,&bilo,acchi,acclo);
   }
   bhi[idx] = bihi;
   blo[idx] = bilo;
}

__global__ void dbl2_bals_qtb
 ( int ncols, int szt, double *Qthi, double *Qtlo,
   double *bhi, double *blo, double *rhi, double *rlo )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int idx = bdx*szt + tdx; // thread tdx computes r[idx]

   double Qjhi;           // register for Q^Thi[idx][j]
   double Qjlo;           // register for Q^Tlo[idx][j]
   double bjhi;           // register for bhi[j]
   double bjlo;           // register for blo[j]
   double rihi = 0.0;     // register for result, rhi[idx]
   double rilo = 0.0;     // register for result, rlo[idx]
   double acchi,acclo;

   int offset = idx*ncols;

   for(int j=0; j<ncols; j++)
   {
      Qjhi = Qthi[offset+j];
      Qjlo = Qtlo[offset+j];
      bjhi = bhi[j];
      bjlo = blo[j];
      // ri = ri + Qj*bj;
      ddg_mul(Qjhi,Qjlo,bjhi,bjlo,&acchi,&acclo);
      ddg_inc(&rihi,&rilo,acchi,acclo);
   }
   rhi[idx] = rihi;
   rlo[idx] = rilo;
}

void GPU_dbl2_bals_head
 ( int nrows, int ncols, int szt, int nbt,
   double **Ahi, double **Alo, double **Qhi, double **Qlo,
   double **Rhi, double **Rlo, double *bhi, double *blo, 
   double *xhi, double *xlo, bool verbose )
{
   double qrtimelapsed_d;
   double houselapsedms,RTvlapsedms,tileRlapsedms,vb2Wlapsedms;
   double WYTlapsedms,QWYTlapsedms,Qaddlapsedms;
   double YWTlapsedms,YWTClapsedms,Raddlapsedms;
   long long int qraddcnt = 0;
   long long int qrmulcnt = 0;
   long long int qrdivcnt = 0;
   long long int sqrtcnt = 0;

   if(verbose) 
      cout << "-> GPU computes the blocked Householder QR ..." << endl;

   GPU_dbl2_blocked_houseqr
      (nrows,ncols,szt,nbt,Ahi,Alo,Qhi,Qlo,Rhi,Rlo,
       &houselapsedms,&RTvlapsedms,&tileRlapsedms,&vb2Wlapsedms,
       &WYTlapsedms,&QWYTlapsedms,&Qaddlapsedms,
       &YWTlapsedms,&YWTClapsedms,&Raddlapsedms,&qrtimelapsed_d,
       &qraddcnt,&qrmulcnt,&qrdivcnt,&sqrtcnt,verbose);

   double bstimelapsed_d;
   double elapsedms,invlapsed,mullapsed,sublapsed;
   long long int bsaddcnt = 0;
   long long int bsmulcnt = 0;
   long long int bsdivcnt = 0;

   if(verbose)
      cout << "-> GPU solves an upper triangular system ..." << endl;

   if(verbose)
   {
      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++)
            cout << "R[" << i << "][" << j << "] : "
                 << Rhi[i][j] << "  " << Rlo[i][j] << endl;

      for(int i=0; i<nrows; i++)
         cout << "b[" << i << "] : "
              << bhi[i] << "  " << blo[i] << endl;
   }
   double **workRhi = new double*[nrows]; // work around because
   double **workRlo = new double*[nrows]; // solver modifies R ...

   for(int i=0; i<nrows; i++)
   {
      workRhi[i] = new double[ncols];
      workRlo[i] = new double[ncols];

      for(int j=0; j<ncols; j++)
      {
         workRhi[i][j] = Rhi[i][j];
         workRlo[i][j] = Rlo[i][j];
      }
   }
   GPU_dbl2_upper_tiled_solver
      (ncols,szt,nbt,workRhi,workRlo,bhi,blo,xhi,xlo,
       &invlapsed,&mullapsed,&sublapsed,&elapsedms,&bstimelapsed_d,
       &bsaddcnt,&bsmulcnt,&bsdivcnt);

   if(verbose)
   {
      cout << "-> after calling the GPU upper solver ..." << endl;
      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++)
            cout << "R[" << i << "][" << j << "] : "
                 << Rhi[i][j] << "  " << Rlo[i][j] << endl;

      for(int i=0; i<nrows; i++)
         cout << "b[" << i << "] : "
               << bhi[i] << "  " << blo[i] << endl;
   }
   for(int i=0; i<nrows; i++)
   {
      free(workRhi[i]);
      free(workRlo[i]);
   }
   free(workRhi); free(workRlo);
}

void GPU_dbl2_bals_tail
 ( int nrows, int ncols, int szt, int nbt, int degp1, int stage,
   double ***mathi, double ***matlo, double **rhshi, double **rhslo,
   double **solhi, double **sollo, bool verbose )
{
   if(verbose)
   {
      cout << "GPU_dbl2_bals_tail input blocks of rhs :" << endl;
      for(int k=0; k<degp1; k++)
      {
         for(int i=0; i<nrows; i++)
            cout << "rhs[" << k << "][" << i << "] : "
                 << rhshi[k][i] << "  " << rhslo[k][i] << endl;
      }
   }

   double *bhi_d;
   double *blo_d;
   const size_t szrhs = nrows*sizeof(double);
   cudaMalloc((void**)&bhi_d,szrhs);
   cudaMalloc((void**)&blo_d,szrhs);

   double *xhi_d;
   double *xlo_d;
   const size_t szsol = ncols*sizeof(double);
   cudaMalloc((void**)&xhi_d,szsol);
   cudaMalloc((void**)&xlo_d,szsol);
   cudaMemcpy(xhi_d,solhi[stage-1],szsol,cudaMemcpyHostToDevice);
   cudaMemcpy(xlo_d,sollo[stage-1],szsol,cudaMemcpyHostToDevice);

   double *Ahi_d;
   double *Alo_d;
   const size_t szmat = nrows*ncols*sizeof(double);
   cudaMalloc((void**)&Ahi_d,szmat);
   cudaMalloc((void**)&Alo_d,szmat);

   double *Ahi_h = new double[szmat];
   double *Alo_h = new double[szmat];

   for(int k=stage; k<degp1; k++)
   {
      if(verbose)
         cout << "GPU_dbl2_bals_tail launches " << nbt
              << " thread blocks in step " << k-stage << endl;

      int idx=0;
      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++)
         {
            Ahi_h[idx]   = mathi[k-stage+1][i][j];
            Alo_h[idx++] = matlo[k-stage+1][i][j];
         }

      cudaMemcpy(bhi_d,rhshi[k],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(blo_d,rhslo[k],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(Ahi_d,Ahi_h,szmat,cudaMemcpyHostToDevice);
      cudaMemcpy(Alo_d,Alo_h,szmat,cudaMemcpyHostToDevice);

      if(verbose)
         cout << "nbt = " << nbt << ", szt = " << szt
              << ", ncols = " << ncols << endl;

      dbl2_bals_tail<<<nbt,szt>>>
          (ncols,szt,Ahi_d,Alo_d,xhi_d,xlo_d,bhi_d,blo_d);
      
      if(verbose)
         cout << "copying block " << k << " of right hand side ..." << endl;

      cudaMemcpy(rhshi[k],bhi_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(rhslo[k],blo_d,szrhs,cudaMemcpyDeviceToHost);
   }
   free(Ahi_h); free(Alo_h);

   if(verbose)
   {
      cout << "GPU_dbl2_bals_tail copied blocks of rhs :" << endl;
      for(int k=0; k<degp1; k++)
      {
         for(int i=0; i<nrows; i++)
            cout << "rhs[" << k << "][" << i << "] : "
                 << rhshi[k][i] << "  " << rhslo[k][i] << endl;
      }
   }
}

void GPU_dbl2_bals_qtb
 ( int ncols, int szt, int nbt,
   double **Qhi, double **Qlo, double *bhi, double *blo, bool verbose )
{
   double *bhi_d;
   double *blo_d;
   const size_t szrhs = ncols*sizeof(double);
   cudaMalloc((void**)&bhi_d,szrhs);
   cudaMalloc((void**)&blo_d,szrhs);
   cudaMemcpy(bhi_d,bhi,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(blo_d,blo,szrhs,cudaMemcpyHostToDevice);

   double *rhi_d;
   double *rlo_d;
   const size_t szsol = ncols*sizeof(double);
   cudaMalloc((void**)&rhi_d,szsol);
   cudaMalloc((void**)&rlo_d,szsol);

   double *Qthi_d;
   double *Qtlo_d;
   const size_t szmat = ncols*ncols*sizeof(double);
   cudaMalloc((void**)&Qthi_d,szmat);
   cudaMalloc((void**)&Qtlo_d,szmat);

   double *Qthi_h = new double[szmat];
   double *Qtlo_h = new double[szmat];

   int idx=0;
   for(int i=0; i<ncols; i++)
      for(int j=0; j<ncols; j++)
      {
         Qthi_h[idx]   = Qhi[j][i];
         Qtlo_h[idx++] = Qlo[j][i];
      }

   cudaMemcpy(Qthi_d,Qthi_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(Qtlo_d,Qtlo_h,szmat,cudaMemcpyHostToDevice);

   dbl2_bals_qtb<<<nbt,szt>>>
      (ncols,szt,Qthi_d,Qtlo_d,bhi_d,blo_d,rhi_d,rlo_d);

   cudaMemcpy(bhi,rhi_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(blo,rlo_d,szrhs,cudaMemcpyDeviceToHost);

   free(Qthi_h); free(Qtlo_h);
}

void GPU_dbl2_bals_solve
 ( int dim, int degp1, int szt, int nbt,
   double ***mathi, double ***matlo, double **Qhi, double **Qlo,
   double **Rhi, double **Rlo, double **rhshi, double **rhslo,
   double **solhi, double **sollo, int vrblvl )
{
   const int nrows = dim;
   const int ncols = dim;
   const bool bvrb = (vrblvl > 0);

   double **Ahi = new double*[nrows];
   double **Alo = new double*[nrows];
   double *bhi = new double[nrows];
   double *blo = new double[nrows];
   double *xhi = new double[ncols];
   double *xlo = new double[ncols];

   double **workRhi = new double*[nrows]; // GPU upper solver changes R
   double **workRlo = new double*[nrows];
   for(int i=0; i<nrows; i++)
   {
      workRhi[i] = new double[ncols];
      workRlo[i] = new double[ncols];
   }
   if(vrblvl)
   {
      cout << "GPU_dbl2_bals_solve blocks of rhs :" << endl;
      for(int k=0; k<degp1; k++)
      {
         for(int i=0; i<nrows; i++)
            cout << "rhs[" << k << "][" << i << "] : "
                 << rhshi[k][i] << "  " << rhslo[k][i] << endl;
      }
   }

   for(int i=0; i<nrows; i++)
   {
      Ahi[i] = new double[ncols];
      Alo[i] = new double[ncols];

      for(int j=0; j<ncols; j++)
      {
         Ahi[i][j] = mathi[0][i][j];
         Alo[i][j] = matlo[0][i][j];
      }
      bhi[i] = rhshi[0][i];
      blo[i] = rhslo[0][i];

      for(int j=0; j<ncols; j++)
      {
         Rhi[i][j] = mathi[0][i][j];
         Rlo[i][j] = matlo[0][i][j];
      }
   }

   GPU_dbl2_bals_head
      (nrows,ncols,szt,nbt,Ahi,Alo,Qhi,Qlo,Rhi,Rlo,bhi,blo,xhi,xlo,bvrb);

   for(int j=0; j<ncols; j++)
   {
      solhi[0][j] = xhi[j];
      sollo[0][j] = xlo[j];
   }
   for(int stage=1; stage<degp1; stage++)
   {
      if(vrblvl > 0)
         cout << "stage " << stage << " in solve tail ..." << endl;

      GPU_dbl2_bals_tail
         (nrows,ncols,szt,nbt,degp1,stage,
          mathi,matlo,rhshi,rhslo,solhi,sollo,bvrb);

      if(vrblvl > 0)
      {
         cout << "blocks of rhs before assignment :" << endl;
         for(int k=0; k<degp1; k++)
         {
            for(int i=0; i<nrows; i++)
               cout << "rhs[" << k << "][" << i << "] : "
                    << rhshi[k][i] << "  " << rhslo[k][i] << endl;
         }
      }

      for(int i=0; i<nrows; i++) 
      {
         cout << "assigning component " << i
              << ", stage = " << stage << endl;
         bhi[i] = rhshi[stage][i];
         blo[i] = rhslo[stage][i];
         cout << "b[" << i << "] : "
              << bhi[i] << "  " << blo[i] << endl;
      }
      double bstimelapsed_d;
      double elapsedms,invlapsed,mullapsed,sublapsed;
      long long int bsaddcnt = 0;
      long long int bsmulcnt = 0;
      long long int bsdivcnt = 0;

      if(vrblvl > 0)
         cout << "-> GPU multiplies rhs with Q^T ..." << endl;

      GPU_dbl2_bals_qtb(ncols,szt,nbt,Qhi,Qlo,bhi,blo,bvrb);

      if(vrblvl > 0)
      {
         for(int i=0; i<nrows; i++)
            cout << "Qtb[" << i << "] : "
                 << bhi[i] << "  " << blo[i] << endl;
      }
      if(vrblvl > 0)
      {
         cout << "-> GPU solves an upper triangular system ..." << endl;
 
         for(int i=0; i<nrows; i++)
            for(int j=0; j<ncols; j++)
               cout << "R[" << i << "][" << j << "] : "
                    << Rhi[i][j] << "  " << Rlo[i][j] << endl;
      }
      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++)
         {
            workRhi[i][j] = Rhi[i][j];
            workRlo[i][j] = Rlo[i][j];
         }

      GPU_dbl2_upper_tiled_solver
         (ncols,szt,nbt,workRhi,workRlo,bhi,blo,xhi,xlo,
          &invlapsed,&mullapsed,&sublapsed,&elapsedms,&bstimelapsed_d,
          &bsaddcnt,&bsmulcnt,&bsdivcnt);

     if(vrblvl > 0)
        for(int i=0; i<ncols; i++)
           cout << "x[" << i << "] : "
                << xhi[i] << "  " << xlo[i] << endl;

      for(int j=0; j<ncols; j++)
      {
         solhi[stage][j] = xhi[j];
         sollo[stage][j] = xlo[j];
      }
   }

   for(int i=0; i<nrows; i++)
   {
      free(Ahi[i]); free(workRhi[i]);
      free(Alo[i]); free(workRlo[i]);
   }
   free(Ahi); free(bhi); free(xhi); free(workRhi);
   free(Alo); free(blo); free(xlo); free(workRlo);
}
