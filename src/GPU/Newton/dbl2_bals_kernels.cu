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
#include "dbl2_tail_kernels.h"
#include "dbl2_bals_kernels.h"
#include "write_dbl2_bstimeflops.h"
#include "write_dbl2_qrtimeflops.h"
#include "dbl_onenorms_host.h"

using namespace std;

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

__global__ void cmplx2_bals_qhb
 ( int ncols, int szt,
   double *QHrehi, double *QHrelo, double *QHimhi, double *QHimlo,
   double *brehi, double *brelo, double *bimhi, double *bimlo,
   double *rrehi, double *rrelo, double *rimhi, double *rimlo )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int idx = bdx*szt + tdx; // thread tdx computes r[idx]

   double Qjrehi;         // registers for real part of Q^H[idx][j]
   double Qjrelo;
   double Qjimhi;         // registers for imaginary part of Q^H[idx][j]
   double Qjimlo;
   double bjrehi;         // register for brehi[j]
   double bjrelo;         // register for brelo[j]
   double bjimhi;         // register for bimhi[j]
   double bjimlo;         // register for bimlo[j]
   double rirehi = 0.0;   // register for result, rrehi[idx]
   double rirelo = 0.0;   // register for result, rrelo[idx]
   double riimhi = 0.0;   // register for result, rimhi[idx]
   double riimlo = 0.0;   // register for result, rimlo[idx]
   double acchi,acclo;

   int offset = idx*ncols;

   for(int j=0; j<ncols; j++)
   {
      Qjrehi = QHrehi[offset+j];
      Qjrelo = QHrelo[offset+j];
      Qjimhi = QHimhi[offset+j];
      Qjimlo = QHimlo[offset+j];
      bjrehi = brehi[j];
      bjrelo = brelo[j];
      bjimhi = bimhi[j];
      bjimlo = bimlo[j];
      // ri = ri + Qj*bj;
      // zre = Qjre*bjre - Qjim*bjim;
      // rire = rire + zre;
      ddg_mul(Qjrehi,Qjrelo,bjrehi,bjrelo,&acchi,&acclo);
      ddg_inc(&rirehi,&rirelo,acchi,acclo);
      ddg_mul(Qjimhi,Qjimlo,bjimhi,bjimlo,&acchi,&acclo);
      ddg_dec(&rirehi,&rirelo,acchi,acclo);
      // zim = Qjre*bjim + Qjim*bjre;
      // riim = riim + zim;
      ddg_mul(Qjrehi,Qjrelo,bjimhi,bjimlo,&acchi,&acclo);
      ddg_inc(&riimhi,&riimlo,acchi,acclo);
      ddg_mul(Qjimhi,Qjimlo,bjrehi,bjrelo,&acchi,&acclo);
      ddg_inc(&riimhi,&riimlo,acchi,acclo);
   }
   rrehi[idx] = rirehi; rrelo[idx] = rirelo;
   rimhi[idx] = riimhi; rimlo[idx] = riimlo;
}

void GPU_dbl2_bals_head
 ( int nrows, int ncols, int szt, int nbt,
   double **Ahi, double **Alo, double **Qhi, double **Qlo,
   double **Rhi, double **Rlo, double *bhi, double *blo, 
   double *xhi, double *xlo,
   double *totqrlapsedms, double *totqtblapsedms, double *totbslapsedms,
   int vrblvl )
{
   double qrtimelapsed_d;
   double houselapsedms,RTvlapsedms,tileRlapsedms,vb2Wlapsedms;
   double WYTlapsedms,QWYTlapsedms,Qaddlapsedms;
   double YWTlapsedms,YWTClapsedms,Raddlapsedms;
   long long int qraddcnt = 0;
   long long int qrmulcnt = 0;
   long long int qrdivcnt = 0;
   long long int sqrtcnt = 0;
   bool verbose = (vrblvl > 1);

   if(vrblvl > 0) 
      cout << "-> GPU computes the blocked Householder QR ..." << endl;

   GPU_dbl2_blocked_houseqr
      (nrows,ncols,szt,nbt,Ahi,Alo,Qhi,Qlo,Rhi,Rlo,
       &houselapsedms,&RTvlapsedms,&tileRlapsedms,&vb2Wlapsedms,
       &WYTlapsedms,&QWYTlapsedms,&Qaddlapsedms,
       &YWTlapsedms,&YWTClapsedms,&Raddlapsedms,&qrtimelapsed_d,
       &qraddcnt,&qrmulcnt,&qrdivcnt,&sqrtcnt,verbose);

   *totqrlapsedms = *totqrlapsedms + houselapsedms + RTvlapsedms
      + tileRlapsedms + vb2Wlapsedms + WYTlapsedms + QWYTlapsedms
      + Qaddlapsedms + YWTlapsedms + YWTClapsedms + Raddlapsedms;

   if(vrblvl > 0)
      write_dbl2_qrtimeflops
         (0,nrows,ncols,houselapsedms,RTvlapsedms,tileRlapsedms,vb2Wlapsedms,
          WYTlapsedms,QWYTlapsedms,Qaddlapsedms,YWTlapsedms,YWTClapsedms,
          Raddlapsedms,qrtimelapsed_d,qraddcnt,qrmulcnt,qrdivcnt,sqrtcnt);

   if(vrblvl > 0) cout << "-> GPU multiplies rhs with Q^T ..." << endl;

   GPU_dbl2_bals_qtb(ncols,szt,nbt,Qhi,Qlo,bhi,blo,totqtblapsedms,vrblvl);

   if(vrblvl > 1)
   {
      for(int i=0; i<nrows; i++)
         cout << "Qtb[" << i << "] : "
              << bhi[i] << "  " << blo[i] << endl;
   }
   double bstimelapsed_d;
   double elapsedms,invlapsed,mullapsed,sublapsed;
   long long int bsaddcnt = 0;
   long long int bsmulcnt = 0;
   long long int bsdivcnt = 0;

   if(vrblvl > 0)
      cout << "-> GPU solves an upper triangular system ..." << endl;

   if(vrblvl > 1)
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

   *totbslapsedms += elapsedms;

   if(vrblvl > 0)
      write_dbl2_bstimeflops
         (szt,nbt,0,invlapsed,mullapsed,sublapsed,elapsedms,bstimelapsed_d,
          bsaddcnt,bsmulcnt,bsdivcnt);

   if(vrblvl > 1)
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

void GPU_cmplx2_bals_head
 ( int nrows, int ncols, int szt, int nbt,
   double **Arehi, double **Arelo, double **Aimhi, double **Aimlo,
   double **Qrehi, double **Qrelo, double **Qimhi, double **Qimlo,
   double **Rrehi, double **Rrelo, double **Rimhi, double **Rimlo, 
   double *brehi, double *brelo, double *bimhi, double *bimlo,
   double *xrehi, double *xrelo, double *ximhi, double *ximlo,
   double *totqrlapsedms, double *totqtblapsedms, double *totbslapsedms,
   int vrblvl )
{
   double qrtimelapsed_d;
   double houselapsedms,RTvlapsedms,tileRlapsedms,vb2Wlapsedms;
   double WYTlapsedms,QWYTlapsedms,Qaddlapsedms;
   double YWTlapsedms,YWTClapsedms,Raddlapsedms;
   long long int qraddcnt = 0;
   long long int qrmulcnt = 0;
   long long int qrdivcnt = 0;
   long long int sqrtcnt = 0;
   bool verbose = (vrblvl > 1);

   if(vrblvl > 0) 
      cout << "-> GPU computes the blocked Householder QR ..." << endl;

   GPU_cmplx2_blocked_houseqr
      (nrows,ncols,szt,nbt,Arehi,Arelo,Aimhi,Aimlo,
       Qrehi,Qrelo,Qimhi,Qimlo,Rrehi,Rrelo,Rimhi,Rimlo,
       &houselapsedms,&RTvlapsedms,&tileRlapsedms,&vb2Wlapsedms,
       &WYTlapsedms,&QWYTlapsedms,&Qaddlapsedms,
       &YWTlapsedms,&YWTClapsedms,&Raddlapsedms,&qrtimelapsed_d,
       &qraddcnt,&qrmulcnt,&qrdivcnt,&sqrtcnt,verbose);

   *totqrlapsedms = *totqrlapsedms + houselapsedms + RTvlapsedms
      + tileRlapsedms + vb2Wlapsedms + WYTlapsedms + QWYTlapsedms
      + Qaddlapsedms + YWTlapsedms + YWTClapsedms + Raddlapsedms;

   if(vrblvl > 0)
      write_dbl2_qrtimeflops
         (1,nrows,ncols,houselapsedms,RTvlapsedms,tileRlapsedms,vb2Wlapsedms,
          WYTlapsedms,QWYTlapsedms,Qaddlapsedms,YWTlapsedms,YWTClapsedms,
          Raddlapsedms,qrtimelapsed_d,qraddcnt,qrmulcnt,qrdivcnt,sqrtcnt);

   if(vrblvl > 0) cout << "-> GPU multiplies rhs with Q^H ..." << endl;

   GPU_cmplx2_bals_qhb
      (ncols,szt,nbt,Qrehi,Qrelo,Qimhi,Qimlo,brehi,brelo,bimhi,bimlo,
       totqtblapsedms,vrblvl);

   if(vrblvl > 1)
   {
      for(int i=0; i<nrows; i++)
         cout << "QHb[" << i << "] : "
              << brehi[i] << "  " << brelo[i] << endl << "  "
              << bimhi[i] << "  " << bimlo[i] << endl;
   }

   double bstimelapsed_d;
   double elapsedms,invlapsed,mullapsed,sublapsed;
   long long int bsaddcnt = 0;
   long long int bsmulcnt = 0;
   long long int bsdivcnt = 0;

   if(vrblvl > 0)
      cout << "-> GPU solves an upper triangular system ..." << endl;

   if(vrblvl > 1)
   {
      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++)
            cout << "R[" << i << "][" << j << "] : "
                 << Rrehi[i][j] << "  " << Rrelo[i][j] << endl << "  "
                 << Rimhi[i][j] << "  " << Rimlo[i][j] << endl;

      for(int i=0; i<nrows; i++)
         cout << "b[" << i << "] : "
              << brehi[i] << "  " << brelo[i] << "  "
              << bimhi[i] << "  " << bimlo[i] << endl;
   }
   double **workRrehi = new double*[nrows]; // work around ...
   double **workRrelo = new double*[nrows];
   double **workRimhi = new double*[nrows];
   double **workRimlo = new double*[nrows];

   for(int i=0; i<nrows; i++)
   {
      workRrehi[i] = new double[ncols]; workRrelo[i] = new double[ncols];
      workRimhi[i] = new double[ncols]; workRimlo[i] = new double[ncols];

      for(int j=0; j<ncols; j++)
      {
         workRrehi[i][j] = Rrehi[i][j]; workRrelo[i][j] = Rrelo[i][j];
         workRimhi[i][j] = Rimhi[i][j]; workRimlo[i][j] = Rimlo[i][j];
      }
   }
   GPU_cmplx2_upper_tiled_solver
      (ncols,szt,nbt,workRrehi,workRrelo,workRimhi,workRimlo,
       brehi,brelo,bimhi,bimlo,xrehi,xrelo,ximhi,ximlo,
       &invlapsed,&mullapsed,&sublapsed,&elapsedms,&bstimelapsed_d,
       &bsaddcnt,&bsmulcnt,&bsdivcnt);

   *totbslapsedms += elapsedms;

   if(vrblvl > 0)
      write_dbl2_bstimeflops
         (szt,nbt,1,invlapsed,mullapsed,sublapsed,elapsedms,bstimelapsed_d,
          bsaddcnt,bsmulcnt,bsdivcnt);

   if(vrblvl > 1)
   {
      cout << "-> after calling the GPU upper solver ..." << endl;
      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++)
            cout << "R[" << i << "][" << j << "] : "
                 << Rrehi[i][j] << "  " << Rrelo[i][j] << endl << "  "
                 << Rimhi[i][j] << "  " << Rimlo[i][j] << endl;

      for(int i=0; i<nrows; i++)
         cout << "b[" << i << "] : "
              << brehi[i] << "  " << brelo[i] << endl << "  "
              << bimhi[i] << "  " << bimlo[i] << endl;
   }
   for(int i=0; i<nrows; i++)
   {
      free(workRrehi[i]); free(workRimhi[i]);
      free(workRrelo[i]); free(workRimlo[i]);
   }
   free(workRrehi); free(workRimhi);
   free(workRrelo); free(workRimlo);
}

void write_dbl2_qtbflops ( int ctype, int ncols, float lapsms )
{
   cout << fixed << setprecision(3);
   cout << "Time spent for Q^T*b : " << lapsms << " milliseconds." << endl;

   long long int flopcnt;
   const long long int longncols2 = ncols*ncols; // to avoid overflow
   if(ctype == 0)
      flopcnt = 20*longncols2 + 23*longncols2;
      // as many + as * in one inner product
   else
      flopcnt = 4*20*longncols2 + 4*23*longncols2;
      // for complex *: 2 ops for +, 6 for *, which is 8 in total

   cout << "    Total number of floating-point operations : "
        << flopcnt << endl;

   long long int bytecnt;

   if(ctype == 0)
      bytecnt = 2*ncols*ncols;
   else
      bytecnt = 4*ncols*ncols;

   cout << "    Total number of bytes : " << bytecnt << endl;

   double intensity = ((double) flopcnt)/bytecnt;
   cout << "     Arithmetic intensity : "
        << scientific << setprecision(3) << intensity
        << " #flops/#bytes" << endl;

   double kernflops = 1000.0*((double) flopcnt)/lapsms;
   // double wallflops = ((double) flopcnt)/timelapsed;
   const int gigacnt = pow(2.0,30);

   cout << "Kernel Time Flops : "
        << scientific << setprecision(3) << kernflops;
   cout << fixed << setprecision(3)
        << " = " << kernflops/gigacnt << " Gigaflops" << endl;
/*
   cout << " Wall Clock Flops : "
        << scientific << setprecision(3) << wallflops;
   cout << fixed << setprecision(3)
        << " = " << wallflops/gigacnt << " Gigaflops" << endl;
 */
}

void GPU_dbl2_bals_qtb
 ( int ncols, int szt, int nbt,
   double **Qhi, double **Qlo, double *bhi, double *blo,
   double *totqtblapsedms, int vrblvl )
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

   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;

   cudaEventRecord(start);
   dbl2_bals_qtb<<<nbt,szt>>>
      (ncols,szt,Qthi_d,Qtlo_d,bhi_d,blo_d,rhi_d,rlo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);

   *totqtblapsedms += milliseconds;

   cudaMemcpy(bhi,rhi_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(blo,rlo_d,szrhs,cudaMemcpyDeviceToHost);

   if(vrblvl > 0) write_dbl2_qtbflops(0,ncols,milliseconds);

   cudaFree(bhi_d); cudaFree(blo_d);
   cudaFree(rhi_d); cudaFree(rlo_d);
   cudaFree(Qthi_d); cudaFree(Qtlo_d);

   free(Qthi_h); free(Qtlo_h);
}

void GPU_cmplx2_bals_qhb
 ( int ncols, int szt, int nbt,
   double **Qrehi, double **Qrelo, double **Qimhi, double **Qimlo,
   double *brehi, double *brelo, double *bimhi, double *bimlo,
   double *totqtblapsedms, int vrblvl )
{
   double *brehi_d;
   double *brelo_d;
   double *bimhi_d;
   double *bimlo_d;
   const size_t szrhs = ncols*sizeof(double);
   cudaMalloc((void**)&brehi_d,szrhs);
   cudaMalloc((void**)&brelo_d,szrhs);
   cudaMalloc((void**)&bimhi_d,szrhs);
   cudaMalloc((void**)&bimlo_d,szrhs);
   cudaMemcpy(brehi_d,brehi,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(brelo_d,brelo,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(bimhi_d,bimhi,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(bimlo_d,bimlo,szrhs,cudaMemcpyHostToDevice);

   double *rrehi_d;
   double *rrelo_d;
   double *rimhi_d;
   double *rimlo_d;
   const size_t szsol = ncols*sizeof(double);
   cudaMalloc((void**)&rrehi_d,szsol);
   cudaMalloc((void**)&rrelo_d,szsol);
   cudaMalloc((void**)&rimhi_d,szsol);
   cudaMalloc((void**)&rimlo_d,szsol);

   double *QHrehi_d;
   double *QHrelo_d;
   double *QHimhi_d;
   double *QHimlo_d;
   const size_t szmat = ncols*ncols*sizeof(double);
   cudaMalloc((void**)&QHrehi_d,szmat);
   cudaMalloc((void**)&QHrelo_d,szmat);
   cudaMalloc((void**)&QHimhi_d,szmat);
   cudaMalloc((void**)&QHimlo_d,szmat);

   double *QHrehi_h = new double[szmat];
   double *QHrelo_h = new double[szmat];
   double *QHimhi_h = new double[szmat];
   double *QHimlo_h = new double[szmat];

   int idx=0;
   for(int i=0; i<ncols; i++)
      for(int j=0; j<ncols; j++)
      {
         QHrehi_h[idx]   = Qrehi[j][i];
         QHrelo_h[idx]   = Qrelo[j][i];
         QHimhi_h[idx]   = - Qimhi[j][i]; // Hermitian transpose !
         QHimlo_h[idx++] = - Qimlo[j][i]; // Hermitian transpose !
      }

   cudaMemcpy(QHrehi_d,QHrehi_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(QHrelo_d,QHrelo_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(QHimhi_d,QHimhi_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(QHimlo_d,QHimlo_h,szmat,cudaMemcpyHostToDevice);

   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;

   cudaEventRecord(start);
   cmplx2_bals_qhb<<<nbt,szt>>>
      (ncols,szt,QHrehi_d,QHrelo_d,QHimhi_d,QHimlo_d,
       brehi_d,brelo_d,bimhi_d,bimlo_d,rrehi_d,rrelo_d,rimhi_d,rimlo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);

   *totqtblapsedms += milliseconds;

   cudaMemcpy(brehi,rrehi_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(brelo,rrelo_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(bimhi,rimhi_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(bimlo,rimlo_d,szrhs,cudaMemcpyDeviceToHost);

   if(vrblvl > 0) write_dbl2_qtbflops(1,ncols,milliseconds);

   cudaFree(brehi_d); cudaFree(brelo_d); cudaFree(bimhi_d); cudaFree(bimlo_d);
   cudaFree(rrehi_d); cudaFree(rrelo_d); cudaFree(rimhi_d); cudaFree(rimlo_d);
   cudaFree(QHrehi_d); cudaFree(QHrelo_d);
   cudaFree(QHimhi_d); cudaFree(QHimlo_d);

   free(QHrehi_h); free(QHimhi_h);
   free(QHrelo_h); free(QHimlo_h);
}

void GPU_dbl2_bals_solve
 ( int dim, int degp1, int szt, int nbt, int tailidx,
   double ***mathi, double ***matlo, double **Qhi, double **Qlo,
   double **Rhi, double **Rlo, double **rhshi, double **rhslo,
   double **solhi, double **sollo,
   bool *zeroQ, bool *noqr, int *upidx, int *bsidx, int *newtail,
   double *totqrlapsedms, double *totqtblapsedms, double *totbslapsedms,
   double *totupdlapsedms, int vrblvl )
{
   const int nrows = dim;
   const int ncols = dim;
   int skipupcnt = 0; // counts the skipped updates
   int skipbscnt = 0; // counts the skipped backsubstitutions
   double prevnorm = 1.0e+99;
   bool firstbs = true; // at the first back substitution
   *newtail = degp1; // in case no back substitution happens

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
   if(vrblvl > 1)
   {
      cout << "GPU_dbl2_bals_solve blocks of rhs :" << endl;
      for(int k=0; k<degp1; k++)
      {
         for(int i=0; i<nrows; i++)
            cout << "rhs[" << k << "][" << i << "] : "
                 << rhshi[k][i] << "  " << rhslo[k][i] << endl;
      }
   }
   double nrm;

   if((*noqr) && (!*zeroQ))
   {
      if(vrblvl > 0) cout << "-> skipping GPU_dbl2_bals_head ..." << endl;
   }
   else
   {
      CPU_dbl_onenorm(nrows,rhshi[0],&nrm);
      if(vrblvl > 0)
         cout << scientific << setprecision(3)
              << "1-norm of b : " << nrm << endl;

      if(nrm < 1.0e-28)
      {
         if(*zeroQ)
         {
            if(vrblvl > 0)
               cout << "-> no skipping GPU_dbl2_bals_head because zeroQ"
                    << endl;

            *noqr = false;
         }
         else
         {
            if(vrblvl > 0)
               cout << "-> skip call to GPU_dbl2_bals_head ..." << endl;

            *noqr = true;

            for(int j=0; j<ncols; j++)
            {
               solhi[0][j] = 0.0; sollo[0][j] = 0.0;
            }
         }
      }
      if(!*noqr)
      {
         if(vrblvl > 0) cout << "-> calling GPU_dbl2_bals_head ..." << endl;

         double **Ahi = new double*[nrows];
         double **Alo = new double*[nrows];

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
            (nrows,ncols,szt,nbt,Ahi,Alo,Qhi,Qlo,Rhi,Rlo,bhi,blo,
             xhi,xlo,totqrlapsedms,totqtblapsedms,totbslapsedms,vrblvl);

         *zeroQ = false;

         for(int j=0; j<ncols; j++)
         {
            solhi[0][j] = xhi[j];
            sollo[0][j] = xlo[j];
         }
         for(int i=0; i<nrows; i++)
         {
            free(Ahi[i]); free(Alo[i]);
         }
         free(Ahi); free(Alo);
      }
   }
   for(int stage=tailidx; stage<degp1; stage++)
   {
      if(vrblvl > 0)
         cout << "stage " << stage << " in solve tail ..." << endl;

      double *xshi = solhi[stage-1];       // solution to do the update with
      CPU_dbl_onenorm(dim,xshi,&nrm);
      if(vrblvl > 0)
         cout << scientific << setprecision(3)
              << "1-norm of x[" << stage-1 << "] : " << nrm << endl;

      if(nrm < 1.0e-28)
      {
         skipupcnt = skipupcnt + 1;

         if(vrblvl > 0)
            cout << "skip update with x[" << stage-1 << "] ..." << endl;
      }
      else
      {
         if(vrblvl > 0)
            cout << "updating with x[" << stage-1 << "] ..." << endl;

         GPU_dbl2_bals_tail
            (nrows,ncols,szt,nbt,degp1,stage,
             mathi,matlo,rhshi,rhslo,solhi,sollo,totupdlapsedms,vrblvl);

         if(vrblvl > 1)
         {
            cout << "blocks of rhs before assignment :" << endl;
            for(int k=0; k<degp1; k++)
            {
               for(int i=0; i<nrows; i++)
                  cout << "rhs[" << k << "][" << i << "] : "
                       << rhshi[k][i] << "  " << rhslo[k][i] << endl;
            }
         }
      }
      for(int i=0; i<nrows; i++) 
      {
         // cout << "assigning component " << i
         //      << ", stage = " << stage << endl;
         bhi[i] = rhshi[stage][i];
         blo[i] = rhslo[stage][i];
         // cout << "b[" << i << "] : "
         //      << bhi[i] << "  " << blo[i] << endl;
      }
      CPU_dbl_onenorm(nrows,bhi,&nrm);
      if(vrblvl > 0)
         cout << scientific << setprecision(3)
              << "1-norm of b[" << stage << "] : " << nrm << endl;

      if((nrm < 1.0e-28) || (nrm > prevnorm))
      {
         skipbscnt = skipbscnt + 1;

         if(vrblvl > 0)
            cout << "-> skip backsubstitution for x[" << stage << "] ..."
                 << endl;

         for(int i=0; i<ncols; i++)
         {
            xhi[i] = 0.0;
            xlo[i] = 0.0;
         }
      }
      else
      {
         // prevnorm = 1.0e+8; // nrm*1.0e+8;

         if(vrblvl > 0)
            cout << "-> run backsubstitution for x[" << stage << "] ..."
                 << endl;

         if(firstbs)
         {
            *newtail = stage;
            firstbs = false;
         }
         double bstimelapsed_d;
         double elapsedms,invlapsed,mullapsed,sublapsed;
         long long int bsaddcnt = 0;
         long long int bsmulcnt = 0;
         long long int bsdivcnt = 0;

         if(vrblvl > 0)
            cout << "-> GPU multiplies rhs with Q^T ..." << endl;

         GPU_dbl2_bals_qtb
            (ncols,szt,nbt,Qhi,Qlo,bhi,blo,totqtblapsedms,vrblvl);

         if(vrblvl > 1)
         {
            for(int i=0; i<nrows; i++)
               cout << "Qtb[" << i << "] : "
                    << bhi[i] << "  " << blo[i] << endl;
         }
         if(vrblvl > 0)
         {
            cout << "-> GPU solves an upper triangular system ..." << endl;
 
            if(vrblvl > 1)
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

         *totbslapsedms += elapsedms;

         if(vrblvl > 0)
            write_dbl2_bstimeflops
               (szt,nbt,0,invlapsed,mullapsed,sublapsed,elapsedms,
                bstimelapsed_d,bsaddcnt,bsmulcnt,bsdivcnt);

        if(vrblvl > 1)
           for(int i=0; i<ncols; i++)
              cout << "x[" << i << "] : "
                   << xhi[i] << "  " << xlo[i] << endl;
      }
      for(int j=0; j<ncols; j++)
      {
         solhi[stage][j] = xhi[j];
         sollo[stage][j] = xlo[j];
      }
   }
   if(vrblvl > 0)
      cout << "*** solve tail skipped " << skipupcnt
           << " updates and " << skipbscnt
           << " backsubstitutions ***" << endl;

   *upidx = skipupcnt;
   *bsidx = skipbscnt;

   for(int i=0; i<nrows; i++)
   {
      free(workRhi[i]); free(workRlo[i]);
   }
   free(bhi); free(xhi); free(workRhi);
   free(blo); free(xlo); free(workRlo);
}

void GPU_cmplx2_bals_solve
 ( int dim, int degp1, int szt, int nbt, int tailidx,
   double ***matrehi, double ***matrelo, double ***matimhi, double ***matimlo,
   double **Qrehi, double **Qrelo, double **Qimhi, double **Qimlo,
   double **Rrehi, double **Rrelo, double **Rimhi, double **Rimlo,
   double **rhsrehi, double **rhsrelo, double **rhsimhi, double **rhsimlo,
   double **solrehi, double **solrelo, double **solimhi, double **solimlo, 
   bool *zeroQ, bool *noqr, int *upidx, int *bsidx, int *newtail,
   double *totqrlapsedms, double *totqtblapsedms, double *totbslapsedms,
   double *totupdlapsedms, int vrblvl )
{
   const int nrows = dim;
   const int ncols = dim;
   int skipupcnt = 0; // counts the skipped updates
   int skipbscnt = 0; // counts the skipped backsubstitutions
   double prevnorm = 1.0e+99;
   bool firstbs = true; // at the first back substitution
   *newtail = degp1; // in case no back substitution happens

   double *brehi = new double[nrows];
   double *brelo = new double[nrows];
   double *bimhi = new double[nrows];
   double *bimlo = new double[nrows];
   double *xrehi = new double[ncols];
   double *xrelo = new double[ncols];
   double *ximhi = new double[ncols];
   double *ximlo = new double[ncols];

   double **workRrehi = new double*[nrows]; // GPU upper solver changes R
   double **workRrelo = new double*[nrows];
   double **workRimhi = new double*[nrows]; 
   double **workRimlo = new double*[nrows]; 

   for(int i=0; i<nrows; i++)
   {
      workRrehi[i] = new double[ncols]; workRrelo[i] = new double[ncols];
      workRimhi[i] = new double[ncols]; workRimlo[i] = new double[ncols];
   }
   if(vrblvl > 1)
   {
      cout << "GPU_cmplx2_bals_solve blocks of rhs :" << endl;
      for(int k=0; k<degp1; k++)
      {
         for(int i=0; i<nrows; i++)
            cout << "rhs[" << k << "][" << i << "] : "
                 << rhsrehi[k][i] << "  " << rhsrelo[k][i] << endl << "  "
                 << rhsimhi[k][i] << "  " << rhsimlo[k][i] << endl;
      }
   }
   double nrm;

   if((*noqr) && (!*zeroQ))
   {
      if(vrblvl > 0) cout << "-> skipping GPU_cmplx2_bals_head ..." << endl;
   }
   else
   {
      CPU_cmplx_onenorm(nrows,rhsrehi[0],rhsimhi[0],&nrm);
      if(vrblvl > 0)
         cout << scientific << setprecision(3)
              << "1-norm of b : " << nrm << endl;

      if(nrm < 1.0e-28)
      {
         if(*zeroQ)
         {
            if(vrblvl > 0)
               cout << "-> no skipping GPU_cmplx2_bals_head because zeroQ"
                    << endl;

            *noqr = false;
         }
         else
         {
            if(vrblvl > 0)
               cout << "-> skip call to GPU_cmplx2_bals_head ..." << endl;

            *noqr = true;

            for(int j=0; j<ncols; j++)
            {
               solrehi[0][j] = 0.0; solrelo[0][j] = 0.0;
               solimhi[0][j] = 0.0; solimlo[0][j] = 0.0;
            }
         }
      }
      if(!*noqr)
      {
         if(vrblvl > 0) cout << "-> calling GPU_cmplx2_bals_head ..." << endl;

         double **Arehi = new double*[nrows];
         double **Arelo = new double*[nrows];
         double **Aimhi = new double*[nrows];
         double **Aimlo = new double*[nrows];

         for(int i=0; i<nrows; i++)
         {
            Arehi[i] = new double[ncols]; Arelo[i] = new double[ncols];
            Aimhi[i] = new double[ncols]; Aimlo[i] = new double[ncols];

            for(int j=0; j<ncols; j++)
            {
               Arehi[i][j] = matrehi[0][i][j]; Arelo[i][j] = matrelo[0][i][j];
               Aimhi[i][j] = matimhi[0][i][j]; Aimlo[i][j] = matimlo[0][i][j];
            }
            brehi[i] = rhsrehi[0][i]; brelo[i] = rhsrelo[0][i];
            bimhi[i] = rhsimhi[0][i]; bimlo[i] = rhsimlo[0][i];

            for(int j=0; j<ncols; j++)
            {
               Rrehi[i][j] = matrehi[0][i][j]; Rrelo[i][j] = matrelo[0][i][j];
               Rimhi[i][j] = matimhi[0][i][j]; Rimlo[i][j] = matimlo[0][i][j];
            }
         }
         GPU_cmplx2_bals_head
            (nrows,ncols,szt,nbt,Arehi,Arelo,Aimhi,Aimlo,
             Qrehi,Qrelo,Qimhi,Qimlo,Rrehi,Rrelo,Rimhi,Rimlo,
             brehi,brelo,bimhi,bimlo,xrehi,xrelo,ximhi,ximlo,
             totqrlapsedms,totqtblapsedms,totbslapsedms,vrblvl);

         *zeroQ = false;

         for(int j=0; j<ncols; j++)
         {
            solrehi[0][j] = xrehi[j]; solrelo[0][j] = xrelo[j];
            solimhi[0][j] = ximhi[j]; solimlo[0][j] = ximlo[j];
         }
         for(int i=0; i<nrows; i++)
         {
            free(Arehi[i]); free(Arelo[i]);
            free(Aimhi[i]); free(Aimlo[i]);
         }
         free(Arehi); free(Arelo);
         free(Aimhi); free(Aimlo);
      }
   }
   for(int stage=tailidx; stage<degp1; stage++)
   {
      if(vrblvl > 0)
         cout << "stage " << stage << " in solve tail ..." << endl;

      double *xrs = solrehi[stage-1];  // solution to do the update with
      double *xis = solimhi[stage-1];

      CPU_cmplx_onenorm(dim,xrs,xis,&nrm);
      if(vrblvl > 0)
         cout << scientific << setprecision(3)
              << "1-norm of x[" << stage-1 << "] : " << nrm << endl;

      if(nrm < 1.0e-28)
      {
         skipupcnt = skipupcnt + 1;

         if(vrblvl > 0)
            cout << "-> skip update with x[" << stage-1 << "] ..." << endl;
      }
      else
      {
         if(vrblvl > 0)
            cout << "-> updating with x[" << stage-1 << "] ..." << endl;

         GPU_cmplx2_bals_tail
            (nrows,ncols,szt,nbt,degp1,stage,
             matrehi,matrelo,matimhi,matimlo,
             rhsrehi,rhsrelo,rhsimhi,rhsimlo,
             solrehi,solrelo,solimhi,solimlo,totupdlapsedms,vrblvl);

         if(vrblvl > 1)
         {
            cout << "blocks of rhs before assignment :" << endl;
            for(int k=0; k<degp1; k++)
            {
               for(int i=0; i<nrows; i++)
                  cout << "rhs[" << k << "][" << i << "] : "
                       << rhsrehi[k][i] << "  " << rhsrelo[k][i] << endl
                       << "  " 
                       << rhsimhi[k][i] << "  " << rhsimlo[k][i] << endl;
            }
         }
      }
      for(int i=0; i<nrows; i++) 
      {
         // cout << "assigning component " << i
         //      << ", stage = " << stage << endl;
         brehi[i] = rhsrehi[stage][i];
         brelo[i] = rhsrelo[stage][i];
         bimhi[i] = rhsimhi[stage][i];
         bimlo[i] = rhsimlo[stage][i];
         // cout << "b[" << i << "] : "
         //      << brehi[i] << "  " << brelo[i] << endl << "  "
         //      << bimhi[i] << "  " << bimlo[i] << endl;
      }
      CPU_cmplx_onenorm(dim,brehi,bimhi,&nrm);
      if(vrblvl > 0)
         cout << scientific << setprecision(3)
              << "1-norm of b[" << stage << "] : " << nrm << endl;

      if((nrm < 1.0e-28) || (nrm > prevnorm))
      {
         skipbscnt = skipbscnt + 1;

         if(vrblvl > 0)
            cout << "-> skip backsubstitutions for x[" << stage << "] ..."
                 << endl;

         for(int i=0; i<ncols; i++)
         {
            xrehi[i] = 0.0; xrelo[i] = 0.0;
            ximhi[i] = 0.0; ximlo[i] = 0.0;
         }
      }
      else
      {
         // prevnorm = 1.0e+8; // nrm*1.0e+8;

         if(vrblvl > 0)
            cout << "-> run backsubstitutions for x[" << stage << "] ..."
                 << endl;

         if(firstbs)
         {
            *newtail = stage;
            firstbs = false;
         }

         double bstimelapsed_d;
         double elapsedms,invlapsed,mullapsed,sublapsed;
         long long int bsaddcnt = 0;
         long long int bsmulcnt = 0;
         long long int bsdivcnt = 0;

         if(vrblvl > 0)
            cout << "-> GPU multiplies rhs with Q^H ..." << endl;

         GPU_cmplx2_bals_qhb
            (ncols,szt,nbt,Qrehi,Qrelo,Qimhi,Qimlo,
                           brehi,brelo,bimhi,bimlo,totqtblapsedms,vrblvl);

         if(vrblvl > 1)
         {
            for(int i=0; i<nrows; i++)
               cout << "QHb[" << i << "] : "
                    << brehi[i] << "  " << brelo[i] << endl << "  "
                    << bimhi[i] << "  " << bimlo[i] << endl;
         }
         if(vrblvl > 0)
         {
            cout << "-> GPU solves an upper triangular system ..." << endl;
    
            if(vrblvl > 1)
               for(int i=0; i<nrows; i++)
                  for(int j=0; j<ncols; j++)
                     cout << "R[" << i << "][" << j << "] : "
                          << Rrehi[i][j] << "  " << Rrelo[i][j] << endl
                          << "  "
                          << Rimhi[i][j] << "  " << Rimlo[i][j] << endl;
         }
         for(int i=0; i<nrows; i++)
            for(int j=0; j<ncols; j++)
            {
               workRrehi[i][j] = Rrehi[i][j]; workRrelo[i][j] = Rrelo[i][j];
               workRimhi[i][j] = Rimhi[i][j]; workRimlo[i][j] = Rimlo[i][j];
            }

         GPU_cmplx2_upper_tiled_solver
            (ncols,szt,nbt,workRrehi,workRrelo,workRimhi,workRimlo,
             brehi,brelo,bimhi,bimlo,xrehi,xrelo,ximhi,ximlo,
             &invlapsed,&mullapsed,&sublapsed,&elapsedms,&bstimelapsed_d,
             &bsaddcnt,&bsmulcnt,&bsdivcnt);

         *totbslapsedms += elapsedms;

         if(vrblvl > 0)
            write_dbl2_bstimeflops
               (szt,nbt,1,invlapsed,mullapsed,sublapsed,elapsedms,
                bstimelapsed_d,bsaddcnt,bsmulcnt,bsdivcnt);

        if(vrblvl > 1)
           for(int i=0; i<ncols; i++)
              cout << "x[" << i << "] : "
                   << xrehi[i] << "  " << xrelo[i] << endl << "  "
                   << ximhi[i] << "  " << ximlo[i] << endl;
      }
      for(int j=0; j<ncols; j++)
      {
         solrehi[stage][j] = xrehi[j]; solrelo[stage][j] = xrelo[j];
         solimhi[stage][j] = ximhi[j]; solimlo[stage][j] = ximlo[j];
      }
   }
   if(vrblvl > 0)
      cout << "*** solve tail skipped " << skipupcnt
           << " updates and " << skipbscnt
           << " backsubstitutions ***" << endl;

   *upidx = skipupcnt;
   *bsidx = skipbscnt;

   for(int i=0; i<nrows; i++)
   {
      free(workRrehi[i]); free(workRrelo[i]);
      free(workRimhi[i]); free(workRimlo[i]);
   }
   free(brehi); free(xrehi); free(workRrehi);
   free(brelo); free(xrelo); free(workRrelo);
   free(bimhi); free(ximhi); free(workRimhi);
   free(bimlo); free(ximlo); free(workRimlo);
}
