// The file dbl_bals_kernels.cu defines the functions with prototypes in
// the file dbl_bals_kernels.h.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#include "dbl_baqr_kernels.h"
#include "dbl_tabs_kernels.h"
#include "dbl_tail_kernels.h"
#include "dbl_bals_kernels.h"
#include "write_dbl_qrtimeflops.h"
#include "write_dbl_bstimeflops.h"
#include "dbl_onenorms_host.h"

using namespace std;

__global__ void dbl_bals_qtb
 ( int ncols, int szt, double *Qt, double *b, double *r )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int idx = bdx*szt + tdx; // thread tdx computes r[idx]

   double Qj;           // register for Q^T[idx][j]
   double bj;           // register for b[j]
   double ri = 0.0;     // register for result, r[idx]

   int offset = idx*ncols;

   for(int j=0; j<ncols; j++)
   {
      Qj = Qt[offset+j];
      bj = b[j];
      ri = ri + Qj*bj;
   }
   r[idx] = ri;
}

__global__ void cmplx_bals_qhb
 ( int ncols, int szt, double *QHre, double *QHim,
   double *bre, double *bim, double *rre, double *rim )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int idx = bdx*szt + tdx; // thread tdx computes r[idx]

   double Qjre;           // register for real part of Q^H[idx][j]
   double Qjim;           // register for imaginary part of Q^H[idx][j]
   double bjre;           // register for bre[j]
   double bjim;           // register for bim[j]
   double rire = 0.0;     // register for result, rre[idx]
   double riim = 0.0;     // register for result, rim[idx]
   double zre,zim;

   int offset = idx*ncols;

   for(int j=0; j<ncols; j++)
   {
      Qjre = QHre[offset+j];
      Qjim = QHim[offset+j];
      bjre = bre[j];
      bjim = bim[j];
      // ri = ri + Qj*bj;
      zre = Qjre*bjre - Qjim*bjim;
      zim = Qjre*bjim + Qjim*bjre;
      rire = rire + zre;
      riim = riim + zim;
   }
   rre[idx] = rire;
   rim[idx] = riim;
}

void GPU_dbl_bals_head
 ( int nrows, int ncols, int szt, int nbt,
   double **A, double **Q, double **R, double *b, double *x,
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

   GPU_dbl_blocked_houseqr
      (nrows,ncols,szt,nbt,A,Q,R,
       &houselapsedms,&RTvlapsedms,&tileRlapsedms,&vb2Wlapsedms,
       &WYTlapsedms,&QWYTlapsedms,&Qaddlapsedms,
       &YWTlapsedms,&YWTClapsedms,&Raddlapsedms,&qrtimelapsed_d,
       &qraddcnt,&qrmulcnt,&qrdivcnt,&sqrtcnt,verbose);

   *totqrlapsedms = *totqrlapsedms + houselapsedms + RTvlapsedms
      + tileRlapsedms + vb2Wlapsedms + WYTlapsedms + QWYTlapsedms
      + Qaddlapsedms + YWTlapsedms + YWTClapsedms + Raddlapsedms;

   if(vrblvl > 0)
      write_dbl_qrtimeflops
         (0,nrows,ncols,
          houselapsedms,RTvlapsedms,tileRlapsedms,vb2Wlapsedms,WYTlapsedms,
          QWYTlapsedms,Qaddlapsedms,YWTlapsedms,YWTClapsedms,Raddlapsedms,
          qrtimelapsed_d,qraddcnt,qrmulcnt,qrdivcnt,sqrtcnt);

   if(vrblvl > 0) cout << "-> GPU multiplies rhs with Q^T ..." << endl;

   GPU_dbl_bals_qtb(ncols,szt,nbt,Q,b,totqtblapsedms,vrblvl);

   if(vrblvl > 1)
   {
      for(int i=0; i<nrows; i++)
         cout << "Qtb[" << i << "] : " << b[i] << endl;
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
            cout << "R[" << i << "][" << j << "] : " << R[i][j] << endl;

      for(int i=0; i<nrows; i++)
         cout << "b[" << i << "] : " << b[i] << endl;
   }
   double **workR = new double*[nrows]; // work around ...
   for(int i=0; i<nrows; i++)
   {
      workR[i] = new double[ncols];
      for(int j=0; j<ncols; j++) workR[i][j] = R[i][j];
   }
   GPU_dbl_upper_tiled_solver
      (ncols,szt,nbt,workR,b,x,
       &invlapsed,&mullapsed,&sublapsed,&elapsedms,&bstimelapsed_d,
       &bsaddcnt,&bsmulcnt,&bsdivcnt);

   *totbslapsedms += elapsedms;

   if(vrblvl > 0)
      write_dbl_bstimeflops
         (szt,nbt,0,invlapsed,mullapsed,sublapsed,elapsedms,bstimelapsed_d,
          bsaddcnt,bsmulcnt,bsdivcnt);

   if(vrblvl > 1)
   {
      cout << "-> after calling the GPU upper solver ..." << endl;
      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++)
            cout << "R[" << i << "][" << j << "] : " << R[i][j] << endl;

      for(int i=0; i<nrows; i++)
         cout << "b[" << i << "] : " << b[i] << endl;
      for(int i=0; i<ncols; i++)
         cout << "x[" << i << "] : " << x[i] << endl;
   }
   for(int i=0; i<nrows; i++) free(workR[i]);
   free(workR);
}

void GPU_cmplx_bals_head
 ( int nrows, int ncols, int szt, int nbt,
   double **Are, double **Aim, double **Qre, double **Qim,
   double **Rre, double **Rim, double *bre, double *bim,
   double *xre, double *xim,
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

   GPU_cmplx_blocked_houseqr
      (nrows,ncols,szt,nbt,Are,Aim,Qre,Qim,Rre,Rim,
       &houselapsedms,&RTvlapsedms,&tileRlapsedms,&vb2Wlapsedms,
       &WYTlapsedms,&QWYTlapsedms,&Qaddlapsedms,
       &YWTlapsedms,&YWTClapsedms,&Raddlapsedms,&qrtimelapsed_d,
       &qraddcnt,&qrmulcnt,&qrdivcnt,&sqrtcnt,verbose);

   *totqrlapsedms = *totqrlapsedms + houselapsedms + RTvlapsedms
      + tileRlapsedms + vb2Wlapsedms + WYTlapsedms + QWYTlapsedms
      + Qaddlapsedms + YWTlapsedms + YWTClapsedms + Raddlapsedms;

   if(vrblvl > 0)
      write_dbl_qrtimeflops
         (1,nrows,ncols,
          houselapsedms,RTvlapsedms,tileRlapsedms,vb2Wlapsedms,WYTlapsedms,
          QWYTlapsedms,Qaddlapsedms,YWTlapsedms,YWTClapsedms,Raddlapsedms,
          qrtimelapsed_d,qraddcnt,qrmulcnt,qrdivcnt,sqrtcnt);

   if(vrblvl > 0)
      cout << "-> GPU multiplies rhs with Q^H ..." << endl;

   GPU_cmplx_bals_qhb
      (ncols,szt,nbt,Qre,Qim,bre,bim,totqtblapsedms,vrblvl);

   if(vrblvl > 1)
   {
      for(int i=0; i<nrows; i++)
         cout << "Qtb[" << i << "] : "
              << bre[i] << "  " << bim[i] << endl;
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
                 << Rre[i][j] << "  " << Rim[i][j] << endl;

      for(int i=0; i<nrows; i++)
         cout << "b[" << i << "] : "
              << bre[i] << "  " << bim[i] << endl;
   }
   double **workRre = new double*[nrows]; // work around ...
   double **workRim = new double*[nrows];

   for(int i=0; i<nrows; i++)
   {
      workRre[i] = new double[ncols];
      workRim[i] = new double[ncols];

      for(int j=0; j<ncols; j++)
      {
         workRre[i][j] = Rre[i][j];
         workRim[i][j] = Rim[i][j];
      }
   }
   if(vrblvl > 1)
   {
      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++)
            cout << "workR[" << i << "][" << j << "] : "
                 << workRre[i][j] << "  " << workRim[i][j] << endl;

      cout << "  ncols = " << ncols
           << "  szt = " << szt
           << "  nbt = " << nbt << endl;
   }
   GPU_cmplx_upper_tiled_solver
      (ncols,szt,nbt,workRre,workRim,bre,bim,xre,xim,
       &invlapsed,&mullapsed,&sublapsed,&elapsedms,&bstimelapsed_d,
       &bsaddcnt,&bsmulcnt,&bsdivcnt);

   *totbslapsedms += elapsedms;

   if(vrblvl > 0)
      write_dbl_bstimeflops
         (szt,nbt,1,invlapsed,mullapsed,sublapsed,elapsedms,bstimelapsed_d,
          bsaddcnt,bsmulcnt,bsdivcnt);

   if(vrblvl > 1)
   {
      cout << "-> GPU_cmplx_bals_head, after GPU upper solver ..." << endl;
      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++)
            cout << "R[" << i << "][" << j << "] : "
                 << Rre[i][j] << "  " << Rim[i][j] << endl;

      for(int i=0; i<nrows; i++)
         cout << "b[" << i << "] : "
              << bre[i] << "  " << bim[i] << endl;
      for(int i=0; i<ncols; i++)
         cout << "x[" << i << "] : "
              << xre[i] << "  " << xim[i] << endl;
   }
   for(int i=0; i<nrows; i++)
   {
      free(workRre[i]); free(workRim[i]);
   }
   free(workRre); free(workRim);
}

void write_dbl_qtbflops ( int ctype, int ncols, float lapsms )
{
   cout << fixed << setprecision(3);
   cout << "Time spent for Q^T*b : " << lapsms << " milliseconds." << endl;

   long long int flopcnt;
   const long long int longncols2 = ncols*ncols; // to avoid overflow
   if(ctype == 0)
      flopcnt = 2*longncols2; // as many + as * in one inner product
   else
      flopcnt = 16*longncols2;
   // for complex: 2 ops for +, 6 for *, which is 8 in total

   cout << "    Total number of floating-point operations : "
        << flopcnt << endl;

   long long int bytecnt;

   if(ctype == 0)
      bytecnt = ncols*ncols;
   else
      bytecnt = 2*ncols*ncols;

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

void GPU_dbl_bals_qtb
 ( int ncols, int szt, int nbt, double **Q, double *b,
   double *totqtblapsedms, int vrblvl )
{
   double *b_d;
   const size_t szrhs = ncols*sizeof(double);
   cudaMalloc((void**)&b_d,szrhs);
   cudaMemcpy(b_d,b,szrhs,cudaMemcpyHostToDevice);

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

   cudaMemcpy(Qt_d,Qt_h,szmat,cudaMemcpyHostToDevice);

   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;

   cudaEventRecord(start);
   dbl_bals_qtb<<<nbt,szt>>>(ncols,szt,Qt_d,b_d,r_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);

   *totqtblapsedms += milliseconds;

   cudaMemcpy(b,r_d,szrhs,cudaMemcpyDeviceToHost);

   if(vrblvl > 0) write_dbl_qtbflops(0,ncols,milliseconds);

   cudaFree(b_d); cudaFree(r_d); cudaFree(Qt_d);

   free(Qt_h);
}

void GPU_cmplx_bals_qhb
 ( int ncols, int szt, int nbt, double **Qre, double **Qim,
   double *bre, double *bim, double *totqtblapsedms, int vrblvl )
{
   double *bre_d;
   double *bim_d;
   const size_t szrhs = ncols*sizeof(double);
   cudaMalloc((void**)&bre_d,szrhs);
   cudaMalloc((void**)&bim_d,szrhs);
   cudaMemcpy(bre_d,bre,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(bim_d,bim,szrhs,cudaMemcpyHostToDevice);

   double *rre_d;
   double *rim_d;
   const size_t szsol = ncols*sizeof(double);
   cudaMalloc((void**)&rre_d,szsol);
   cudaMalloc((void**)&rim_d,szsol);

   double *QHre_d;
   double *QHim_d;
   const size_t szmat = ncols*ncols*sizeof(double);
   cudaMalloc((void**)&QHre_d,szmat);
   cudaMalloc((void**)&QHim_d,szmat);

   double *QHre_h = new double[szmat];
   double *QHim_h = new double[szmat];

   int idx=0;
   for(int i=0; i<ncols; i++)
      for(int j=0; j<ncols; j++)
      {
         QHre_h[idx]   = Qre[j][i];
         QHim_h[idx++] = - Qim[j][i]; // Hermitian transpose !
      }

   cudaMemcpy(QHre_d,QHre_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(QHim_d,QHim_h,szmat,cudaMemcpyHostToDevice);

   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   float milliseconds;

   cudaEventRecord(start);
   cmplx_bals_qhb<<<nbt,szt>>>
      (ncols,szt,QHre_d,QHim_d,bre_d,bim_d,rre_d,rim_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);

   *totqtblapsedms += milliseconds;

   cudaMemcpy(bre,rre_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(bim,rim_d,szrhs,cudaMemcpyDeviceToHost);

   if(vrblvl > 0) write_dbl_qtbflops(1,ncols,milliseconds);

   cudaFree(bre_d); cudaFree(bim_d); cudaFree(rre_d); cudaFree(rim_d);
   cudaFree(QHre_d); cudaFree(QHim_d);

   free(QHre_h); free(QHim_h);
}

void GPU_dbl_bals_solve
 ( int dim, int degp1, int szt, int nbt, int tailidx,
   double ***mat, double **Q, double **R, double **rhs, double **sol,
   bool *zeroQ, bool *noqr, int *upidx, int *bsidx, int *newtail,
   double *totqrlapsedms, double *totqtblapsedms, double *totbslapsedms,
   double *totupdlapsedms, int vrblvl )
{
   const int nrows = dim;
   const int ncols = dim;
   int skipupcnt = 0; // counts skipped updates
   int skipbscnt = 0; // counts skipped back substitutions
   double prevnorm = 1.0e+99;
   bool firstbs = true; // at the first back substitution
   *newtail = degp1; // in case no back substitution happens

   double *b = new double[nrows];
   double *x = new double[ncols];

   double **workR = new double*[nrows]; // GPU upper solver changes R
   for(int i=0; i<nrows; i++) workR[i] = new double[ncols];

   if(vrblvl > 1)
   {
      cout << "GPU_dbl_bals_solve blocks of rhs :" << endl;
      for(int k=0; k<degp1; k++)
      {
         for(int i=0; i<nrows; i++)
            cout << "rhs[" << k << "][" << i << "] : " << rhs[k][i] << endl;
      }
   }
   double nrm;

   if((*noqr) && (!*zeroQ))
   {
      if(vrblvl > 0) cout << "-> skipping GPU_dbl_bals_head ..." << endl;
   }
   else
   {
      CPU_dbl_onenorm(nrows,rhs[0],&nrm);
      if(vrblvl > 0)
         cout << scientific << setprecision(3)
              << "1-norm of b : " << nrm << endl;

      if(nrm + 1.0 == 1.0)
      {
         if(*zeroQ)
         {
            if(vrblvl > 0)
               cout << "-> no skipping of GPU_dbl_bals_head because zeroQ"
                    << endl;

            *noqr = false;
         }
         else
         {
            if(vrblvl > 0)
               cout << "-> skip call to GPU_dbl_bals_head ..." << endl;
   
            *noqr = true;

            for(int j=0; j<ncols; j++) sol[0][j] = 0.0;
         }
      }
      if(!*noqr)
      {
         if(vrblvl > 0) cout << "-> calling GPU_dbl_bals_head ..." << endl;

         double **A = new double*[nrows];

         for(int i=0; i<nrows; i++)
         {
            A[i] = new double[ncols];
            for(int j=0; j<ncols; j++) A[i][j] = mat[0][i][j];
            b[i] = rhs[0][i];
            for(int j=0; j<ncols; j++) R[i][j] = mat[0][i][j];
         }
         GPU_dbl_bals_head
            (nrows,ncols,szt,nbt,A,Q,R,b,x,
             totqrlapsedms,totqtblapsedms,totbslapsedms,vrblvl);

         *zeroQ = false;

         if(vrblvl > 0)
         {
            CPU_dbl_onenorm(ncols,x,&nrm);
            cout << scientific << setprecision(3)
                 << "1-norm of x : " << nrm << endl;
         }
         for(int j=0; j<ncols; j++) sol[0][j] = x[j];

         for(int i=0; i<nrows; i++) free(A[i]);
         free(A);
      }
   }
   for(int stage=tailidx; stage<degp1; stage++)
   {
      if(vrblvl > 0)
         cout << "stage " << stage << " in solve tail ..." << endl;

      double *xs = sol[stage-1];       // solution to do the update with
      CPU_dbl_onenorm(dim,xs,&nrm);
      if(vrblvl > 0)
         cout << scientific << setprecision(3)
              << "1-norm of x[" << stage-1 << "] : " << nrm << endl;

      if(nrm < 1.0e-15)
      {
         skipupcnt = skipupcnt + 1;

         if(vrblvl > 0)
            cout << "-> skip update with x[" << stage-1 << "] ..." << endl;
      }
      else
      {
         if(vrblvl > 0)
            cout << "-> updating with x[" << stage-1 << "] ..." << endl;

         GPU_dbl_bals_tail
            (nrows,ncols,szt,nbt,degp1,stage,mat,rhs,sol,
             totupdlapsedms,vrblvl);

         if(vrblvl > 1)
         {
            cout << "blocks of rhs before assignment :" << endl;
            for(int k=0; k<degp1; k++)
            {
               for(int i=0; i<nrows; i++)
                  cout << "rhs[" << k << "][" << i << "] : "
                       << rhs[k][i] << endl;
            }
         }
      }
      for(int i=0; i<nrows; i++) 
      {
         // cout << "assigning component " << i
         //      << ", stage = " << stage << endl;
         b[i] = rhs[stage][i];
         // cout << "b[" << i << "] : " << b[i] << endl;
      }
      CPU_dbl_onenorm(nrows,b,&nrm);
      if(vrblvl > 0)
         cout << scientific << setprecision(3)
              << "1-norm of b[" << stage << "] : " << nrm << endl;

      if((nrm < 1.0e-15) || (nrm > prevnorm))
      {
         skipbscnt = skipbscnt + 1;

         if(vrblvl > 0)
            cout << "skip backsubstitution for x[" << stage << "] ..." << endl;

         for(int i=0; i<ncols; i++) x[i] = 0.0;
      }
      else
      {
         // prevnorm = 1.0e+8; // nrm*1.0e+8;

         if(vrblvl > 0)
            cout << "run backsubstitution for x[" << stage << "] ..." << endl;

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

         GPU_dbl_bals_qtb(ncols,szt,nbt,Q,b,totqtblapsedms,vrblvl);

         if(vrblvl > 1)
         {
            for(int i=0; i<nrows; i++)
               cout << "Qtb[" << i << "] : " << b[i] << endl;
         }
         if(vrblvl > 0)
         {
            cout << "-> GPU solves an upper triangular system ..." << endl;
 
            if(vrblvl > 1)
               for(int i=0; i<nrows; i++)
                  for(int j=0; j<ncols; j++)
                     cout << "R[" << i << "][" << j << "] : "
                          << R[i][j] << endl;
         }
         for(int i=0; i<nrows; i++)
            for(int j=0; j<ncols; j++) workR[i][j] = R[i][j];

         GPU_dbl_upper_tiled_solver
            (ncols,szt,nbt,workR,b,x,
             &invlapsed,&mullapsed,&sublapsed,&elapsedms,&bstimelapsed_d,
             &bsaddcnt,&bsmulcnt,&bsdivcnt);

         *totbslapsedms += elapsedms;

         if(vrblvl > 0)
            write_dbl_bstimeflops
               (szt,nbt,0,invlapsed,mullapsed,sublapsed,elapsedms,
                bstimelapsed_d,bsaddcnt,bsmulcnt,bsdivcnt);

         if(vrblvl > 1)
            for(int i=0; i<ncols; i++)
               cout << "x[" << i << "] : " << x[i] << endl;
      }
      for(int j=0; j<ncols; j++) sol[stage][j] = x[j];
   }
   if(vrblvl > 0)
      cout << "*** solve tail skipped " << skipupcnt
           << " updates and " << skipbscnt
           << " backsubstitutions ***" << endl;

   *upidx = skipupcnt;
   *bsidx = skipbscnt;

   for(int i=0; i<nrows; i++) free(workR[i]);

   free(b); free(x); free(workR);
}

void GPU_cmplx_bals_solve
 ( int dim, int degp1, int szt, int nbt, int tailidx,
   double ***matre, double ***matim, double **Qre, double **Qim,
   double **Rre, double **Rim, double **rhsre, double **rhsim,
   double **solre, double **solim,
   bool *zeroQ, bool *noqr, int *upidx, int *bsidx, int *newtail,
   double *totqrlapsedms, double *totqtblapsedms, double *totbslapsedms,
   double *totupdlapsedms, int vrblvl )
{
   const int nrows = dim;
   const int ncols = dim;
   int skipupcnt = 0; // counts skipped updates
   int skipbscnt = 0; // counts skipped backsubstitutions
   double prevnorm = 1.0e+99;
   bool firstbs = true; // at the first back substitution
   *newtail = degp1; // in case no back substitution happens

   double *bre = new double[nrows];
   double *bim = new double[nrows];
   double *xre = new double[ncols];
   double *xim = new double[ncols];

   double **workRre = new double*[nrows]; // GPU upper solver changes R
   double **workRim = new double*[nrows]; 

   for(int i=0; i<nrows; i++)
   {
      workRre[i] = new double[ncols];
      workRim[i] = new double[ncols];
   }
   if(vrblvl > 1)
   {
      cout << "GPU_cmplx_bals_solve blocks of rhs :" << endl;
      for(int k=0; k<degp1; k++)
      {
         for(int i=0; i<nrows; i++)
            cout << "rhs[" << k << "][" << i << "] : "
                 << rhsre[k][i] << "  " << rhsim[k][i] << endl;
      }
   }
   double nrm;

   if((*noqr) && (!*zeroQ))
   {
      if(vrblvl > 0) cout << "skipping GPU_cmplx_bals_head ..." << endl;
   }
   else
   {
      CPU_cmplx_onenorm(nrows,rhsre[0],rhsim[0],&nrm);
      if(vrblvl > 0)
         cout << scientific << setprecision(3)
              << "1-norm of b : " << nrm << endl;

      if(nrm + 1.0 == 1.0)
      {
         if(*zeroQ)
         {
            if(vrblvl > 0)
               cout << "no skipping of GPU_cmplx_bals_head because zeroQ"
                    << endl;

            *noqr = false;
         }
         else
         {
            if(vrblvl > 0)
               cout << "skip call to GPU_cmplx_bals_head ..." << endl;

            *noqr = true;
  
            for(int j=0; j<ncols; j++)
            {
               solre[0][j] = 0.0; solim[0][j] = 0.0;
            }
         }
      }
      if(!*noqr)
      {
         if(vrblvl > 0) cout << "calling GPU_cmplx_bals_head ..." << endl;

         double **Are = new double*[nrows];
         double **Aim = new double*[nrows];

         for(int i=0; i<nrows; i++)
         {
            Are[i] = new double[ncols]; Aim[i] = new double[ncols];
   
            for(int j=0; j<ncols; j++)
            {
               Are[i][j] = matre[0][i][j]; Aim[i][j] = matim[0][i][j];
            }
            bre[i] = rhsre[0][i]; bim[i] = rhsim[0][i];

            for(int j=0; j<ncols; j++)
            {
               Rre[i][j] = matre[0][i][j]; Rim[i][j] = matim[0][i][j];
            }
         }
         GPU_cmplx_bals_head
            (nrows,ncols,szt,nbt,Are,Aim,Qre,Qim,Rre,Rim,bre,bim,
             xre,xim,totqrlapsedms,totqtblapsedms,totbslapsedms,vrblvl);

         *zeroQ = false;

         for(int j=0; j<ncols; j++)
         {
            solre[0][j] = xre[j];
            solim[0][j] = xim[j];
         }
         for(int i=0; i<nrows; i++)
         {
            free(Are[i]); free(Aim[i]);
         }
         free(Are); free(Aim);
      }
   }
   for(int stage=tailidx; stage<degp1; stage++)
   {
      if(vrblvl > 0)
         cout << "stage " << stage << " in solve tail ..." << endl;

      double *xrs = solre[stage-1];       // solution to do the update with
      double *xis = solim[stage-1];

      CPU_cmplx_onenorm(dim,xrs,xis,&nrm);
      if(vrblvl > 0)
         cout << scientific << setprecision(3)
              << "1-norm of x[" << stage-1 << "] : " << nrm << endl;

      if(nrm < 1.0e-15)
      {
         skipupcnt = skipupcnt + 1;

         if(vrblvl > 0)
            cout << "skip update with x[" << stage-1 << "] ..." << endl;
      }
      else
      {
         if(vrblvl > 0)
            cout << "updating with x[" << stage-1 << "] ..." << endl;

         GPU_cmplx_bals_tail
            (nrows,ncols,szt,nbt,degp1,stage,matre,matim,
             rhsre,rhsim,solre,solim,totupdlapsedms,vrblvl);

         if(vrblvl > 1)
         {
            cout << "blocks of rhs before assignment :" << endl;
            for(int k=0; k<degp1; k++)
            {
               for(int i=0; i<nrows; i++)
                  cout << "rhs[" << k << "][" << i << "] : "
                       << rhsre[k][i] << "  " << rhsim[k][i] << endl;
            }
         }
      }
      for(int i=0; i<nrows; i++) 
      {
         // cout << "assigning component " << i
         //      << ", stage = " << stage << endl;
         bre[i] = rhsre[stage][i];
         bim[i] = rhsim[stage][i];
         // cout << "b[" << i << "] : "
         //      << bre[i] << "  " << bim[i] << endl;
      }
      CPU_cmplx_onenorm(dim,bre,bim,&nrm);
      if(vrblvl > 0)
         cout << scientific << setprecision(3)
              << "1-norm of b[" << stage << "] : " << nrm << endl;

      if((nrm < 1.0e-15) || (nrm > prevnorm))
      {
         skipbscnt = skipbscnt + 1;

         if(vrblvl > 0)
            cout << "-> skip backsubstitutions for x[" << stage << "] ..."
                 << endl;

         for(int i=0; i<ncols; i++)
         {
            xre[i] = 0.0;
            xim[i] = 0.0;
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

         GPU_cmplx_bals_qhb
            (ncols,szt,nbt,Qre,Qim,bre,bim,totqtblapsedms,vrblvl);

         if(vrblvl > 1)
         {
            for(int i=0; i<nrows; i++)
               cout << "QHb[" << i << "] : "
                    << bre[i] << "  " << bim[i] << endl;
         }
         if(vrblvl > 0)
         {
            cout << "-> GPU solves an upper triangular system ..." << endl;
 
            if(vrblvl > 1)
               for(int i=0; i<nrows; i++)
                  for(int j=0; j<ncols; j++)
                     cout << "R[" << i << "][" << j << "] : "
                          << Rre[i][j] << "  " << Rim[i][j] << endl;
         }
         for(int i=0; i<nrows; i++)
            for(int j=0; j<ncols; j++)
            {
               workRre[i][j] = Rre[i][j]; workRim[i][j] = Rim[i][j];
            }

         GPU_cmplx_upper_tiled_solver
            (ncols,szt,nbt,workRre,workRim,bre,bim,xre,xim,
             &invlapsed,&mullapsed,&sublapsed,&elapsedms,&bstimelapsed_d,
             &bsaddcnt,&bsmulcnt,&bsdivcnt);

         *totbslapsedms += elapsedms;

         if(vrblvl > 0)
            write_dbl_bstimeflops
               (szt,nbt,1,invlapsed,mullapsed,sublapsed,elapsedms,
                bstimelapsed_d,bsaddcnt,bsmulcnt,bsdivcnt);

         if(vrblvl > 1)
            for(int i=0; i<ncols; i++)
               cout << "x[" << i << "] : " << xre[i] << "  " << xim[i] << endl;
      }
      for(int j=0; j<ncols; j++)
      {
         solre[stage][j] = xre[j];
         solim[stage][j] = xim[j];
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
      free(workRre[i]); free(workRim[i]);
   }
   free(bre); free(xre); free(workRre);
   free(bim); free(xim); free(workRim);
}
