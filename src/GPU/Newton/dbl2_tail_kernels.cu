// The file dbl2_tail_kernels.cu defines the functions with prototypes in
// the file dbl2_tail_kernels.h.

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
#include "dbl_bals_flopcounts.h"

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

__global__ void cmplx2_bals_tail
 ( int ncols, int szt,
   double *Arehi, double *Arelo, double *Aimhi, double *Aimlo,
   double *xrehi, double *xrelo, double *ximhi, double *ximlo, 
   double *brehi, double *brelo, double *bimhi, double *bimlo )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int idx = bdx*szt + tdx; // thread tdx updates b[idx]

   double Ajrehi;             // register for Arehi[idx][j]
   double Ajrelo;             // register for Arelo[idx][j]
   double Ajimhi;             // register for Aimhi[idx][j]
   double Ajimlo;             // register for Aimlo[idx][j]
   double xjrehi;             // register for xrehi[j]
   double xjrelo;             // register for xrelo[j]
   double xjimhi;             // register for ximhi[j]
   double xjimlo;             // register for ximlo[j]
   double birehi = brehi[idx];  // register for brehi[idx]
   double birelo = brelo[idx];  // register for brelo[idx]
   double biimhi = bimhi[idx];  // register for bimhi[idx]
   double biimlo = bimlo[idx];  // register for bimlo[idx]
   double acchi,acclo;

   int offset = idx*ncols;

   for(int j=0; j<ncols; j++)
   {
      Ajrehi = Arehi[offset+j];
      Ajrelo = Arelo[offset+j];
      Ajimhi = Aimhi[offset+j];
      Ajimlo = Aimlo[offset+j];
      xjrehi = xrehi[j];
      xjrelo = xrelo[j];
      xjimhi = ximhi[j];
      xjimlo = ximlo[j];
      // bi = bi - Aj*xj;
      // zre = Ajre*xjre - Ajim*xjim;
      // bire = bire - zre;
      ddg_mul(Ajrehi,Ajrelo,xjrehi,xjrelo,&acchi,&acclo);
      ddg_dec(&birehi,&birelo,acchi,acclo);
      ddg_mul(Ajimhi,Ajimlo,xjimhi,xjimlo,&acchi,&acclo);
      ddg_inc(&birehi,&birelo,acchi,acclo);
      // zim = Ajre*xjim + Ajim*xjre;
      // biim = biim - zim;
      ddg_mul(Ajrehi,Ajrelo,xjimhi,xjimlo,&acchi,&acclo);
      ddg_dec(&biimhi,&biimlo,acchi,acclo);
      ddg_mul(Ajimhi,Ajimlo,xjrehi,xjrelo,&acchi,&acclo);
      ddg_dec(&biimhi,&biimlo,acchi,acclo);
   }
   brehi[idx] = birehi;
   brelo[idx] = birelo;
   bimhi[idx] = biimhi;
   bimlo[idx] = biimlo;
}

void write_dbl2_balsflops ( int ctype, int ncols, float lapsms )
{
   cout << fixed << setprecision(3);
   cout << "Time spent for b = b - A*x: " << lapsms
        << " milliseconds." << endl;

   long long int flopcnt;
   if(ctype == 0)
      flopcnt = 20*ncols*ncols + 23*ncols*ncols;
      // as many + as * in one inner product
   else
      flopcnt = 4*20*ncols*ncols + 4*23*ncols*ncols;
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

void GPU_dbl2_bals_tail
 ( int nrows, int ncols, int szt, int nbt, int degp1, int stage,
   double ***mathi, double ***matlo, double **rhshi, double **rhslo,
   double **solhi, double **sollo, double *totupdlapsedms, int vrblvl )
{
   if(vrblvl > 1)
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
      if(vrblvl > 1)
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

      if(vrblvl > 1)
         cout << "nbt = " << nbt << ", szt = " << szt
              << ", ncols = " << ncols << endl;

      cudaEvent_t start,stop;       // to measure time spent by kernels 
      cudaEventCreate(&start);
      cudaEventCreate(&stop);
      float milliseconds;

      cudaEventRecord(start);
      dbl2_bals_tail<<<nbt,szt>>>
          (ncols,szt,Ahi_d,Alo_d,xhi_d,xlo_d,bhi_d,blo_d);
      cudaEventRecord(stop);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&milliseconds,start,stop);

      *totupdlapsedms += milliseconds;

      if(vrblvl > 0) write_dbl2_balsflops(0,ncols,milliseconds);
      
      if(vrblvl > 1)
         cout << "copying block " << k << " of right hand side ..." << endl;

      cudaMemcpy(rhshi[k],bhi_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(rhslo[k],blo_d,szrhs,cudaMemcpyDeviceToHost);
   }
   free(Ahi_h); free(Alo_h);

   if(vrblvl > 1)
   {
      cout << "GPU_dbl2_bals_tail copied blocks of rhs :" << endl;
      for(int k=0; k<degp1; k++)
      {
         for(int i=0; i<nrows; i++)
            cout << "rhs[" << k << "][" << i << "] : "
                 << rhshi[k][i] << "  " << rhslo[k][i] << endl;
      }
   }
   cudaFree(bhi_d); cudaFree(blo_d);
   cudaFree(xhi_d); cudaFree(xlo_d);
   cudaFree(Ahi_d); cudaFree(Alo_d);
}

void GPU_cmplx2_bals_tail
 ( int nrows, int ncols, int szt, int nbt, int degp1, int stage,
   double ***matrehi, double ***matrelo, double ***matimhi, double ***matimlo,
   double **rhsrehi, double **rhsrelo, double **rhsimhi, double **rhsimlo,
   double **solrehi, double **solrelo, double **solimhi, double **solimlo,
   double *totupdlapsedms, int vrblvl )
{
   if(vrblvl > 1)
   {
      cout << "GPU_cmplx2_bals_tail input blocks of rhs :" << endl;
      for(int k=0; k<degp1; k++)
      {
         for(int i=0; i<nrows; i++)
            cout << "rhs[" << k << "][" << i << "] : "
                 << rhsrehi[k][i] << "  " << rhsrelo[k][i] << endl << "  "
                 << rhsimhi[k][i] << "  " << rhsimlo[k][i] << endl;
      }
   }
   double *brehi_d;
   double *brelo_d;
   double *bimhi_d;
   double *bimlo_d;
   const size_t szrhs = nrows*sizeof(double);
   cudaMalloc((void**)&brehi_d,szrhs);
   cudaMalloc((void**)&brelo_d,szrhs);
   cudaMalloc((void**)&bimhi_d,szrhs);
   cudaMalloc((void**)&bimlo_d,szrhs);

   double *xrehi_d;
   double *xrelo_d;
   double *ximhi_d;
   double *ximlo_d;
   const size_t szsol = ncols*sizeof(double);
   cudaMalloc((void**)&xrehi_d,szsol);
   cudaMalloc((void**)&xrelo_d,szsol);
   cudaMalloc((void**)&ximhi_d,szsol);
   cudaMalloc((void**)&ximlo_d,szsol);
   cudaMemcpy(xrehi_d,solrehi[stage-1],szsol,cudaMemcpyHostToDevice);
   cudaMemcpy(xrelo_d,solrelo[stage-1],szsol,cudaMemcpyHostToDevice);
   cudaMemcpy(ximhi_d,solimhi[stage-1],szsol,cudaMemcpyHostToDevice);
   cudaMemcpy(ximlo_d,solimlo[stage-1],szsol,cudaMemcpyHostToDevice);

   double *Arehi_d;
   double *Arelo_d;
   double *Aimhi_d;
   double *Aimlo_d;
   const size_t szmat = nrows*ncols*sizeof(double);
   cudaMalloc((void**)&Arehi_d,szmat);
   cudaMalloc((void**)&Arelo_d,szmat);
   cudaMalloc((void**)&Aimhi_d,szmat);
   cudaMalloc((void**)&Aimlo_d,szmat);

   double *Arehi_h = new double[szmat];
   double *Arelo_h = new double[szmat];
   double *Aimhi_h = new double[szmat];
   double *Aimlo_h = new double[szmat];

   for(int k=stage; k<degp1; k++)
   {
      if(vrblvl > 1)
         cout << "GPU_cmplx2_bals_tail launches " << nbt
              << " thread blocks in step " << k-stage << endl;

      int idx=0;
      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++)
         {
            Arehi_h[idx]   = matrehi[k-stage+1][i][j];
            Arelo_h[idx]   = matrelo[k-stage+1][i][j];
            Aimhi_h[idx]   = matimhi[k-stage+1][i][j];
            Aimlo_h[idx++] = matimlo[k-stage+1][i][j];
         }
      
      cudaMemcpy(brehi_d,rhsrehi[k],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(brelo_d,rhsrelo[k],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(bimhi_d,rhsimhi[k],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(bimlo_d,rhsimlo[k],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(Arehi_d,Arehi_h,szmat,cudaMemcpyHostToDevice);
      cudaMemcpy(Arelo_d,Arelo_h,szmat,cudaMemcpyHostToDevice);
      cudaMemcpy(Aimhi_d,Aimhi_h,szmat,cudaMemcpyHostToDevice);
      cudaMemcpy(Aimlo_d,Aimlo_h,szmat,cudaMemcpyHostToDevice);

      if(vrblvl > 1)
         cout << "nbt = " << nbt << ", szt = " << szt
              << ", ncols = " << ncols << endl;

      cudaEvent_t start,stop;       // to measure time spent by kernels 
      cudaEventCreate(&start);
      cudaEventCreate(&stop);
      float milliseconds;

      cudaEventRecord(start);
      cmplx2_bals_tail<<<nbt,szt>>>
         (ncols,szt,Arehi_d,Arelo_d,Aimhi_d,Aimlo_d,
          xrehi_d,xrelo_d,ximhi_d,ximlo_d,brehi_d,brelo_d,bimhi_d,bimlo_d);
      cudaEventRecord(stop);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&milliseconds,start,stop);

      *totupdlapsedms += milliseconds;

      if(vrblvl > 0) write_dbl2_balsflops(1,ncols,milliseconds);
      
      if(vrblvl > 1)
         cout << "copying block " << k << " of right hand side ..." << endl;

      cudaMemcpy(rhsrehi[k],brehi_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(rhsrelo[k],brelo_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(rhsimhi[k],bimhi_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(rhsimlo[k],bimlo_d,szrhs,cudaMemcpyDeviceToHost);
   }
   free(Arehi_h); free(Aimhi_h); free(Arelo_h); free(Aimlo_h);

   if(vrblvl > 1)
   {
      cout << "GPU_cmplx2_bals_tail copied blocks of rhs :" << endl;
      for(int k=0; k<degp1; k++)
      {
         for(int i=0; i<nrows; i++)
            cout << "rhs[" << k << "][" << i << "] : "
                 << rhsrehi[k][i] << "  " << rhsrelo[k][i] << endl << "  "
                 << rhsimhi[k][i] << "  " << rhsimlo[k][i] << endl;
      }
   }
   cudaFree(brehi_d); cudaFree(brelo_d); cudaFree(bimhi_d); cudaFree(bimlo_d);
   cudaFree(xrehi_d); cudaFree(xrelo_d); cudaFree(ximhi_d); cudaFree(ximlo_d);
   cudaFree(Arehi_d); cudaFree(Arelo_d); cudaFree(Aimhi_d); cudaFree(Aimlo_d);
}

void GPU_dbl2_linear_residue
 ( int dim, int degp1, int szt, int nbt, int tailidx,
   double ***mathi, double ***matlo,
   double **rhshi, double **rhslo, double **solhi, double **sollo,
   double **resvechi, double **resveclo, double *resmaxhi, double *resmaxlo,
   double *lapms, long long int *add, long long int *mul, int vrblvl )
{
   double *rhi_d;
   double *rlo_d;
   const size_t szrhs = dim*sizeof(double);
   cudaMalloc((void**)&rhi_d,szrhs);
   cudaMalloc((void**)&rlo_d,szrhs);

   double *xhi_d;
   double *xlo_d;
   const size_t szsol = dim*sizeof(double);
   cudaMalloc((void**)&xhi_d,szsol);
   cudaMalloc((void**)&xlo_d,szsol);

   double *Ahi_d;
   double *Alo_d;
   const size_t szmat = dim*dim*sizeof(double);
   cudaMalloc((void**)&Ahi_d,szmat);
   cudaMalloc((void**)&Alo_d,szmat);

   double *Ahi_h = new double[dim*dim];
   double *Alo_h = new double[dim*dim];

   *add = 0; // initialize number of additions
   *mul = 0; // initialize number of multiplications

   if(vrblvl > 0)
      cout << "GPU_dbl2_linear_residue for deg+1 : " << degp1 << endl;

   for(int i=tailidx; i<degp1; i++)  // compute i-th residual vector
   {
      cudaMemcpy(rhi_d,rhshi[i],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(rlo_d,rhslo[i],szrhs,cudaMemcpyHostToDevice);

      for(int j=0; j<=(i-tailidx); j++)  // multiply mat[j] with sol[i-j]
      {
         int idx=0;
         for(int i1=0; i1<dim; i1++)
            for(int j1=0; j1<dim; j1++)
            {
               Ahi_h[idx]   = mathi[j][i1][j1];
               Alo_h[idx++] = matlo[j][i1][j1];
            }
      
         cudaMemcpy(Ahi_d,Ahi_h,szmat,cudaMemcpyHostToDevice);
         cudaMemcpy(Alo_d,Alo_h,szmat,cudaMemcpyHostToDevice);
         cudaMemcpy(xhi_d,solhi[i-j],szsol,cudaMemcpyHostToDevice);
         cudaMemcpy(xlo_d,sollo[i-j],szsol,cudaMemcpyHostToDevice);
         cudaMemcpy(Ahi_d,Ahi_h,szmat,cudaMemcpyHostToDevice);
         cudaMemcpy(Alo_d,Alo_h,szmat,cudaMemcpyHostToDevice);

         if(vrblvl > 1)
            cout << "GPU_dbl2_linear_residue launches " << nbt
                 << " thread blocks in step " << i << ", " << j << endl;

         cudaEvent_t start,stop;       // to measure time spent by kernels 
         cudaEventCreate(&start);
         cudaEventCreate(&stop);
         float milliseconds;

         cudaEventRecord(start);
         dbl2_bals_tail<<<nbt,szt>>>
            (dim,szt,Ahi_d,Alo_d,xhi_d,xlo_d,rhi_d,rlo_d);
         cudaEventRecord(stop);
         cudaEventSynchronize(stop);
         cudaEventElapsedTime(&milliseconds,start,stop);
         *lapms += milliseconds;
         flopcount_dbl_bals_tail(dim,add,mul);

         if(vrblvl > 0) write_dbl2_balsflops(0,dim,milliseconds);
      }
      cudaMemcpy(resvechi[i],rhi_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(resveclo[i],rlo_d,szrhs,cudaMemcpyDeviceToHost);
   }
   if(vrblvl > 1)
   {
      for(int i=tailidx; i<degp1; i++) 
      {
         cout << "Solution vector " << i << " :" << endl;
         for(int j=0; j<dim; j++)
            cout << solhi[i][j] << "  " << sollo[i][j] << endl;
         cout << "Residual vector " << i << " :" << endl;
         for(int j=0; j<dim; j++)
            cout << resvechi[i][j] << "  " << resveclo[i][j] << endl;
      }
   }
   *resmaxhi = 0.0; *resmaxlo = 0.0;

   for(int i=tailidx; i<degp1; i++)
   {
      double *rihi = resvechi[i];
      double *rilo = resveclo[i];

      for(int j=0; j<dim; j++)
         if(abs(rihi[j]) > *resmaxhi)
         {
            *resmaxhi = abs(rihi[j]);
            *resmaxlo = abs(rilo[j]);
         }
   }
   free(Ahi_h); free(Alo_h);

   cudaFree(rhi_d); cudaFree(rlo_d);
   cudaFree(xhi_d); cudaFree(xlo_d);
   cudaFree(Ahi_d); cudaFree(Alo_d);
}

void GPU_cmplx2_linear_residue
 ( int dim, int degp1, int szt, int nbt, int tailidx,
   double ***matrehi, double ***matrelo, double ***matimhi, double ***matimlo,
   double **rhsrehi, double **rhsrelo, double **rhsimhi, double **rhsimlo, 
   double **solrehi, double **solrelo, double **solimhi, double **solimlo,
   double **resvecrehi, double **resvecrelo,
   double **resvecimhi, double **resvecimlo,
   double *resmaxhi, double *resmaxlo,
   double *lapms, long long int *add, long long int *mul, int vrblvl )
{
   double *rrehi_d;
   double *rrelo_d;
   double *rimhi_d;
   double *rimlo_d;
   const size_t szrhs = dim*sizeof(double);
   cudaMalloc((void**)&rrehi_d,szrhs);
   cudaMalloc((void**)&rrelo_d,szrhs);
   cudaMalloc((void**)&rimhi_d,szrhs);
   cudaMalloc((void**)&rimlo_d,szrhs);

   double *xrehi_d;
   double *xrelo_d;
   double *ximhi_d;
   double *ximlo_d;
   const size_t szsol = dim*sizeof(double);
   cudaMalloc((void**)&xrehi_d,szsol);
   cudaMalloc((void**)&xrelo_d,szsol);
   cudaMalloc((void**)&ximhi_d,szsol);
   cudaMalloc((void**)&ximlo_d,szsol);

   double *Arehi_d;
   double *Arelo_d;
   double *Aimhi_d;
   double *Aimlo_d;
   const size_t szmat = dim*dim*sizeof(double);
   cudaMalloc((void**)&Arehi_d,szmat);
   cudaMalloc((void**)&Arelo_d,szmat);
   cudaMalloc((void**)&Aimhi_d,szmat);
   cudaMalloc((void**)&Aimlo_d,szmat);

   double *Arehi_h = new double[dim*dim];
   double *Arelo_h = new double[dim*dim];
   double *Aimhi_h = new double[dim*dim];
   double *Aimlo_h = new double[dim*dim];

   *add = 0; // initialize number of additions
   *mul = 0; // initialize number of multiplications

   if(vrblvl > 0)
      cout << "GPU_cmplx2_linear_residue for deg+1 : " << degp1 << endl;

   for(int i=tailidx; i<degp1; i++)  // compute i-th residual vector
   {
      cudaMemcpy(rrehi_d,rhsrehi[i],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(rrelo_d,rhsrelo[i],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(rimhi_d,rhsimhi[i],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(rimlo_d,rhsimlo[i],szrhs,cudaMemcpyHostToDevice);

      for(int j=0; j<=(i-tailidx); j++)  // multiply mat[j] with sol[i-j]
      {
         int idx=0;
         for(int i1=0; i1<dim; i1++)
            for(int j1=0; j1<dim; j1++)
            {
               Arehi_h[idx]   = matrehi[j][i1][j1];
               Arelo_h[idx]   = matrelo[j][i1][j1];
               Aimhi_h[idx]   = matimhi[j][i1][j1];
               Aimlo_h[idx++] = matimlo[j][i1][j1];
            }
      
         cudaMemcpy(Arehi_d,Arehi_h,szmat,cudaMemcpyHostToDevice);
         cudaMemcpy(Arelo_d,Arelo_h,szmat,cudaMemcpyHostToDevice);
         cudaMemcpy(Aimhi_d,Aimhi_h,szmat,cudaMemcpyHostToDevice);
         cudaMemcpy(Aimlo_d,Aimlo_h,szmat,cudaMemcpyHostToDevice);
         cudaMemcpy(xrehi_d,solrehi[i-j],szsol,cudaMemcpyHostToDevice);
         cudaMemcpy(xrelo_d,solrelo[i-j],szsol,cudaMemcpyHostToDevice);
         cudaMemcpy(ximhi_d,solimhi[i-j],szsol,cudaMemcpyHostToDevice);
         cudaMemcpy(ximlo_d,solimlo[i-j],szsol,cudaMemcpyHostToDevice);
         cudaMemcpy(Arehi_d,Arehi_h,szmat,cudaMemcpyHostToDevice);
         cudaMemcpy(Arelo_d,Arelo_h,szmat,cudaMemcpyHostToDevice);
         cudaMemcpy(Aimhi_d,Aimhi_h,szmat,cudaMemcpyHostToDevice);
         cudaMemcpy(Aimlo_d,Aimlo_h,szmat,cudaMemcpyHostToDevice);

         if(vrblvl > 1)
            cout << "GPU_cmplx2_linear_residue launches " << nbt
                 << " thread blocks in step " << i << ", " << j << endl;

         cudaEvent_t start,stop;       // to measure time spent by kernels 
         cudaEventCreate(&start);
         cudaEventCreate(&stop);
         float milliseconds;

         cudaEventRecord(start);
         cmplx2_bals_tail<<<nbt,szt>>>
            (dim,szt,Arehi_d,Arelo_d,Aimhi_d,Aimlo_d,
                     xrehi_d,xrelo_d,ximhi_d,ximlo_d,
                     rrehi_d,rrelo_d,rimhi_d,rimlo_d);
         cudaEventRecord(stop);
         cudaEventSynchronize(stop);
         cudaEventElapsedTime(&milliseconds,start,stop);
         *lapms += milliseconds;
         flopcount_cmplx_bals_tail(dim,add,mul);

         if(vrblvl > 0) write_dbl2_balsflops(1,dim,milliseconds);
      }
      cudaMemcpy(resvecrehi[i],rrehi_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(resvecrelo[i],rrelo_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(resvecimhi[i],rimhi_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(resvecimlo[i],rimlo_d,szrhs,cudaMemcpyDeviceToHost);
   }
   if(vrblvl > 1)
   {
      for(int i=tailidx; i<degp1; i++)
      {
         cout << "Solution vector " << i << " :" << endl;
         for(int j=0; j<dim; j++)
            cout << solrehi[i][j] << "  " << solrelo[i][j] << endl << "  "
                 << solimhi[i][j] << "  " << solimlo[i][j] << endl;

         cout << "Residual vector " << i << " :" << endl;
         for(int j=0; j<dim; j++)
            cout << resvecrehi[i][j] << "  "
                 << resvecrelo[i][j] << endl << "  "
                 << resvecimhi[i][j] << "  "
                 << resvecimlo[i][j] << endl;
      }
   }
   *resmaxhi = 0.0; *resmaxlo = 0.0;

   for(int i=tailidx; i<degp1; i++)
   {
      double *rirehi = resvecrehi[i];
      double *rirelo = resvecrelo[i];
      double *riimhi = resvecimhi[i];
      double *riimlo = resvecimlo[i];

      for(int j=0; j<dim; j++)
         if(abs(rirehi[j]) + abs(riimhi[j]) > *resmaxhi)
         {
            *resmaxhi = abs(rirehi[j]) + abs(riimhi[j]);
            *resmaxlo = abs(rirelo[j]) + abs(riimlo[j]);
         }
   }
   free(Arehi_h); free(Arelo_h); free(Aimhi_h); free(Aimlo_h);

   cudaFree(rrehi_d); cudaFree(rrelo_d);
   cudaFree(rimhi_d); cudaFree(rimlo_d);
   cudaFree(xrehi_d); cudaFree(xrelo_d);
   cudaFree(ximhi_d); cudaFree(ximlo_d);
   cudaFree(Arehi_d); cudaFree(Arelo_d);
   cudaFree(Aimhi_d); cudaFree(Aimlo_d);
}
