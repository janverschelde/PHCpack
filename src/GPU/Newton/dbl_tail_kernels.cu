// The file dbl_tail_kernels.cu defines the functions with prototypes in
// the file dbl_tail_kernels.h.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#include "dbl_tail_kernels.h"
#include "dbl_bals_flopcounts.h"

using namespace std;

__global__ void dbl_bals_tail
 ( int ncols, int szt, double *A, double *x, double *b )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int idx = bdx*szt + tdx; // thread tdx updates b[idx]

   double Aj;           // register for A[idx][j]
   double xj;           // register for x[j]
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

__global__ void cmplx_bals_tail
 ( int ncols, int szt, double *Are, double *Aim,
   double *xre, double *xim, double *bre, double *bim )
{
   const int bdx = blockIdx.x;
   const int tdx = threadIdx.x;
   const int idx = bdx*szt + tdx; // thread tdx updates b[idx]

   double Ajre;             // register for Are[idx][j]
   double Ajim;             // register for Aim[idx][j]
   double xjre;             // register for xre[j]
   double xjim;             // register for xim[j]
   double bire = bre[idx];  // register for bre[idx]
   double biim = bim[idx];  // register for bim[idx]
   double zre,zim;

   int offset = idx*ncols;

   for(int j=0; j<ncols; j++)
   {
      Ajre = Are[offset+j];
      Ajim = Aim[offset+j];
      xjre = xre[j];
      xjim = xim[j];
      // bi = bi - Aj*xj;
      zre = Ajre*xjre - Ajim*xjim;
      zim = Ajre*xjim + Ajim*xjre;
      bire = bire - zre;
      biim = biim - zim;
   }
   bre[idx] = bire;
   bim[idx] = biim;
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
   cudaMemcpy(x_d,sol[stage-1],szsol,cudaMemcpyHostToDevice);

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
         for(int j=0; j<ncols; j++) A_h[idx++] = mat[k-stage+1][i][j];
      
      cudaMemcpy(b_d,rhs[k],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(A_d,A_h,szmat,cudaMemcpyHostToDevice);

      if(verbose)
         cout << "nbt = " << nbt << ", szt = " << szt
              << ", ncols = " << ncols << endl;

      dbl_bals_tail<<<nbt,szt>>>(ncols,szt,A_d,x_d,b_d);
      
      if(verbose)
         cout << "copying block " << k << " of right hand side ..." << endl;

      cudaMemcpy(rhs[k],b_d,szrhs,cudaMemcpyDeviceToHost);
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
   cudaFree(b_d); cudaFree(x_d); cudaFree(A_d);
}

void GPU_cmplx_bals_tail
 ( int nrows, int ncols, int szt, int nbt, int degp1, int stage,
   double ***matre, double ***matim, double **rhsre, double **rhsim,
   double **solre, double **solim, bool verbose )
{
   if(verbose)
   {
      cout << "GPU_cmplx_bals_tail input blocks of rhs :" << endl;
      for(int k=0; k<degp1; k++)
      {
         for(int i=0; i<nrows; i++)
            cout << "rhs[" << k << "][" << i << "] : "
                 << rhsre[k][i] << "  " << rhsim[k][i] << endl;
      }
   }
   double *bre_d;
   double *bim_d;
   const size_t szrhs = nrows*sizeof(double);
   cudaMalloc((void**)&bre_d,szrhs);
   cudaMalloc((void**)&bim_d,szrhs);

   double *xre_d;
   double *xim_d;
   const size_t szsol = ncols*sizeof(double);
   cudaMalloc((void**)&xre_d,szsol);
   cudaMalloc((void**)&xim_d,szsol);
   cudaMemcpy(xre_d,solre[stage-1],szsol,cudaMemcpyHostToDevice);
   cudaMemcpy(xim_d,solim[stage-1],szsol,cudaMemcpyHostToDevice);

   if(verbose)
   {
      cout << "GPU_cmplx_bals_tail solution x :" << endl;
      for(int i=0; i<ncols; i++)
         cout << solre[stage-1][i] << "  " << solim[stage-1][i] << endl;
   }
   double *Are_d;
   double *Aim_d;
   const size_t szmat = nrows*ncols*sizeof(double);
   cudaMalloc((void**)&Are_d,szmat);
   cudaMalloc((void**)&Aim_d,szmat);

   double *Are_h = new double[szmat];
   double *Aim_h = new double[szmat];

   for(int k=stage; k<degp1; k++)
   {
      if(verbose)
         cout << "GPU_cmplx_bals_tail launches " << nbt
              << " thread blocks in step " << k-stage << endl;

      int idx=0;
      for(int i=0; i<nrows; i++)
         for(int j=0; j<ncols; j++)
         {
            Are_h[idx]   = matre[k-stage+1][i][j];
            Aim_h[idx++] = matim[k-stage+1][i][j];
         }
      
      cudaMemcpy(bre_d,rhsre[k],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(bim_d,rhsim[k],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(Are_d,Are_h,szmat,cudaMemcpyHostToDevice);
      cudaMemcpy(Aim_d,Aim_h,szmat,cudaMemcpyHostToDevice);

      if(verbose)
         cout << "nbt = " << nbt << ", szt = " << szt
              << ", ncols = " << ncols << endl;

      cmplx_bals_tail<<<nbt,szt>>>
         (ncols,szt,Are_d,Aim_d,xre_d,xim_d,bre_d,bim_d);
      
      if(verbose)
         cout << "copying block " << k << " of right hand side ..." << endl;

      cudaMemcpy(rhsre[k],bre_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(rhsim[k],bim_d,szrhs,cudaMemcpyDeviceToHost);
   }
   free(Are_h); free(Aim_h);

   if(verbose)
   {
      cout << "GPU_cmplx_bals_tail copied blocks of rhs :" << endl;
      for(int k=0; k<degp1; k++)
      {
         for(int i=0; i<nrows; i++)
            cout << "rhs[" << k << "][" << i << "] : "
                 << rhsre[k][i] << "  " << rhsim[k][i] << endl;
      }
   }
   cudaFree(bre_d); cudaFree(bim_d);
   cudaFree(xre_d); cudaFree(xim_d);
   cudaFree(Are_d); cudaFree(Aim_d);
}

void GPU_dbl_linear_residue
 ( int dim, int degp1, int szt, int nbt,
   double ***mat, double **rhs, double **sol,
   double **resvec, double *resmax,
   double *lapms, long long int *add, long long int *mul,
   int vrblvl )
{
   double *r_d;
   const size_t szrhs = dim*sizeof(double);
   cudaMalloc((void**)&r_d,szrhs);

   double *x_d;
   const size_t szsol = dim*sizeof(double);
   cudaMalloc((void**)&x_d,szsol);

   double *A_d;
   const size_t szmat = dim*dim*sizeof(double);
   cudaMalloc((void**)&A_d,szmat);

   double *A_h = new double[dim*dim];

   *add = 0; // initialize number of additions
   *mul = 0; // initialize number of multiplications

   for(int i=0; i<degp1; i++)  // compute i-th residual vector
   {
      cudaMemcpy(r_d,rhs[i],szrhs,cudaMemcpyHostToDevice);

      for(int j=0; j<=i; j++)  // multiply mat[j] with sol[i-j]
      {
         int idx=0;
         for(int i1=0; i1<dim; i1++)
            for(int j1=0; j1<dim; j1++) A_h[idx++] = mat[j][i1][j1];
      
         cudaMemcpy(A_d,A_h,szmat,cudaMemcpyHostToDevice);
         cudaMemcpy(x_d,sol[i-j],szsol,cudaMemcpyHostToDevice);
         cudaMemcpy(A_d,A_h,szmat,cudaMemcpyHostToDevice);

         if(vrblvl > 0)
            cout << "GPU_dbl_linear_residue launches " << nbt
                 << " thread blocks in step " << i << ", " << j << endl;

         cudaEvent_t start,stop;       // to measure time spent by kernels 
         cudaEventCreate(&start);
         cudaEventCreate(&stop);
         float milliseconds;

         cudaEventRecord(start);
         dbl_bals_tail<<<nbt,szt>>>(dim,szt,A_d,x_d,r_d);
         cudaEventRecord(stop);
         cudaEventSynchronize(stop);
         cudaEventElapsedTime(&milliseconds,start,stop);
         *lapms += milliseconds;
         flopcount_dbl_bals_tail(dim,add,mul);
      }
      cudaMemcpy(resvec[i],r_d,szrhs,cudaMemcpyDeviceToHost);
   }
   if(vrblvl > 1)
   {
      for(int i=0; i<degp1; i++) 
      {
         cout << "Solution vector " << i << " :" << endl;
         for(int j=0; j<dim; j++) cout << sol[i][j] << endl;
         cout << "Residual vector " << i << " :" << endl;
         for(int j=0; j<dim; j++) cout << resvec[i][j] << endl;
      }
   }
   *resmax = 0.0;
   for(int i=0; i<degp1; i++)
   {
      double *ri = resvec[i];
      for(int j=0; j<dim; j++)
         if(abs(ri[j]) > *resmax) *resmax = abs(ri[j]);
   }
   free(A_h);

   cudaFree(r_d); cudaFree(x_d); cudaFree(A_d);
}

void GPU_cmplx_linear_residue
 ( int dim, int degp1, int szt, int nbt,
   double ***matre, double ***matim, double **rhsre, double **rhsim,
   double **solre, double **solim,
   double **resvecre, double **resvecim, double *resmax,
   double *lapms, long long int *add, long long int *mul,
   int vrblvl )
{
   double *rre_d;
   double *rim_d;
   const size_t szrhs = dim*sizeof(double);
   cudaMalloc((void**)&rre_d,szrhs);
   cudaMalloc((void**)&rim_d,szrhs);

   double *xre_d;
   double *xim_d;
   const size_t szsol = dim*sizeof(double);
   cudaMalloc((void**)&xre_d,szsol);
   cudaMalloc((void**)&xim_d,szsol);

   double *Are_d;
   double *Aim_d;
   const size_t szmat = dim*dim*sizeof(double);
   cudaMalloc((void**)&Are_d,szmat);
   cudaMalloc((void**)&Aim_d,szmat);

   double *Are_h = new double[dim*dim];
   double *Aim_h = new double[dim*dim];

   *add = 0; // initialize number of additions
   *mul = 0; // initialize number of multiplications

   for(int i=0; i<degp1; i++)  // compute i-th residual vector
   {
      cudaMemcpy(rre_d,rhsre[i],szrhs,cudaMemcpyHostToDevice);
      cudaMemcpy(rim_d,rhsim[i],szrhs,cudaMemcpyHostToDevice);

      for(int j=0; j<=i; j++)  // multiply mat[j] with sol[i-j]
      {
         int idx=0;
         for(int i1=0; i1<dim; i1++)
            for(int j1=0; j1<dim; j1++)
            {
               Are_h[idx]   = matre[j][i1][j1];
               Aim_h[idx++] = matim[j][i1][j1];
            }
      
         cudaMemcpy(Are_d,Are_h,szmat,cudaMemcpyHostToDevice);
         cudaMemcpy(Aim_d,Aim_h,szmat,cudaMemcpyHostToDevice);
         cudaMemcpy(xre_d,solre[i-j],szsol,cudaMemcpyHostToDevice);
         cudaMemcpy(xim_d,solim[i-j],szsol,cudaMemcpyHostToDevice);
         cudaMemcpy(Are_d,Are_h,szmat,cudaMemcpyHostToDevice);
         cudaMemcpy(Aim_d,Aim_h,szmat,cudaMemcpyHostToDevice);

         if(vrblvl > 0)
            cout << "GPU_cmplx_linear_residue launches " << nbt
                 << " thread blocks in step " << i << ", " << j << endl;

         cudaEvent_t start,stop;       // to measure time spent by kernels 
         cudaEventCreate(&start);
         cudaEventCreate(&stop);
         float milliseconds;

         cudaEventRecord(start);
         cmplx_bals_tail<<<nbt,szt>>>
            (dim,szt,Are_d,Aim_d,xre_d,xim_d,rre_d,rim_d);
         cudaEventRecord(stop);
         cudaEventSynchronize(stop);
         cudaEventElapsedTime(&milliseconds,start,stop);
         *lapms += milliseconds;
         flopcount_cmplx_bals_tail(dim,add,mul);
      }
      cudaMemcpy(resvecre[i],rre_d,szrhs,cudaMemcpyDeviceToHost);
      cudaMemcpy(resvecim[i],rim_d,szrhs,cudaMemcpyDeviceToHost);
   }
   if(vrblvl > 1)
   {
      for(int i=0; i<degp1; i++) 
      {
         cout << "Solution vector " << i << " :" << endl;
         for(int j=0; j<dim; j++)
            cout << solre[i][j] << "  " << solim[i][j] << endl;

         cout << "Residual vector " << i << " :" << endl;
         for(int j=0; j<dim; j++)
            cout << resvecre[i][j] << "  " << resvecim[i][j] << endl;
      }
   }
   *resmax = 0.0;
   for(int i=0; i<degp1; i++)
   {
      double *rire = resvecre[i];
      double *riim = resvecim[i];

      for(int j=0; j<dim; j++)
         if(abs(rire[j]) + abs(riim[j]) > *resmax)
            *resmax = abs(rire[j]) + abs(riim[j]);
   }
   free(Are_h); free(Aim_h);

   cudaFree(rre_d); cudaFree(rim_d);
   cudaFree(xre_d); cudaFree(xim_d);
   cudaFree(Are_d); cudaFree(Aim_d);
}
