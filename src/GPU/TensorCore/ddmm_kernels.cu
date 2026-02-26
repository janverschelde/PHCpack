/* Definitions the kernels for a double double matrix multiplication. */

#include <iostream>
#include "double_double_gpufun.cu"

using namespace std;

__global__ void ddmm
 ( int nrows, int ncols, int dim,
   double *Ahi, double *Alo, double *Bhi, double *Blo,
   double *Chi, double *Clo )
{
   const int idx = blockIdx.x*blockDim.x + threadIdx.x;
   // element computed in C
   const int row = idx / ncols;
   const int col = idx % ncols;

   double prdhi = 0.0;
   double prdlo = 0.0;

   double a_hi,a_lo,b_hi,b_lo,acchi,acclo;

   for(int k=0; k<dim; k++)
   {
      a_hi = Ahi[row*dim+k]; a_lo = Alo[row*dim+k];
      b_hi = Bhi[col*dim+k]; b_lo = Blo[col*dim+k];

      ddg_mul(a_hi, a_lo, b_hi, b_lo, &acchi, &acclo);
      ddg_inc(&prdhi, &prdlo, acchi, acclo);
   }
   Chi[idx] = prdhi;
   Clo[idx] = prdlo;
}

void GPU_dd_matmatmul
 ( int nrows, int ncols, int dim,
   double **Ahi, double **Alo, double **BThi, double **BTlo,
   double **Chi, double **Clo, float* milliseconds )
{
   const int Adim = nrows*dim;
   const int Bdim = ncols*dim;
   const int Cdim = nrows*ncols;

   double *Ahi_h = new double[Adim];
   double *Alo_h = new double[Adim];
   double *Bhi_h = new double[Bdim];
   double *Blo_h = new double[Bdim];
   double *Chi_h = new double[Cdim];
   double *Clo_h = new double[Cdim];

   int idx = 0;
   for(int i=0; i<nrows; i++)
      for(int j=0; j<dim; j++)
      {
         Ahi_h[idx]   = Ahi[i][j]; 
         Alo_h[idx++] = Alo[i][j]; 
      }
   idx = 0;
   for(int i=0; i<ncols; i++)
      for(int j=0; j<dim; j++)
      {
         Bhi_h[idx]   = BThi[i][j]; 
         Blo_h[idx++] = BTlo[i][j]; 
      }

   double *Ahi_d;
   double *Alo_d;
   double *Bhi_d;
   double *Blo_d;
   double *Chi_d;
   double *Clo_d;

   const size_t Asize = Adim*sizeof(double);
   const size_t Bsize = Bdim*sizeof(double);
   const size_t Csize = Cdim*sizeof(double);
   cudaMalloc((void**)&Ahi_d, Asize);
   cudaMalloc((void**)&Alo_d, Asize);
   cudaMalloc((void**)&Bhi_d, Bsize);
   cudaMalloc((void**)&Blo_d, Bsize);
   cudaMalloc((void**)&Chi_d, Csize);
   cudaMalloc((void**)&Clo_d, Csize);
   cudaMemcpy(Ahi_d, Ahi_h, Asize, cudaMemcpyHostToDevice);
   cudaMemcpy(Alo_d, Alo_h, Asize, cudaMemcpyHostToDevice);
   cudaMemcpy(Bhi_d, Bhi_h, Bsize, cudaMemcpyHostToDevice);
   cudaMemcpy(Blo_d, Blo_h, Bsize, cudaMemcpyHostToDevice);

   const int sizeblock = 64;
   const int nbrblocks = (int) ceil(Cdim/((double) sizeblock));

   cout << "-> launching "
        << nbrblocks << " blocks of "
        << sizeblock << " threads ..." << endl;

   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);

   cudaEventRecord(start);
   ddmm<<<nbrblocks,sizeblock>>>
     (nrows, ncols, dim, Ahi_d, Alo_d, Bhi_d, Blo_d, Chi_d, Clo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(milliseconds, start, stop);

   cudaMemcpy(Chi_h, Chi_d, Csize, cudaMemcpyDeviceToHost);
   cudaMemcpy(Clo_h, Clo_d, Csize, cudaMemcpyDeviceToHost);

   idx = 0;
   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
      {
         Chi[i][j] = Chi_h[idx]; 
         Clo[i][j] = Clo_h[idx++]; 
      }

   free(Ahi_h); free(Alo_h);
   free(Bhi_h); free(Blo_h);
   free(Chi_h); free(Clo_h);
}
