// The file dbl_matrices_kernels.cu defines the kernels with prototypes
// in dbl_matrices_kernels.h.

#include <iostream>
#include <cmath>
#include "dbl_matrices_kernels.h"

__global__ void dbl_convolutions ( double *x, double *y, int deg1 )
{
   const int j = blockIdx.x;     // convolution of j-th series in x and y
   const int k = threadIdx.x;    // thread k computes k-th coefficient 
   const int offset = j*deg1+k;  // position of the k-th coefficient

   __shared__ double xv[d_shmemsize];
   __shared__ double yv[d_shmemsize];
   __shared__ double zv[d_shmemsize];

   int idx = deg1+k;

   xv[k] = x[offset];
   yv[k] = 0.0;         // padded with zeros
   yv[idx] = y[offset];

   zv[k] = xv[0]*yv[idx];

   for(int i=1; i<deg1; i++)
   {
      idx = deg1 + k - i;
      zv[k] = zv[k] + xv[i]*yv[idx];
   }
   x[offset] = zv[k];
}

__global__ void dbl_additions ( double *x, int lag, int deg1 )
{
   const int j = blockIdx.x; 
   const int k = threadIdx.x;
   const int x_offset = j*deg1+k;
   int y_offset;

   __shared__ double xv[d_shmemsize];
   __shared__ double yv[d_shmemsize];

   xv[k] = x[x_offset];

   if(j < lag)
   {
      y_offset = x_offset + lag*deg1;
      yv[k] = x[y_offset];
      xv[k] += yv[k];
      x[x_offset] = xv[k]; // store for next block in next round
   }
}

void GPU_dbl_inner_product
 ( int BS, int dim, int deg, double **x, double **y, double *z,
   int mode )
{
   const int deg1 = deg+1;         // coefficient series length

   double* x_d;                    // x_d is x_h on the device
   double* y_d;                    // y_d is y_h on the device
   double* z_d;                    // z_d is z_h on the device

   size_t szdeg = deg1*sizeof(double);
   size_t szdim = dim*szdeg;

   double* x_h = new double[dim*deg1];
   double* y_h = new double[dim*deg1];
   int ix = 0;
   for(int i=0; i<dim; i++)
      for(int j=0; j<deg1; j++)
      {
         x_h[ix] = x[i][j]; y_h[ix++] = y[i][j];
      }

   cudaMalloc((void**)&x_d,szdim);
   cudaMalloc((void**)&y_d,szdim);
   cudaMalloc((void**)&z_d,szdeg);
   cudaMemcpy(x_d,x_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(y_d,y_h,szdim,cudaMemcpyHostToDevice);

   if(BS == deg1) dbl_convolutions<<<dim,BS>>>(x_d,y_d,deg1);

   if(mode==1) // do all additions on the host
   {
      cudaMemcpy(x_h,x_d,szdim,cudaMemcpyDeviceToHost);

      for(int i=0; i<deg1; i++) z[i] = 0.0;

      int ix=0;
      for(int j=0; j<dim; j++)
         for(int i=0; i<deg1; i++) z[i] = z[i] + x_h[ix++];
   }
   else
   {
      double logdim = log2((double)dim);
      double ceil_logdim = ceil(logdim);
      int ceil2log = int(ceil_logdim);
      std::cout << "log(" << dim << ") : " << logdim << std::endl;
      std::cout << "ceil(log(" << dim << ")) : " << ceil_logdim << std::endl;
      std::cout << "ceil2log : " << ceil2log << std::endl;

      if(BS == deg1)
      {
          int lag = dim/2;
          for(int L=0; L<ceil2log; L++)
          {
              dbl_additions<<<lag,BS>>>(x_d,lag,deg1);
              lag = lag/2;
          }
      }
      cudaMemcpy(z,x_d,szdeg,cudaMemcpyDeviceToHost);
   }
   free(x_h); free(y_h);
}
