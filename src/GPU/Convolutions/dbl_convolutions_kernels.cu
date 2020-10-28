// The file dbl_convolutions_kernels.cu defines kernels with prototypes
// in dbl_convolution_kernels.h.

#include "dbl_convolutions_kernels.h"

__global__ void dbl_convolute
 ( double *x, double *y, double *z, int dim )
{
   int k = threadIdx.x;                 // thread k computes z[k]

   __shared__ double xv[d_shmemsize];
   __shared__ double yv[d_shmemsize];
   __shared__ double zv[d_shmemsize];

   xv[k] = x[k];
   yv[k] = y[k];

   zv[k] = xv[0]*yv[k];

   for(int i=1; i<=k; i++) zv[k] = zv[k] + xv[i]*yv[k-i];

   __syncthreads();

   z[k] = zv[k];
}

void GPU_dbl_product
 ( double *x_h, double *y_h, double *z_h, int deg, int freq, int BS )
{
   const int dim = deg+1;            // length of all vectors
   double* x_d;                      // x_d is x_h on the device
   double* y_d;                      // y_d is y_h on the device
   double* z_d;                      // z_d is z_h on the device
   size_t size = dim*sizeof(double); // number of bytes for each vector

   cudaMalloc((void**)&x_d,size);
   cudaMalloc((void**)&y_d,size);
   cudaMalloc((void**)&z_d,size);
   cudaMemcpy(x_d,x_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(y_d,y_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(z_d,z_h,size,cudaMemcpyHostToDevice);

   if(dim == BS)
   {
      for(int i=0; i<freq; i++)
         dbl_convolute<<<1,BS>>>(x_d,y_d,z_d,dim);
   }

   cudaMemcpy(z_h,z_d,size,cudaMemcpyDeviceToHost);
}
