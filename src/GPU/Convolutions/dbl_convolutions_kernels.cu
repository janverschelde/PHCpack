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

__global__ void cmplx_convolute
 ( double *xre, double *xim, double *yre, double *yim,
   double *zre, double *zim, int dim )
{
   int k = threadIdx.x;       // thread k computes zre[k] and zim[k]

   __shared__ double xvre[d_shmemsize];
   __shared__ double xvim[d_shmemsize];
   __shared__ double yvre[d_shmemsize];
   __shared__ double yvim[d_shmemsize];
   __shared__ double zvre[d_shmemsize];
   __shared__ double zvim[d_shmemsize];

   double xr,xi,yr,yi,zr,zi;

   xvre[k] = xre[k]; xvim[k] = xim[k];
   yvre[k] = yre[k]; yvim[k] = yim[k];

   xr = xvre[0]; xi = xvim[0];    // z[k] = x[0]*y[k]
   yr = yvre[k]; yi = yvim[k];
   zr = xr*yr - xi*yi;
   zi = xr*yi + xi*yr;
   zvre[k] = zr;
   zvim[k] = zi;

   for(int i=1; i<=k; i++) // z[k] = z[k] + x[i]*y[k-i]
   {
      xr = xvre[i];   xi = xvim[i];
      yr = yvre[k-i]; yi = yvim[k-i];
      zr = xr*yr - xi*yi;
      zi = xr*yi + xi*yr;
      zvre[k] += zr;
      zvim[k] += zi;
   }
   __syncthreads();

   zre[k] = zvre[k];
   zim[k] = zvim[k];
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

void GPU_cmplx_product
 ( double *xre_h, double *xim_h, double *yre_h, double *yim_h,
   double *zre_h, double *zim_h, int deg, int freq, int BS )
{
   const int dim = deg+1;            // length of all vectors
   double* xre_d;                    // xre_d is xre_h on the device
   double* xim_d;                    // xim_d is xim_h on the device
   double* yre_d;                    // yre_d is yre_h on the device
   double* yim_d;                    // yim_d is yim_h on the device
   double* zre_d;                    // zre_d is zre_h on the device
   double* zim_d;                    // zim_d is zim_h on the device
   size_t size = dim*sizeof(double); // number of bytes for each vector

   cudaMalloc((void**)&xre_d,size);
   cudaMalloc((void**)&xim_d,size);
   cudaMalloc((void**)&yre_d,size);
   cudaMalloc((void**)&yim_d,size);
   cudaMalloc((void**)&zre_d,size);
   cudaMalloc((void**)&zim_d,size);
   cudaMemcpy(xre_d,xre_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(xim_d,xim_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yre_d,yre_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yim_d,yim_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(zre_d,zre_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(zim_d,zim_h,size,cudaMemcpyHostToDevice);

   if(dim == BS)
   {
      for(int i=0; i<freq; i++)
         cmplx_convolute<<<1,BS>>>(xre_d,xim_d,yre_d,yim_d,zre_d,zim_d,dim);
   }

   cudaMemcpy(zre_h,zre_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zim_h,zim_d,size,cudaMemcpyDeviceToHost);
}
