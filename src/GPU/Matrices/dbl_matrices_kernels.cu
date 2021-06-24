// The file dbl_matrices_kernels.cu defines the kernels with prototypes
// in dbl_matrices_kernels.h.

#include "dbl_matrices_kernels.h"

__global__ void real_convolutions ( double *x, double *y, int deg1 )
{
   int j = blockIdx.x;     // convolution of j-th series in x and y
   int k = threadIdx.x;    // thread k computes k-th coefficient in product
   int offset = j*deg1+k;  // position of the k-th coefficient of the series

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
   x[offset+k] = zv[k];
}

void GPU_real_inner_product
 ( int BS, int dim, int deg, double **x, double **y, double *z )
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

   real_convolutions<<<dim,BS>>>(x_d,y_d,deg1);
}
