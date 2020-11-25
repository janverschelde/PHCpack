// The file dbl2_sqrt_kernels.cu contains the definitions of the
// functions with prototypes in dbl2_sqrt_kernels.h.

#include "double_double_gpufun.cu"
#include "dbl2_sqrt_kernels.h"

__global__ void dbl2_sqrt
 ( double *hi, double *lo, int max_steps )
{
   double x_hi,x_lo,z_hi,z_lo;

   int k = threadIdx.x;

   if(k == 0)
   {
      x_hi = *hi; x_lo = *lo;

      for(int i=0; i<max_steps; i++)
      {
         ddg_sqr(x_hi,x_lo,&z_hi,&z_lo);    // z = x*x
         ddg_inc(&z_hi,&z_lo,*hi,*lo);      // x^2 + input number
         ddg_div(z_hi,z_lo,x_hi,x_lo,&z_hi,&z_lo);
         ddg_mlt_d(&z_hi,&z_lo,0.5);        // z = z/(2*x)
         x_hi = z_hi; x_lo = z_lo;
      }
      *hi = x_hi; *lo = x_lo;
   }
}

void GPU_dbl2_sqrt ( double *hi, double *lo, int maxstp )
{
   double* hi_d;                 // high part on device
   double* lo_d;                 // low part on device

   size_t size = sizeof(double);
   cudaMalloc((void**)&hi_d,size);
   cudaMalloc((void**)&lo_d,size);
   cudaMemcpy(hi_d,hi,size,cudaMemcpyHostToDevice);
   cudaMemcpy(lo_d,lo,size,cudaMemcpyHostToDevice);

   dbl2_sqrt<<<1,1>>>(hi_d,lo_d,maxstp);

   cudaMemcpy(hi,hi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(lo,lo_d,size,cudaMemcpyDeviceToHost);
}
