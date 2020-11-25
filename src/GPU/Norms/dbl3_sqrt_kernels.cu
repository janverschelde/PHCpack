// The file dbl3_sqrt_kernels.cu contains the definitions of the
// functions with prototypes in dbl3_sqrt_kernels.h.

#include "double_double_gpufun.cu"
#include "triple_double_gpufun.cu"
#include "dbl3_sqrt_kernels.h"

__global__ void dbl3_sqrt
 ( double *hi, double *mi, double *lo, int max_steps )
{
   double x_hi,x_mi,x_lo,z_hi,z_mi,z_lo;

   int k = threadIdx.x;

   if(k == 0)
   {
      x_hi = *hi; x_mi = *mi; x_lo = *lo;

      for(int i=0; i<max_steps; i++)
      {
         tdg_sqr(x_hi,x_mi,x_lo,&z_hi,&z_mi,&z_lo); // z = x*x
         tdg_inc(&z_hi,&z_mi,&z_lo,*hi,*mi,*lo);    // x^2 + input number
         tdg_div(z_hi,z_mi,z_lo,x_hi,x_mi,x_lo,&z_hi,&z_mi,&z_lo);
         tdg_mul_td_d(z_hi,z_mi,z_lo,0.5,&z_hi,&z_mi,&z_lo); // z = z/(2*x)
         tdg_copy(z_hi,z_mi,z_lo,&x_hi,&x_mi,&x_lo);
      }
      *hi = x_hi; *mi = x_mi; *lo = x_lo;
   }
}

void GPU_dbl3_sqrt
 ( double *hi, double *mi, double *lo, int maxstp )
{
   double* hi_d;                 // high part on device
   double* mi_d;                 // middle part on device
   double* lo_d;                 // low part on device

   size_t size = sizeof(double);
   cudaMalloc((void**)&hi_d,size);
   cudaMalloc((void**)&mi_d,size);
   cudaMalloc((void**)&lo_d,size);
   cudaMemcpy(hi_d,hi,size,cudaMemcpyHostToDevice);
   cudaMemcpy(mi_d,mi,size,cudaMemcpyHostToDevice);
   cudaMemcpy(lo_d,lo,size,cudaMemcpyHostToDevice);

   dbl3_sqrt<<<1,1>>>(hi_d,mi_d,lo_d,maxstp);

   cudaMemcpy(hi,hi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(mi,mi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(lo,lo_d,size,cudaMemcpyDeviceToHost);
}
