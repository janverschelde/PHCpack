// The file dbl4_sqrt_kernels.cu contains the definitions of the
// functions with prototypes in dbl4_sqrt_kernels.h.

#include "double_double_gpufun.cu"
#include "quad_double_gpufun.cu"
#include "dbl4_sqrt_kernels.h"

__global__ void dbl4_sqrt
 ( double *hihi, double *lohi, double *hilo, double *lolo, int max_steps )
{
   double x_hihi,x_lohi,x_hilo,x_lolo;
   double z_hihi,z_lohi,z_hilo,z_lolo;

   int k = threadIdx.x;

   if(k == 0)
   {
      x_hihi = *hihi; x_lohi = *lohi; x_hilo = *hilo; x_lolo = *lolo;

      for(int i=0; i<max_steps; i++)
      {
         qdg_sqr( x_hihi, x_lohi, x_hilo, x_lolo,
                 &z_hihi,&z_lohi,&z_hilo,&z_lolo);    // z = x*x
         qdg_inc(&z_hihi,&z_lohi,&z_hilo,&z_lolo,
                   *hihi,  *lohi,  *hilo,  *lolo);    // x^2 + input number
         qdg_div( z_hihi, z_lohi, z_hilo, z_lolo,
                  x_hihi, x_lohi, x_hilo, x_lolo,
                 &z_hihi,&z_lohi,&z_hilo,&z_lolo);
         qdg_mul_qd_d( z_hihi, z_lohi, z_hilo, z_lolo, 0.5,
                      &z_hihi,&z_lohi,&z_hilo,&z_lolo);    // z = z/(2*x)
         qdg_copy( z_hihi, z_lohi, z_hilo, z_lolo,
                  &x_hihi,&x_lohi,&x_hilo,&x_lolo);
      }
      *hihi = x_hihi; *lohi = x_lohi; *hilo = x_hilo; *lolo = x_lolo;
   }
}

void GPU_dbl4_sqrt
 ( double *hihi, double *lohi, double *hilo, double *lolo, int maxstp )
{
   double* hihi_d;                 // highest part on device
   double* lohi_d;                 // second highest part on device
   double* hilo_d;                 // second lowest part on device
   double* lolo_d;                 // lowest part on device

   size_t size = sizeof(double);
   cudaMalloc((void**)&hihi_d,size);
   cudaMalloc((void**)&lohi_d,size);
   cudaMalloc((void**)&hilo_d,size);
   cudaMalloc((void**)&lolo_d,size);
   cudaMemcpy(hihi_d,hihi,size,cudaMemcpyHostToDevice);
   cudaMemcpy(lohi_d,lohi,size,cudaMemcpyHostToDevice);
   cudaMemcpy(hilo_d,hilo,size,cudaMemcpyHostToDevice);
   cudaMemcpy(lolo_d,lolo,size,cudaMemcpyHostToDevice);

   dbl4_sqrt<<<1,1>>>(hihi_d,lohi_d,hilo_d,lolo_d,maxstp);

   cudaMemcpy(hihi,hihi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(lohi,lohi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(hilo,hilo_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(lolo,lolo_d,size,cudaMemcpyDeviceToHost);
}
