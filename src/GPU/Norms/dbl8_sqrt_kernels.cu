// The file dbl8_sqrt_kernels.cu contains the definitions of the
// functions with prototypes in dbl8_sqrt_kernels.h.

#include "double_double_gpufun.cu"
#include "quad_double_gpufun.cu"
#include "octo_double_gpufun.cu"
#include "dbl8_sqrt_kernels.h"

__global__ void dbl8_sqrt
 ( double *hihihi, double *lohihi, double *hilohi, double *lolohi,
   double *hihilo, double *lohilo, double *hilolo, double *lololo,
   int max_steps )
{
   double x_hihihi,x_lohihi,x_hilohi,x_lolohi;
   double x_hihilo,x_lohilo,x_hilolo,x_lololo;
   double z_hihihi,z_lohihi,z_hilohi,z_lolohi;
   double z_hihilo,z_lohilo,z_hilolo,z_lololo;

   int k = threadIdx.x;

   if(k == 0)
   {
      x_hihihi = *hihihi; x_lohihi = *lohihi;
      x_hilohi = *hilohi; x_lolohi = *lolohi;
      x_hihilo = *hihilo; x_lohilo = *lohilo;
      x_hilolo = *hilolo; x_lololo = *lololo;

      for(int i=0; i<max_steps; i++)
      {
         odg_sqr( x_hihihi, x_lohihi, x_hilohi, x_lolohi,
                  x_hihilo, x_lohilo, x_hilolo, x_lololo,
                 &z_hihihi,&z_lohihi,&z_hilohi,&z_lolohi,
                 &z_hihilo,&z_lohilo,&z_hilolo,&z_lololo);    // z = x*x
         odg_inc(&z_hihihi,&z_lohihi,&z_hilohi,&z_lolohi,
                 &z_hihilo,&z_lohilo,&z_hilolo,&z_lololo,
                   *hihihi,  *lohihi,  *hilohi,  *lolohi,
                   *hihilo,  *lohilo,  *hilolo,  *lololo);    // x^2 + input
         odg_div( z_hihihi, z_lohihi, z_hilohi, z_lolohi,
                  z_hihilo, z_lohilo, z_hilolo, z_lololo,
                  x_hihihi, x_lohihi, x_hilohi, x_lolohi,
                  x_hihilo, x_lohilo, x_hilolo, x_lololo,
                 &z_hihihi,&z_lohihi,&z_hilohi,&z_lolohi,
                 &z_hihilo,&z_lohilo,&z_hilolo,&z_lololo);
         odg_mul_od_d( z_hihihi, z_lohihi, z_hilohi, z_lolohi,
                       z_hihilo, z_lohilo, z_hilolo, z_lololo, 0.5,
                      &z_hihihi,&z_lohihi,&z_hilohi,&z_lolohi,
                      &z_hihilo,&z_lohilo,&z_hilolo,&z_lololo); // z = z/(2*x)
         odg_copy( z_hihihi, z_lohihi, z_hilohi, z_lolohi,
                   z_hihilo, z_lohilo, z_hilolo, z_lololo,
                  &x_hihihi,&x_lohihi,&x_hilohi,&x_lolohi,
                  &x_hihilo,&x_lohilo,&x_hilolo,&x_lololo);
      }
      *hihihi = x_hihihi; *lohihi = x_lohihi;
      *hilohi = x_hilohi; *lolohi = x_lolohi;
      *hihilo = x_hihilo; *lohilo = x_lohilo;
      *hilolo = x_hilolo; *lololo = x_lololo;
   }
}

void GPU_dbl8_sqrt
 ( double *hihihi, double *lohihi, double *hilohi, double *lolohi,
   double *hihilo, double *lohilo, double *hilolo, double *lololo,
   int maxstp )
{
   double* hihihi_d;               // highest part on device
   double* lohihi_d;               // second highest part on device
   double* hilohi_d;               // third highest part on device
   double* lolohi_d;               // fourth highest part on device
   double* hihilo_d;               // fourth lowest part on device
   double* lohilo_d;               // third lowest part on device
   double* hilolo_d;               // second lowest part on device
   double* lololo_d;               // lowest part on device

   size_t size = sizeof(double);
   cudaMalloc((void**)&hihihi_d,size);
   cudaMalloc((void**)&lohihi_d,size);
   cudaMalloc((void**)&hilohi_d,size);
   cudaMalloc((void**)&lolohi_d,size);
   cudaMalloc((void**)&hihilo_d,size);
   cudaMalloc((void**)&lohilo_d,size);
   cudaMalloc((void**)&hilolo_d,size);
   cudaMalloc((void**)&lololo_d,size);
   cudaMemcpy(hihihi_d,hihihi,size,cudaMemcpyHostToDevice);
   cudaMemcpy(lohihi_d,lohihi,size,cudaMemcpyHostToDevice);
   cudaMemcpy(hilohi_d,hilohi,size,cudaMemcpyHostToDevice);
   cudaMemcpy(lolohi_d,lolohi,size,cudaMemcpyHostToDevice);
   cudaMemcpy(hihilo_d,hihilo,size,cudaMemcpyHostToDevice);
   cudaMemcpy(lohilo_d,lohilo,size,cudaMemcpyHostToDevice);
   cudaMemcpy(hilolo_d,hilolo,size,cudaMemcpyHostToDevice);
   cudaMemcpy(lololo_d,lololo,size,cudaMemcpyHostToDevice);

   dbl8_sqrt<<<1,1>>>(hihihi_d,lohihi_d,hilohi_d,lolohi_d,
                      hihilo_d,lohilo_d,hilolo_d,lololo_d,maxstp);

   cudaMemcpy(hihihi,hihihi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(lohihi,lohihi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(hilohi,hilohi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(lolohi,lolohi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(hihilo,hihilo_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(lohilo,lohilo_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(hilolo,hilolo_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(lololo,lololo_d,size,cudaMemcpyDeviceToHost);
}
