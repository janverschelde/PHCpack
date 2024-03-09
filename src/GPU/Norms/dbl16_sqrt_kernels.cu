// The file dbl16_sqrt_kernels.cu contains the definitions of the
// functions with prototypes in dbl16_sqrt_kernels.h.

#include "double_double_gpufun.cu"
#include "quad_double_gpufun.cu"
#include "octo_double_gpufun.cu"
#include "hexa_double_gpufun.cu"
#include "dbl16_sqrt_kernels.h"

__global__ void dbl16_sqrt
 ( double *hihihihi, double *lohihihi, double *hilohihi, double *lolohihi,
   double *hihilohi, double *lohilohi, double *hilolohi, double *lololohi,
   double *hihihilo, double *lohihilo, double *hilohilo, double *lolohilo,
   double *hihilolo, double *lohilolo, double *hilololo, double *lolololo,
   int max_steps )
{
   double x_hihihihi,x_lohihihi,x_hilohihi,x_lolohihi;
   double x_hihilohi,x_lohilohi,x_hilolohi,x_lololohi;
   double x_hihihilo,x_lohihilo,x_hilohilo,x_lolohilo;
   double x_hihilolo,x_lohilolo,x_hilololo,x_lolololo;
   double z_hihihihi,z_lohihihi,z_hilohihi,z_lolohihi;
   double z_hihilohi,z_lohilohi,z_hilolohi,z_lololohi;
   double z_hihihilo,z_lohihilo,z_hilohilo,z_lolohilo;
   double z_hihilolo,z_lohilolo,z_hilololo,z_lolololo;

   int k = threadIdx.x;

   if(k == 0)
   {
      x_hihihihi = *hihihihi; x_lohihihi = *lohihihi;
      x_hilohihi = *hilohihi; x_lolohihi = *lolohihi;
      x_hihilohi = *hihilohi; x_lohilohi = *lohilohi;
      x_hilolohi = *hilolohi; x_lololohi = *lololohi;
      x_hihihilo = *hihihilo; x_lohihilo = *lohihilo;
      x_hilohilo = *hilohilo; x_lolohilo = *lolohilo;
      x_hihilolo = *hihilolo; x_lohilolo = *lohilolo;
      x_hilololo = *hilololo; x_lolololo = *lolololo;

      for(int i=0; i<max_steps; i++)
      {
         hdg_sqr( x_hihihihi, x_lohihihi, x_hilohihi, x_lolohihi,
                  x_hihilohi, x_lohilohi, x_hilolohi, x_lololohi,
                  x_hihihilo, x_lohihilo, x_hilohilo, x_lolohilo,
                  x_hihilolo, x_lohilolo, x_hilololo, x_lolololo,
                 &z_hihihihi,&z_lohihihi,&z_hilohihi,&z_lolohihi,
                 &z_hihilohi,&z_lohilohi,&z_hilolohi,&z_lololohi,
                 &z_hihihilo,&z_lohihilo,&z_hilohilo,&z_lolohilo,
                 &z_hihilolo,&z_lohilolo,&z_hilololo,&z_lolololo);    // z = x*x
         hdg_inc(&z_hihihihi,&z_lohihihi,&z_hilohihi,&z_lolohihi,
                 &z_hihilohi,&z_lohilohi,&z_hilolohi,&z_lololohi,
                 &z_hihihilo,&z_lohihilo,&z_hilohilo,&z_lolohilo,
                 &z_hihilolo,&z_lohilolo,&z_hilololo,&z_lolololo,
                   *hihihihi,  *lohihihi,  *hilohihi,  *lolohihi,
                   *hihilohi,  *lohilohi,  *hilolohi,  *lololohi,
                   *hihihilo,  *lohihilo,  *hilohilo,  *lolohilo,
                   *hihilolo,  *lohilolo,  *hilololo,  *lolololo);    // x^2 + input
         hdg_div( z_hihihihi, z_lohihihi, z_hilohihi, z_lolohihi,
                  z_hihilohi, z_lohilohi, z_hilolohi, z_lololohi,
                  z_hihihilo, z_lohihilo, z_hilohilo, z_lolohilo,
                  z_hihilolo, z_lohilolo, z_hilololo, z_lolololo,
                  x_hihihihi, x_lohihihi, x_hilohihi, x_lolohihi,
                  x_hihilohi, x_lohilohi, x_hilolohi, x_lololohi,
                  x_hihihilo, x_lohihilo, x_hilohilo, x_lolohilo,
                  x_hihilolo, x_lohilolo, x_hilololo, x_lolololo,
                 &z_hihihihi,&z_lohihihi,&z_hilohihi,&z_lolohihi,
                 &z_hihilohi,&z_lohilohi,&z_hilolohi,&z_lololohi,
                 &z_hihihilo,&z_lohihilo,&z_hilohilo,&z_lolohilo,
                 &z_hihilolo,&z_lohilolo,&z_hilololo,&z_lolololo);
         hdg_mul_hd_d( z_hihihihi, z_lohihihi, z_hilohihi, z_lolohihi,
                       z_hihilohi, z_lohilohi, z_hilolohi, z_lololohi,
                       z_hihihilo, z_lohihilo, z_hilohilo, z_lolohilo,
                       z_hihilolo, z_lohilolo, z_hilololo, z_lolololo, 0.5,
                      &z_hihihihi,&z_lohihihi,&z_hilohihi,&z_lolohihi,
                      &z_hihilohi,&z_lohilohi,&z_hilolohi,&z_lololohi,
                      &z_hihihilo,&z_lohihilo,&z_hilohilo,&z_lolohilo,
                      &z_hihilolo,&z_lohilolo,&z_hilololo,&z_lolololo); // z = z/(2*x)
         hdg_copy( z_hihihihi, z_lohihihi, z_hilohihi, z_lolohihi,
                   z_hihilohi, z_lohilohi, z_hilolohi, z_lololohi,
                   z_hihihilo, z_lohihilo, z_hilohilo, z_lolohilo,
                   z_hihilolo, z_lohilolo, z_hilololo, z_lolololo,
                  &x_hihihihi,&x_lohihihi,&x_hilohihi,&x_lolohihi,
                  &x_hihilohi,&x_lohilohi,&x_hilolohi,&x_lololohi,
                  &x_hihihilo,&x_lohihilo,&x_hilohilo,&x_lolohilo,
                  &x_hihilolo,&x_lohilolo,&x_hilololo,&x_lolololo);
      }
      *hihihihi = x_hihihihi; *lohihihi = x_lohihihi;
      *hilohihi = x_hilohihi; *lolohihi = x_lolohihi;
      *hihilohi = x_hihilohi; *lohilohi = x_lohilohi;
      *hilolohi = x_hilolohi; *lololohi = x_lololohi;
      *hihihilo = x_hihihilo; *lohihilo = x_lohihilo;
      *hilohilo = x_hilohilo; *lolohilo = x_lolohilo;
      *hihilolo = x_hihilolo; *lohilolo = x_lohilolo;
      *hilololo = x_hilololo; *lolololo = x_lolololo;
   }
}

void GPU_dbl16_sqrt
 ( double *hihihihi, double *lohihihi, double *hilohihi, double *lolohihi,
   double *hihilohi, double *lohilohi, double *hilolohi, double *lololohi,
   double *hihihilo, double *lohihilo, double *hilohilo, double *lolohilo,
   double *hihilolo, double *lohilolo, double *hilololo, double *lolololo,
   int maxstp )
{
   double* hihihihi_d;               // highest part on device
   double* lohihihi_d;               // second highest part on device
   double* hilohihi_d;               // third highest part on device
   double* lolohihi_d;               // fourth highest part on device
   double* hihilohi_d;               // fifth highest part on device
   double* lohilohi_d;               // sixth highst part on device
   double* hilolohi_d;               // seventh highest part on device
   double* lololohi_d;               // eigth highest part on device
   double* hihihilo_d;               // eigth lowest part on device
   double* lohihilo_d;               // seventh lowest part on device
   double* hilohilo_d;               // sixth lowst part on device
   double* lolohilo_d;               // fifth lowest part on device
   double* hihilolo_d;               // fourth lowest part on device
   double* lohilolo_d;               // third lowest part on device
   double* hilololo_d;               // second lowst part on device
   double* lolololo_d;               // lowest part on device

   size_t size = sizeof(double);
   cudaMalloc((void**)&hihihihi_d,size);
   cudaMalloc((void**)&lohihihi_d,size);
   cudaMalloc((void**)&hilohihi_d,size);
   cudaMalloc((void**)&lolohihi_d,size);
   cudaMalloc((void**)&hihilohi_d,size);
   cudaMalloc((void**)&lohilohi_d,size);
   cudaMalloc((void**)&hilolohi_d,size);
   cudaMalloc((void**)&lololohi_d,size);
   cudaMalloc((void**)&hihihilo_d,size);
   cudaMalloc((void**)&lohihilo_d,size);
   cudaMalloc((void**)&hilohilo_d,size);
   cudaMalloc((void**)&lolohilo_d,size);
   cudaMalloc((void**)&hihilolo_d,size);
   cudaMalloc((void**)&lohilolo_d,size);
   cudaMalloc((void**)&hilololo_d,size);
   cudaMalloc((void**)&lolololo_d,size);
   cudaMemcpy(hihihihi_d,hihihihi,size,cudaMemcpyHostToDevice);
   cudaMemcpy(lohihihi_d,lohihihi,size,cudaMemcpyHostToDevice);
   cudaMemcpy(hilohihi_d,hilohihi,size,cudaMemcpyHostToDevice);
   cudaMemcpy(lolohihi_d,lolohihi,size,cudaMemcpyHostToDevice);
   cudaMemcpy(hihilohi_d,hihilohi,size,cudaMemcpyHostToDevice);
   cudaMemcpy(lohilohi_d,lohilohi,size,cudaMemcpyHostToDevice);
   cudaMemcpy(hilolohi_d,hilolohi,size,cudaMemcpyHostToDevice);
   cudaMemcpy(lololohi_d,lololohi,size,cudaMemcpyHostToDevice);
   cudaMemcpy(hihihilo_d,hihihilo,size,cudaMemcpyHostToDevice);
   cudaMemcpy(lohihilo_d,lohihilo,size,cudaMemcpyHostToDevice);
   cudaMemcpy(hilohilo_d,hilohilo,size,cudaMemcpyHostToDevice);
   cudaMemcpy(lolohilo_d,lolohilo,size,cudaMemcpyHostToDevice);
   cudaMemcpy(hihilolo_d,hihilolo,size,cudaMemcpyHostToDevice);
   cudaMemcpy(lohilolo_d,lohilolo,size,cudaMemcpyHostToDevice);
   cudaMemcpy(hilololo_d,hilololo,size,cudaMemcpyHostToDevice);
   cudaMemcpy(lolololo_d,lolololo,size,cudaMemcpyHostToDevice);

   dbl16_sqrt<<<1,1>>>(hihihihi_d,lohihihi_d,hilohihi_d,lolohihi_d,
                       hihilohi_d,lohilohi_d,hilolohi_d,lololohi_d,
                       hihihilo_d,lohihilo_d,hilohilo_d,lolohilo_d,
                       hihilolo_d,lohilolo_d,hilololo_d,lolololo_d,maxstp);

   cudaMemcpy(hihihihi,hihihihi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(lohihihi,lohihihi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(hilohihi,hilohihi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(lolohihi,lolohihi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(hihilohi,hihilohi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(lohilohi,lohilohi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(hilolohi,hilolohi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(lololohi,lololohi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(hihihilo,hihihilo_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(lohihilo,lohihilo_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(hilohilo,hilohilo_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(lolohilo,lolohilo_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(hihilolo,hihilolo_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(lohilolo,lohilolo_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(hilololo,hilololo_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(lolololo,lolololo_d,size,cudaMemcpyDeviceToHost);
}
