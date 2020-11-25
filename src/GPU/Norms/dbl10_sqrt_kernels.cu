// The file dbl10_sqrt_kernels.cu contains the definitions of the
// functions with prototypes in dbl10_sqrt_kernels.h.

#include "double_double_gpufun.cu"
#include "quad_double_gpufun.cu"
#include "octo_double_gpufun.cu"
#include "deca_double_gpufun.cu"
#include "dbl10_sqrt_kernels.h"

__global__ void dbl10_sqrt
 ( double *rtb, double *rix, double *rmi, double *rrg, double *rpk,
   double *ltb, double *lix, double *lmi, double *lrg, double *lpk,
   int max_steps )
{
   double x_rtb,x_rix,x_rmi,x_rrg,x_rpk;
   double x_ltb,x_lix,x_lmi,x_lrg,x_lpk;
   double z_rtb,z_rix,z_rmi,z_rrg,z_rpk;
   double z_ltb,z_lix,z_lmi,z_lrg,z_lpk;

   int k = threadIdx.x;

   if(k == 0)
   {
      x_rtb = *rtb; x_rix = *rix; x_rmi = *rmi; x_rrg = *rrg; x_rpk = *rpk;
      x_ltb = *ltb; x_lix = *lix; x_lmi = *lmi; x_lrg = *lrg; x_lpk = *lpk;

      for(int i=0; i<max_steps; i++)
      {
         dag_sqr( x_rtb, x_rix, x_rmi, x_rrg, x_rpk,
                  x_ltb, x_lix, x_lmi, x_lrg, x_lpk,
                 &z_rtb,&z_rix,&z_rmi,&z_rrg,&z_rpk,
                 &z_ltb,&z_lix,&z_lmi,&z_lrg,&z_lpk);    // z = x*x
         dag_inc(&z_rtb,&z_rix,&z_rmi,&z_rrg,&z_rpk,
                 &z_ltb,&z_lix,&z_lmi,&z_lrg,&z_lpk,
                   *rtb,  *rix,  *rmi,  *rrg,  *rpk,
                   *ltb,  *lix,  *lmi,  *lrg,  *lpk);  // x^2 + input number
         dag_div( z_rtb, z_rix, z_rmi, z_rrg, z_rpk,
                  z_ltb, z_lix, z_lmi, z_lrg, z_lpk,
                  x_rtb, x_rix, x_rmi, x_rrg, x_rpk,
                  x_ltb, x_lix, x_lmi, x_lrg, x_lpk,
                 &z_rtb,&z_rix,&z_rmi,&z_rrg,&z_rpk,
                 &z_ltb,&z_lix,&z_lmi,&z_lrg,&z_lpk);
         dag_mul_da_d( z_rtb, z_rix, z_rmi, z_rrg, z_rpk,
                       z_ltb, z_lix, z_lmi, z_lrg, z_lpk, 0.5,
                      &z_rtb,&z_rix,&z_rmi,&z_rrg,&z_rpk,
                      &z_ltb,&z_lix,&z_lmi,&z_lrg,&z_lpk);    // z = z/(2*x)
         dag_copy( z_rtb, z_rix, z_rmi, z_rrg, z_rpk,
                   z_ltb, z_lix, z_lmi, z_lrg, z_lpk,
                  &x_rtb,&x_rix,&x_rmi,&x_rrg,&x_rpk,
                  &x_ltb,&x_lix,&x_lmi,&x_lrg,&x_lpk);
      }
      *rtb = x_rtb; *rix = x_rix; *rmi = x_rmi; *rrg = x_rrg; *rpk = x_rpk;
      *ltb = x_ltb; *lix = x_lix; *lmi = x_lmi; *lrg = x_lrg; *lpk = x_lpk;
   }
}

void GPU_dbl10_sqrt
 ( double *rtb, double *rix, double *rmi, double *rrg, double *rpk,
   double *ltb, double *lix, double *lmi, double *lrg, double *lpk,
   int maxstp )
{
   double* rtb_d;                 // highest part on device
   double* rix_d;                 // second highest part on device
   double* rmi_d;                 // third highest part on device
   double* rrg_d;                 // fourth highest part on device
   double* rpk_d;                 // fifth highest part on device
   double* ltb_d;                 // fifth lowest part on device
   double* lix_d;                 // fourth lowest part on device
   double* lmi_d;                 // third lowest part on device
   double* lrg_d;                 // second lowest part on device
   double* lpk_d;                 // lowest part on device

   size_t size = sizeof(double);
   cudaMalloc((void**)&rtb_d,size);
   cudaMalloc((void**)&rix_d,size);
   cudaMalloc((void**)&rmi_d,size);
   cudaMalloc((void**)&rrg_d,size);
   cudaMalloc((void**)&rpk_d,size);
   cudaMalloc((void**)&ltb_d,size);
   cudaMalloc((void**)&lix_d,size);
   cudaMalloc((void**)&lmi_d,size);
   cudaMalloc((void**)&lrg_d,size);
   cudaMalloc((void**)&lpk_d,size);
   cudaMemcpy(rtb_d,rtb,size,cudaMemcpyHostToDevice);
   cudaMemcpy(rix_d,rix,size,cudaMemcpyHostToDevice);
   cudaMemcpy(rmi_d,rmi,size,cudaMemcpyHostToDevice);
   cudaMemcpy(rrg_d,rrg,size,cudaMemcpyHostToDevice);
   cudaMemcpy(rpk_d,rpk,size,cudaMemcpyHostToDevice);
   cudaMemcpy(ltb_d,ltb,size,cudaMemcpyHostToDevice);
   cudaMemcpy(lix_d,lix,size,cudaMemcpyHostToDevice);
   cudaMemcpy(lmi_d,lmi,size,cudaMemcpyHostToDevice);
   cudaMemcpy(lrg_d,lrg,size,cudaMemcpyHostToDevice);
   cudaMemcpy(lpk_d,lpk,size,cudaMemcpyHostToDevice);

   dbl10_sqrt<<<1,1>>>(rtb_d,rix_d,rmi_d,rrg_d,rpk_d,
                       ltb_d,lix_d,lmi_d,lrg_d,lpk_d,maxstp);

   cudaMemcpy(rtb,rtb_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(rix,rix_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(rmi,rmi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(rrg,rrg_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(rpk,rpk_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(ltb,ltb_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(lix,lix_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(lmi,lmi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(lrg,lrg_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(lpk,lpk_d,size,cudaMemcpyDeviceToHost);
}
