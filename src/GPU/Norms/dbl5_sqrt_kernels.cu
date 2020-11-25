// The file dbl5_sqrt_kernels.cu contains the definitions of the
// functions with prototypes in dbl5_sqrt_kernels.h.

#include "double_double_gpufun.cu"
#include "penta_double_gpufun.cu"
#include "dbl5_sqrt_kernels.h"

__global__ void dbl5_sqrt
 ( double *tb, double *ix, double *mi, double *rg, double *pk,
   int max_steps )
{
   double x_tb,x_ix,x_mi,x_rg,x_pk;
   double z_tb,z_ix,z_mi,z_rg,z_pk;

   int k = threadIdx.x;

   if(k == 0)
   {
      x_tb = *tb; x_ix = *ix; x_mi = *mi; x_rg = *rg; x_pk = *pk;

      for(int i=0; i<max_steps; i++)
      {
         pdg_sqr( x_tb, x_ix, x_mi, x_rg, x_pk,
                 &z_tb,&z_ix,&z_mi,&z_rg,&z_pk);    // z = x*x
         pdg_inc(&z_tb,&z_ix,&z_mi,&z_rg,&z_pk,
                   *tb,  *ix,  *mi,  *rg,  *pk);    // x^2 + input number
         pdg_div( z_tb, z_ix, z_mi, z_rg, z_pk,
                  x_tb, x_ix, x_mi, x_rg, x_pk,
                 &z_tb,&z_ix,&z_mi,&z_rg,&z_pk);
         pdg_mul_pd_d( z_tb, z_ix, z_mi, z_rg, z_pk, 0.5,
                      &z_tb,&z_ix,&z_mi,&z_rg,&z_pk);    // z = z/(2*x)
         pdg_copy( z_tb, z_ix, z_mi, z_rg, z_pk,
                  &x_tb,&x_ix,&x_mi,&x_rg,&x_pk);
      }
      *tb = x_tb; *ix = x_ix; *mi = x_mi; *rg = x_rg; *pk = x_pk;
   }
}

void GPU_dbl5_sqrt
 ( double *tb, double *ix, double *mi, double *rg, double *pk,
   int maxstp )
{
   double* tb_d;                 // highest part on device
   double* ix_d;                 // second highest part on device
   double* mi_d;                 // middle part on device
   double* rg_d;                 // second lowest part on device
   double* pk_d;                 // lowest part on device

   size_t size = sizeof(double);
   cudaMalloc((void**)&tb_d,size);
   cudaMalloc((void**)&ix_d,size);
   cudaMalloc((void**)&mi_d,size);
   cudaMalloc((void**)&rg_d,size);
   cudaMalloc((void**)&pk_d,size);
   cudaMemcpy(tb_d,tb,size,cudaMemcpyHostToDevice);
   cudaMemcpy(ix_d,ix,size,cudaMemcpyHostToDevice);
   cudaMemcpy(mi_d,mi,size,cudaMemcpyHostToDevice);
   cudaMemcpy(rg_d,rg,size,cudaMemcpyHostToDevice);
   cudaMemcpy(pk_d,pk,size,cudaMemcpyHostToDevice);

   dbl5_sqrt<<<1,1>>>(tb_d,ix_d,mi_d,rg_d,pk_d,maxstp);

   cudaMemcpy(tb,tb_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(ix,ix_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(mi,mi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(rg,rg_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(pk,pk_d,size,cudaMemcpyDeviceToHost);
}
