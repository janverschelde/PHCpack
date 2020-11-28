// The file dbl10_convolutions_kernels.cu defines kernels with prototypes
// in dbl10_convolutions_kernels.h.

#ifdef gpufun
#include "double_double_gpufun.cu"
#include "quad_double_gpufun.cu"
#include "octo_double_gpufun.cu"
#include "deca_double_gpufun.cu"
#endif
#include "dbl10_convolutions_kernels.h"

__global__ void dbl10_increment
 ( double *xrtb, double *xrix, double *xrmi, double *xrrg, double *xrpk,
   double *xltb, double *xlix, double *xlmi, double *xlrg, double *xlpk,
   double *yrtb, double *yrix, double *yrmi, double *yrrg, double *yrpk,
   double *yltb, double *ylix, double *ylmi, double *ylrg, double *ylpk,
   double *zrtb, double *zrix, double *zrmi, double *zrrg, double *zrpk,
   double *zltb, double *zlix, double *zlmi, double *zlrg, double *zlpk,
   int dim )
{
   int k = threadIdx.x;                 // thread k computes z[k]

   __shared__ double xvrtb[da_shmemsize];
   __shared__ double xvrix[da_shmemsize];
   __shared__ double xvrmi[da_shmemsize];
   __shared__ double xvrrg[da_shmemsize];
   __shared__ double xvrpk[da_shmemsize];
   __shared__ double xvltb[da_shmemsize];
   __shared__ double xvlix[da_shmemsize];
   __shared__ double xvlmi[da_shmemsize];
   __shared__ double xvlrg[da_shmemsize];
   __shared__ double xvlpk[da_shmemsize];
   __shared__ double yvrtb[da_shmemsize];
   __shared__ double yvrix[da_shmemsize];
   __shared__ double yvrmi[da_shmemsize];
   __shared__ double yvrrg[da_shmemsize];
   __shared__ double yvrpk[da_shmemsize];
   __shared__ double yvltb[da_shmemsize];
   __shared__ double yvlix[da_shmemsize];
   __shared__ double yvlmi[da_shmemsize];
   __shared__ double yvlrg[da_shmemsize];
   __shared__ double yvlpk[da_shmemsize];
   __shared__ double zvrtb[da_shmemsize];
   __shared__ double zvrix[da_shmemsize];
   __shared__ double zvrmi[da_shmemsize];
   __shared__ double zvrrg[da_shmemsize];
   __shared__ double zvrpk[da_shmemsize];
   __shared__ double zvltb[da_shmemsize];
   __shared__ double zvlix[da_shmemsize];
   __shared__ double zvlmi[da_shmemsize];
   __shared__ double zvlrg[da_shmemsize];
   __shared__ double zvlpk[da_shmemsize];

   xvrtb[k] = xrtb[k]; xvrix[k] = xrix[k]; xvrmi[k] = xrmi[k];
   xvrrg[k] = xrrg[k]; xvrpk[k] = xrpk[k];
   xvltb[k] = xltb[k]; xvlix[k] = xlix[k]; xvlmi[k] = xlmi[k];
   xvlrg[k] = xlrg[k]; xvlpk[k] = xlpk[k];
   yvrtb[k] = yrtb[k]; yvrix[k] = yrix[k]; yvrmi[k] = yrmi[k];
   yvrrg[k] = yrrg[k]; yvrpk[k] = yrpk[k];
   yvltb[k] = yltb[k]; yvlix[k] = ylix[k]; yvlmi[k] = ylmi[k];
   yvlrg[k] = ylrg[k]; yvlpk[k] = ylpk[k];

   __syncthreads();

   dag_add(xvrtb[k],xvrix[k],xvrmi[k],xvrrg[k],xvrpk[k],
           xvltb[k],xvlix[k],xvlmi[k],xvlrg[k],xvlpk[k],
           yvrtb[k],yvrix[k],yvrmi[k],yvrrg[k],yvrpk[k],
           yvltb[k],yvlix[k],yvlmi[k],yvlrg[k],yvlpk[k],
           &zvrtb[k],&zvrix[k],&zvrmi[k],&zvrrg[k],&zvrpk[k],
           &zvltb[k],&zvlix[k],&zvlmi[k],&zvlrg[k],&zvlpk[k]);

   __syncthreads();

   zrtb[k] = zvrtb[k]; zrix[k] = zvrix[k]; zrmi[k] = zvrmi[k];
   zrrg[k] = zvrrg[k]; zrpk[k] = zvrpk[k];
   zltb[k] = zvltb[k]; zlix[k] = zvlix[k]; zlmi[k] = zvlmi[k];
   zlrg[k] = zvlrg[k]; zlpk[k] = zvlpk[k];
}

__global__ void dbl10_decrement
 ( double *xrtb, double *xrix, double *xrmi, double *xrrg, double *xrpk,
   double *xltb, double *xlix, double *xlmi, double *xlrg, double *xlpk,
   double *yrtb, double *yrix, double *yrmi, double *yrrg, double *yrpk,
   double *yltb, double *ylix, double *ylmi, double *ylrg, double *ylpk,
   double *zrtb, double *zrix, double *zrmi, double *zrrg, double *zrpk,
   double *zltb, double *zlix, double *zlmi, double *zlrg, double *zlpk,
   int dim )
{
   int k = threadIdx.x;                 // thread k computes z[k]

   __shared__ double xvrtb[da_shmemsize];
   __shared__ double xvrix[da_shmemsize];
   __shared__ double xvrmi[da_shmemsize];
   __shared__ double xvrrg[da_shmemsize];
   __shared__ double xvrpk[da_shmemsize];
   __shared__ double xvltb[da_shmemsize];
   __shared__ double xvlix[da_shmemsize];
   __shared__ double xvlmi[da_shmemsize];
   __shared__ double xvlrg[da_shmemsize];
   __shared__ double xvlpk[da_shmemsize];
   __shared__ double yvrtb[da_shmemsize];
   __shared__ double yvrix[da_shmemsize];
   __shared__ double yvrmi[da_shmemsize];
   __shared__ double yvrrg[da_shmemsize];
   __shared__ double yvrpk[da_shmemsize];
   __shared__ double yvltb[da_shmemsize];
   __shared__ double yvlix[da_shmemsize];
   __shared__ double yvlmi[da_shmemsize];
   __shared__ double yvlrg[da_shmemsize];
   __shared__ double yvlpk[da_shmemsize];
   __shared__ double zvrtb[da_shmemsize];
   __shared__ double zvrix[da_shmemsize];
   __shared__ double zvrmi[da_shmemsize];
   __shared__ double zvrrg[da_shmemsize];
   __shared__ double zvrpk[da_shmemsize];
   __shared__ double zvltb[da_shmemsize];
   __shared__ double zvlix[da_shmemsize];
   __shared__ double zvlmi[da_shmemsize];
   __shared__ double zvlrg[da_shmemsize];
   __shared__ double zvlpk[da_shmemsize];

   xvrtb[k] = xrtb[k]; xvrix[k] = xrix[k]; xvrmi[k] = xrmi[k];
   xvrrg[k] = xrrg[k]; xvrpk[k] = xrpk[k];
   xvltb[k] = xltb[k]; xvlix[k] = xlix[k]; xvlmi[k] = xlmi[k];
   xvlrg[k] = xlrg[k]; xvlpk[k] = xlpk[k];
   yvrtb[k] = yrtb[k]; yvrix[k] = yrix[k]; yvrmi[k] = yrmi[k];
   yvrrg[k] = yrrg[k]; yvrpk[k] = yrpk[k];
   yvltb[k] = yltb[k]; yvlix[k] = ylix[k]; yvlmi[k] = ylmi[k];
   yvlrg[k] = ylrg[k]; yvlpk[k] = ylpk[k];

   __syncthreads();

   dag_sub(xvrtb[k],xvrix[k],xvrmi[k],xvrrg[k],xvrpk[k],
           xvltb[k],xvlix[k],xvlmi[k],xvlrg[k],xvlpk[k],
           yvrtb[k],yvrix[k],yvrmi[k],yvrrg[k],yvrpk[k],
           yvltb[k],yvlix[k],yvlmi[k],yvlrg[k],yvlpk[k],
           &zvrtb[k],&zvrix[k],&zvrmi[k],&zvrrg[k],&zvrpk[k],
           &zvltb[k],&zvlix[k],&zvlmi[k],&zvlrg[k],&zvlpk[k]);

   __syncthreads();

   zrtb[k] = zvrtb[k]; zrix[k] = zvrix[k]; zrmi[k] = zvrmi[k];
   zrrg[k] = zvrrg[k]; zrpk[k] = zvrpk[k];
   zltb[k] = zvltb[k]; zlix[k] = zvlix[k]; zlmi[k] = zvlmi[k];
   zlrg[k] = zvlrg[k]; zlpk[k] = zvlpk[k];
}

__global__ void dbl10_convolute
 ( double *xrtb, double *xrix, double *xrmi, double *xrrg, double *xrpk,
   double *xltb, double *xlix, double *xlmi, double *xlrg, double *xlpk,
   double *yrtb, double *yrix, double *yrmi, double *yrrg, double *yrpk,
   double *yltb, double *ylix, double *ylmi, double *ylrg, double *ylpk,
   double *zrtb, double *zrix, double *zrmi, double *zrrg, double *zrpk,
   double *zltb, double *zlix, double *zlmi, double *zlrg, double *zlpk,
   int dim )
{
   int k = threadIdx.x;                 // thread k computes z[k]

   __shared__ double xvrtb[da_shmemsize];
   __shared__ double xvrix[da_shmemsize];
   __shared__ double xvrmi[da_shmemsize];
   __shared__ double xvrrg[da_shmemsize];
   __shared__ double xvrpk[da_shmemsize];
   __shared__ double xvltb[da_shmemsize];
   __shared__ double xvlix[da_shmemsize];
   __shared__ double xvlmi[da_shmemsize];
   __shared__ double xvlrg[da_shmemsize];
   __shared__ double xvlpk[da_shmemsize];
   __shared__ double yvrtb[da_shmemsize];
   __shared__ double yvrix[da_shmemsize];
   __shared__ double yvrmi[da_shmemsize];
   __shared__ double yvrrg[da_shmemsize];
   __shared__ double yvrpk[da_shmemsize];
   __shared__ double yvltb[da_shmemsize];
   __shared__ double yvlix[da_shmemsize];
   __shared__ double yvlmi[da_shmemsize];
   __shared__ double yvlrg[da_shmemsize];
   __shared__ double yvlpk[da_shmemsize];
   __shared__ double zvrtb[da_shmemsize];
   __shared__ double zvrix[da_shmemsize];
   __shared__ double zvrmi[da_shmemsize];
   __shared__ double zvrrg[da_shmemsize];
   __shared__ double zvrpk[da_shmemsize];
   __shared__ double zvltb[da_shmemsize];
   __shared__ double zvlix[da_shmemsize];
   __shared__ double zvlmi[da_shmemsize];
   __shared__ double zvlrg[da_shmemsize];
   __shared__ double zvlpk[da_shmemsize];

   double prdrtb,prdrix,prdrmi,prdrrg,prdrpk;
   double prdltb,prdlix,prdlmi,prdlrg,prdlpk;

   xvrtb[k] = xrtb[k]; xvrix[k] = xrix[k]; xvrmi[k] = xrmi[k];
   xvrrg[k] = xrrg[k]; xvrpk[k] = xrpk[k];
   xvltb[k] = xltb[k]; xvlix[k] = xlix[k]; xvlmi[k] = xlmi[k];
   xvlrg[k] = xlrg[k]; xvlpk[k] = xlpk[k];
   yvrtb[k] = yrtb[k]; yvrix[k] = yrix[k]; yvrmi[k] = yrmi[k];
   yvrrg[k] = yrrg[k]; yvrpk[k] = yrpk[k];
   yvltb[k] = yltb[k]; yvlix[k] = ylix[k]; yvlmi[k] = ylmi[k];
   yvlrg[k] = ylrg[k]; yvlpk[k] = ylpk[k];

   __syncthreads();

   // zv[k] = xv[0]*yv[k];
   dag_mul(xvrtb[0],xvrix[0],xvrmi[0],xvrrg[0],xvrpk[0],
           xvltb[0],xvlix[0],xvlmi[0],xvlrg[0],xvlpk[0],
           yvrtb[k],yvrix[k],yvrmi[k],yvrrg[k],yvrpk[k],
           yvltb[k],yvlix[k],yvlmi[k],yvlrg[k],yvlpk[k],
           &zvrtb[k],&zvrix[k],&zvrmi[k],&zvrrg[k],&zvrpk[k],
           &zvltb[k],&zvlix[k],&zvlmi[k],&zvlrg[k],&zvlpk[k]);

   for(int i=1; i<=k; i++) // zv[k] = zv[k] + xv[i]*yv[k-i];
   {
      dag_mul(xvrtb[i],xvrix[i],xvrmi[i],xvrrg[i],xvrpk[i],
              xvltb[i],xvlix[i],xvlmi[i],xvlrg[i],xvlpk[i],
              yvrtb[k-i],yvrix[k-i],yvrmi[k-i],yvrrg[k-i],yvrpk[k-i],
              yvltb[k-i],yvlix[k-i],yvlmi[k-i],yvlrg[k-i],yvlpk[k-i],
              &prdrtb,&prdrix,&prdrmi,&prdrrg,&prdrpk,
              &prdltb,&prdlix,&prdlmi,&prdlrg,&prdlpk);
      dag_inc(&zvrtb[k],&zvrix[k],&zvrmi[k],&zvrrg[k],&zvrpk[k],
              &zvltb[k],&zvlix[k],&zvlmi[k],&zvlrg[k],&zvlpk[k],
              prdrtb,prdrix,prdrmi,prdrrg,prdrpk,
              prdltb,prdlix,prdlmi,prdlrg,prdlpk);
   }

   __syncthreads();

   zrtb[k] = zvrtb[k]; zrix[k] = zvrix[k]; zrmi[k] = zvrmi[k];
   zrrg[k] = zvrrg[k]; zrpk[k] = zvrpk[k];
   zltb[k] = zvltb[k]; zlix[k] = zvlix[k]; zlmi[k] = zvlmi[k];
   zlrg[k] = zvlrg[k]; zlpk[k] = zvlpk[k];
}

__global__ void dbl10_padded_convolute
 ( double *xrtb, double *xrix, double *xrmi, double *xrrg, double *xrpk,
   double *xltb, double *xlix, double *xlmi, double *xlrg, double *xlpk,
   double *yrtb, double *yrix, double *yrmi, double *yrrg, double *yrpk,
   double *yltb, double *ylix, double *ylmi, double *ylrg, double *ylpk,
   double *zrtb, double *zrix, double *zrmi, double *zrrg, double *zrpk,
   double *zltb, double *zlix, double *zlmi, double *zlrg, double *zlpk,
   int dim )
{
   int k = threadIdx.x;                 // thread k computes z[k]

   __shared__ double xvrtb[da_shmemsize];
   __shared__ double xvrix[da_shmemsize];
   __shared__ double xvrmi[da_shmemsize];
   __shared__ double xvrrg[da_shmemsize];
   __shared__ double xvrpk[da_shmemsize];
   __shared__ double xvltb[da_shmemsize];
   __shared__ double xvlix[da_shmemsize];
   __shared__ double xvlmi[da_shmemsize];
   __shared__ double xvlrg[da_shmemsize];
   __shared__ double xvlpk[da_shmemsize];
   __shared__ double yvrtb[da_shmemsize];
   __shared__ double yvrix[da_shmemsize];
   __shared__ double yvrmi[da_shmemsize];
   __shared__ double yvrrg[da_shmemsize];
   __shared__ double yvrpk[da_shmemsize];
   __shared__ double yvltb[da_shmemsize];
   __shared__ double yvlix[da_shmemsize];
   __shared__ double yvlmi[da_shmemsize];
   __shared__ double yvlrg[da_shmemsize];
   __shared__ double yvlpk[da_shmemsize];
   __shared__ double zvrtb[da_shmemsize];
   __shared__ double zvrix[da_shmemsize];
   __shared__ double zvrmi[da_shmemsize];
   __shared__ double zvrrg[da_shmemsize];
   __shared__ double zvrpk[da_shmemsize];
   __shared__ double zvltb[da_shmemsize];
   __shared__ double zvlix[da_shmemsize];
   __shared__ double zvlmi[da_shmemsize];
   __shared__ double zvlrg[da_shmemsize];
   __shared__ double zvlpk[da_shmemsize];

   double prdrtb,prdrix,prdrmi,prdrrg,prdrpk;
   double prdltb,prdlix,prdlmi,prdlrg,prdlpk;
   int idx = dim+k;

   xvrtb[k] = xrtb[k]; xvrix[k] = xrix[k]; xvrmi[k] = xrmi[k];
   xvrrg[k] = xrrg[k]; xvrpk[k] = xrpk[k];
   xvltb[k] = xltb[k]; xvlix[k] = xlix[k]; xvlmi[k] = xlmi[k];
   xvlrg[k] = xlrg[k]; xvlpk[k] = xlpk[k];
   yvrtb[k] = 0.0; yvrix[k] = 0.0; yvrmi[k] = 0.0;
   yvrrg[k] = 0.0; yvrpk[k] = 0.0;
   yvltb[k] = 0.0; yvlix[k] = 0.0; yvlmi[k] = 0.0;
   yvlrg[k] = 0.0; yvlpk[k] = 0.0;
   yvrtb[idx] = yrtb[k]; yvrix[idx] = yrix[k]; yvrmi[idx] = yrmi[k];
   yvrrg[idx] = yrrg[k]; yvrpk[idx] = yrpk[k];
   yvltb[idx] = yltb[k]; yvlix[idx] = ylix[k]; yvlmi[idx] = ylmi[k];
   yvlrg[idx] = ylrg[k]; yvlpk[idx] = ylpk[k];

   __syncthreads();

   // zv[k] = xv[0]*yv[k];
   dag_mul(xvrtb[0],xvrix[0],xvrmi[0],xvrrg[0],xvrpk[0],
           xvltb[0],xvlix[0],xvlmi[0],xvlrg[0],xvlpk[0],
           yvrtb[idx],yvrix[idx],yvrmi[idx],yvrrg[idx],yvrpk[idx],
           yvltb[idx],yvlix[idx],yvlmi[idx],yvlrg[idx],yvlpk[idx],
           &zvrtb[k],&zvrix[k],&zvrmi[k],&zvrrg[k],&zvrpk[k],
           &zvltb[k],&zvlix[k],&zvlmi[k],&zvlrg[k],&zvlpk[k]);
   __syncthreads();

   for(int i=1; i<dim; i++) // zv[k] = zv[k] + xv[i]*yv[k-i];
   {
      int idx = dim + k - i;
      dag_mul(xvrtb[i],xvrix[i],xvrmi[i],xvrrg[i],xvrpk[i],
              xvltb[i],xvlix[i],xvlmi[i],xvlrg[i],xvlpk[i],
              yvrtb[idx],yvrix[idx],yvrmi[idx],yvrrg[idx],yvrpk[idx],
              yvltb[idx],yvlix[idx],yvlmi[idx],yvlrg[idx],yvlpk[idx],
              &prdrtb,&prdrix,&prdrmi,&prdrrg,&prdrpk,
              &prdltb,&prdlix,&prdlmi,&prdlrg,&prdlpk);
      __syncthreads();
      dag_inc(&zvrtb[k],&zvrix[k],&zvrmi[k],&zvrrg[k],&zvrpk[k],
              &zvltb[k],&zvlix[k],&zvlmi[k],&zvlrg[k],&zvlpk[k],
              prdrtb,prdrix,prdrmi,prdrrg,prdrpk,
              prdltb,prdlix,prdlmi,prdlrg,prdlpk);
      __syncthreads();
   }
   zrtb[k] = zvrtb[k]; zrix[k] = zvrix[k]; zrmi[k] = zvrmi[k];
   zrrg[k] = zvrrg[k]; zrpk[k] = zvrpk[k];
   zltb[k] = zvltb[k]; zlix[k] = zvlix[k]; zlmi[k] = zvlmi[k];
   zlrg[k] = zvlrg[k]; zlpk[k] = zvlpk[k];
}

__global__ void cmplx10_convolute
 ( double *xrertb, double *xrerix, double *xrermi, double *xrerrg,
   double *xrerpk, double *xreltb, double *xrelix, double *xrelmi,
   double *xrelrg, double *xrelpk,
   double *ximrtb, double *ximrix, double *ximrmi, double *ximrrg,
   double *ximrpk, double *ximltb, double *ximlix, double *ximlmi,
   double *ximlrg, double *ximlpk,
   double *yrertb, double *yrerix, double *yrermi, double *yrerrg,
   double *yrerpk, double *yreltb, double *yrelix, double *yrelmi,
   double *yrelrg, double *yrelpk,
   double *yimrtb, double *yimrix, double *yimrmi, double *yimrrg,
   double *yimrpk, double *yimltb, double *yimlix, double *yimlmi,
   double *yimlrg, double *yimlpk,
   double *zrertb, double *zrerix, double *zrermi, double *zrerrg,
   double *zrerpk, double *zreltb, double *zrelix, double *zrelmi,
   double *zrelrg, double *zrelpk,
   double *zimrtb, double *zimrix, double *zimrmi, double *zimrrg,
   double *zimrpk, double *zimltb, double *zimlix, double *zimlmi,
   double *zimlrg, double *zimlpk, int dim )
{
   int k = threadIdx.x;       // thread k computes zre[k] and zim[k]

   __shared__ double xvrertb[da_shmemsize];
   __shared__ double xvrerix[da_shmemsize];
   __shared__ double xvrermi[da_shmemsize];
   __shared__ double xvrerrg[da_shmemsize];
   __shared__ double xvrerpk[da_shmemsize];
   __shared__ double xvreltb[da_shmemsize];
   __shared__ double xvrelix[da_shmemsize];
   __shared__ double xvrelmi[da_shmemsize];
   __shared__ double xvrelrg[da_shmemsize];
   __shared__ double xvrelpk[da_shmemsize];
   __shared__ double xvimrtb[da_shmemsize];
   __shared__ double xvimrix[da_shmemsize];
   __shared__ double xvimrmi[da_shmemsize];
   __shared__ double xvimrrg[da_shmemsize];
   __shared__ double xvimrpk[da_shmemsize];
   __shared__ double xvimltb[da_shmemsize];
   __shared__ double xvimlix[da_shmemsize];
   __shared__ double xvimlmi[da_shmemsize];
   __shared__ double xvimlrg[da_shmemsize];
   __shared__ double xvimlpk[da_shmemsize];
   __shared__ double yvrertb[da_shmemsize];
   __shared__ double yvrerix[da_shmemsize];
   __shared__ double yvrermi[da_shmemsize];
   __shared__ double yvrerrg[da_shmemsize];
   __shared__ double yvrerpk[da_shmemsize];
   __shared__ double yvreltb[da_shmemsize];
   __shared__ double yvrelix[da_shmemsize];
   __shared__ double yvrelmi[da_shmemsize];
   __shared__ double yvrelrg[da_shmemsize];
   __shared__ double yvrelpk[da_shmemsize];
   __shared__ double yvimrtb[da_shmemsize];
   __shared__ double yvimrix[da_shmemsize];
   __shared__ double yvimrmi[da_shmemsize];
   __shared__ double yvimrrg[da_shmemsize];
   __shared__ double yvimrpk[da_shmemsize];
   __shared__ double yvimltb[da_shmemsize];
   __shared__ double yvimlix[da_shmemsize];
   __shared__ double yvimlmi[da_shmemsize];
   __shared__ double yvimlrg[da_shmemsize];
   __shared__ double yvimlpk[da_shmemsize];
   __shared__ double zvrertb[da_shmemsize];
   __shared__ double zvrerix[da_shmemsize];
   __shared__ double zvrermi[da_shmemsize];
   __shared__ double zvrerrg[da_shmemsize];
   __shared__ double zvrerpk[da_shmemsize];
   __shared__ double zvreltb[da_shmemsize];
   __shared__ double zvrelix[da_shmemsize];
   __shared__ double zvrelmi[da_shmemsize];
   __shared__ double zvrelrg[da_shmemsize];
   __shared__ double zvrelpk[da_shmemsize];
   __shared__ double zvimrtb[da_shmemsize];
   __shared__ double zvimrix[da_shmemsize];
   __shared__ double zvimrmi[da_shmemsize];
   __shared__ double zvimrrg[da_shmemsize];
   __shared__ double zvimrpk[da_shmemsize];
   __shared__ double zvimltb[da_shmemsize];
   __shared__ double zvimlix[da_shmemsize];
   __shared__ double zvimlmi[da_shmemsize];
   __shared__ double zvimlrg[da_shmemsize];
   __shared__ double zvimlpk[da_shmemsize];

   double xrrtb,xirtb,yrrtb,yirtb,zrrtb,zirtb,accrtb;
   double xrltb,xiltb,yrltb,yiltb,zrltb,ziltb,accltb;
   double xrrix,xirix,yrrix,yirix,zrrix,zirix,accrix;
   double xrlix,xilix,yrlix,yilix,zrlix,zilix,acclix;
   double xrrmi,xirmi,yrrmi,yirmi,zrrmi,zirmi,accrmi;
   double xrlmi,xilmi,yrlmi,yilmi,zrlmi,zilmi,acclmi;
   double xrrrg,xirrg,yrrrg,yirrg,zrrrg,zirrg,accrrg;
   double xrlrg,xilrg,yrlrg,yilrg,zrlrg,zilrg,acclrg;
   double xrrpk,xirpk,yrrpk,yirpk,zrrpk,zirpk,accrpk;
   double xrlpk,xilpk,yrlpk,yilpk,zrlpk,zilpk,acclpk;
   int idx;

   xvrertb[k] = xrertb[k]; xvimrtb[k] = ximrtb[k];
   xvreltb[k] = xreltb[k]; xvimltb[k] = ximltb[k];
   xvrerix[k] = xrerix[k]; xvimrix[k] = ximrix[k];
   xvrelix[k] = xrelix[k]; xvimlix[k] = ximlix[k];
   xvrermi[k] = xrermi[k]; xvimrmi[k] = ximrmi[k];
   xvrelmi[k] = xrelmi[k]; xvimlmi[k] = ximlmi[k];
   xvrerrg[k] = xrerrg[k]; xvimrrg[k] = ximrrg[k];
   xvrelrg[k] = xrelrg[k]; xvimlrg[k] = ximlrg[k];
   xvrerpk[k] = xrerpk[k]; xvimrpk[k] = ximrpk[k];
   xvrelpk[k] = xrelpk[k]; xvimlpk[k] = ximlpk[k];
   yvrertb[k] = yrertb[k]; yvimrtb[k] = yimrtb[k];
   yvreltb[k] = yreltb[k]; yvimltb[k] = yimltb[k];
   yvrerix[k] = yrerix[k]; yvimrix[k] = yimrix[k];
   yvrelix[k] = yrelix[k]; yvimlix[k] = yimlix[k];
   yvrermi[k] = yrermi[k]; yvimrmi[k] = yimrmi[k];
   yvrelmi[k] = yrelmi[k]; yvimlmi[k] = yimlmi[k];
   yvrerrg[k] = yrerrg[k]; yvimrrg[k] = yimrrg[k];
   yvrelrg[k] = yrelrg[k]; yvimlrg[k] = yimlrg[k];
   yvrerpk[k] = yrerpk[k]; yvimrpk[k] = yimrpk[k];
   yvrelpk[k] = yrelpk[k]; yvimlpk[k] = yimlpk[k];

   __syncthreads();

   // z[k] = x[0]*y[k]
   xrrtb = xvrertb[0]; xrrix = xvrerix[0]; xrrmi = xvrermi[0];
   xrrrg = xvrerrg[0]; xrrpk = xvrerpk[0];
   xrltb = xvreltb[0]; xrlix = xvrelix[0]; xrlmi = xvrelmi[0];
   xrlrg = xvrelrg[0]; xrlpk = xvrelpk[0];
   xirtb = xvimrtb[0]; xirix = xvimrix[0]; xirmi = xvimrmi[0];
   xirrg = xvimrrg[0]; xirpk = xvimrpk[0];
   xiltb = xvimltb[0]; xilix = xvimlix[0]; xilmi = xvimlmi[0];
   xilrg = xvimlrg[0]; xilpk = xvimlpk[0];
   yrrtb = yvrertb[k]; yrrix = yvrerix[k]; yrrmi = yvrermi[k];
   yrrrg = yvrerrg[k]; yrrpk = yvrerpk[k];
   yrltb = yvreltb[k]; yrlix = yvrelix[k]; yrlmi = yvrelmi[k];
   yrlrg = yvrelrg[k]; yrlpk = yvrelpk[k];
   yirtb = yvimrtb[k]; yirix = yvimrix[k]; yirmi = yvimrmi[k];
   yirrg = yvimrrg[k]; yirpk = yvimrpk[k];
   yiltb = yvimltb[k]; yilix = yvimlix[k]; yilmi = yvimlmi[k];
   yilrg = yvimlrg[k]; yilpk = yvimlpk[k];

   dag_mul(xrrtb,xrrix,xrrmi,xrrrg,xrrpk,xrltb,xrlix,xrlmi,xrlrg,xrlpk,
           yrrtb,yrrix,yrrmi,yrrrg,yrrpk,yrltb,yrlix,yrlmi,yrlrg,yrlpk,
           &zrrtb,&zrrix,&zrrmi,&zrrrg,&zrrpk,
           &zrltb,&zrlix,&zrlmi,&zrlrg,&zrlpk);         // zr = xr*yr
   dag_mul(xirtb,xirix,xirmi,xirrg,xirpk,xiltb,xilix,xilmi,xilrg,xilpk,
           yirtb,yirix,yirmi,yirrg,yirpk,yiltb,yilix,yilmi,yilrg,yilpk,
           &accrtb,&accrix,&accrmi,&accrrg,&accrpk,
           &accltb,&acclix,&acclmi,&acclrg,&acclpk);    // acc = xi*yi
   dag_minus(&accrtb,&accrix,&accrmi,&accrrg,&accrpk,
             &accltb,&acclix,&acclmi,&acclrg,&acclpk);
   dag_inc(&zrrtb,&zrrix,&zrrmi,&zrrrg,&zrrpk,
           &zrltb,&zrlix,&zrlmi,&zrlrg,&zrlpk,
           accrtb,accrix,accrmi,accrrg,accrpk,
           accltb,acclix,acclmi,acclrg,acclpk);         // zr = xr*yr - xi*yi
   dag_mul(xrrtb,xrrix,xrrmi,xrrrg,xrrpk,xrltb,xrlix,xrlmi,xrlrg,xrlpk,
           yirtb,yirix,yirmi,yirrg,yirpk,yiltb,yilix,yilmi,yilrg,yilpk,
           &zirtb,&zirix,&zirmi,&zirrg,&zirpk,
           &ziltb,&zilix,&zilmi,&zilrg,&zilpk);         // zi = xr*yi
   dag_mul(xirtb,xirix,xirmi,xirrg,xirpk,xiltb,xilix,xilmi,xilrg,xilpk,
           yrrtb,yrrix,yrrmi,yrrrg,yrrpk,yrltb,yrlix,yrlmi,yrlrg,yrlpk,
           &accrtb,&accrix,&accrmi,&accrrg,&accrpk,
           &accltb,&acclix,&acclmi,&acclrg,&acclpk);    // acc = xi*yr
   dag_inc(&zirtb,&zirix,&zirmi,&zirrg,&zirpk,
           &ziltb,&zilix,&zilmi,&zilrg,&zilpk,
           accrtb,accrix,accrmi,accrrg,accrpk,
           accltb,acclix,acclmi,acclrg,acclpk);         // zr = xr*yr + xi*yi

   zvrertb[k] = zrrtb; zvrerix[k] = zrrix; zvrermi[k] = zrrmi;
   zvrerrg[k] = zrrrg; zvrerpk[k] = zrrpk;
   zvreltb[k] = zrltb; zvrelix[k] = zrlix; zvrelmi[k] = zrlmi;
   zvrelrg[k] = zrlrg; zvrelpk[k] = zrlpk;
   zvimrtb[k] = zirtb; zvimrix[k] = zirix; zvimrmi[k] = zirmi;
   zvimrrg[k] = zirrg; zvimrpk[k] = zirpk;
   zvimltb[k] = ziltb; zvimlix[k] = zilix; zvimlmi[k] = zilmi;
   zvimlrg[k] = zilrg; zvimlpk[k] = zilpk;

   for(int i=1; i<=k; i++) // z[k] = z[k] + x[i]*y[k-i]
   {
      xrrtb = xvrertb[i]; xrrix = xvrerix[i]; xrrmi = xvrermi[i];
      xrrrg = xvrerrg[i]; xrrpk = xvrerpk[i];
      xrltb = xvreltb[i]; xrlix = xvrelix[i]; xrlmi = xvrelmi[i];
      xrlrg = xvrelrg[i]; xrlpk = xvrelpk[i];
      xirtb = xvimrtb[i]; xirix = xvimrix[i]; xirmi = xvimrmi[i];
      xirrg = xvimrrg[i]; xirpk = xvimrpk[i];
      xiltb = xvimltb[i]; xilix = xvimlix[i]; xilmi = xvimlmi[i];
      xilrg = xvimlrg[i]; xilpk = xvimlpk[i];
      idx = k-i;
      yrrtb = yvrertb[idx]; yrrix = yvrerix[idx]; yrrmi = yvrermi[idx];
      yrrrg = yvrerrg[idx]; yrrpk = yvrerpk[idx];
      yrltb = yvreltb[idx]; yrlix = yvrelix[idx]; yrlmi = yvrelmi[idx];
      yrlrg = yvrelrg[idx]; yrlpk = yvrelpk[idx];
      yirtb = yvimrtb[idx]; yirix = yvimrix[idx]; yirmi = yvimrmi[idx];
      yirrg = yvimrrg[idx]; yirpk = yvimrpk[idx];
      yiltb = yvimltb[idx]; yilix = yvimlix[idx]; yilmi = yvimlmi[idx];
      yilrg = yvimlrg[idx]; yilpk = yvimlpk[idx];

      dag_mul(xrrtb,xrrix,xrrmi,xrrrg,xrrpk,xrltb,xrlix,xrlmi,xrlrg,xrlpk,
              yrrtb,yrrix,yrrmi,yrrrg,yrrpk, yrltb,yrlix,yrlmi,yrlrg,yrlpk,
              &zrrtb,&zrrix,&zrrmi,&zrrrg,&zrrpk,
              &zrltb,&zrlix,&zrlmi,&zrlrg,&zrlpk);        // zr = xr*yr
      dag_mul(xirtb,xirix,xirmi,xirrg,xirpk,xiltb,xilix,xilmi,xilrg,xilpk,
              yirtb,yirix,yirmi,yirrg,yirpk,yiltb,yilix,yilmi,yilrg,yilpk,
              &accrtb,&accrix,&accrmi,&accrrg,&accrpk,
              &accltb,&acclix,&acclmi,&acclrg,&acclpk);   // xi*yi
      dag_minus(&accrtb,&accrix,&accrmi,&accrrg,&accrpk,
                &accltb,&acclix,&acclmi,&acclrg,&acclpk);
      dag_inc(&zrrtb,&zrrix,&zrrmi,&zrrrg,&zrrpk,
              &zrltb,&zrlix,&zrlmi,&zrlrg,&zrlpk,
              accrtb,accrix,accrmi,accrrg,accrpk,
              accltb,acclix,acclmi,acclrg,acclpk);     // zr = xr*yr - xi*yi
      dag_mul(xrrtb,xrrix,xrrmi,xrrrg,xrrpk,xrltb,xrlix,xrlmi,xrlrg,xrlpk,
              yirtb,yirix,yirmi,yirrg,yirpk,yiltb,yilix,yilmi,yilrg,yilpk,
              &zirtb,&zirix,&zirmi,&zirrg,&zirpk,
              &ziltb,&zilix,&zilmi,&zilrg,&zilpk);        // zi = xr*yi
      dag_mul(xirtb,xirix,xirmi,xirrg,xirpk,xiltb,xilix,xilmi,xilrg,xilpk,
              yrrtb,yrrix,yrrmi,yrrrg,yrrpk,yrltb,yrlix,yrlmi,yrlrg,yrlpk,
              &accrtb,&accrix,&accrmi,&accrrg,&accrpk,
              &accltb,&acclix,&acclmi,&acclrg,&acclpk);   // xi*yr
      dag_inc(&zirtb,&zirix,&zirmi,&zirrg,&zirpk,
              &ziltb,&zilix,&zilmi,&zilrg,&zilpk,
              accrtb,accrix,accrmi,accrrg,accrpk,
              accltb,acclix,acclmi,acclrg,acclpk);     // zr = xr*yi + xi*yr
      // zvre[k] += zr; zvim[k] += zi
      dag_inc(&zvrertb[k],&zvrerix[k],&zvrermi[k],&zvrerrg[k],&zvrerpk[k],
              &zvreltb[k],&zvrelix[k],&zvrelmi[k],&zvrelrg[k],&zvrelpk[k],
              zrrtb,zrrix,zrrmi,zrrrg,zrrpk,zrltb,zrlix,zrlmi,zrlrg,zrlpk);
      dag_inc(&zvimrtb[k],&zvimrix[k],&zvimrmi[k],&zvimrrg[k],&zvimrpk[k],
              &zvimltb[k],&zvimlix[k],&zvimlmi[k],&zvimlrg[k],&zvimlpk[k],
              zirtb,zirix,zirmi,zirrg,zirpk,ziltb,zilix,zilmi,zilrg,zilpk);
   }

   __syncthreads();

   zrertb[k] = zvrertb[k]; zrerix[k] = zvrerix[k]; zrermi[k] = zvrermi[k];
   zrerrg[k] = zvrerrg[k]; zrerpk[k] = zvrerpk[k];
   zreltb[k] = zvreltb[k]; zrelix[k] = zvrelix[k]; zrelmi[k] = zvrelmi[k];
   zrelrg[k] = zvrelrg[k]; zrelpk[k] = zvrelpk[k];
   zimrtb[k] = zvimrtb[k]; zimrix[k] = zvimrix[k]; zimrmi[k] = zvimrmi[k];
   zimrrg[k] = zvimrrg[k]; zimrpk[k] = zvimrpk[k];
   zimltb[k] = zvimltb[k]; zimlix[k] = zvimlix[k]; zimlmi[k] = zvimlmi[k];
   zimlrg[k] = zvimlrg[k]; zimlpk[k] = zvimlpk[k];
}

__global__ void cmplx10_padded_convolute
 ( double *xrertb, double *xrerix, double *xrermi, double *xrerrg,
   double *xrerpk, double *xreltb, double *xrelix, double *xrelmi,
   double *xrelrg, double *xrelpk,
   double *ximrtb, double *ximrix, double *ximrmi, double *ximrrg,
   double *ximrpk, double *ximltb, double *ximlix, double *ximlmi,
   double *ximlrg, double *ximlpk,
   double *yrertb, double *yrerix, double *yrermi, double *yrerrg,
   double *yrerpk, double *yreltb, double *yrelix, double *yrelmi,
   double *yrelrg, double *yrelpk,
   double *yimrtb, double *yimrix, double *yimrmi, double *yimrrg,
   double *yimrpk, double *yimltb, double *yimlix, double *yimlmi,
   double *yimlrg, double *yimlpk,
   double *zrertb, double *zrerix, double *zrermi, double *zrerrg,
   double *zrerpk, double *zreltb, double *zrelix, double *zrelmi,
   double *zrelrg, double *zrelpk,
   double *zimrtb, double *zimrix, double *zimrmi, double *zimrrg,
   double *zimrpk, double *zimltb, double *zimlix, double *zimlmi,
   double *zimlrg, double *zimlpk, int dim )
{
   int k = threadIdx.x;       // thread k computes zre[k] and zim[k]

   __shared__ double xvrertb[da_shmemsize];
   __shared__ double xvrerix[da_shmemsize];
   __shared__ double xvrermi[da_shmemsize];
   __shared__ double xvrerrg[da_shmemsize];
   __shared__ double xvrerpk[da_shmemsize];
   __shared__ double xvreltb[da_shmemsize];
   __shared__ double xvrelix[da_shmemsize];
   __shared__ double xvrelmi[da_shmemsize];
   __shared__ double xvrelrg[da_shmemsize];
   __shared__ double xvrelpk[da_shmemsize];
   __shared__ double xvimrtb[da_shmemsize];
   __shared__ double xvimrix[da_shmemsize];
   __shared__ double xvimrmi[da_shmemsize];
   __shared__ double xvimrrg[da_shmemsize];
   __shared__ double xvimrpk[da_shmemsize];
   __shared__ double xvimltb[da_shmemsize];
   __shared__ double xvimlix[da_shmemsize];
   __shared__ double xvimlmi[da_shmemsize];
   __shared__ double xvimlrg[da_shmemsize];
   __shared__ double xvimlpk[da_shmemsize];
   __shared__ double yvrertb[da_shmemsize];
   __shared__ double yvrerix[da_shmemsize];
   __shared__ double yvrermi[da_shmemsize];
   __shared__ double yvrerrg[da_shmemsize];
   __shared__ double yvrerpk[da_shmemsize];
   __shared__ double yvreltb[da_shmemsize];
   __shared__ double yvrelix[da_shmemsize];
   __shared__ double yvrelmi[da_shmemsize];
   __shared__ double yvrelrg[da_shmemsize];
   __shared__ double yvrelpk[da_shmemsize];
   __shared__ double yvimrtb[da_shmemsize];
   __shared__ double yvimrix[da_shmemsize];
   __shared__ double yvimrmi[da_shmemsize];
   __shared__ double yvimrrg[da_shmemsize];
   __shared__ double yvimrpk[da_shmemsize];
   __shared__ double yvimltb[da_shmemsize];
   __shared__ double yvimlix[da_shmemsize];
   __shared__ double yvimlmi[da_shmemsize];
   __shared__ double yvimlrg[da_shmemsize];
   __shared__ double yvimlpk[da_shmemsize];
   __shared__ double zvrertb[da_shmemsize];
   __shared__ double zvrerix[da_shmemsize];
   __shared__ double zvrermi[da_shmemsize];
   __shared__ double zvrerrg[da_shmemsize];
   __shared__ double zvrerpk[da_shmemsize];
   __shared__ double zvreltb[da_shmemsize];
   __shared__ double zvrelix[da_shmemsize];
   __shared__ double zvrelmi[da_shmemsize];
   __shared__ double zvrelrg[da_shmemsize];
   __shared__ double zvrelpk[da_shmemsize];
   __shared__ double zvimrtb[da_shmemsize];
   __shared__ double zvimrix[da_shmemsize];
   __shared__ double zvimrmi[da_shmemsize];
   __shared__ double zvimrrg[da_shmemsize];
   __shared__ double zvimrpk[da_shmemsize];
   __shared__ double zvimltb[da_shmemsize];
   __shared__ double zvimlix[da_shmemsize];
   __shared__ double zvimlmi[da_shmemsize];
   __shared__ double zvimlrg[da_shmemsize];
   __shared__ double zvimlpk[da_shmemsize];

   double xrrtb,xirtb,yrrtb,yirtb,zrrtb,zirtb,accrtb;
   double xrltb,xiltb,yrltb,yiltb,zrltb,ziltb,accltb;
   double xrrix,xirix,yrrix,yirix,zrrix,zirix,accrix;
   double xrlix,xilix,yrlix,yilix,zrlix,zilix,acclix;
   double xrrmi,xirmi,yrrmi,yirmi,zrrmi,zirmi,accrmi;
   double xrlmi,xilmi,yrlmi,yilmi,zrlmi,zilmi,acclmi;
   double xrrrg,xirrg,yrrrg,yirrg,zrrrg,zirrg,accrrg;
   double xrlrg,xilrg,yrlrg,yilrg,zrlrg,zilrg,acclrg;
   double xrrpk,xirpk,yrrpk,yirpk,zrrpk,zirpk,accrpk;
   double xrlpk,xilpk,yrlpk,yilpk,zrlpk,zilpk,acclpk;
   int idx = dim+k;

   xvrertb[k] = xrertb[k]; xvimrtb[k] = ximrtb[k];
   xvreltb[k] = xreltb[k]; xvimltb[k] = ximltb[k];
   xvrerix[k] = xrerix[k]; xvimrix[k] = ximrix[k];
   xvrelix[k] = xrelix[k]; xvimlix[k] = ximlix[k];
   xvrermi[k] = xrermi[k]; xvimrmi[k] = ximrmi[k];
   xvrelmi[k] = xrelmi[k]; xvimlmi[k] = ximlmi[k];
   xvrerrg[k] = xrerrg[k]; xvimrrg[k] = ximrrg[k];
   xvrelrg[k] = xrelrg[k]; xvimlrg[k] = ximlrg[k];
   xvrerpk[k] = xrerpk[k]; xvimrpk[k] = ximrpk[k];
   xvrelpk[k] = xrelpk[k]; xvimlpk[k] = ximlpk[k];
   yvrertb[k] = 0.0; yvimrtb[k] = 0.0;
   yvreltb[k] = 0.0; yvimltb[k] = 0.0;
   yvrerix[k] = 0.0; yvimrix[k] = 0.0;
   yvrelix[k] = 0.0; yvimlix[k] = 0.0;
   yvrermi[k] = 0.0; yvimrmi[k] = 0.0;
   yvrelmi[k] = 0.0; yvimlmi[k] = 0.0;
   yvrerrg[k] = 0.0; yvimrrg[k] = 0.0;
   yvrelrg[k] = 0.0; yvimlrg[k] = 0.0;
   yvrerpk[k] = 0.0; yvimrpk[k] = 0.0;
   yvrelpk[k] = 0.0; yvimlpk[k] = 0.0;
   yvrertb[idx] = yrertb[k]; yvimrtb[idx] = yimrtb[k];
   yvreltb[idx] = yreltb[k]; yvimltb[idx] = yimltb[k];
   yvrerix[idx] = yrerix[k]; yvimrix[idx] = yimrix[k];
   yvrelix[idx] = yrelix[k]; yvimlix[idx] = yimlix[k];
   yvrermi[idx] = yrermi[k]; yvimrmi[idx] = yimrmi[k];
   yvrelmi[idx] = yrelmi[k]; yvimlmi[idx] = yimlmi[k];
   yvrerrg[idx] = yrerrg[k]; yvimrrg[idx] = yimrrg[k];
   yvrelrg[idx] = yrelrg[k]; yvimlrg[idx] = yimlrg[k];
   yvrerpk[idx] = yrerpk[k]; yvimrpk[idx] = yimrpk[k];
   yvrelpk[idx] = yrelpk[k]; yvimlpk[idx] = yimlpk[k];

   __syncthreads();

   // z[k] = x[0]*y[k]
   xrrtb = xvrertb[0]; xrrix = xvrerix[0]; xrrmi = xvrermi[0];
   xrrrg = xvrerrg[0]; xrrpk = xvrerpk[0];
   xrltb = xvreltb[0]; xrlix = xvrelix[0]; xrlmi = xvrelmi[0];
   xrlrg = xvrelrg[0]; xrlpk = xvrelpk[0];
   xirtb = xvimrtb[0]; xirix = xvimrix[0]; xirmi = xvimrmi[0];
   xirrg = xvimrrg[0]; xirpk = xvimrpk[0];
   xiltb = xvimltb[0]; xilix = xvimlix[0]; xilmi = xvimlmi[0];
   xilrg = xvimlrg[0]; xilpk = xvimlpk[0];
   yrrtb = yvrertb[idx]; yrrix = yvrerix[idx]; yrrmi = yvrermi[idx];
   yrrrg = yvrerrg[idx]; yrrpk = yvrerpk[idx];
   yrltb = yvreltb[idx]; yrlix = yvrelix[idx]; yrlmi = yvrelmi[idx];
   yrlrg = yvrelrg[idx]; yrlpk = yvrelpk[idx];
   yirtb = yvimrtb[idx]; yirix = yvimrix[idx]; yirmi = yvimrmi[idx];
   yirrg = yvimrrg[idx]; yirpk = yvimrpk[idx];
   yiltb = yvimltb[idx]; yilix = yvimlix[idx]; yilmi = yvimlmi[idx];
   yilrg = yvimlrg[idx]; yilpk = yvimlpk[idx];

   dag_mul(xrrtb,xrrix,xrrmi,xrrrg,xrrpk,xrltb,xrlix,xrlmi,xrlrg,xrlpk,
           yrrtb,yrrix,yrrmi,yrrrg,yrrpk,yrltb,yrlix,yrlmi,yrlrg,yrlpk,
           &zrrtb,&zrrix,&zrrmi,&zrrrg,&zrrpk,
           &zrltb,&zrlix,&zrlmi,&zrlrg,&zrlpk);         // zr = xr*yr
   __syncthreads();
   dag_mul(xirtb,xirix,xirmi,xirrg,xirpk,xiltb,xilix,xilmi,xilrg,xilpk,
           yirtb,yirix,yirmi,yirrg,yirpk,yiltb,yilix,yilmi,yilrg,yilpk,
           &accrtb,&accrix,&accrmi,&accrrg,&accrpk,
           &accltb,&acclix,&acclmi,&acclrg,&acclpk);    // acc = xi*yi
   __syncthreads();
   dag_minus(&accrtb,&accrix,&accrmi,&accrrg,&accrpk,
             &accltb,&acclix,&acclmi,&acclrg,&acclpk);
   __syncthreads();
   dag_inc(&zrrtb,&zrrix,&zrrmi,&zrrrg,&zrrpk,
           &zrltb,&zrlix,&zrlmi,&zrlrg,&zrlpk,
           accrtb,accrix,accrmi,accrrg,accrpk,
           accltb,acclix,acclmi,acclrg,acclpk);         // zr = xr*yr - xi*yi
   __syncthreads();
   dag_mul(xrrtb,xrrix,xrrmi,xrrrg,xrrpk,xrltb,xrlix,xrlmi,xrlrg,xrlpk,
           yirtb,yirix,yirmi,yirrg,yirpk,yiltb,yilix,yilmi,yilrg,yilpk,
           &zirtb,&zirix,&zirmi,&zirrg,&zirpk,
           &ziltb,&zilix,&zilmi,&zilrg,&zilpk);         // zi = xr*yi
   __syncthreads();
   dag_mul(xirtb,xirix,xirmi,xirrg,xirpk,xiltb,xilix,xilmi,xilrg,xilpk,
           yrrtb,yrrix,yrrmi,yrrrg,yrrpk,yrltb,yrlix,yrlmi,yrlrg,yrlpk,
           &accrtb,&accrix,&accrmi,&accrrg,&accrpk,
           &accltb,&acclix,&acclmi,&acclrg,&acclpk);    // acc = xi*yr
   __syncthreads();
   dag_inc(&zirtb,&zirix,&zirmi,&zirrg,&zirpk,
           &ziltb,&zilix,&zilmi,&zilrg,&zilpk,
           accrtb,accrix,accrmi,accrrg,accrpk,
           accltb,acclix,acclmi,acclrg,acclpk);         // zr = xr*yr + xi*yi
   __syncthreads();

   zvrertb[k] = zrrtb; zvrerix[k] = zrrix; zvrermi[k] = zrrmi;
   zvrerrg[k] = zrrrg; zvrerpk[k] = zrrpk;
   zvreltb[k] = zrltb; zvrelix[k] = zrlix; zvrelmi[k] = zrlmi;
   zvrelrg[k] = zrlrg; zvrelpk[k] = zrlpk;
   zvimrtb[k] = zirtb; zvimrix[k] = zirix; zvimrmi[k] = zirmi;
   zvimrrg[k] = zirrg; zvimrpk[k] = zirpk;
   zvimltb[k] = ziltb; zvimlix[k] = zilix; zvimlmi[k] = zilmi;
   zvimlrg[k] = zilrg; zvimlpk[k] = zilpk;

   for(int i=1; i<dim; i++) // z[k] = z[k] + x[i]*y[k-i]
   {
      idx = dim + k - i;
      xrrtb = xvrertb[i]; xrrix = xvrerix[i]; xrrmi = xvrermi[i];
      xrrrg = xvrerrg[i]; xrrpk = xvrerpk[i];
      xrltb = xvreltb[i]; xrlix = xvrelix[i]; xrlmi = xvrelmi[i];
      xrlrg = xvrelrg[i]; xrlpk = xvrelpk[i];
      xirtb = xvimrtb[i]; xirix = xvimrix[i]; xirmi = xvimrmi[i];
      xirrg = xvimrrg[i]; xirpk = xvimrpk[i];
      xiltb = xvimltb[i]; xilix = xvimlix[i]; xilmi = xvimlmi[i];
      xilrg = xvimlrg[i]; xilpk = xvimlpk[i];
      yrrtb = yvrertb[idx]; yrrix = yvrerix[idx]; yrrmi = yvrermi[idx];
      yrrrg = yvrerrg[idx]; yrrpk = yvrerpk[idx];
      yrltb = yvreltb[idx]; yrlix = yvrelix[idx]; yrlmi = yvrelmi[idx];
      yrlrg = yvrelrg[idx]; yrlpk = yvrelpk[idx];
      yirtb = yvimrtb[idx]; yirix = yvimrix[idx]; yirmi = yvimrmi[idx];
      yirrg = yvimrrg[idx]; yirpk = yvimrpk[idx];
      yiltb = yvimltb[idx]; yilix = yvimlix[idx]; yilmi = yvimlmi[idx];
      yilrg = yvimlrg[idx]; yilpk = yvimlpk[idx];

      dag_mul(xrrtb,xrrix,xrrmi,xrrrg,xrrpk,xrltb,xrlix,xrlmi,xrlrg,xrlpk,
              yrrtb,yrrix,yrrmi,yrrrg,yrrpk, yrltb,yrlix,yrlmi,yrlrg,yrlpk,
              &zrrtb,&zrrix,&zrrmi,&zrrrg,&zrrpk,
              &zrltb,&zrlix,&zrlmi,&zrlrg,&zrlpk);        // zr = xr*yr
      __syncthreads();
      dag_mul(xirtb,xirix,xirmi,xirrg,xirpk,xiltb,xilix,xilmi,xilrg,xilpk,
              yirtb,yirix,yirmi,yirrg,yirpk,yiltb,yilix,yilmi,yilrg,yilpk,
              &accrtb,&accrix,&accrmi,&accrrg,&accrpk,
              &accltb,&acclix,&acclmi,&acclrg,&acclpk);   // xi*yi
      __syncthreads();
      dag_minus(&accrtb,&accrix,&accrmi,&accrrg,&accrpk,
                &accltb,&acclix,&acclmi,&acclrg,&acclpk);
      __syncthreads();
      dag_inc(&zrrtb,&zrrix,&zrrmi,&zrrrg,&zrrpk,
              &zrltb,&zrlix,&zrlmi,&zrlrg,&zrlpk,
              accrtb,accrix,accrmi,accrrg,accrpk,
              accltb,acclix,acclmi,acclrg,acclpk);     // zr = xr*yr - xi*yi
      __syncthreads();
      dag_mul(xrrtb,xrrix,xrrmi,xrrrg,xrrpk,xrltb,xrlix,xrlmi,xrlrg,xrlpk,
              yirtb,yirix,yirmi,yirrg,yirpk,yiltb,yilix,yilmi,yilrg,yilpk,
              &zirtb,&zirix,&zirmi,&zirrg,&zirpk,
              &ziltb,&zilix,&zilmi,&zilrg,&zilpk);        // zi = xr*yi
      __syncthreads();
      dag_mul(xirtb,xirix,xirmi,xirrg,xirpk,xiltb,xilix,xilmi,xilrg,xilpk,
              yrrtb,yrrix,yrrmi,yrrrg,yrrpk,yrltb,yrlix,yrlmi,yrlrg,yrlpk,
              &accrtb,&accrix,&accrmi,&accrrg,&accrpk,
              &accltb,&acclix,&acclmi,&acclrg,&acclpk);   // xi*yr
      __syncthreads();
      dag_inc(&zirtb,&zirix,&zirmi,&zirrg,&zirpk,
              &ziltb,&zilix,&zilmi,&zilrg,&zilpk,
              accrtb,accrix,accrmi,accrrg,accrpk,
              accltb,acclix,acclmi,acclrg,acclpk);     // zr = xr*yi + xi*yr
      __syncthreads();
      // zvre[k] += zr; zvim[k] += zi
      dag_inc(&zvrertb[k],&zvrerix[k],&zvrermi[k],&zvrerrg[k],&zvrerpk[k],
              &zvreltb[k],&zvrelix[k],&zvrelmi[k],&zvrelrg[k],&zvrelpk[k],
              zrrtb,zrrix,zrrmi,zrrrg,zrrpk,zrltb,zrlix,zrlmi,zrlrg,zrlpk);
      __syncthreads();
      dag_inc(&zvimrtb[k],&zvimrix[k],&zvimrmi[k],&zvimrrg[k],&zvimrpk[k],
              &zvimltb[k],&zvimlix[k],&zvimlmi[k],&zvimlrg[k],&zvimlpk[k],
              zirtb,zirix,zirmi,zirrg,zirpk,ziltb,zilix,zilmi,zilrg,zilpk);
      __syncthreads();
   }
   zrertb[k] = zvrertb[k]; zrerix[k] = zvrerix[k]; zrermi[k] = zvrermi[k];
   zrerrg[k] = zvrerrg[k]; zrerpk[k] = zvrerpk[k];
   zreltb[k] = zvreltb[k]; zrelix[k] = zvrelix[k]; zrelmi[k] = zvrelmi[k];
   zrelrg[k] = zvrelrg[k]; zrelpk[k] = zvrelpk[k];
   zimrtb[k] = zvimrtb[k]; zimrix[k] = zvimrix[k]; zimrmi[k] = zvimrmi[k];
   zimrrg[k] = zvimrrg[k]; zimrpk[k] = zvimrpk[k];
   zimltb[k] = zvimltb[k]; zimlix[k] = zvimlix[k]; zimlmi[k] = zvimlmi[k];
   zimlrg[k] = zvimlrg[k]; zimlpk[k] = zvimlpk[k];
}

void GPU_dbl10_product
 ( double *xrtb_h, double *xrix_h, double *xrmi_h, double *xrrg_h,
   double *xrpk_h, double *xltb_h, double *xlix_h, double *xlmi_h,
   double *xlrg_h, double *xlpk_h,
   double *yrtb_h, double *yrix_h, double *yrmi_h, double *yrrg_h,
   double *yrpk_h, double *yltb_h, double *ylix_h, double *ylmi_h,
   double *ylrg_h, double *ylpk_h,
   double *zrtb_h, double *zrix_h, double *zrmi_h, double *zrrg_h,
   double *zrpk_h, double *zltb_h, double *zlix_h, double *zlmi_h,
   double *zlrg_h, double *zlpk_h, int deg, int freq, int BS, int padded )
{
   const int dim = deg+1;            // length of all vectors
   double* xrtb_d;                   // xrtb_d is xrtb_h on the device
   double* xrix_d;                   // xrix_d is xrix_h on the device
   double* xrmi_d;                   // xrmi_d is xrmi_h on the device
   double* xrrg_d;                   // xrrg_d is xrrg_h on the device
   double* xrpk_d;                   // xrpk_d is xrpk_h on the device
   double* xltb_d;                   // xltb_d is xltb_h on the device
   double* xlix_d;                   // xlix_d is xlix_h on the device
   double* xlmi_d;                   // xlmi_d is xlmi_h on the device
   double* xlrg_d;                   // xlrg_d is xlrg_h on the device
   double* xlpk_d;                   // xlpk_d is xlpk_h on the device
   double* yrtb_d;                   // yrtb_d is yrtb_h on the device
   double* yrix_d;                   // yrix_d is yrix_h on the device
   double* yrmi_d;                   // yrmi_d is yrmi_h on the device
   double* yrrg_d;                   // yrrg_d is yrrg_h on the device
   double* yrpk_d;                   // yrpk_d is yrpk_h on the device
   double* yltb_d;                   // yltb_d is yltb_h on the device
   double* ylix_d;                   // ylix_d is ylix_h on the device
   double* ylmi_d;                   // ylmi_d is ylmi_h on the device
   double* ylrg_d;                   // ylrg_d is ylrg_h on the device
   double* ylpk_d;                   // ylpk_d is ylpk_h on the device
   double* zrtb_d;                   // zrtb_d is zrtb_h on the device
   double* zrix_d;                   // zrix_d is zrix_h on the device
   double* zrmi_d;                   // zrmi_d is zrmi_h on the device
   double* zrrg_d;                   // zrrg_d is zrrg_h on the device
   double* zrpk_d;                   // zrpk_d is zrpk_h on the device
   double* zltb_d;                   // zltb_d is zltb_h on the device
   double* zlix_d;                   // zlix_d is zlix_h on the device
   double* zlmi_d;                   // zlmi_d is zlmi_h on the device
   double* zlrg_d;                   // zlrg_d is zlrg_h on the device
   double* zlpk_d;                   // zlpk_d is zlpk_h on the device
   size_t size = dim*sizeof(double); // number of bytes for each vector

   cudaMalloc((void**)&xrtb_d,size);
   cudaMalloc((void**)&xrix_d,size);
   cudaMalloc((void**)&xrmi_d,size);
   cudaMalloc((void**)&xrrg_d,size);
   cudaMalloc((void**)&xrpk_d,size);
   cudaMalloc((void**)&xltb_d,size);
   cudaMalloc((void**)&xlix_d,size);
   cudaMalloc((void**)&xlmi_d,size);
   cudaMalloc((void**)&xlrg_d,size);
   cudaMalloc((void**)&xlpk_d,size);
   cudaMalloc((void**)&yrtb_d,size);
   cudaMalloc((void**)&yrix_d,size);
   cudaMalloc((void**)&yrmi_d,size);
   cudaMalloc((void**)&yrrg_d,size);
   cudaMalloc((void**)&yrpk_d,size);
   cudaMalloc((void**)&yltb_d,size);
   cudaMalloc((void**)&ylix_d,size);
   cudaMalloc((void**)&ylmi_d,size);
   cudaMalloc((void**)&ylrg_d,size);
   cudaMalloc((void**)&ylpk_d,size);
   cudaMalloc((void**)&zrtb_d,size);
   cudaMalloc((void**)&zrix_d,size);
   cudaMalloc((void**)&zrmi_d,size);
   cudaMalloc((void**)&zrrg_d,size);
   cudaMalloc((void**)&zrpk_d,size);
   cudaMalloc((void**)&zltb_d,size);
   cudaMalloc((void**)&zlix_d,size);
   cudaMalloc((void**)&zlmi_d,size);
   cudaMalloc((void**)&zlrg_d,size);
   cudaMalloc((void**)&zlpk_d,size);
   cudaMemcpy(xrtb_d,xrtb_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(xrix_d,xrix_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(xrmi_d,xrmi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(xrrg_d,xrrg_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(xrpk_d,xrpk_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(xltb_d,xltb_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(xlix_d,xlix_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(xlmi_d,xlmi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(xlrg_d,xlrg_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(xlpk_d,xlpk_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yrtb_d,yrtb_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yrix_d,yrix_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yrmi_d,yrmi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yrrg_d,yrrg_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yrpk_d,yrpk_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yltb_d,yltb_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(ylix_d,ylix_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(ylmi_d,ylmi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(ylrg_d,ylrg_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(ylpk_d,ylpk_h,size,cudaMemcpyHostToDevice);

   if(dim == BS)
   {
      if(padded == 1)
      {
         for(int i=0; i<freq; i++)
            dbl10_padded_convolute<<<1,BS>>>
               (xrtb_d,xrix_d,xrmi_d,xrrg_d,xrpk_d,
                xltb_d,xlix_d,xlmi_d,xlrg_d,xlpk_d,
                yrtb_d,yrix_d,yrmi_d,yrrg_d,yrpk_d,
                yltb_d,ylix_d,ylmi_d,ylrg_d,ylpk_d,
                zrtb_d,zrix_d,zrmi_d,zrrg_d,zrpk_d,
                zltb_d,zlix_d,zlmi_d,zlrg_d,zlpk_d,dim);
      }
      else
      {
         for(int i=0; i<freq; i++)
            dbl10_convolute<<<1,BS>>>
               (xrtb_d,xrix_d,xrmi_d,xrrg_d,xrpk_d,
                xltb_d,xlix_d,xlmi_d,xlrg_d,xlpk_d,
                yrtb_d,yrix_d,yrmi_d,yrrg_d,yrpk_d,
                yltb_d,ylix_d,ylmi_d,ylrg_d,ylpk_d,
                zrtb_d,zrix_d,zrmi_d,zrrg_d,zrpk_d,
                zltb_d,zlix_d,zlmi_d,zlrg_d,zlpk_d,dim);
      }
   }
   cudaMemcpy(zrtb_h,zrtb_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zrix_h,zrix_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zrmi_h,zrmi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zrrg_h,zrrg_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zrpk_h,zrpk_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zltb_h,zltb_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zlix_h,zlix_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zlmi_h,zlmi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zlrg_h,zlrg_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zlpk_h,zlpk_d,size,cudaMemcpyDeviceToHost);
}

void GPU_cmplx10_product
 ( double *xrertb_h, double *xrerix_h, double *xrermi_h, double *xrerrg_h,
   double *xrerpk_h, double *xreltb_h, double *xrelix_h, double *xrelmi_h,
   double *xrelrg_h, double *xrelpk_h,
   double *ximrtb_h, double *ximrix_h, double *ximrmi_h, double *ximrrg_h,
   double *ximrpk_h, double *ximltb_h, double *ximlix_h, double *ximlmi_h,
   double *ximlrg_h, double *ximlpk_h,
   double *yrertb_h, double *yrerix_h, double *yrermi_h, double *yrerrg_h,
   double *yrerpk_h, double *yreltb_h, double *yrelix_h, double *yrelmi_h,
   double *yrelrg_h, double *yrelpk_h,
   double *yimrtb_h, double *yimrix_h, double *yimrmi_h, double *yimrrg_h,
   double *yimrpk_h, double *yimltb_h, double *yimlix_h, double *yimlmi_h,
   double *yimlrg_h, double *yimlpk_h,
   double *zrertb_h, double *zrerix_h, double *zrermi_h, double *zrerrg_h,
   double *zrerpk_h, double *zreltb_h, double *zrelix_h, double *zrelmi_h,
   double *zrelrg_h, double *zrelpk_h,
   double *zimrtb_h, double *zimrix_h, double *zimrmi_h, double *zimrrg_h,
   double *zimrpk_h, double *zimltb_h, double *zimlix_h, double *zimlmi_h,
   double *zimlrg_h, double *zimlpk_h, int deg, int freq, int BS, int mode )
{
   const int dim = deg+1;            // length of all vectors
   double* xrertb_d;                 // xrertb_d is xrertb_h on the device
   double* xrerix_d;                 // xrerix_d is xrerix_h on the device
   double* xrermi_d;                 // xrermi_d is xrermi_h on the device
   double* xrerrg_d;                 // xrerrg_d is xrerrg_h on the device
   double* xrerpk_d;                 // xrerpk_d is xrerpk_h on the device
   double* xreltb_d;                 // xreltb_d is xreltb_h on the device
   double* xrelix_d;                 // xrelix_d is xrelix_h on the device
   double* xrelmi_d;                 // xrelmi_d is xrelmi_h on the device
   double* xrelrg_d;                 // xrelrg_d is xrelrg_h on the device
   double* xrelpk_d;                 // xrelpk_d is xrelpk_h on the device
   double* ximrtb_d;                 // ximrtb_d is ximrtb_h on the device
   double* ximrix_d;                 // ximrix_d is ximrix_h on the device
   double* ximrmi_d;                 // ximrmi_d is ximrmi_h on the device
   double* ximrrg_d;                 // ximrrg_d is ximrrg_h on the device
   double* ximrpk_d;                 // ximrpk_d is ximrpk_h on the device
   double* ximltb_d;                 // ximltb_d is ximltb_h on the device
   double* ximlix_d;                 // ximlix_d is ximlix_h on the device
   double* ximlmi_d;                 // ximlmi_d is ximlmi_h on the device
   double* ximlrg_d;                 // ximlrg_d is ximlrg_h on the device
   double* ximlpk_d;                 // ximlpk_d is ximlpk_h on the device
   double* yrertb_d;                 // yrertb_d is yrertb_h on the device
   double* yrerix_d;                 // yrerix_d is yrerix_h on the device
   double* yrermi_d;                 // yrermi_d is yrermi_h on the device
   double* yrerrg_d;                 // yrerrg_d is yrerrg_h on the device
   double* yrerpk_d;                 // yrerpk_d is yrerpk_h on the device
   double* yreltb_d;                 // yreltb_d is yreltb_h on the device
   double* yrelix_d;                 // yrelix_d is yrelix_h on the device
   double* yrelmi_d;                 // yrelmi_d is yrelmi_h on the device
   double* yrelrg_d;                 // yrelrg_d is yrelrg_h on the device
   double* yrelpk_d;                 // yrelpk_d is yrelpk_h on the device
   double* yimrtb_d;                 // yimrtb_d is yimrtb_h on the device
   double* yimrix_d;                 // yimrix_d is yimrix_h on the device
   double* yimrmi_d;                 // yimrmi_d is yimrmi_h on the device
   double* yimrrg_d;                 // yimrrg_d is yimrrg_h on the device
   double* yimrpk_d;                 // yimrpk_d is yimrpk_h on the device
   double* yimltb_d;                 // yimltb_d is yimltb_h on the device
   double* yimlix_d;                 // yimlix_d is yimlix_h on the device
   double* yimlmi_d;                 // yimlmi_d is yimlmi_h on the device
   double* yimlrg_d;                 // yimlrg_d is yimlrg_h on the device
   double* yimlpk_d;                 // yimlpk_d is yimlpk_h on the device
   double* zrertb_d;                 // zrertb_d is zrertb_h on the device
   double* zrerix_d;                 // zrerix_d is zrerix_h on the device
   double* zrermi_d;                 // zrermi_d is zrermi_h on the device
   double* zrerrg_d;                 // zrerrg_d is zrerrg_h on the device
   double* zrerpk_d;                 // zrerpk_d is zrerpk_h on the device
   double* zreltb_d;                 // zreltb_d is zreltb_h on the device
   double* zrelix_d;                 // zrelix_d is zrelix_h on the device
   double* zrelmi_d;                 // zrelmi_d is zrelmi_h on the device
   double* zrelrg_d;                 // zrelrg_d is zrelrg_h on the device
   double* zrelpk_d;                 // zrelpk_d is zrelpk_h on the device
   double* zimrtb_d;                 // zimrtb_d is zimrtb_h on the device
   double* zimrix_d;                 // zimrix_d is zimrix_h on the device
   double* zimrmi_d;                 // zimrmi_d is zimrmi_h on the device
   double* zimrrg_d;                 // zimrrg_d is zimrrg_h on the device
   double* zimrpk_d;                 // zimrpk_d is zimrpk_h on the device
   double* zimltb_d;                 // zimltb_d is zimltb_h on the device
   double* zimlix_d;                 // zimlix_d is zimlix_h on the device
   double* zimlmi_d;                 // zimlmi_d is zimlmi_h on the device
   double* zimlrg_d;                 // zimlrg_d is zimlrg_h on the device
   double* zimlpk_d;                 // zimlpk_d is zimlpk_h on the device
   double* accrtb_d;                 // accumulates highest doubles
   double* accrix_d;                 // accumulates second highest doubles
   double* accrmi_d;                 // accumulates third highest doubles
   double* accrrg_d;                 // accumulates fourth highest doubles
   double* accrpk_d;                 // accumulates fifth highest doubles
   double* accltb_d;                 // accumulates fifth lowest doubles
   double* acclix_d;                 // accumulates fourth lowest doubles
   double* acclmi_d;                 // accumulates third lowest doubles
   double* acclrg_d;                 // accumulates second lowest doubles
   double* acclpk_d;                 // accumulates lowest doubles
   size_t size = dim*sizeof(double); // number of bytes for each vector

   cudaMalloc((void**)&xrertb_d,size);
   cudaMalloc((void**)&xrerix_d,size);
   cudaMalloc((void**)&xrermi_d,size);
   cudaMalloc((void**)&xrerrg_d,size);
   cudaMalloc((void**)&xrerpk_d,size);
   cudaMalloc((void**)&xreltb_d,size);
   cudaMalloc((void**)&xrelix_d,size);
   cudaMalloc((void**)&xrelmi_d,size);
   cudaMalloc((void**)&xrelrg_d,size);
   cudaMalloc((void**)&xrelpk_d,size);
   cudaMalloc((void**)&ximrtb_d,size);
   cudaMalloc((void**)&ximrix_d,size);
   cudaMalloc((void**)&ximrmi_d,size);
   cudaMalloc((void**)&ximrrg_d,size);
   cudaMalloc((void**)&ximrpk_d,size);
   cudaMalloc((void**)&ximltb_d,size);
   cudaMalloc((void**)&ximlix_d,size);
   cudaMalloc((void**)&ximlmi_d,size);
   cudaMalloc((void**)&ximlrg_d,size);
   cudaMalloc((void**)&ximlpk_d,size);
   cudaMalloc((void**)&yrertb_d,size);
   cudaMalloc((void**)&yrerix_d,size);
   cudaMalloc((void**)&yrermi_d,size);
   cudaMalloc((void**)&yrerrg_d,size);
   cudaMalloc((void**)&yrerpk_d,size);
   cudaMalloc((void**)&yreltb_d,size);
   cudaMalloc((void**)&yrelix_d,size);
   cudaMalloc((void**)&yrelmi_d,size);
   cudaMalloc((void**)&yrelrg_d,size);
   cudaMalloc((void**)&yrelpk_d,size);
   cudaMalloc((void**)&yimrtb_d,size);
   cudaMalloc((void**)&yimrix_d,size);
   cudaMalloc((void**)&yimrmi_d,size);
   cudaMalloc((void**)&yimrrg_d,size);
   cudaMalloc((void**)&yimrpk_d,size);
   cudaMalloc((void**)&yimltb_d,size);
   cudaMalloc((void**)&yimlix_d,size);
   cudaMalloc((void**)&yimlmi_d,size);
   cudaMalloc((void**)&yimlrg_d,size);
   cudaMalloc((void**)&yimlpk_d,size);
   cudaMalloc((void**)&zrertb_d,size);
   cudaMalloc((void**)&zrerix_d,size);
   cudaMalloc((void**)&zrermi_d,size);
   cudaMalloc((void**)&zrerrg_d,size);
   cudaMalloc((void**)&zrerpk_d,size);
   cudaMalloc((void**)&zreltb_d,size);
   cudaMalloc((void**)&zrelix_d,size);
   cudaMalloc((void**)&zrelmi_d,size);
   cudaMalloc((void**)&zrelrg_d,size);
   cudaMalloc((void**)&zrelpk_d,size);
   cudaMalloc((void**)&zimrtb_d,size);
   cudaMalloc((void**)&zimrix_d,size);
   cudaMalloc((void**)&zimrmi_d,size);
   cudaMalloc((void**)&zimrrg_d,size);
   cudaMalloc((void**)&zimrpk_d,size);
   cudaMalloc((void**)&zimltb_d,size);
   cudaMalloc((void**)&zimlix_d,size);
   cudaMalloc((void**)&zimlmi_d,size);
   cudaMalloc((void**)&zimlrg_d,size);
   cudaMalloc((void**)&zimlpk_d,size);
   cudaMalloc((void**)&accrtb_d,size);
   cudaMalloc((void**)&accrix_d,size);
   cudaMalloc((void**)&accrmi_d,size);
   cudaMalloc((void**)&accrrg_d,size);
   cudaMalloc((void**)&accrpk_d,size);
   cudaMalloc((void**)&accltb_d,size);
   cudaMalloc((void**)&acclix_d,size);
   cudaMalloc((void**)&acclmi_d,size);
   cudaMalloc((void**)&acclrg_d,size);
   cudaMalloc((void**)&acclpk_d,size);
   cudaMemcpy(xrertb_d,xrertb_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(xrerix_d,xrerix_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(xrermi_d,xrermi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(xrerrg_d,xrerrg_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(xrerpk_d,xrerpk_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(xreltb_d,xreltb_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(xrelix_d,xrelix_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(xrelmi_d,xrelmi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(xrelrg_d,xrelrg_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(xrelpk_d,xrelpk_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(ximrtb_d,ximrtb_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(ximrix_d,ximrix_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(ximrmi_d,ximrmi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(ximrrg_d,ximrrg_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(ximrpk_d,ximrpk_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(ximltb_d,ximltb_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(ximlix_d,ximlix_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(ximlmi_d,ximlmi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(ximlrg_d,ximlrg_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(ximlpk_d,ximlpk_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yrertb_d,yrertb_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yrerix_d,yrerix_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yrermi_d,yrermi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yrerrg_d,yrerrg_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yrerpk_d,yrerpk_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yreltb_d,yreltb_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yrelix_d,yrelix_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yrelmi_d,yrelmi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yrelrg_d,yrelrg_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yrelpk_d,yrelpk_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yimrtb_d,yimrtb_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yimrix_d,yimrix_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yimrmi_d,yimrmi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yimrrg_d,yimrrg_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yimrpk_d,yimrpk_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yimltb_d,yimltb_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yimlix_d,yimlix_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yimlmi_d,yimlmi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yimlrg_d,yimlrg_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yimlpk_d,yimlpk_h,size,cudaMemcpyHostToDevice);

   if(dim == BS)
   {
      if(mode == 1)
      {
         for(int i=0; i<freq; i++)
            cmplx10_padded_convolute<<<1,BS>>>
               (xrertb_d,xrerix_d,xrermi_d,xrerrg_d,xrerpk_d,
                xreltb_d,xrelix_d,xrelmi_d,xrelrg_d,xrelpk_d,
                ximrtb_d,ximrix_d,ximrmi_d,ximrrg_d,ximrpk_d,
                ximltb_d,ximlix_d,ximlmi_d,ximlrg_d,ximlpk_d,
                yrertb_d,yrerix_d,yrermi_d,yrerrg_d,yrerpk_d,
                yreltb_d,yrelix_d,yrelmi_d,yrelrg_d,yrelpk_d,
                yimrtb_d,yimrix_d,yimrmi_d,yimrrg_d,yimrpk_d,
                yimltb_d,yimlix_d,yimlmi_d,yimlrg_d,yimlpk_d,
                zrertb_d,zrerix_d,zrermi_d,zrerrg_d,zrerpk_d,
                zreltb_d,zrelix_d,zrelmi_d,zrelrg_d,zrelpk_d,
                zimrtb_d,zimrix_d,zimrmi_d,zimrrg_d,zimrpk_d,
                zimltb_d,zimlix_d,zimlmi_d,zimlrg_d,zimlpk_d,dim);
      }
      else if(mode == 1)
      {
         for(int i=0; i<freq; i++)
         {
            dbl10_padded_convolute<<<1,BS>>>
               (xrertb_d,xrerix_d,xrermi_d,xrerrg_d,xrerpk_d,
                xreltb_d,xrelix_d,xrelmi_d,xrelrg_d,xrelpk_d,
                yrertb_d,yrerix_d,yrermi_d,yrerrg_d,yrerpk_d,
                yreltb_d,yrelix_d,yrelmi_d,yrelrg_d,yrelpk_d,
                zrertb_d,zrerix_d,zrermi_d,zrerrg_d,zrerpk_d,
                zreltb_d,zrelix_d,zrelmi_d,zrelrg_d,zrelpk_d,dim);
            dbl10_padded_convolute<<<1,BS>>>
               (ximrtb_d,ximrix_d,ximrmi_d,ximrrg_d,ximrpk_d,
                ximltb_d,ximlix_d,ximlmi_d,ximlrg_d,ximlpk_d,
                yimrtb_d,yimrix_d,yimrmi_d,yimrrg_d,yimrpk_d,
                yimltb_d,yimlix_d,yimlmi_d,yimlrg_d,yimlpk_d,
                accrtb_d,accrix_d,accrmi_d,accrrg_d,zimrpk_d,
                accltb_d,acclix_d,acclmi_d,acclrg_d,zimlpk_d,dim);
            dbl10_decrement<<<1,BS>>>
               (zrertb_d,zrerix_d,zrermi_d,zrerrg_d,zrerpk_d,
                zreltb_d,zrelix_d,zrelmi_d,zrelrg_d,zrelpk_d,
                accrtb_d,accrix_d,accrmi_d,accrrg_d,accrpk_d,
                accltb_d,acclix_d,acclmi_d,acclrg_d,acclpk_d,
                zrertb_d,zrerix_d,zrermi_d,zrerrg_d,zrerpk_d,
                zreltb_d,zrelix_d,zrelmi_d,zrelrg_d,zrelpk_d,dim);
            dbl10_padded_convolute<<<1,BS>>>
               (xrertb_d,xrerix_d,xrermi_d,xrerrg_d,xrerpk_d,
                xreltb_d,xrelix_d,xrelmi_d,xrelrg_d,xrelpk_d,
                yimrtb_d,yimrix_d,yimrmi_d,yimrrg_d,yimrpk_d,
                yimltb_d,yimlix_d,yimlmi_d,yimlrg_d,yimlpk_d,
                zimrtb_d,zimrix_d,zimrmi_d,zimrrg_d,zimrpk_d,
                zimltb_d,zimlix_d,zimlmi_d,zimlrg_d,zimlpk_d,dim);
            dbl10_padded_convolute<<<1,BS>>>
               (ximrtb_d,ximrix_d,ximrmi_d,ximrrg_d,ximrpk_d,
                ximltb_d,ximlix_d,ximlmi_d,ximlrg_d,ximlpk_d,
                yrertb_d,yrerix_d,yrermi_d,yrerrg_d,yrerpk_d,
                yreltb_d,yrelix_d,yrelmi_d,yrelrg_d,yrelpk_d,
                accrtb_d,accrix_d,accrmi_d,accrrg_d,accrpk_d,
                accltb_d,acclix_d,acclmi_d,acclrg_d,acclpk_d,dim);
            dbl10_increment<<<1,BS>>>
               (zimrtb_d,zimrix_d,zimrmi_d,zimrrg_d,zimrpk_d,
                zimltb_d,zimlix_d,zimlmi_d,zimlrg_d,zimlpk_d,
                accrtb_d,accrix_d,accrmi_d,accrrg_d,accrpk_d,
                accltb_d,acclix_d,acclmi_d,acclrg_d,acclpk_d,
                zimrtb_d,zimrix_d,zimrmi_d,zimrrg_d,zimrpk_d,
                zimltb_d,zimlix_d,zimlmi_d,zimlrg_d,zimlpk_d,dim);
         }
      }
      else
      {
         for(int i=0; i<freq; i++)
            cmplx10_convolute<<<1,BS>>>
               (xrertb_d,xrerix_d,xrermi_d,xrerrg_d,xrerpk_d,
                xreltb_d,xrelix_d,xrelmi_d,xrelrg_d,xrelpk_d,
                ximrtb_d,ximrix_d,ximrmi_d,ximrrg_d,ximrpk_d,
                ximltb_d,ximlix_d,ximlmi_d,ximlrg_d,ximlpk_d,
                yrertb_d,yrerix_d,yrermi_d,yrerrg_d,yrerpk_d,
                yreltb_d,yrelix_d,yrelmi_d,yrelrg_d,yrelpk_d,
                yimrtb_d,yimrix_d,yimrmi_d,yimrrg_d,yimrpk_d,
                yimltb_d,yimlix_d,yimlmi_d,yimlrg_d,yimlpk_d,
                zrertb_d,zrerix_d,zrermi_d,zrerrg_d,zrerpk_d,
                zreltb_d,zrelix_d,zrelmi_d,zrelrg_d,zrelpk_d,
                zimrtb_d,zimrix_d,zimrmi_d,zimrrg_d,zimrpk_d,
                zimltb_d,zimlix_d,zimlmi_d,zimlrg_d,zimlpk_d,dim);
      }
   }
   cudaMemcpy(zrertb_h,zrertb_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zrerix_h,zrerix_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zrermi_h,zrermi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zrerrg_h,zrerrg_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zrerpk_h,zrerpk_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zreltb_h,zreltb_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zrelix_h,zrelix_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zrelmi_h,zrelmi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zrelrg_h,zrelrg_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zrelpk_h,zrelpk_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zimrtb_h,zimrtb_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zimrix_h,zimrix_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zimrmi_h,zimrmi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zimrrg_h,zimrrg_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zimrpk_h,zimrpk_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zimltb_h,zimltb_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zimlix_h,zimlix_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zimlmi_h,zimlmi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zimlrg_h,zimlrg_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zimlpk_h,zimlpk_d,size,cudaMemcpyDeviceToHost);
}
