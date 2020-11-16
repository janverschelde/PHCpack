// The file dbl10_monomials_kernels.cu defines the kernels specified
// in dbl10_monomials_kernels.h.

/* The algorithm to compute forward, backward, and cross products
 * (denoted respectively by arrays f, b, and c)
 * for a monomial cff*x[0]*x[1]* .. *x[n-1] goes as follows:
 *
 * f[0] := cff*x[0]
 * for i from 1 to n-1 do f[i] := f[i-1]*x[i]
 * if n > 2 then
 *    b[0] := x[n-1]*x[n-2]
 *    for i from 1 to n-3 do b[i] := b[i-1]*x[n-2-i]
 *    b[n-3] := b[n-3]*cff
 *    if n = 3 then
 *       c[0] = f[0]*x[2]
 *    else
 *       for i from 0 to n-4 do c[i] := f[i]*b[n-4-i]
 *       c[n-3] := f[n-3]*x[n-1]
 *
 * Compared to the evaluation and differentiation of a product of variables,
 * (without coefficient cff), two extra multiplications must be done,
 * but this is better than n+1 multiplications with cff afterwards. */

#include "double_double_gpufun.cu"
#include "quad_double_gpufun.cu"
#include "octo_double_gpufun.cu"
#include "deca_double_gpufun.cu"
#include "dbl10_convolutions_kernels.h"
#include "dbl10_monomials_kernels.h"

__device__ void dbl10_convolute
 ( double *xrtb, double *xrix, double *xrmi, double *xrrg, double *xrpk,
   double *xltb, double *xlix, double *xlmi, double *xlrg, double *xlpk,
   double *yrtb, double *yrix, double *yrmi, double *yrrg, double *yrpk,
   double *yltb, double *ylix, double *ylmi, double *ylrg, double *ylpk,
   double *zrtb, double *zrix, double *zrmi, double *zrrg, double *zrpk,
   double *zltb, double *zlix, double *zlmi, double *zlrg, double *zlpk,
   int dim, int k )
{
   double prdrtb,prdrix,prdrmi,prdrrg,prdrpk;
   double prdltb,prdlix,prdlmi,prdlrg,prdlpk;
   int idx;

   // zv[k] = xv[0]*yv[k];
   dag_mul(xrtb[0],xrix[0],xrmi[0],xrrg[0],xrpk[0],
           xltb[0],xlix[0],xlmi[0],xlrg[0],xlpk[0],
           yrtb[k],yrix[k],yrmi[k],yrrg[k],yrpk[k],
           yltb[k],ylix[k],ylmi[k],ylrg[k],ylpk[k],
           &zrtb[k],&zrix[k],&zrmi[k],&zrrg[k],&zrpk[k],
           &zltb[k],&zlix[k],&zlmi[k],&zlrg[k],&zlpk[k]);

   for(int i=1; i<=k; i++) // zv[k] = zv[k] + xv[i]*yv[k-i];
   {
      idx = k-i;
      dag_mul(xrtb[i],xrix[i],xrmi[i],xrrg[i],xrpk[i],
              xltb[i],xlix[i],xlmi[i],xlrg[i],xlpk[i],
              yrtb[idx],yrix[idx],yrmi[idx],yrrg[idx],yrpk[idx],
              yltb[idx],ylix[idx],ylmi[idx],ylrg[idx],ylpk[idx],
              &prdrtb,&prdrix,&prdrmi,&prdrrg,&prdrpk,
              &prdltb,&prdlix,&prdlmi,&prdlrg,&prdlpk);
      dag_inc(&zrtb[k],&zrix[k],&zrmi[k],&zrrg[k],&zrpk[k],
              &zltb[k],&zlix[k],&zlmi[k],&zlrg[k],&zlpk[k],
              prdrtb,prdrix,prdrmi,prdrrg,prdrpk,
              prdltb,prdlix,prdlmi,prdlrg,prdlpk);
   }
}

__device__ void cmplx10_convolute
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
   double *zimlrg, double *zimlpk, int dim, int k )
{
}

__global__ void GPU_dbl10_speel
 ( int nvr, int deg, int *idx,
   double *cffrtb, double *cffrix, double *cffrmi, double *cffrrg,
   double *cffrpk, double *cffltb, double *cfflix, double *cfflmi,
   double *cfflrg, double *cfflpk,
   double *inputrtb, double *inputrix, double *inputrmi, double *inputrrg,
   double *inputrpk, double *inputltb, double *inputlix, double *inputlmi,
   double *inputlrg, double *inputlpk,
   double *forwardrtb, double *forwardrix, double *forwardrmi,
   double *forwardrrg, double *forwardrpk,
   double *forwardltb, double *forwardlix, double *forwardlmi,
   double *forwardlrg, double *forwardlpk,
   double *backwardrtb, double *backwardrix, double *backwardrmi,
   double *backwardrrg, double *backwardrpk,
   double *backwardltb, double *backwardlix, double *backwardlmi,
   double *backwardlrg, double *backwardlpk,
   double *crossrtb, double *crossrix, double *crossrmi, double *crossrrg,
   double *crossrpk, double *crossltb, double *crosslix, double *crosslmi,
   double *crosslrg, double *crosslpk )
{
   const int k = threadIdx.x;
   const int deg1 = deg+1;
   int ix1,ix2;

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
  
   xvrtb[k] = cffrtb[k]; xvrix[k] = cffrix[k]; xvrmi[k] = cffrmi[k]; 
   xvrrg[k] = cffrrg[k]; xvrpk[k] = cffrpk[k];
   xvltb[k] = cffltb[k]; xvlix[k] = cfflix[k]; xvlmi[k] = cfflmi[k]; 
   xvlrg[k] = cfflrg[k]; xvlpk[k] = cfflpk[k];
   ix1 = idx[0]*deg1+k;
   yvrtb[k] = inputrtb[ix1]; yvrix[k] = inputrix[ix1];
   yvrmi[k] = inputrmi[ix1]; yvrrg[k] = inputrrg[ix1];
   yvrpk[k] = inputrpk[ix1]; 
   yvltb[k] = inputltb[ix1]; yvlix[k] = inputlix[ix1];
   yvlmi[k] = inputlmi[ix1]; yvlrg[k] = inputlrg[ix1];
   yvlpk[k] = inputlpk[ix1]; 
   __syncthreads();
   dbl10_convolute(xvrtb,xvrix,xvrmi,xvrrg,xvrpk,
                   xvltb,xvlix,xvlmi,xvlrg,xvlpk,
                   yvrtb,yvrix,yvrmi,yvrrg,yvrpk,
                   yvltb,yvlix,yvlmi,yvlrg,yvlpk,
                   zvrtb,zvrix,zvrmi,zvrrg,zvrpk,
                   zvltb,zvlix,zvlmi,zvlrg,zvlpk,deg1,k);
   __syncthreads();
   forwardrtb[k] = zvrtb[k]; forwardrix[k] = zvrix[k];
   forwardrmi[k] = zvrmi[k]; forwardrrg[k] = zvrrg[k];
   forwardrpk[k] = zvrpk[k];
   forwardltb[k] = zvltb[k]; forwardlix[k] = zvlix[k];
   forwardlmi[k] = zvlmi[k]; forwardlrg[k] = zvlrg[k];
   forwardlpk[k] = zvlpk[k];                            // f[0] = cff*x[0]

   for(int i=1; i<nvr; i++)
   {
      xvrtb[k] = zvrtb[k]; xvrix[k] = zvrix[k]; xvrmi[k] = zvrmi[k];
      xvrrg[k] = zvrrg[k]; xvrpk[k] = zvrpk[k];
      xvltb[k] = zvltb[k]; xvlix[k] = zvlix[k]; xvlmi[k] = zvlmi[k];
      xvlrg[k] = zvlrg[k]; xvlpk[k] = zvlpk[k];
      ix2 = idx[i]*deg1+k;
      yvrtb[k] = inputrtb[ix2]; yvrix[k] = inputrix[ix2];
      yvrmi[k] = inputrmi[ix2]; yvrrg[k] = inputrrg[ix2];
      yvrpk[k] = inputrpk[ix2];
      yvltb[k] = inputltb[ix2]; yvlix[k] = inputlix[ix2];
      yvlmi[k] = inputlmi[ix2]; yvlrg[k] = inputlrg[ix2];
      yvlpk[k] = inputlpk[ix2];
      __syncthreads();
      dbl10_convolute(xvrtb,xvrix,xvrmi,xvrrg,xvrpk,
                      xvltb,xvlix,xvlmi,xvlrg,xvlpk,
                      yvrtb,yvrix,yvrmi,yvrrg,yvrpk,
                      yvltb,yvlix,yvlmi,yvlrg,yvlpk,
                      zvrtb,zvrix,zvrmi,zvrrg,zvrpk,
                      zvltb,zvlix,zvlmi,zvlrg,zvlpk,deg1,k);
      __syncthreads();
      ix1 = i*deg1+k;
      forwardrtb[ix1] = zvrtb[k]; forwardrix[ix1] = zvrix[k]; 
      forwardrmi[ix1] = zvrmi[k]; forwardrrg[ix1] = zvrrg[k]; 
      forwardrpk[ix1] = zvrpk[k];
      forwardltb[ix1] = zvltb[k]; forwardlix[ix1] = zvlix[k]; 
      forwardlmi[ix1] = zvlmi[k]; forwardlrg[ix1] = zvlrg[k]; 
      forwardlpk[ix1] = zvlpk[k];                     // f[i] = f[i-1]*x[i]
   }
   if(nvr > 2)
   {
      ix1 = idx[nvr-1]*deg1+k;
      xvrtb[k] = inputrtb[ix1]; xvrix[k] = inputrix[ix1];
      xvrmi[k] = inputrmi[ix1]; xvrrg[k] = inputrrg[ix1];
      xvrpk[k] = inputrpk[ix1];
      xvltb[k] = inputltb[ix1]; xvlix[k] = inputlix[ix1];
      xvlmi[k] = inputlmi[ix1]; xvlrg[k] = inputlrg[ix1];
      xvlpk[k] = inputlpk[ix1];
      ix2 = idx[nvr-2]*deg1+k;
      yvrtb[k] = inputrtb[ix2]; yvrix[k] = inputrix[ix2];
      yvrmi[k] = inputrmi[ix2]; yvrrg[k] = inputrrg[ix2];
      yvrpk[k] = inputrpk[ix2];
      yvltb[k] = inputltb[ix2]; yvlix[k] = inputlix[ix2];
      yvlmi[k] = inputlmi[ix2]; yvlrg[k] = inputlrg[ix2];
      yvlpk[k] = inputlpk[ix2];
      __syncthreads();
      dbl10_convolute(xvrtb,xvrix,xvrmi,xvrrg,xvrpk,
                      xvltb,xvlix,xvlmi,xvlrg,xvlpk,
                      yvrtb,yvrix,yvrmi,yvrrg,yvrpk,
                      yvltb,yvlix,yvlmi,yvlrg,yvlpk,
                      zvrtb,zvrix,zvrmi,zvrrg,zvrpk,
                      zvltb,zvlix,zvlmi,zvlrg,zvlpk,deg1,k);
      __syncthreads();
      backwardrtb[k] = zvrtb[k]; backwardrix[k] = zvrix[k];
      backwardrmi[k] = zvrmi[k]; backwardrrg[k] = zvrrg[k];
      backwardrpk[k] = zvrpk[k]; 
      backwardltb[k] = zvltb[k]; backwardlix[k] = zvlix[k];
      backwardlmi[k] = zvlmi[k]; backwardlrg[k] = zvlrg[k];
      backwardlpk[k] = zvlpk[k];                     // b[0] = x[n-1]*x[n-2]

      for(int i=1; i<nvr-2; i++)
      {
         xvrtb[k] = zvrtb[k]; xvrix[k] = zvrix[k]; xvrmi[k] = zvrmi[k]; 
         xvrrg[k] = zvrrg[k]; xvrpk[k] = zvrpk[k];
         xvltb[k] = zvltb[k]; xvlix[k] = zvlix[k]; xvlmi[k] = zvlmi[k]; 
         xvlrg[k] = zvlrg[k]; xvlpk[k] = zvlpk[k];
         ix2 = idx[nvr-2-i]*deg1+k;
         yvrtb[k] = inputrtb[ix2]; yvrix[k] = inputrix[ix2];
         yvrmi[k] = inputrmi[ix2]; yvrrg[k] = inputrrg[ix2];
         yvrpk[k] = inputrpk[ix2];
         yvltb[k] = inputltb[ix2]; yvlix[k] = inputlix[ix2];
         yvlmi[k] = inputlmi[ix2]; yvlrg[k] = inputlrg[ix2];
         yvlpk[k] = inputlpk[ix2];
         __syncthreads();
         dbl10_convolute(xvrtb,xvrix,xvrmi,xvrrg,xvrpk,
                         xvltb,xvlix,xvlmi,xvlrg,xvlpk,
                         yvrtb,yvrix,yvrmi,yvrrg,yvrpk,
                         yvltb,yvlix,yvlmi,yvlrg,yvlpk,
                         zvrtb,zvrix,zvrmi,zvrrg,zvrpk,
                         zvltb,zvlix,zvlmi,zvlrg,zvlpk,deg1,k);
         __syncthreads();
         ix1 = i*deg1+k;
         backwardrtb[ix1] = zvrtb[k]; backwardrix[ix1] = zvrix[k];
         backwardrmi[ix1] = zvrmi[k]; backwardrrg[ix1] = zvrrg[k];
         backwardrpk[ix1] = zvrpk[k]; 
         backwardltb[ix1] = zvltb[k]; backwardlix[ix1] = zvlix[k];
         backwardlmi[ix1] = zvlmi[k]; backwardlrg[ix1] = zvlrg[k];
         backwardlpk[ix1] = zvlpk[k];            // b[i] = b[i-1]*x[n-2-i]
      }
      xvrtb[k] = zvrtb[k];  xvrix[k] = zvrix[k];  xvrmi[k] = zvrmi[k];
      xvrrg[k] = zvrrg[k];  xvrpk[k] = zvrpk[k];
      xvltb[k] = zvltb[k];  xvlix[k] = zvlix[k];  xvlmi[k] = zvlmi[k];
      xvlrg[k] = zvlrg[k];  xvlpk[k] = zvlpk[k];
      yvrtb[k] = cffrtb[k]; yvrix[k] = cffrix[k]; yvrmi[k] = cffrmi[k];
      yvrrg[k] = cffrrg[k]; yvrpk[k] = cffrpk[k];
      yvltb[k] = cffltb[k]; yvlix[k] = cfflix[k]; yvlmi[k] = cfflmi[k];
      yvlrg[k] = cfflrg[k]; yvlpk[k] = cfflpk[k];
      __syncthreads();
      dbl10_convolute(xvrtb,xvrix,xvrmi,xvrrg,xvrpk,
                      xvltb,xvlix,xvlmi,xvlrg,xvlpk,
                      yvrtb,yvrix,yvrmi,yvrrg,yvrpk,
                      yvltb,yvlix,yvlmi,yvlrg,yvlpk,
                      zvrtb,zvrix,zvrmi,zvrrg,zvrpk,
                      zvltb,zvlix,zvlmi,zvlrg,zvlpk,deg1,k);
      __syncthreads();
      ix2 = (nvr-3)*deg1+k;
      backwardrtb[ix2] = zvrtb[k]; backwardrix[ix2] = zvrix[k];
      backwardrmi[ix2] = zvrmi[k]; backwardrrg[ix2] = zvrrg[k];
      backwardrpk[ix2] = zvrpk[k];
      backwardltb[ix2] = zvltb[k]; backwardlix[ix2] = zvlix[k];
      backwardlmi[ix2] = zvlmi[k]; backwardlrg[ix2] = zvlrg[k];
      backwardlpk[ix2] = zvlpk[k];                   // b[n-3] = b[n-3]*cff

      if(nvr == 3)
      {
         xvrtb[k] = forwardrtb[k]; xvrix[k] = forwardrix[k];
         xvrmi[k] = forwardrmi[k]; xvrrg[k] = forwardrrg[k];
         xvrpk[k] = forwardrpk[k];
         xvltb[k] = forwardltb[k]; xvlix[k] = forwardlix[k];
         xvlmi[k] = forwardlmi[k]; xvlrg[k] = forwardlrg[k];
         xvlpk[k] = forwardlpk[k];
         ix2 = idx[2]*deg1+k;
         yvrtb[k] = inputrtb[ix2]; yvrix[k] = inputrix[ix2];
         yvrmi[k] = inputrmi[ix2]; yvrrg[k] = inputrrg[ix2]; 
         yvrpk[k] = inputrpk[ix2];
         yvltb[k] = inputltb[ix2]; yvlix[k] = inputlix[ix2];
         yvlmi[k] = inputlmi[ix2]; yvlrg[k] = inputlrg[ix2]; 
         yvlpk[k] = inputlpk[ix2];
         __syncthreads();
         dbl10_convolute(xvrtb,xvrix,xvrmi,xvrrg,xvrpk,
                         xvltb,xvlix,xvlmi,xvlrg,xvlpk,
                         yvrtb,yvrix,yvrmi,yvrrg,yvrpk,
                         yvltb,yvlix,yvlmi,yvlrg,yvlpk,
                         zvrtb,zvrix,zvrmi,zvrrg,zvrpk,
                         zvltb,zvlix,zvlmi,zvlrg,zvlpk,deg1,k);
         __syncthreads();
         crossrtb[k] = zvrtb[k]; crossrix[k] = zvrix[k];
         crossrmi[k] = zvrmi[k]; crossrrg[k] = zvrrg[k];
         crossrpk[k] = zvrpk[k]; 
         crossltb[k] = zvltb[k]; crosslix[k] = zvlix[k];
         crosslmi[k] = zvlmi[k]; crosslrg[k] = zvlrg[k];
         crosslpk[k] = zvlpk[k];                       // c[0] = f[0]*x[2]
      }
      else
      {
         for(int i=0; i<nvr-3; i++)
         {
            ix1 = i*deg1+k; 
            xvrtb[k] = forwardrtb[ix1]; xvrix[k] = forwardrix[ix1];
            xvrmi[k] = forwardrmi[ix1]; xvrrg[k] = forwardrrg[ix1];
            xvrpk[k] = forwardrpk[ix1];
            xvltb[k] = forwardltb[ix1]; xvlix[k] = forwardlix[ix1];
            xvlmi[k] = forwardlmi[ix1]; xvlrg[k] = forwardlrg[ix1];
            xvlpk[k] = forwardlpk[ix1];
            ix2 = (nvr-4-i)*deg1+k;
            yvrtb[k] = backwardrtb[ix2]; yvrix[k] = backwardrix[ix2];
            yvrmi[k] = backwardrmi[ix2]; yvrrg[k] = backwardrrg[ix2];
            yvrpk[k] = backwardrpk[ix2];
            yvltb[k] = backwardltb[ix2]; yvlix[k] = backwardlix[ix2];
            yvlmi[k] = backwardlmi[ix2]; yvlrg[k] = backwardlrg[ix2];
            yvlpk[k] = backwardlpk[ix2];
            __syncthreads();
            dbl10_convolute(xvrtb,xvrix,xvrmi,xvrrg,xvrpk,
                            xvltb,xvlix,xvlmi,xvlrg,xvlpk,
                            yvrtb,yvrix,yvrmi,yvrrg,yvrpk,
                            yvltb,yvlix,yvlmi,yvlrg,yvlpk,
                            zvrtb,zvrix,zvrmi,zvrrg,zvrpk,
                            zvltb,zvlix,zvlmi,zvlrg,zvlpk,deg1,k);
            __syncthreads();
            crossrtb[ix1] = zvrtb[k]; crossrix[ix1] = zvrix[k];
            crossrmi[ix1] = zvrmi[k]; crossrrg[ix1] = zvrrg[k];
            crossrpk[ix1] = zvrpk[k];
            crossltb[ix1] = zvltb[k]; crosslix[ix1] = zvlix[k];
            crosslmi[ix1] = zvlmi[k]; crosslrg[ix1] = zvlrg[k];
            crosslpk[ix1] = zvlpk[k];            // c[i] = f[i]*b[n-4-i]
         }
         ix1 = (nvr-3)*deg1+k;
         xvrtb[k] = forwardrtb[ix1]; xvrix[k] = forwardrix[ix1];
         xvrmi[k] = forwardrmi[ix1]; xvrrg[k] = forwardrrg[ix1];
         xvrpk[k] = forwardrpk[ix1];
         xvltb[k] = forwardltb[ix1]; xvlix[k] = forwardlix[ix1];
         xvlmi[k] = forwardlmi[ix1]; xvlrg[k] = forwardlrg[ix1];
         xvlpk[k] = forwardlpk[ix1];
         ix2 = idx[nvr-1]*deg1+k;
         yvrtb[k] = inputrtb[ix2]; yvrix[k] = inputrix[ix2];
         yvrmi[k] = inputrmi[ix2]; yvrrg[k] = inputrrg[ix2];
         yvrpk[k] = inputrpk[ix2];
         yvltb[k] = inputltb[ix2]; yvlix[k] = inputlix[ix2];
         yvlmi[k] = inputlmi[ix2]; yvlrg[k] = inputlrg[ix2];
         yvlpk[k] = inputlpk[ix2];
         __syncthreads();
         dbl10_convolute(xvrtb,xvrix,xvrmi,xvrrg,xvrpk,
                         xvltb,xvlix,xvlmi,xvlrg,xvlpk,
                         yvrtb,yvrix,yvrmi,yvrrg,yvrpk,
                         yvltb,yvlix,yvlmi,yvlrg,yvlpk,
                         zvrtb,zvrix,zvrmi,zvrrg,zvrpk,
                         zvltb,zvlix,zvlmi,zvlrg,zvlpk,deg1,k);
         __syncthreads();
         crossrtb[ix1] = zvrtb[k]; crossrix[ix1] = zvrix[k];
         crossrmi[ix1] = zvrmi[k]; crossrrg[ix1] = zvrrg[k];
         crossrpk[ix1] = zvrpk[k];
         crossltb[ix1] = zvltb[k]; crosslix[ix1] = zvlix[k];
         crosslmi[ix1] = zvlmi[k]; crosslrg[ix1] = zvlrg[k];
         crosslpk[ix1] = zvlpk[k];                 // c[n-3] = f[n-3]*x[n-1]
      }
   }
}

__global__ void GPU_cmplx10_speel
 ( int nvr, int deg, int *idx,
   double *cffrertb, double *cffrerix, double *cffrermi, double *cffrerrg,
   double *cffrerpk, double *cffreltb, double *cffrelix, double *cffrelmi,
   double *cffrelrg, double *cffrelpk,
   double *cffimrtb, double *cffimrix, double *cffimrmi, double *cffimrrg,
   double *cffimrpk, double *cffimltb, double *cffimlix, double *cffimlmi,
   double *cffimlrg, double *cffimlpk,
   double *inputrertb, double *inputrerix, double *inputrermi,
   double *inputrerrg, double *inputrerpk,
   double *inputreltb, double *inputrelix, double *inputrelmi,
   double *inputrelrg, double *inputrelpk,
   double *inputimrtb, double *inputimrix, double *inputimrmi,
   double *inputimrrg, double *inputimrpk,
   double *inputimltb, double *inputimlix, double *inputimlmi,
   double *inputimlrg, double *inputimlpk,
   double *forwardrertb, double *forwardrerix, double *forwardrermi,
   double *forwardrerrg, double *forwardrerpk,
   double *forwardreltb, double *forwardrelix, double *forwardrelmi,
   double *forwardrelrg, double *forwardrelpk,
   double *forwardimrtb, double *forwardimrix, double *forwardimrmi,
   double *forwardimrrg, double *forwardimrpk,
   double *forwardimltb, double *forwardimlix, double *forwardimlmi,
   double *forwardimlrg, double *forwardimlpk,
   double *backwardrertb, double *backwardrerix, double *backwardrermi,
   double *backwardrerrg, double *backwardrerpk,
   double *backwardreltb, double *backwardrelix, double *backwardrelmi,
   double *backwardrelrg, double *backwardrelpk,
   double *backwardimrtb, double *backwardimrix, double *backwardimrmi,
   double *backwardimrrg, double *backwardimrpk,
   double *backwardimltb, double *backwardimlix, double *backwardimlmi,
   double *backwardimlrg, double *backwardimlpk,
   double *crossrertb, double *crossrerix, double *crossrermi,
   double *crossrerrg, double *crossrerpk,
   double *crossreltb, double *crossrelix, double *crossrelmi,
   double *crossrelrg, double *crossrelpk,
   double *crossimrtb, double *crossimrix, double *crossimrmi,
   double *crossimrrg, double *crossimrpk,
   double *crossimltb, double *crossimlix, double *crossimlmi,
   double *crossimlrg, double *crossimlpk )
{
}

void GPU_dbl10_evaldiff
 ( int BS, int dim, int nvr, int deg, int *idx,
   double *cffrtb, double *cffrix, double *cffrmi, double *cffrrg,
   double *cffrpk, double *cffltb, double *cfflix, double *cfflmi,
   double *cfflrg, double *cfflpk,
   double **inputrtb, double **inputrix, double **inputrmi,
   double **inputrrg, double **inputrpk,
   double **inputltb, double **inputlix, double **inputlmi,
   double **inputlrg, double **inputlpk,
   double **outputrtb, double **outputrix, double **outputrmi,
   double **outputrrg, double **outputrpk,
   double **outputltb, double **outputlix, double **outputlmi,
   double **outputlrg, double **outputlpk )
{
   const int deg1 = deg+1;      // length of all vectors
   double *inputrtb_d;          // inputrtb_d is inputrtb on the device
   double *inputrix_d;          // inputrix_d is inputrix on the device
   double *inputrmi_d;          // inputrmi_d is inputrmi on the device
   double *inputrrg_d;          // inputrrg_d is inputrrg on the device
   double *inputrpk_d;          // inputrpk_d is inputrpk on the device
   double *inputltb_d;          // inputltb_d is inputltb on the device
   double *inputlix_d;          // inputlix_d is inputlix on the device
   double *inputlmi_d;          // inputlmi_d is inputlmi on the device
   double *inputlrg_d;          // inputlrg_d is inputlrg on the device
   double *inputlpk_d;          // inputlpk_d is inputlpk on the device
   double *forwardrtb_d;        // highest forward products on the device
   double *forwardrix_d;        // second highest forward products
   double *forwardrmi_d;        // third highest forward products
   double *forwardrrg_d;        // fourth highest forward products
   double *forwardrpk_d;        // fifth highest forward products
   double *forwardltb_d;        // fifth lowest forward products
   double *forwardlix_d;        // fourth lowest forward products
   double *forwardlmi_d;        // third lowest forward products
   double *forwardlrg_d;        // second lowest forward products
   double *forwardlpk_d;        // lowest forward products
   double *backwardrtb_d;       // highest backward products on the device
   double *backwardrix_d;       // second highest backward products
   double *backwardrmi_d;       // third highest backward products
   double *backwardrrg_d;       // fourth highest backward products
   double *backwardrpk_d;       // fifth highest backward products
   double *backwardltb_d;       // fifth lowest backward products
   double *backwardlix_d;       // fourth lowest backward products
   double *backwardlmi_d;       // third lowest backward products
   double *backwardlrg_d;       // second lowest backward products
   double *backwardlpk_d;       // lowest backward products
   double *crossrtb_d;          // highest cross products on the device
   double *crossrix_d;          // second highest cross products
   double *crossrmi_d;          // third highest cross products
   double *crossrrg_d;          // fourth highest cross products
   double *crossrpk_d;          // fifth highest cross products
   double *crossltb_d;          // fifth lowest cross products
   double *crosslix_d;          // fourth lowest cross products
   double *crosslmi_d;          // third lowest cross products
   double *crosslrg_d;          // second lowest cross products
   double *crosslpk_d;          // lowest cross products
   double *cffrtb_d;            // cffrtb_d is cffrtb on device
   double *cffrix_d;            // cffrix_d is cffrix on device
   double *cffrmi_d;            // cffrmi_d is cffrmi on device
   double *cffrrg_d;            // cffrrg_d is cffrrg on device
   double *cffrpk_d;            // cffrpk_d is cffrpk on device
   double *cffltb_d;            // cffltb_d is cffltb on device
   double *cfflix_d;            // cfflix_d is cfflix on device
   double *cfflmi_d;            // cfflmi_d is cfflmi on device
   double *cfflrg_d;            // cfflrg_d is cfflrg on device
   double *cfflpk_d;            // cfflpk_d is cfflpk on device
   int *idx_d;                  // idx_d is idx on device

   size_t szcff = deg1*sizeof(double);
   size_t szdim = dim*(deg1)*sizeof(double);
   size_t sznvr = nvr*(deg1)*sizeof(double);
   size_t sznvr2 = (nvr-2)*(deg1)*sizeof(double);
   size_t szidx = nvr*sizeof(int);

   cudaMalloc((void**)&idx_d,szidx);
   cudaMalloc((void**)&cffrtb_d,szcff);
   cudaMalloc((void**)&cffrix_d,szcff);
   cudaMalloc((void**)&cffrmi_d,szcff);
   cudaMalloc((void**)&cffrrg_d,szcff);
   cudaMalloc((void**)&cffrpk_d,szcff);
   cudaMalloc((void**)&cffltb_d,szcff);
   cudaMalloc((void**)&cfflix_d,szcff);
   cudaMalloc((void**)&cfflmi_d,szcff);
   cudaMalloc((void**)&cfflrg_d,szcff);
   cudaMalloc((void**)&cfflpk_d,szcff);
   cudaMalloc((void**)&inputrtb_d,szdim);
   cudaMalloc((void**)&inputrix_d,szdim);
   cudaMalloc((void**)&inputrmi_d,szdim);
   cudaMalloc((void**)&inputrrg_d,szdim);
   cudaMalloc((void**)&inputrpk_d,szdim);
   cudaMalloc((void**)&inputltb_d,szdim);
   cudaMalloc((void**)&inputlix_d,szdim);
   cudaMalloc((void**)&inputlmi_d,szdim);
   cudaMalloc((void**)&inputlrg_d,szdim);
   cudaMalloc((void**)&inputlpk_d,szdim);
   cudaMalloc((void**)&forwardrtb_d,sznvr);
   cudaMalloc((void**)&forwardrix_d,sznvr);
   cudaMalloc((void**)&forwardrmi_d,sznvr);
   cudaMalloc((void**)&forwardrrg_d,sznvr);
   cudaMalloc((void**)&forwardrpk_d,sznvr);
   cudaMalloc((void**)&forwardltb_d,sznvr);
   cudaMalloc((void**)&forwardlix_d,sznvr);
   cudaMalloc((void**)&forwardlmi_d,sznvr);
   cudaMalloc((void**)&forwardlrg_d,sznvr);
   cudaMalloc((void**)&forwardlpk_d,sznvr);
   cudaMalloc((void**)&backwardrtb_d,sznvr2);
   cudaMalloc((void**)&backwardrix_d,sznvr2);
   cudaMalloc((void**)&backwardrmi_d,sznvr2);
   cudaMalloc((void**)&backwardrrg_d,sznvr2);
   cudaMalloc((void**)&backwardrpk_d,sznvr2);
   cudaMalloc((void**)&backwardltb_d,sznvr2);
   cudaMalloc((void**)&backwardlix_d,sznvr2);
   cudaMalloc((void**)&backwardlmi_d,sznvr2);
   cudaMalloc((void**)&backwardlrg_d,sznvr2);
   cudaMalloc((void**)&backwardlpk_d,sznvr2);
   cudaMalloc((void**)&crossrtb_d,sznvr2);
   cudaMalloc((void**)&crossrix_d,sznvr2);
   cudaMalloc((void**)&crossrmi_d,sznvr2);
   cudaMalloc((void**)&crossrrg_d,sznvr2);
   cudaMalloc((void**)&crossrpk_d,sznvr2);
   cudaMalloc((void**)&crossltb_d,sznvr2);
   cudaMalloc((void**)&crosslix_d,sznvr2);
   cudaMalloc((void**)&crosslmi_d,sznvr2);
   cudaMalloc((void**)&crosslrg_d,sznvr2);
   cudaMalloc((void**)&crosslpk_d,sznvr2);

   double *inputrtb_h = new double[dim*(deg1)];
   double *inputrix_h = new double[dim*(deg1)];
   double *inputrmi_h = new double[dim*(deg1)];
   double *inputrrg_h = new double[dim*(deg1)];
   double *inputrpk_h = new double[dim*(deg1)];
   double *inputltb_h = new double[dim*(deg1)];
   double *inputlix_h = new double[dim*(deg1)];
   double *inputlmi_h = new double[dim*(deg1)];
   double *inputlrg_h = new double[dim*(deg1)];
   double *inputlpk_h = new double[dim*(deg1)];
   int ix = 0;
   for(int i=0; i<dim; i++)
      for(int j=0; j<deg1; j++)
      {
         inputrtb_h[ix] = inputrtb[i][j];
         inputrix_h[ix] = inputrix[i][j];
         inputrmi_h[ix] = inputrmi[i][j];
         inputrrg_h[ix] = inputrrg[i][j];
         inputrpk_h[ix] = inputrpk[i][j];
         inputltb_h[ix] = inputltb[i][j];
         inputlix_h[ix] = inputlix[i][j];
         inputlmi_h[ix] = inputlmi[i][j];
         inputlrg_h[ix] = inputlrg[i][j];
         inputlpk_h[ix++] = inputlpk[i][j];
      }

   cudaMemcpy(idx_d,idx,szidx,cudaMemcpyHostToDevice);
   cudaMemcpy(cffrtb_d,cffrtb,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffrix_d,cffrix,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffrmi_d,cffrmi,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffrrg_d,cffrrg,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffrpk_d,cffrpk,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffltb_d,cffltb,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cfflix_d,cfflix,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cfflmi_d,cfflmi,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cfflrg_d,cfflrg,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cfflpk_d,cfflpk,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(inputrtb_d,inputrtb_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputrix_d,inputrix_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputrmi_d,inputrmi_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputrrg_d,inputrrg_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputrpk_d,inputrpk_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputltb_d,inputltb_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputlix_d,inputlix_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputlmi_d,inputlmi_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputlrg_d,inputlrg_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputlpk_d,inputlpk_h,szdim,cudaMemcpyHostToDevice);

   if(BS == deg1)
   {
      GPU_dbl10_speel<<<1,BS>>>
         (nvr,deg,idx_d,
          cffrtb_d,cffrix_d,cffrmi_d,cffrrg_d,cffrpk_d,
          cffltb_d,cfflix_d,cfflmi_d,cfflrg_d,cfflpk_d,
          inputrtb_d,inputrix_d,inputrmi_d,inputrrg_d,inputrpk_d,
          inputltb_d,inputlix_d,inputlmi_d,inputlrg_d,inputlpk_d,
          forwardrtb_d,forwardrix_d,forwardrmi_d,forwardrrg_d,forwardrpk_d,
          forwardltb_d,forwardlix_d,forwardlmi_d,forwardlrg_d,forwardlpk_d,
          backwardrtb_d,backwardrix_d,backwardrmi_d,backwardrrg_d,
          backwardlpk_d,backwardltb_d,backwardlix_d,backwardlmi_d,
          backwardlrg_d,backwardlpk_d,
          crossrtb_d,crossrix_d,crossrmi_d,crossrrg_d,crossrpk_d,
          crossltb_d,crosslix_d,crosslmi_d,crosslrg_d,crosslpk_d);
   }
   double *forwardrtb_h = new double[(deg1)*nvr];
   double *forwardrix_h = new double[(deg1)*nvr];
   double *forwardrmi_h = new double[(deg1)*nvr];
   double *forwardrrg_h = new double[(deg1)*nvr];
   double *forwardrpk_h = new double[(deg1)*nvr];
   double *forwardltb_h = new double[(deg1)*nvr];
   double *forwardlix_h = new double[(deg1)*nvr];
   double *forwardlmi_h = new double[(deg1)*nvr];
   double *forwardlrg_h = new double[(deg1)*nvr];
   double *forwardlpk_h = new double[(deg1)*nvr];
   double *backwardrtb_h = new double[(deg1)*(nvr-2)];
   double *backwardrix_h = new double[(deg1)*(nvr-2)];
   double *backwardrmi_h = new double[(deg1)*(nvr-2)];
   double *backwardrrg_h = new double[(deg1)*(nvr-2)];
   double *backwardrpk_h = new double[(deg1)*(nvr-2)];
   double *backwardltb_h = new double[(deg1)*(nvr-2)];
   double *backwardlix_h = new double[(deg1)*(nvr-2)];
   double *backwardlmi_h = new double[(deg1)*(nvr-2)];
   double *backwardlrg_h = new double[(deg1)*(nvr-2)];
   double *backwardlpk_h = new double[(deg1)*(nvr-2)];
   double *crossrtb_h = new double[(deg1)*(nvr-2)];
   double *crossrix_h = new double[(deg1)*(nvr-2)];
   double *crossrmi_h = new double[(deg1)*(nvr-2)];
   double *crossrrg_h = new double[(deg1)*(nvr-2)];
   double *crossrpk_h = new double[(deg1)*(nvr-2)];
   double *crossltb_h = new double[(deg1)*(nvr-2)];
   double *crosslix_h = new double[(deg1)*(nvr-2)];
   double *crosslmi_h = new double[(deg1)*(nvr-2)];
   double *crosslrg_h = new double[(deg1)*(nvr-2)];
   double *crosslpk_h = new double[(deg1)*(nvr-2)];
  
   cudaMemcpy(forwardrtb_h,forwardrtb_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardrix_h,forwardrix_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardrmi_h,forwardrmi_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardrrg_h,forwardrrg_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardrpk_h,forwardrpk_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardltb_h,forwardltb_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardlix_h,forwardlix_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardlmi_h,forwardlmi_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardlrg_h,forwardlrg_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardlpk_h,forwardlpk_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardrtb_h,backwardrtb_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardrix_h,backwardrix_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardrmi_h,backwardrmi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardrrg_h,backwardrrg_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardrpk_h,backwardrpk_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardltb_h,backwardltb_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardlix_h,backwardlix_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardlmi_h,backwardlmi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardlrg_h,backwardlrg_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardlpk_h,backwardlpk_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossrtb_h,crossrtb_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossrix_h,crossrix_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossrmi_h,crossrmi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossrrg_h,crossrrg_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossrpk_h,crossrpk_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossltb_h,crossltb_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crosslix_h,crosslix_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crosslmi_h,crosslmi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crosslrg_h,crosslrg_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crosslpk_h,crosslpk_d,sznvr2,cudaMemcpyDeviceToHost);

   int offset = (nvr-1)*deg1;            // assign value of the monomial
   for(int i=0; i<deg1; i++)
   {
      outputrtb[dim][i] = forwardrtb_h[offset+i];
      outputrix[dim][i] = forwardrix_h[offset+i];
      outputrmi[dim][i] = forwardrmi_h[offset+i];
      outputrrg[dim][i] = forwardrrg_h[offset+i];
      outputrpk[dim][i] = forwardrpk_h[offset+i];
      outputltb[dim][i] = forwardltb_h[offset+i];
      outputlix[dim][i] = forwardlix_h[offset+i];
      outputlmi[dim][i] = forwardlmi_h[offset+i];
      outputlrg[dim][i] = forwardlrg_h[offset+i];
      outputlpk[dim][i] = forwardlpk_h[offset+i];
   }
   ix = idx[nvr-1];                      // derivative with respect to x[n-1]
   offset = (nvr-2)*deg1;
   for(int i=0; i<deg1; i++)
   {
      outputrtb[ix][i] = forwardrtb_h[offset+i];
      outputrix[ix][i] = forwardrix_h[offset+i];
      outputrmi[ix][i] = forwardrmi_h[offset+i];
      outputrrg[ix][i] = forwardrrg_h[offset+i];
      outputrpk[ix][i] = forwardrpk_h[offset+i];
      outputltb[ix][i] = forwardltb_h[offset+i];
      outputlix[ix][i] = forwardlix_h[offset+i];
      outputlmi[ix][i] = forwardlmi_h[offset+i];
      outputlrg[ix][i] = forwardlrg_h[offset+i];
      outputlpk[ix][i] = forwardlpk_h[offset+i];
   }
   ix = idx[0];                          // derivative with respect to x[0]
   offset = (nvr-3)*deg1;
   for(int i=0; i<deg1; i++)
   {
      outputrtb[ix][i] = backwardrtb_h[offset+i];
      outputrix[ix][i] = backwardrix_h[offset+i];
      outputrmi[ix][i] = backwardrmi_h[offset+i];
      outputrrg[ix][i] = backwardrrg_h[offset+i];
      outputrpk[ix][i] = backwardrpk_h[offset+i];
      outputltb[ix][i] = backwardltb_h[offset+i];
      outputlix[ix][i] = backwardlix_h[offset+i];
      outputlmi[ix][i] = backwardlmi_h[offset+i];
      outputlrg[ix][i] = backwardlrg_h[offset+i];
      outputlpk[ix][i] = backwardlpk_h[offset+i];
   }
   for(int k=1; k<nvr-1; k++)            // derivative with respect to x[k]
   {
      ix = idx[k]; offset = (k-1)*deg1;
      for(int i=0; i<deg1; i++)
      {
         outputrtb[ix][i] = crossrtb_h[offset+i];
         outputrix[ix][i] = crossrix_h[offset+i];
         outputrmi[ix][i] = crossrmi_h[offset+i];
         outputrrg[ix][i] = crossrrg_h[offset+i];
         outputrpk[ix][i] = crossrpk_h[offset+i];
         outputltb[ix][i] = crossltb_h[offset+i];
         outputlix[ix][i] = crosslix_h[offset+i];
         outputlmi[ix][i] = crosslmi_h[offset+i];
         outputlrg[ix][i] = crosslrg_h[offset+i];
         outputlpk[ix][i] = crosslpk_h[offset+i];
      }
   }
}

void GPU_cmplx10_evaldiff
 ( int BS, int dim, int nvr, int deg, int *idx,
   double *cffretb, double *cffreix, double *cffremi,
   double *cffrerg, double *cffrepk,
   double *cffimtb, double *cffimix, double *cffimmi,
   double *cffimrg, double *cffimpk,
   double **inputretb, double **inputreix, double **inputremi,
   double **inputrerg, double **inputrepk,
   double **inputimtb, double **inputimix, double **inputimmi,
   double **inputimrg, double **inputimpk,
   double **outputretb, double **outputreix, double **outputremi,
   double **outputrerg, double **outputrepk,
   double **outputimtb, double **outputimix, double **outputimmi,
   double **outputimrg, double **outputimpk )
{
}
