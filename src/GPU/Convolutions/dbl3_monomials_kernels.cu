// The file dbl3_monomials_kernels.cu defines the kernels specified
// in dbl3_monomials_kernels.h.

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
#include "triple_double_gpufun.cu"
#include "dbl3_convolutions_kernels.h"
#include "dbl3_monomials_kernels.h"

__device__ void dbl3_convolute
 ( double *xhi, double *xmi, double *xlo,
   double *yhi, double *ymi, double *ylo,
   double *zhi, double *zmi, double *zlo, int dim, int k )
{
   double prdhi,prdmi,prdlo;

   // z[k] = x[0]*y[k];
   tdg_mul(xhi[0],xmi[0],xlo[0],yhi[k],ymi[k],ylo[k],&zhi[k],&zmi[k],&zlo[k]);

   for(int i=1; i<=k; i++) // z[k] = z[k] + x[i]*y[k-i];
   {
      tdg_mul(xhi[i],xmi[i],xlo[i],yhi[k-i],ymi[k-i],ylo[k-i],
              &prdhi,&prdmi,&prdlo);
      tdg_inc(&zhi[k],&zmi[k],&zlo[k],prdhi,prdmi,prdlo);
   }
}

__device__ void cmplx3_convolute
 ( double *xrehi, double *xremi, double *xrelo,
   double *ximhi, double *ximmi, double *ximlo,
   double *yrehi, double *yremi, double *yrelo,
   double *yimhi, double *yimmi, double *yimlo,
   double *zrehi, double *zremi, double *zrelo,
   double *zimhi, double *zimmi, double *zimlo, int dim, int k )
{
   double xrhi,xihi,yrhi,yihi,zrhi,zihi,acchi;
   double xrmi,ximi,yrmi,yimi,zrmi,zimi,accmi;
   double xrlo,xilo,yrlo,yilo,zrlo,zilo,acclo;

   // z[k] = x[0]*y[k]
   xrhi = xrehi[0]; xrmi = xremi[0]; xrlo = xrelo[0];
   xihi = ximhi[0]; ximi = ximmi[0]; xilo = ximlo[0];
   yrhi = yrehi[k]; yrmi = yremi[k]; yrlo = yrelo[k];
   yihi = yimhi[k]; yimi = yimmi[k]; yilo = yimlo[k];

   tdg_mul(xrhi,xrmi,xrlo,yrhi,yrmi,yrlo,&zrhi,&zrmi,&zrlo);     // zr = xr*yr
   tdg_mul(xihi,ximi,xilo,yihi,yimi,yilo,&acchi,&accmi,&acclo); // acc = xi*yi
   tdg_minus(&acchi,&accmi,&acclo);
   tdg_inc(&zrhi,&zrmi,&zrlo,acchi,accmi,acclo);         // zr = xr*yr - xi*yi
   tdg_mul(xrhi,xrmi,xrlo,yihi,yimi,yilo,&zihi,&zimi,&zilo);     // zi = xr*yi
   tdg_mul(xihi,ximi,xilo,yrhi,yrmi,yrlo,&acchi,&accmi,&acclo); // acc = xi*yr
   tdg_inc(&zihi,&zimi,&zilo,acchi,accmi,acclo);         // zr = xr*yr + xi*yi

   zrehi[k] = zrhi; zremi[k] = zrmi; zrelo[k] = zrlo;
   zimhi[k] = zihi; zimmi[k] = zimi; zimlo[k] = zilo;

   for(int i=1; i<=k; i++) // z[k] = z[k] + x[i]*y[k-i]
   {
      xrhi = xrehi[i]; xrmi = xremi[i]; xrlo = xrelo[i];
      xihi = ximhi[i]; ximi = ximmi[i]; xilo = ximlo[i];
      yrhi = yrehi[k-i]; yrmi = yremi[k-i]; yrlo = yrelo[k-i];
      yihi = yimhi[k-i]; yimi = yimmi[k-i]; yilo = yimlo[k-i];

      tdg_mul(xrhi,xrmi,xrlo,yrhi,yrmi,yrlo,&zrhi,&zrmi,&zrlo); // zr = xr*yr
      tdg_mul(xihi,ximi,xilo,yihi,yimi,yilo,&acchi,&accmi,&acclo);   // xi*yi
      tdg_minus(&acchi,&accmi,&acclo);
      tdg_inc(&zrhi,&zrmi,&zrlo,acchi,accmi,acclo);     // zr = xr*yr - xi*yi
      tdg_mul(xrhi,xrmi,xrlo,yihi,yimi,yilo,&zihi,&zimi,&zilo); // zi = xr*yi
      tdg_mul(xihi,ximi,xilo,yrhi,yrmi,yrlo,&acchi,&accmi,&acclo);   // xi*yr
      tdg_inc(&zihi,&zimi,&zilo,acchi,accmi,acclo);     // zr = xr*yr + xi*yi
      // zvre[k] += zr; zvim[k] += zi
      tdg_inc(&zrehi[k],&zremi[k],&zrelo[k],zrhi,zrmi,zrlo);
      tdg_inc(&zimhi[k],&zimmi[k],&zimlo[k],zihi,zimi,zilo);
   }
}

__global__ void GPU_dbl3_speel
 ( int nvr, int deg, int *idx, double *cffhi, double *cffmi, double *cfflo,
   double *inputhi, double *inputmi, double *inputlo,
   double *forwardhi, double *forwardmi, double *forwardlo,
   double *backwardhi, double *backwardmi, double *backwardlo,
   double *crosshi, double *crossmi, double *crosslo )
{
   const int k = threadIdx.x;
   const int deg1 = deg+1;
   int ix1,ix2;

   __shared__ double xvhi[td_shmemsize];
   __shared__ double xvmi[td_shmemsize];
   __shared__ double xvlo[td_shmemsize];
   __shared__ double yvhi[td_shmemsize];
   __shared__ double yvmi[td_shmemsize];
   __shared__ double yvlo[td_shmemsize];
   __shared__ double zvhi[td_shmemsize];
   __shared__ double zvmi[td_shmemsize];
   __shared__ double zvlo[td_shmemsize];
  
   xvhi[k] = cffhi[k]; xvmi[k] = cffmi[k]; xvlo[k] = cfflo[k];
   ix1 = idx[0]*deg1+k;
   yvhi[k] = inputhi[ix1]; yvmi[k] = inputmi[ix1]; yvlo[k] = inputlo[ix1]; 
   __syncthreads();
   dbl3_convolute(xvhi,xvmi,xvlo,yvhi,yvmi,yvlo,zvhi,zvmi,zvlo,deg1,k);
   forwardhi[k] = zvhi[k];
   forwardmi[k] = zvmi[k];
   forwardlo[k] = zvlo[k];                            // f[0] = cff*x[0]

   for(int i=1; i<nvr; i++)
   {
      xvhi[k] = zvhi[k]; xvmi[k] = zvmi[k]; xvlo[k] = zvlo[k];
      ix2 = idx[i]*deg1+k;
      yvhi[k] = inputhi[ix2]; yvmi[k] = inputmi[ix2]; yvlo[k] = inputlo[ix2];
      __syncthreads();
      dbl3_convolute(xvhi,xvmi,xvlo,yvhi,yvmi,yvlo,zvhi,zvmi,zvlo,deg1,k);
      forwardhi[i*deg1+k] = zvhi[k]; 
      forwardmi[i*deg1+k] = zvmi[k]; 
      forwardlo[i*deg1+k] = zvlo[k];                  // f[i] = f[i-1]*x[i]
   }
   if(nvr > 2)
   {
      ix1 = idx[nvr-1]*deg1+k;
      xvhi[k] = inputhi[ix1]; xvmi[k] = inputmi[ix1]; xvlo[k] = inputlo[ix1];
      ix2 = idx[nvr-2]*deg1+k;
      yvhi[k] = inputhi[ix2]; yvmi[k] = inputmi[ix2]; yvlo[k] = inputlo[ix2];
      __syncthreads();
      dbl3_convolute(xvhi,xvmi,xvlo,yvhi,yvmi,yvlo,zvhi,zvmi,zvlo,deg1,k);
      backwardhi[k] = zvhi[k];
      backwardmi[k] = zvmi[k];
      backwardlo[k] = zvlo[k];                       // b[0] = x[n-1]*x[n-2]
      for(int i=1; i<nvr-2; i++)
      {
         xvhi[k] = zvhi[k]; xvmi[k] = zvmi[k]; xvlo[k] = zvlo[k];
         ix2 = idx[nvr-2-i]*deg1+k;
         yvhi[k] = inputhi[ix2]; yvmi[k] = inputmi[ix2];
         yvlo[k] = inputlo[ix2];
         __syncthreads();
         dbl3_convolute(xvhi,xvmi,xvlo,yvhi,yvmi,yvlo,zvhi,zvmi,zvlo,deg1,k);
         backwardhi[i*deg1+k] = zvhi[k];
         backwardmi[i*deg1+k] = zvmi[k];
         backwardlo[i*deg1+k] = zvlo[k];             // b[i] = b[i-1]*x[n-2-i]
      }
      xvhi[k] = zvhi[k];  xvmi[k] = zvmi[k];  xvlo[k] = zvlo[k];
      yvhi[k] = cffhi[k]; yvmi[k] = cffmi[k]; yvlo[k] = cfflo[k];
      __syncthreads();
      dbl3_convolute(xvhi,xvmi,xvlo,yvhi,yvmi,yvlo,zvhi,zvmi,zvlo,deg1,k);
      ix2 = (nvr-3)*deg1+k;
      backwardhi[ix2] = zvhi[k];
      backwardmi[ix2] = zvmi[k];
      backwardlo[ix2] = zvlo[k];                    // b[n-3] = b[n-3]*cff

      if(nvr == 3)
      {
         xvhi[k] = forwardhi[k]; xvmi[k] = forwardmi[k];
         xvlo[k] = forwardlo[k];
         ix2 = idx[2]*deg1+k;
         yvhi[k] = inputhi[ix2]; yvmi[k] = inputmi[ix2]; 
         yvlo[k] = inputlo[ix2];
         __syncthreads();
         dbl3_convolute(xvhi,xvmi,xvlo,yvhi,yvmi,yvlo,zvhi,zvmi,zvlo,deg1,k);
         crosshi[k] = zvhi[k];
         crossmi[k] = zvmi[k];
         crosslo[k] = zvlo[k];                      // c[0] = f[0]*x[2]
      }
      else
      {
         for(int i=0; i<nvr-3; i++)
         {
            ix1 = i*deg1+k; 
            xvhi[k] = forwardhi[ix1]; xvmi[k] = forwardmi[ix1];
            xvlo[k] = forwardlo[ix1];
            ix2 = (nvr-4-i)*deg1+k;
            yvhi[k] = backwardhi[ix2]; yvmi[k] = backwardmi[ix2];
            yvlo[k] = backwardlo[ix2];
            __syncthreads();
            dbl3_convolute
               (xvhi,xvmi,xvlo,yvhi,yvmi,yvlo,zvhi,zvmi,zvlo,deg1,k);
            crosshi[i*deg1+k] = zvhi[k];
            crossmi[i*deg1+k] = zvmi[k];
            crosslo[i*deg1+k] = zvlo[k];            // c[i] = f[i]*b[n-4-i]
         }
         ix1 = (nvr-3)*deg1+k;
         xvhi[k] = forwardhi[ix1]; xvmi[k] = forwardmi[ix1];
         xvlo[k] = forwardlo[ix1];
         ix2 = idx[nvr-1]*deg1+k;
         yvhi[k] = inputhi[ix2]; yvmi[k] = inputmi[ix2];
         yvlo[k] = inputlo[ix2];
         __syncthreads();
         dbl3_convolute(xvhi,xvmi,xvlo,yvhi,yvmi,yvlo,zvhi,zvmi,zvlo,deg1,k);
         crosshi[(nvr-3)*deg1+k] = zvhi[k];
         crossmi[(nvr-3)*deg1+k] = zvmi[k];
         crosslo[(nvr-3)*deg1+k] = zvlo[k];         // c[n-3] = f[n-3]*x[n-1]
      }
   }
}

__global__ void GPU_cmplx3_speel
 ( int nvr, int deg, int *idx,
   double *cffrehi, double *cffremi, double *cffrelo,
   double *cffimhi, double *cffimmi, double *cffimlo,
   double *inputrehi, double *inputremi, double *inputrelo,
   double *inputimhi, double *inputimmi, double *inputimlo,
   double *forwardrehi, double *forwardremi, double *forwardrelo,
   double *forwardimhi, double *forwardimmi, double *forwardimlo,
   double *backwardrehi, double *backwardremi, double *backwardrelo,
   double *backwardimhi, double *backwardimmi, double *backwardimlo,
   double *crossrehi, double *crossremi, double *crossrelo,
   double *crossimhi, double *crossimmi, double *crossimlo )
{
   const int k = threadIdx.x;
   const int deg1 = deg+1;
   int ix1,ix2;

   __shared__ double xvrehi[td_shmemsize];
   __shared__ double xvremi[td_shmemsize];
   __shared__ double xvrelo[td_shmemsize];
   __shared__ double xvimhi[td_shmemsize];
   __shared__ double xvimmi[td_shmemsize];
   __shared__ double xvimlo[td_shmemsize];
   __shared__ double yvrehi[td_shmemsize];
   __shared__ double yvremi[td_shmemsize];
   __shared__ double yvrelo[td_shmemsize];
   __shared__ double yvimhi[td_shmemsize];
   __shared__ double yvimmi[td_shmemsize];
   __shared__ double yvimlo[td_shmemsize];
   __shared__ double zvrehi[td_shmemsize];
   __shared__ double zvremi[td_shmemsize];
   __shared__ double zvrelo[td_shmemsize];
   __shared__ double zvimhi[td_shmemsize];
   __shared__ double zvimmi[td_shmemsize];
   __shared__ double zvimlo[td_shmemsize];

   xvrehi[k] = cffrehi[k]; xvremi[k] = cffremi[k]; xvrelo[k] = cffrelo[k];
   xvimhi[k] = cffimhi[k]; xvimmi[k] = cffimmi[k]; xvimlo[k] = cffimlo[k];
   ix1 = idx[0]*deg1+k;
   yvrehi[k] = inputrehi[ix1]; yvremi[k] = inputremi[ix1];
   yvrelo[k] = inputrelo[ix1];
   yvimhi[k] = inputimhi[ix1]; yvimmi[k] = inputimmi[ix1];
   yvimlo[k] = inputimlo[ix1];
   __syncthreads();                                       // f[0] = cff*x[0] 
   cmplx3_convolute(xvrehi,xvremi,xvrelo,xvimhi,xvimmi,xvimlo,
                    yvrehi,yvremi,yvrelo,yvimhi,yvimmi,yvimlo,
                    zvrehi,zvremi,zvrelo,zvimhi,zvimmi,zvimlo,deg1,k);
   forwardrehi[k] = zvrehi[k]; forwardremi[k] = zvremi[k];
   forwardrelo[k] = zvrelo[k];
   forwardimhi[k] = zvimhi[k]; forwardimmi[k] = zvimmi[k];
   forwardimlo[k] = zvimlo[k];

   for(int i=1; i<nvr; i++)
   {
      xvrehi[k] = zvrehi[k]; xvremi[k] = zvremi[k]; xvrelo[k] = zvrelo[k];
      xvimhi[k] = zvimhi[k]; xvimmi[k] = zvimmi[k]; xvimlo[k] = zvimlo[k];
      ix2 = idx[i]*deg1+k;
      yvrehi[k] = inputrehi[ix2]; yvremi[k] = inputremi[ix2];
      yvrelo[k] = inputrelo[ix2];
      yvimhi[k] = inputimhi[ix2]; yvimmi[k] = inputimmi[ix2];
      yvimlo[k] = inputimlo[ix2];
      __syncthreads();                                 // f[i] = f[i-i]*x[i]
      cmplx3_convolute(xvrehi,xvremi,xvrelo,xvimhi,xvimmi,xvimlo,
                       yvrehi,yvremi,yvrelo,yvimhi,yvimmi,yvimlo,
                       zvrehi,zvremi,zvrelo,zvimhi,zvimmi,zvimlo,deg1,k);
      ix1 = i*deg1+k;                                   
      forwardrehi[ix1] = zvrehi[k]; forwardremi[ix1] = zvremi[k];
      forwardrelo[ix1] = zvrelo[k];
      forwardimhi[ix1] = zvimhi[k]; forwardimmi[ix1] = zvimmi[k];
      forwardimlo[ix1] = zvimlo[k]; 
   }
   if(nvr > 2)
   {
      ix1 = idx[nvr-1]*deg1+k;
      xvrehi[k] = inputrehi[ix1]; xvremi[k] = inputremi[ix1];
      xvrelo[k] = inputrelo[ix1];
      xvimhi[k] = inputimhi[ix1]; xvimmi[k] = inputimmi[ix1];
      xvimlo[k] = inputimlo[ix1];
      ix2 = idx[nvr-2]*deg1+k;
      yvrehi[k] = inputrehi[ix2]; yvremi[k] = inputremi[ix2];
      yvrelo[k] = inputrelo[ix2];
      yvimhi[k] = inputimhi[ix2]; yvimmi[k] = inputimmi[ix2];
      yvimlo[k] = inputimlo[ix2];
      __syncthreads();                               // b[0] = x[n-1]*x[n-2]
      cmplx3_convolute(xvrehi,xvremi,xvrelo,xvimhi,xvimmi,xvimlo,
                       yvrehi,yvremi,yvrelo,yvimhi,yvimmi,yvimlo,
                       zvrehi,zvremi,zvrelo,zvimhi,zvimmi,zvimlo,deg1,k);
      backwardrehi[k] = zvrehi[k]; backwardremi[k] = zvremi[k];
      backwardrelo[k] = zvrelo[k];
      backwardimhi[k] = zvimhi[k]; backwardimmi[k] = zvimmi[k];
      backwardimlo[k] = zvimlo[k];

      for(int i=1; i<nvr-2; i++)
      {
         xvrehi[k] = zvrehi[k]; xvremi[k] = zvremi[k]; xvrelo[k] = zvrelo[k];
         xvimhi[k] = zvimhi[k]; xvimmi[k] = zvimmi[k]; xvimlo[k] = zvimlo[k];
         ix2 = idx[nvr-2-i]*deg1+k;
         yvrehi[k] = inputrehi[ix2]; yvremi[k] = inputremi[ix2];
         yvrelo[k] = inputrelo[ix2];
         yvimhi[k] = inputimhi[ix2]; yvimmi[k] = inputimmi[ix2];
         yvimlo[k] = inputimlo[ix2];
         __syncthreads();                           // b[i] = b[i]*x[n-2-i]
         cmplx3_convolute(xvrehi,xvremi,xvrelo,xvimhi,xvimmi,xvimlo,
                          yvrehi,yvremi,yvrelo,yvimhi,yvimmi,yvimlo,
                          zvrehi,zvremi,zvrelo,zvimhi,zvimmi,zvimlo,deg1,k);
         ix1 = i*deg1+k;
         backwardrehi[ix1] = zvrehi[k]; backwardremi[ix1] = zvremi[k];
         backwardrelo[ix1] = zvrelo[k];
         backwardimhi[ix1] = zvimhi[k]; backwardimmi[ix1] = zvimmi[k];
         backwardimlo[ix1] = zvimlo[k];
      }
      xvrehi[k] = zvrehi[k]; xvremi[k] = zvremi[k]; xvrelo[k] = zvrelo[k];
      xvimhi[k] = zvimhi[k]; xvimmi[k] = zvimmi[k]; xvimlo[k] = zvimlo[k];
      yvrehi[k] = cffrehi[k]; yvremi[k] = cffremi[k]; yvrelo[k] = cffrelo[k];
      yvimhi[k] = cffimhi[k]; yvimmi[k] = cffimmi[k]; yvimlo[k] = cffimlo[k];
      __syncthreads();                               // b[n-3] = b[n-3]*cff
      cmplx3_convolute(xvrehi,xvremi,xvrelo,xvimhi,xvimmi,xvimlo,
                       yvrehi,yvremi,yvrelo,yvimhi,yvimmi,yvimlo,
                       zvrehi,zvremi,zvrelo,zvimhi,zvimmi,zvimlo,deg1,k);
      ix1 = (nvr-3)*deg1+k;
      backwardrehi[ix1] = zvrehi[k]; backwardremi[ix1] = zvremi[k];
      backwardrelo[ix1] = zvrelo[k];
      backwardimhi[ix1] = zvimhi[k]; backwardimmi[ix1] = zvimmi[k];
      backwardimlo[ix1] = zvimlo[k];

      if(nvr == 3)
      {
         xvrehi[k] = forwardrehi[k]; xvremi[k] = forwardremi[k];
         xvrelo[k] = forwardrelo[k];
         xvimhi[k] = forwardimhi[k]; xvimmi[k] = forwardimmi[k];
         xvimlo[k] = forwardimlo[k];
         ix2 = idx[2]*deg1+k;
         yvrehi[k] = inputrehi[ix2]; yvremi[k] = inputremi[ix2];
         yvrelo[k] = inputrelo[ix2];
         yvimhi[k] = inputimhi[ix2]; yvimmi[k] = inputimmi[ix2];
         yvimlo[k] = inputimlo[ix2];
         __syncthreads();                               // c[0] = f[0]*x[2]
         cmplx3_convolute(xvrehi,xvremi,xvrelo,xvimhi,xvimmi,xvimlo,
                          yvrehi,yvremi,yvrelo,yvimhi,yvimmi,yvimlo,
                          zvrehi,zvremi,zvrelo,zvimhi,zvimmi,zvimlo,deg1,k);
         crossrehi[k] = zvrehi[k]; crossremi[k] = zvremi[k];
         crossrelo[k] = zvrelo[k];
         crossimhi[k] = zvimhi[k]; crossimmi[k] = zvimmi[k];
         crossimlo[k] = zvimlo[k];
      }
      else
      {
         for(int i=0; i<nvr-3; i++)
         {
            ix1 = i*deg1+k;   
            xvrehi[k] = forwardrehi[ix1]; xvremi[k] = forwardremi[ix1];
            xvrelo[k] = forwardrelo[ix1];
            xvimhi[k] = forwardimhi[ix1]; xvimmi[k] = forwardimmi[ix1];
            xvimlo[k] = forwardimlo[ix1];
            ix2 = (nvr-4-i)*deg1+k;
            yvrehi[k] = backwardrehi[ix2]; yvremi[k] = backwardremi[ix2];
            yvrelo[k] = backwardrelo[ix2];
            yvimhi[k] = backwardimhi[ix2]; yvimmi[k] = backwardimmi[ix2];
            yvimlo[k] = backwardimlo[ix2];
            __syncthreads();                        // c[i] = f[i]*b[n-4-i]
            cmplx3_convolute(xvrehi,xvremi,xvrelo,xvimhi,xvimmi,xvimlo,
                             yvrehi,yvremi,yvrelo,yvimhi,yvimmi,yvimlo,
                             zvrehi,zvremi,zvrelo,zvimhi,zvimmi,zvimlo,deg1,k);
            ix1 = i*deg1+k;
            crossrehi[ix1] = zvrehi[k]; crossremi[ix1] = zvremi[k];
            crossrelo[ix1] = zvrelo[k];
            crossimhi[ix1] = zvimhi[k]; crossimmi[ix1] = zvimmi[k];
            crossimlo[ix1] = zvimlo[k];
         }
         ix1 = (nvr-3)*deg1+k;
         xvrehi[k] = forwardrehi[ix1]; xvremi[k] = forwardremi[ix1];
         xvrelo[k] = forwardrelo[ix1];
         xvimhi[k] = forwardimhi[ix1]; xvimmi[k] = forwardimmi[ix1];
         xvimlo[k] = forwardimlo[ix1];
         ix2 = idx[nvr-1]*deg1+k;
         yvrehi[k] = inputrehi[ix2]; yvremi[k] = inputremi[ix2];
         yvrelo[k] = inputrelo[ix2];
         yvimhi[k] = inputimhi[ix2]; yvimmi[k] = inputimmi[ix2];
         yvimlo[k] = inputimlo[ix2];
         __syncthreads();                         // c[n-3] = f[n-3]*x[n-1]
         cmplx3_convolute(xvrehi,xvremi,xvrelo,xvimhi,xvimmi,xvimlo,
                          yvrehi,yvremi,yvrelo,yvimhi,yvimmi,yvimlo,
                          zvrehi,zvremi,zvrelo,zvimhi,zvimmi,zvimlo,deg1,k);
         ix1 = (nvr-3)*deg1+k;
         crossrehi[ix1] = zvrehi[k]; crossremi[ix1] = zvremi[k];
         crossrelo[ix1] = zvrelo[k];
         crossimhi[ix1] = zvimhi[k]; crossimmi[ix1] = zvimmi[k];
         crossimlo[ix1] = zvimlo[k];
      }
   }
}

void GPU_dbl3_evaldiff
 ( int BS, int dim, int nvr, int deg, int *idx,
   double *cffhi, double *cffmi, double *cfflo,
   double **inputhi, double **inputmi, double **inputlo,
   double **outputhi, double **outputmi, double **outputlo )
{
   const int deg1 = deg+1;            // length of all vectors
   double *inputhi_d;                 // inputhi_d is input on the device
   double *inputmi_d;                 // inputmi_d is input on the device
   double *inputlo_d;                 // inputlo_d is input on the device
   double *forwardhi_d;               // high forward products on the device
   double *forwardmi_d;               // middle forward products on the device
   double *forwardlo_d;               // low forward products on the device
   double *backwardhi_d;              // high backward products on the device
   double *backwardmi_d;              // middle backward products on the device
   double *backwardlo_d;              // low backward products on the device
   double *crosshi_d;                 // high cross products on the device
   double *crossmi_d;                 // middle cross products on the device
   double *crosslo_d;                 // low cross products on the device
   double *cffhi_d;                   // cffhi_d is cffhi on device
   double *cffmi_d;                   // cffmi_d is cffmi on device
   double *cfflo_d;                   // cfflo_d is cfflo on device
   int *idx_d;                        // idx_d is idx on device

   size_t szcff = deg1*sizeof(double);
   size_t szdim = dim*(deg1)*sizeof(double);
   size_t sznvr = nvr*(deg1)*sizeof(double);
   size_t sznvr2 = (nvr-2)*(deg1)*sizeof(double);
   size_t szidx = nvr*sizeof(int);

   cudaMalloc((void**)&idx_d,szidx);
   cudaMalloc((void**)&cffhi_d,szcff);
   cudaMalloc((void**)&cffmi_d,szcff);
   cudaMalloc((void**)&cfflo_d,szcff);
   cudaMalloc((void**)&inputhi_d,szdim);
   cudaMalloc((void**)&inputmi_d,szdim);
   cudaMalloc((void**)&inputlo_d,szdim);
   cudaMalloc((void**)&forwardhi_d,sznvr);
   cudaMalloc((void**)&forwardmi_d,sznvr);
   cudaMalloc((void**)&forwardlo_d,sznvr);
   cudaMalloc((void**)&backwardhi_d,sznvr2);
   cudaMalloc((void**)&backwardmi_d,sznvr2);
   cudaMalloc((void**)&backwardlo_d,sznvr2);
   cudaMalloc((void**)&crosshi_d,sznvr2);
   cudaMalloc((void**)&crossmi_d,sznvr2);
   cudaMalloc((void**)&crosslo_d,sznvr2);

   double *inputhi_h = new double[dim*(deg1)];
   double *inputmi_h = new double[dim*(deg1)];
   double *inputlo_h = new double[dim*(deg1)];
   int ix = 0;
   for(int i=0; i<dim; i++)
      for(int j=0; j<deg1; j++)
      {
         inputhi_h[ix] = inputhi[i][j];
         inputmi_h[ix] = inputmi[i][j];
         inputlo_h[ix++] = inputlo[i][j];
      }

   cudaMemcpy(idx_d,idx,szidx,cudaMemcpyHostToDevice);
   cudaMemcpy(cffhi_d,cffhi,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffmi_d,cffmi,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cfflo_d,cfflo,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(inputhi_d,inputhi_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputmi_d,inputmi_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputlo_d,inputlo_h,szdim,cudaMemcpyHostToDevice);

   if(BS = deg1)
   {
      GPU_dbl3_speel<<<1,BS>>>
         (nvr,deg,idx_d,cffhi_d,cffmi_d,cfflo_d,inputhi_d,inputmi_d,
          inputlo_d,forwardhi_d,forwardmi_d,forwardlo_d,backwardhi_d,
          backwardmi_d,backwardlo_d,crosshi_d,crossmi_d,crosslo_d);
   }
   double *forwardhi_h = new double[(deg1)*nvr];
   double *forwardmi_h = new double[(deg1)*nvr];
   double *forwardlo_h = new double[(deg1)*nvr];
   double *backwardhi_h = new double[(deg1)*(nvr-2)];
   double *backwardmi_h = new double[(deg1)*(nvr-2)];
   double *backwardlo_h = new double[(deg1)*(nvr-2)];
   double *crosshi_h = new double[(deg1)*(nvr-2)];
   double *crossmi_h = new double[(deg1)*(nvr-2)];
   double *crosslo_h = new double[(deg1)*(nvr-2)];
  
   cudaMemcpy(forwardhi_h,forwardhi_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardmi_h,forwardmi_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardlo_h,forwardlo_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardhi_h,backwardhi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardmi_h,backwardmi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardlo_h,backwardlo_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crosshi_h,crosshi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossmi_h,crossmi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crosslo_h,crosslo_d,sznvr2,cudaMemcpyDeviceToHost);

   int offset = (nvr-1)*deg1;            // assign value of the monomial
   for(int i=0; i<deg1; i++)
   {
      outputhi[dim][i] = forwardhi_h[offset+i];
      outputmi[dim][i] = forwardmi_h[offset+i];
      outputlo[dim][i] = forwardlo_h[offset+i];
   }
   ix = idx[nvr-1];                      // derivative with respect to x[n-1]
   offset = (nvr-2)*deg1;
   for(int i=0; i<deg1; i++)
   {
      outputhi[ix][i] = forwardhi_h[offset+i];
      outputmi[ix][i] = forwardmi_h[offset+i];
      outputlo[ix][i] = forwardlo_h[offset+i];
   }
   ix = idx[0];                          // derivative with respect to x[0]
   offset = (nvr-3)*deg1;
   for(int i=0; i<deg1; i++)
   {
      outputhi[ix][i] = backwardhi_h[offset+i];
      outputmi[ix][i] = backwardmi_h[offset+i];
      outputlo[ix][i] = backwardlo_h[offset+i];
   }
   for(int k=1; k<nvr-1; k++)            // derivative with respect to x[k]
   {
      ix = idx[k]; offset = (k-1)*deg1;
      for(int i=0; i<deg1; i++)
      {
         outputhi[ix][i] = crosshi_h[offset+i];
         outputmi[ix][i] = crossmi_h[offset+i];
         outputlo[ix][i] = crosslo_h[offset+i];
      }
   }
}

void GPU_cmplx3_evaldiff
 ( int BS, int dim, int nvr, int deg, int *idx,
   double *cffrehi, double *cffremi, double *cffrelo,
   double *cffimhi, double *cffimmi, double *cffimlo,
   double **inputrehi, double **inputremi, double **inputrelo,
   double **inputimhi, double **inputimmi, double **inputimlo,
   double **outputrehi, double **outputremi, double **outputrelo,
   double **outputimhi, double **outputimmi, double **outputimlo )
{
   const int deg1 = deg+1;          // length of all vectors
   double *inputrehi_d;             // inputrehi_d is inputrehi on the device
   double *inputremi_d;             // inputremi_d is inputrehi on the device
   double *inputrelo_d;             // inputrelo_d is inputrelo on the device
   double *inputimhi_d;             // inputimhi_d is inputrehi on the device
   double *inputimmi_d;             // inputimmi_d is inputrehi on the device
   double *inputimlo_d;             // inputimlo_d is inputrelo on the device
   double *forwardrehi_d;
   double *forwardremi_d;
   double *forwardrelo_d;
   double *forwardimhi_d;
   double *forwardimmi_d;
   double *forwardimlo_d;
   double *backwardrehi_d;
   double *backwardremi_d;
   double *backwardrelo_d;
   double *backwardimhi_d;
   double *backwardimmi_d;
   double *backwardimlo_d;
   double *crossrehi_d;
   double *crossremi_d;
   double *crossrelo_d;
   double *crossimhi_d;
   double *crossimmi_d;
   double *crossimlo_d;
   double *cffrehi_d;               // cffrehi_d is cffrehi on the device
   double *cffremi_d;               // cffremi_d is cffrehi on the device
   double *cffrelo_d;               // cffrelo_d is cffrelo on the device
   double *cffimhi_d;               // cffimhi_d is cffimhi on the device
   double *cffimmi_d;               // cffimmi_d is cffimhi on the device
   double *cffimlo_d;               // cffimlo_d is cffimlo on the device
   int *idx_d;                      // idx_d is idx on the device

   size_t szdim = dim*(deg1)*sizeof(double);
   size_t sznvr = nvr*(deg1)*sizeof(double);
   size_t sznvr2 = (nvr-2)*(deg1)*sizeof(double);
   size_t szidx = nvr*sizeof(int);
   size_t szcff = deg1*sizeof(double);

   cudaMalloc((void**)&idx_d,szidx);
   cudaMalloc((void**)&cffrehi_d,szcff);
   cudaMalloc((void**)&cffremi_d,szcff);
   cudaMalloc((void**)&cffrelo_d,szcff);
   cudaMalloc((void**)&cffimhi_d,szcff);
   cudaMalloc((void**)&cffimmi_d,szcff);
   cudaMalloc((void**)&cffimlo_d,szcff);
   cudaMalloc((void**)&inputrehi_d,szdim);
   cudaMalloc((void**)&inputremi_d,szdim);
   cudaMalloc((void**)&inputrelo_d,szdim);
   cudaMalloc((void**)&inputimhi_d,szdim);
   cudaMalloc((void**)&inputimmi_d,szdim);
   cudaMalloc((void**)&inputimlo_d,szdim);
   cudaMalloc((void**)&forwardrehi_d,sznvr);
   cudaMalloc((void**)&forwardremi_d,sznvr);
   cudaMalloc((void**)&forwardrelo_d,sznvr);
   cudaMalloc((void**)&forwardimhi_d,sznvr);
   cudaMalloc((void**)&forwardimmi_d,sznvr);
   cudaMalloc((void**)&forwardimlo_d,sznvr);
   cudaMalloc((void**)&backwardrehi_d,sznvr2);
   cudaMalloc((void**)&backwardremi_d,sznvr2);
   cudaMalloc((void**)&backwardrelo_d,sznvr2);
   cudaMalloc((void**)&backwardimhi_d,sznvr2);
   cudaMalloc((void**)&backwardimmi_d,sznvr2);
   cudaMalloc((void**)&backwardimlo_d,sznvr2);
   cudaMalloc((void**)&crossrehi_d,sznvr2);
   cudaMalloc((void**)&crossremi_d,sznvr2);
   cudaMalloc((void**)&crossrelo_d,sznvr2);
   cudaMalloc((void**)&crossimhi_d,sznvr2);
   cudaMalloc((void**)&crossimmi_d,sznvr2);
   cudaMalloc((void**)&crossimlo_d,sznvr2);

   double *inputrehi_h = new double[dim*(deg1)];
   double *inputremi_h = new double[dim*(deg1)];
   double *inputrelo_h = new double[dim*(deg1)];
   double *inputimhi_h = new double[dim*(deg1)];
   double *inputimmi_h = new double[dim*(deg1)];
   double *inputimlo_h = new double[dim*(deg1)];
   int ix = 0;
   for(int i=0; i<dim; i++)
      for(int j=0; j<deg1; j++)
      {
         inputrehi_h[ix] = inputrehi[i][j];
         inputremi_h[ix] = inputremi[i][j];
         inputrelo_h[ix] = inputrelo[i][j];
         inputimhi_h[ix] = inputimhi[i][j];
         inputimmi_h[ix] = inputimmi[i][j];
         inputimlo_h[ix++] = inputimlo[i][j];
      }

   cudaMemcpy(idx_d,idx,szidx,cudaMemcpyHostToDevice);
   cudaMemcpy(cffrehi_d,cffrehi,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffremi_d,cffremi,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffrelo_d,cffrelo,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffimhi_d,cffimhi,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffimmi_d,cffimmi,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffimlo_d,cffimlo,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(inputrehi_d,inputrehi_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputremi_d,inputremi_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputrelo_d,inputrelo_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputimhi_d,inputimhi_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputimmi_d,inputimmi_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputimlo_d,inputimlo_h,szdim,cudaMemcpyHostToDevice);

   if(BS = deg1)
   {
      GPU_cmplx3_speel<<<1,BS>>>(nvr,deg,idx_d,
         cffrehi_d,cffremi_d,cffrelo_d,cffimhi_d,cffimmi_d,cffimlo_d,
         inputrehi_d,inputremi_d,inputrelo_d,
         inputimhi_d,inputimmi_d,inputimlo_d,
         forwardrehi_d,forwardremi_d,forwardrelo_d,
         forwardimhi_d,forwardimmi_d,forwardimlo_d,
         backwardrehi_d,backwardremi_d,backwardrelo_d,
         backwardimhi_d,backwardimmi_d,backwardimlo_d,
         crossrehi_d,crossremi_d,crossrelo_d,
         crossimhi_d,crossimmi_d,crossimlo_d);
   }
   double *forwardrehi_h = new double[(deg1)*nvr];
   double *forwardremi_h = new double[(deg1)*nvr];
   double *forwardrelo_h = new double[(deg1)*nvr];
   double *forwardimhi_h = new double[(deg1)*nvr];
   double *forwardimmi_h = new double[(deg1)*nvr];
   double *forwardimlo_h = new double[(deg1)*nvr];
   double *backwardrehi_h = new double[(deg1)*(nvr-2)];
   double *backwardremi_h = new double[(deg1)*(nvr-2)];
   double *backwardrelo_h = new double[(deg1)*(nvr-2)];
   double *backwardimhi_h = new double[(deg1)*(nvr-2)];
   double *backwardimmi_h = new double[(deg1)*(nvr-2)];
   double *backwardimlo_h = new double[(deg1)*(nvr-2)];
   double *crossrehi_h = new double[(deg1)*(nvr-2)];
   double *crossremi_h = new double[(deg1)*(nvr-2)];
   double *crossrelo_h = new double[(deg1)*(nvr-2)];
   double *crossimhi_h = new double[(deg1)*(nvr-2)];
   double *crossimmi_h = new double[(deg1)*(nvr-2)];
   double *crossimlo_h = new double[(deg1)*(nvr-2)];
  
   cudaMemcpy(forwardrehi_h,forwardrehi_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardremi_h,forwardremi_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardrelo_h,forwardrelo_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardimhi_h,forwardimhi_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardimmi_h,forwardimmi_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardimlo_h,forwardimlo_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardrehi_h,backwardrehi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardremi_h,backwardremi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardrelo_h,backwardrelo_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardimhi_h,backwardimhi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardimmi_h,backwardimmi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardimlo_h,backwardimlo_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossrehi_h,crossrehi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossremi_h,crossremi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossrelo_h,crossrelo_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossimhi_h,crossimhi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossimmi_h,crossimmi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossimlo_h,crossimlo_d,sznvr2,cudaMemcpyDeviceToHost);

   int offset = (nvr-1)*deg1;
   for(int i=0; i<deg1; i++)   // assign value of the monomial
   {
      outputrehi[dim][i] = forwardrehi_h[offset+i];
      outputremi[dim][i] = forwardremi_h[offset+i];
      outputrelo[dim][i] = forwardrelo_h[offset+i];
      outputimhi[dim][i] = forwardimhi_h[offset+i];
      outputimmi[dim][i] = forwardimmi_h[offset+i];
      outputimlo[dim][i] = forwardimlo_h[offset+i];
   }
   ix = idx[nvr-1];
   offset = (nvr-2)*deg1;
   for(int i=0; i<deg1; i++)  // derivative with respect to x[n-1]
   {
      outputrehi[ix][i] = forwardrehi_h[offset+i];
      outputremi[ix][i] = forwardremi_h[offset+i];
      outputrelo[ix][i] = forwardrelo_h[offset+i];
      outputimhi[ix][i] = forwardimhi_h[offset+i];
      outputimmi[ix][i] = forwardimmi_h[offset+i];
      outputimlo[ix][i] = forwardimlo_h[offset+i];
   }
   ix = idx[0]; 
   offset = (nvr-3)*deg1;
   for(int i=0; i<deg1; i++)   // derivative with respect to x[0]
   {
      outputrehi[ix][i] = backwardrehi_h[offset+i];
      outputremi[ix][i] = backwardremi_h[offset+i];
      outputrelo[ix][i] = backwardrelo_h[offset+i];
      outputimhi[ix][i] = backwardimhi_h[offset+i];
      outputimmi[ix][i] = backwardimmi_h[offset+i];
      outputimlo[ix][i] = backwardimlo_h[offset+i];
   }
   for(int k=1; k<nvr-1; k++)  // derivative with respect to x[k]
   {
      ix = idx[k]; offset = (k-1)*deg1;
      for(int i=0; i<deg1; i++)
      {
         outputrehi[ix][i] = crossrehi_h[offset+i];
         outputremi[ix][i] = crossremi_h[offset+i];
         outputrelo[ix][i] = crossrelo_h[offset+i];
         outputimhi[ix][i] = crossimhi_h[offset+i];
         outputimmi[ix][i] = crossimmi_h[offset+i];
         outputimlo[ix][i] = crossimlo_h[offset+i];
      }
   }
}
