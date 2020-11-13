// The file dbl5_monomials_kernels.cu defines the kernels specified
// in dbl5_monomials_kernels.h.

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
#include "penta_double_gpufun.cu"
#include "dbl5_convolutions_kernels.h"
#include "dbl5_monomials_kernels.h"

__device__ void dbl5_convolute
 ( double *xtb, double *xix, double *xmi, double *xrg, double *xpk,
   double *ytb, double *yix, double *ymi, double *yrg, double *ypk,
   double *ztb, double *zix, double *zmi, double *zrg, double *zpk,
   int dim, int k )
{
   double prdtb,prdix,prdmi,prdrg,prdpk;
   int idx;

   // zv[k] = xv[0]*yv[k];
   pdg_mul(xtb[0],xix[0],xmi[0],xrg[0],xpk[0],
           ytb[k],yix[k],ymi[k],yrg[k],ypk[k],
           &ztb[k],&zix[k],&zmi[k],&zrg[k],&zpk[k]);

   for(int i=1; i<=k; i++) // zv[k] = zv[k] + xv[i]*yv[k-i];
   {
      idx = k-i;
      pdg_mul(xtb[i],xix[i],xmi[i],xrg[i],xpk[i],
              ytb[idx],yix[idx],ymi[idx],yrg[idx],ypk[idx],
              &prdtb,&prdix,&prdmi,&prdrg,&prdpk);
      pdg_inc(&ztb[k],&zix[k],&zmi[k],&zrg[k],&zpk[k],
              prdtb,prdix,prdmi,prdrg,prdpk);
   }
}

__device__ void cmplx5_convolute
 ( double *xretb, double *xreix, double *xremi, double *xrerg, double *xrepk,
   double *ximtb, double *ximix, double *ximmi, double *ximrg, double *ximpk,
   double *yretb, double *yreix, double *yremi, double *yrerg, double *yrepk,
   double *yimtb, double *yimix, double *yimmi, double *yimrg, double *yimpk,
   double *zretb, double *zreix, double *zremi, double *zrerg, double *zrepk,
   double *zimtb, double *zimix, double *zimmi, double *zimrg, double *zimpk,
   int dim, int k )
{
   double xrtb,xitb,yrtb,yitb,zrtb,zitb,acctb;
   double xrix,xiix,yrix,yiix,zrix,ziix,accix;
   double xrmi,ximi,yrmi,yimi,zrmi,zimi,accmi;
   double xrrg,xirg,yrrg,yirg,zrrg,zirg,accrg;
   double xrpk,xipk,yrpk,yipk,zrpk,zipk,accpk;
   int idx;

   // z[k] = x[0]*y[k]
   xrtb = xretb[0]; xrix = xreix[0]; xrmi = xremi[0];
   xrrg = xrerg[0]; xrpk = xrepk[0];
   xitb = ximtb[0]; xiix = ximix[0]; ximi = ximmi[0];
   xirg = ximrg[0]; xipk = ximpk[0];
   yrtb = yretb[k]; yrix = yreix[k]; yrmi = yremi[k];
   yrrg = yrerg[k]; yrpk = yrepk[k];
   yitb = yimtb[k]; yiix = yimix[k]; yimi = yimmi[k];
   yirg = yimrg[k]; yipk = yimpk[k];

   pdg_mul(xrtb,xrix,xrmi,xrrg,xrpk,yrtb,yrix,yrmi,yrrg,yrpk,
           &zrtb,&zrix,&zrmi,&zrrg,&zrpk);         // zr = xr*yr
   pdg_mul(xitb,xiix,ximi,xirg,xipk,yitb,yiix,yimi,yirg,yipk,
           &acctb,&accix,&accmi,&accrg,&accpk);    // acc = xi*yi
   pdg_minus(&acctb,&accix,&accmi,&accrg,&accpk);
   pdg_inc(&zrtb,&zrix,&zrmi,&zrrg,&zrpk,
           acctb,accix,accmi,accrg,accpk);         // zr = xr*yr - xi*yi
   pdg_mul(xrtb,xrix,xrmi,xrrg,xrpk,yitb,yiix,yimi,yirg,yipk,
           &zitb,&ziix,&zimi,&zirg,&zipk);         // zi = xr*yi
   pdg_mul(xitb,xiix,ximi,xirg,xipk,yrtb,yrix,yrmi,yrrg,yrpk,
           &acctb,&accix,&accmi,&accrg,&accpk);    // acc = xi*yr
   pdg_inc(&zitb,&ziix,&zimi,&zirg,&zipk,
           acctb,accix,accmi,accrg,accpk);         // zr = xr*yr + xi*yi

   zretb[k] = zrtb; zreix[k] = zrix; zremi[k] = zrmi;
   zrerg[k] = zrrg; zrepk[k] = zrpk;
   zimtb[k] = zitb; zimix[k] = ziix; zimmi[k] = zimi;
   zimrg[k] = zirg; zimpk[k] = zipk;

   for(int i=1; i<=k; i++) // z[k] = z[k] + x[i]*y[k-i]
   {
      xrtb = xretb[i]; xrix = xreix[i]; xrmi = xremi[i];
      xrrg = xrerg[i]; xrpk = xrepk[i];
      xitb = ximtb[i]; xiix = ximix[i]; ximi = ximmi[i];
      xirg = ximrg[i]; xipk = ximpk[i];
      idx = k-i;
      yrtb = yretb[idx]; yrix = yreix[idx]; yrmi = yremi[idx];
      yrrg = yrerg[idx]; yrpk = yrepk[idx];
      yitb = yimtb[idx]; yiix = yimix[idx]; yimi = yimmi[idx];
      yirg = yimrg[idx]; yipk = yimpk[idx];

      pdg_mul(xrtb,xrix,xrmi,xrrg,xrpk,yrtb,yrix,yrmi,yrrg,yrpk,
              &zrtb,&zrix,&zrmi,&zrrg,&zrpk);        // zr = xr*yr
      pdg_mul(xitb,xiix,ximi,xirg,xipk,yitb,yiix,yimi,yirg,yipk,
              &acctb,&accix,&accmi,&accrg,&accpk);   // xi*yi
      pdg_minus(&acctb,&accix,&accmi,&accrg,&accpk);
      pdg_inc(&zrtb,&zrix,&zrmi,&zrrg,&zrpk,
              acctb,accix,accmi,accrg,accpk);        // zr = xr*yr - xi*yi
      pdg_mul(xrtb,xrix,xrmi,xrrg,xrpk,yitb,yiix,yimi,yirg,yipk,
              &zitb,&ziix,&zimi,&zirg,&zipk);        // zi = xr*yi
      pdg_mul(xitb,xiix,ximi,xirg,xipk,yrtb,yrix,yrmi,yrrg,yrpk,
              &acctb,&accix,&accmi,&accrg,&accpk);   // xi*yr
      pdg_inc(&zitb,&ziix,&zimi,&zirg,&zipk,
              acctb,accix,accmi,accrg,accpk);        // zr = xr*yr + xi*yi
      // zre[k] += zr; zim[k] += zi
      pdg_inc(&zretb[k],&zreix[k],&zremi[k],&zrerg[k],&zrepk[k],
              zrtb,zrix,zrmi,zrrg,zrpk);
      pdg_inc(&zimtb[k],&zimix[k],&zimmi[k],&zimrg[k],&zimpk[k],
              zitb,ziix,zimi,zirg,zipk);
   }
}

__global__ void GPU_dbl5_speel
 ( int nvr, int deg, int *idx,
   double *cfftb, double *cffix, double *cffmi, double *cffrg, double *cffpk,
   double *inputtb, double *inputix, double *inputmi,
   double *inputrg, double *inputpk,
   double *forwardtb, double *forwardix, double *forwardmi,
   double *forwardrg, double *forwardpk,
   double *backwardtb, double *backwardix, double *backwardmi,
   double *backwardrg, double *backwardpk,
   double *crosstb, double *crossix, double *crossmi,
   double *crossrg, double *crosspk )
{
   const int k = threadIdx.x;
   const int deg1 = deg+1;
   int ix1,ix2;

   __shared__ double xvtb[pd_shmemsize];
   __shared__ double xvix[pd_shmemsize];
   __shared__ double xvmi[pd_shmemsize];
   __shared__ double xvrg[pd_shmemsize];
   __shared__ double xvpk[pd_shmemsize];
   __shared__ double yvtb[pd_shmemsize];
   __shared__ double yvix[pd_shmemsize];
   __shared__ double yvmi[pd_shmemsize];
   __shared__ double yvrg[pd_shmemsize];
   __shared__ double yvpk[pd_shmemsize];
   __shared__ double zvtb[pd_shmemsize];
   __shared__ double zvix[pd_shmemsize];
   __shared__ double zvmi[pd_shmemsize];
   __shared__ double zvrg[pd_shmemsize];
   __shared__ double zvpk[pd_shmemsize];
  
   xvtb[k] = cfftb[k]; xvix[k] = cffix[k]; xvmi[k] = cffmi[k]; 
   xvrg[k] = cffrg[k]; xvpk[k] = cffpk[k];
   ix1 = idx[0]*deg1+k;
   yvtb[k] = inputtb[ix1]; yvix[k] = inputix[ix1]; yvmi[k] = inputmi[ix1];
   yvrg[k] = inputrg[ix1]; yvpk[k] = inputpk[ix1]; 
   __syncthreads();
   dbl5_convolute(xvtb,xvix,xvmi,xvrg,xvpk,
                  yvtb,yvix,yvmi,yvrg,yvpk,
                  zvtb,zvix,zvmi,zvrg,zvpk,deg1,k);
   __syncthreads();
   forwardtb[k] = zvtb[k]; forwardix[k] = zvix[k];
   forwardmi[k] = zvmi[k]; forwardrg[k] = zvrg[k];
   forwardpk[k] = zvpk[k];                            // f[0] = cff*x[0]

   for(int i=1; i<nvr; i++)
   {
      xvtb[k] = zvtb[k]; xvix[k] = zvix[k]; xvmi[k] = zvmi[k];
      xvrg[k] = zvrg[k]; xvpk[k] = zvpk[k];
      ix2 = idx[i]*deg1+k;
      yvtb[k] = inputtb[ix2]; yvix[k] = inputix[ix2]; yvmi[k] = inputmi[ix2];
      yvrg[k] = inputrg[ix2]; yvpk[k] = inputpk[ix2];
      __syncthreads();
      dbl5_convolute(xvtb,xvix,xvmi,xvrg,xvpk,
                     yvtb,yvix,yvmi,yvrg,yvpk,
                     zvtb,zvix,zvmi,zvrg,zvpk,deg1,k);
      __syncthreads();
      ix1 = i*deg1+k;
      forwardtb[ix1] = zvtb[k]; forwardix[ix1] = zvix[k]; 
      forwardmi[ix1] = zvmi[k]; forwardrg[ix1] = zvrg[k]; 
      forwardpk[ix1] = zvpk[k];                        // f[i] = f[i-1]*x[i]
   }
   if(nvr > 2)
   {
      ix1 = idx[nvr-1]*deg1+k;
      xvtb[k] = inputtb[ix1]; xvix[k] = inputix[ix1]; xvmi[k] = inputmi[ix1];
      xvrg[k] = inputrg[ix1]; xvpk[k] = inputpk[ix1];
      ix2 = idx[nvr-2]*deg1+k;
      yvtb[k] = inputtb[ix2]; yvix[k] = inputix[ix2]; yvmi[k] = inputmi[ix2];
      yvrg[k] = inputrg[ix2]; yvpk[k] = inputpk[ix2];
      __syncthreads();
      dbl5_convolute(xvtb,xvix,xvmi,xvrg,xvpk,
                     yvtb,yvix,yvmi,yvrg,yvpk,
                     zvtb,zvix,zvmi,zvrg,zvpk,deg1,k);
      __syncthreads();
      backwardtb[k] = zvtb[k]; backwardix[k] = zvix[k];
      backwardmi[k] = zvmi[k]; backwardrg[k] = zvrg[k];
      backwardpk[k] = zvpk[k];                       // b[0] = x[n-1]*x[n-2]
      for(int i=1; i<nvr-2; i++)
      {
         xvtb[k] = zvtb[k]; xvix[k] = zvix[k]; xvmi[k] = zvmi[k]; 
         xvrg[k] = zvrg[k]; xvpk[k] = zvpk[k];
         ix2 = idx[nvr-2-i]*deg1+k;
         yvtb[k] = inputtb[ix2]; yvix[k] = inputix[ix2];
         yvmi[k] = inputmi[ix2]; yvrg[k] = inputrg[ix2];
         yvpk[k] = inputpk[ix2];
         __syncthreads();
         dbl5_convolute(xvtb,xvix,xvmi,xvrg,xvpk,
                        yvtb,yvix,yvmi,yvrg,yvpk,
                        zvtb,zvix,zvmi,zvrg,zvpk,deg1,k);
         __syncthreads();
         ix1 = i*deg1+k;
         backwardtb[ix1] = zvtb[k]; backwardix[ix1] = zvix[k];
         backwardmi[ix1] = zvmi[k]; backwardrg[ix1] = zvrg[k];
         backwardpk[ix1] = zvpk[k];             // b[i] = b[i-1]*x[n-2-i]
      }
      xvtb[k] = zvtb[k];  xvix[k] = zvix[k];  xvmi[k] = zvmi[k];
      xvrg[k] = zvrg[k];  xvpk[k] = zvpk[k];
      yvtb[k] = cfftb[k]; yvix[k] = cffix[k]; yvmi[k] = cffmi[k];
      yvrg[k] = cffrg[k]; yvpk[k] = cffpk[k];
      __syncthreads();
      dbl5_convolute(xvtb,xvix,xvmi,xvrg,xvpk,
                     yvtb,yvix,yvmi,yvrg,yvpk,
                     zvtb,zvix,zvmi,zvrg,zvpk,deg1,k);
      __syncthreads();
      ix2 = (nvr-3)*deg1+k;
      backwardtb[ix2] = zvtb[k]; backwardix[ix2] = zvix[k];
      backwardmi[ix2] = zvmi[k]; backwardrg[ix2] = zvrg[k];
      backwardpk[ix2] = zvpk[k];                    // b[n-3] = b[n-3]*cff

      if(nvr == 3)
      {
         xvtb[k] = forwardtb[k]; xvix[k] = forwardix[k];
         xvmi[k] = forwardmi[k]; xvrg[k] = forwardrg[k];
         xvpk[k] = forwardpk[k];
         ix2 = idx[2]*deg1+k;
         yvtb[k] = inputtb[ix2]; yvix[k] = inputix[ix2];
         yvmi[k] = inputmi[ix2]; yvrg[k] = inputrg[ix2]; 
         yvpk[k] = inputpk[ix2];
         __syncthreads();
         dbl5_convolute(xvtb,xvix,xvmi,xvrg,xvpk,
                        yvtb,yvix,yvmi,yvrg,yvpk,
                        zvtb,zvix,zvmi,zvrg,zvpk,deg1,k);
         __syncthreads();
         crosstb[k] = zvtb[k]; crossix[k] = zvix[k];
         crossmi[k] = zvmi[k]; crossrg[k] = zvrg[k];
         crosspk[k] = zvpk[k];                      // c[0] = f[0]*x[2]
      }
      else
      {
         for(int i=0; i<nvr-3; i++)
         {
            ix1 = i*deg1+k; 
            xvtb[k] = forwardtb[ix1]; xvix[k] = forwardix[ix1];
            xvmi[k] = forwardmi[ix1]; xvrg[k] = forwardrg[ix1];
            xvpk[k] = forwardpk[ix1];
            ix2 = (nvr-4-i)*deg1+k;
            yvtb[k] = backwardtb[ix2]; yvix[k] = backwardix[ix2];
            yvmi[k] = backwardmi[ix2]; yvrg[k] = backwardrg[ix2];
            yvpk[k] = backwardpk[ix2];
            __syncthreads();
            dbl5_convolute(xvtb,xvix,xvmi,xvrg,xvpk,
                           yvtb,yvix,yvmi,yvrg,yvpk,
                           zvtb,zvix,zvmi,zvrg,zvpk,deg1,k);
            __syncthreads();
            crosstb[ix1] = zvtb[k]; crossix[ix1] = zvix[k];
            crossmi[ix1] = zvmi[k]; crossrg[ix1] = zvrg[k];
            crosspk[ix1] = zvpk[k];            // c[i] = f[i]*b[n-4-i]
         }
         ix1 = (nvr-3)*deg1+k;
         xvtb[k] = forwardtb[ix1]; xvix[k] = forwardix[ix1];
         xvmi[k] = forwardmi[ix1]; xvrg[k] = forwardrg[ix1];
         xvpk[k] = forwardpk[ix1];
         ix2 = idx[nvr-1]*deg1+k;
         yvtb[k] = inputtb[ix2]; yvix[k] = inputix[ix2];
         yvmi[k] = inputmi[ix2]; yvrg[k] = inputrg[ix2];
         yvpk[k] = inputpk[ix2];
         __syncthreads();
         dbl5_convolute(xvtb,xvix,xvmi,xvrg,xvpk,
                        yvtb,yvix,yvmi,yvrg,yvpk,
                        zvtb,zvix,zvmi,zvrg,zvpk,deg1,k);
         __syncthreads();
         crosstb[ix1] = zvtb[k]; crossix[ix1] = zvix[k];
         crossmi[ix1] = zvmi[k]; crossrg[ix1] = zvrg[k];
         crosspk[ix1] = zvpk[k];                   // c[n-3] = f[n-3]*x[n-1]
      }
   }
}

__global__ void GPU_cmplx5_speel
 ( int nvr, int deg, int *idx,
   double *cffretb, double *cffreix, double *cffremi,
   double *cffrerg, double *cffrepk,
   double *cffimtb, double *cffimix, double *cffimmi,
   double *cffimrg, double *cffimpk,
   double *inputretb, double *inputreix, double *inputremi,
   double *inputrerg, double *inputrepk,
   double *inputimtb, double *inputimix, double *inputimmi,
   double *inputimrg, double *inputimpk,
   double *forwardretb, double *forwardreix, double *forwardremi,
   double *forwardrerg, double *forwardrepk,
   double *forwardimtb, double *forwardimix, double *forwardimmi,
   double *forwardimrg, double *forwardimpk,
   double *backwardretb, double *backwardreix, double *backwardremi,
   double *backwardrerg, double *backwardrepk,
   double *backwardimtb, double *backwardimix, double *backwardimmi,
   double *backwardimrg, double *backwardimpk,
   double *crossretb, double *crossreix, double *crossremi,
   double *crossrerg, double *crossrepk,
   double *crossimtb, double *crossimix, double *crossimmi,
   double *crossimrg, double *crossimpk )
{
   const int k = threadIdx.x;
   const int deg1 = deg+1;
   int ix1,ix2;

   __shared__ double xvretb[pd_shmemsize];
   __shared__ double xvreix[pd_shmemsize];
   __shared__ double xvremi[pd_shmemsize];
   __shared__ double xvrerg[pd_shmemsize];
   __shared__ double xvrepk[pd_shmemsize];
   __shared__ double xvimtb[pd_shmemsize];
   __shared__ double xvimix[pd_shmemsize];
   __shared__ double xvimmi[pd_shmemsize];
   __shared__ double xvimrg[pd_shmemsize];
   __shared__ double xvimpk[pd_shmemsize];
   __shared__ double yvretb[pd_shmemsize];
   __shared__ double yvreix[pd_shmemsize];
   __shared__ double yvremi[pd_shmemsize];
   __shared__ double yvrerg[pd_shmemsize];
   __shared__ double yvrepk[pd_shmemsize];
   __shared__ double yvimtb[pd_shmemsize];
   __shared__ double yvimix[pd_shmemsize];
   __shared__ double yvimmi[pd_shmemsize];
   __shared__ double yvimrg[pd_shmemsize];
   __shared__ double yvimpk[pd_shmemsize];
   __shared__ double zvretb[pd_shmemsize];
   __shared__ double zvreix[pd_shmemsize];
   __shared__ double zvremi[pd_shmemsize];
   __shared__ double zvrerg[pd_shmemsize];
   __shared__ double zvrepk[pd_shmemsize];
   __shared__ double zvimtb[pd_shmemsize];
   __shared__ double zvimix[pd_shmemsize];
   __shared__ double zvimmi[pd_shmemsize];
   __shared__ double zvimrg[pd_shmemsize];
   __shared__ double zvimpk[pd_shmemsize];

   xvretb[k] = cffretb[k]; xvreix[k] = cffreix[k]; xvremi[k] = cffremi[k];
   xvrerg[k] = cffrerg[k]; xvrepk[k] = cffrepk[k];
   xvimtb[k] = cffimtb[k]; xvimix[k] = cffimix[k]; xvimmi[k] = cffimmi[k];
   xvimrg[k] = cffimrg[k]; xvimpk[k] = cffimpk[k];
   ix1 = idx[0]*deg1+k;
   yvretb[k] = inputretb[ix1]; yvreix[k] = inputreix[ix1];
   yvremi[k] = inputremi[ix1]; yvrerg[k] = inputrerg[ix1];
   yvrepk[k] = inputrepk[ix1];
   yvimtb[k] = inputimtb[ix1]; yvimix[k] = inputimix[ix1];
   yvimmi[k] = inputimmi[ix1]; yvimrg[k] = inputimrg[ix1];
   yvimpk[k] = inputimpk[ix1];
   __syncthreads();                                      // f[0] = cff*x[0] 
   cmplx5_convolute(xvretb,xvreix,xvremi,xvrerg,xvrepk,
                    xvimtb,xvimix,xvimmi,xvimrg,xvimpk,
                    yvretb,yvreix,yvremi,yvrerg,yvrepk,
                    yvimtb,yvimix,yvimmi,yvimrg,yvimpk,
                    zvretb,zvreix,zvremi,zvrerg,zvrepk,
                    zvimtb,zvimix,zvimmi,zvimrg,zvimpk,deg1,k);
   __syncthreads();
   forwardretb[k] = zvretb[k]; forwardreix[k] = zvreix[k];
   forwardremi[k] = zvremi[k]; forwardrerg[k] = zvrerg[k];
   forwardrepk[k] = zvrepk[k];
   forwardimtb[k] = zvimtb[k]; forwardimix[k] = zvimix[k];
   forwardimmi[k] = zvimmi[k]; forwardimrg[k] = zvimrg[k];
   forwardimpk[k] = zvimpk[k];

   for(int i=1; i<nvr; i++)
   {
      xvretb[k] = zvretb[k]; xvreix[k] = zvreix[k];
      xvremi[k] = zvremi[k]; xvrerg[k] = zvrerg[k];
      xvrepk[k] = zvrepk[k];
      xvimtb[k] = zvimtb[k]; xvimix[k] = zvimix[k];
      xvimmi[k] = zvimmi[k]; xvimrg[k] = zvimrg[k];
      xvimpk[k] = zvimpk[k];
      ix2 = idx[i]*deg1+k;
      yvretb[k] = inputretb[ix2]; yvreix[k] = inputreix[ix2];
      yvremi[k] = inputremi[ix2]; yvrerg[k] = inputrerg[ix2];
      yvrepk[k] = inputrepk[ix2];
      yvimtb[k] = inputimtb[ix2]; yvimix[k] = inputimix[ix2];
      yvimmi[k] = inputimmi[ix2]; yvimrg[k] = inputimrg[ix2];
      yvimpk[k] = inputimpk[ix2];
      __syncthreads();                                 // f[i] = f[i-i]*x[i]
      cmplx5_convolute(xvretb,xvreix,xvremi,xvrerg,xvrepk,
                       xvimtb,xvimix,xvimmi,xvimrg,xvimpk,
                       yvretb,yvreix,yvremi,yvrerg,yvrepk,
                       yvimtb,yvimix,yvimmi,yvimrg,yvimpk,
                       zvretb,zvreix,zvremi,zvrerg,zvrepk,
                       zvimtb,zvimix,zvimmi,zvimrg,zvimpk,deg1,k);
      __syncthreads();
      ix1 = i*deg1+k;                                   
      forwardretb[ix1] = zvretb[k]; forwardreix[ix1] = zvreix[k];
      forwardremi[ix1] = zvremi[k]; forwardrerg[ix1] = zvrerg[k];
      forwardrepk[ix1] = zvrepk[k];
      forwardimtb[ix1] = zvimtb[k]; forwardimix[ix1] = zvimix[k];
      forwardimmi[ix1] = zvimmi[k]; forwardimrg[ix1] = zvimrg[k];
      forwardimpk[ix1] = zvimpk[k]; 
   }
   if(nvr > 2)
   {
      ix1 = idx[nvr-1]*deg1+k;
      xvretb[k] = inputretb[ix1]; xvreix[k] = inputreix[ix1];
      xvremi[k] = inputremi[ix1]; xvrerg[k] = inputrerg[ix1];
      xvrepk[k] = inputrepk[ix1];
      xvimtb[k] = inputimtb[ix1]; xvimix[k] = inputimix[ix1];
      xvimmi[k] = inputimmi[ix1]; xvimrg[k] = inputimrg[ix1];
      xvimpk[k] = inputimpk[ix1];
      ix2 = idx[nvr-2]*deg1+k;
      yvretb[k] = inputretb[ix2]; yvreix[k] = inputreix[ix2];
      yvremi[k] = inputremi[ix2]; yvrerg[k] = inputrerg[ix2];
      yvrepk[k] = inputrepk[ix2];
      yvimtb[k] = inputimtb[ix2]; yvimix[k] = inputimix[ix2];
      yvimmi[k] = inputimmi[ix2]; yvimrg[k] = inputimrg[ix2];
      yvimpk[k] = inputimpk[ix2];
      __syncthreads();                               // b[0] = x[n-1]*x[n-2]
      cmplx5_convolute(xvretb,xvreix,xvremi,xvrerg,xvrepk,
                       xvimtb,xvimix,xvimmi,xvimrg,xvimpk,
                       yvretb,yvreix,yvremi,yvrerg,yvrepk,
                       yvimtb,yvimix,yvimmi,yvimrg,yvimpk,
                       zvretb,zvreix,zvremi,zvrerg,zvrepk,
                       zvimtb,zvimix,zvimmi,zvimrg,zvimpk,deg1,k);
      __syncthreads();
      backwardretb[k] = zvretb[k]; backwardreix[k] = zvreix[k];
      backwardremi[k] = zvremi[k]; backwardrerg[k] = zvrerg[k];
      backwardrepk[k] = zvrepk[k];
      backwardimtb[k] = zvimtb[k]; backwardimix[k] = zvimix[k];
      backwardimmi[k] = zvimmi[k]; backwardimrg[k] = zvimrg[k];
      backwardimpk[k] = zvimpk[k];

      for(int i=1; i<nvr-2; i++)
      {
         xvretb[k] = zvretb[k]; xvreix[k] = zvreix[k];
         xvremi[k] = zvremi[k]; xvrerg[k] = zvrerg[k];
         xvrepk[k] = zvrepk[k];
         xvimtb[k] = zvimtb[k]; xvimix[k] = zvimix[k];
         xvimmi[k] = zvimmi[k]; xvimrg[k] = zvimrg[k];
         xvimpk[k] = zvimpk[k];
         ix2 = idx[nvr-2-i]*deg1+k;
         yvretb[k] = inputretb[ix2]; yvreix[k] = inputreix[ix2];
         yvremi[k] = inputremi[ix2]; yvrerg[k] = inputrerg[ix2];
         yvrepk[k] = inputrepk[ix2];
         yvimtb[k] = inputimtb[ix2]; yvimix[k] = inputimix[ix2];
         yvimmi[k] = inputimmi[ix2]; yvimrg[k] = inputimrg[ix2];
         yvimpk[k] = inputimpk[ix2];
         __syncthreads();                           // b[i] = b[i]*x[n-2-i]
         cmplx5_convolute(xvretb,xvreix,xvremi,xvrerg,xvrepk,
                          xvimtb,xvimix,xvimmi,xvimrg,xvimpk,
                          yvretb,yvreix,yvremi,yvrerg,yvrepk,
                          yvimtb,yvimix,yvimmi,yvimrg,yvimpk,
                          zvretb,zvreix,zvremi,zvrerg,zvrepk,
                          zvimtb,zvimix,zvimmi,zvimrg,zvimpk,deg1,k);
         __syncthreads();
         ix1 = i*deg1+k;
         backwardretb[ix1] = zvretb[k]; backwardreix[ix1] = zvreix[k];
         backwardremi[ix1] = zvremi[k]; backwardrerg[ix1] = zvrerg[k];
         backwardrepk[ix1] = zvrepk[k];
         backwardimtb[ix1] = zvimtb[k]; backwardimix[ix1] = zvimix[k];
         backwardimmi[ix1] = zvimmi[k]; backwardimrg[ix1] = zvimrg[k];
         backwardimpk[ix1] = zvimpk[k];
      }
      xvretb[k] = zvretb[k]; xvreix[k] = zvreix[k]; xvremi[k] = zvremi[k];
      xvrerg[k] = zvrerg[k]; xvrepk[k] = zvrepk[k];
      xvimtb[k] = zvimtb[k]; xvimix[k] = zvimix[k]; xvimmi[k] = zvimmi[k];
      xvimrg[k] = zvimrg[k]; xvimpk[k] = zvimpk[k];
      yvretb[k] = cffretb[k]; yvreix[k] = cffreix[k]; yvremi[k] = cffremi[k];
      yvrerg[k] = cffrerg[k]; yvrepk[k] = cffrepk[k];
      yvimtb[k] = cffimtb[k]; yvimix[k] = cffimix[k]; yvimmi[k] = cffimmi[k];
      yvimrg[k] = cffimrg[k]; yvimpk[k] = cffimpk[k];
      __syncthreads();                               // b[n-3] = b[n-3]*cff
      cmplx5_convolute(xvretb,xvreix,xvremi,xvrerg,xvrepk,
                       xvimtb,xvimix,xvimmi,xvimrg,xvimpk,
                       yvretb,yvreix,yvremi,yvrerg,yvrepk,
                       yvimtb,yvimix,yvimmi,yvimrg,yvimpk,
                       zvretb,zvreix,zvremi,zvrerg,zvrepk,
                       zvimtb,zvimix,zvimmi,zvimrg,zvimpk,deg1,k);
      __syncthreads();
      ix1 = (nvr-3)*deg1+k;
      backwardretb[ix1] = zvretb[k]; backwardreix[ix1] = zvreix[k];
      backwardremi[ix1] = zvremi[k]; backwardrerg[ix1] = zvrerg[k];
      backwardrepk[ix1] = zvrepk[k];
      backwardimtb[ix1] = zvimtb[k]; backwardimix[ix1] = zvimix[k];
      backwardimmi[ix1] = zvimmi[k]; backwardimrg[ix1] = zvimrg[k];
      backwardimpk[ix1] = zvimpk[k];

      if(nvr == 3)
      {
         xvretb[k] = forwardretb[k]; xvreix[k] = forwardreix[k];
         xvremi[k] = forwardremi[k]; xvrerg[k] = forwardrerg[k];
         xvrepk[k] = forwardrepk[k];
         xvimtb[k] = forwardimtb[k]; xvimix[k] = forwardimix[k];
         xvimmi[k] = forwardimmi[k]; xvimrg[k] = forwardimrg[k];
         xvimpk[k] = forwardimpk[k];
         ix2 = idx[2]*deg1+k;
         yvretb[k] = inputretb[ix2]; yvreix[k] = inputreix[ix2];
         yvremi[k] = inputremi[ix2]; yvrerg[k] = inputrerg[ix2];
         yvrepk[k] = inputrepk[ix2];
         yvimtb[k] = inputimtb[ix2]; yvimix[k] = inputimix[ix2];
         yvimmi[k] = inputimmi[ix2]; yvimrg[k] = inputimrg[ix2];
         yvimpk[k] = inputimpk[ix2];
         __syncthreads();                               // c[0] = f[0]*x[2]
         cmplx5_convolute(xvretb,xvreix,xvremi,xvrerg,xvrepk,
                          xvimtb,xvimix,xvimmi,xvimrg,xvimpk,
                          yvretb,yvreix,yvremi,yvrerg,yvrepk,
                          yvimtb,yvimix,yvimmi,yvimrg,yvimpk,
                          zvretb,zvreix,zvremi,zvrerg,zvrepk,
                          zvimtb,zvimix,zvimmi,zvimrg,zvimpk,deg1,k);
         __syncthreads();
         crossretb[k] = zvretb[k]; crossreix[k] = zvreix[k];
         crossremi[k] = zvremi[k]; crossrerg[k] = zvrerg[k];
         crossrepk[k] = zvrepk[k];
         crossimtb[k] = zvimtb[k]; crossimix[k] = zvimix[k];
         crossimmi[k] = zvimmi[k]; crossimrg[k] = zvimrg[k];
         crossimpk[k] = zvimpk[k];
      }
      else
      {
         for(int i=0; i<nvr-3; i++)
         {
            ix1 = i*deg1+k;   
            xvretb[k] = forwardretb[ix1]; xvreix[k] = forwardreix[ix1];
            xvremi[k] = forwardremi[ix1]; xvrerg[k] = forwardrerg[ix1];
            xvrepk[k] = forwardrepk[ix1];
            xvimtb[k] = forwardimtb[ix1]; xvimix[k] = forwardimix[ix1];
            xvimmi[k] = forwardimmi[ix1]; xvimrg[k] = forwardimrg[ix1];
            xvimpk[k] = forwardimpk[ix1];
            ix2 = (nvr-4-i)*deg1+k;
            yvretb[k] = backwardretb[ix2]; yvreix[k] = backwardreix[ix2];
            yvremi[k] = backwardremi[ix2]; yvrerg[k] = backwardrerg[ix2];
            yvrepk[k] = backwardrepk[ix2];
            yvimtb[k] = backwardimtb[ix2]; yvimix[k] = backwardimix[ix2];
            yvimmi[k] = backwardimmi[ix2]; yvimrg[k] = backwardimrg[ix2];
            yvimpk[k] = backwardimpk[ix2];
            __syncthreads();                        // c[i] = f[i]*b[n-4-i]
            cmplx5_convolute(xvretb,xvreix,xvremi,xvrerg,xvrepk,
                             xvimtb,xvimix,xvimmi,xvimrg,xvimpk,
                             yvretb,yvreix,yvremi,yvrerg,yvrepk,
                             yvimtb,yvimix,yvimmi,yvimrg,yvimpk,
                             zvretb,zvreix,zvremi,zvrerg,zvrepk,
                             zvimtb,zvimix,zvimmi,zvimrg,zvimpk,deg1,k);
            __syncthreads();
            ix1 = i*deg1+k;
            crossretb[ix1] = zvretb[k]; crossreix[ix1] = zvreix[k];
            crossremi[ix1] = zvremi[k]; crossrerg[ix1] = zvrerg[k];
            crossrepk[ix1] = zvrepk[k];
            crossimtb[ix1] = zvimtb[k]; crossimix[ix1] = zvimix[k];
            crossimmi[ix1] = zvimmi[k]; crossimrg[ix1] = zvimrg[k];
            crossimpk[ix1] = zvimpk[k];
         }
         ix1 = (nvr-3)*deg1+k;
         xvretb[k] = forwardretb[ix1]; xvreix[k] = forwardreix[ix1];
         xvremi[k] = forwardremi[ix1]; xvrerg[k] = forwardrerg[ix1];
         xvrepk[k] = forwardrepk[ix1];
         xvimtb[k] = forwardimtb[ix1]; xvimix[k] = forwardimix[ix1];
         xvimmi[k] = forwardimmi[ix1]; xvimrg[k] = forwardimrg[ix1];
         xvimpk[k] = forwardimpk[ix1];
         ix2 = idx[nvr-1]*deg1+k;
         yvretb[k] = inputretb[ix2]; yvreix[k] = inputreix[ix2];
         yvremi[k] = inputremi[ix2]; yvrerg[k] = inputrerg[ix2];
         yvrepk[k] = inputrepk[ix2];
         yvimtb[k] = inputimtb[ix2]; yvimix[k] = inputimix[ix2];
         yvimmi[k] = inputimmi[ix2]; yvimrg[k] = inputimrg[ix2];
         yvimpk[k] = inputimpk[ix2];
         __syncthreads();                         // c[n-3] = f[n-3]*x[n-1]
         cmplx5_convolute(xvretb,xvreix,xvremi,xvrerg,xvrepk,
                          xvimtb,xvimix,xvimmi,xvimrg,xvimpk,
                          yvretb,yvreix,yvremi,yvrerg,yvrepk,
                          yvimtb,yvimix,yvimmi,yvimrg,yvimpk,
                          zvretb,zvreix,zvremi,zvrerg,zvrepk,
                          zvimtb,zvimix,zvimmi,zvimrg,zvimpk,deg1,k);
         __syncthreads();
         ix1 = (nvr-3)*deg1+k;
         crossretb[ix1] = zvretb[k]; crossreix[ix1] = zvreix[k];
         crossremi[ix1] = zvremi[k]; crossrerg[ix1] = zvrerg[k];
         crossrepk[ix1] = zvrepk[k];
         crossimtb[ix1] = zvimtb[k]; crossimix[ix1] = zvimix[k];
         crossimmi[ix1] = zvimmi[k]; crossimrg[ix1] = zvimrg[k];
         crossimpk[ix1] = zvimpk[k];
      }
   }
}

void GPU_dbl5_evaldiff
 ( int BS, int dim, int nvr, int deg, int *idx,
   double *cfftb, double *cffix, double *cffmi, double *cffrg, double *cffpk,
   double **inputtb, double **inputix, double **inputmi,
   double **inputrg, double **inputpk,
   double **outputtb, double **outputix, double **outputmi,
   double **outputrg, double **outputpk )
{
   const int deg1 = deg+1;      // length of all vectors
   double *inputtb_d;           // inputtb_d is inputtb on the device
   double *inputix_d;           // inputix_d is inputix on the device
   double *inputmi_d;           // inputmi_d is inputmi on the device
   double *inputrg_d;           // inputrg_d is inputrg on the device
   double *inputpk_d;           // inputpk_d is inputpk on the device
   double *forwardtb_d;         // highest forward products on the device
   double *forwardix_d;         // second highest forward products
   double *forwardmi_d;         // middle forward products
   double *forwardrg_d;         // second lowest forward products
   double *forwardpk_d;         // lowest forward products
   double *backwardtb_d;        // highest backward products on the device
   double *backwardix_d;        // second highest backward products
   double *backwardmi_d;        // middle backward products
   double *backwardrg_d;        // second lowest backward products
   double *backwardpk_d;        // lowest backward products
   double *crosstb_d;           // highest cross products on the device
   double *crossix_d;           // second highest cross products
   double *crossmi_d;           // middle cross products
   double *crossrg_d;           // second lowest cross products
   double *crosspk_d;           // lowest cross products
   double *cfftb_d;             // cfftb_d is cfftb on device
   double *cffix_d;             // cffix_d is cffix on device
   double *cffmi_d;             // cffmi_d is cffmi on device
   double *cffrg_d;             // cffrg_d is cffrg on device
   double *cffpk_d;             // cffpk_d is cffpk on device
   int *idx_d;                  // idx_d is idx on device

   size_t szcff = deg1*sizeof(double);
   size_t szdim = dim*(deg1)*sizeof(double);
   size_t sznvr = nvr*(deg1)*sizeof(double);
   size_t sznvr2 = (nvr-2)*(deg1)*sizeof(double);
   size_t szidx = nvr*sizeof(int);

   cudaMalloc((void**)&idx_d,szidx);
   cudaMalloc((void**)&cfftb_d,szcff);
   cudaMalloc((void**)&cffix_d,szcff);
   cudaMalloc((void**)&cffmi_d,szcff);
   cudaMalloc((void**)&cffrg_d,szcff);
   cudaMalloc((void**)&cffpk_d,szcff);
   cudaMalloc((void**)&inputtb_d,szdim);
   cudaMalloc((void**)&inputix_d,szdim);
   cudaMalloc((void**)&inputmi_d,szdim);
   cudaMalloc((void**)&inputrg_d,szdim);
   cudaMalloc((void**)&inputpk_d,szdim);
   cudaMalloc((void**)&forwardtb_d,sznvr);
   cudaMalloc((void**)&forwardix_d,sznvr);
   cudaMalloc((void**)&forwardmi_d,sznvr);
   cudaMalloc((void**)&forwardrg_d,sznvr);
   cudaMalloc((void**)&forwardpk_d,sznvr);
   cudaMalloc((void**)&backwardtb_d,sznvr2);
   cudaMalloc((void**)&backwardix_d,sznvr2);
   cudaMalloc((void**)&backwardmi_d,sznvr2);
   cudaMalloc((void**)&backwardrg_d,sznvr2);
   cudaMalloc((void**)&backwardpk_d,sznvr2);
   cudaMalloc((void**)&crosstb_d,sznvr2);
   cudaMalloc((void**)&crossix_d,sznvr2);
   cudaMalloc((void**)&crossmi_d,sznvr2);
   cudaMalloc((void**)&crossrg_d,sznvr2);
   cudaMalloc((void**)&crosspk_d,sznvr2);

   double *inputtb_h = new double[dim*(deg1)];
   double *inputix_h = new double[dim*(deg1)];
   double *inputmi_h = new double[dim*(deg1)];
   double *inputrg_h = new double[dim*(deg1)];
   double *inputpk_h = new double[dim*(deg1)];
   int ix = 0;
   for(int i=0; i<dim; i++)
      for(int j=0; j<deg1; j++)
      {
         inputtb_h[ix] = inputtb[i][j];
         inputix_h[ix] = inputix[i][j];
         inputmi_h[ix] = inputmi[i][j];
         inputrg_h[ix] = inputrg[i][j];
         inputpk_h[ix++] = inputpk[i][j];
      }

   cudaMemcpy(idx_d,idx,szidx,cudaMemcpyHostToDevice);
   cudaMemcpy(cfftb_d,cfftb,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffix_d,cffix,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffmi_d,cffmi,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffrg_d,cffrg,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffpk_d,cffpk,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(inputtb_d,inputtb_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputix_d,inputix_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputmi_d,inputmi_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputrg_d,inputrg_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputpk_d,inputpk_h,szdim,cudaMemcpyHostToDevice);

   if(BS == deg1)
   {
      GPU_dbl5_speel<<<1,BS>>>
         (nvr,deg,idx_d,cfftb_d,cffix_d,cffmi_d,cffrg_d,cffpk_d,
          inputtb_d,inputix_d,inputmi_d,inputrg_d,inputpk_d,
          forwardtb_d,forwardix_d,forwardmi_d,forwardrg_d,forwardpk_d,
          backwardtb_d,backwardix_d,backwardmi_d,backwardrg_d,backwardpk_d,
          crosstb_d,crossix_d,crossmi_d,crossrg_d,crosspk_d);
   }
   double *forwardtb_h = new double[(deg1)*nvr];
   double *forwardix_h = new double[(deg1)*nvr];
   double *forwardmi_h = new double[(deg1)*nvr];
   double *forwardrg_h = new double[(deg1)*nvr];
   double *forwardpk_h = new double[(deg1)*nvr];
   double *backwardtb_h = new double[(deg1)*(nvr-2)];
   double *backwardix_h = new double[(deg1)*(nvr-2)];
   double *backwardmi_h = new double[(deg1)*(nvr-2)];
   double *backwardrg_h = new double[(deg1)*(nvr-2)];
   double *backwardpk_h = new double[(deg1)*(nvr-2)];
   double *crosstb_h = new double[(deg1)*(nvr-2)];
   double *crossix_h = new double[(deg1)*(nvr-2)];
   double *crossmi_h = new double[(deg1)*(nvr-2)];
   double *crossrg_h = new double[(deg1)*(nvr-2)];
   double *crosspk_h = new double[(deg1)*(nvr-2)];
  
   cudaMemcpy(forwardtb_h,forwardtb_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardix_h,forwardix_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardmi_h,forwardmi_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardrg_h,forwardrg_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardpk_h,forwardpk_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardtb_h,backwardtb_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardix_h,backwardix_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardmi_h,backwardmi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardrg_h,backwardrg_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardpk_h,backwardpk_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crosstb_h,crosstb_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossix_h,crossix_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossmi_h,crossmi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossrg_h,crossrg_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crosspk_h,crosspk_d,sznvr2,cudaMemcpyDeviceToHost);

   int offset = (nvr-1)*deg1;            // assign value of the monomial
   for(int i=0; i<deg1; i++)
   {
      outputtb[dim][i] = forwardtb_h[offset+i];
      outputix[dim][i] = forwardix_h[offset+i];
      outputmi[dim][i] = forwardmi_h[offset+i];
      outputrg[dim][i] = forwardrg_h[offset+i];
      outputpk[dim][i] = forwardpk_h[offset+i];
   }
   ix = idx[nvr-1];                      // derivative with respect to x[n-1]
   offset = (nvr-2)*deg1;
   for(int i=0; i<deg1; i++)
   {
      outputtb[ix][i] = forwardtb_h[offset+i];
      outputix[ix][i] = forwardix_h[offset+i];
      outputmi[ix][i] = forwardmi_h[offset+i];
      outputrg[ix][i] = forwardrg_h[offset+i];
      outputpk[ix][i] = forwardpk_h[offset+i];
   }
   ix = idx[0];                          // derivative with respect to x[0]
   offset = (nvr-3)*deg1;
   for(int i=0; i<deg1; i++)
   {
      outputtb[ix][i] = backwardtb_h[offset+i];
      outputix[ix][i] = backwardix_h[offset+i];
      outputmi[ix][i] = backwardmi_h[offset+i];
      outputrg[ix][i] = backwardrg_h[offset+i];
      outputpk[ix][i] = backwardpk_h[offset+i];
   }
   for(int k=1; k<nvr-1; k++)            // derivative with respect to x[k]
   {
      ix = idx[k]; offset = (k-1)*deg1;
      for(int i=0; i<deg1; i++)
      {
         outputtb[ix][i] = crosstb_h[offset+i];
         outputix[ix][i] = crossix_h[offset+i];
         outputmi[ix][i] = crossmi_h[offset+i];
         outputrg[ix][i] = crossrg_h[offset+i];
         outputpk[ix][i] = crosspk_h[offset+i];
      }
   }
}

void GPU_cmplx5_evaldiff
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
   const int deg1 = deg+1;          // length of all vectors
   double *inputretb_d;             // inputretb_d is inputretb on the device
   double *inputreix_d;             // inputreix_d is inputreix on the device
   double *inputremi_d;             // inputremi_d is inputremi on the device
   double *inputrerg_d;             // inputrerg_d is inputrerg on the device
   double *inputrepk_d;             // inputrepk_d is inputrepk on the device
   double *inputimtb_d;             // inputimtb_d is inputretb on the device
   double *inputimix_d;             // inputimix_d is inputreix on the device
   double *inputimmi_d;             // inputimmi_d is inputremi on the device
   double *inputimrg_d;             // inputimrg_d is inputrerg on the device
   double *inputimpk_d;             // inputimpk_d is inputrepk on the device
   double *forwardretb_d;
   double *forwardreix_d;
   double *forwardremi_d;
   double *forwardrerg_d;
   double *forwardrepk_d;
   double *forwardimtb_d;
   double *forwardimix_d;
   double *forwardimmi_d;
   double *forwardimrg_d;
   double *forwardimpk_d;
   double *backwardretb_d;
   double *backwardreix_d;
   double *backwardremi_d;
   double *backwardrerg_d;
   double *backwardrepk_d;
   double *backwardimtb_d;
   double *backwardimix_d;
   double *backwardimmi_d;
   double *backwardimrg_d;
   double *backwardimpk_d;
   double *crossretb_d;
   double *crossreix_d;
   double *crossremi_d;
   double *crossrerg_d;
   double *crossrepk_d;
   double *crossimtb_d;
   double *crossimix_d;
   double *crossimmi_d;
   double *crossimrg_d;
   double *crossimpk_d;
   double *cffretb_d;               // cffretb_d is cffretb on the device
   double *cffreix_d;               // cffreix_d is cffreix on the device
   double *cffremi_d;               // cffremi_d is cffremi on the device
   double *cffrerg_d;               // cffremi_d is cffrerg on the device
   double *cffrepk_d;               // cffrepk_d is cffrepk on the device
   double *cffimtb_d;               // cffimtb_d is cffimtb on the device
   double *cffimix_d;               // cffimix_d is cffimix on the device
   double *cffimmi_d;               // cffimmi_d is cffimmi on the device
   double *cffimrg_d;               // cffimrg_d is cffimrg on the device
   double *cffimpk_d;               // cffimpk_d is cffimpk on the device
   int *idx_d;                      // idx_d is idx on the device

   size_t szdim = dim*(deg1)*sizeof(double);
   size_t sznvr = nvr*(deg1)*sizeof(double);
   size_t sznvr2 = (nvr-2)*(deg1)*sizeof(double);
   size_t szidx = nvr*sizeof(int);
   size_t szcff = deg1*sizeof(double);

   cudaMalloc((void**)&idx_d,szidx);
   cudaMalloc((void**)&cffretb_d,szcff);
   cudaMalloc((void**)&cffreix_d,szcff);
   cudaMalloc((void**)&cffremi_d,szcff);
   cudaMalloc((void**)&cffrerg_d,szcff);
   cudaMalloc((void**)&cffrepk_d,szcff);
   cudaMalloc((void**)&cffimtb_d,szcff);
   cudaMalloc((void**)&cffimix_d,szcff);
   cudaMalloc((void**)&cffimmi_d,szcff);
   cudaMalloc((void**)&cffimrg_d,szcff);
   cudaMalloc((void**)&cffimpk_d,szcff);
   cudaMalloc((void**)&inputretb_d,szdim);
   cudaMalloc((void**)&inputreix_d,szdim);
   cudaMalloc((void**)&inputremi_d,szdim);
   cudaMalloc((void**)&inputrerg_d,szdim);
   cudaMalloc((void**)&inputrepk_d,szdim);
   cudaMalloc((void**)&inputimtb_d,szdim);
   cudaMalloc((void**)&inputimix_d,szdim);
   cudaMalloc((void**)&inputimmi_d,szdim);
   cudaMalloc((void**)&inputimrg_d,szdim);
   cudaMalloc((void**)&inputimpk_d,szdim);
   cudaMalloc((void**)&forwardretb_d,sznvr);
   cudaMalloc((void**)&forwardreix_d,sznvr);
   cudaMalloc((void**)&forwardremi_d,sznvr);
   cudaMalloc((void**)&forwardrerg_d,sznvr);
   cudaMalloc((void**)&forwardrepk_d,sznvr);
   cudaMalloc((void**)&forwardimtb_d,sznvr);
   cudaMalloc((void**)&forwardimix_d,sznvr);
   cudaMalloc((void**)&forwardimmi_d,sznvr);
   cudaMalloc((void**)&forwardimrg_d,sznvr);
   cudaMalloc((void**)&forwardimpk_d,sznvr);
   cudaMalloc((void**)&backwardretb_d,sznvr2);
   cudaMalloc((void**)&backwardreix_d,sznvr2);
   cudaMalloc((void**)&backwardremi_d,sznvr2);
   cudaMalloc((void**)&backwardrerg_d,sznvr2);
   cudaMalloc((void**)&backwardrepk_d,sznvr2);
   cudaMalloc((void**)&backwardimtb_d,sznvr2);
   cudaMalloc((void**)&backwardimix_d,sznvr2);
   cudaMalloc((void**)&backwardimmi_d,sznvr2);
   cudaMalloc((void**)&backwardimrg_d,sznvr2);
   cudaMalloc((void**)&backwardimpk_d,sznvr2);
   cudaMalloc((void**)&crossretb_d,sznvr2);
   cudaMalloc((void**)&crossreix_d,sznvr2);
   cudaMalloc((void**)&crossremi_d,sznvr2);
   cudaMalloc((void**)&crossrerg_d,sznvr2);
   cudaMalloc((void**)&crossrepk_d,sznvr2);
   cudaMalloc((void**)&crossimtb_d,sznvr2);
   cudaMalloc((void**)&crossimix_d,sznvr2);
   cudaMalloc((void**)&crossimmi_d,sznvr2);
   cudaMalloc((void**)&crossimrg_d,sznvr2);
   cudaMalloc((void**)&crossimpk_d,sznvr2);

   double *inputretb_h = new double[dim*(deg1)];
   double *inputreix_h = new double[dim*(deg1)];
   double *inputremi_h = new double[dim*(deg1)];
   double *inputrerg_h = new double[dim*(deg1)];
   double *inputrepk_h = new double[dim*(deg1)];
   double *inputimtb_h = new double[dim*(deg1)];
   double *inputimix_h = new double[dim*(deg1)];
   double *inputimmi_h = new double[dim*(deg1)];
   double *inputimrg_h = new double[dim*(deg1)];
   double *inputimpk_h = new double[dim*(deg1)];
   int ix = 0;
   for(int i=0; i<dim; i++)
      for(int j=0; j<deg1; j++)
      {
         inputretb_h[ix] = inputretb[i][j];
         inputreix_h[ix] = inputreix[i][j];
         inputremi_h[ix] = inputremi[i][j];
         inputrerg_h[ix] = inputrerg[i][j];
         inputrepk_h[ix] = inputrepk[i][j];
         inputimtb_h[ix] = inputimtb[i][j];
         inputimix_h[ix] = inputimix[i][j];
         inputimmi_h[ix] = inputimmi[i][j];
         inputimrg_h[ix] = inputimrg[i][j];
         inputimpk_h[ix++] = inputimpk[i][j];
      }

   cudaMemcpy(idx_d,idx,szidx,cudaMemcpyHostToDevice);
   cudaMemcpy(cffretb_d,cffretb,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffreix_d,cffreix,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffremi_d,cffremi,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffrerg_d,cffrerg,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffrepk_d,cffrepk,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffimtb_d,cffimtb,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffimix_d,cffimix,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffimmi_d,cffimmi,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffimrg_d,cffimrg,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(cffimpk_d,cffimpk,szcff,cudaMemcpyHostToDevice);
   cudaMemcpy(inputretb_d,inputretb_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputreix_d,inputreix_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputremi_d,inputremi_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputrerg_d,inputrerg_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputrepk_d,inputrepk_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputimtb_d,inputimtb_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputimix_d,inputimix_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputimmi_d,inputimmi_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputimrg_d,inputimrg_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(inputimpk_d,inputimpk_h,szdim,cudaMemcpyHostToDevice);

   if(BS == deg1)
   {
      GPU_cmplx5_speel<<<1,BS>>>(nvr,deg,idx_d,
         cffretb_d,cffreix_d,cffremi_d,cffrerg_d,cffrepk_d,
         cffimtb_d,cffimix_d,cffimmi_d,cffimrg_d,cffimpk_d,
         inputretb_d,inputreix_d,inputremi_d,inputrerg_d,inputrepk_d,
         inputimtb_d,inputimix_d,inputimmi_d,inputimrg_d,inputimpk_d,
         forwardretb_d,forwardreix_d,forwardremi_d,
         forwardrerg_d,forwardrepk_d,
         forwardimtb_d,forwardimix_d,forwardimmi_d,
         forwardimrg_d,forwardimpk_d,
         backwardretb_d,backwardreix_d,backwardremi_d,
         backwardrerg_d,backwardrepk_d,
         backwardimtb_d,backwardimix_d,backwardimmi_d,
         backwardimrg_d,backwardimpk_d,
         crossretb_d,crossreix_d,crossremi_d,crossrerg_d,crossrepk_d,
         crossimtb_d,crossimix_d,crossimmi_d,crossimrg_d,crossimpk_d);
   }
   double *forwardretb_h = new double[(deg1)*nvr];
   double *forwardreix_h = new double[(deg1)*nvr];
   double *forwardremi_h = new double[(deg1)*nvr];
   double *forwardrerg_h = new double[(deg1)*nvr];
   double *forwardrepk_h = new double[(deg1)*nvr];
   double *forwardimtb_h = new double[(deg1)*nvr];
   double *forwardimix_h = new double[(deg1)*nvr];
   double *forwardimmi_h = new double[(deg1)*nvr];
   double *forwardimrg_h = new double[(deg1)*nvr];
   double *forwardimpk_h = new double[(deg1)*nvr];
   double *backwardretb_h = new double[(deg1)*(nvr-2)];
   double *backwardreix_h = new double[(deg1)*(nvr-2)];
   double *backwardremi_h = new double[(deg1)*(nvr-2)];
   double *backwardrerg_h = new double[(deg1)*(nvr-2)];
   double *backwardrepk_h = new double[(deg1)*(nvr-2)];
   double *backwardimtb_h = new double[(deg1)*(nvr-2)];
   double *backwardimix_h = new double[(deg1)*(nvr-2)];
   double *backwardimmi_h = new double[(deg1)*(nvr-2)];
   double *backwardimrg_h = new double[(deg1)*(nvr-2)];
   double *backwardimpk_h = new double[(deg1)*(nvr-2)];
   double *crossretb_h = new double[(deg1)*(nvr-2)];
   double *crossreix_h = new double[(deg1)*(nvr-2)];
   double *crossremi_h = new double[(deg1)*(nvr-2)];
   double *crossrerg_h = new double[(deg1)*(nvr-2)];
   double *crossrepk_h = new double[(deg1)*(nvr-2)];
   double *crossimtb_h = new double[(deg1)*(nvr-2)];
   double *crossimix_h = new double[(deg1)*(nvr-2)];
   double *crossimmi_h = new double[(deg1)*(nvr-2)];
   double *crossimrg_h = new double[(deg1)*(nvr-2)];
   double *crossimpk_h = new double[(deg1)*(nvr-2)];
  
   cudaMemcpy(forwardretb_h,forwardretb_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardreix_h,forwardreix_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardremi_h,forwardremi_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardrerg_h,forwardrerg_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardrepk_h,forwardrepk_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardimtb_h,forwardimtb_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardimix_h,forwardimix_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardimmi_h,forwardimmi_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardimrg_h,forwardimrg_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(forwardimpk_h,forwardimpk_d,sznvr,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardretb_h,backwardretb_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardreix_h,backwardreix_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardremi_h,backwardremi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardrerg_h,backwardrerg_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardrepk_h,backwardrepk_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardimtb_h,backwardimtb_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardimix_h,backwardimix_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardimmi_h,backwardimmi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardimrg_h,backwardimrg_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(backwardimpk_h,backwardimpk_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossretb_h,crossretb_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossreix_h,crossreix_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossremi_h,crossremi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossrerg_h,crossrerg_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossrepk_h,crossrepk_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossimtb_h,crossimtb_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossimix_h,crossimix_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossimmi_h,crossimmi_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossimrg_h,crossimrg_d,sznvr2,cudaMemcpyDeviceToHost);
   cudaMemcpy(crossimpk_h,crossimpk_d,sznvr2,cudaMemcpyDeviceToHost);

   int offset = (nvr-1)*deg1;
   for(int i=0; i<deg1; i++)   // assign value of the monomial
   {
      outputretb[dim][i] = forwardretb_h[offset+i];
      outputreix[dim][i] = forwardreix_h[offset+i];
      outputremi[dim][i] = forwardremi_h[offset+i];
      outputrerg[dim][i] = forwardrerg_h[offset+i];
      outputrepk[dim][i] = forwardrepk_h[offset+i];
      outputimtb[dim][i] = forwardimtb_h[offset+i];
      outputimix[dim][i] = forwardimix_h[offset+i];
      outputimmi[dim][i] = forwardimmi_h[offset+i];
      outputimrg[dim][i] = forwardimrg_h[offset+i];
      outputimpk[dim][i] = forwardimpk_h[offset+i];
   }
   ix = idx[nvr-1];
   offset = (nvr-2)*deg1;
   for(int i=0; i<deg1; i++)  // derivative with respect to x[n-1]
   {
      outputretb[ix][i] = forwardretb_h[offset+i];
      outputreix[ix][i] = forwardreix_h[offset+i];
      outputremi[ix][i] = forwardremi_h[offset+i];
      outputrerg[ix][i] = forwardrerg_h[offset+i];
      outputrepk[ix][i] = forwardrepk_h[offset+i];
      outputimtb[ix][i] = forwardimtb_h[offset+i];
      outputimix[ix][i] = forwardimix_h[offset+i];
      outputimmi[ix][i] = forwardimmi_h[offset+i];
      outputimrg[ix][i] = forwardimrg_h[offset+i];
      outputimpk[ix][i] = forwardimpk_h[offset+i];
   }
   ix = idx[0]; 
   offset = (nvr-3)*deg1;
   for(int i=0; i<deg1; i++)   // derivative with respect to x[0]
   {
      outputretb[ix][i] = backwardretb_h[offset+i];
      outputreix[ix][i] = backwardreix_h[offset+i];
      outputremi[ix][i] = backwardremi_h[offset+i];
      outputrerg[ix][i] = backwardrerg_h[offset+i];
      outputrepk[ix][i] = backwardrepk_h[offset+i];
      outputimtb[ix][i] = backwardimtb_h[offset+i];
      outputimix[ix][i] = backwardimix_h[offset+i];
      outputimmi[ix][i] = backwardimmi_h[offset+i];
      outputimrg[ix][i] = backwardimrg_h[offset+i];
      outputimpk[ix][i] = backwardimpk_h[offset+i];
   }
   for(int k=1; k<nvr-1; k++)  // derivative with respect to x[k]
   {
      ix = idx[k]; offset = (k-1)*deg1;
      for(int i=0; i<deg1; i++)
      {
         outputretb[ix][i] = crossretb_h[offset+i];
         outputreix[ix][i] = crossreix_h[offset+i];
         outputremi[ix][i] = crossremi_h[offset+i];
         outputrerg[ix][i] = crossrerg_h[offset+i];
         outputrepk[ix][i] = crossrepk_h[offset+i];
         outputimtb[ix][i] = crossimtb_h[offset+i];
         outputimix[ix][i] = crossimix_h[offset+i];
         outputimmi[ix][i] = crossimmi_h[offset+i];
         outputimrg[ix][i] = crossimrg_h[offset+i];
         outputimpk[ix][i] = crossimpk_h[offset+i];
      }
   }
}
