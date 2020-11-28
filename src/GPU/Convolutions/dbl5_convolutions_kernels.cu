// The file dbl5_convolutions_kernels.cu defines kernels with prototypes
// in dbl5_convolutions_kernels.h.

#ifdef gpufun
#include "double_double_gpufun.cu"
#include "penta_double_gpufun.cu"
#endif
#include "dbl5_convolutions_kernels.h"

__global__ void dbl5_increment
 ( double *xtb, double *xix, double *xmi, double *xrg, double *xpk,
   double *ytb, double *yix, double *ymi, double *yrg, double *ypk,
   double *ztb, double *zix, double *zmi, double *zrg, double *zpk, int dim )
{
   int k = threadIdx.x;                 // thread k computes z[k]

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

   xvtb[k] = xtb[k]; xvix[k] = xix[k]; xvmi[k] = xmi[k];
   xvrg[k] = xrg[k]; xvpk[k] = xpk[k];
   yvtb[k] = ytb[k]; yvix[k] = yix[k]; yvmi[k] = ymi[k];
   yvrg[k] = yrg[k]; yvpk[k] = ypk[k];

   pdg_add(xvtb[k],xvix[k],xvmi[k],xvrg[k],xvpk[k],
           yvtb[k],yvix[k],yvmi[k],yvrg[k],yvpk[k],
           &zvtb[k],&zvix[k],&zvmi[k],&zvrg[k],&zvpk[k]);

   // __syncthreads();

   ztb[k] = zvtb[k]; zix[k] = zvix[k]; zmi[k] = zvmi[k];
   zrg[k] = zvrg[k]; zpk[k] = zvpk[k];
}

__global__ void dbl5_decrement
 ( double *xtb, double *xix, double *xmi, double *xrg, double *xpk,
   double *ytb, double *yix, double *ymi, double *yrg, double *ypk,
   double *ztb, double *zix, double *zmi, double *zrg, double *zpk, int dim )
{
   int k = threadIdx.x;                 // thread k computes z[k]

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

   xvtb[k] = xtb[k]; xvix[k] = xix[k]; xvmi[k] = xmi[k];
   xvrg[k] = xrg[k]; xvpk[k] = xpk[k];
   yvtb[k] = ytb[k]; yvix[k] = yix[k]; yvmi[k] = ymi[k];
   yvrg[k] = yrg[k]; yvpk[k] = ypk[k];

   pdg_sub(xvtb[k],xvix[k],xvmi[k],xvrg[k],xvpk[k],
           yvtb[k],yvix[k],yvmi[k],yvrg[k],yvpk[k],
           &zvtb[k],&zvix[k],&zvmi[k],&zvrg[k],&zvpk[k]);

   // __syncthreads();

   ztb[k] = zvtb[k]; zix[k] = zvix[k]; zmi[k] = zvmi[k];
   zrg[k] = zvrg[k]; zpk[k] = zvpk[k];
}

__global__ void dbl5_convolute
 ( double *xtb, double *xix, double *xmi, double *xrg, double *xpk,
   double *ytb, double *yix, double *ymi, double *yrg, double *ypk,
   double *ztb, double *zix, double *zmi, double *zrg, double *zpk, int dim )
{
   int k = threadIdx.x;                 // thread k computes z[k]

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

   double prdtb,prdix,prdmi,prdrg,prdpk;
   int idx;

   xvtb[k] = xtb[k]; xvix[k] = xix[k]; xvmi[k] = xmi[k];
   xvrg[k] = xrg[k]; xvpk[k] = xpk[k];
   yvtb[k] = ytb[k]; yvix[k] = yix[k]; yvmi[k] = ymi[k];
   yvrg[k] = yrg[k]; yvpk[k] = ypk[k];

   __syncthreads();

   // zv[k] = xv[0]*yv[k];
   pdg_mul( xvtb[0], xvix[0], xvmi[0], xvrg[0], xvpk[0],
            yvtb[k], yvix[k], yvmi[k], yvrg[k], yvpk[k],
           &zvtb[k],&zvix[k],&zvmi[k],&zvrg[k],&zvpk[k]);

   // __syncthreads();

   for(int i=1; i<=k; i++) // zv[k] = zv[k] + xv[i]*yv[k-i];
   {
      idx = k-i;
      pdg_mul(  xvtb[i],  xvix[i],  xvmi[i],  xvrg[i],  xvpk[i],
                yvtb[idx],yvix[idx],yvmi[idx],yvrg[idx],yvpk[idx],
              &prdtb,   &prdix,   &prdmi,   &prdrg,   &prdpk);
      pdg_inc(&zvtb[k],&zvix[k],&zvmi[k],&zvrg[k],&zvpk[k],
              prdtb,   prdix,   prdmi,   prdrg,   prdpk);
   }

   // __syncthreads();

   ztb[k] = zvtb[k]; zix[k] = zvix[k]; zmi[k] = zvmi[k];
   zrg[k] = zvrg[k]; zpk[k] = zvpk[k];
}

__global__ void dbl5_padded_convolute
 ( double *xtb, double *xix, double *xmi, double *xrg, double *xpk,
   double *ytb, double *yix, double *ymi, double *yrg, double *ypk,
   double *ztb, double *zix, double *zmi, double *zrg, double *zpk, int dim )
{
   int k = threadIdx.x;                 // thread k computes z[k]

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

   double prdtb,prdix,prdmi,prdrg,prdpk;
   int idx = dim+k;

   xvtb[k] = xtb[k]; xvix[k] = xix[k]; xvmi[k] = xmi[k];
   xvrg[k] = xrg[k]; xvpk[k] = xpk[k];
   yvtb[k] = 0.0; yvix[k] = 0.0; yvmi[k] = 0.0;
   yvrg[k] = 0.0; yvpk[k] = 0.0;
   yvtb[idx] = ytb[k]; yvix[idx] = yix[k]; yvmi[idx] = ymi[k];
   yvrg[idx] = yrg[k]; yvpk[idx] = ypk[k];

   // zv[k] = xv[0]*yv[k];
   pdg_mul( xvtb[0], xvix[0], xvmi[0], xvrg[0], xvpk[0],
            yvtb[idx], yvix[idx], yvmi[idx], yvrg[idx], yvpk[idx],
           &zvtb[k],&zvix[k],&zvmi[k],&zvrg[k],&zvpk[k]);

   for(int i=1; i<dim; i++) // zv[k] = zv[k] + xv[i]*yv[k-i];
   {
      idx = dim+k-i;
      pdg_mul(  xvtb[i],  xvix[i],  xvmi[i],  xvrg[i],  xvpk[i],
                yvtb[idx],yvix[idx],yvmi[idx],yvrg[idx],yvpk[idx],
              &prdtb,   &prdix,   &prdmi,   &prdrg,   &prdpk);
      pdg_inc(&zvtb[k],&zvix[k],&zvmi[k],&zvrg[k],&zvpk[k],
              prdtb,   prdix,   prdmi,   prdrg,   prdpk);
   }

   ztb[k] = zvtb[k]; zix[k] = zvix[k]; zmi[k] = zvmi[k];
   zrg[k] = zvrg[k]; zpk[k] = zvpk[k];
}

__global__ void cmplx5_convolute
 ( double *xretb, double *xreix, double *xremi, double *xrerg, double *xrepk,
   double *ximtb, double *ximix, double *ximmi, double *ximrg, double *ximpk,
   double *yretb, double *yreix, double *yremi, double *yrerg, double *yrepk,
   double *yimtb, double *yimix, double *yimmi, double *yimrg, double *yimpk,
   double *zretb, double *zreix, double *zremi, double *zrerg, double *zrepk,
   double *zimtb, double *zimix, double *zimmi, double *zimrg, double *zimpk,
   int dim )
{
   int k = threadIdx.x;       // thread k computes zre[k] and zim[k]

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

   double xrtb,xitb,yrtb,yitb,zrtb,zitb,acctb;
   double xrix,xiix,yrix,yiix,zrix,ziix,accix;
   double xrmi,ximi,yrmi,yimi,zrmi,zimi,accmi;
   double xrrg,xirg,yrrg,yirg,zrrg,zirg,accrg;
   double xrpk,xipk,yrpk,yipk,zrpk,zipk,accpk;

   xvretb[k] = xretb[k]; xvimtb[k] = ximtb[k];
   xvreix[k] = xreix[k]; xvimix[k] = ximix[k];
   xvremi[k] = xremi[k]; xvimmi[k] = ximmi[k];
   xvrerg[k] = xrerg[k]; xvimrg[k] = ximrg[k];
   xvrepk[k] = xrepk[k]; xvimpk[k] = ximpk[k];
   yvretb[k] = yretb[k]; yvimtb[k] = yimtb[k];
   yvreix[k] = yreix[k]; yvimix[k] = yimix[k];
   yvremi[k] = yremi[k]; yvimmi[k] = yimmi[k];
   yvrerg[k] = yrerg[k]; yvimrg[k] = yimrg[k];
   yvrepk[k] = yrepk[k]; yvimpk[k] = yimpk[k];

   __syncthreads();

   // z[k] = x[0]*y[k]
   xrtb = xvretb[0]; xrix = xvreix[0]; xrmi = xvremi[0];
   xrrg = xvrerg[0]; xrpk = xvrepk[0];
   xitb = xvimtb[0]; xiix = xvimix[0]; ximi = xvimmi[0];
   xirg = xvimrg[0]; xipk = xvimpk[0];
   yrtb = yvretb[k]; yrix = yvreix[k]; yrmi = yvremi[k];
   yrrg = yvrerg[k]; yrpk = yvrepk[k];
   yitb = yvimtb[k]; yiix = yvimix[k]; yimi = yvimmi[k];
   yirg = yvimrg[k]; yipk = yvimpk[k];

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

   zvretb[k] = zrtb; zvreix[k] = zrix; zvremi[k] = zrmi;
   zvrerg[k] = zrrg; zvrepk[k] = zrpk;
   zvimtb[k] = zitb; zvimix[k] = ziix; zvimmi[k] = zimi;
   zvimrg[k] = zirg; zvimpk[k] = zipk;

   for(int i=1; i<=k; i++) // z[k] = z[k] + x[i]*y[k-i]
   {
      xrtb = xvretb[i]; xrix = xvreix[i]; xrmi = xvremi[i];
      xrrg = xvrerg[i]; xrpk = xvrepk[i];
      xitb = xvimtb[i]; xiix = xvimix[i]; ximi = xvimmi[i];
      xirg = xvimrg[i]; xipk = xvimpk[i];
      yrtb = yvretb[k-i]; yrix = yvreix[k-i]; yrmi = yvremi[k-i];
      yrrg = yvrerg[k-i]; yrpk = yvrepk[k-i];
      yitb = yvimtb[k-i]; yiix = yvimix[k-i]; yimi = yvimmi[k-i];
      yirg = yvimrg[k-i]; yipk = yvimpk[k-i];

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
      // zvre[k] += zr; zvim[k] += zi
      pdg_inc(&zvretb[k],&zvreix[k],&zvremi[k],&zvrerg[k],&zvrepk[k],
              zrtb,zrix,zrmi,zrrg,zrpk);
      pdg_inc(&zvimtb[k],&zvimix[k],&zvimmi[k],&zvimrg[k],&zvimpk[k],
              zitb,ziix,zimi,zirg,zipk);
   }

   __syncthreads();

   zretb[k] = zvretb[k]; zreix[k] = zvreix[k]; zremi[k] = zvremi[k];
   zrerg[k] = zvrerg[k]; zrepk[k] = zvrepk[k];
   zimtb[k] = zvimtb[k]; zimix[k] = zvimix[k]; zimmi[k] = zvimmi[k];
   zimrg[k] = zvimrg[k]; zimpk[k] = zvimpk[k];
}

__global__ void cmplx5_padded_convolute
 ( double *xretb, double *xreix, double *xremi, double *xrerg, double *xrepk,
   double *ximtb, double *ximix, double *ximmi, double *ximrg, double *ximpk,
   double *yretb, double *yreix, double *yremi, double *yrerg, double *yrepk,
   double *yimtb, double *yimix, double *yimmi, double *yimrg, double *yimpk,
   double *zretb, double *zreix, double *zremi, double *zrerg, double *zrepk,
   double *zimtb, double *zimix, double *zimmi, double *zimrg, double *zimpk,
   int dim )
{
   int k = threadIdx.x;       // thread k computes zre[k] and zim[k]

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

   double xrtb,xitb,yrtb,yitb,zrtb,zitb,acctb;
   double xrix,xiix,yrix,yiix,zrix,ziix,accix;
   double xrmi,ximi,yrmi,yimi,zrmi,zimi,accmi;
   double xrrg,xirg,yrrg,yirg,zrrg,zirg,accrg;
   double xrpk,xipk,yrpk,yipk,zrpk,zipk,accpk;
   int idx = dim+k;

   xvretb[k] = xretb[k]; xvimtb[k] = ximtb[k];
   xvreix[k] = xreix[k]; xvimix[k] = ximix[k];
   xvremi[k] = xremi[k]; xvimmi[k] = ximmi[k];
   xvrerg[k] = xrerg[k]; xvimrg[k] = ximrg[k];
   xvrepk[k] = xrepk[k]; xvimpk[k] = ximpk[k];
   yvretb[k] = 0.0; yvimtb[k] = 0.0;
   yvreix[k] = 0.0; yvimix[k] = 0.0;
   yvremi[k] = 0.0; yvimmi[k] = 0.0;
   yvrerg[k] = 0.0; yvimrg[k] = 0.0;
   yvrepk[k] = 0.0; yvimpk[k] = 0.0;
   yvretb[idx] = yretb[k]; yvimtb[idx] = yimtb[k];
   yvreix[idx] = yreix[k]; yvimix[idx] = yimix[k];
   yvremi[idx] = yremi[k]; yvimmi[idx] = yimmi[k];
   yvrerg[idx] = yrerg[k]; yvimrg[idx] = yimrg[k];
   yvrepk[idx] = yrepk[k]; yvimpk[idx] = yimpk[k];

   __syncthreads();

   // z[k] = x[0]*y[k]
   xrtb = xvretb[0]; xrix = xvreix[0]; xrmi = xvremi[0];
   xrrg = xvrerg[0]; xrpk = xvrepk[0];
   xitb = xvimtb[0]; xiix = xvimix[0]; ximi = xvimmi[0];
   xirg = xvimrg[0]; xipk = xvimpk[0];
   yrtb = yvretb[idx]; yrix = yvreix[idx]; yrmi = yvremi[idx];
   yrrg = yvrerg[idx]; yrpk = yvrepk[idx];
   yitb = yvimtb[idx]; yiix = yvimix[idx]; yimi = yvimmi[idx];
   yirg = yvimrg[idx]; yipk = yvimpk[idx];

   pdg_mul(xrtb,xrix,xrmi,xrrg,xrpk,yrtb,yrix,yrmi,yrrg,yrpk,
           &zrtb,&zrix,&zrmi,&zrrg,&zrpk);         // zr = xr*yr
   __syncthreads();
   pdg_mul(xitb,xiix,ximi,xirg,xipk,yitb,yiix,yimi,yirg,yipk,
           &acctb,&accix,&accmi,&accrg,&accpk);    // acc = xi*yi
   __syncthreads();
   pdg_minus(&acctb,&accix,&accmi,&accrg,&accpk);
   pdg_inc(&zrtb,&zrix,&zrmi,&zrrg,&zrpk,
           acctb,accix,accmi,accrg,accpk);         // zr = xr*yr - xi*yi
   __syncthreads();
   pdg_mul(xrtb,xrix,xrmi,xrrg,xrpk,yitb,yiix,yimi,yirg,yipk,
           &zitb,&ziix,&zimi,&zirg,&zipk);         // zi = xr*yi
   __syncthreads();
   pdg_mul(xitb,xiix,ximi,xirg,xipk,yrtb,yrix,yrmi,yrrg,yrpk,
           &acctb,&accix,&accmi,&accrg,&accpk);    // acc = xi*yr
   __syncthreads();
   pdg_inc(&zitb,&ziix,&zimi,&zirg,&zipk,
           acctb,accix,accmi,accrg,accpk);         // zr = xr*yr + xi*yi
   __syncthreads();

   zvretb[k] = zrtb; zvreix[k] = zrix; zvremi[k] = zrmi;
   zvrerg[k] = zrrg; zvrepk[k] = zrpk;
   zvimtb[k] = zitb; zvimix[k] = ziix; zvimmi[k] = zimi;
   zvimrg[k] = zirg; zvimpk[k] = zipk;

   for(int i=1; i<dim; i++) // z[k] = z[k] + x[i]*y[k-i]
   {
      idx = dim + k - i;
      xrtb = xvretb[i]; xrix = xvreix[i]; xrmi = xvremi[i];
      xrrg = xvrerg[i]; xrpk = xvrepk[i];
      xitb = xvimtb[i]; xiix = xvimix[i]; ximi = xvimmi[i];
      xirg = xvimrg[i]; xipk = xvimpk[i];
      yrtb = yvretb[idx]; yrix = yvreix[idx]; yrmi = yvremi[idx];
      yrrg = yvrerg[idx]; yrpk = yvrepk[idx];
      yitb = yvimtb[idx]; yiix = yvimix[idx]; yimi = yvimmi[idx];
      yirg = yvimrg[idx]; yipk = yvimpk[idx];

      pdg_mul(xrtb,xrix,xrmi,xrrg,xrpk,yrtb,yrix,yrmi,yrrg,yrpk,
              &zrtb,&zrix,&zrmi,&zrrg,&zrpk);        // zr = xr*yr
      __syncthreads();
      pdg_mul(xitb,xiix,ximi,xirg,xipk,yitb,yiix,yimi,yirg,yipk,
              &acctb,&accix,&accmi,&accrg,&accpk);   // xi*yi
      __syncthreads();
      pdg_minus(&acctb,&accix,&accmi,&accrg,&accpk);
      pdg_inc(&zrtb,&zrix,&zrmi,&zrrg,&zrpk,
              acctb,accix,accmi,accrg,accpk);        // zr = xr*yr - xi*yi
      __syncthreads();
      pdg_mul(xrtb,xrix,xrmi,xrrg,xrpk,yitb,yiix,yimi,yirg,yipk,
              &zitb,&ziix,&zimi,&zirg,&zipk);        // zi = xr*yi
      __syncthreads();
      pdg_mul(xitb,xiix,ximi,xirg,xipk,yrtb,yrix,yrmi,yrrg,yrpk,
              &acctb,&accix,&accmi,&accrg,&accpk);   // xi*yr
      __syncthreads();
      pdg_inc(&zitb,&ziix,&zimi,&zirg,&zipk,
              acctb,accix,accmi,accrg,accpk);        // zr = xr*yr + xi*yi
      __syncthreads();
      // zvre[k] += zr; zvim[k] += zi
      pdg_inc(&zvretb[k],&zvreix[k],&zvremi[k],&zvrerg[k],&zvrepk[k],
              zrtb,zrix,zrmi,zrrg,zrpk);
      __syncthreads();
      pdg_inc(&zvimtb[k],&zvimix[k],&zvimmi[k],&zvimrg[k],&zvimpk[k],
              zitb,ziix,zimi,zirg,zipk);
      __syncthreads();
   }
   zretb[k] = zvretb[k]; zreix[k] = zvreix[k]; zremi[k] = zvremi[k];
   zrerg[k] = zvrerg[k]; zrepk[k] = zvrepk[k];
   zimtb[k] = zvimtb[k]; zimix[k] = zvimix[k]; zimmi[k] = zvimmi[k];
   zimrg[k] = zvimrg[k]; zimpk[k] = zvimpk[k];
}

void GPU_dbl5_product
 ( double *xtb_h, double *xix_h, double *xmi_h, double *xrg_h, double *xpk_h,
   double *ytb_h, double *yix_h, double *ymi_h, double *yrg_h, double *ypk_h,
   double *ztb_h, double *zix_h, double *zmi_h, double *zrg_h, double *zpk_h,
   int deg, int freq, int BS, int padded )
{
   const int dim = deg+1;            // length of all vectors
   double* xtb_d;                    // xtb_d is xtb_h on the device
   double* xix_d;                    // xix_d is xix_h on the device
   double* xmi_d;                    // xmi_d is xmi_h on the device
   double* xrg_d;                    // xrg_d is xrg_h on the device
   double* xpk_d;                    // xpk_d is xpk_h on the device
   double* ytb_d;                    // ytb_d is ytb_h on the device
   double* yix_d;                    // yix_d is yix_h on the device
   double* ymi_d;                    // ymi_d is ymi_h on the device
   double* yrg_d;                    // yrg_d is yrg_h on the device
   double* ypk_d;                    // ypk_d is ypk_h on the device
   double* ztb_d;                    // ztb_d is ztb_h on the device
   double* zix_d;                    // zix_d is zix_h on the device
   double* zmi_d;                    // zmi_d is zmi_h on the device
   double* zrg_d;                    // zrg_d is zrg_h on the device
   double* zpk_d;                    // zpk_d is zpk_h on the device
   size_t size = dim*sizeof(double); // number of bytes for each vector

   cudaMalloc((void**)&xtb_d,size);
   cudaMalloc((void**)&xix_d,size);
   cudaMalloc((void**)&xmi_d,size);
   cudaMalloc((void**)&xrg_d,size);
   cudaMalloc((void**)&xpk_d,size);
   cudaMalloc((void**)&ytb_d,size);
   cudaMalloc((void**)&yix_d,size);
   cudaMalloc((void**)&ymi_d,size);
   cudaMalloc((void**)&yrg_d,size);
   cudaMalloc((void**)&ypk_d,size);
   cudaMalloc((void**)&ztb_d,size);
   cudaMalloc((void**)&zix_d,size);
   cudaMalloc((void**)&zmi_d,size);
   cudaMalloc((void**)&zrg_d,size);
   cudaMalloc((void**)&zpk_d,size);
   cudaMemcpy(xtb_d,xtb_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(xix_d,xix_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(xmi_d,xmi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(xrg_d,xrg_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(xpk_d,xpk_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(ytb_d,ytb_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yix_d,yix_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(ymi_d,ymi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yrg_d,yrg_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(ypk_d,ypk_h,size,cudaMemcpyHostToDevice);

   if(dim == BS)
   {
      if(padded == 1)
      {
         for(int i=0; i<freq; i++)
            dbl5_padded_convolute<<<1,BS>>>
               (xtb_d,xix_d,xmi_d,xrg_d,xpk_d,
                ytb_d,yix_d,ymi_d,yrg_d,ypk_d,
                ztb_d,zix_d,zmi_d,zrg_d,zpk_d,dim);
      }
      else
      {
         for(int i=0; i<freq; i++)
            dbl5_convolute<<<1,BS>>>
               (xtb_d,xix_d,xmi_d,xrg_d,xpk_d,
                ytb_d,yix_d,ymi_d,yrg_d,ypk_d,
                ztb_d,zix_d,zmi_d,zrg_d,zpk_d,dim);
      }
   }
   cudaMemcpy(ztb_h,ztb_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zix_h,zix_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zmi_h,zmi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zrg_h,zrg_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zpk_h,zpk_d,size,cudaMemcpyDeviceToHost);
}

void GPU_cmplx5_product
 ( double *xretb_h, double *xreix_h, double *xremi_h, double *xrerg_h,
   double *xrepk_h, double *ximtb_h, double *ximix_h, double *ximmi_h,
   double *ximrg_h, double *ximpk_h, double *yretb_h, double *yreix_h,
   double *yremi_h, double *yrerg_h, double *yrepk_h, double *yimtb_h,
   double *yimix_h, double *yimmi_h, double *yimrg_h, double *yimpk_h,
   double *zretb_h, double *zreix_h, double *zremi_h, double *zrerg_h,
   double *zrepk_h, double *zimtb_h, double *zimix_h, double *zimmi_h,
   double *zimrg_h, double *zimpk_h, int deg, int freq, int BS,
   int mode )
{
   const int dim = deg+1;            // length of all vectors
   double* xretb_d;                  // xretb_d is xretb_h on the device
   double* xreix_d;                  // xreix_d is xreix_h on the device
   double* xremi_d;                  // xremi_d is xremi_h on the device
   double* xrerg_d;                  // xrerg_d is xrerg_h on the device
   double* xrepk_d;                  // xrepk_d is xrepk_h on the device
   double* ximtb_d;                  // ximtb_d is ximtb_h on the device
   double* ximix_d;                  // ximix_d is ximix_h on the device
   double* ximmi_d;                  // ximmi_d is ximmi_h on the device
   double* ximrg_d;                  // ximrg_d is ximrg_h on the device
   double* ximpk_d;                  // ximpk_d is ximpk_h on the device
   double* yretb_d;                  // yretb_d is yretb_h on the device
   double* yreix_d;                  // yreix_d is yreix_h on the device
   double* yremi_d;                  // yremi_d is yremi_h on the device
   double* yrerg_d;                  // yrerg_d is yrerg_h on the device
   double* yrepk_d;                  // yrepk_d is yrepk_h on the device
   double* yimtb_d;                  // yimtb_d is yimtb_h on the device
   double* yimix_d;                  // yimix_d is yimix_h on the device
   double* yimmi_d;                  // yimmi_d is yimmi_h on the device
   double* yimrg_d;                  // yimrg_d is yimrg_h on the device
   double* yimpk_d;                  // yimpk_d is yimpk_h on the device
   double* zretb_d;                  // zretb_d is zretb_h on the device
   double* zreix_d;                  // zreix_d is zreix_h on the device
   double* zremi_d;                  // zremi_d is zremi_h on the device
   double* zrerg_d;                  // zrerg_d is zrerg_h on the device
   double* zrepk_d;                  // zrepk_d is zrepk_h on the device
   double* zimtb_d;                  // zimtb_d is zimtb_h on the device
   double* zimix_d;                  // zimix_d is zimix_h on the device
   double* zimmi_d;                  // zimmi_d is zimmi_h on the device
   double* zimrg_d;                  // zimrg_d is zimrg_h on the device
   double* zimpk_d;                  // zimpk_d is zimpk_h on the device
   double* acctb_d;                  // accumulates highest doubles
   double* accix_d;                  // accumulates second highest doubles
   double* accmi_d;                  // accumulates middle doubles
   double* accrg_d;                  // accumulates second lowest doubles
   double* accpk_d;                  // accumulates lowest doubles
   size_t size = dim*sizeof(double); // number of bytes for each vector

   cudaMalloc((void**)&xretb_d,size);
   cudaMalloc((void**)&xreix_d,size);
   cudaMalloc((void**)&xremi_d,size);
   cudaMalloc((void**)&xrerg_d,size);
   cudaMalloc((void**)&xrepk_d,size);
   cudaMalloc((void**)&ximtb_d,size);
   cudaMalloc((void**)&ximix_d,size);
   cudaMalloc((void**)&ximmi_d,size);
   cudaMalloc((void**)&ximrg_d,size);
   cudaMalloc((void**)&ximpk_d,size);
   cudaMalloc((void**)&yretb_d,size);
   cudaMalloc((void**)&yreix_d,size);
   cudaMalloc((void**)&yremi_d,size);
   cudaMalloc((void**)&yrerg_d,size);
   cudaMalloc((void**)&yrepk_d,size);
   cudaMalloc((void**)&yimtb_d,size);
   cudaMalloc((void**)&yimix_d,size);
   cudaMalloc((void**)&yimmi_d,size);
   cudaMalloc((void**)&yimrg_d,size);
   cudaMalloc((void**)&yimpk_d,size);
   cudaMalloc((void**)&zretb_d,size);
   cudaMalloc((void**)&zreix_d,size);
   cudaMalloc((void**)&zremi_d,size);
   cudaMalloc((void**)&zrerg_d,size);
   cudaMalloc((void**)&zrepk_d,size);
   cudaMalloc((void**)&zimtb_d,size);
   cudaMalloc((void**)&zimix_d,size);
   cudaMalloc((void**)&zimmi_d,size);
   cudaMalloc((void**)&zimrg_d,size);
   cudaMalloc((void**)&zimpk_d,size);
   cudaMalloc((void**)&acctb_d,size);
   cudaMalloc((void**)&accix_d,size);
   cudaMalloc((void**)&accmi_d,size);
   cudaMalloc((void**)&accrg_d,size);
   cudaMalloc((void**)&accpk_d,size);
   cudaMemcpy(xretb_d,xretb_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(xreix_d,xreix_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(xremi_d,xremi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(xrerg_d,xrerg_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(xrepk_d,xrepk_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(ximtb_d,ximtb_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(ximix_d,ximix_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(ximmi_d,ximmi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(ximrg_d,ximrg_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(ximpk_d,ximpk_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yretb_d,yretb_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yreix_d,yreix_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yremi_d,yremi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yrerg_d,yrerg_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yrepk_d,yrepk_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yimtb_d,yimtb_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yimix_d,yimix_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yimmi_d,yimmi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yimrg_d,yimrg_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yimpk_d,yimpk_h,size,cudaMemcpyHostToDevice);

   if(dim == BS)
   {
      if(mode == 2)
      {
         for(int i=0; i<freq; i++)
            cmplx5_padded_convolute<<<1,BS>>>
               (xretb_d,xreix_d,xremi_d,xrerg_d,xrepk_d,
                ximtb_d,ximix_d,ximmi_d,ximrg_d,ximpk_d,
                yretb_d,yreix_d,yremi_d,yrerg_d,yrepk_d,
                yimtb_d,yimix_d,yimmi_d,yimrg_d,yimpk_d,
                zretb_d,zreix_d,zremi_d,zrerg_d,zrepk_d,
                zimtb_d,zimix_d,zimmi_d,zimrg_d,zimpk_d,dim);
      }
      else if(mode == 1)
      {
         for(int i=0; i<freq; i++)
         {
            dbl5_padded_convolute<<<1,BS>>>
               (xretb_d,xreix_d,xremi_d,xrerg_d,xrepk_d,
                yretb_d,yreix_d,yremi_d,yrerg_d,yrepk_d,
                zretb_d,zreix_d,zremi_d,zrerg_d,zrepk_d,dim);
            dbl5_padded_convolute<<<1,BS>>>
               (ximtb_d,ximix_d,ximmi_d,ximrg_d,ximpk_d,
                yimtb_d,yimix_d,yimmi_d,yimrg_d,yimpk_d,
                acctb_d,accix_d,accmi_d,accrg_d,zimpk_d,dim);
            dbl5_decrement<<<1,BS>>>
               (zretb_d,zreix_d,zremi_d,zrerg_d,zrepk_d,
                acctb_d,accix_d,accmi_d,accrg_d,accpk_d,
                zretb_d,zreix_d,zremi_d,zrerg_d,zrepk_d,dim);
            dbl5_padded_convolute<<<1,BS>>>
               (xretb_d,xreix_d,xremi_d,xrerg_d,xrepk_d,
                yimtb_d,yimix_d,yimmi_d,yimrg_d,yimpk_d,
                zimtb_d,zimix_d,zimmi_d,zimrg_d,zimpk_d,dim);
            dbl5_padded_convolute<<<1,BS>>>
               (ximtb_d,ximix_d,ximmi_d,ximrg_d,ximpk_d,
                yretb_d,yreix_d,yremi_d,yrerg_d,yrepk_d,
                acctb_d,accix_d,accmi_d,accrg_d,accpk_d,dim);
            dbl5_increment<<<1,BS>>>
               (zimtb_d,zimix_d,zimmi_d,zimrg_d,zimpk_d,
                acctb_d,accix_d,accmi_d,accrg_d,accpk_d,
                zimtb_d,zimix_d,zimmi_d,zimrg_d,zimpk_d,dim);
         }
      }
      else
      {
         for(int i=0; i<freq; i++)
            cmplx5_convolute<<<1,BS>>>
               (xretb_d,xreix_d,xremi_d,xrerg_d,xrepk_d,
                ximtb_d,ximix_d,ximmi_d,ximrg_d,ximpk_d,
                yretb_d,yreix_d,yremi_d,yrerg_d,yrepk_d,
                yimtb_d,yimix_d,yimmi_d,yimrg_d,yimpk_d,
                zretb_d,zreix_d,zremi_d,zrerg_d,zrepk_d,
                zimtb_d,zimix_d,zimmi_d,zimrg_d,zimpk_d,dim);
      }
   }

   cudaMemcpy(zretb_h,zretb_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zreix_h,zreix_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zremi_h,zremi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zrerg_h,zrerg_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zrepk_h,zrepk_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zimtb_h,zimtb_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zimix_h,zimix_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zimmi_h,zimmi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zimrg_h,zimrg_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zimpk_h,zimpk_d,size,cudaMemcpyDeviceToHost);
}
