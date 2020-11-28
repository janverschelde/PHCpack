// The file dbl3_convolutions_kernels.cu defines kernels with prototypes
// in dbl3_convolutions_kernels.h.

#ifdef gpufun
#include "double_double_gpufun.cu"
#include "triple_double_gpufun.cu"
#endif
#include "dbl3_convolutions_kernels.h"

__global__ void dbl3_increment
 ( double *xhi, double *xmi, double *xlo,
   double *yhi, double *ymi, double *ylo,
   double *zhi, double *zmi, double *zlo, int dim )
{
   int k = threadIdx.x;                 // thread k computes z[k]

   __shared__ double xvhi[td_shmemsize];
   __shared__ double xvmi[td_shmemsize];
   __shared__ double xvlo[td_shmemsize];
   __shared__ double yvhi[td_shmemsize];
   __shared__ double yvmi[td_shmemsize];
   __shared__ double yvlo[td_shmemsize];
   __shared__ double zvhi[td_shmemsize];
   __shared__ double zvmi[td_shmemsize];
   __shared__ double zvlo[td_shmemsize];

   xvhi[k] = xhi[k]; xvmi[k] = xmi[k]; xvlo[k] = xlo[k];
   yvhi[k] = yhi[k]; yvmi[k] = ymi[k]; yvlo[k] = ylo[k];

   tdg_add(xvhi[k],xvmi[k],xvlo[k],yvhi[k],yvmi[k],yvlo[k],
           &zvhi[k],&zvmi[k],&zvlo[k]);

   __syncthreads();

   zhi[k] = zvhi[k]; zmi[k] = zvmi[k]; zlo[k] = zvlo[k];
}

__global__ void dbl3_decrement
 ( double *xhi, double *xmi, double *xlo,
   double *yhi, double *ymi, double *ylo,
   double *zhi, double *zmi, double *zlo, int dim )
{
   int k = threadIdx.x;                 // thread k computes z[k]

   __shared__ double xvhi[td_shmemsize];
   __shared__ double xvmi[td_shmemsize];
   __shared__ double xvlo[td_shmemsize];
   __shared__ double yvhi[td_shmemsize];
   __shared__ double yvmi[td_shmemsize];
   __shared__ double yvlo[td_shmemsize];
   __shared__ double zvhi[td_shmemsize];
   __shared__ double zvmi[td_shmemsize];
   __shared__ double zvlo[td_shmemsize];

   xvhi[k] = xhi[k]; xvmi[k] = xmi[k]; xvlo[k] = xlo[k];
   yvhi[k] = yhi[k]; yvmi[k] = ymi[k]; yvlo[k] = ylo[k];

   tdg_sub(xvhi[k],xvmi[k],xvlo[k],yvhi[k],yvmi[k],yvlo[k],
           &zvhi[k],&zvmi[k],&zvlo[k]);

   __syncthreads();

   zhi[k] = zvhi[k]; zmi[k] = zvmi[k]; zlo[k] = zvlo[k];
}

__global__ void dbl3_convolute
 ( double *xhi, double *xmi, double *xlo,
   double *yhi, double *ymi, double *ylo,
   double *zhi, double *zmi, double *zlo, int dim )
{
   int k = threadIdx.x;                 // thread k computes z[k]

   __shared__ double xvhi[td_shmemsize];
   __shared__ double xvmi[td_shmemsize];
   __shared__ double xvlo[td_shmemsize];
   __shared__ double yvhi[td_shmemsize];
   __shared__ double yvmi[td_shmemsize];
   __shared__ double yvlo[td_shmemsize];
   __shared__ double zvhi[td_shmemsize];
   __shared__ double zvmi[td_shmemsize];
   __shared__ double zvlo[td_shmemsize];

   double prdhi,prdmi,prdlo;

   xvhi[k] = xhi[k]; xvmi[k] = xmi[k]; xvlo[k] = xlo[k];
   yvhi[k] = yhi[k]; yvmi[k] = ymi[k]; yvlo[k] = ylo[k];
   __syncthreads();

   // zv[k] = xv[0]*yv[k];
   tdg_mul(xvhi[0],xvmi[0],xvlo[0],yvhi[k],yvmi[k],yvlo[k],
           &zvhi[k],&zvmi[k],&zvlo[k]);
   __syncthreads();

   for(int i=1; i<=k; i++) // zv[k] = zv[k] + xv[i]*yv[k-i];
   {
      tdg_mul(xvhi[i],xvmi[i],xvlo[i],yvhi[k-i],yvmi[k-i],yvlo[k-i],
              &prdhi,&prdmi,&prdlo);
      __syncthreads();
      tdg_inc(&zvhi[k],&zvmi[k],&zvlo[k],prdhi,prdmi,prdlo);
      __syncthreads();
   }
   zhi[k] = zvhi[k]; zmi[k] = zvmi[k]; zlo[k] = zvlo[k];
}

__global__ void dbl3_padded_convolute
 ( double *xhi, double *xmi, double *xlo,
   double *yhi, double *ymi, double *ylo,
   double *zhi, double *zmi, double *zlo, int dim )
{
   int k = threadIdx.x;                 // thread k computes z[k]

   __shared__ double xvhi[td_shmemsize];
   __shared__ double xvmi[td_shmemsize];
   __shared__ double xvlo[td_shmemsize];
   __shared__ double yvhi[td_shmemsize];
   __shared__ double yvmi[td_shmemsize];
   __shared__ double yvlo[td_shmemsize];
   __shared__ double zvhi[td_shmemsize];
   __shared__ double zvmi[td_shmemsize];
   __shared__ double zvlo[td_shmemsize];

   double prdhi,prdmi,prdlo;
   int idx = dim+k;

   xvhi[k] = xhi[k];   xvmi[k] = xmi[k];   xvlo[k] = xlo[k];
   yvhi[k] = 0.0;      yvmi[k] = 0.0;      yvlo[k] = 0.0;
   yvhi[idx] = yhi[k]; yvmi[idx] = ymi[k]; yvlo[idx] = ylo[k];

   __syncthreads();

   // zv[k] = xv[0]*yv[k];
   tdg_mul(xvhi[0],xvmi[0],xvlo[0],yvhi[idx],yvmi[idx],yvlo[idx],
           &zvhi[k],&zvmi[k],&zvlo[k]);
   __syncthreads();

   for(int i=1; i<dim; i++) // zv[k] = zv[k] + xv[i]*yv[k-i];
   {
      int idx = dim + k - i;
      tdg_mul(xvhi[i],xvmi[i],xvlo[i],yvhi[idx],yvmi[idx],yvlo[idx],
              &prdhi,&prdmi,&prdlo);
      __syncthreads();
      tdg_inc(&zvhi[k],&zvmi[k],&zvlo[k],prdhi,prdmi,prdlo);
      __syncthreads();
   }

   zhi[k] = zvhi[k]; zmi[k] = zvmi[k]; zlo[k] = zvlo[k];
}

__global__ void cmplx3_convolute
 ( double *xrehi, double *xremi, double *xrelo,
   double *ximhi, double *ximmi, double *ximlo,
   double *yrehi, double *yremi, double *yrelo,
   double *yimhi, double *yimmi, double *yimlo,
   double *zrehi, double *zremi, double *zrelo,
   double *zimhi, double *zimmi, double *zimlo, int dim )
{
   int k = threadIdx.x;       // thread k computes zre[k] and zim[k]

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

   double xrhi,xihi,yrhi,yihi,zrhi,zihi,acchi;
   double xrmi,ximi,yrmi,yimi,zrmi,zimi,accmi;
   double xrlo,xilo,yrlo,yilo,zrlo,zilo,acclo;

   xvrehi[k] = xrehi[k]; xvimhi[k] = ximhi[k];
   xvremi[k] = xremi[k]; xvimmi[k] = ximmi[k];
   xvrelo[k] = xrelo[k]; xvimlo[k] = ximlo[k];
   yvrehi[k] = yrehi[k]; yvimhi[k] = yimhi[k];
   yvremi[k] = yremi[k]; yvimmi[k] = yimmi[k];
   yvrelo[k] = yrelo[k]; yvimlo[k] = yimlo[k];

   __syncthreads();

   // z[k] = x[0]*y[k]
   xrhi = xvrehi[0]; xrmi = xvremi[0]; xrlo = xvrelo[0];
   xihi = xvimhi[0]; ximi = xvimmi[0]; xilo = xvimlo[0];
   yrhi = yvrehi[k]; yrmi = yvremi[k]; yrlo = yvrelo[k];
   yihi = yvimhi[k]; yimi = yvimmi[k]; yilo = yvimlo[k];

   tdg_mul(xrhi,xrmi,xrlo,yrhi,yrmi,yrlo,&zrhi,&zrmi,&zrlo);     // zr = xr*yr
   tdg_mul(xihi,ximi,xilo,yihi,yimi,yilo,&acchi,&accmi,&acclo); // acc = xi*yi
   tdg_minus(&acchi,&accmi,&acclo);
   tdg_inc(&zrhi,&zrmi,&zrlo,acchi,accmi,acclo);         // zr = xr*yr - xi*yi
   tdg_mul(xrhi,xrmi,xrlo,yihi,yimi,yilo,&zihi,&zimi,&zilo);     // zi = xr*yi
   tdg_mul(xihi,ximi,xilo,yrhi,yrmi,yrlo,&acchi,&accmi,&acclo); // acc = xi*yr
   tdg_inc(&zihi,&zimi,&zilo,acchi,accmi,acclo);         // zr = xr*yr + xi*yi

   zvrehi[k] = zrhi; zvremi[k] = zrmi; zvrelo[k] = zrlo;
   zvimhi[k] = zihi; zvimmi[k] = zimi; zvimlo[k] = zilo;

   for(int i=1; i<=k; i++) // z[k] = z[k] + x[i]*y[k-i]
   {
      xrhi = xvrehi[i]; xrmi = xvremi[i]; xrlo = xvrelo[i];
      xihi = xvimhi[i]; ximi = xvimmi[i]; xilo = xvimlo[i];
      yrhi = yvrehi[k-i]; yrmi = yvremi[k-i]; yrlo = yvrelo[k-i];
      yihi = yvimhi[k-i]; yimi = yvimmi[k-i]; yilo = yvimlo[k-i];

      tdg_mul(xrhi,xrmi,xrlo,yrhi,yrmi,yrlo,&zrhi,&zrmi,&zrlo); // zr = xr*yr
      tdg_mul(xihi,ximi,xilo,yihi,yimi,yilo,&acchi,&accmi,&acclo);   // xi*yi
      tdg_minus(&acchi,&accmi,&acclo);
      tdg_inc(&zrhi,&zrmi,&zrlo,acchi,accmi,acclo);     // zr = xr*yr - xi*yi
      tdg_mul(xrhi,xrmi,xrlo,yihi,yimi,yilo,&zihi,&zimi,&zilo); // zi = xr*yi
      tdg_mul(xihi,ximi,xilo,yrhi,yrmi,yrlo,&acchi,&accmi,&acclo);   // xi*yr
      tdg_inc(&zihi,&zimi,&zilo,acchi,accmi,acclo);     // zr = xr*yr + xi*yi
      // zvre[k] += zr; zvim[k] += zi
      tdg_inc(&zvrehi[k],&zvremi[k],&zvrelo[k],zrhi,zrmi,zrlo);
      tdg_inc(&zvimhi[k],&zvimmi[k],&zvimlo[k],zihi,zimi,zilo);
   }

   __syncthreads();

   zrehi[k] = zvrehi[k]; zremi[k] = zvremi[k]; zrelo[k] = zvrelo[k];
   zimhi[k] = zvimhi[k]; zimmi[k] = zvimmi[k]; zimlo[k] = zvimlo[k];
}

__global__ void cmplx3_looped_convolute
 ( double *xrehi, double *xremi, double *xrelo,
   double *ximhi, double *ximmi, double *ximlo,
   double *yrehi, double *yremi, double *yrelo,
   double *yimhi, double *yimmi, double *yimlo,
   double *zrehi, double *zremi, double *zrelo,
   double *zimhi, double *zimmi, double *zimlo, int dim )
{
   int k = threadIdx.x;       // thread k computes zre[k] and zim[dim-1-k]

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

   const int dim1 = dim-1;
   int idx;
   double xrhi,xihi,yrhi,yihi,zvhi,acchi;
   double xrmi,ximi,yrmi,yimi,zvmi,accmi;
   double xrlo,xilo,yrlo,yilo,zvlo,acclo;

   xvrehi[k] = xrehi[k]; xvimhi[k] = ximhi[k];
   xvremi[k] = xremi[k]; xvimmi[k] = ximmi[k];
   xvrelo[k] = xrelo[k]; xvimlo[k] = ximlo[k];
   yvrehi[k] = yrehi[k]; yvimhi[k] = yimhi[k];
   yvremi[k] = yremi[k]; yvimmi[k] = yimmi[k];
   yvrelo[k] = yrelo[k]; yvimlo[k] = yimlo[k];

   __syncthreads();

   // z[k] = x[0]*y[k]
   xrhi = xvrehi[0]; xrmi = xvremi[0]; xrlo = xvrelo[0];
   xihi = xvimhi[0]; ximi = xvimmi[0]; xilo = xvimlo[0];
   yrhi = yvrehi[k]; yrmi = yvremi[k]; yrlo = yvrelo[k];
   yihi = yvimhi[k]; yimi = yvimmi[k]; yilo = yvimlo[k];

   tdg_mul(xrhi,xrmi,xrlo,yrhi,yrmi,yrlo,&zvhi,&zvmi,&zvlo);     // zr = xr*yr
   tdg_mul(xihi,ximi,xilo,yihi,yimi,yilo,&acchi,&accmi,&acclo); // acc = xi*yi
   tdg_minus(&acchi,&accmi,&acclo);
   tdg_inc(&zvhi,&zvmi,&zvlo,acchi,accmi,acclo);         // zr = xr*yr - xi*yi

   zvrehi[k] = zvhi; zvremi[k] = zvmi; zvrelo[k] = zvlo;

   for(int i=1; i<=k; i++) // z[k] = z[k] + x[i]*y[k-i]
   {
      idx = k-i;
      xrhi = xvrehi[i]; xrmi = xvremi[i]; xrlo = xvrelo[i];
      xihi = xvimhi[i]; ximi = xvimmi[i]; xilo = xvimlo[i];
      yrhi = yvrehi[idx]; yrmi = yvremi[idx]; yrlo = yvrelo[idx];
      yihi = yvimhi[idx]; yimi = yvimmi[idx]; yilo = yvimlo[idx];

      tdg_mul(xrhi,xrmi,xrlo,yrhi,yrmi,yrlo,&zvhi,&zvmi,&zvlo); // zr = xr*yr
      tdg_mul(xihi,ximi,xilo,yihi,yimi,yilo,&acchi,&accmi,&acclo);   // xi*yi
      tdg_minus(&acchi,&accmi,&acclo);
      tdg_inc(&zvhi,&zvmi,&zvlo,acchi,accmi,acclo);     // zr = xr*yr - xi*yi
      // zvim[k] += zi
      tdg_inc(&zvrehi[k],&zvremi[k],&zvrelo[k],zvhi,zvmi,zvlo);
   }
                                                          // z[k] = x[0]*y[k]
   idx = dim1-k;
   xrhi = xvrehi[0]; xrmi = xvremi[0]; xrlo = xvrelo[0];
   xihi = xvimhi[0]; ximi = xvimmi[0]; xilo = xvimlo[0];
   yrhi = yvrehi[idx]; yrmi = yvremi[idx]; yrlo = yvrelo[idx];
   yihi = yvimhi[idx]; yimi = yvimmi[idx]; yilo = yvimlo[idx];

   tdg_mul(xrhi,xrmi,xrlo,yihi,yimi,yilo,&zvhi,&zvmi,&zvlo);     // zi = xr*yi
   tdg_mul(xihi,ximi,xilo,yrhi,yrmi,yrlo,&acchi,&accmi,&acclo); // acc = xi*yr
   tdg_inc(&zvhi,&zvmi,&zvlo,acchi,accmi,acclo);         // zi = xr*yi + xi*yr

   zvimhi[k] = zvhi; zvimmi[k] = zvmi; zvimlo[k] = zvlo;

   for(int i=1; i<dim-k; i++) // z[k] = z[k] + x[i]*y[k-i]
   {
      idx = dim1-k-i;
      xrhi = xvrehi[i]; xrmi = xvremi[i]; xrlo = xvrelo[i];
      xihi = xvimhi[i]; ximi = xvimmi[i]; xilo = xvimlo[i];
      yrhi = yvrehi[idx]; yrmi = yvremi[idx]; yrlo = yvrelo[idx];
      yihi = yvimhi[idx]; yimi = yvimmi[idx]; yilo = yvimlo[idx];

      tdg_mul(xrhi,xrmi,xrlo,yihi,yimi,yilo,&zvhi,&zvmi,&zvlo); // zi = xr*yi
      tdg_mul(xihi,ximi,xilo,yrhi,yrmi,yrlo,&acchi,&accmi,&acclo);   // xi*yr
      tdg_inc(&zvhi,&zvmi,&zvlo,acchi,accmi,acclo);     // zi = xr*yi + xi*yr
      // zvim[k] += zi
      tdg_inc(&zvimhi[k],&zvimmi[k],&zvimlo[k],zvhi,zvmi,zvlo);
   }
   __syncthreads();

   zrehi[k] = zvrehi[k]; zremi[k] = zvremi[k]; zrelo[k] = zvrelo[k];
   zimhi[k] = zvimhi[dim1-k]; zimmi[k] = zvimmi[dim1-k];
   zimlo[k] = zvimlo[dim1-k];
}

__global__ void cmplx3_padded_convolute
 ( double *xrehi, double *xremi, double *xrelo,
   double *ximhi, double *ximmi, double *ximlo,
   double *yrehi, double *yremi, double *yrelo,
   double *yimhi, double *yimmi, double *yimlo,
   double *zrehi, double *zremi, double *zrelo,
   double *zimhi, double *zimmi, double *zimlo, int dim )
{
   int k = threadIdx.x;       // thread k computes zre[k] and zim[k]

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

   double xrhi,xihi,yrhi,yihi,zrhi,zihi,acchi;
   double xrmi,ximi,yrmi,yimi,zrmi,zimi,accmi;
   double xrlo,xilo,yrlo,yilo,zrlo,zilo,acclo;
   int idx = dim+k;

   xvrehi[k] = xrehi[k];   xvimhi[k] = ximhi[k];
   xvremi[k] = xremi[k];   xvimmi[k] = ximmi[k];
   xvrelo[k] = xrelo[k];   xvimlo[k] = ximlo[k];
   yvrehi[k] = 0.0;        yvimhi[k] = 0.0;
   yvremi[k] = 0.0;        yvimmi[k] = 0.0;
   yvrelo[k] = 0.0;        yvimlo[k] = 0.0;
   yvrehi[idx] = yrehi[k]; yvimhi[idx] = yimhi[k];
   yvremi[idx] = yremi[k]; yvimmi[idx] = yimmi[k];
   yvrelo[idx] = yrelo[k]; yvimlo[idx] = yimlo[k];

   // z[k] = x[0]*y[k]
   xrhi = xvrehi[0]; xrmi = xvremi[0]; xrlo = xvrelo[0];
   xihi = xvimhi[0]; ximi = xvimmi[0]; xilo = xvimlo[0];
   yrhi = yvrehi[idx]; yrmi = yvremi[idx]; yrlo = yvrelo[idx];
   yihi = yvimhi[idx]; yimi = yvimmi[idx]; yilo = yvimlo[idx];

   tdg_mul(xrhi,xrmi,xrlo,yrhi,yrmi,yrlo,&zrhi,&zrmi,&zrlo);     // zr = xr*yr
   __syncthreads();
   tdg_mul(xihi,ximi,xilo,yihi,yimi,yilo,&acchi,&accmi,&acclo); // acc = xi*yi
   __syncthreads();
   tdg_minus(&acchi,&accmi,&acclo);
   tdg_inc(&zrhi,&zrmi,&zrlo,acchi,accmi,acclo);         // zr = xr*yr - xi*yi
   tdg_mul(xrhi,xrmi,xrlo,yihi,yimi,yilo,&zihi,&zimi,&zilo);     // zi = xr*yi
   __syncthreads();
   tdg_mul(xihi,ximi,xilo,yrhi,yrmi,yrlo,&acchi,&accmi,&acclo); // acc = xi*yr
   __syncthreads();
   tdg_inc(&zihi,&zimi,&zilo,acchi,accmi,acclo);         // zi = xr*yi + xi*yr

   zvrehi[k] = zrhi; zvremi[k] = zrmi; zvrelo[k] = zrlo;
   zvimhi[k] = zihi; zvimmi[k] = zimi; zvimlo[k] = zilo;

   for(int i=1; i<dim; i++) // z[k] = z[k] + x[i]*y[k-i]
   {
      idx = dim + k - i;
      xrhi = xvrehi[i]; xrmi = xvremi[i]; xrlo = xvrelo[i];
      xihi = xvimhi[i]; ximi = xvimmi[i]; xilo = xvimlo[i];
      yrhi = yvrehi[idx]; yrmi = yvremi[idx]; yrlo = yvrelo[idx];
      yihi = yvimhi[idx]; yimi = yvimmi[idx]; yilo = yvimlo[idx];

      tdg_mul(xrhi,xrmi,xrlo,yrhi,yrmi,yrlo,&zrhi,&zrmi,&zrlo); // zr = xr*yr
      __syncthreads();
      tdg_mul(xihi,ximi,xilo,yihi,yimi,yilo,&acchi,&accmi,&acclo);   // xi*yi
      __syncthreads();
      tdg_minus(&acchi,&accmi,&acclo);
      tdg_inc(&zrhi,&zrmi,&zrlo,acchi,accmi,acclo);     // zr = xr*yr - xi*yi
      tdg_mul(xrhi,xrmi,xrlo,yihi,yimi,yilo,&zihi,&zimi,&zilo); // zi = xr*yi
      __syncthreads();
      tdg_mul(xihi,ximi,xilo,yrhi,yrmi,yrlo,&acchi,&accmi,&acclo);   // xi*yr
      __syncthreads();
      tdg_inc(&zihi,&zimi,&zilo,acchi,accmi,acclo);     // zr = xr*yr + xi*yi
      // zvre[k] += zr; zvim[k] += zi
      tdg_inc(&zvrehi[k],&zvremi[k],&zvrelo[k],zrhi,zrmi,zrlo);
      tdg_inc(&zvimhi[k],&zvimmi[k],&zvimlo[k],zihi,zimi,zilo);
   }
   zrehi[k] = zvrehi[k]; zremi[k] = zvremi[k]; zrelo[k] = zvrelo[k];
   zimhi[k] = zvimhi[k]; zimmi[k] = zvimmi[k]; zimlo[k] = zvimlo[k];
}

void GPU_dbl3_product
 ( double *xhi_h, double *xmi_h, double *xlo_h,
   double *yhi_h, double *ymi_h, double *ylo_h,
   double *zhi_h, double *zmi_h, double *zlo_h, int deg, int freq, int BS,
   int padded )
{
   const int dim = deg+1;            // length of all vectors
   double* xhi_d;                    // xhi_d is xhi_h on the device
   double* xmi_d;                    // xmi_d is xmi_h on the device
   double* xlo_d;                    // xlo_d is xlo_h on the device
   double* yhi_d;                    // yhi_d is yhi_h on the device
   double* ymi_d;                    // ymi_d is ymi_h on the device
   double* ylo_d;                    // ylo_d is ylo_h on the device
   double* zhi_d;                    // zhi_d is zhi_h on the device
   double* zmi_d;                    // zmi_d is zmi_h on the device
   double* zlo_d;                    // zlo_d is zlo_h on the device
   size_t size = dim*sizeof(double); // number of bytes for each vector

   cudaMalloc((void**)&xhi_d,size);
   cudaMalloc((void**)&xmi_d,size);
   cudaMalloc((void**)&xlo_d,size);
   cudaMalloc((void**)&yhi_d,size);
   cudaMalloc((void**)&ymi_d,size);
   cudaMalloc((void**)&ylo_d,size);
   cudaMalloc((void**)&zhi_d,size);
   cudaMalloc((void**)&zmi_d,size);
   cudaMalloc((void**)&zlo_d,size);
   cudaMemcpy(xhi_d,xhi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(xmi_d,xmi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(xlo_d,xlo_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yhi_d,yhi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(ymi_d,ymi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(ylo_d,ylo_h,size,cudaMemcpyHostToDevice);

   if(dim == BS)
   {
      if(padded == 1)
      {
         for(int i=0; i<freq; i++)
            dbl3_padded_convolute<<<1,BS>>>
               (xhi_d,xmi_d,xlo_d,yhi_d,ymi_d,ylo_d,zhi_d,zmi_d,zlo_d,dim);
      }
      else
      {
         for(int i=0; i<freq; i++)
            dbl3_convolute<<<1,BS>>>
               (xhi_d,xmi_d,xlo_d,yhi_d,ymi_d,ylo_d,zhi_d,zmi_d,zlo_d,dim);
      }
   }
   cudaMemcpy(zhi_h,zhi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zmi_h,zmi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zlo_h,zlo_d,size,cudaMemcpyDeviceToHost);
}

void GPU_cmplx3_product
 ( double *xrehi_h, double *xremi_h, double *xrelo_h,
   double *ximhi_h, double *ximmi_h, double *ximlo_h,
   double *yrehi_h, double *yremi_h, double *yrelo_h,
   double *yimhi_h, double *yimmi_h, double *yimlo_h,
   double *zrehi_h, double *zremi_h, double *zrelo_h,
   double *zimhi_h, double *zimmi_h, double *zimlo_h,
   int deg, int freq, int BS, int mode )
{
   const int dim = deg+1;            // length of all vectors
   double* xrehi_d;                  // xrehi_d is xrehi_h on the device
   double* xremi_d;                  // xremi_d is xremi_h on the device
   double* xrelo_d;                  // xrelo_d is xrelo_h on the device
   double* ximhi_d;                  // ximhi_d is ximhi_h on the device
   double* ximmi_d;                  // ximmi_d is ximmi_h on the device
   double* ximlo_d;                  // ximlo_d is ximlo_h on the device
   double* yrehi_d;                  // yrehi_d is yrehi_h on the device
   double* yremi_d;                  // yremi_d is yremi_h on the device
   double* yrelo_d;                  // yrelo_d is yrelo_h on the device
   double* yimhi_d;                  // yimhi_d is yimhi_h on the device
   double* yimmi_d;                  // yimmi_d is yimmi_h on the device
   double* yimlo_d;                  // yimlo_d is yimlo_h on the device
   double* zrehi_d;                  // zrehi_d is zrehi_h on the device
   double* zremi_d;                  // zremi_d is zremi_h on the device
   double* zrelo_d;                  // zrelo_d is zrelo_h on the device
   double* zimhi_d;                  // zimhi_d is zimhi_h on the device
   double* zimmi_d;                  // zimmi_d is zimmi_h on the device
   double* zimlo_d;                  // zimlo_d is zimlo_h on the device
   double* acchi_d;                  // accumulates high parts on the device
   double* accmi_d;                  // accumulates middle parts
   double* acclo_d;                  // accumulates low parts 
   size_t size = dim*sizeof(double); // number of bytes for each vector

   cudaMalloc((void**)&xrehi_d,size);
   cudaMalloc((void**)&xremi_d,size);
   cudaMalloc((void**)&xrelo_d,size);
   cudaMalloc((void**)&ximhi_d,size);
   cudaMalloc((void**)&ximmi_d,size);
   cudaMalloc((void**)&ximlo_d,size);
   cudaMalloc((void**)&yrehi_d,size);
   cudaMalloc((void**)&yremi_d,size);
   cudaMalloc((void**)&yrelo_d,size);
   cudaMalloc((void**)&yimhi_d,size);
   cudaMalloc((void**)&yimmi_d,size);
   cudaMalloc((void**)&yimlo_d,size);
   cudaMalloc((void**)&zrehi_d,size);
   cudaMalloc((void**)&zremi_d,size);
   cudaMalloc((void**)&zrelo_d,size);
   cudaMalloc((void**)&zimhi_d,size);
   cudaMalloc((void**)&zimmi_d,size);
   cudaMalloc((void**)&zimlo_d,size);
   cudaMalloc((void**)&acchi_d,size);
   cudaMalloc((void**)&accmi_d,size);
   cudaMalloc((void**)&acclo_d,size);
   cudaMemcpy(xrehi_d,xrehi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(xremi_d,xremi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(xrelo_d,xrelo_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(ximhi_d,ximhi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(ximmi_d,ximmi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(ximlo_d,ximlo_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yrehi_d,yrehi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yremi_d,yremi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yrelo_d,yrelo_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yimhi_d,yimhi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yimmi_d,yimmi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yimlo_d,yimlo_h,size,cudaMemcpyHostToDevice);

   if(dim == BS)
   {
      if(mode == 2)
      {
         for(int i=0; i<freq; i++)
         {
            cmplx3_padded_convolute<<<1,BS>>>
               (xrehi_d,xremi_d,xrelo_d,ximhi_d,ximmi_d,ximlo_d,
                yrehi_d,yremi_d,yrelo_d,yimhi_d,yimmi_d,yimlo_d,
                zrehi_d,zremi_d,zrelo_d,zimhi_d,zimmi_d,zimlo_d,dim);
         }
      }
      else if(mode == 2)
      {
         for(int i=0; i<freq; i++)
         {
            dbl3_padded_convolute<<<1,BS>>>
               (xrehi_d,xremi_d,xrelo_d,yrehi_d,yremi_d,yrelo_d,
                zrehi_d,zremi_d,zrelo_d,dim);
            dbl3_padded_convolute<<<1,BS>>>
               (ximhi_d,ximmi_d,ximlo_d,yimhi_d,yimmi_d,yimlo_d,
                acchi_d,accmi_d,acclo_d,dim);
            dbl3_decrement<<<1,BS>>>
               (zrehi_d,zremi_d,zrelo_d,acchi_d,accmi_d,acclo_d,
                zrehi_d,zremi_d,zrelo_d,dim);
            dbl3_padded_convolute<<<1,BS>>>
               (xrehi_d,xremi_d,xrelo_d,yimhi_d,yimmi_d,yimlo_d,
                zimhi_d,zimmi_d,zimlo_d,dim);
            dbl3_padded_convolute<<<1,BS>>>
               (ximhi_d,ximmi_d,ximlo_d,yrehi_d,yremi_d,yrelo_d,
                acchi_d,accmi_d,acclo_d,dim);
            dbl3_increment<<<1,BS>>>
               (zimhi_d,zimmi_d,zimlo_d,acchi_d,accmi_d,acclo_d,
                zimhi_d,zimmi_d,zimlo_d,dim);
         }
      }
      else if(mode == 1)
      {
         for(int i=0; i<freq; i++)
            cmplx3_looped_convolute<<<1,BS>>>
               (xrehi_d,xremi_d,xrelo_d,ximhi_d,ximmi_d,ximlo_d,
                yrehi_d,yremi_d,yrelo_d,yimhi_d,yimmi_d,yimlo_d,
                zrehi_d,zremi_d,zrelo_d,zimhi_d,zimmi_d,zimlo_d,dim);
      }
      else
      {
         for(int i=0; i<freq; i++)
            cmplx3_convolute<<<1,BS>>>
               (xrehi_d,xremi_d,xrelo_d,ximhi_d,ximmi_d,ximlo_d,
                yrehi_d,yremi_d,yrelo_d,yimhi_d,yimmi_d,yimlo_d,
                zrehi_d,zremi_d,zrelo_d,zimhi_d,zimmi_d,zimlo_d,dim);
      }

   }

   cudaMemcpy(zrehi_h,zrehi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zremi_h,zremi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zrelo_h,zrelo_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zimhi_h,zimhi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zimmi_h,zimmi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zimlo_h,zimlo_d,size,cudaMemcpyDeviceToHost);
}
