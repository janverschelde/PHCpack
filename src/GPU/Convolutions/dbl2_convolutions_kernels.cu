// The file dbl2_convolutions_kernels.cu defines kernels with prototypes
// in dbl2_convolutions_kernels.h.

#ifdef gpufun
#include "double_double_gpufun.cu"
#endif
#include "dbl2_convolutions_kernels.h"

__global__ void dbl2_increment
 ( double *xhi, double *xlo, double *yhi, double *ylo,
   double *zhi, double *zlo, int dim )
{
   int k = threadIdx.x;                 // thread k computes z[k]

   __shared__ double xvhi[dd_shmemsize];
   __shared__ double xvlo[dd_shmemsize];
   __shared__ double yvhi[dd_shmemsize];
   __shared__ double yvlo[dd_shmemsize];
   __shared__ double zvhi[dd_shmemsize];
   __shared__ double zvlo[dd_shmemsize];

   xvhi[k] = xhi[k]; xvlo[k] = xlo[k];
   yvhi[k] = yhi[k]; yvlo[k] = ylo[k];

   ddg_add(xvhi[k],xvlo[k],yvhi[k],yvlo[k],&zvhi[k],&zvlo[k]);

   __syncthreads();

   zhi[k] = zvhi[k];
   zlo[k] = zvlo[k];
}

__global__ void dbl2_decrement
 ( double *xhi, double *xlo, double *yhi, double *ylo,
   double *zhi, double *zlo, int dim )
{
   int k = threadIdx.x;                 // thread k computes z[k]

   __shared__ double xvhi[dd_shmemsize];
   __shared__ double xvlo[dd_shmemsize];
   __shared__ double yvhi[dd_shmemsize];
   __shared__ double yvlo[dd_shmemsize];
   __shared__ double zvhi[dd_shmemsize];
   __shared__ double zvlo[dd_shmemsize];

   xvhi[k] = xhi[k]; xvlo[k] = xlo[k];
   yvhi[k] = yhi[k]; yvlo[k] = ylo[k];

   ddg_sub(xvhi[k],xvlo[k],yvhi[k],yvlo[k],&zvhi[k],&zvlo[k]);

   __syncthreads();

   zhi[k] = zvhi[k];
   zlo[k] = zvlo[k];
}

__global__ void dbl2_convolute
 ( double *xhi, double *xlo, double *yhi, double *ylo,
   double *zhi, double *zlo, int dim )
{
   int k = threadIdx.x;                 // thread k computes z[k]

   __shared__ double xvhi[dd_shmemsize];
   __shared__ double xvlo[dd_shmemsize];
   __shared__ double yvhi[dd_shmemsize];
   __shared__ double yvlo[dd_shmemsize];
   __shared__ double zvhi[dd_shmemsize];
   __shared__ double zvlo[dd_shmemsize];

   double prdhi,prdlo;

   xvhi[k] = xhi[k]; xvlo[k] = xlo[k];
   yvhi[k] = yhi[k]; yvlo[k] = ylo[k];

   __syncthreads();

   // zv[k] = xv[0]*yv[k];
   ddg_mul(xvhi[0],xvlo[0],yvhi[k],yvlo[k],&zvhi[k],&zvlo[k]);

   for(int i=1; i<=k; i++) // zv[k] = zv[k] + xv[i]*yv[k-i];
   {
      ddg_mul(xvhi[i],xvlo[i],yvhi[k-i],yvlo[k-i],&prdhi,&prdlo);
      ddg_inc(&zvhi[k],&zvlo[k],prdhi,prdlo);
   }

   __syncthreads();

   zhi[k] = zvhi[k];
   zlo[k] = zvlo[k];
}

__global__ void dbl2_padded_convolute
 ( double *xhi, double *xlo, double *yhi, double *ylo,
   double *zhi, double *zlo, int dim )
{
   int k = threadIdx.x;                 // thread k computes z[k]

   __shared__ double xvhi[dd_shmemsize];
   __shared__ double xvlo[dd_shmemsize];
   __shared__ double yvhi[dd_shmemsize];
   __shared__ double yvlo[dd_shmemsize];
   __shared__ double zvhi[dd_shmemsize];
   __shared__ double zvlo[dd_shmemsize];

   double prdhi,prdlo;
   int idx = dim+k;

   xvhi[k] = xhi[k];   xvlo[k] = xlo[k];
   yvhi[k] = 0.0;      yvlo[k] = 0.0;      // padding with zeros
   yvhi[idx] = yhi[k]; yvlo[idx] = ylo[k];

   __syncthreads();

   // zv[k] = xv[0]*yv[k];
   ddg_mul(xvhi[0],xvlo[0],yvhi[idx],yvlo[idx],&zvhi[k],&zvlo[k]);

   for(int i=1; i<dim; i++) // zv[k] = zv[k] + xv[i]*yv[k-i];
   {
      idx = dim + k - i;
      ddg_mul(xvhi[i],xvlo[i],yvhi[idx],yvlo[idx],&prdhi,&prdlo);
      ddg_inc(&zvhi[k],&zvlo[k],prdhi,prdlo);
   }
   __syncthreads();

   zhi[k] = zvhi[k]; zlo[k] = zvlo[k];
}

__global__ void cmplx2_convolute
 ( double *xrehi, double *xrelo, double *ximhi, double *ximlo,
   double *yrehi, double *yrelo, double *yimhi, double *yimlo,
   double *zrehi, double *zrelo, double *zimhi, double *zimlo, int dim )
{
   int k = threadIdx.x;       // thread k computes zre[k] and zim[k]

   __shared__ double xvrehi[dd_shmemsize];
   __shared__ double xvrelo[dd_shmemsize];
   __shared__ double xvimhi[dd_shmemsize];
   __shared__ double xvimlo[dd_shmemsize];
   __shared__ double yvrehi[dd_shmemsize];
   __shared__ double yvrelo[dd_shmemsize];
   __shared__ double yvimhi[dd_shmemsize];
   __shared__ double yvimlo[dd_shmemsize];
   __shared__ double zvrehi[dd_shmemsize];
   __shared__ double zvrelo[dd_shmemsize];
   __shared__ double zvimhi[dd_shmemsize];
   __shared__ double zvimlo[dd_shmemsize];

   double xrhi,xihi,yrhi,yihi,zrhi,zihi,acchi;
   double xrlo,xilo,yrlo,yilo,zrlo,zilo,acclo;

   xvrehi[k] = xrehi[k]; xvimhi[k] = ximhi[k];
   xvrelo[k] = xrelo[k]; xvimlo[k] = ximlo[k];
   yvrehi[k] = yrehi[k]; yvimhi[k] = yimhi[k];
   yvrelo[k] = yrelo[k]; yvimlo[k] = yimlo[k];

   __syncthreads();

   // z[k] = x[0]*y[k]
   xrhi = xvrehi[0]; xrlo = xvrelo[0]; xihi = xvimhi[0]; xilo = xvimlo[0];
   yrhi = yvrehi[k]; yrlo = yvrelo[k]; yihi = yvimhi[k]; yilo = yvimlo[k];

   ddg_mul(xrhi,xrlo,yrhi,yrlo,&zrhi,&zrlo);   // zr = xr*yr
   ddg_mul(xihi,xilo,yihi,yilo,&acchi,&acclo); // acc = xi*yi
   ddg_dec(&zrhi,&zrlo,acchi,acclo);           // zr = xr*yr - xi*yi
   ddg_mul(xrhi,xrlo,yihi,yilo,&zihi,&zilo);   // zi = xr*yi
   ddg_mul(xihi,xilo,yrhi,yrlo,&acchi,&acclo); // acc = xi*yr
   ddg_inc(&zihi,&zilo,acchi,acclo);           // zi = xr*yi + xi*yr

   zvrehi[k] = zrhi; zvrelo[k] = zrlo;
   zvimhi[k] = zihi; zvimlo[k] = zilo;

   for(int i=1; i<=k; i++) // z[k] = z[k] + x[i]*y[k-i]
   {
      xrhi = xvrehi[i]; xrlo = xvrelo[i];
      xihi = xvimhi[i]; xilo = xvimlo[i];
      yrhi = yvrehi[k-i]; yrlo = yvrelo[k-i];
      yihi = yvimhi[k-i]; yilo = yvimlo[k-i];

      ddg_mul(xrhi,xrlo,yrhi,yrlo,&zrhi,&zrlo);   // zr = xr*yr
      ddg_mul(xihi,xilo,yihi,yilo,&acchi,&acclo); // acc = xi*yi
      ddg_dec(&zrhi,&zrlo,acchi,acclo);           // zr = xr*yr - xi*yi
      ddg_mul(xrhi,xrlo,yihi,yilo,&zihi,&zilo);   // zi = xr*yi
      ddg_mul(xihi,xilo,yrhi,yrlo,&acchi,&acclo); // acc = xi*yr
      ddg_inc(&zihi,&zilo,acchi,acclo);           // zi = xr*yi + xi*yr

      ddg_inc(&zvrehi[k],&zvrelo[k],zrhi,zrlo);   // zvre[k] += zr;
      ddg_inc(&zvimhi[k],&zvimlo[k],zihi,zilo);   // zvim[k] += zi;
   }

   __syncthreads();

   zrehi[k] = zvrehi[k]; zrelo[k] = zvrelo[k];
   zimhi[k] = zvimhi[k]; zimlo[k] = zvimlo[k];
}

__global__ void cmplx2_padded_convolute
 ( double *xrehi, double *xrelo, double *ximhi, double *ximlo,
   double *yrehi, double *yrelo, double *yimhi, double *yimlo,
   double *zrehi, double *zrelo, double *zimhi, double *zimlo, int dim )
{
   int k = threadIdx.x;       // thread k computes zre[k] and zim[k]

   __shared__ double xvrehi[dd_shmemsize];
   __shared__ double xvrelo[dd_shmemsize];
   __shared__ double xvimhi[dd_shmemsize];
   __shared__ double xvimlo[dd_shmemsize];
   __shared__ double yvrehi[dd_shmemsize];
   __shared__ double yvrelo[dd_shmemsize];
   __shared__ double yvimhi[dd_shmemsize];
   __shared__ double yvimlo[dd_shmemsize];
   __shared__ double zvrehi[dd_shmemsize];
   __shared__ double zvrelo[dd_shmemsize];
   __shared__ double zvimhi[dd_shmemsize];
   __shared__ double zvimlo[dd_shmemsize];

   double xrhi,xihi,yrhi,yihi,zrhi,zihi,acchi;
   double xrlo,xilo,yrlo,yilo,zrlo,zilo,acclo;
   int idx = dim+k;

   xvrehi[k] = xrehi[k];   xvimhi[k] = ximhi[k];
   xvrelo[k] = xrelo[k];   xvimlo[k] = ximlo[k];
   yvrehi[k] = 0.0;        yvimhi[k] = 0.0;
   yvrelo[k] = 0.0;        yvimlo[k] = 0.0;
   yvrehi[idx] = yrehi[k]; yvimhi[idx] = yimhi[k];
   yvrelo[idx] = yrelo[k]; yvimlo[idx] = yimlo[k];

   __syncthreads();

   // z[k] = x[0]*y[k]
   xrhi = xvrehi[0]; xrlo = xvrelo[0];
   xihi = xvimhi[0]; xilo = xvimlo[0];
   yrhi = yvrehi[idx]; yrlo = yvrelo[idx];
   yihi = yvimhi[idx]; yilo = yvimlo[idx];

   ddg_mul(xrhi,xrlo,yrhi,yrlo,&zrhi,&zrlo);   // zr = xr*yr
   __syncthreads();
   ddg_mul(xihi,xilo,yihi,yilo,&acchi,&acclo); // acc = xi*yi
   __syncthreads();
   ddg_dec(&zrhi,&zrlo,acchi,acclo);           // zr = xr*yr - xi*yi
   __syncthreads();
   ddg_mul(xrhi,xrlo,yihi,yilo,&zihi,&zilo);   // zi = xr*yi
   __syncthreads();
   ddg_mul(xihi,xilo,yrhi,yrlo,&acchi,&acclo); // acc = xi*yr
   __syncthreads();
   ddg_inc(&zihi,&zilo,acchi,acclo);           // zi = xr*yi + xi*yr
   __syncthreads();
   zvrehi[k] = zrhi; zvrelo[k] = zrlo;
   zvimhi[k] = zihi; zvimlo[k] = zilo;

   for(int i=1; i<dim; i++) // z[k] = z[k] + x[i]*y[k-i]
   {
      idx = dim + k - i;
      xrhi = xvrehi[i]; xrlo = xvrelo[i];
      xihi = xvimhi[i]; xilo = xvimlo[i];
      yrhi = yvrehi[idx]; yrlo = yvrelo[idx];
      yihi = yvimhi[idx]; yilo = yvimlo[idx];

      ddg_mul(xrhi,xrlo,yrhi,yrlo,&zrhi,&zrlo);   // zr = xr*yr
      __syncthreads();
      ddg_mul(xihi,xilo,yihi,yilo,&acchi,&acclo); // acc = xi*yi
      __syncthreads();
      ddg_dec(&zrhi,&zrlo,acchi,acclo);           // zr = xr*yr - xi*yi
      __syncthreads();
      ddg_mul(xrhi,xrlo,yihi,yilo,&zihi,&zilo);   // zi = xr*yi
      __syncthreads();
      ddg_mul(xihi,xilo,yrhi,yrlo,&acchi,&acclo); // acc = xi*yr
      __syncthreads();
      ddg_inc(&zihi,&zilo,acchi,acclo);           // zi = xr*yi + xi*yr

      __syncthreads();
      ddg_inc(&zvrehi[k],&zvrelo[k],zrhi,zrlo);   // zvre[k] += zr;
      __syncthreads();
      ddg_inc(&zvimhi[k],&zvimlo[k],zihi,zilo);   // zvim[k] += zi;
   }
   __syncthreads();
   zrehi[k] = zvrehi[k]; zrelo[k] = zvrelo[k];
   zimhi[k] = zvimhi[k]; zimlo[k] = zvimlo[k];
}

__global__ void cmplx2_looped_convolute
 ( double *xrehi, double *xrelo, double *ximhi, double *ximlo,
   double *yrehi, double *yrelo, double *yimhi, double *yimlo,
   double *zrehi, double *zrelo, double *zimhi, double *zimlo, int dim )
{
   int k = threadIdx.x;       // thread k computes zre[k] and zim[dim-1-k]

   __shared__ double xvrehi[dd_shmemsize];
   __shared__ double xvrelo[dd_shmemsize];
   __shared__ double xvimhi[dd_shmemsize];
   __shared__ double xvimlo[dd_shmemsize];
   __shared__ double yvrehi[dd_shmemsize];
   __shared__ double yvrelo[dd_shmemsize];
   __shared__ double yvimhi[dd_shmemsize];
   __shared__ double yvimlo[dd_shmemsize];
   __shared__ double zvrehi[dd_shmemsize];
   __shared__ double zvrelo[dd_shmemsize];
   __shared__ double zvimhi[dd_shmemsize];
   __shared__ double zvimlo[dd_shmemsize];

   const int dim1 = dim-1;
   int idx;
   double xrhi,xihi,yrhi,yihi,zvhi,acchi;
   double xrlo,xilo,yrlo,yilo,zvlo,acclo;

   xvrehi[k] = xrehi[k]; xvimhi[k] = ximhi[k];
   xvrelo[k] = xrelo[k]; xvimlo[k] = ximlo[k];
   yvrehi[k] = yrehi[k]; yvimhi[k] = yimhi[k];
   yvrelo[k] = yrelo[k]; yvimlo[k] = yimlo[k];

   __syncthreads();

   // z[k] = x[0]*y[k]
   xrhi = xvrehi[0]; xrlo = xvrelo[0]; xihi = xvimhi[0]; xilo = xvimlo[0];
   yrhi = yvrehi[k]; yrlo = yvrelo[k]; yihi = yvimhi[k]; yilo = yvimlo[k];

   ddg_mul(xrhi,xrlo,yrhi,yrlo,&zvhi,&zvlo);   // zr = xr*yr
   ddg_mul(xihi,xilo,yihi,yilo,&acchi,&acclo); // acc = xi*yi
   ddg_dec(&zvhi,&zvlo,acchi,acclo);           // zr = xr*yr - xi*yi

   zvrehi[k] = zvhi; zvrelo[k] = zvlo;

   for(int i=1; i<=k; i++) // z[k] = z[k] + x[i]*y[k-i]
   {
      idx = k-i;
      xrhi = xvrehi[i]; xrlo = xvrelo[i];
      xihi = xvimhi[i]; xilo = xvimlo[i];
      yrhi = yvrehi[idx]; yrlo = yvrelo[idx];
      yihi = yvimhi[idx]; yilo = yvimlo[idx];

      ddg_mul(xrhi,xrlo,yrhi,yrlo,&zvhi,&zvlo);   // zr = xr*yr
      ddg_mul(xihi,xilo,yihi,yilo,&acchi,&acclo); // acc = xi*yi
      ddg_dec(&zvhi,&zvlo,acchi,acclo);           // zr = xr*yr - xi*yi

      ddg_inc(&zvrehi[k],&zvrelo[k],zvhi,zvlo);   // zvre[k] += zr;
   }

   // z[k] = x[0]*y[k]
   xrhi = xvrehi[0]; xrlo = xvrelo[0];
   xihi = xvimhi[0]; xilo = xvimlo[0];
   idx = dim1-k;
   yrhi = yvrehi[idx]; yrlo = yvrelo[idx];
   yihi = yvimhi[idx]; yilo = yvimlo[idx];

   ddg_mul(xrhi,xrlo,yihi,yilo,&zvhi,&zvlo);   // zi = xr*yi
   ddg_mul(xihi,xilo,yrhi,yrlo,&acchi,&acclo); // acc = xi*yr
   ddg_inc(&zvhi,&zvlo,acchi,acclo);           // zi = xr*yi + xi*yr

   zvimhi[k] = zvhi; zvimlo[k] = zvlo;

   for(int i=1; i<dim-k; i++) // z[k] = z[k] + x[i]*y[k-i]
   {
      idx = dim1-k-i;
      xrhi = xvrehi[i]; xrlo = xvrelo[i];
      xihi = xvimhi[i]; xilo = xvimlo[i];
      yrhi = yvrehi[idx]; yrlo = yvrelo[idx];
      yihi = yvimhi[idx]; yilo = yvimlo[idx];

      ddg_mul(xrhi,xrlo,yihi,yilo,&zvhi,&zvlo);   // zi = xr*yi
      ddg_mul(xihi,xilo,yrhi,yrlo,&acchi,&acclo); // acc = xi*yr
      ddg_inc(&zvhi,&zvlo,acchi,acclo);           // zr = xr*yi + xi*yr

      ddg_inc(&zvimhi[k],&zvimlo[k],zvhi,zvlo);   // zvim[k] += zi;
   }

   __syncthreads();

   zrehi[k] = zvrehi[k]; zrelo[k] = zvrelo[k];
   zimhi[k] = zvimhi[dim1-k];
   zimlo[k] = zvimlo[dim1-k];
}

void GPU_dbl2_product
 ( double *xhi_h, double *xlo_h, double *yhi_h, double *ylo_h,
   double *zhi_h, double *zlo_h, int deg, int freq, int BS, int padded )
{
   const int dim = deg+1;            // length of all vectors
   double* xhi_d;                    // xhi_d is xhi_h on the device
   double* xlo_d;                    // xlo_d is xlo_h on the device
   double* yhi_d;                    // yhi_d is yhi_h on the device
   double* ylo_d;                    // ylo_d is ylo_h on the device
   double* zhi_d;                    // zhi_d is zhi_h on the device
   double* zlo_d;                    // zlo_d is zlo_h on the device
   size_t size = dim*sizeof(double); // number of bytes for each vector

   cudaMalloc((void**)&xhi_d,size);
   cudaMalloc((void**)&xlo_d,size);
   cudaMalloc((void**)&yhi_d,size);
   cudaMalloc((void**)&ylo_d,size);
   cudaMalloc((void**)&zhi_d,size);
   cudaMalloc((void**)&zlo_d,size);
   cudaMemcpy(xhi_d,xhi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(xlo_d,xlo_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yhi_d,yhi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(ylo_d,ylo_h,size,cudaMemcpyHostToDevice);

   if(dim == BS)
   {
      if(padded == 1)
      {
         for(int i=0; i<freq; i++)
            dbl2_padded_convolute<<<1,BS>>>
               (xhi_d,xlo_d,yhi_d,ylo_d,zhi_d,zlo_d,dim);
      }
      else
      {
         for(int i=0; i<freq; i++)
            dbl2_convolute<<<1,BS>>>
               (xhi_d,xlo_d,yhi_d,ylo_d,zhi_d,zlo_d,dim);
      }
   }
   cudaMemcpy(zhi_h,zhi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zlo_h,zlo_d,size,cudaMemcpyDeviceToHost);
}

void GPU_cmplx2_product
 ( double *xrehi_h, double *xrelo_h, double *ximhi_h, double *ximlo_h,
   double *yrehi_h, double *yrelo_h, double *yimhi_h, double *yimlo_h,
   double *zrehi_h, double *zrelo_h, double *zimhi_h, double *zimlo_h,
   int deg, int freq, int BS, int mode )
{
   const int dim = deg+1;            // length of all vectors
   double* xrehi_d;                  // xrehi_d is xrehi_h on the device
   double* xrelo_d;                  // xrelo_d is xrelo_h on the device
   double* ximhi_d;                  // ximhi_d is ximhi_h on the device
   double* ximlo_d;                  // ximlo_d is ximlo_h on the device
   double* yrehi_d;                  // yrehi_d is yrehi_h on the device
   double* yrelo_d;                  // yrelo_d is yrelo_h on the device
   double* yimhi_d;                  // yimhi_d is yimhi_h on the device
   double* yimlo_d;                  // yimlo_d is yimlo_h on the device
   double* zrehi_d;                  // zrehi_d is zrehi_h on the device
   double* zrelo_d;                  // zrelo_d is zrelo_h on the device
   double* zimhi_d;                  // zimhi_d is zimhi_h on the device
   double* zimlo_d;                  // zimlo_d is zimlo_h on the device
   double* acchi_d;                  // accumulates high parts
   double* acclo_d;                  // accumulates low parts
   size_t size = dim*sizeof(double); // number of bytes for each vector

   cudaMalloc((void**)&xrehi_d,size);
   cudaMalloc((void**)&xrelo_d,size);
   cudaMalloc((void**)&ximhi_d,size);
   cudaMalloc((void**)&ximlo_d,size);
   cudaMalloc((void**)&yrehi_d,size);
   cudaMalloc((void**)&yrelo_d,size);
   cudaMalloc((void**)&yimhi_d,size);
   cudaMalloc((void**)&yimlo_d,size);
   cudaMalloc((void**)&zrehi_d,size);
   cudaMalloc((void**)&zrelo_d,size);
   cudaMalloc((void**)&zimhi_d,size);
   cudaMalloc((void**)&zimlo_d,size);
   cudaMalloc((void**)&acchi_d,size);
   cudaMalloc((void**)&acclo_d,size);
   cudaMemcpy(xrehi_d,xrehi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(xrelo_d,xrelo_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(ximhi_d,ximhi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(ximlo_d,ximlo_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yrehi_d,yrehi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yrelo_d,yrelo_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yimhi_d,yimhi_h,size,cudaMemcpyHostToDevice);
   cudaMemcpy(yimlo_d,yimlo_h,size,cudaMemcpyHostToDevice);

   if(dim == BS)
   {
      if(mode == 3)
      {
         for(int i=0; i<freq; i++)
         {
            cmplx2_padded_convolute<<<1,BS>>>
               (xrehi_d,xrelo_d,ximhi_d,ximlo_d,
                yrehi_d,yrelo_d,yimhi_d,yimlo_d,
                zrehi_d,zrelo_d,zimhi_d,zimlo_d,dim);
         }
      }
      else if(mode == 2)
      {
         for(int i=0; i<freq; i++)
         {
            dbl2_padded_convolute<<<1,BS>>>
               (xrehi_d,xrelo_d,yrehi_d,yrelo_d,zrehi_d,zrelo_d,dim);
            dbl2_padded_convolute<<<1,BS>>>
               (ximhi_d,ximlo_d,yimhi_d,yimlo_d,acchi_d,acclo_d,dim);
            dbl2_decrement<<<1,BS>>>
               (zrehi_d,zrelo_d,acchi_d,acclo_d,zrehi_d,zrelo_d,dim);
            dbl2_padded_convolute<<<1,BS>>>
               (xrehi_d,xrelo_d,yimhi_d,yimlo_d,zimhi_d,zimlo_d,dim);
            dbl2_padded_convolute<<<1,BS>>>
               (ximhi_d,ximlo_d,yrehi_d,yrelo_d,acchi_d,acclo_d,dim);
            dbl2_increment<<<1,BS>>>
               (zimhi_d,zimlo_d,acchi_d,acclo_d,zimhi_d,zimlo_d,dim);
         }
      }
      else if(mode == 1)
      {
         for(int i=0; i<freq; i++)
            cmplx2_looped_convolute<<<1,BS>>>
               (xrehi_d,xrelo_d,ximhi_d,ximlo_d,
                yrehi_d,yrelo_d,yimhi_d,yimlo_d,
                zrehi_d,zrelo_d,zimhi_d,zimlo_d,dim);
      }
      else
      {
         for(int i=0; i<freq; i++)
            cmplx2_convolute<<<1,BS>>>
               (xrehi_d,xrelo_d,ximhi_d,ximlo_d,
                yrehi_d,yrelo_d,yimhi_d,yimlo_d,
                zrehi_d,zrelo_d,zimhi_d,zimlo_d,dim);
      }
   }
   cudaMemcpy(zrehi_h,zrehi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zrelo_h,zrelo_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zimhi_h,zimhi_d,size,cudaMemcpyDeviceToHost);
   cudaMemcpy(zimlo_h,zimlo_d,size,cudaMemcpyDeviceToHost);
}
