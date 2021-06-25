// The file dbl_matrices_kernels.cu defines the kernels with prototypes
// in dbl_matrices_kernels.h.

#include <iostream>
#include <cmath>
#include "dbl_matrices_kernels.h"

__global__ void dbl_convolutions ( double *x, double *y, int deg1 )
{
   const int j = blockIdx.x;     // convolution of j-th series in x and y
   const int k = threadIdx.x;    // thread k computes k-th coefficient 
   const int offset = j*deg1+k;  // position of the k-th coefficient

   __shared__ double xv[d_shmemsize];
   __shared__ double yv[d_shmemsize];
   __shared__ double zv[d_shmemsize];

   int idx = deg1+k;

   xv[k] = x[offset];
   yv[k] = 0.0;         // padded with zeros
   yv[idx] = y[offset];

   zv[k] = xv[0]*yv[idx];

   for(int i=1; i<deg1; i++)
   {
      idx = deg1 + k - i;
      zv[k] = zv[k] + xv[i]*yv[idx];
   }
   x[offset] = zv[k];
}

__global__ void cmplx_convolutions
 ( double *xre, double *xim, double *yre, double *yim, int deg1 )
{
   const int j = blockIdx.x;     // convolution of j-th series in x and y
   const int k = threadIdx.x;    // thread k computes k-th coefficient 
   const int offset = j*deg1+k;  // position of the k-th coefficient

   __shared__ double xvre[d_shmemsize];
   __shared__ double xvim[d_shmemsize];
   __shared__ double yvre[d_shmemsize];
   __shared__ double yvim[d_shmemsize];
   __shared__ double zvre[d_shmemsize];
   __shared__ double zvim[d_shmemsize];

   double xr,xi,yr,yi,zr,zi;

   int idx = deg1+k;

   xvre[k] = xre[offset];
   xvim[k] = xim[offset];
   yvre[k] = 0.0;         // padded with zeros
   yvim[k] = 0.0;
   yvre[idx] = yre[offset];
   yvim[idx] = yim[offset];

   xr = xvre[0];   xi = xvim[0];    // zv[k] = xv[0]*yv[idx];
   yr = yvre[idx]; yi = yvim[idx];
   zr = xr*yr - xi*yi;
   zi = xr*yi + xi*yr;
   zvre[k] = zr;
   zvim[k] = zi;

   for(int i=1; i<deg1; i++)
   {
      idx = deg1 + k - i;
      xr = xvre[i];   xi = xvim[i];   // zv[k] = zv[k] + xv[i]*yv[idx];
      yr = yvre[idx]; yi = yvim[idx];
      zr = xr*yr - xi*yi;
      zi = xr*yi + xi*yr;
      zvre[k] += zr;
      zvim[k] += zi;
   }
   xre[offset] = zvre[k];
   xim[offset] = zvim[k];
}

__global__ void dbl_additions ( double *x, int lag, int shf, int deg1 )
{
   const int j = blockIdx.x; 
   const int k = threadIdx.x;
   const int x_offset = shf+j*deg1+k;
   int y_offset;

   __shared__ double xv[d_shmemsize];
   __shared__ double yv[d_shmemsize];

   xv[k] = x[x_offset];

   if(j < lag)
   {
      y_offset = x_offset + lag*deg1;
      yv[k] = x[y_offset];
      xv[k] += yv[k];
      x[x_offset] = xv[k]; // store for next block in next round
   }
}

__global__ void cmplx_additions
 ( double *xre, double *xim, int lag, int shf, int deg1 )
{
   const int j = blockIdx.x; 
   const int k = threadIdx.x;
   const int x_offset = shf+j*deg1+k;
   int y_offset;

   __shared__ double xvre[d_shmemsize];
   __shared__ double xvim[d_shmemsize];
   __shared__ double yvre[d_shmemsize];
   __shared__ double yvim[d_shmemsize];

   xvre[k] = xre[x_offset];
   xvim[k] = xim[x_offset];

   if(j < lag)
   {
      y_offset = x_offset + lag*deg1;
      yvre[k] = xre[y_offset];
      yvim[k] = xim[y_offset];
      xvre[k] += yvre[k];
      xvim[k] += yvim[k];
      xre[x_offset] = xvre[k]; // store for next block in next round
      xim[x_offset] = xvim[k];
   }
}

void GPU_dbl_inner_product
 ( int BS, int dim, int deg, double **x, double **y, double *z,
   int mode, bool verbose )
{
   const int deg1 = deg+1;         // coefficient series length

   double* x_d;                    // x_d is x_h on the device
   double* y_d;                    // y_d is y_h on the device
   double* z_d;                    // z_d is z_h on the device

   size_t szdeg = deg1*sizeof(double);
   size_t szdim = dim*szdeg;

   double* x_h = new double[dim*deg1];
   double* y_h = new double[dim*deg1];
   int ix = 0;
   for(int i=0; i<dim; i++)
      for(int j=0; j<deg1; j++)
      {
         x_h[ix] = x[i][j]; y_h[ix++] = y[i][j];
      }

   cudaMalloc((void**)&x_d,szdim);
   cudaMalloc((void**)&y_d,szdim);
   cudaMalloc((void**)&z_d,szdeg);
   cudaMemcpy(x_d,x_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(y_d,y_h,szdim,cudaMemcpyHostToDevice);

   if(BS == deg1) dbl_convolutions<<<dim,BS>>>(x_d,y_d,deg1);

   if(mode==1) // do all additions on the host
   {
      cudaMemcpy(x_h,x_d,szdim,cudaMemcpyDeviceToHost);

      for(int i=0; i<deg1; i++) z[i] = 0.0;

      int ix=0;
      for(int j=0; j<dim; j++)
         for(int i=0; i<deg1; i++) z[i] = z[i] + x_h[ix++];
   }
   else
   {
      using namespace std;

      double logdim = log2((double)dim);
      double ceil_logdim = ceil(logdim);
      int ceil2log = int(ceil_logdim);

      if(verbose)
      {
         cout << "log(" << dim << ") : " << logdim << endl;
         cout << "ceil(log(" << dim << ")) : " << ceil_logdim << endl;
         cout << "ceil2log : " << ceil2log << endl;
      }
      if(BS == deg1)
      {
         int restshift = (dim % 2 == 0 ? 0 : deg1);
         int lag = dim/2;

         for(int L=0; L<ceil2log; L++)
         {
            if(L == ceil2log-1) restshift = 0; // no shift at end

            if(verbose)
               cout << "restshift : " << restshift
                    << "  lag : " << lag << "  L : " << L << endl;

            dbl_additions<<<lag,BS>>>(x_d,lag,restshift,deg1);

            if(restshift == 0)                // no shift left
            {
               restshift = (lag % 2 == 0 ? 0 : deg1);
               lag = (lag == 1 ? 1 : lag/2);
            }
            else                              // have shift left
            {
               if(lag % 2 == 0)               // if even lag,
                  lag = lag/2;                // then keep the shift
               else
               {
                  restshift = 0;              // set shift to zero
                  lag = (lag == 1 ? 1 : lag/2);
                  if(L < ceil2log-2) lag = lag + 1;
               }
            }
         }
      }
      cudaMemcpy(z,x_d,szdeg,cudaMemcpyDeviceToHost);
   }
   free(x_h); free(y_h);
}

void GPU_cmplx_inner_product
 ( int BS, int dim, int deg,
   double **xre, double **xim, double **yre, double **yim,
   double *zre, double *zim, int mode, bool verbose )
{
   const int deg1 = deg+1;         // coefficient series length

   double* xre_d;                  // xre_d is xre_h on the device
   double* xim_d;                  // xim_d is xim_h on the device
   double* yre_d;                  // yre_d is yre_h on the device
   double* yim_d;                  // yim_d is yim_h on the device
   double* zre_d;                  // zre_d is zre_h on the device
   double* zim_d;                  // zim_d is zim_h on the device

   size_t szdeg = deg1*sizeof(double);
   size_t szdim = dim*szdeg;

   double* xre_h = new double[dim*deg1];
   double* xim_h = new double[dim*deg1];
   double* yre_h = new double[dim*deg1];
   double* yim_h = new double[dim*deg1];
   int ix = 0;
   for(int i=0; i<dim; i++)
      for(int j=0; j<deg1; j++)
      {
         xre_h[ix] = xre[i][j];
         xim_h[ix] = xim[i][j];
         yre_h[ix] = yre[i][j];
         yim_h[ix++] = yim[i][j];
      }

   cudaMalloc((void**)&xre_d,szdim);
   cudaMalloc((void**)&xim_d,szdim);
   cudaMalloc((void**)&yre_d,szdim);
   cudaMalloc((void**)&yim_d,szdim);
   cudaMalloc((void**)&zre_d,szdeg);
   cudaMalloc((void**)&zim_d,szdeg);
   cudaMemcpy(xre_d,xre_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(xim_d,xim_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(yre_d,yre_h,szdim,cudaMemcpyHostToDevice);
   cudaMemcpy(yim_d,yim_h,szdim,cudaMemcpyHostToDevice);

   if(BS == deg1)
      cmplx_convolutions<<<dim,BS>>>(xre_d,xim_d,yre_d,yim_d,deg1);

   if(mode==1) // do all additions on the host
   {
      cudaMemcpy(xre_h,xre_d,szdim,cudaMemcpyDeviceToHost);
      cudaMemcpy(xim_h,xim_d,szdim,cudaMemcpyDeviceToHost);

      for(int i=0; i<deg1; i++)
      {
         zre[i] = 0.0;
         zim[i] = 0.0;
      }
      int ix=0;
      for(int j=0; j<dim; j++)
         for(int i=0; i<deg1; i++)
         {
            zre[i] = zre[i] + xre_h[ix];
            zim[i] = zim[i] + xim_h[ix++];
         }
   }
   else
   {
      using namespace std;

      double logdim = log2((double)dim);
      double ceil_logdim = ceil(logdim);
      int ceil2log = int(ceil_logdim);

      if(verbose)
      {
         cout << "log(" << dim << ") : " << logdim << endl;
         cout << "ceil(log(" << dim << ")) : " << ceil_logdim << endl;
         cout << "ceil2log : " << ceil2log << endl;
      }
      if(BS == deg1)
      {
         int restshift = (dim % 2 == 0 ? 0 : deg1);
         int lag = dim/2;

         for(int L=0; L<ceil2log; L++)
         {
            if(L == ceil2log-1) restshift = 0; // no shift at end

            if(verbose)
               cout << "restshift : " << restshift
                    << "  lag : " << lag << "  L : " << L << endl;

            cmplx_additions<<<lag,BS>>>(xre_d,xim_d,lag,restshift,deg1);

            if(restshift == 0)                // no shift left
            {
               restshift = (lag % 2 == 0 ? 0 : deg1);
               lag = (lag == 1 ? 1 : lag/2);
            }
            else                              // have shift left
            {
               if(lag % 2 == 0)               // if even lag,
                  lag = lag/2;                // then keep the shift
               else
               {
                  restshift = 0;              // set shift to zero
                  lag = (lag == 1 ? 1 : lag/2);
                  if(L < ceil2log-2) lag = lag + 1;
               }
            }
         }
      }
      cudaMemcpy(zre,xre_d,szdeg,cudaMemcpyDeviceToHost);
      cudaMemcpy(zim,xim_d,szdeg,cudaMemcpyDeviceToHost);
   }
   free(xre_h); free(yre_h);
   free(xim_h); free(yim_h);
}

void GPU_dbl_matrix_vector_product
 ( int BS, int rows, int cols, int deg, double ***A, double **x,
   double **y, int mode, bool verbose )
{
   const int deg1 = deg+1;         // coefficient series length

   double* row_d;                  // row_d is a row of A
   double* x_d;                    // x_d is x_h on the device

   size_t szdeg = deg1*sizeof(double);
   size_t szcol = cols*szdeg;

   double* x_h = new double[cols*deg1];
   int ix = 0;
   for(int i=0; i<cols; i++)
      for(int j=0; j<deg1; j++) x_h[ix++] = x[i][j];

   cudaMalloc((void**)&x_d,szcol);
   cudaMemcpy(x_d,x_h,szcol,cudaMemcpyHostToDevice);

   cudaMalloc((void**)&row_d,szcol);
   double* row_h = new double[cols*deg1];

   for(int i=0; i<rows; i++)
   {
      using namespace std;

      if(verbose) cout << "Computing row " << i << " ..." << endl;

      ix = 0;
      for(int j=0; j<cols; j++)
         for(int k=0; k<deg1; k++) row_h[ix++] = A[i][j][k];

      cudaMemcpy(row_d,row_h,szcol,cudaMemcpyHostToDevice);

      if(verbose) cout << "  ... launching the kernel ..." << endl;

      if(BS == deg1) dbl_convolutions<<<cols,BS>>>(row_d,x_d,deg1);

      if(mode == 1) // do all additions on the host
      {
         cudaMemcpy(row_h,row_d,szcol,cudaMemcpyDeviceToHost);

         for(int k=0; k<deg1; k++) y[i][k] = 0.0;

         ix=0;
         for(int j=0; j<cols; j++)
            for(int k=0; k<deg1; k++) y[i][k] = y[i][k] + row_h[ix++];
      }
      else          // launch log2(cols) summation kernels
      {
         using namespace std;

         int dim = cols;
         double logdim = log2((double)dim);
         double ceil_logdim = ceil(logdim);
         int ceil2log = int(ceil_logdim);

         if(verbose)
         {
            cout << "log(" << dim << ") : " << logdim << endl;
            cout << "ceil(log(" << dim << ")) : " << ceil_logdim << endl;
            cout << "ceil2log : " << ceil2log << endl;
         }
         if(BS == deg1)
         {
            int restshift = (dim % 2 == 0 ? 0 : deg1);
            int lag = dim/2;

            for(int L=0; L<ceil2log; L++)
            {
               if(L == ceil2log-1) restshift = 0; // no shift at end

               if(verbose)
                  cout << "restshift : " << restshift
                       << "  lag : " << lag << "  L : " << L << endl;

               dbl_additions<<<lag,BS>>>(row_d,lag,restshift,deg1);

               if(restshift == 0)                // no shift left
               {
                  restshift = (lag % 2 == 0 ? 0 : deg1);
                  lag = (lag == 1 ? 1 : lag/2);
               }
               else                              // have shift left
               {
                  if(lag % 2 == 0)               // if even lag,
                     lag = lag/2;                // then keep the shift
                  else
                  {
                     restshift = 0;              // set shift to zero
                     lag = (lag == 1 ? 1 : lag/2);
                     if(L < ceil2log-2) lag = lag + 1;
                  }
               }
            }
            cudaMemcpy(y[i],row_d,szdeg,cudaMemcpyDeviceToHost);
         }
      }
   }
}

void GPU_cmplx_matrix_vector_product
 ( int BS, int rows, int cols, int deg, double ***Are, double ***Aim,
   double **xre, double **xim, double **yre, double **yim,
   int mode, bool verbose )
{
   const int deg1 = deg+1;         // coefficient series length

   double* rowre_d;                // row_d is a row of Are
   double* rowim_d;                // row_d is a row of Aim
   double* xre_d;                  // xre_d is xre_h on the device
   double* xim_d;                  // xim_d is xim_h on the device
   // double* y_d;                    // y_d is y_h on the device

   size_t szdeg = deg1*sizeof(double);
   size_t szcol = cols*szdeg;

   double* xre_h = new double[cols*deg1];
   double* xim_h = new double[cols*deg1];
   int ix = 0;
   for(int i=0; i<cols; i++)
      for(int j=0; j<deg1; j++)
      {
         xre_h[ix] = xre[i][j];
         xim_h[ix++] = xim[i][j];
      }

   cudaMalloc((void**)&xre_d,szcol);
   cudaMalloc((void**)&xim_d,szcol);
   cudaMemcpy(xre_d,xre_h,szcol,cudaMemcpyHostToDevice);
   cudaMemcpy(xim_d,xim_h,szcol,cudaMemcpyHostToDevice);

   cudaMalloc((void**)&rowre_d,szcol);
   cudaMalloc((void**)&rowim_d,szcol);
   double* rowre_h = new double[cols*deg1];
   double* rowim_h = new double[cols*deg1];

   for(int i=0; i<rows; i++)
   {
      using namespace std;

      if(verbose) cout << "Computing row " << i << " ..." << endl;

      ix = 0;
      for(int j=0; j<cols; j++)
         for(int k=0; k<deg1; k++)
         {
            rowre_h[ix] = Are[i][j][k];
            rowim_h[ix++] = Aim[i][j][k];
         }

      cudaMemcpy(rowre_d,rowre_h,szcol,cudaMemcpyHostToDevice);
      cudaMemcpy(rowim_d,rowim_h,szcol,cudaMemcpyHostToDevice);

      if(verbose) cout << "  ... launching the kernel ..." << endl;

      if(BS == deg1)
         cmplx_convolutions<<<cols,BS>>>(rowre_d,rowim_d,xre_d,xim_d,deg1);

      cudaMemcpy(rowre_h,rowre_d,szcol,cudaMemcpyDeviceToHost);
      cudaMemcpy(rowim_h,rowim_d,szcol,cudaMemcpyDeviceToHost);

      if(mode == 1) // do all additions on the host
      {
         for(int k=0; k<deg1; k++)
         {
            yre[i][k] = 0.0;
            yim[i][k] = 0.0;
         }
         ix=0;
         for(int j=0; j<cols; j++)
            for(int k=0; k<deg1; k++)
            {
               yre[i][k] = yre[i][k] + rowre_h[ix];
               yim[i][k] = yim[i][k] + rowim_h[ix++];
            }
      }
      else          // launch log2(cols) summation kernels
      {
         using namespace std;

         int dim = cols;
         double logdim = log2((double)dim);
         double ceil_logdim = ceil(logdim);
         int ceil2log = int(ceil_logdim);

         if(verbose)
         {
            cout << "log(" << dim << ") : " << logdim << endl;
            cout << "ceil(log(" << dim << ")) : " << ceil_logdim << endl;
            cout << "ceil2log : " << ceil2log << endl;
         }
         if(BS == deg1)
         {
            int restshift = (dim % 2 == 0 ? 0 : deg1);
            int lag = dim/2;

            for(int L=0; L<ceil2log; L++)
            {
               if(L == ceil2log-1) restshift = 0; // no shift at end

               if(verbose)
                  cout << "restshift : " << restshift
                       << "  lag : " << lag << "  L : " << L << endl;

               cmplx_additions<<<lag,BS>>>
                  (rowre_d,rowim_d,lag,restshift,deg1);

               if(restshift == 0)                // no shift left
               {
                  restshift = (lag % 2 == 0 ? 0 : deg1);
                  lag = (lag == 1 ? 1 : lag/2);
               }
               else                              // have shift left
               {
                  if(lag % 2 == 0)               // if even lag,
                     lag = lag/2;                // then keep the shift
                  else
                  {
                     restshift = 0;              // set shift to zero
                     lag = (lag == 1 ? 1 : lag/2);
                     if(L < ceil2log-2) lag = lag + 1;
                  }
               }
            }
            cudaMemcpy(yre[i],rowre_d,szdeg,cudaMemcpyDeviceToHost);
            cudaMemcpy(yim[i],rowim_d,szdeg,cudaMemcpyDeviceToHost);
         }
      }
   }
}
