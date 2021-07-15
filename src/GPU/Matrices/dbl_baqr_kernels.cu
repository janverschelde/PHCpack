/* The file dbl_baqr_kernels.cu defines the functions with prototypes in
 * the file dbl_baqr_kernels.h. */

#include <iostream>
#ifdef winwalltime
#include "wingettimeofday.h"
#else
#include <sys/time.h>
#endif
#include "dbl_baqr_kernels.h"

__global__ void dbl_small_house
 ( double x0, double *x1, int dim, int dimLog2, double *v, double *beta )
{
   int j = threadIdx.x;

   __shared__ double shv[d_shmemsize];
   __shared__ double prd[d_shmemsize];

   bool stopflag = false;
   double mu,v0,v0p2;

   shv[j] = x1[j];              // reading of vector into shared memory
   prd[j] = shv[j]*shv[j];

   v[j+1] = shv[j];             // copies x to v, for a zero beta
   if(j == 0) v[0] = 1.0;

   __syncthreads();
   int powTwo = 1;                          // sum reduction
   for(int k=0; k < dimLog2; k++)
   {
      if((j%(powTwo*2)) == 0)
         if(j+powTwo < dim) prd[j] = prd[j] + prd[j+powTwo];
      powTwo = powTwo*2;
      __syncthreads();
   }
   // thread 0 computes the sqrt of the inner product, others wait
   if(j == 0)
   {
      if(prd[0] == 0.0)                    // prd[0] is sigma of house
      {
         *beta = 0.0; stopflag = true;
      }
   }
   __syncthreads();
   if(stopflag) return;                    // case when sigma is zero
   if(j == 0)
   {
      mu = sqrt(x0*x0 + prd[0]);
      if(x0 <= 0.0)
         v0 = x0 - mu;
      else
         v0 = -prd[0]/(x0 + mu);

      v0p2 = v0*v0;
      *beta = 2.0*v0p2/(prd[0] + v0p2);
   }
   __syncthreads();
   shv[j] = shv[j]/v0;
   v[j+1] = shv[j];
   if(j == 0) v[0] = 1.0;
}

void GPU_dbl_blocked_qr
 ( int nrows, int ncols, int szt, int nbt,
   double **A, double **Q, double **R,
   double *lapms, double *walltimesec )
{
   const int dim = nrows*ncols;   // total number of doubles
   double *A_h = new double[dim]; // matrix A on the host
   double *A_d;                   // matrix on the device

   int ix=0;                      // copy the columns of A to A_h
   for(int j=0; j<ncols; j++)   
      for(int i=0; i<nrows; i++) A_h[ix++] = A[i][j];

   const size_t sznum = dim*sizeof(double);
   cudaMalloc((void**)&A_d,sznum);
   cudaMemcpy(A_d,A,sznum,cudaMemcpyHostToDevice);

   free(A_h);
}
