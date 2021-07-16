/* The file dbl_baqr_kernels.cu defines the functions with prototypes in
 * the file dbl_baqr_kernels.h. */

#include <iostream>
#include <iomanip>
#ifdef winwalltime
#include "wingettimeofday.h"
#else
#include <sys/time.h>
#endif
#include "dbl_baqr_kernels.h"

using namespace std;

__global__ void dbl_small_house
 ( double *x0, double *x1, int dim, int dimLog2, double *v, double *beta )
{
   int j = threadIdx.x;

   __shared__ double shv[d_shmemsize];
   __shared__ double prd[d_shmemsize];

   bool stopflag = false;
   double mu,v0,v0p2;

   shv[j] = x1[j];              // reading of vector into shared memory
   prd[j] = shv[j]*shv[j];      // for the 2-norm computation

   v[j+1] = shv[j];             // copies x to v, in case beta is zero
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
   if(j == 0)                              // thread zero sets beta
   {
      mu = sqrt((*x0)*(*x0) + prd[0]);
      if(*x0 <= 0.0)
         v0 = *x0 - mu;
      else
         v0 = -prd[0]/(*x0 + mu);

      v0p2 = v0*v0;
      *beta = 2.0*v0p2/(prd[0] + v0p2);
      prd[0] = v0;                         // v0 needed for normalization
   }
   __syncthreads();
   shv[j] = shv[j]/prd[0];
   v[j+1] = shv[j];
   if(j == 0) v[0] = 1.0;
}

void GPU_dbl_blocked_houseqr
 ( int nrows, int ncols, int szt, int nbt,
   double **A, double **Q, double **R,
   double *lapms, double *walltimesec, bool verbose )
{
   const int dim = nrows*ncols;     // total number of doubles
   double *A_h = new double[dim];   // matrix A on the host
   double *A_d;                     // matrix on the device
   double *v_h = new double[nrows]; // Householder vector on host
   double *v_d;                     // Householder vector on device
   double *x0_d;                    // first element for house on device
   double beta_h;                   // beta on the host
   double *beta_d;                  // beta on the device

   int ix=0;                        // copy the columns of A to A_h
   for(int j=0; j<ncols; j++)   
      for(int i=0; i<nrows; i++) A_h[ix++] = A[i][j];

   const size_t sznum = dim*sizeof(double);
   cudaMalloc((void**)&A_d,sznum);
   cudaMemcpy(A_d,A_h,sznum,cudaMemcpyHostToDevice);
   const size_t szhouse = nrows*sizeof(double);
   cudaMalloc((void**)&v_d,szhouse);
   cudaMalloc((void**)&x0_d,sizeof(double));
   cudaMemcpy(x0_d,&A_h[0],sizeof(double),cudaMemcpyHostToDevice);
   cudaMalloc((void**)&beta_d,sizeof(double));

   const int nrows1 = nrows-1;
   const int nrLog2 = ceil(log2((double) nrows1));

   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   *lapms = 0.0;
   float milliseconds;
   struct timeval begintime,endtime; // wall clock time of computations

   gettimeofday(&begintime,0);

   cudaEventRecord(start);
   dbl_small_house<<<1,nrows1>>>(x0_d,&A_d[1],nrows1,nrLog2,v_d,beta_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;

   gettimeofday(&endtime,0);
   long seconds = endtime.tv_sec - begintime.tv_sec;
   long microseconds = endtime.tv_usec - begintime.tv_usec;
   *walltimesec = seconds + microseconds*1.0e-6;

   cudaMemcpy(&beta_h,beta_d,sizeof(double),cudaMemcpyDeviceToHost);

   if(verbose)
      cout << scientific << setprecision(16)
           << "The first beta : " << beta_h << endl;

   free(A_h); free(v_h);
}
