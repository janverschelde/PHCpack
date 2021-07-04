/* The file dbl2_tabs_kernels.cu defines the functions specified in
 * the file dbl2_tabs_kernels.h. */

#include <iostream>
#ifdef gpufun
#include "double_double_gpufun.cu"
#endif
#include "dbl2_tabs_kernels.h"

using namespace std;

__global__ void dbl2_small_invert_upper 
( int dim, double *Uhi, double *Ulo, double *invUhi, double *invUlo )
{
   const int k = threadIdx.x; // thread k computes k-th column of inverse

   __shared__ double Ucolhi[dd_shmemsize];
   __shared__ double Ucollo[dd_shmemsize];
   __shared__ double invUrowshi[dd_shmemsize];
   __shared__ double invUrowslo[dd_shmemsize];

   double rhshi,rhslo,xvalhi,xvallo,acchi,acclo;

   int colidx = dim*(dim-1);          // start with the last column

   Ucolhi[k] = Uhi[colidx+k];         // load the last column
   Ucollo[k] = Ulo[colidx+k];
   rhshi = ((double) int(k == dim-1));  // right hand side for each thread
   rhslo = 0.0;
   int rowidx = (dim - 1)*dim + k;      // the row index in the inverse

   __syncthreads();
   // invUrows[rowidx] = rhs/Ucol[k]; // last row of the inverse
   ddg_div(rhshi,rhslo,Ucolhi[k],Ucollo[k],
           &invUrowshi[rowidx],&invUrowslo[rowidx]);

   for(int i=dim-2; i>=0; i--)        // compute row with index i
   {
      rhshi = ((double) int(k == i)); // set rhs for i-th unit vector
      rhslo = 0.0;

      for(int j=i+1; j<dim; j++)
      {
         colidx = dim*j;              // need column j of U

         Ucolhi[k] = Uhi[colidx+k];
         Ucollo[k] = Ulo[colidx+k];

         rowidx = j*dim + k;          // need solution value

         xvalhi = invUrowshi[rowidx];
         xvallo = invUrowslo[rowidx];

         __syncthreads();
         // rhs = rhs - Ucol[i]*xval; // update right hand side
         ddg_mul(Ucolhi[i],Ucollo[i],xvalhi,xvallo,&acchi,&acclo);
         ddg_dec(&rhshi,&rhslo,acchi,acclo);
      }
      rowidx = i*dim + k;             // save in i-th row of inverse

      colidx = dim*i;                 // need column i of U
      Ucolhi[k] = Uhi[colidx+k];
      Ucollo[k] = Ulo[colidx+k];

      __syncthreads();
      // invUrows[rowidx] = rhs/Ucol[i];
      ddg_div(rhshi,rhslo,Ucolhi[i],Ucollo[i],
              &invUrowshi[rowidx],&invUrowslo[rowidx]);
   }
   rowidx = 0;
   for(int i=0; i<dim; i++)
   {
      __syncthreads();
      invUhi[rowidx+k] = invUrowshi[rowidx+k];
      invUlo[rowidx+k] = invUrowslo[rowidx+k];
      rowidx = rowidx + dim;
   }
}

__global__ void dbl2_medium_invert_upper
 ( int dim, double *Uhi, double *Ulo, double *invUhi, double *invUlo)
{
   const int k = threadIdx.x;  // thread k computes k-th column of inverse

   __shared__ double Ucolhi[dd_shmemsize];      // one column of U
   __shared__ double Ucollo[dd_shmemsize];      // one column of U
   __shared__ double invUrowhi[dd_shmemsize];   // one row of invU
   __shared__ double invUrowlo[dd_shmemsize];   // one row of invU

   double rhshi,rhslo,xvalhi,xvallo,acchi,acclo;

   int colidx = dim*(dim-1);           // start with the last column

   Ucolhi[k] = Uhi[colidx+k];          // load the last column
   Ucollo[k] = Ulo[colidx+k];
   rhshi = ((double) int(k == dim-1)); // right hand side for each thread
   rhslo = 0.0;
   int rowidx = (dim - 1)*dim + k;     // the row index in the inverse

   // invUrow[k] = rhs/Ucol[k];          // last row of the inverse
   ddg_div(rhshi,rhslo,Ucolhi[k],Ucollo[k],&invUrowhi[k],&invUrowlo[k]);
   invUhi[rowidx] = invUrowhi[k];     // store the last row into invU
   invUlo[rowidx] = invUrowlo[k]; 

   for(int i=dim-2; i>=0; i--)        // compute row with index i
   {
      rhshi = ((double) int(k == i)); // set rhs for i-th unit vector
      rhslo = 0.0;

      for(int j=i+1; j<dim; j++)
      {
         colidx = dim*j;              // need column j of U
         Ucolhi[k] = Uhi[colidx+k];
         Ucollo[k] = Ulo[colidx+k];

         rowidx = j*dim + k;            // need solution value
         invUrowhi[k] = invUhi[rowidx]; // load invU row into invUrow
         invUrowlo[k] = invUlo[rowidx];
         xvalhi = invUrowhi[k];
         xvallo = invUrowlo[k];

         __syncthreads();
         // rhs = rhs - Ucol[i]*xval;    // update right hand side
         ddg_mul(Ucolhi[i],Ucollo[i],xvalhi,xvallo,&acchi,&acclo);
         ddg_dec(&rhshi,&rhslo,acchi,acclo);
      }
      colidx = dim*i;                 // need column i of U
      Ucolhi[k] = Uhi[colidx+k];
      Ucollo[k] = Ulo[colidx+k];
      rowidx = i*dim + k;             // save in i-th row of inverse

      __syncthreads();
      // invUrow[k] = rhs/Ucol[i];
      ddg_div(rhshi,rhslo,Ucolhi[i],Ucollo[i],&invUrowhi[k],&invUrowlo[k]);
      invUhi[rowidx] = invUrowhi[k];
      invUlo[rowidx] = invUrowlo[k];
   }
}

void GPU_dbl2_upper_inverse
 ( int dim, double **Uhi, double **Ulo, double **invUhi, double **invUlo )
{
   const int szU = dim*dim;

   double *Uhi_h = new double[szU];     // Uhi_h stores the columns of Uhi
   double *Ulo_h = new double[szU];     // Ulo_h stores the columns of Ulo 
   double *Uhi_d;                       // Uhi_d is Uhi_h on the device
   double *Ulo_d;                       // Ulo_d is Ulo_h on the device
   double *invUhi_h = new double[szU];  // high doubles of the inverse
   double *invUlo_h = new double[szU];  // low doubles of the inverse
   double *invUhi_d;                    // invUhi_d is invUhi_h on the device
   double *invUlo_d;                    // invUlo_d is invUlo_h on the device

   int ix = 0;
   for(int j=0; j<dim; j++)
      for(int i=0; i<dim; i++)
      {
         Uhi_h[ix] = Uhi[i][j];
         Ulo_h[ix++] = Ulo[i][j];
      }

   // only for debugging
   // test_dbl2_small_invert_upper(dim,Uhi_h,Ulo_h,invUhi,invUlo_h);

   size_t szmat = szU*sizeof(double);
   cudaMalloc((void**)&Uhi_d,szmat);
   cudaMalloc((void**)&Ulo_d,szmat);
   cudaMalloc((void**)&invUhi_d,szmat);
   cudaMalloc((void**)&invUlo_d,szmat);
   cudaMemcpy(Uhi_d,Uhi_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(Ulo_d,Ulo_h,szmat,cudaMemcpyHostToDevice);

   if(dim <= 16)
      dbl2_small_invert_upper<<<1,dim>>>(dim,Uhi_d,Ulo_d,invUhi_d,invUlo_d);
   else
      dbl2_medium_invert_upper<<<1,dim>>>(dim,Uhi_d,Ulo_d,invUhi_d,invUlo_d);

   cudaMemcpy(invUhi_h,invUhi_d,szmat,cudaMemcpyDeviceToHost);
   cudaMemcpy(invUlo_h,invUlo_d,szmat,cudaMemcpyDeviceToHost);

   ix = 0;
   for(int i=0; i<dim; i++)
      for(int j=0; j<dim; j++)
      {
         invUhi[i][j] = invUhi_h[ix];
         invUlo[i][j] = invUlo_h[ix++];
      }

   free(Uhi_h); free(invUhi_h);
   free(Ulo_h); free(invUlo_h);
}
