/* The file dbl_tabs_kernels.cu defines the functions with prototypes in
 * the file dbl_tabs_kernels.h. */

#include "dbl_tabs_kernels.h"

__global__ void dbl_invert_upper ( int dim, double *U, double *invU )
{
   const int k = threadIdx.x; // thread k computes k-th column of inverse

   __shared__ double Ucol[d_shmemsize];
   __shared__ double invUrows[d_shmemsize];

   double rhs,xval;

   int colidx = dim*(dim-1);          // start with the last column

   Ucol[k] = U[colidx+k];             // load the last column

   rhs = ((double) int(k == dim-1));  // right hand side for each thread

   int rowidx = dim - 1 + k*dim;      // the row index in the inverse

   invUrows[rowidx] = rhs/Ucol[k];    // last row of the inverse

   for(int i=dim-2; i>=0; i--)        // compute row with index i
   {
      rhs = ((double) int(k == i));   // set rhs for i-th unit vector

      for(int j=i+1; j<dim; j++)
      {
         colidx = dim*j;              // need column j of U
         Ucol[k] = U[colidx+k];
         rowidx = j + k*dim;          // need solution value
         xval = invUrows[rowidx];
         __syncthreads();
         rhs = rhs - Ucol[j]*xval;    // update right hand side
      }
      rowidx = i + k*dim;             // save in i-th row of inverse
      colidx = dim*i;
      Ucol[k] = U[colidx+k];
      __syncthreads();
      invUrows[rowidx] = rhs/Ucol[i];
   }
   rowidx = 0;
   for(int i=0; i<dim; i++)
   {
      invU[rowidx+k] = invUrows[rowidx+k];
      rowidx = rowidx + dim;
   }
}

void GPU_dbl_upper_inverse ( int dim, double **U, double **invU )
{
   const int szU = dim*dim;

   double *U_h = new double[szU];     // U_h stores the columns of U 
   double *U_d;                       // U_d is U_h on the device
   double *invU_h = new double[szU];  // the columns of the inverse
   double *invU_d;                    // invU_d is invU_h on the device

   int ix = 0;
   for(int j=0; j<dim; j++)
      for(int i=0; i<dim; i++) U_h[ix++] = U[i][j];

   size_t szmat = szU*sizeof(double);
   cudaMalloc((void**)&U_d,szmat);
   cudaMalloc((void**)&invU_d,szmat);
   cudaMemcpy(U_d,U_h,szmat,cudaMemcpyHostToDevice);

   dbl_invert_upper<<<1,dim>>>(dim,U_d,invU_d);

   cudaMemcpy(invU_h,invU_d,szmat,cudaMemcpyDeviceToHost);

   ix = 0;
   for(int i=0; i<dim; i++)
      for(int j=0; j<dim; j++) invU[i][j] = invU_h[ix++];

   free(U_h); free(invU_h);
}
