/* The file dbl_tabs_kernels.cu defines the functions with prototypes in
 * the file dbl_tabs_kernels.h. */

#include <iostream>
#include "dbl_tabs_kernels.h"

using namespace std;

__global__ void dbl_small_invert_upper ( int dim, double *U, double *invU )
{
   const int k = threadIdx.x; // thread k computes k-th column of inverse

   __shared__ double Ucol[d_shmemsize];
   __shared__ double invUrows[d_shmemsize];

   double rhs,xval;

   int colidx = dim*(dim-1);          // start with the last column

   Ucol[k] = U[colidx+k];             // load the last column
   rhs = ((double) int(k == dim-1));  // right hand side for each thread
   int rowidx = (dim - 1)*dim + k;    // the row index in the inverse

   __syncthreads();
   invUrows[rowidx] = rhs/Ucol[k];    // last row of the inverse

   for(int i=dim-2; i>=0; i--)        // compute row with index i
   {
      rhs = ((double) int(k == i));   // set rhs for i-th unit vector

      for(int j=i+1; j<dim; j++)
      {
         colidx = dim*j;              // need column j of U

         Ucol[k] = U[colidx+k];

         rowidx = j*dim + k;          // need solution value

         xval = invUrows[rowidx];

         __syncthreads();
         rhs = rhs - Ucol[i]*xval;    // update right hand side
      }
      rowidx = i*dim + k;             // save in i-th row of inverse

      colidx = dim*i;                 // need column i of U
      Ucol[k] = U[colidx+k];

      __syncthreads();
      invUrows[rowidx] = rhs/Ucol[i];
   }
   rowidx = 0;
   for(int i=0; i<dim; i++)
   {
      __syncthreads();
      invU[rowidx+k] = invUrows[rowidx+k];
      rowidx = rowidx + dim;
   }
}

void test_dbl_small_invert_upper ( int dim, double *U, double *invU )
{
   double *Ucol = new double[d_shmemsize];
   double *invUrows = new double[d_shmemsize];

   for(int i=0; i<dim*dim; i++)
      cout << "U[" << i << "] : " << U[i] << endl;

   for(int k=0; k<dim; k++)  // compute k-th column of inverse
   {
      double rhs,xval;

      cout << "******** step " << k << " *********" << endl;

      int colidx = dim*(dim-1);          // start with the last column

      Ucol[k] = U[colidx+k];             // load the last column
      rhs = ((double) int(k == dim-1));  // right hand side for each thread
      int rowidx = (dim - 1)*dim + k;    // the row index in the inverse

      cout << "  k : " << k
           << "  i : " << dim-1
           << "  rhs : " << rhs << endl;
      cout << "  k : " << k
           << "  i : " << dim-1
           << "  assigning to rowidx : " << rowidx << endl;
      invUrows[rowidx] = rhs/Ucol[k];    // last row of the inverse
      cout << "  value : " << invUrows[rowidx] << endl;
    
      cout << "            entering loop" << endl;

      for(int i=dim-2; i>=0; i--)        // compute row with index i
      {
         rhs = ((double) int(k == i));   // set rhs for i-th unit vector
 
         cout << "  k : " << k
              << "  i : " << i
              << "  rhs : " << rhs << endl;

         for(int j=i+1; j<dim; j++)
         {
            colidx = dim*j;              // need column j of U
            cout << "i : " << i << "  j : " << j << "  k : " << k
                 << "  colidx : " << colidx << endl;
            // Ucol[k] = U[colidx+k];
            for(int L=0; L<dim; L++) Ucol[L] = U[colidx+L];
            rowidx = j*dim + k;          // need solution value
            xval = invUrows[rowidx];
            cout << "Ucol[" << i << "] : " << Ucol[i] << endl;
            cout << "xval : " << xval << endl;
            rhs = rhs - Ucol[i]*xval;    // update right hand side
         }
         rowidx = i*dim + k;             // save in i-th row of inverse

         colidx = dim*i;                 // need column i of U
         // Ucol[k] = U[colidx+k];
         for(int L=0; L<dim; L++) Ucol[L] = U[colidx+L];
         cout << "  k : " << k
              << "  i : " << i
              << "  assigning to rowidx : " << rowidx << endl;
         invUrows[rowidx] = rhs/Ucol[i];
         cout << "  value : " << invUrows[rowidx] << endl;
      }
      rowidx = 0;
      for(int i=0; i<dim; i++)
      {
         invU[rowidx+k] = invUrows[rowidx+k];
         rowidx = rowidx + dim;
      }
   }
   for(int i=0; i<dim*dim; i++)
      cout << "invUrows[" << i << "] : " << invUrows[i] << endl;
}

__global__ void dbl_medium_invert_upper ( int dim, double *U, double *invU )
{
   const int k = threadIdx.x;  // thread k computes k-th column of inverse

   __shared__ double Ucol[d_shmemsize];      // one column of U
   __shared__ double invUrow[d_shmemsize];   // one row of invU

   double rhs,xval;

   int colidx = dim*(dim-1);          // start with the last column

   Ucol[k] = U[colidx+k];             // load the last column
   rhs = ((double) int(k == dim-1));  // right hand side for each thread
   int rowidx = (dim - 1)*dim + k;    // the row index in the inverse

   invUrow[k] = rhs/Ucol[k];          // last row of the inverse
   invU[rowidx] = invUrow[k];         // store the last row into invU

   for(int i=dim-2; i>=0; i--)        // compute row with index i
   {
      rhs = ((double) int(k == i));   // set rhs for i-th unit vector

      for(int j=i+1; j<dim; j++)
      {
         colidx = dim*j;              // need column j of U
         Ucol[k] = U[colidx+k];

         rowidx = j*dim + k;          // need solution value
         invUrow[k] = invU[rowidx];   // load invU row into invUrow
         xval = invUrow[k];

         __syncthreads();
         rhs = rhs - Ucol[i]*xval;    // update right hand side
      }
      colidx = dim*i;                 // need column i of U
      Ucol[k] = U[colidx+k];
      rowidx = i*dim + k;             // save in i-th row of inverse

      __syncthreads();
      invUrow[k] = rhs/Ucol[i];
      invU[rowidx] = invUrow[k];
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

   // test_small_invert_upper(dim,U_h,invU_h); // only for debugging

   size_t szmat = szU*sizeof(double);
   cudaMalloc((void**)&U_d,szmat);
   cudaMalloc((void**)&invU_d,szmat);
   cudaMemcpy(U_d,U_h,szmat,cudaMemcpyHostToDevice);

   if(dim <= 16)
      dbl_small_invert_upper<<<1,dim>>>(dim,U_d,invU_d);
   else
      dbl_medium_invert_upper<<<1,dim>>>(dim,U_d,invU_d);

   cudaMemcpy(invU_h,invU_d,szmat,cudaMemcpyDeviceToHost);

   ix = 0;
   for(int i=0; i<dim; i++)
      for(int j=0; j<dim; j++) invU[i][j] = invU_h[ix++];

   free(U_h); free(invU_h);
}
