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

__global__ void cmplx_small_invert_upper
 ( int dim, double *Ure, double *Uim, double *invUre, double *invUim )
{
   const int k = threadIdx.x; // thread k computes k-th column of inverse

   __shared__ double Ucolre[d_shmemsize];
   __shared__ double Ucolim[d_shmemsize];
   __shared__ double invUrowsre[d_shmemsize];
   __shared__ double invUrowsim[d_shmemsize];

   double rhsre,rhsim,xvalre,xvalim,accre,accim,det;

   int colidx = dim*(dim-1);           // start with the last column

   Ucolre[k] = Ure[colidx+k];          // load the last column
   Ucolim[k] = Uim[colidx+k];
   rhsre = ((double) int(k == dim-1)); // right hand side for each thread
   rhsim = 0.0;
   int rowidx = (dim - 1)*dim + k;     // the row index in the inverse

   __syncthreads();
   // invUrows[rowidx] = rhs/Ucol[k];  // last row of the inverse
   det = Ucolre[k]*Ucolre[k] + Ucolim[k]*Ucolim[k];
   accre = Ucolre[k]/det;
   accim = -Ucolim[k]/det;
   invUrowsre[rowidx] = rhsre*accre - rhsim*accim;
   invUrowsim[rowidx] = rhsim*accre + rhsre*accim;

   for(int i=dim-2; i>=0; i--)        // compute row with index i
   {
      rhsre = ((double) int(k == i)); // set rhs for i-th unit vector
      rhsim = 0.0;

      for(int j=i+1; j<dim; j++)
      {
         colidx = dim*j;              // need column j of U

         Ucolre[k] = Ure[colidx+k];
         Ucolim[k] = Uim[colidx+k];

         rowidx = j*dim + k;          // need solution value

         xvalre = invUrowsre[rowidx];
         xvalim = invUrowsim[rowidx];

         __syncthreads();
         // rhs = rhs - Ucol[i]*xval; // update right hand side
         accre = Ucolre[i]*xvalre - Ucolim[i]*xvalim;
         accim = Ucolim[i]*xvalre + Ucolre[i]*xvalim;
         rhsre = rhsre - accre;
         rhsim = rhsim - accim;
      }
      rowidx = i*dim + k;             // save in i-th row of inverse

      colidx = dim*i;                 // need column i of U
      Ucolre[k] = Ure[colidx+k];
      Ucolim[k] = Uim[colidx+k];

      __syncthreads();
      // invUrows[rowidx] = rhs/Ucol[i];
      det = Ucolre[i]*Ucolre[i] + Ucolim[i]*Ucolim[i];
      accre = Ucolre[i]/det;
      accim = -Ucolim[i]/det;
      invUrowsre[rowidx] = rhsre*accre - rhsim*accim;
      invUrowsim[rowidx] = rhsim*accre + rhsre*accim;
   }
   rowidx = 0;
   for(int i=0; i<dim; i++)
   {
      __syncthreads();
      invUre[rowidx+k] = invUrowsre[rowidx+k];
      invUim[rowidx+k] = invUrowsim[rowidx+k];
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

__global__ void cmplx_medium_invert_upper
 ( int dim, double *Ure, double *Uim, double *invUre, double *invUim )
{
   const int k = threadIdx.x;  // thread k computes k-th column of inverse

   __shared__ double Ucolre[d_shmemsize];    // one column of U
   __shared__ double Ucolim[d_shmemsize];    // imaginary parts of U
   __shared__ double invUrowre[d_shmemsize]; // one row of invU
   __shared__ double invUrowim[d_shmemsize]; // imaginary parts of invU

   double rhsre,rhsim,xvalre,xvalim,det,accre,accim;

   int colidx = dim*(dim-1);            // start with the last column

   Ucolre[k] = Ure[colidx+k];           // load the last column
   Ucolim[k] = Uim[colidx+k]; 
   rhsre = ((double) int(k == dim-1));  // right hand side for each thread
   rhsim = 0.0;
   int rowidx = (dim - 1)*dim + k;      // the row index in the inverse

   // invUrow[k] = rhs/Ucol[k];         // last row of the inverse
   det = Ucolre[k]*Ucolre[k] + Ucolim[k]*Ucolim[k];
   accre = Ucolre[k]/det;
   accim = -Ucolim[k]/det;
   invUrowre[k] = rhsre*accre - rhsim*accim;
   invUrowim[k] = rhsim*accre + rhsre*accim;
   invUre[rowidx] = invUrowre[k];       // store the last row into invU
   invUim[rowidx] = invUrowim[k];

   for(int i=dim-2; i>=0; i--)          // compute row with index i
   {
      rhsre = ((double) int(k == i));   // set rhs for i-th unit vector
      rhsim = 0.0;

      for(int j=i+1; j<dim; j++)
      {
         colidx = dim*j;                // need column j of U
         Ucolre[k] = Ure[colidx+k];
         Ucolim[k] = Uim[colidx+k];

         rowidx = j*dim + k;            // need solution value
         invUrowre[k] = invUre[rowidx]; // load invU row into invUrow
         invUrowim[k] = invUim[rowidx];
         xvalre = invUrowre[k];
         xvalim = invUrowim[k];

         __syncthreads();
         // rhs = rhs - Ucol[i]*xval;   // update right hand side
         accre = Ucolre[i]*xvalre - Ucolim[i]*xvalim;
         accim = Ucolim[i]*xvalre + Ucolre[i]*xvalim;
         rhsre = rhsre - accre;
         rhsim = rhsim - accim;
      }
      colidx = dim*i;                   // need column i of U
      Ucolre[k] = Ure[colidx+k];
      Ucolim[k] = Uim[colidx+k];
      rowidx = i*dim + k;               // save in i-th row of inverse

      __syncthreads();
      // invUrow[k] = rhs/Ucol[i];
      det = Ucolre[i]*Ucolre[i] + Ucolim[i]*Ucolim[i];
      accre = Ucolre[i]/det;
      accim = -Ucolim[i]/det;
      invUrowre[k] = rhsre*accre - rhsim*accim;
      invUrowim[k] = rhsim*accre + rhsre*accim;
      invUre[rowidx] = invUrowre[k];
      invUim[rowidx] = invUrowim[k];
   }
}

__global__ void  dbl_invert_tiles ( int dim, double *U, double *invU )
{
   const int B = blockIdx.x;   // block index
   const int k = threadIdx.x;  // thread k computes k-th column of inverse
   const int offset = dim*dim*B; // offset in U and invU

   __shared__ double Ucol[d_shmemsize];      // one column of U
   __shared__ double invUrow[d_shmemsize];   // one row of invU

   double rhs,xval;

   int colidx = offset + dim*(dim-1); // start with the last column

   Ucol[k] = U[colidx+k];             // load the last column
   rhs = ((double) int(k == dim-1));  // right hand side for each thread
   int rowidx = offset + (dim - 1)*dim + k; // row index in the inverse

   invUrow[k] = rhs/Ucol[k];          // last row of the inverse
   invU[rowidx] = invUrow[k];         // store the last row into invU

   for(int i=dim-2; i>=0; i--)        // compute row with index i
   {
      rhs = ((double) int(k == i));   // set rhs for i-th unit vector

      for(int j=i+1; j<dim; j++)
      {
         colidx = offset + dim*j;     // need column j of U
         Ucol[k] = U[colidx+k];

         rowidx = offset + j*dim + k; // need solution value
         invUrow[k] = invU[rowidx];   // load invU row into invUrow
         xval = invUrow[k];

         __syncthreads();
         rhs = rhs - Ucol[i]*xval;    // update right hand side
      }
      colidx = offset + dim*i;        // need column i of U
      Ucol[k] = U[colidx+k];
      rowidx = offset + i*dim + k;    // save in i-th row of inverse

      __syncthreads();
      invUrow[k] = rhs/Ucol[i];
      invU[rowidx] = invUrow[k];
   }
}

__global__ void dbl_multiply_inverse
 ( int dim, int idx, double *invU, double *w )
{
   const int k = threadIdx.x;     // thread k computes k-th product
   const int rhsoff = dim*idx;    // offset for the right hand size
   const int offset = dim*rhsoff; // offset for diagonal tile

   __shared__ double work[d_shmemsize];      // copy of w

   work[k] = w[rhsoff+k];

   double result = 0.0; // each thread stores its product in result
   double coeff;

   for(int j=0; j<dim; j++)  // column j of the inverse diagonal tile
   {
      coeff = invU[offset+k*dim+j]; // thread k does row k
      result = result + coeff*work[j];
   }
   w[rhsoff+k] = result;
}

__global__ void dbl_back_substitute
 ( int dim, int idx, double *U, double *w )
{
   const int B = blockIdx.x;     // block index
   const int k = threadIdx.x;    // thread k computes k-th product
   const int offset = B*dim*dim; // numbers to skip

   __shared__ double wrk[d_shmemsize];    // copy of w
   __shared__ double sol[d_shmemsize];    // solution to update with

   wrk[k] = w[B*dim+k];    // block B updates B-th slice of w
   sol[k] = w[idx*dim+k];  // solution that is back substituted

   double result = 0.0; // each thread stores its product in result
   double coeff;

   for(int j=0; j<dim; j++)  // column j of the inverse diagonal tile
   {
      coeff = U[offset+k*dim+j];
      result = result + coeff*sol[j];
   }
   wrk[k] = wrk[k] - result; // subtract product
   w[B*dim+k] = wrk[k];
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

void GPU_cmplx_upper_inverse
 ( int dim, double **Ure, double **Uim, double **invUre, double **invUim )
{
   const int szU = dim*dim;

   double *Ure_h = new double[szU];    // Ure_h stores real U 
   double *Uim_h = new double[szU];    // Uim_h stores imaginary U 
   double *Ure_d;                      // Ure_d is Ure_h on the device
   double *Uim_d;                      // Uim_d is Uim_h on the device
   double *invUre_h = new double[szU]; // real parts of the inverse
   double *invUim_h = new double[szU]; // imaginary parts of the inverse
   double *invUre_d;                   // invUre_d is invUre_h on the device
   double *invUim_d;                   // invUim_d is invUim_h on the device

   int ix = 0;
   for(int j=0; j<dim; j++)
      for(int i=0; i<dim; i++)
      {
         Ure_h[ix]   = Ure[i][j];
         Uim_h[ix++] = Uim[i][j];
      }

   size_t szmat = szU*sizeof(double);
   cudaMalloc((void**)&Ure_d,szmat);
   cudaMalloc((void**)&Uim_d,szmat);
   cudaMalloc((void**)&invUre_d,szmat);
   cudaMalloc((void**)&invUim_d,szmat);
   cudaMemcpy(Ure_d,Ure_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(Uim_d,Uim_h,szmat,cudaMemcpyHostToDevice);

   if(dim <= 16)
      cmplx_small_invert_upper<<<1,dim>>>(dim,Ure_d,Uim_d,invUre_d,invUim_d);
   else
      cmplx_medium_invert_upper<<<1,dim>>>(dim,Ure_d,Uim_d,invUre_d,invUim_d);

   cudaMemcpy(invUre_h,invUre_d,szmat,cudaMemcpyDeviceToHost);
   cudaMemcpy(invUim_h,invUim_d,szmat,cudaMemcpyDeviceToHost);

   ix = 0;
   for(int i=0; i<dim; i++)
      for(int j=0; j<dim; j++)
      {
         invUre[i][j] = invUre_h[ix];
         invUim[i][j] = invUim_h[ix++];
      }

   free(Ure_h); free(invUre_h);
   free(Uim_h); free(invUim_h);
}

void GPU_dbl_upper_tiled_solver
 ( int dim, int szt, int nbt, double **U, double *b, double *x )
{
   const int nbr = nbt*szt*szt;   // number of doubles on diagonal tiles
   double *D_h = new double[nbr];    // the diagonal tiles on the host
   double *D_d;                      // diagonal tiles on the device
   double *invD_h = new double[nbr]; // inverse of diagonal tiles on host 
   double *invD_d;                   // invD_d is invD_h on device
   int offset;
   int ix = 0;

   for(int k=0; k<nbt; k++) // copy columns of the k-th tile
   {
      offset = k*szt;
      for(int j=0; j<szt; j++)
         for(int i=0; i<szt; i++) D_h[ix++] = U[offset+i][offset+j];
   }
   const size_t sznum = nbr*sizeof(double);
   cudaMalloc((void**)&D_d,sznum);
   cudaMalloc((void**)&invD_d,sznum);
   cudaMemcpy(D_d,D_h,sznum,cudaMemcpyHostToDevice);

   dbl_invert_tiles<<<nbt,szt>>>(szt,D_d,invD_d);

   double *rhs_d;                    // right hand side on device
   const size_t szrhs = dim*sizeof(double);
   cudaMalloc((void**)&rhs_d,szrhs);
   cudaMemcpy(rhs_d,b,szrhs,cudaMemcpyHostToDevice);

   dbl_multiply_inverse<<<1,szt>>>(szt,nbt-1,invD_d,rhs_d);

   int nbrUcol = (nbt-1)*szt*szt;         // #doubles in column of U
   double *Ucol_h = new double[nbrUcol];  // column of U on host
   double *Ucol_d;
   const size_t szUcol = nbrUcol*sizeof(double);
   cudaMalloc((void**)&Ucol_d,szUcol);

   int coloff,rowoff;

   for(int k=nbt-1; k>0; k--)      // update with solution tile k
   {
      coloff = k*szt;      // column offset to update with solution tile k
      ix = 0;
      for(int L=0; L<k; L++)       // copy k tiles of U
      {
         rowoff = L*szt;           // row offset for update data
         for(int i=0; i<szt; i++)
            for(int j=0; j<szt; j++) Ucol_h[ix++] = U[rowoff+i][coloff+j];
      }
      cudaMemcpy(Ucol_d,Ucol_h,nbrUcol*sizeof(double),
                 cudaMemcpyHostToDevice);

      dbl_back_substitute<<<k,szt>>>(szt,k,Ucol_d,rhs_d);

      // (k-1)-th solution tile is ready for inverse multiplication
      dbl_multiply_inverse<<<1,szt>>>(szt,k-1,invD_d,rhs_d);

      nbrUcol = nbrUcol - szt*szt; // one tile less used in update
   }
   cudaMemcpy(x,rhs_d,szrhs,cudaMemcpyDeviceToHost);

   // copy of invD_d is needed only for testing purposes
   cudaMemcpy(invD_h,invD_d,sznum,cudaMemcpyDeviceToHost);

   ix = 0;
   for(int k=0; k<nbt; k++) // copy rows of the inverse of the k-th tile
   {
      offset = k*szt;
      for(int i=0; i<szt; i++)
         for(int j=0; j<szt; j++) U[offset+i][offset+j] = invD_h[ix++];
   }
   free(D_h); free(invD_h); free(Ucol_h);
}
