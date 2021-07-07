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

__global__ void  dbl2_invert_tiles
 ( int dim, double *Uhi, double *Ulo, double *invUhi, double *invUlo )
{
   const int B = blockIdx.x;   // block index
   const int k = threadIdx.x;  // thread k computes k-th column of inverse
   const int offset = dim*dim*B; // offset in U and invU

   __shared__ double Ucolhi[dd_shmemsize];      // one column of U
   __shared__ double Ucollo[dd_shmemsize];
   __shared__ double invUrowhi[dd_shmemsize];   // one row of invU
   __shared__ double invUrowlo[dd_shmemsize]; 

   double rhshi,rhslo,xvalhi,xvallo,acchi,acclo;

   int colidx = offset + dim*(dim-1); // start with the last column

   Ucolhi[k] = Uhi[colidx+k];         // load the last column
   Ucollo[k] = Ulo[colidx+k];
   rhshi = ((double) int(k == dim-1));  // right hand side for each thread
   rhslo = 0.0;
   int rowidx = offset + (dim - 1)*dim + k; // row index in the inverse

   // invUrow[k] = rhs/Ucol[k];       // last row of the inverse
   ddg_div(rhshi,rhslo,Ucolhi[k],Ucollo[k],&invUrowhi[k],&invUrowlo[k]);
   invUhi[rowidx] = invUrowhi[k];     // store the last row into invU
   invUlo[rowidx] = invUrowlo[k];

   for(int i=dim-2; i>=0; i--)        // compute row with index i
   {
      rhshi = ((double) int(k == i));   // set rhs for i-th unit vector
      rhslo = 0.0;

      for(int j=i+1; j<dim; j++)
      {
         colidx = offset + dim*j;     // need column j of U
         Ucolhi[k] = Uhi[colidx+k];
         Ucollo[k] = Ulo[colidx+k];

         rowidx = offset + j*dim + k; // need solution value
         invUrowhi[k] = invUhi[rowidx]; // load invU row into invUrow
         invUrowlo[k] = invUlo[rowidx]; // load invU row into invUrow
         xvalhi = invUrowhi[k];
         xvallo = invUrowlo[k];

         __syncthreads();
         // rhs = rhs - Ucol[i]*xval;    // update right hand side
         ddg_mul(Ucolhi[i],Ucollo[i],xvalhi,xvallo,&acchi,&acclo);
         ddg_dec(&rhshi,&rhslo,acchi,acclo);
      }
      colidx = offset + dim*i;        // need column i of U
      Ucolhi[k] = Uhi[colidx+k];
      Ucollo[k] = Ulo[colidx+k];
      rowidx = offset + i*dim + k;    // save in i-th row of inverse

      __syncthreads();
      // invUrow[k] = rhs/Ucol[i];
      ddg_div(rhshi,rhslo,Ucolhi[i],Ucollo[i],&invUrowhi[k],&invUrowlo[k]);
      invUhi[rowidx] = invUrowhi[k];
      invUlo[rowidx] = invUrowlo[k];
   }
}

__global__ void dbl2_multiply_inverse
 ( int dim, int idx, double *invUhi, double *invUlo,
   double *whi, double *wlo )
{
   const int k = threadIdx.x;     // thread k computes k-th product
   const int rhsoff = dim*idx;    // offset for the right hand size
   const int offset = dim*rhsoff; // offset for diagonal tile

   __shared__ double workhi[dd_shmemsize];      // copy of w
   __shared__ double worklo[dd_shmemsize];      // copy of w

   workhi[k] = whi[rhsoff+k];
   worklo[k] = wlo[rhsoff+k];

   double resulthi = 0.0; // each thread stores its product in result
   double resultlo = 0.0;
   double coeffhi,coefflo,acchi,acclo;

   for(int j=0; j<dim; j++)  // column j of the inverse diagonal tile
   {
      coeffhi = invUhi[offset+k*dim+j]; // thread k does row k
      coefflo = invUlo[offset+k*dim+j];
      // result = result + coeff*work[j];
      ddg_mul(coeffhi,coefflo,workhi[j],worklo[j],&acchi,&acclo);
      ddg_inc(&resulthi,&resultlo,acchi,acclo);
   }
   whi[rhsoff+k] = resulthi;
   wlo[rhsoff+k] = resultlo;
}

__global__ void dbl2_back_substitute
 ( int dim, int idx, double *Uhi, double *Ulo, double *whi, double *wlo )
{
   const int B = blockIdx.x;     // block index
   const int k = threadIdx.x;    // thread k computes k-th product
   const int offset = B*dim*dim; // numbers to skip

   __shared__ double wrkhi[dd_shmemsize];   // copy of w
   __shared__ double wrklo[dd_shmemsize]; 
   __shared__ double solhi[dd_shmemsize];    // solution to update with
   __shared__ double sollo[dd_shmemsize];

   wrkhi[k] = whi[B*dim+k];    // block B updates B-th slice of w
   wrklo[k] = wlo[B*dim+k];
   solhi[k] = whi[idx*dim+k];  // solution that is back substituted
   sollo[k] = wlo[idx*dim+k];

   double resulthi = 0.0; // each thread stores its product in result
   double resultlo = 0.0;
   double coeffhi,coefflo,acchi,acclo;

   for(int j=0; j<dim; j++)  // column j of the inverse diagonal tile
   {
      coeffhi = Uhi[offset+k*dim+j];
      coefflo = Ulo[offset+k*dim+j];
      // result = result + coeff*sol[j];
      ddg_mul(coeffhi,coefflo,solhi[j],sollo[j],&acchi,&acclo);
      ddg_inc(&resulthi,&resultlo,acchi,acclo);
   }
   // wrk[k] = wrk[k] - result; // subtract product
   ddg_dec(&wrkhi[k],&wrklo[k],resulthi,resultlo);
   whi[B*dim+k] = wrkhi[k];
   wlo[B*dim+k] = wrklo[k];
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
         Uhi_h[ix]   = Uhi[i][j];
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

void GPU_dbl2_upper_tiled_solver
 ( int dim, int szt, int nbt, double **Uhi, double **Ulo,
   double *bhi, double *blo, double *xhi, double *xlo )
{
   const int nbr = nbt*szt*szt;   // number of doubles on diagonal tiles
   double *Dhi_h = new double[nbr];    // the diagonal tiles on the host
   double *Dlo_h = new double[nbr];    // low doubles of diagonal tiles
   double *Dhi_d;                      // diagonal tiles on the device
   double *Dlo_d;                      // low doubles of diagonal tiles
   double *invDhi_h = new double[nbr]; // inverse of diagonal tiles on host 
   double *invDlo_h = new double[nbr]; // low doubles of inverse tiles
   double *invDhi_d;                   // invDhi_d is invDhi_h on device
   double *invDlo_d;                   // invDlo_d is invDlo_h on device
   int offset;
   int ix = 0;

   for(int k=0; k<nbt; k++) // copy columns of the k-th tile
   {
      offset = k*szt;
      for(int j=0; j<szt; j++)
         for(int i=0; i<szt; i++)
         {
            Dhi_h[ix]   = Uhi[offset+i][offset+j];
            Dlo_h[ix++] = Ulo[offset+i][offset+j];
         }
   }
   const size_t sznum = nbr*sizeof(double);
   cudaMalloc((void**)&Dhi_d,sznum);
   cudaMalloc((void**)&Dlo_d,sznum);
   cudaMalloc((void**)&invDhi_d,sznum);
   cudaMalloc((void**)&invDlo_d,sznum);
   cudaMemcpy(Dhi_d,Dhi_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Dlo_d,Dlo_h,sznum,cudaMemcpyHostToDevice);

   dbl2_invert_tiles<<<nbt,szt>>>(szt,Dhi_d,Dlo_d,invDhi_d,invDlo_d);

   double *rhshi_d;                    // right hand side on device
   double *rhslo_d;
   const size_t szrhs = dim*sizeof(double);
   cudaMalloc((void**)&rhshi_d,szrhs);
   cudaMalloc((void**)&rhslo_d,szrhs);
   cudaMemcpy(rhshi_d,bhi,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(rhslo_d,blo,szrhs,cudaMemcpyHostToDevice);

   dbl2_multiply_inverse<<<1,szt>>>
      (szt,nbt-1,invDhi_d,invDlo_d,rhshi_d,rhslo_d);

   int nbrUcol = (nbt-1)*szt*szt;           // #doubles in column of U
   double *Ucolhi_h = new double[nbrUcol];  // column of U on host
   double *Ucollo_h = new double[nbrUcol];  // column of U on host
   double *Ucolhi_d;
   double *Ucollo_d;
   const size_t szUcol = nbrUcol*sizeof(double);
   cudaMalloc((void**)&Ucolhi_d,szUcol);
   cudaMalloc((void**)&Ucollo_d,szUcol);

   int coloff,rowoff;

   for(int k=nbt-1; k>0; k--)      // update with solution tile k
   {
      coloff = k*szt;      // column offset to update with solution tile k
      ix = 0;
      for(int L=0; L<k; L++)       // copy k tiles of U
      {
         rowoff = L*szt;           // row offset for update data
         for(int i=0; i<szt; i++)
            for(int j=0; j<szt; j++)
            {
               Ucolhi_h[ix]   = Uhi[rowoff+i][coloff+j];
               Ucollo_h[ix++] = Ulo[rowoff+i][coloff+j];
            }
      }
      cudaMemcpy(Ucolhi_d,Ucolhi_h,nbrUcol*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(Ucollo_d,Ucollo_h,nbrUcol*sizeof(double),
                 cudaMemcpyHostToDevice);

      dbl2_back_substitute<<<k,szt>>>
         (szt,k,Ucolhi_d,Ucollo_d,rhshi_d,rhslo_d);

      // (k-1)-th solution tile is ready for inverse multiplication
      dbl2_multiply_inverse<<<1,szt>>>
         (szt,k-1,invDhi_d,invDlo_d,rhshi_d,rhslo_d);

      nbrUcol = nbrUcol - szt*szt; // one tile less used in update
   }
   cudaMemcpy(xhi,rhshi_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(xlo,rhslo_d,szrhs,cudaMemcpyDeviceToHost);

   // copy of invD_d is needed only for testing purposes
   cudaMemcpy(invDhi_h,invDhi_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(invDlo_h,invDlo_d,sznum,cudaMemcpyDeviceToHost);

   ix = 0;
   for(int k=0; k<nbt; k++) // copy rows of the inverse of the k-th tile
   {
      offset = k*szt;
      for(int i=0; i<szt; i++)
         for(int j=0; j<szt; j++)
         {
            Uhi[offset+i][offset+j] = invDhi_h[ix];
            Ulo[offset+i][offset+j] = invDlo_h[ix++];
         }
   }
   free(Dhi_h); free(invDhi_h); free(Ucolhi_h);
   free(Dlo_h); free(invDlo_h); free(Ucollo_h);
}
