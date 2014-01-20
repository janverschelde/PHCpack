// Kernels to compute the Gram matrix of a sequence of complex vectors
// are defined below.  Because this code prepares for the revision of the
// modified Gram-Schmidt method, a double indexed block structure is not
// applied, instead, kernels are launced with a double loop.

#include <gqd.cu>
#include <gqd_type.h>
#include "complex.h"

#define d  0 
#define dd 1
#define qd 2

#ifdef precision
#define p precision
#else
#define p 0
#endif

#if(p == 0)
typedef double T;
#define shmemsize 256
#elif(p == 1)
typedef gdd_real T;
#define shmemsize 128
#else
typedef gqd_real T;
#define shmemsize 64
#endif

#define maxrounds 32

using namespace std;

__global__ void small_normalize_and_reduce
 ( complex<T>* v, int rows, int rowsLog2, int cols, int pivot )
{
   int block = blockIdx.x+pivot; // column for reduction w.r.t. pivot
   int j = threadIdx.x;
   __shared__ complex<T> piv[shmemsize];    // contains pivot column
   __shared__ complex<T> shv[2][shmemsize]; // for the reduction
   __shared__ T prd[shmemsize];             // for norm of the pivot
   int powTwo;

   piv[j] = v[pivot*rows+j];
   prd[j] = piv[j].real*piv[j].real + piv[j].imag*piv[j].imag;
   __syncthreads();

   powTwo = 1;                            // sum for the norm
   for(int k=0; k < rowsLog2; k++)
   {
      if((j%(powTwo*2)) == 0)
         if(j+powTwo < rows) prd[j] = prd[j] + prd[j+powTwo];
      powTwo = powTwo*2;
      __syncthreads();
   }
   if(j == 0) prd[0] = sqrt(prd[0]);

   __syncthreads();

   piv[j] = piv[j]/prd[0];
   if(block == pivot) v[pivot*rows+j] = piv[j];
   __syncthreads();
   if(block != pivot)
   {
      shv[0][j] = v[block*rows+j];
      shv[1][j] = piv[j].adj()*shv[0][j];
      __syncthreads();
      powTwo = 1;                          // sum for inner product
      for(int k=0; k < rowsLog2; k++)
      {
         if((j%(powTwo*2)) == 0)
            if(j+powTwo < rows) shv[1][j] = shv[1][j] + shv[1][j+powTwo];
         powTwo = powTwo*2;
         __syncthreads();
      }
      shv[0][j] = shv[0][j] - shv[1][0]*piv[j];
      __syncthreads();
      v[block*rows+j] = shv[0][j];
   }
}

__global__ void small_QR_normalize_and_reduce
 ( complex<T>* v, complex<T>* R,
   int dimR, int rows, int rowsLog2, int cols, int pivot )
{
   int b = blockIdx.x;
   int j = threadIdx.x;
   int block = b+pivot;    // column for reduction w.r.t. pivot
   int i = block*rows + j;
   int L = pivot*rows + j;
   __shared__ complex<T> piv[shmemsize];    // contains pivot column
   __shared__ complex<T> shv[2][shmemsize]; // for the reduction
   __shared__ T prd[shmemsize];             // for norm of the pivot
   int powTwo;
   int indR = (dimR-1) - (pivot*(pivot+1))/2 - (b*(b+1))/2 - b*(pivot+1);

   piv[j] = v[L];
   prd[j] = piv[j].real*piv[j].real + piv[j].imag*piv[j].imag;
   __syncthreads();

   powTwo = 1;                            // sum for the norm
   for(int k=0; k < rowsLog2; k++)
   {
      if((j%(powTwo*2)) == 0)
         if(j+powTwo < rows) prd[j] = prd[j] + prd[j+powTwo];
      powTwo = powTwo*2;
      __syncthreads();
   }
   if(j == 0) prd[0] = sqrt(prd[0]);
   __syncthreads();
   piv[j] = piv[j]/prd[0];
   if(block == pivot) 
   {
       v[L] = piv[j];
       if(j == 0)
       {
          R[indR].init(0.0,0.0);
          R[indR].real = prd[0];
       }
   }
   __syncthreads();
   if(block != pivot)
   {
      shv[0][j] = v[i];
      shv[1][j] = piv[j].adj()*shv[0][j];
      __syncthreads();
      powTwo = 1;                          // sum for inner product
      for(int k=0; k < rowsLog2; k++)
      {
         if((j%(powTwo*2)) == 0)
            if(j+powTwo < rows) shv[1][j] = shv[1][j] + shv[1][j+powTwo];
         powTwo = powTwo*2;
         __syncthreads();
      }
      shv[0][j] = shv[0][j] - shv[1][0]*piv[j];
      __syncthreads();
      v[i] = shv[0][j];
      if(j == 0) R[indR] = shv[1][0];
   }
}

__global__ void large_normalize
 ( complex<T>* v, int rows, int pivot, int BS, T *pivnorm )
{
   int block = blockIdx.x;
   int j = threadIdx.x;
   int ind = block*BS + j;
   __shared__ complex<T> piv[shmemsize];    // contains pivot column

   if(ind < rows)
   {
      piv[j] = v[pivot*rows + ind];
      piv[j] = piv[j]/(*pivnorm);
      v[pivot*rows + ind] = piv[j];
   }
}
__global__ void update_solution
 ( complex<T>* x_old, complex<T>* x_d, int rows, int BS )
{
   int block = blockIdx.x;
   int j = threadIdx.x;
   int ind = block*BS + j;
   __shared__ complex<T> piv[shmemsize];    // contains pivot column

   if(ind < rows)
   {
      x_old[ind] = x_old[ind] + x_d[ind];
   }
}

__global__ void large_normalize_and_reduce
 ( complex<T>* v, int rows, int rowsLog2, int cols, int pivot,
   int rnd, int rndLog2, int BS, int BSLog2, T* pivnorm )
{
   int block = blockIdx.x+pivot; // column for reduction w.r.t. pivot
   int j = threadIdx.x;
   int vBSind = 0;
   int powTwo;

   __shared__ complex<T> piv[shmemsize];    // contains pivot column
   __shared__ complex<T> shv[shmemsize];    // for the reduction
   __shared__ complex<T> sums[maxrounds];   // partial sums in rounds
   __shared__ T newpivnorm;                 // norm of the pivot

   for(int i=0; i<rnd; i++)
   {
      if(vBSind+j>=rows)          // exclude extra threads in last round
         shv[j].init(0.0,0.0);
      else
      {
         piv[j] = v[pivot*rows+j+vBSind];
         shv[j].real = piv[j].real*piv[j].real + piv[j].imag*piv[j].imag;
      }
      __syncthreads();
      powTwo = 1;                       // partial sums for the norm
      for(int k=0; k < BSLog2; k++)
      {
         if((j%(powTwo*2)) == 0)
            if(j+powTwo < BS)
               shv[j].real = shv[j].real + shv[j+powTwo].real;
         powTwo = powTwo*2;
         __syncthreads();
      }
      if(j == 0) sums[i].real = shv[0].real;
      __syncthreads();
      vBSind = vBSind + BS;
   }
   __syncthreads();
   powTwo = 1;                          // sum reduction for the norm
   for(int k=0; k < rndLog2; k++)
   {
      if((j%(powTwo*2)) == 0)
         if(j+powTwo < rnd)
            sums[j].real = sums[j].real + sums[j+powTwo].real;
      powTwo = powTwo*2;
      __syncthreads();
   }
   if(j == 0) newpivnorm = sqrt(sums[0].real);
   __syncthreads();
   vBSind = 0;
   for(int i=0; i<rnd; i++)  // normalize and partial sums for inner product
   {
      if(vBSind+j<rows)          // exclude extra threads in last round
      {
         if((block == pivot) && (pivot > 0)) // delayed normalization
         {
            piv[j] = v[(pivot-1)*rows+j+vBSind];
            piv[j] = piv[j]/(*pivnorm);
            v[(pivot-1)*rows+j+vBSind] = piv[j];
         }
         else  // every other block applies normalization to pivot column
         {
            piv[j] = v[pivot*rows+j+vBSind];
            piv[j] = piv[j]/newpivnorm;
         }
      }
      if(block != pivot)         // nonpivot blocks make inner product
      {
         if(vBSind+j>=rows)
            shv[j].init(0.0,0.0);
         else
         {
            shv[j] = v[block*rows+j+vBSind];
            shv[j] = piv[j].adj()*shv[j];
         }
         __syncthreads();
         powTwo = 1;             // partial sums for inner product
         for(int k=0; k < BSLog2; k++)
         {
            if((j%(powTwo*2)) == 0)
               if(j+powTwo < BS) shv[j] = shv[j] + shv[j+powTwo];
            powTwo = powTwo*2;
            __syncthreads();
         }
         if(j == 0) sums[i] = shv[0];
         __syncthreads();
      }
      __syncthreads();
      vBSind = vBSind + BS;
   }
   if(block == pivot) *pivnorm = newpivnorm;
   if(block != pivot)
   {
      powTwo = 1;                          // sum reduction for inner product
      for(int k=0; k < rndLog2; k++)
      {
         if((j%(powTwo*2)) == 0)
            if(j+powTwo < rnd)
               sums[j] = sums[j] + sums[j+powTwo];
         powTwo = powTwo*2;
         __syncthreads();
      }
      vBSind = 0;
      for(int i=0; i<rnd; i++)            // perform reduction
      {
         if(vBSind+j < rows)
         {
            piv[j] = v[pivot*rows+j+vBSind];
            shv[j] = v[block*rows+j+vBSind];
            shv[j] = shv[j] - sums[0]*piv[j]/newpivnorm;
            v[block*rows+j+vBSind] = shv[j];
         }
         __syncthreads();
         vBSind = vBSind + BS;
      }
   }
}

__global__ void chandra_evaluate_and_differentiate
 ( complex<T>* v, complex<T>* R, complex<T>* x, int dimR, int rows, int rowsLog2, int cols,
   int rnd, int rndLog2, int BS, int BSLog2, T *pivnorm )
{
   int b = blockIdx.x;
   int j = threadIdx.x;
   int block = b; // column for reduction w.r.t. pivot
   int vBSind = 0;
   int powTwo;
   complex<T> par;
   par.init(33.0,0.0);
   complex<T> local_tmp;
   local_tmp.init(64.0,0.0);
   par = par/local_tmp;

   //__shared__ complex<T> piv[shmemsize];    // contains pivot column
   __shared__ complex<T> shv[shmemsize];    // for the reduction
   __shared__ complex<T> sums[maxrounds];   // partial sums in rounds
   __shared__ complex<T> tmp;   // partial sums in rounds
   __shared__ T newpivnorm;                 // norm of the pivot

   for(int i=0; i<rnd; i++)
   {
      if(vBSind+j>=rows-1)          // exclude extra threads in last round
         shv[j].init(0.0,0.0);
      else
      {
         int x_ind = j+vBSind;
         local_tmp.init(x_ind+b+2,0);
         shv[j] = R[x_ind]/local_tmp;
      }
      __syncthreads();
      powTwo = 1;                       // partial sums for the norm
      for(int k=0; k < BSLog2; k++)
      {
         if((j%(powTwo*2)) == 0)
            if(j+powTwo < BS)
               //shv[j].real = shv[j].real + shv[j+powTwo].real;
               shv[j] = shv[j] + shv[j+powTwo];
         powTwo = powTwo*2;
         __syncthreads();
      }
      if(j == 0){ 
         //sums[i].real = shv[0].real;
         sums[i] = shv[0];
      }
      __syncthreads();
      vBSind = vBSind + BS;
   }
   __syncthreads();
   powTwo = 1;                          // sum reduction for the norm
   for(int k=0; k < rndLog2; k++)
   {
      if((j%(powTwo*2)) == 0)
         if(j+powTwo < rnd)
            sums[j] = sums[j] + sums[j+powTwo];
      powTwo = powTwo*2;
      __syncthreads();
   }
   //x[j] = shv[j];

   if(j==0){
      sums[0] = sums[0]*(b+1);
      sums[0].real = sums[0].real + 1;
      local_tmp = sums[0];
      local_tmp = local_tmp * par;
      int coef = 2*rows;
      local_tmp.real = local_tmp.real - coef;
      local_tmp = R[b] * local_tmp;
      local_tmp.real = local_tmp.real + coef;
      v[rows*rows+b] = local_tmp;


      local_tmp = sums[0];
      if(b!=rows-1){
         complex<T> local_tmp2;
         local_tmp2.init(2.0,0.0);
         local_tmp = local_tmp + R[b]/local_tmp2;
      }

      local_tmp = local_tmp * par;
      local_tmp.real = coef - local_tmp.real;
      local_tmp.imag = 0 - local_tmp.imag;

      v[rows*b+b] = local_tmp;

      tmp = R[b]*par*(-b-1);
   }

   __syncthreads();

   local_tmp = tmp;

   vBSind = 0;
   
   for(int i=0; i<rnd; i++)
   {
      int x_ind = j+vBSind;
      if(x_ind<(rows-1) && x_ind != b){
         complex<T> local_tmp2;
         local_tmp2.init(x_ind+b+2, 0);
         v[rows*x_ind+b] = local_tmp/local_tmp2;
         //v[rows*x_ind+b].init(b, x_ind+b+2);
      }
      vBSind = vBSind + BS;
   }
   if(j==0 and b != rows-1){
      v[rows*(rows-1)+b].init(0.0,0.0);
   }
}

__global__ void large_QR_normalize_and_reduce
 ( complex<T>* v, complex<T>* R, int dimR, int rows, int rowsLog2, int cols,
   int pivot, int rnd, int rndLog2, int BS, int BSLog2, T *pivnorm )
{
   int b = blockIdx.x;
   int j = threadIdx.x;
   int block = b+pivot; // column for reduction w.r.t. pivot
   int vBSind = 0;
   int powTwo;
   int indR = (dimR-1) - (pivot*(pivot+1))/2 - (b*(b+1))/2 - b*(pivot+1);

   __shared__ complex<T> piv[shmemsize];    // contains pivot column
   __shared__ complex<T> shv[shmemsize];    // for the reduction
   __shared__ complex<T> sums[maxrounds];   // partial sums in rounds
   __shared__ T newpivnorm;                 // norm of the pivot

   for(int i=0; i<rnd; i++)
   {
      if(vBSind+j>=rows)          // exclude extra threads in last round
         shv[j].init(0.0,0.0);
      else
      {
         piv[j] = v[pivot*rows+j+vBSind];
         shv[j].real = piv[j].real*piv[j].real + piv[j].imag*piv[j].imag;
      }
      __syncthreads();
      powTwo = 1;                       // partial sums for the norm
      for(int k=0; k < BSLog2; k++)
      {
         if((j%(powTwo*2)) == 0)
            if(j+powTwo < BS)
               shv[j].real = shv[j].real + shv[j+powTwo].real;
         powTwo = powTwo*2;
         __syncthreads();
      }
      if(j == 0) sums[i].real = shv[0].real;
      __syncthreads();
      vBSind = vBSind + BS;
   }
   __syncthreads();
   powTwo = 1;                          // sum reduction for the norm
   for(int k=0; k < rndLog2; k++)
   {
      if((j%(powTwo*2)) == 0)
         if(j+powTwo < rnd)
            sums[j].real = sums[j].real + sums[j+powTwo].real;
      powTwo = powTwo*2;
      __syncthreads();
   }
   if(j == 0)
   {
      newpivnorm = sqrt(sums[0].real);
      R[indR].init(0.0,0.0);
      R[indR].real = newpivnorm;
   }
   __syncthreads();
   vBSind = 0;
   for(int i=0; i<rnd; i++)  // normalize and partial sums for inner product
   {
      if(vBSind+j<rows)          // exclude extra threads in last round
      {
         if((block == pivot) && (pivot > 0)) // delayed normalization
         {
            piv[j] = v[(pivot-1)*rows+j+vBSind];
            piv[j] = piv[j]/(*pivnorm);
            v[(pivot-1)*rows+j+vBSind] = piv[j];
         }
         else  // every other block applies normalization to pivot column
         {
            piv[j] = v[pivot*rows+j+vBSind];
            piv[j] = piv[j]/newpivnorm;
         }
      }
      if(block != pivot)         // nonpivot blocks make inner product
      {
         if(vBSind+j>=rows)
            shv[j].init(0.0,0.0);
         else
         {
            shv[j] = v[block*rows+j+vBSind];
            shv[j] = piv[j].adj()*shv[j];
         }
         __syncthreads();
         powTwo = 1;             // partial sums for inner product
         for(int k=0; k < BSLog2; k++)
         {
            if((j%(powTwo*2)) == 0)
               if(j+powTwo < BS) shv[j] = shv[j] + shv[j+powTwo];
            powTwo = powTwo*2;
            __syncthreads();
         }
         if(j == 0) sums[i] = shv[0];
         __syncthreads();
      }
      __syncthreads();
      vBSind = vBSind + BS;
   }
   if(block == pivot) *pivnorm = newpivnorm;
   if(block != pivot)
   {
      powTwo = 1;                          // sum reduction for inner product
      for(int k=0; k < rndLog2; k++)
      {
         if((j%(powTwo*2)) == 0)
            if(j+powTwo < rnd)
               sums[j] = sums[j] + sums[j+powTwo];
         powTwo = powTwo*2;
         __syncthreads();
      }
      vBSind = 0;
      for(int i=0; i<rnd; i++)            // perform reduction
      {
         if(vBSind+j < rows)
         {
            piv[j] = v[pivot*rows+j+vBSind];
            shv[j] = v[block*rows+j+vBSind];
            shv[j] = shv[j] - sums[0]*piv[j]/newpivnorm;
            v[block*rows+j+vBSind] = shv[j];
         }
         __syncthreads();
         vBSind = vBSind + BS;
      }
      if(j == 0) R[indR] = sums[0];
   }
}

__global__ void small_backsubstitution
 ( complex<T>* R, complex<T>* x, int dim )
{
   int j = threadIdx.x;
   __shared__ complex<T> sol[shmemsize];
   __shared__ complex<T> Rcl[shmemsize];
   int ind;
   complex<T> update;

   update = R[j];
   __syncthreads();
   for(int k=dim-1; k>=0; k--)  // compute k-th component of solution
   {
      if(j < k+1)
      {
         ind = (dim - k)*(dim + 3 + k)/2 + j; 
         __syncthreads();
         Rcl[j] = R[ind];
         __syncthreads();
         if(j == k) sol[j] = update/Rcl[j]; // all other threads wait
         __syncthreads();
         if(j < k) update = update - sol[k]*Rcl[j];  // update
      }
      __syncthreads();
   }
   x[j] = sol[j];
}

__global__ void large_backsubstitution
 ( complex<T>* R, complex<T>* x, int dim, int rnd, int pivot, int BS )
{
   int j = threadIdx.x;
   int b = blockIdx.x;
   __shared__ complex<T> sol[shmemsize];
   __shared__ complex<T> Rcl[shmemsize];
   int ind;
   complex<T> update;
   int offset = pivot*BS;

   if(pivot == rnd-1)
      update = R[offset+j];
   else
      update = x[offset+j];
   __syncthreads();
   for(int k=BS-1; k>=0; k--)  // compute k-th component of solution
   {
      if(j < k+1)
      {
         ind = (dim - k - offset)*(dim + 3 + k + offset)/2 + j; 
         __syncthreads();
         Rcl[j] = R[ind+offset];
         __syncthreads();
         if(j == k) sol[j] = update/Rcl[j]; // all other threads wait
         __syncthreads();
         if(j < k) update = update - sol[k]*Rcl[j];  // update
      }
      __syncthreads();
   }
   if(b == 0) x[offset+j] = sol[j];
   if(b != 0)
   {
      int block_offset = b*BS;
      if(pivot == rnd-1)
         update = R[offset-block_offset+j];
      else
         update = x[offset-block_offset+j];
      for(int k=BS-1; k>=0; k--)  // continue updates
      {
         ind = (dim - k - offset)*(dim + 3 + k + offset)/2 + j; 
         __syncthreads();
         Rcl[j] = R[ind+offset-block_offset];
         __syncthreads();
         update = update - sol[k]*Rcl[j];  // update
      }
      __syncthreads();
      x[offset-block_offset+j] = update;
   }
}

void GPU_mgs2 ( complex<T>* v_h, int rows, int cols, int freq, int BS )
{
   int rowsLog2 = ceil(log2((double) rows)); // ceil for sum reduction
   complex<T>* v_d; // memory allocation on device and data transfer
   size_t size = rows*cols*sizeof(complex<T>);
   cudaMalloc((void**)&v_d,size);
   cudaMemcpy(v_d,v_h,size,cudaMemcpyHostToDevice);
   T* pivnrm; // norm of the pivot column
   cudaMalloc((void**)&pivnrm,sizeof(T));

   if(rows == BS)
      for(int piv=0; piv<cols; piv++)
         small_normalize_and_reduce<<<BS-piv,BS>>>(v_d,rows,rowsLog2,cols,piv);
   else
   {
      int rf = ceil(((double) rows)/BS);
      int rfLog2 = ceil(log2((double) rf));
      int BSLog2 = ceil(log2((double) BS));
      for(int piv=0; piv<cols; piv++)
         large_normalize_and_reduce<<<cols-piv,BS>>>
            (v_d,rows,rowsLog2,cols,piv,rf,rfLog2,BS,BSLog2,pivnrm);
      large_normalize<<<rf,BS>>>(v_d,rows,cols-1,BS,pivnrm);
   }
   cudaMemcpy(v_h,v_d,size,cudaMemcpyDeviceToHost);
}

/*void GPU_mgs2qr
 ( complex<T>* v_h, complex<T>* R_h,
   int dimR, int rows, int cols, int freq, int BS )
{
   int rowsLog2 = ceil(log2((double) rows)); // ceil for sum reduction
   complex<T>* v_d; // memory allocation on device and data transfer
   size_t v_size = rows*cols*sizeof(complex<T>);
   cudaMalloc((void**)&v_d,v_size);
   cudaMemcpy(v_d,v_h,v_size,cudaMemcpyHostToDevice);
   complex<T> *R_d;
   size_t R_size = dimR*sizeof(complex<T>);
   cudaMalloc((void**)&R_d,R_size);
   T* pivnrm; // norm of the pivot column
   cudaMalloc((void**)&pivnrm,sizeof(T));

   if(rows == BS)
      for(int piv=0; piv<cols; piv++)
         small_QR_normalize_and_reduce<<<BS-piv,BS>>>
            (v_d,R_d,dimR,rows,rowsLog2,cols,piv);
   else
   {
      int rf = ceil(((double) rows)/BS);
      int rfLog2 = ceil(log2((double) rf));
      int BSLog2 = ceil(log2((double) BS));
      for(int piv=0; piv<cols; piv++)
         large_QR_normalize_and_reduce<<<cols-piv,BS>>>
            (v_d,R_d,dimR,rows,rowsLog2,cols,piv,rf,rfLog2,BS,BSLog2,pivnrm);
      large_normalize<<<rf,BS>>>(v_d,rows,cols-1,BS,pivnrm);
   }

   cudaMemcpy(v_h,v_d,v_size,cudaMemcpyDeviceToHost);
   cudaMemcpy(R_h,R_d,R_size,cudaMemcpyDeviceToHost);
}*/

void GPU_mgs2qrls
 ( complex<T>* v_h, complex<T>* R_h, complex<T> *x_h,
   int dimR, int rows, int cols, int freq, int BS )
{
   int rowsLog2 = ceil(log2((double) rows)); // ceil for sum reduction
   complex<T>* v_d; // memory allocation on device and data transfer
   size_t v_size = rows*cols*sizeof(complex<T>);
   cudaMalloc((void**)&v_d,v_size);
   //cudaMemcpy(v_d,v_h,v_size,cudaMemcpyHostToDevice);
   complex<T> *R_d;
   size_t R_size = dimR*sizeof(complex<T>);
   cudaMalloc((void**)&R_d,R_size);
   complex<T> *x_d;
   size_t x_size = rows*sizeof(complex<T>);
   cudaMalloc((void**)&x_d,x_size);

   complex<T> *x_old;
   cudaMalloc((void**)&x_old,x_size);
   cudaMemcpy(x_old,x_h,x_size,cudaMemcpyHostToDevice);

   T* pivnrm; // norm of the pivot column
   cudaMalloc((void**)&pivnrm,sizeof(T));

   for(int i=0; i<freq; i++){
      int rf = ceil(((double) rows)/BS);
      int rfLog2 = ceil(log2((double) rf));
      int BSLog2 = ceil(log2((double) BS));
      chandra_evaluate_and_differentiate<<<rows,BS>>>
         (v_d,x_old,x_d,dimR,rows,rowsLog2,cols,rf,rfLog2,BS,BSLog2,pivnrm);
      for(int piv=0; piv<cols-1; piv++)
         large_QR_normalize_and_reduce<<<cols-piv,BS>>>
            (v_d,R_d,dimR,rows,rowsLog2,cols,piv,rf,rfLog2,BS,BSLog2,pivnrm);
      large_normalize<<<rf,BS>>>(v_d,rows,cols-2,BS,pivnrm);
      for(int piv=rf-1; piv>=0; piv--)
         large_backsubstitution<<<piv+1,BS>>>(R_d,x_d,rows,rf,piv,BS);
      update_solution<<<rf,BS>>>(x_old, x_d, rows, BS);
   }
   cudaMemcpy(x_h,x_old,x_size,cudaMemcpyDeviceToHost);
}
