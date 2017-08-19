#include "large_norm.cu"
#include "large_row_reduce.cu"

__global__ void mgs_large_row_backsubstitution_kernel
 ( GT* R, GT* x, int dim, int rnd, int pivot, int BS )
{
   int tidx = threadIdx.x;
   int b = blockIdx.x;
   __shared__ GT sol[BS_QR_Back];
   __shared__ GT Rcl[BS_QR_Back];
   int ind;
   GT update;
   int offset = pivot*BS;

   if(pivot == rnd-1)
      update = R[offset+tidx];
   else
      update = x[offset+tidx];
   __syncthreads();
   for(int k=BS-1; k>=0; k--)  // compute k-th component of solution
   {
      if(tidx < k+1)
      {
         ind = (dim - k - offset)*(dim + 3 + k + offset)/2 + tidx;
         Rcl[tidx] = R[ind+offset];
      }
      if(tidx == k) sol[tidx] = update/Rcl[tidx]; // all other threads wait
      __syncthreads();
      if(tidx < k) update = update - sol[k]*Rcl[tidx];// update
      __syncthreads();
   }
   if(b == 0) x[offset+tidx] = sol[tidx];
   if(b != 0)
   {
      int block_offset = b*BS;
      if(pivot == rnd-1)
         update = R[offset-block_offset+tidx];
      else
         update = x[offset-block_offset+tidx];
      for(int k=BS-1; k>=0; k--)  // continue updates
      {
         ind = (dim - k - offset)*(dim + 3 + k + offset)/2 + tidx;
         __syncthreads();
         Rcl[tidx] = R[ind+offset-block_offset];
         __syncthreads();
         update = update - sol[k]*Rcl[tidx];  // update
      }
      __syncthreads();
      x[offset-block_offset+tidx] = update;
   }
}

void mgs_large_old ( GT* V, GT* R, GT* sol, int rows, int cols )
{
   int BS = min(BS_QR,rows);
   // int BS = 32;

   int rowsLog2 = log2ceil(rows); // ceil for sum reduction
   // int dimR = cols*(cols+1)/2;

   T* pivnrm; // norm of the pivot column
   cudaMalloc((void**)&pivnrm,sizeof(T));

   T* sums_global; // norm of the pivot column
   cudaMalloc((void**)&sums_global,maxrounds*sizeof(T));

   int rf = ceil(((double) rows)/BS);
   int rfLog2 = log2ceil(rf);
   int BSLog2 = log2ceil(BS);
   int lastBSLog2 = log2ceil(rows-BS*(rf-1));

   /*
     std::cout << "BS     = " << BS << std::endl;
     std::cout << "rf     = " << rf << std::endl;
     std::cout << "rfLog2 = " << rfLog2 << std::endl;
     std::cout << "BSLog2 = " << BSLog2 << std::endl;
     std::cout << "lastBSLog2 = " << lastBSLog2 << std::endl;
    */

   for(int piv=0; piv<cols-1; piv++)
   {
      mgs_large_normalize_kernel1<<<rf,BS>>>
         (V,R,rows,rowsLog2,cols,piv,rf,rfLog2,BS,BSLog2,pivnrm,
          lastBSLog2,sums_global);
      mgs_large_normalize_kernel2<<<rf,BS>>>
         (V,R,rows,rowsLog2,cols,piv,rf,rfLog2,BS,BSLog2,pivnrm,
          lastBSLog2,sums_global);
      // XXX BS should be greater than maxround
      mgs_large_row_reduce_kernel<<<cols-piv-1,BS>>>
         (V,R,cols,rows,rowsLog2, piv,rf,rfLog2,BS,BSLog2,pivnrm,lastBSLog2);
   }
   // BS = BS_QR_Back;

   rf = ceil(((double) (cols-1))/BS);
   for(int piv=rf-1; piv>=0; piv--)
   {
      mgs_large_row_backsubstitution_kernel<<<piv+1,BS>>>
         (R,sol,cols-1,rf,piv,BS);
   }
}
