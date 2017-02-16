#ifndef __PATH_GPU_MGS_LARGE_BACKSUB_CU__
#define __PATH_GPU_MGS_LARGE_BACKSUB_CU__

__global__ void mgs_large_backsubstitution_kernel1
 ( GT* R, GT* x, int dim, int pivot, int BS0, int BS )
{
   int tidx = threadIdx.x;
   __shared__ GT sol[BS_QR_Back];
   __shared__ GT Rcl[BS_QR_Back];
   __shared__ GT update[BS_QR_Back];
   // GT update;
   int offset = pivot*BS;

   update[tidx] = R[offset+tidx];
   // __syncthreads();
   for(int k=BS0-1; k>=0; k--) // compute k-th component of solution
   {
      if(tidx < k+1)
      {
         int ind = (dim - k - offset)*(dim + 3 + k + offset)/2 + tidx;
         Rcl[tidx] = R[ind+offset];
      }
      if(tidx == k)
         sol[tidx] = update[tidx]/Rcl[tidx]; // all other threads wait
      __syncthreads();
      if(tidx < k) update[tidx] = update[tidx] - sol[k]*Rcl[tidx]; // update
      // __syncthreads();
   }
   x[offset+tidx] = sol[tidx];
}

__global__ void mgs_large_backsubstitution_kernel2
 ( GT* R, GT* x, int dim, int pivot, int BS0, int BS )
{
   int tidx = threadIdx.x;
   int b = blockIdx.x+1;
   __shared__ GT sol[BS_QR_Back];
   __shared__ GT Rcl[BS_QR_Back];
   int ind;
   GT update;
   int offset = pivot*BS;//+BS0;

   if(tidx < BS0)
   {
      sol[tidx] = x[offset+tidx];
   }
   if(BS>32)
   {
      __syncthreads();
   }
   int block_offset = b*BS;
   update = R[offset-block_offset+tidx];

   for(int k=BS0-1; k>=0; k--)  // continue updates
   {
      ind = (dim - k - offset)*(dim + 3 + k + offset)/2 + tidx;
      Rcl[tidx] = R[ind+offset-block_offset];
      update = update - sol[k]*Rcl[tidx];  // update
   }
   R[offset-block_offset+tidx] = update;
}


void mgs_large_backsubstitution ( GT* R, GT* sol, int rows, int cols )
{
   int BS = BS_QR_Back;

   int rf = ceil(((double) (cols-1))/BS);
   // std::cout << "rf = " << rf << std::endl;
   int BS_col = cols-1 - BS*(rf-1);

   // std::cout << "BS_col = " << BS_col << std::endl;

   for(int piv=rf-1; piv>=0; piv--) 
   {
      // std::cout<< "piv = " << piv << std::endl;
      if(piv==rf-2) BS_col = BS;
      mgs_large_backsubstitution_kernel1<<<1,BS_col>>>
         (R,sol,cols-1,piv,BS_col,BS);
      if(piv==0) break;
      mgs_large_backsubstitution_kernel2<<<piv,BS>>>
        (R,sol,cols-1,piv,BS_col,BS);
   }
}

#endif /*__PATH_GPU_MGS_LARGE_BACKSUB_CU__*/
