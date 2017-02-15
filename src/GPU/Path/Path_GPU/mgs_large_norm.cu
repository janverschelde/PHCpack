// kernels to compute the norm of a large vector

#ifndef __PATH_GPU_MGS_LARGE_NORM_CU__
#define __PATH_GPU_MGS_LARGE_NORM_CU__

__global__ void mgs_large_normalize_kernel1
 ( GT* v, GT* R, int rows, int rowsLog2, int cols, int pivot,
   int rnd, int rndLog2, int BS, int BSLog2, T *pivnorm,
   int lastBSLog2, T* sums_global )
{
   int b = blockIdx.x;
   int j = threadIdx.x;
   int i = b*BS + j;

   __shared__ T shv[shmemsize/2]; // for the reduction

   if(b == rnd-1)
   {
      BSLog2 = lastBSLog2;
   }
   BSLog2 -= 1;

   if(i < rows)
   {
      GT tmp_piv = v[pivot*rows+i];
      shv[j] = tmp_piv.real*tmp_piv.real + tmp_piv.imag*tmp_piv.imag;
   }
   int half_size = 1 << (BSLog2);// sum for the norm

   if(half_size > 16) 
   {
      __syncthreads();
   }
   if(i + half_size < rows && j < half_size)
   {
      shv[j] = shv[j] + shv[j+half_size];
   }
   for(int k=0; k < BSLog2; k++)
   {
      half_size /= 2;
      if(half_size > 16)
      {
         __syncthreads();
      }
      if(j < half_size)
      {
         shv[j] = shv[j] + shv[j+half_size];
      }
   }
   if(j == 0) sums_global[b] = shv[0];
}

__global__ void mgs_large_normalize_kernel2
 ( GT* v, GT* R, int rows, int rowsLog2, int cols, int pivot,
   int rnd, int rndLog2, int BS, int BSLog2, T *pivnorm,
   int lastBSLog2, T* sums_global )
{
   int b = blockIdx.x;
   int j = threadIdx.x;
   int i = b*BS + j;

   __shared__ T sums[maxrounds];// partial sums in rounds
   // maxrounds 32 reduce sync

   T newpivnorm;// norm of the pivot

   if(j < rnd)
   {
      sums[j] = sums_global[j];
   }
   if(rndLog2 > 0)
   {
      int powTwo = 1<<(rndLog2-1); // sum for the norm

      if(powTwo>16)
      {
         __syncthreads();
      }
      if(j + powTwo < rnd)
      {
         sums[j] = sums[j] + sums[j+powTwo];
      }
      for(int k=0; k < rndLog2-1; k++)
      {
         if(powTwo>16)
         {
            __syncthreads();
         }
         powTwo = powTwo/2;
         if(j < powTwo)
         {
            sums[j] = sums[j] + sums[j+powTwo];
         }
      }
   }
   __syncthreads();

   newpivnorm = sqrt(sums[0]);

   if(i<rows)          // exclude extra threads in last round
   {
      v[pivot*rows+i] /= newpivnorm;
   }
   if(i == 0)
   {
      int indR = r_pos(pivot, pivot+b, cols);
      // (dimR-1) - (pivot*(pivot+1))/2 - (b*(b+1))/2 - b*(pivot+1);
      R[indR].init_imag();
      R[indR].real = newpivnorm;
   }
}

#endif /*__PATH_GPU_MGS_LARGE_NORM_CU__*/
