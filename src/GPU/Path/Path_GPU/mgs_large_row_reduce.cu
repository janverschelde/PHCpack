#ifndef __PATH_GPU_MGS_LARGE_ROW_REDUCE_CU__
#define __PATH_GPU_MGS_LARGE_ROW_REDUCE_CU__

__global__ void mgs_large_row_reduce_kernel
 ( GT* v, GT* R, int cols, int rows, int rowsLog2,
   int pivot, int rnd, int rndLog2, int BS, int BSLog2, 
   T *pivnorm, int lastBSLog2, int piv_end = 0 )
{
   int b = blockIdx.x + 1 + piv_end;
   int j = threadIdx.x;
   int block = b+pivot; // column for reduction w.r.t. pivot
   int vBSind = 0;
   int powTwo;

   __shared__ GT piv[BS_QR];      // contains pivot column
   __shared__ GT shv[BS_QR];      // for the reduction
   __shared__ GT sums[maxrounds]; // partial sums in rounds

   vBSind = 0;
   for(int i=0; i<rnd; i++) // normalize and partial sums for inner product
   {
      if(vBSind+j<rows)
      {
         piv[j] = v[pivot*rows+j+vBSind];
         shv[j] = v[block*rows+j+vBSind];
         shv[j] = piv[j].adj()*shv[j];
      }
      if(i==rnd-1)
      {
         BSLog2 = lastBSLog2;
      }
      powTwo = 1<<(BSLog2-1); // sum for the norm

      if(powTwo > 16)
      {
         __syncthreads();
      }
      if(j+powTwo<BS && vBSind+j+powTwo<rows)
      {
         shv[j] = shv[j] + shv[j+powTwo];
      }
      for(int k=0; k<BSLog2-1; k++)
      {
         powTwo = powTwo/2;
         if(powTwo > 16)
         {
            __syncthreads();
         }
         if(j < powTwo)
         {
            shv[j] = shv[j] + shv[j+powTwo];
         }
         // __syncthreads();
      }
      if(j == 0) sums[i] = shv[0];
      __syncthreads();  // avoid shv[0] is changed in next round
      vBSind = vBSind + BS;
   }
   if(rndLog2 > 0)
   {
      powTwo = 1<<(rndLog2-1); // sum for the norm

      if(j + powTwo < rnd)
      {
         sums[j] = sums[j] + sums[j+powTwo];
      }
      for(int k=0; k < rndLog2-1; k++)
      {
         powTwo = powTwo/2;
         // Maxround < 32, so it is not necessary to sync.
         if(powTwo > 16)
         {
            __syncthreads();
         }
         if(j < powTwo)
         {
            sums[j] = sums[j] + sums[j+powTwo];
         }
      }
   }
   __syncthreads();

   vBSind = 0;
   for(int i=0; i<rnd; i++)            // perform reduction
   {
      if(vBSind+j < rows)
      {
         piv[j] = v[pivot*rows+j+vBSind];
         shv[j] = v[block*rows+j+vBSind];
         shv[j] = shv[j] - sums[0]*piv[j];
         v[block*rows+j+vBSind] = shv[j];
      }
      vBSind = vBSind + BS;
   }
   if(j == 0)
   {
      int indR = r_pos(pivot, pivot+b, cols);
      // (dimR-1) - (pivot*(pivot+1))/2 - (b*(b+1))/2 - b*(pivot+1);
      R[indR] = sums[0];
   }
}

#endif /*__PATH_GPU_MGS_LARGE_ROW_REDUCE_CU__*/
