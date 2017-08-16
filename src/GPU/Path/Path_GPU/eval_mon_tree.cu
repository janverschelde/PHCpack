// Sum block, mulitithread gpu sum unroll, for test
template <unsigned int n_th>
__global__ void eval_mon_tree_kernel
 ( GT* workspace_d, GT* x_d, GT* workspace_coef, int* mon_pos_start,
   unsigned short* pos_d, int n_mon )
{
   // GT* x_d, unsigned short* pos_d, GT* workspace_d, int n_mon,
   //                         int dim, int workspace_size_int){
   __shared__ GT x_sh[shmemsize];
   int BS = blockDim.x;
   int bidx = (gridDim.x*blockIdx.y+blockIdx.x)*BS;
   // int idx = bidx + threadIdx.x;
   int tidx = threadIdx.x;
   // int tidx2 = tidx + BS;

   int midx = tidx/n_th; // monomial index
   int midx_global =  midx + bidx/n_th;

   if(midx_global < n_mon)
   {
      // int sys_idx = blockIdx.z;
      // GT* x_d_tmp = x_d + sys_idx*dim;
      // GT* workspace_d_tmp = workspace_d + sys_idx*workspace_size_int;

      int pidx = tidx - midx*n_th; // thread index in monomial
      int pidx2 = pidx + n_th;
      // int xidx0 = midx_global*(n_th*2+1);
      int xidx0 = mon_pos_start[midx_global];
      int xidx1 = xidx0+ pidx+1; // pos index 1
      int xidx2 = xidx1 + n_th;  // pos index 2

      // int* pos = pos_d;// + BS*2*blockIdx.x;
      int n = pos_d[xidx0];

      GT* x_level = x_sh + midx*n_th*2;

      // Load to Shared Memory
      GT x0, x1;
      x0 = x_d[pos_d[xidx1]];
      if(pidx2 < n)
      {
         x1 = x_d[pos_d[xidx2]];
         x_level[pidx] = x0*x1;
      }
      else
      {
         x_level[pidx] = x0;
      }
      if(n_th > 32)
      {
         __syncthreads();
      }
      // Up
      GT* x_last_level;

      if(n_th > 256)
      {
         x_last_level = x_level; x_level += 512;
         if(pidx < 256)
         {
             x_level[pidx] = x_last_level[pidx] * x_last_level[pidx+256];
         }
         __syncthreads();
      }
      if(n_th > 128)
      {
         x_last_level = x_level; x_level += 256;
         if(pidx < 128)
         {
            x_level[pidx] = x_last_level[pidx] * x_last_level[pidx+128];
         }
         __syncthreads();
      }
      if(n_th > 64)
      {
         x_last_level = x_level; x_level += 128;
         if(pidx < 64)
         { 
            x_level[pidx] = x_last_level[pidx] * x_last_level[pidx+64];
         }
         __syncthreads();
      }
      if(n_th > 32)
      {
         x_last_level = x_level; x_level += 64;
         if(pidx < 32)
         {
            x_level[pidx] = x_last_level[pidx] * x_last_level[pidx+32];
         }
      }
      if(n_th > 16)
      {
         x_last_level = x_level; x_level += 32;
         if(pidx < 16)
         {
            x_level[pidx] = x_last_level[pidx] * x_last_level[pidx+16];
         }
      }
      if(n_th > 8)
      {
         x_last_level = x_level; x_level += 16;
         if(pidx <  8)
         {
            x_level[pidx] = x_last_level[pidx] * x_last_level[pidx+ 8];
         }
      }
      if(n_th > 4)
      {
         x_last_level = x_level; x_level +=  8;
         if(pidx <  4)
         {
            x_level[pidx] = x_last_level[pidx] * x_last_level[pidx+ 4];
         }
      }
      if(n_th > 2)
      {
         x_last_level = x_level; x_level +=  4;
         if(pidx <  2)
         {
            x_level[pidx] = x_last_level[pidx] * x_last_level[pidx+ 2];
         }
      }
      if(n_th > 1)
      {
         x_last_level = x_level; x_level +=  2;
         if(pidx <  1)
         {
            x_level[pidx] = x_last_level[pidx] * x_last_level[pidx+ 1];
         }
      }
      if(pidx < 3)
      {
         x_last_level[pidx] *= workspace_coef[midx_global];
      }
      // Common Factor to be multiplied
      if(pidx == 0)
      {
         workspace_d[xidx0] = x_level[0];
      }
      // Down
      if(n_th > 2)
      {
         x_level = x_last_level;
         x_last_level -= 4;
         if(pidx <  4)
         {
            int new_pidx = (pidx&1) ^ 1; // XOR to swith the term in last level
            x_last_level[pidx]*=x_level[new_pidx];
         }
      }
      if(n_th > 4)
      {
         x_level = x_last_level;
         x_last_level -= 8;
         if(pidx <  8)
         {
            int new_pidx = (pidx&3) ^ 2; // XOR to swith the term in last level
            x_last_level[pidx]*=x_level[new_pidx];
         }
      }
      if(n_th > 8)
      {
         x_level = x_last_level;
         x_last_level -= 16;
         if(pidx <  16)
         {
            int new_pidx = (pidx & 7) ^ 4;
            // XOR to swith the term in last level
            x_last_level[pidx]*=x_level[new_pidx];
         }
      }
      if(n_th > 16)
      {
         x_level = x_last_level;
         x_last_level -= 32;
         if(pidx < 32)
         {
            int new_pidx = (pidx & 15) ^ 8;
            // XOR to swith the term in last level
            x_last_level[pidx]*=x_level[new_pidx];
         }
         __syncthreads();
      }
      if(n_th > 32)
      {
         x_level = x_last_level;
         x_last_level -= 64;
         if(pidx < 64)
         {
            int new_pidx = (pidx & 31) ^ 16;
            // XOR to swith the term in last level
            x_last_level[pidx]*=x_level[new_pidx];
         }
         __syncthreads();
      }
      if(n_th > 64)
      {
         x_level = x_last_level;
         x_last_level -= 128;
         if(pidx < 128)
         {
            int new_pidx = (pidx & 63) ^ 32;
            // XOR to swith the term in last level
            x_last_level[pidx]*=x_level[new_pidx];
         }
         __syncthreads();
      }
      if(n_th > 128)
      {
         x_level = x_last_level;
         x_last_level -= 256;
         if(pidx < 256)
         {
            int new_pidx = (pidx & 127) ^ 64;
            // XOR to swith the term in last level
            x_last_level[pidx]*=x_level[new_pidx];
         }
         __syncthreads();
      }
      if(n_th > 256)
      {
         x_level = x_last_level;
         x_last_level -= 512;
         if(pidx < 512)
         {
            int new_pidx = (pidx & 255) ^ 128;
            // XOR to swith the term in last level
            x_last_level[pidx]*=x_level[new_pidx];
         }
         __syncthreads();
      }
      int new_tidx = pidx ^ (n_th/2); // XOR to swith the term in last level
      GT x_upper = x_last_level[new_tidx];
      if(pidx2 < n)
      {
         workspace_d[xidx1] = x1*x_upper;
         workspace_d[xidx2] = x0*x_upper;
      }
      else
      {
         workspace_d[xidx1] = x_upper;
      }
   }
}

// Sum block, mulitithread gpu sum unroll, for test
template <unsigned int n_th>
__global__ void eval_mon_tree_n_kernel
 ( GT* workspace_d, GT* x_d, GT* workspace_coef, int* mon_pos_start,
   unsigned short* pos_d, int n_mon)
   // int dim, int workspace_size_int){
{
   __shared__ GT x_sh[shmemsize];
   int BS = blockDim.x;
   int bidx = (gridDim.x*blockIdx.y+blockIdx.x)*BS;
   // int idx = bidx + threadIdx.x;
   int tidx = threadIdx.x;
   // int tidx2 = tidx + BS;

   int midx = tidx/n_th; // monomial index
   int midx_global =  midx + bidx/n_th;

   if(midx_global < n_mon)
   {
      // int sys_idx = blockIdx.z;
      // GT* x_d_tmp = x_d + sys_idx*dim;
      // GT* workspace_d_tmp = workspace_d + sys_idx*workspace_size_int;

      int pidx = tidx - midx*n_th; // thread index in monomial
      // int xidx0 = midx_global*(n_th*2+1);
      int xidx0 = mon_pos_start[midx_global];

      // int* pos = pos_d;// + BS*2*blockIdx.x;
      int n = pos_d[xidx0];

      GT* x_level = x_sh + midx*n_th*2;

      // Load to Shared Memory
      unsigned short* pos_tmp = pos_d+xidx0+1;
      GT* workspace_tmp = workspace_d+xidx0+1;

      int pos_idx = pidx;
      GT tmp = x_d[pos_tmp[pos_idx]];
      workspace_tmp[pos_idx] = GT(1,0);
      pos_idx += n_th;
      while(pos_idx < n)
      {
         workspace_tmp[pos_idx] = tmp;
         tmp *= x_d[pos_tmp[pos_idx]];
         pos_idx += n_th;
      }
      x_level[pidx] = tmp;

      // Up
      GT* x_last_level;

      if(n_th > 256)
      {
         __syncthreads();
         x_last_level = x_level; x_level += 512;
         if(pidx < 256)
         {
            x_level[pidx] = x_last_level[pidx] * x_last_level[pidx+256];
         }
      }
      if(n_th > 128)
      {
         __syncthreads();
         x_last_level = x_level; x_level += 256;
         if(pidx < 128)
         {
            x_level[pidx] = x_last_level[pidx] * x_last_level[pidx+128];
         }
      }
      if(n_th > 64)
      {
         __syncthreads();
         x_last_level = x_level; x_level += 128;
         if(pidx < 64)
         {
            x_level[pidx] = x_last_level[pidx] * x_last_level[pidx+64];
         }
      }
      if(n_th > 32)
      {
         __syncthreads();
         x_last_level = x_level; x_level += 64;
         if(pidx < 32)
         {
            x_level[pidx] = x_last_level[pidx] * x_last_level[pidx+32];
         }
      }
      if(n_th > 16)
      {
         x_last_level = x_level; x_level += 32;
         if(pidx < 16)
         {
            x_level[pidx] = x_last_level[pidx] * x_last_level[pidx+16];
         }
      }
      if(n_th > 8)
      {
         x_last_level = x_level; x_level += 16;
         if(pidx <  8)
         {
            x_level[pidx] = x_last_level[pidx] * x_last_level[pidx+ 8];
         }
      }
      if(n_th > 4)
      {
         x_last_level = x_level; x_level +=  8;
         if(pidx <  4)
         {
            x_level[pidx] = x_last_level[pidx] * x_last_level[pidx+ 4];
         }
      }
      if(n_th > 2)
      {
         x_last_level = x_level; x_level +=  4;
         if(pidx <  2)
         {
            x_level[pidx] = x_last_level[pidx] * x_last_level[pidx+ 2];
         }
      }
      if(n_th > 1)
      {
         x_last_level = x_level; x_level +=  2;
         if(pidx == 0)
         { 
            x_level[0] = x_last_level[0] * x_last_level[1];
         }
      }
      if(pidx < 3)
      {
         x_last_level[pidx] *= workspace_coef[midx_global];
      }
      // Common Factor to be multiplied
      if(pidx == 0)
      {
         workspace_d[xidx0] = x_level[0];
      }
      // Down
      if(n_th > 2)
      {
         x_level = x_last_level;
         x_last_level -= 4;
         if(pidx <  4)
         {
            int new_pidx = (pidx&1) ^ 1; // XOR to swith the term in last level
            x_last_level[pidx]*=x_level[new_pidx];
         }
      }
      if(n_th > 4)
      {
         x_level = x_last_level;
         x_last_level -= 8;
         if(pidx <  8)
         {
            int new_pidx = (pidx&3) ^ 2; // XOR to swith the term in last level
            x_last_level[pidx]*=x_level[new_pidx];
         }
      }
      if(n_th > 8)
      {
         x_level = x_last_level;
         x_last_level -= 16;
         if(pidx <  16)
         {
            int new_pidx = (pidx & 7) ^ 4;
            // XOR to swith the term in last level
            x_last_level[pidx]*=x_level[new_pidx];
         }
      }
      if(n_th > 16)
      {
         x_level = x_last_level;
         x_last_level -= 32;
         if(pidx < 32)
         {
            int new_pidx = (pidx & 15) ^ 8;
            // XOR to swith the term in last level
            x_last_level[pidx]*=x_level[new_pidx];
         }
         __syncthreads();
      }
      if(n_th > 32)
      {
         x_level = x_last_level;
         x_last_level -= 64;
         if(pidx < 64)
         {
            int new_pidx = (pidx & 31) ^ 16;
            // XOR to swith the term in last level
            x_last_level[pidx]*=x_level[new_pidx];
         }
         __syncthreads();
      }
      if(n_th > 64)
      {
         x_level = x_last_level;
         x_last_level -= 128;
         if(pidx < 128)
         {
            int new_pidx = (pidx & 63) ^ 32;
            // XOR to swith the term in last level
            x_last_level[pidx]*=x_level[new_pidx];
         }
         __syncthreads();
      }
      if(n_th > 128)
      {
         x_level = x_last_level;
         x_last_level -= 256;
         if(pidx < 256)
         {
            int new_pidx = (pidx & 127) ^ 64;
            // XOR to swith the term in last level
            x_last_level[pidx]*=x_level[new_pidx];
         }
         __syncthreads();
      }
      if(n_th > 256)
      {
         x_level = x_last_level;
         x_last_level -= 512;
         if(pidx < 512)
         {
            int new_pidx = (pidx & 255) ^ 128;
            // XOR to swith the term in last level
            x_last_level[pidx]*=x_level[new_pidx];
         }
         __syncthreads();
      }
      // Copy Derivative
      // bidx *= 2;
      // tidx += 1; // for test, keep tidx in memory
      int new_tidx = pidx ^ (n_th/2); // XOR to swith the term in last level
      tmp = x_last_level[new_tidx];

      pos_idx = pidx + (n-1)/n_th*n_th;
      // int tmp_idx = pos_idx;

      if(pos_idx < n)
      {
         workspace_tmp[pos_idx] *= tmp;
         tmp *= x_d[pos_tmp[pos_idx]];
      }
      pos_idx -= n_th;
      // GT tmp = x_d[pos_tmp[pos_idx]];
      // workspace_tmp[pos_idx] = tmp;

      while(pos_idx >= n_th)
      {
         workspace_tmp[pos_idx] *= tmp;
         tmp *= x_d[pos_tmp[pos_idx]];
         pos_idx -= n_th;
      }
      workspace_tmp[pos_idx] *= tmp;
   }
}

// Mon evalutation and differentiation on GPU, level 1 x0*x1, x1, x0
__global__ void eval_mon_tree_2_kernel
 ( GT* workspace_d, GT* x_d, GT* workspace_coef, int* mon_pos_start,
   unsigned short* pos_d, int n_mon )
{
   int idx = blockIdx.x*blockDim.x + threadIdx.x;

   if(idx < n_mon)
   {
      // int sys_idx = blockIdx.z;
      // GT* x_d_tmp = x_d + sys_idx*dim;
      // GT* workspace_d_tmp = workspace_d + sys_idx*workspace_size_int;

      // int xidx = idx*3;
      int xidx = mon_pos_start[idx];

      // Load to Shared Memory
      GT x0, x1, tmp;

      x1  = x_d[pos_d[xidx+2]];
      tmp = workspace_coef[idx];
      x0  = x_d[pos_d[xidx+1]]*tmp;

      workspace_d[xidx+2] = x0;
      workspace_d[xidx]   = x0*x1;
      workspace_d[xidx+1] = x1*tmp;
   }
}

// Sum block, mulitithread gpu sum unroll, for test
__global__ void eval_mon_tree_4_kernel
 ( GT* workspace_d, GT* x_d, GT* workspace_coef, int* mon_pos_start,
   unsigned short* pos_d, int n_mon )
{
   // GT* x_d, unsigned short* pos_d, GT* workspace_d, int n_mon,
   //                         int dim, int workspace_size_int){
   __shared__ GT x_sh[shmemsize];
   int BS = blockDim.x;
   int bidx = (gridDim.x*blockIdx.y+blockIdx.x)*BS;
   // int idx = bidx + threadIdx.x;
   int tidx = threadIdx.x;
   // int tidx2 = tidx + BS;

   int midx = tidx/2; // monomial index
   int midx_global =  midx + bidx/2;

   if(midx_global < n_mon)
   {
      // int sys_idx = blockIdx.z;
      // GT* x_d_tmp = x_d + sys_idx*dim;
      // GT* workspace_d_tmp = workspace_d + sys_idx*workspace_size_int;

      int pidx = tidx - midx*2; // thread index in monomial
      int pidx2 = pidx + 2;
      // int xidx0 = midx_global*(n_th*2+1);
      int xidx0 = mon_pos_start[midx_global];
      int xidx1 = xidx0+ pidx+1; // pos index 1
      int xidx2 = xidx1 + 2;  // pos index 2

      // int* pos = pos_d; // + BS*2*blockIdx.x;
      int n = pos_d[xidx0];

      GT* x_level = x_sh + midx*3;

      // Load to Shared Memory
      GT x0, x1;
      x0 = x_d[pos_d[xidx1]];
      if(pidx2 < n)
      {
         x1 = x_d[pos_d[xidx2]];
         x_level[pidx] = x0*x1;
      }
      else
      {
         x_level[pidx] = x0;
      }
      if(pidx == 0)
      {
         x_level[2] = x_level[0] * x_level[1];
      }
      GT coef_tmp = workspace_coef[midx_global];

      x_level[pidx] *= coef_tmp;

      if(pidx == 0)
      {
         workspace_d[xidx0] = x_level[2]*coef_tmp;
      }

      // Copy Derivative
      int new_tidx = pidx ^ 1; // XOR to swith the term in last level
      GT x_upper = x_level[new_tidx];
      if(pidx2 < n)
      {
         workspace_d[xidx1] = x1*x_upper;
         workspace_d[xidx2] = x0*x_upper;
      }
      else
      {
         workspace_d[xidx1] = x_upper;
      }
   }
}

// Sum block, mulitithread gpu sum unroll, for test
__global__ void eval_mon_tree_4n_kernel
 ( GT* workspace_d, GT* x_d, GT* workspace_coef, int* mon_pos_start,
   unsigned short* pos_d, int n_mon )
{
   // GT* x_d, unsigned short* pos_d, GT* workspace_d, int n_mon,
   //                        int dim, int workspace_size_int){
   __shared__ GT x_sh[shmemsize];
   int BS = blockDim.x;
   int bidx = (gridDim.x*blockIdx.y+blockIdx.x)*BS;
   // int idx = bidx + threadIdx.x;
   int tidx = threadIdx.x;
   // int tidx2 = tidx + BS;

   int midx = tidx/2; // monomial index
   int midx_global =  midx + bidx/2;

   if(midx_global < n_mon)
   {
      //int sys_idx = blockIdx.z;
      //GT* x_d_tmp = x_d + sys_idx*dim;
      //GT* workspace_d_tmp = workspace_d + sys_idx*workspace_size_int;

      int pidx = tidx - midx*2; // thread index in monomial
      //int pidx2 = pidx + 2;
      //int xidx0 = midx_global*(n_th*2+1);
      int xidx0 = mon_pos_start[midx_global];
      //int xidx1 = xidx0+ pidx+1; // pos index 1
      //int xidx2 = xidx1 + 2;  // pos index 2

      //int* pos = pos_d;// + BS*2*blockIdx.x;
      int n = pos_d[xidx0];

      GT* x_level = x_sh + midx*3;

      // Load to Shared Memory
      unsigned short* pos_tmp = pos_d+xidx0+1;
      GT* workspace_tmp = workspace_d+xidx0+1;

      int pos_idx = pidx;
      GT tmp = x_d[pos_tmp[pos_idx]];
      workspace_tmp[pos_idx] = GT(1,0);
      pos_idx += 2;
      while(pos_idx < n)
      {
         workspace_tmp[pos_idx] = tmp;
         tmp *= x_d[pos_tmp[pos_idx]];
         pos_idx += 2;
      }
      x_level[pidx] = tmp;

      if(pidx == 0)
      {
         x_level[2] = x_level[0] * x_level[1];
      }
      GT coef_tmp = workspace_coef[midx_global];
      x_level[pidx] *= coef_tmp;

      if(pidx == 0)
      {
         workspace_d[xidx0] = x_level[2]*coef_tmp;
      }
      // Copy Derivative
      int new_tidx = pidx ^ 1; // XOR to swith the term in last level
      tmp = x_level[new_tidx];
      pos_idx = pidx + (n-1)/2*2;
      // int tmp_idx = pos_idx;
      if(pos_idx < n)
      {
         workspace_tmp[pos_idx] *= tmp;
         tmp *= x_d[pos_tmp[pos_idx]];
      }
      pos_idx -= 2;
      // GT tmp = x_d[pos_tmp[pos_idx]];
      // workspace_tmp[pos_idx] = tmp;

      while(pos_idx >= 2)
      {
         workspace_tmp[pos_idx] *= tmp;
         tmp *= x_d[pos_tmp[pos_idx]];
         pos_idx -= 2;
      }
      workspace_tmp[pos_idx] *= tmp;
   }
}

