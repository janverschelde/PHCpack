#ifndef __PATH_GPU_EVAL_MON_SINGLE_CU_
#define __PATH_GPU_EVAL_MON_SINGLE_CU_

// Mon evalutaion and differentiation on GPU
__global__ void eval_mon_single_kernel
 ( GT* workspace_mon, GT* x, GT*workspace_coef, int* mon_pos_start,
   unsigned short* mon_pos, int n_mon) 
{
   int idx = (gridDim.x*blockIdx.y+blockIdx.x)*blockDim.x + threadIdx.x;

   // int tidx = threadIdx.x;
   if(idx < n_mon)
   {
      int tmp_start = mon_pos_start[idx];
      GT* deri = workspace_mon + tmp_start;
      unsigned short* pos = mon_pos + tmp_start;

      GT tmp = workspace_coef[idx];
      deri[1] = tmp;
      deri[0] = x[pos[1]]*tmp;
   }
}

// Mon evalutaion and differentiation on GPU
__global__ void eval_mon_single_kernel
 ( GT* workspace_mon, GT* x, GT*workspace_coef,
   unsigned short* mon_pos, int n_mon )
{
   int idx = (gridDim.x*blockIdx.y+blockIdx.x)*blockDim.x + threadIdx.x;

   // int tidx = threadIdx.x;
   if(idx < n_mon) 
   {
      int tmp_start = idx*2;
      GT* deri = workspace_mon + tmp_start;
      unsigned short* pos = mon_pos + tmp_start;

      GT tmp = workspace_coef[idx];
      deri[1] = tmp;
      deri[0] = x[pos[1]]*tmp;
   }
}

__global__ void eval_mon_single_with_base_kernel
 ( GT* workspace_mon, GT* x, GT*workspace_coef, int* mon_pos_start,
   unsigned short* mon_pos, unsigned short* mon_exp, int n_mon )
{
   int idx = (gridDim.x*blockIdx.y+blockIdx.x)*blockDim.x + threadIdx.x;

   // int tidx = threadIdx.x;
   if(idx < n_mon)
   {
      int tmp_start = mon_pos_start[idx];
      GT* deri = workspace_mon + tmp_start;
      unsigned short* pos = mon_pos + tmp_start;
      unsigned short* exp = mon_exp + tmp_start;

      GT tmp = workspace_coef[idx];
      deri[1] = tmp*GT(exp[1],0.0);
      deri[0] = x[pos[1]]*tmp;
   }
}

// Mon evalutaion and differentiation on GPU
__global__ void eval_mon_single_kernel
 ( GT* workspace_mon, GT* x, GT*workspace_coef, int* mon_pos_start,
   unsigned short* mon_pos, int n_mon, int workspace_size )
{
   int idx = (gridDim.x*blockIdx.y+blockIdx.x)*blockDim.x + threadIdx.x;

   int path_idx = blockIdx.z;
   workspace_mon += workspace_size*path_idx;
   x += workspace_size*path_idx;
   workspace_coef += workspace_size*path_idx;

   // int tidx = threadIdx.x;
   if(idx < n_mon)
   {
      int tmp_start = mon_pos_start[idx];
      GT* deri = workspace_mon + tmp_start;
      unsigned short* pos = mon_pos + tmp_start;

      GT tmp = workspace_coef[idx];
      deri[1] = tmp;
      deri[0] = x[pos[1]]*tmp;
   }
}

// Mon evalutaion and differentiation on GPU
__global__ void eval_mon_single_kernel
 ( GT* workspace_mon, GT* x, GT*workspace_coef, int* mon_pos_start,
   unsigned short* mon_pos, int n_mon, int workspace_size, int* x_t_idx,
   int dim )
{
   int idx = (gridDim.x*blockIdx.y+blockIdx.x)*blockDim.x + threadIdx.x;

   int path_idx = blockIdx.z;
   workspace_mon += workspace_size*path_idx;
   x += workspace_size*path_idx + x_t_idx[path_idx]*dim;
   workspace_coef += workspace_size*path_idx;

   // int tidx = threadIdx.x;
   if(idx < n_mon)
   {
      int tmp_start = mon_pos_start[idx];
      GT* deri = workspace_mon + tmp_start;
      unsigned short* pos = mon_pos + tmp_start;

      GT tmp = workspace_coef[idx];
      deri[1] = tmp;
      deri[0] = x[pos[1]]*tmp;
   }
}

// Mon evalutaion and differentiation on GPU
__global__ void eval_mon_single_kernel
 ( GT* workspace_mon, GT* x, GT*workspace_coef, int* mon_pos_start,
   unsigned short* mon_pos, int n_mon, int workspace_size, int* x_t_idx,
   int dim, int* path_idx_mult )
{
   int idx = (gridDim.x*blockIdx.y+blockIdx.x)*blockDim.x + threadIdx.x;

   int path_idx = path_idx_mult[blockIdx.z];
   workspace_mon += workspace_size*path_idx;
   x += workspace_size*path_idx + x_t_idx[path_idx]*dim;
   workspace_coef += workspace_size*path_idx;

   // int tidx = threadIdx.x;
   if(idx < n_mon)
   {
      int tmp_start = 2*idx; //mon_pos_start[idx];
      GT* deri = workspace_mon + tmp_start;
      unsigned short* pos = mon_pos + tmp_start;

      GT tmp = workspace_coef[idx];
      deri[1] = tmp;
      deri[0] = x[pos[1]]*tmp;
   }
}

#endif /*__PATH_GPU_EVAL_MON_SINGLE_CU_*/
