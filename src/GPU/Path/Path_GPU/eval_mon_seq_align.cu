#include "eval_mon_single.cu"

__global__ void eval_mon_seq_align_kernel
 ( GT* workspace_mon, GT* x, GT* workspace_coef,
   int* mon_pos_start_block, unsigned short* mon_pos_block, int n_mon );

void eval_mon_seq_align ( GPUWorkspace& workspace, const GPUInst& inst )
{
   dim3 mon_single_grid = get_grid(inst.n_mon_single,inst.mon_level0_BS,1);

   eval_mon_single_kernel<<<mon_single_grid, inst.mon_level0_BS>>>
      (workspace.mon, workspace.x, workspace.coef,
       inst.mon_single_pos_block, inst.n_mon_single);

   dim3 mon_block_grid = get_grid(inst.n_mon_block, BS_Mon_Align, 1);

   eval_mon_seq_align_kernel<<<mon_block_grid, BS_Mon_Align>>>
      (workspace.mon, workspace.x, workspace.coef + inst.n_mon_single,
       inst.mon_pos_start_block, inst.mon_pos_block, inst.n_mon_block);
}

__global__ void eval_mon_seq_align_kernel
 ( GT* workspace_mon, GT* x, GT* workspace_coef,
   int* mon_pos_start_block, unsigned short* mon_pos_block, int n_mon )
{
   int idx = (gridDim.x*blockIdx.y+blockIdx.x)*blockDim.x + threadIdx.x;
   int widx = idx/32;
   int wtidx = idx - widx*32;

   if(idx < n_mon)
   {
      int tmp_block_start = mon_pos_start_block[widx];
      GT* deri = workspace_mon + tmp_block_start;
      unsigned short* pos = mon_pos_block + tmp_block_start;

      int n_var = pos[wtidx];

      GT tmp = x[pos[32+wtidx]];

      GT* deri_tmp = deri + 32;
      deri_tmp[32+wtidx] = tmp;

      for(int i=2; i<n_var; i++)
      {
         tmp *= x[pos[i*32+wtidx]];
         deri_tmp[i*32+wtidx] = tmp;
      }
      tmp = workspace_coef[idx];

      for(int i=n_var; i>1; i--)
      {
         deri[i*32+wtidx] *= tmp;
         tmp *= x[pos[i*32+wtidx]];
      }
      deri[32+wtidx] = tmp;
      deri[wtidx] = x[pos[32+wtidx]]*tmp;
   }
}

__global__ void eval_mon_seq_align_kernel
 ( GT* workspace_mon, GT* x, GT* workspace_coef, int* mon_pos_start_block,
   unsigned short* mon_pos_block, int n_mon, int workspace_size,
   int* x_t_idx, int dim );

__global__ void eval_mon_seq_align_kernel
 ( GT* workspace_mon, GT* x, GT* workspace_coef, int* mon_pos_start_block,
   unsigned short* mon_pos_block, int n_mon, int workspace_size, 
   int* x_t_idx, int dim, int* path_idx_mult );

void eval_mon_seq_align2 ( GPUWorkspace& workspace, const GPUInst& inst )
{
   int n_path = workspace.n_path_continuous;
   std::cout << "n_mon_single = " << inst.n_mon_single << std::endl;
   dim3 mon_single_grid
      = get_grid(inst.n_mon_single,inst.mon_level0_BS,n_path);

   eval_mon_single_kernel<<<mon_single_grid, inst.mon_level0_BS>>>
      (workspace.mon, workspace.x_array, workspace.coef, inst.mon_pos_start,
       inst.mon_pos_block, inst.n_mon_single, workspace.workspace_size,
       workspace.x_t_idx_mult, inst.dim, workspace.path_idx);

   dim3 mon_block_grid = get_grid(inst.n_mon_block, BS_Mon_Align, n_path);
   // int NB_mon = (inst.n_mon_block-1)/BS_Mon_Align + 1;
   eval_mon_seq_align_kernel<<<mon_block_grid, BS_Mon_Align>>>
     (workspace.mon, workspace.x_array, workspace.coef + inst.n_mon_single,
      inst.mon_pos_start_block, inst.mon_pos_block, inst.n_mon_block,
      workspace.workspace_size, workspace.x_t_idx_mult, inst.dim,
      workspace.path_idx);

   /*
     eval_mon_seq_align_block_kernel<<<inst.NB_mon_block, inst.BS_mon_block>>>
        (workspace.mon+inst.n_mon_level[0]*2, workspace.x,
         workspace.coef + inst.n_mon_level[0], inst.mon_pos_start_block,
         inst.mon_pos_block, inst.n_mon_block);
    */
}

__global__ void eval_mon_seq_align_kernel
 ( GT* workspace_mon, GT* x, GT* workspace_coef, int* mon_pos_start_block,
   unsigned short* mon_pos_block, int n_mon, int workspace_size, 
   int* x_t_idx, int dim, int* path_idx_mult )
{
   int idx = (gridDim.x*blockIdx.y+blockIdx.x)*blockDim.x + threadIdx.x;
   int widx = idx/32;
   int wtidx = idx - widx*32;

   int path_idx = path_idx_mult[blockIdx.z];
   workspace_mon += workspace_size*path_idx;
   x += workspace_size*path_idx+x_t_idx[path_idx]*dim;
   workspace_coef += workspace_size*path_idx;

   if(idx < n_mon)
   {
      int tmp_block_start = mon_pos_start_block[widx];
      GT* deri = workspace_mon + tmp_block_start;
      unsigned short* pos = mon_pos_block + tmp_block_start;

      int n_var = pos[wtidx];

      GT tmp = x[pos[32+wtidx]];

      GT* deri_tmp = deri + 32;
      deri_tmp[32+wtidx] = tmp;

      for(int i=2; i<n_var; i++)
      {
         tmp *= x[pos[i*32+wtidx]];
         deri_tmp[i*32+wtidx] = tmp;
      }
      tmp = workspace_coef[idx];

      for(int i=n_var; i>1; i--)
      {
         deri[i*32+wtidx] *= tmp;
         tmp *= x[pos[i*32+wtidx]];
      }
      deri[32+wtidx] = tmp;
      // deri[wtidx] = x[pos[32+wtidx]]*tmp;
      deri[wtidx] = x[pos[32+wtidx]]*tmp;
   }
}

__global__ void eval_mon_seq_align_kernel
 ( GT* workspace_mon, GT* x, GT* workspace_coef, int* mon_pos_start_block,
   unsigned short* mon_pos_block, int n_mon, int workspace_size, 
   int* x_t_idx, int dim )
{
   int idx = (gridDim.x*blockIdx.y+blockIdx.x)*blockDim.x + threadIdx.x;
   int widx = idx/32;
   int wtidx = idx - widx*32;

   int path_idx = blockIdx.z;
   workspace_mon += workspace_size*path_idx;
   x += workspace_size*path_idx+x_t_idx[path_idx]*dim;
   workspace_coef += workspace_size*path_idx;

   if(idx < n_mon)
   {
      int tmp_block_start = mon_pos_start_block[widx];
      GT* deri = workspace_mon + tmp_block_start;
      unsigned short* pos = mon_pos_block + tmp_block_start;

      int n_var = pos[wtidx];

      GT tmp = x[pos[32+wtidx]];

      GT* deri_tmp = deri + 32;
      deri_tmp[32+wtidx] = tmp;

      for(int i=2; i<n_var; i++)
      {
         tmp *= x[pos[i*32+wtidx]];
         deri_tmp[i*32+wtidx] = tmp;
      }
      tmp = workspace_coef[idx];

      for(int i=n_var; i>1; i--) 
      {
         deri[i*32+wtidx] *= tmp;
         tmp *= x[pos[i*32+wtidx]];
      }
      deri[32+wtidx] = tmp;
      // deri[wtidx] = x[pos[32+wtidx]]*tmp;
      deri[wtidx] = x[pos[32+wtidx]]*tmp;
   }
}

__global__ void eval_mon_seq_align_kernel
 ( GT* workspace_mon, GT* x, GT* workspace_coef, int* mon_pos_start_block,
   unsigned short* mon_pos_block, int n_mon, int workspace_size )
{
   int idx = (gridDim.x*blockIdx.y+blockIdx.x)*blockDim.x + threadIdx.x;
   int widx = idx/32;
   int wtidx = idx - widx*32;

   int path_idx = blockIdx.z;
   workspace_mon += workspace_size*path_idx;
   x += workspace_size*path_idx;
   workspace_coef += workspace_size*path_idx;

   if(idx < n_mon)
   {
      int tmp_block_start = mon_pos_start_block[widx];
      GT* deri = workspace_mon + tmp_block_start;
      unsigned short* pos = mon_pos_block + tmp_block_start;

      int n_var = pos[wtidx];

      GT tmp = x[pos[32+wtidx]];

      GT* deri_tmp = deri + 32;
      deri_tmp[32+wtidx] = tmp;

      for(int i=2; i<n_var; i++)
      {
         tmp *= x[pos[i*32+wtidx]];
         deri_tmp[i*32+wtidx] = tmp;
      }
      tmp = workspace_coef[idx];

      for(int i=n_var; i>1; i--)
      {
         deri[i*32+wtidx] *= tmp;
         tmp *= x[pos[i*32+wtidx]];
      }
      deri[32+wtidx] = tmp;
      deri[wtidx] = x[pos[32+wtidx]]*tmp;
   }
}
