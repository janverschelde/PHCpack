// Monomial evaluation and differentiation on GPU
#include "eval_mon_single.cu"

__global__ void eval_mon_seq_kernel
 ( GT* workspace_mon, GT* x, GT* workspace_coef,
   int* mon_pos_start, unsigned short* mon_pos, int n_mon) 
{
   int idx = (gridDim.x*blockIdx.y+blockIdx.x)*blockDim.x + threadIdx.x;

   // int tidx = threadIdx.x;
   if(idx < n_mon)
   {
      int tmp_start = mon_pos_start[idx];
      GT* deri = workspace_mon + tmp_start;
      unsigned short* pos = mon_pos + tmp_start;

      int n_var = pos[0];

      GT tmp = x[pos[1]];

      GT* deri_tmp = deri + 1;
      deri_tmp[1] = tmp;

      for(int i=2; i<n_var; i++)
      {
         tmp *= x[pos[i]];
         deri_tmp[i] = tmp;
      }
      tmp = workspace_coef[idx];

      for(int i=n_var; i>1; i--) 
      {
         deri[i] *= tmp;
         tmp *= x[pos[i]];
      }
      deri[1] = tmp;
      deri[0] = x[pos[1]]*tmp;
   }
}

__global__ void eval_mon_seq_with_base_kernel
 ( GT* workspace_mon, GT* x, GT* workspace_coef, int* mon_pos_start,
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

      int n_var = pos[0];

      GT tmp = x[pos[1]];

      GT* deri_tmp = deri + 1;
      deri_tmp[1] = tmp;

      for(int i=2; i<n_var; i++)
      {
         tmp *= x[pos[i]];
         deri_tmp[i] = tmp;
      }
      tmp = workspace_coef[idx];

      for(int i=n_var; i>1; i--) 
      {
         deri[i] *= tmp*GT(exp[i],0.0);
         tmp *= x[pos[i]];
      }
      deri[1] = tmp*GT(exp[1],0.0);
      deri[0] = x[pos[1]]*tmp;
   }
}

void eval_mon_seq_without_base ( GPUWorkspace& workspace, const GPUInst& inst )
{
   eval_mon_single_kernel<<<inst.mon_level_grid[0], inst.mon_level0_BS>>>
      (workspace.mon, workspace.x, workspace.coef, inst.mon_pos_start,
       inst.mon_pos, inst.n_mon_level[0]);

   eval_mon_seq_kernel<<<inst.mon_global_grid, inst.mon_global_BS>>>
      (workspace.mon, workspace.x, workspace.coef + inst.n_mon_level[0],
       inst.mon_pos_start + inst.n_mon_level[0], inst.mon_pos,
       inst.n_mon_global);
}

void eval_mon_seq_with_base ( GPUWorkspace& workspace, const GPUInst& inst )
{
   eval_mon_single_with_base_kernel<<<inst.mon_level_grid[0],
                                      inst.mon_level0_BS>>>
      (workspace.mon, workspace.x, workspace.coef, inst.mon_pos_start,
       inst.mon_pos,  inst.mon_exp, inst.n_mon_level[0]);

   eval_mon_seq_with_base_kernel<<<inst.mon_global_grid, inst.mon_global_BS>>>
      (workspace.mon, workspace.x, workspace.coef + inst.n_mon_level[0], 
       inst.mon_pos_start + inst.n_mon_level[0], inst.mon_pos, inst.mon_exp,
       inst.n_mon_global);
}

void eval_mon_seq ( GPUWorkspace& workspace, const GPUInst& inst )
{
   if(workspace.deg_table == NULL)
   {
      eval_mon_seq_without_base(workspace, inst);
   }
   else
   {
      eval_mon_seq_with_base(workspace, inst);
   }
}

__global__ void eval_mon_seq_kernel
 ( GT* workspace_mon, GT* x, GT* workspace_coef, int* mon_pos_start,
   unsigned short* mon_pos, int n_mon, int workspace_size )
{
   int idx = (gridDim.x*blockIdx.y+blockIdx.x)*blockDim.x + threadIdx.x;

   int sys_idx = blockIdx.z;
   workspace_mon += workspace_size*sys_idx;
   x += workspace_size*sys_idx;
   workspace_coef += workspace_size*sys_idx;

   // int tidx = threadIdx.x;
   if(idx < n_mon)
   {
      int tmp_start = mon_pos_start[idx];
      GT* deri = workspace_mon + tmp_start;
      unsigned short* pos = mon_pos + tmp_start;

      int n_var = pos[0];

      GT tmp = x[pos[1]];

      GT* deri_tmp = deri + 1;
      deri_tmp[1] = tmp;

      for(int i=2; i<n_var; i++)
      {
         tmp *= x[pos[i]];
         deri_tmp[i] = tmp;
      }
      tmp = workspace_coef[idx];

      for(int i=n_var; i>1; i--)
      {
         deri[i] *= tmp;
         tmp *= x[pos[i]];
      }
      deri[1] = tmp;
      deri[0] = x[pos[1]]*tmp;
   }
}

void eval_mon_seq_mult ( GPUWorkspace& workspace, const GPUInst& inst )
{
   eval_mon_single_kernel<<<inst.mon_level_grid[0], inst.mon_level0_BS>>>
      (workspace.mon, workspace.x, workspace.coef, inst.mon_pos_start,
       inst.mon_pos, inst.n_mon_level[0], workspace.workspace_size);

   eval_mon_seq_kernel<<<inst.mon_global_grid, inst.mon_global_BS>>>
      (workspace.mon, workspace.x, workspace.coef + inst.n_mon_level[0],
       inst.mon_pos_start + inst.n_mon_level[0], inst.mon_pos,
       inst.n_mon_global, workspace.workspace_size);
}
