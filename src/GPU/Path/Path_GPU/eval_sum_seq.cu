#ifndef __PATH_GPU_EVAL_SUM_SEQ_CU_
#define __PATH_GPU_EVAL_SUM_SEQ_CU_

// Mon evaluation and differentiation on GPU

__global__ void eval_sum_seq_kernel
 ( GT* workspace_matrix, GT* workspace_sum, 
   int* sum_pos, int* sum_pos_start, int n_sum )
{
   int idx = (gridDim.x*blockIdx.y+blockIdx.x)*blockDim.x + threadIdx.x;

   if(idx < n_sum)
   {
      int* pos = sum_pos + sum_pos_start[idx];
      int n_var = *pos++;

      GT tmp = workspace_sum[*pos++];

      for(int i=1; i<n_var; i++)
      {
         tmp += workspace_sum[*pos++];
      }
      workspace_matrix[*pos] = tmp;
   }
}

// Mon evaluation and differentiation on GPU

__global__ void eval_sum_seq_kernel
 ( GT* workspace_matrix, GT* workspace_sum, int* sum_pos,
   int* sum_pos_start, int n_sum, int workspace_size)
{
   int idx = (gridDim.x*blockIdx.y+blockIdx.x)*blockDim.x + threadIdx.x;

   int path_idx = blockIdx.z;
   workspace_matrix += workspace_size*path_idx;
   workspace_sum += workspace_size*path_idx;

   if(idx < n_sum) 
   {
      int* pos = sum_pos + sum_pos_start[idx];
      int n_var = *pos++;

      GT tmp = workspace_sum[*pos++];

      for(int i=1; i<n_var; i++)
      {
         tmp += workspace_sum[*pos++];
      }
      workspace_matrix[*pos] = tmp;
   }
}

// Mon evaluation and differentiation on GPU

__global__ void eval_sum_seq_kernel
 ( GT* workspace_matrix, GT* workspace_sum, int* sum_pos,
   int* sum_pos_start, int n_sum, int workspace_size, int* path_idx_mult )
{
   int idx = (gridDim.x*blockIdx.y+blockIdx.x)*blockDim.x + threadIdx.x;

   int path_idx = path_idx_mult[blockIdx.z];
   workspace_matrix += workspace_size*path_idx;
   workspace_sum += workspace_size*path_idx;

   if(idx < n_sum)
   {
      int* pos = sum_pos + sum_pos_start[idx];
      int n_var = *pos++;

      GT tmp = workspace_sum[*pos++];

      for(int i=1; i<n_var; i++)
      {
         tmp += workspace_sum[*pos++];
      }
      workspace_matrix[*pos] = tmp;
   }
}

void eval_sum_seq ( GPUWorkspace& workspace, const GPUInst& inst )
{
   dim3 sum_grid = get_grid(inst.n_sum,inst.sum_BS,workspace.n_path_continuous);
   //eval_sum_seq_kernel<<<sum_grid, inst.sum_BS>>>(workspace.matrix, \
   //   workspace.sum, inst.sum_pos, inst.sum_pos_start, inst.n_sum, \
   // workspace.workspace_size, workspace.path_idx);

   eval_sum_seq_kernel<<<sum_grid, inst.sum_BS>>>
      (workspace.matrix, workspace.sum, inst.sum_pos,
       inst.sum_pos_start, inst.n_sum);
}

#endif /*__PATH_GPU_EVAL_SUM_SEQ_CU_*/
