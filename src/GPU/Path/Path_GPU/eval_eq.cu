// #include "eval_mon_single.cu"

__global__ void eval_coef_mult_kernel
 ( GT* workspace_coef, const GT* coef_orig, int n_path,
   GT* t, GT* one_minor_t )
{
   // __shared__ GT div_diff_sh[shmemsize];
   int bidx = (gridDim.x*blockIdx.y+blockIdx.x)*blockDim.x;
   int tidx = threadIdx.x;
   int path_idx = bidx + tidx;

   if(path_idx<n_path)
   {
      int coef_idx = blockIdx.z;
      workspace_coef += coef_idx*n_path;
      t += path_idx;
      one_minor_t += path_idx;

      /*int path_idx = blockIdx.z;
        x_predictor += path_idx*np_predictor*dim;
        t_predictor += path_idx*np_predictor;
        x_new += path_idx*dim;*/

      // workspace_coef[idx] = coef_orig[idx];
      // XXX align coef later (*t)*coef_orig[idx] 
      // + (*one_minor_t)*coef_orig[idx+n_coef]
      workspace_coef[path_idx] = (*t)*coef_orig[2*coef_idx]
                               + (*one_minor_t)*coef_orig[2*coef_idx+1];
   }
}

__global__ void eval_coef_mult_kernel
 ( GT* workspace_coef, const GT* coef_orig, int n_path,
   GT* t, GT* one_minor_t, int workspace_size, int* x_t_idx )
{
   // __shared__ GT div_diff_sh[shmemsize];
   int bidx = (gridDim.x*blockIdx.y+blockIdx.x)*blockDim.x;
   int tidx = threadIdx.x;
   int path_idx = bidx + tidx;

   if(path_idx<n_path)
   {
      int coef_idx = blockIdx.z;
      workspace_coef += coef_idx*n_path;
      t += path_idx*workspace_size + x_t_idx[path_idx];
      one_minor_t += path_idx*workspace_size;

      /*int path_idx = blockIdx.z;
        x_predictor += path_idx*np_predictor*dim;
        t_predictor += path_idx*np_predictor;
        x_new += path_idx*dim;*/

      // workspace_coef[idx] = coef_orig[idx];
      // XXX align coef later (*t)*coef_orig[idx]
      // + (*one_minor_t)*coef_orig[idx+n_coef]
      workspace_coef[path_idx] = (*t)*coef_orig[2*coef_idx]
                               + (*one_minor_t)*coef_orig[2*coef_idx+1];
   }
}

// Monomial evaluation and differentiation on GPU
__global__ void eval_mon_single_mult_kernel
 ( GT* workspace_mon, GT* x, GT*workspace_coef, int* mon_pos_start,
   unsigned short* mon_pos, int n_path, int workspace_size )
{
   int path_idx = (gridDim.x*blockIdx.y+blockIdx.x)*blockDim.x + threadIdx.x;

   if(path_idx < n_path)
   {
      int mon_idx = blockIdx.z;
      workspace_mon += path_idx;
      x += path_idx;
      workspace_coef += path_idx;

      int tmp_start = mon_pos_start[mon_idx];
      GT* deri = workspace_mon + tmp_start*n_path;
      unsigned short* pos = mon_pos + tmp_start;

      GT tmp = workspace_coef[mon_idx*n_path];
      deri[n_path] = tmp;
      deri[0] = x[pos[1]*n_path]*tmp;

      // deri[n_path] = GT(2.0,2.0);
      // deri[0] = GT(3.0,3.0);
   }
}

__global__ void eval_mon_seq_mult_kernel
 ( GT* workspace_mon, GT* x, GT* workspace_coef, int* mon_pos_start,
   unsigned short* mon_pos, int n_path, int workspace_size )
{
   __shared__ int pos[512];

   int t_idx = threadIdx.x;
   int BS = blockDim.x;
   int path_idx = (gridDim.x*blockIdx.y+blockIdx.x)*BS + t_idx;

   int mon_idx = blockIdx.z;
   workspace_mon += path_idx;
   x += path_idx;
   workspace_coef += path_idx;

   int tmp_start = mon_pos_start[mon_idx];
   unsigned short* pos_tmp = mon_pos + tmp_start;

   int n_var = pos_tmp[0];
   int rnd = (n_var-1)/BS+1;
   if(rnd>1)
   {
      for(int rnd_idx=0; rnd_idx<rnd-1; rnd_idx++)
      {
         pos[t_idx+1] = pos_tmp[t_idx+1]*n_path;
         t_idx += BS;
      }
   }
   if(t_idx<n_var)
   {
      pos[t_idx+1] = pos_tmp[t_idx+1]*n_path;
   }
   __syncthreads();

   if(path_idx < n_path)
   {
      GT* deri = workspace_mon + tmp_start*n_path;
      GT tmp = x[pos[1]];

      GT* deri_tmp = deri + n_path;
      deri_tmp[n_path] = tmp;

      for(int i=2; i<n_var; i++)
      {
         tmp *= x[pos[i]];
         deri_tmp[i*n_path] = tmp;
      }
      tmp = workspace_coef[mon_idx*n_path];

      for(int i=n_var; i>1; i--) 
      {
         deri[i*n_path] *= tmp;
         tmp *= x[pos[i]];
      }
      deri[n_path] = tmp;
      deri[0] = x[pos[1]]*tmp;
   }
}

__global__ void eval_eq_kernel
 ( GT* workspace_mon, GT* x_mult, GT* matrix_mult,
   int* mon_pos_start, unsigned short* mon_pos )
{
   __shared__ GT eq[11];
   __shared__ GT x[10];

   int t_idx = threadIdx.x;
   int BS = blockDim.x;
   // int path_idx = (gridDim.x*blockIdx.y+blockIdx.x)*BS + t_idx;

   int mon_idx = t_idx;
   // workspace_mon += path_idx;
   // x += path_idx;
   // workspace_coef += path_idx;

   int tmp_start = mon_pos_start[mon_idx];
   unsigned short* pos_tmp = mon_pos + tmp_start;

   x[t_idx] = x_mult[t_idx];
   eq[t_idx] = GT(0.0,0.0);

   // if(path_idx < n_path){

   GT* deri = workspace_mon;
   unsigned short* pos = mon_pos[tmp_start];

   GT tmp = x[pos[1]];

   GT* deri_tmp = deri + n_path;
   deri_tmp[n_path] = tmp;

   for(int i=2; i<n_var; i++)
   {
      tmp *= x[pos[i]];
      deri_tmp[i*n_path] = tmp;
   }

   tmp = workspace_coef[mon_idx*n_path];

   for(int i=n_var; i>1; i--)
   {
      eq[pos[i]]+= deri[i*n_path] * tmp;
      tmp *= x[pos[i]];
   }
   eq[pos[0]] += tmp;

   matrix_mult[t_idx] = eq[t_idx];
}

// Monomial evaluation and differentiation on GPU
__global__ void eval_sum_seq_mult_kernel
 ( GT* workspace_matrix, GT* workspace_sum, int* sum_pos, 
   int* sum_pos_start, int n_path )
{
   __shared__ int pos[512];
   int t_idx = threadIdx.x;
   int BS = blockDim.x;
   int path_idx = (gridDim.x*blockIdx.y+blockIdx.x)*BS + t_idx;
   int sum_idx = blockIdx.z;
   workspace_matrix += path_idx;
   workspace_sum += path_idx;

   int* pos_tmp = sum_pos + sum_pos_start[sum_idx];
   int n_var = *pos_tmp;

   int rnd = n_var/BS+1;
   if(rnd>1)
   {
      for(int rnd_idx=0; rnd_idx<rnd-1; rnd_idx++)
      {
         pos[t_idx] = pos_tmp[t_idx+1]*n_path;
         t_idx += BS;
      }
   }
   if(t_idx<=n_var)
   {
      pos[t_idx] = pos_tmp[t_idx+1]*n_path;
   }
   __syncthreads();

   if(path_idx < n_path)
   {
      GT tmp = workspace_sum[pos[0]];

      for(int i=1; i<n_var; i++) 
      {
         tmp += workspace_sum[pos[i]];
      }
      workspace_matrix[pos[n_var]] = tmp;
   }
}

// Monomial evaluation and differentiation on GPU
/*
  __global__ void eval_sum_seq_mult_transpose_kernel
  ( GT* workspace_matrix, GT* workspace_sum, int* sum_pos,
    int* sum_pos_start, int n_path, int n_sum, int sum_trunk )
 {
    __shared__ GT tile[16][17];
    int pos_output[16];

    int t_idx = threadIdx.x;

    int path_trunk_idx = gridDim.x*blockIdx.y+blockIdx.x;
    int path_trunk = blockDim.x;
    int path_idx = path_trunk_idx*path_trunk + t_idx;

    int sum_trunk_idx = blockIdx.z;
    int sum_start_idx = sum_trunk_idx*sum_trunk;
    int sum_end_idx = sum_start_idx + sum_trunk;
    if(n_sum < sum_end_idx)
    {
       sum_end_idx = n_sum;
    }
    if(path_idx < n_path)
    {
       workspace_matrix += path_idx;
       workspace_sum += path_idx;
       for(int sum_idx=sum_start_idx; sum_idx<sum_end_idx; sum_idx++)
       {
          int* pos = sum_pos + sum_pos_start[sum_idx];
          int n_var = *pos++;

          GT tmp = workspace_sum[(*pos++)*n_path];

          for(int i=1; i<n_var; i++) 
          {
             tmp += workspace_sum[(*pos++)*n_path];
          }
          workspace_matrix[(*pos)*n_path] = tmp;
          if(t_idx == 0)
          {
             pos_output[sum_idx-sum_start_idx] = *pos;
          }
       }
    }
    int n_sum_remain = sum_end_idx - sum_start_idx;
    if(t_idx < n_sum_remain)
    {
       // postponed instruction summations are ordered by number of terms,
       // need to be changed.
    }
 }
 */

__global__ void eval_mult_init_zero
 ( GT* x_array, GT* t_array, GT* one_minor_t, GT* x_mult,
   GT* t_mult, GT* one_minor_t_mult, int workspace_size,
   int* path_idx_mult, int* x_t_idx_mult, int n_path, int dim,
   GT* matrix_mult, int n_sum_zero, int* sum_zeros )
{
   __shared__ int zero_pos[512];

   int t_idx = threadIdx.x;
   int BS = blockDim.x;
   int eval_idx = (gridDim.x*blockIdx.y+blockIdx.x)*BS+t_idx;
   if(eval_idx < n_path)
   {
      int path_idx = path_idx_mult[eval_idx];

      GT* t = t_array + workspace_size*path_idx + x_t_idx_mult[path_idx];
      one_minor_t += workspace_size*path_idx;

      t_mult[eval_idx] = *t;
      one_minor_t_mult[eval_idx] = *one_minor_t;
      x_mult += eval_idx;
      GT* x = x_array + workspace_size*path_idx + x_t_idx_mult[path_idx]*dim;

      for(int var_idx=0; var_idx<dim; var_idx++)
      {
         x_mult[var_idx*n_path] = x[var_idx];
      }
   }
   int rnd = (n_sum_zero-1)/BS+1;
   if(rnd>1)
   {
      for(int rnd_idx=0; rnd_idx<rnd-1; rnd_idx++)
      {
         zero_pos[t_idx] = sum_zeros[t_idx]*n_path;
         t_idx += BS;
      }
   }
   if(t_idx<n_sum_zero)
   {
      zero_pos[t_idx] = sum_zeros[t_idx]*n_path;
   }
   __syncthreads();

   if(eval_idx < n_path)
   {
      matrix_mult += eval_idx;
      for(int zero_idx=0; zero_idx<n_sum_zero; zero_idx++)
      {
         int m_idx = zero_pos[zero_idx];
         matrix_mult[m_idx] = GT(0.0,0.0);
      }
   }
}

__global__ void eval_mult_init
 ( GT* x_array, GT* t_array, GT* one_minor_t, GT* x_mult, GT* t_mult,
   GT* one_minor_t_mult, int workspace_size, int* path_idx_mult,
   int* x_t_idx_mult, int n_path, int dim )
{
   int eval_idx = (gridDim.x*blockIdx.y+blockIdx.x)*blockDim.x + threadIdx.x;
   if(eval_idx < n_path)
   {
      int path_idx = path_idx_mult[eval_idx];

      GT* t = t_array + workspace_size*path_idx + x_t_idx_mult[path_idx];
      one_minor_t += workspace_size*path_idx;

      t_mult[eval_idx] = *t;
      one_minor_t_mult[eval_idx] = *one_minor_t;
      x_mult += eval_idx;
      GT* x = x_array + workspace_size*path_idx + x_t_idx_mult[path_idx]*dim;

      for(int var_idx=0; var_idx<dim; var_idx++)
      {
         x_mult[var_idx*n_path] = x[var_idx];
      }
   }
}

__global__ void eval_mult_end
 ( GT* matrix_mult, GT* matrix, int matrix_dim, 
   int workspace_size, int* path_idx_mult, int n_path )
{
   int eval_idx = (gridDim.x*blockIdx.y+blockIdx.x)*blockDim.x + threadIdx.x;
   if(eval_idx < n_path)
   {
      int path_idx = path_idx_mult[eval_idx];
      matrix_mult += eval_idx;
      matrix += workspace_size*path_idx;

      for(int var_idx=0; var_idx<matrix_dim; var_idx++)
      {
         matrix[var_idx] = matrix_mult[var_idx*n_path];
      }
   }
}

__global__ void eval_mult_end
 ( GT* matrix_mult, GT* matrix, int matrix_dim, int workspace_size,
   int* path_idx_mult, int n_path, int trunk )
{
   int eval_idx = (gridDim.x*blockIdx.y+blockIdx.x)*blockDim.x + threadIdx.x;
   if(eval_idx < n_path)
   {
      int trunk_id = blockIdx.z;
      int path_idx = path_idx_mult[eval_idx];
      matrix_mult += eval_idx;
      matrix += workspace_size*path_idx;

      int start = trunk_id * trunk;
      int end = start + trunk;
      if(end>matrix_dim)
      {
         end = matrix_dim;
      }
      for(int var_idx=start; var_idx<end; var_idx++)
      {
         matrix[var_idx] = matrix_mult[var_idx*n_path];
      }
   }
}

__global__ void eval_mult_end2
 ( GT* matrix_mult, GT* matrix, int matrix_dim, int workspace_size,
   int* path_idx_mult, int n_path, int var_trunk )
{
   __shared__ GT tile[16][17];

   int t_idx = threadIdx.x;
   int sys_trunk= blockDim.x;
   int var_trunk_idx = blockIdx.z;
   int sys_trunk_idx = gridDim.x*blockIdx.y+blockIdx.x;

   int n_path_remain = n_path - sys_trunk_idx*sys_trunk;

   if(t_idx < n_path_remain)
   {
      int start = var_trunk_idx * var_trunk;
      int end = start + var_trunk;
      if(end > matrix_dim)
      {
         end = matrix_dim;
      }
      matrix_mult += sys_trunk_idx*sys_trunk+t_idx;
      for(int var_idx=start; var_idx<end; var_idx++)
      {
         tile[var_idx-start][t_idx] = matrix_mult[var_idx*n_path];
      }
   }
   __syncthreads();

   int n_var_remain = matrix_dim - var_trunk_idx*var_trunk;

   if(t_idx < n_var_remain)
   {
      int start = sys_trunk_idx * sys_trunk;
      int end = start + sys_trunk;
      if(end > n_path)
      {
         end = n_path;
      }
      matrix += var_trunk_idx*var_trunk+t_idx;
      for(int sys_idx=start; sys_idx<end; sys_idx++)
      {
         int path_idx = path_idx_mult[sys_idx];
         matrix[path_idx*workspace_size] = tile[t_idx][sys_idx-start];
      }
   }
}

void eval_mult ( GPUWorkspace& workspace, const GPUInst& inst )
{
   // std::cout << "workspace.n_path_continuous = " 
   // << workspace.n_path_continuous << std::endl;
   int n_path = workspace.n_path_continuous;

   workspace.mon_mult = workspace.coef_mult + workspace.n_coef*n_path;
   workspace.sum_mult = workspace.mon_mult - workspace.n_constant*n_path;

   dim3 init_grid = get_grid(n_path, inst.coef_BS, 1);

   eval_mult_init_zero<<<init_grid, inst.coef_BS>>>
      (workspace.x_array, workspace.t_array, workspace.one_minor_t,
       workspace.x_mult, workspace.newton_t_mult, workspace.one_minor_t_mult,
       workspace.workspace_size, workspace.path_idx, workspace.x_t_idx_mult,
       n_path, inst.dim, workspace.matrix_mult, inst.n_sum_zero,
       inst.sum_zeros);

   dim3 coef_grid = get_grid(n_path, inst.coef_BS, inst.n_coef);
   eval_coef_mult_kernel<<<coef_grid, inst.coef_BS>>>
      (workspace.coef_mult, inst.coef, n_path, workspace.newton_t_mult,
       workspace.one_minor_t_mult);

   /*
     eval_coef_mult_kernel<<<coef_grid, inst.coef_BS>>>
        (workspace.coef_mult, inst.coef, n_path, workspace.t_array,
         workspace.one_minor_t, workspace.workspace_size,
         workspace.x_t_idx_mult);
    */

   /*
     std::cout << "Coef Part" << std::endl;
     CT* gpu_coef_mult = workspace.get_workspace_mult();
     for(int path_idx=0; path_idx<n_path; path_idx++)
     {
        // gpu_workspace_all[path_idx] = workspace.get_workspace(path_idx);
        // gpu_matrix_all[path_idx] = workspace.get_matrix(path_idx);
        std::cout << "path_idx = " << path_idx << std::endl;
        for(int coef_idx=0; coef_idx<10; coef_idx++)
        {
           std::cout << path_idx << " " << coef_idx << " "
                     << gpu_coef_mult[path_idx + n_path*coef_idx];
        }
     }
    */

   dim3 mon_level_grid0
      = get_grid(n_path, inst.mon_global_BS, inst.n_mon_level[0]);

   eval_mon_single_mult_kernel<<<mon_level_grid0, inst.mon_global_BS>>>
      (workspace.mon_mult, workspace.x_mult, workspace.coef_mult,
       inst.mon_pos_start, inst.mon_pos, n_path, workspace.workspace_size);

   /*
     std::cout << "Mon Part" << std::endl;
     CT* gpu_mon_mult = workspace.get_mon_mult();
     for(int path_idx=0; path_idx<n_path; path_idx++)
     {
        // gpu_workspace_all[path_idx] = workspace.get_workspace(path_idx);
        // gpu_matrix_all[path_idx] = workspace.get_matrix(path_idx);
        std::cout << "path_idx = " << path_idx << std::endl;
        for(int coef_idx=0; coef_idx<20; coef_idx++)
        {
           std::cout << path_idx << " " << coef_idx << " "
                     << gpu_mon_mult[path_idx + n_path*coef_idx];
        }
     }
    */
    // std::cout << "n_mon_global = " << inst.mon_global_BS << std::endl;
    dim3 mon_level_grid
       = get_grid(n_path, inst.mon_global_BS, inst.n_mon_global);
    eval_mon_seq_mult_kernel<<<mon_level_grid, inst.mon_global_BS>>>
       (workspace.mon_mult, workspace.x_mult,
        workspace.coef_mult+inst.n_mon_level[0]*n_path,
        inst.mon_pos_start+inst.n_mon_level[0], inst.mon_pos, 
        n_path, workspace.workspace_size);

    // std::cout << "n_path = " << n_path << std::endl;
    std::cout << "inst.n_sum = " << inst.n_sum << std::endl;
    dim3 sum_grid = get_grid(n_path,inst.sum_BS,inst.n_sum);
    eval_sum_seq_mult_kernel<<<sum_grid, inst.sum_BS>>>
       (workspace.matrix_mult, workspace.sum_mult, inst.sum_pos,
        inst.sum_pos_start, n_path);

    /*
      int sum_trunk = 16;
      dim3 sum_grid2 = get_grid(n_path,16,(inst.n_sum-1)/sum_trunk+1);
      eval_sum_seq_mult_transpose_kernel<<<sum_grid2, 16>>>
         (workspace.matrix_mult, workspace.sum_mult, inst.sum_pos,
          inst.sum_pos_start, n_path, inst.n_sum, sum_trunk);
     */

    /*
      std::cout << "Matrix Part" << std::endl;
      CT* gpu_matrix_mult = workspace.get_matrix_mult();
      for(int path_idx=0; path_idx<n_path; path_idx++)
      {
         // gpu_workspace_all[path_idx] = workspace.get_workspace(path_idx);
         // gpu_matrix_all[path_idx] = workspace.get_matrix(path_idx);
         std::cout << "path_idx = " << path_idx << std::endl;
         for(int coef_idx=workspace.dim*workspace.n_eq;
                 coef_idx<(workspace.dim+1)*workspace.n_eq; coef_idx++)
         {
            std::cout << path_idx << " " << coef_idx << " "
                      << gpu_matrix_mult[path_idx + n_path*coef_idx];
         }
      }
     */

    /*
      int matrix_dim = (inst.n_eq+1)*inst.dim;

      dim3 end_grid0 = get_grid(n_path, 32, 1);
      eval_mult_end<<<end_grid0, 32>>>
         (workspace.matrix_mult, workspace.matrix, workspace.n_matrix,
          workspace.workspace_size, workspace.path_idx, n_path);
     */

    int transpose_trunk = 16;
    dim3 end_grid
       = get_grid(n_path, 16, (workspace.n_matrix-1)/transpose_trunk+1);

    /*
      eval_mult_end<<<end_grid, 16>>>
         (workspace.matrix_mult, workspace.matrix, workspace.n_matrix,
          workspace.workspace_size, workspace.path_idx, n_path,
          transpose_trunk);
     */

    eval_mult_end2<<<end_grid, 16>>>
       (workspace.matrix_mult, workspace.matrix, workspace.n_matrix,
        workspace.workspace_size, workspace.path_idx, n_path, transpose_trunk);
}

void eval_eq ( GPUWorkspace& workspace, const GPUInst& inst )
{
   // std::cout << "workspace.n_path_continuous = "
   // << workspace.n_path_continuous << std::endl;
   int n_path = workspace.n_path_continuous;

   /*
     workspace.mon_mult = workspace.coef_mult + workspace.n_coef*n_path;
     workspace.sum_mult = workspace.mon_mult - workspace.n_constant*n_path;

     dim3 init_grid = get_grid(n_path, inst.coef_BS, 1);

     eval_mult_init_zero<<<init_grid, inst.coef_BS>>>
        (workspace.x_array, workspace.t_array, workspace.one_minor_t,
         workspace.x_mult, workspace.newton_t_mult,
         workspace.one_minor_t_mult, workspace.workspace_size,
         workspace.path_idx, workspace.x_t_idx_mult, n_path, inst.dim,
         workspace.matrix_mult, inst.n_sum_zero, inst.sum_zeros);

     dim3 coef_grid = get_grid(n_path, inst.coef_BS, inst.n_coef);
     eval_coef_mult_kernel<<<coef_grid, inst.coef_BS>>>
        (workspace.coef_mult, inst.coef, n_path, workspace.newton_t_mult,
         workspace.one_minor_t_mult);
    */

   /*
     eval_coef_mult_kernel<<<coef_grid, inst.coef_BS>>>
        (workspace.coef_mult, inst.coef, n_path, workspace.t_array,
         workspace.one_minor_t, workspace.workspace_size,
         workspace.x_t_idx_mult);
    */

   /*
     std::cout << "Coef Part" << std::endl;
     CT* gpu_coef_mult = workspace.get_workspace_mult();
     for(int path_idx=0; path_idx<n_path; path_idx++)
     {
        // gpu_workspace_all[path_idx] = workspace.get_workspace(path_idx);
        // gpu_matrix_all[path_idx] = workspace.get_matrix(path_idx);
        std::cout << "path_idx = " << path_idx << std::endl;
        for(int coef_idx=0; coef_idx<10; coef_idx++)
        {
           std::cout << path_idx << " " << coef_idx << " "
                     << gpu_coef_mult[path_idx + n_path*coef_idx];
        }
     }
    */

   dim3 mon_level_grid0
      = get_grid(n_path, inst.mon_global_BS, inst.n_mon_level[0]);

   eval_mon_single_mult_kernel<<<1, 10>>>
      (workspace.mon_mult, workspace.x_mult, workspace.matrix_mult,
       inst.mon_pos_start, inst.mon_pos);

   std::cout << "Mon Part" << std::endl;
   // CT* gpu_mon_mult = workspace.get_mon_mult();
   GT** gpu_matrix_all = new GT*[n_path];
   for(int path_idx=0; path_idx<n_path; path_idx++)
   {
      std::cout << "path_idx = " << path_idx << std::endl;
      // gpu_workspace_all[path_idx] = workspace.get_workspace(path_idx);
      gpu_matrix_all[path_idx] = workspace.get_matrix(path_idx);
      for(int var_idx=0; var_idx<10; var_idx++)
      {
         std::cout << path_idx << " " << var_idx << " " 
                   << gpu_matrix_all[path_idx][var_idx];
      }
   }

   // std::cout << "n_mon_global = " << inst.mon_global_BS << std::endl;

   /*
     dim3 mon_level_grid
        = get_grid(n_path, inst.mon_global_BS, inst.n_mon_global);
     eval_mon_seq_mult_kernel<<<mon_level_grid, inst.mon_global_BS>>>
        (workspace.mon_mult, workspace.x_mult,
         workspace.coef_mult+inst.n_mon_level[0]*n_path,
         inst.mon_pos_start+inst.n_mon_level[0], inst.mon_pos,
         n_path, workspace.workspace_size);

     // std::cout << "n_path = " << n_path << std::endl;
     std::cout << "inst.n_sum = " << inst.n_sum << std::endl;
     dim3 sum_grid = get_grid(n_path,inst.sum_BS,inst.n_sum);
     eval_sum_seq_mult_kernel<<<sum_grid, inst.sum_BS>>>
        (workspace.matrix_mult, workspace.sum_mult, inst.sum_pos,
         inst.sum_pos_start, n_path);
    */

   /*
     int sum_trunk = 16;
     dim3 sum_grid2 = get_grid(n_path,16,(inst.n_sum-1)/sum_trunk+1);
     eval_sum_seq_mult_transpose_kernel<<<sum_grid2, 16>>>
        (workspace.matrix_mult, workspace.sum_mult, inst.sum_pos,
         inst.sum_pos_start, n_path, inst.n_sum, sum_trunk);
    */

   /*
     std::cout << "Matrix Part" << std::endl;
     CT* gpu_matrix_mult = workspace.get_matrix_mult();
     for(int path_idx=0; path_idx<n_path; path_idx++)
     {
        // gpu_workspace_all[path_idx] = workspace.get_workspace(path_idx);
        // gpu_matrix_all[path_idx] = workspace.get_matrix(path_idx);
        std::cout << "path_idx = " << path_idx << std::endl;
        for(int coef_idx=workspace.dim*workspace.n_eq;
                coef_idx<(workspace.dim+1)*workspace.n_eq; coef_idx++)
        {
           std::cout << path_idx << " " << coef_idx << " "
                     << gpu_matrix_mult[path_idx + n_path*coef_idx];
        }
     }
    */

   /*
     int matrix_dim = (inst.n_eq+1)*inst.dim;

     dim3 end_grid0 = get_grid(n_path, 32, 1);
     eval_mult_end<<<end_grid0, 32>>>
        (workspace.matrix_mult, workspace.matrix, workspace.n_matrix,
         workspace.workspace_size, workspace.path_idx, n_path);
    */

   /*
     int transpose_trunk = 16;
     dim3 end_grid
        = get_grid(n_path, 16, (workspace.n_matrix-1)/transpose_trunk+1);
    */

   /*
     eval_mult_end<<<end_grid, 16>>>
        (workspace.matrix_mult, workspace.matrix, workspace.n_matrix,
         workspace.workspace_size, workspace.path_idx, n_path,
         transpose_trunk);
    */

   /*
     eval_mult_end2<<<end_grid, 16>>>
        (workspace.matrix_mult, workspace.matrix, workspace.n_matrix,
         workspace.workspace_size, workspace.path_idx, n_path,
         transpose_trunk);
    */
}

__global__ void eval_mult_eq
 ( GT* matrix_mult, GT* x_mult, GT* workspace_eq, int* eq_pos_start,
   int* mon_pos_start_eq, GT* coef_eq, GT* t_mult, GT* one_minor_t_mult,
   unsigned short* mon_pos_eq, int n_path, int dim, int n_eq )
{
   __shared__ GT eq[11][16];
   __shared__ int mon_pos[11];

   int t_idx = threadIdx.x;
   int path_idx = (gridDim.x*blockIdx.y+blockIdx.x)*blockDim.x + t_idx;
   int eq_idx = blockIdx.z;

   // init constant
   int coef_pos;
   if(path_idx < n_path)
   {
      coef_pos =  2*eq_pos_start[eq_idx];
      eq[dim][t_idx] = t_mult[path_idx]*coef_eq[coef_pos++]
                     + one_minor_t_mult[path_idx]*coef_eq[coef_pos++];
      // init 0 in Jacobian
      for(int var_idx=0; var_idx<dim; var_idx++)
      {
         eq[var_idx][t_idx] = GT(0.0,0.0);
      }
      x_mult += path_idx;
   }
   int* mon_pos_start = mon_pos_start_eq+eq_pos_start[eq_idx];
   int n_mon = *mon_pos_start++;

   for(int mon_idx=0; mon_idx<n_mon; mon_idx++)
   {
      int tmp_start = *mon_pos_start++;
      unsigned short* mon_pos_tmp = mon_pos_eq+tmp_start;

      int n_var = mon_pos_tmp[0];
      if(t_idx<n_var)
      {
         mon_pos[t_idx+1] = mon_pos_tmp[t_idx+1];
      }
      if(path_idx < n_path)
      {
         GT* deri = workspace_eq + tmp_start*n_path + path_idx;
         GT tmp = x_mult[mon_pos[1]*n_path];

         deri[n_path] = tmp;

         for(int var_idx=2; var_idx<n_var; var_idx++)
         {
            tmp *= x_mult[mon_pos[var_idx]*n_path];
            deri[var_idx*n_path] =tmp;
         }
	 // compute coef
	 tmp = t_mult[path_idx]*coef_eq[coef_pos++]
             + one_minor_t_mult[path_idx]*coef_eq[coef_pos++];

         deri -= n_path;
         for(int var_idx=n_var; var_idx>1; var_idx--) 
         {
            int tmp_var = mon_pos[var_idx];
            eq[tmp_var][t_idx] += deri[var_idx*n_path]*tmp;
            tmp *= x_mult[tmp_var*n_path];
         }
         eq[mon_pos[1]][t_idx] +=  tmp;
         eq[dim][t_idx] += x_mult[mon_pos[1]*n_path]*tmp;
      }
   }
   if(path_idx < n_path)
   {
      for(int var_idx=0; var_idx<=dim; var_idx++)
      {
         int output_idx = var_idx*n_eq+eq_idx;
         matrix_mult[output_idx*n_path+path_idx] = eq[var_idx][t_idx];
      }
   }
}

void eval_mult_eq ( GPUWorkspace& workspace, const GPUInst& inst )
{
   int n_path = workspace.n_path_continuous;
   int eval_mult_eq_BS = 16;

   dim3 init_grid = get_grid(n_path, inst.coef_BS, 1);

   eval_mult_init<<<init_grid, inst.coef_BS>>>
      (workspace.x_array, workspace.t_array, workspace.one_minor_t,
       workspace.x_mult, workspace.newton_t_mult, workspace.one_minor_t_mult,
       workspace.workspace_size, workspace.path_idx, workspace.x_t_idx_mult,
       n_path, inst.dim);

   dim3 eq_grid = get_grid(n_path, eval_mult_eq_BS, inst.n_eq);
   eval_mult_eq<<<eq_grid, eval_mult_eq_BS>>>
      (workspace.matrix_mult, workspace.x_mult, workspace.workspace_eq,
       inst.eq_pos_start, inst.mon_pos_start_eq, inst.coef_eq,
       workspace.t_mult, workspace.one_minor_t_mult,  inst.mon_pos_eq,
       n_path, inst.dim, inst.n_eq);

   int transpose_trunk = 16;
   dim3 end_grid 
      = get_grid(n_path, 16, (workspace.n_matrix-1)/transpose_trunk+1);

   /*eval_mult_end<<<end_grid, 16>>>
        (workspace.matrix_mult, workspace.matrix, workspace.n_matrix,
         workspace.workspace_size, workspace.path_idx, n_path,
         transpose_trunk);*/

   eval_mult_end2<<<end_grid, 16>>>
      (workspace.matrix_mult, workspace.matrix, workspace.n_matrix,
       workspace.workspace_size, workspace.path_idx, n_path, transpose_trunk);
}
