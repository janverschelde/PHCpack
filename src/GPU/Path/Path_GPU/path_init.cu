__global__ void path_mult_x_init_kernel
 ( GT* x_mult_horizontal, GT* x_array_mult, int workspace_size,
   int dim, int n_predictor )
{
   int idx = threadIdx.x;
   // int idx = bidx + tidx;
   int path_idx = blockIdx.x;
   // x_array += path_idx*workspace_size;
   x_array_mult += path_idx*dim*(n_predictor+1);
   x_mult_horizontal += path_idx*dim;

   if(idx < dim)
   {
     // x_array[idx] = x_mult_horizontal[idx];
     x_array_mult[idx] = x_mult_horizontal[idx];
   }
}

__global__ void path_mult_x_finish_kernel
 ( GT* x_mult_horizontal, GT* x_array_mult, int* x_t_idx_mult,
   int workspace_size, int dim, int n_predictor )
{
   int idx = threadIdx.x;
   // int idx = bidx + tidx;
   int path_idx = blockIdx.x;
   int last_x_t_idx = x_t_idx_mult[path_idx]-1;
   if(last_x_t_idx<0)
   {
      last_x_t_idx = n_predictor;
   }
   // x_array += path_idx*workspace_size + last_x_t_idx*dim;
   x_array_mult += path_idx*dim*(n_predictor+1) + last_x_t_idx*dim;
   x_mult_horizontal += path_idx*dim;

   if(idx < dim)
   {
      x_mult_horizontal[idx] = x_array_mult[idx];
   }
}

__global__ void path_mult_init_kernel
 ( GT* t_mult, GT* t_last_mult, GT* delta_t_mult, GT* t_array_mult,
   GT* one_minor_t_mult, GT* alpha, int* path_success, int* newton_success,
   int* end_range, int* n_success, double* max_f_val_last, int* n_point_mult,
   int* x_t_idx_mult, int workspace_size, int* path_idx_mult, int n_path,
   GT* matrix_mult, int n_sum_zero, int* sum_zeros, int n_predictor )
{
   __shared__ int zero_pos[512];

   int BS = blockDim.x;
   int t_idx = threadIdx.x;
   int path_idx = t_idx+blockIdx.x*BS;

   if(path_idx<n_path)
   {
      n_point_mult[path_idx] = 1;
      x_t_idx_mult[path_idx] = 1;

      delta_t_mult[path_idx] = GT(MAX_DELTA_T,0.0);
      t_last_mult[path_idx] = GT(0.0,0.0);
      t_mult[path_idx] = delta_t_mult[path_idx];
      // GT* t = t_array + workspace_size*path_idx;
      // *t =  GT(0.0,0.0);
      // *(t+1) = t_mult[path_idx];
      GT* t = t_array_mult + (n_predictor+1)*path_idx;
      *t =  GT(0.0,0.0);
      *(t+1) = t_mult[path_idx];
      one_minor_t_mult[path_idx] = (*alpha)*(GT(1.0,0) - t_mult[path_idx]);

      end_range[path_idx] = 0;
      n_success[path_idx] = 0;
      path_success[path_idx] = 0;
      path_idx_mult[path_idx] = path_idx;
      max_f_val_last[path_idx] = 1E10;
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

   if(path_idx < n_path)
   {
      matrix_mult += path_idx;
      for(int zero_idx=0; zero_idx<n_sum_zero; zero_idx++)
      {
         int m_idx = zero_pos[zero_idx];
         matrix_mult[m_idx] = GT(0.0,0.0);
      }
   }
}
