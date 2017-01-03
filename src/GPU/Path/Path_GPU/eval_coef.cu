// Mon block evalutation and differentiation on GPU
// template <unsigned int n_th>

__global__ void eval_coef_kernel
 ( GT* workspace_coef, const GT* coef_orig, int n_coef, 
   GT* t, GT* one_minor_t, int workspace_size )
{
   // __shared__ GT div_diff_sh[shmemsize];

   int bidx = (gridDim.x*blockIdx.y+blockIdx.x)*blockDim.x;
   int tidx = threadIdx.x;
   int idx = bidx + tidx;

   int path_idx = blockIdx.z;
   // workspace_coef += workspace_size*path_idx;

   /*int path_idx = blockIdx.z;
     x_predictor += path_idx*np_predictor*dim;
     t_predictor += path_idx*np_predictor;
     x_new += path_idx*dim;*/

   if(idx < n_coef) 
   {
      // workspace_coef[idx] = coef_orig[idx];
      // workspace_coef[idx] = (*t)*coef_orig[2*idx]
      //                     + (*one_minor_t)*coef_orig[2*idx+1];
      workspace_coef[idx] = (*t)*coef_orig[idx]
                          + (*one_minor_t)*coef_orig[idx+n_coef];
   }
}

__global__ void eval_coef_kernel
 ( GT* workspace_coef, const GT* coef_orig, int n_coef,
   GT* t, GT* one_minor_t, int workspace_size, int* x_t_idx )
{
   //__shared__ GT div_diff_sh[shmemsize];
   int bidx = (gridDim.x*blockIdx.y+blockIdx.x)*blockDim.x;
   int tidx = threadIdx.x;
   int idx = bidx + tidx;

   int path_idx = blockIdx.z;
   workspace_coef += workspace_size*path_idx;
   t += workspace_size*path_idx + x_t_idx[path_idx];
   one_minor_t += workspace_size*path_idx;

   /*int path_idx = blockIdx.z;
     x_predictor += path_idx*np_predictor*dim;
     t_predictor += path_idx*np_predictor;
     x_new += path_idx*dim;*/

   if(idx < n_coef)
   {
      // workspace_coef[idx] = coef_orig[idx];
      // workspace_coef[idx] = (*t)*coef_orig[2*idx]
      //                     + (*one_minor_t)*coef_orig[2*idx+1];
      workspace_coef[idx] = (*t)*coef_orig[idx]
                          + (*one_minor_t)*coef_orig[idx+n_coef];
   }
}

__global__ void eval_coef_kernel
 ( GT* workspace_coef, const GT* coef_orig, int n_coef, GT* t,
   GT* one_minor_t, int workspace_size, int* x_t_idx, int* path_idx_mult)
{
   //__shared__ GT div_diff_sh[shmemsize];
   int bidx = (gridDim.x*blockIdx.y+blockIdx.x)*blockDim.x;
   int tidx = threadIdx.x;
   int idx = bidx + tidx;

   int path_idx = path_idx_mult[blockIdx.z];
   workspace_coef += workspace_size*path_idx;
   t += workspace_size*path_idx + x_t_idx[path_idx];
   one_minor_t += workspace_size*path_idx;

   /*int path_idx = blockIdx.z;
     x_predictor += path_idx*np_predictor*dim;
     t_predictor += path_idx*np_predictor;
     x_new += path_idx*dim;*/

   if(idx < n_coef)
   {
      // workspace_coef[idx] = coef_orig[idx];
      // workspace_coef[idx] = (*t)*coef_orig[2*idx]
      //                     + (*one_minor_t)*coef_orig[2*idx+1];
      workspace_coef[idx] = (*t)*coef_orig[idx]
                          + (*one_minor_t)*coef_orig[idx+n_coef];
   }
}

// template <unsigned int n_th>
__global__ void eval_coef_copy_kernel
 ( GT* workspace_coef, const GT* coef_orig, int n_coef )
{
   //__shared__ GT div_diff_sh[shmemsize];
   int bidx = (gridDim.x*blockIdx.y+blockIdx.x)*blockDim.x;
   int tidx = threadIdx.x;
   int idx = bidx + tidx;

   // int path_idx = blockIdx.z;
   // workspace_coef += workspace_size*path_idx;

   /*int path_idx = blockIdx.z;
     x_predictor += path_idx*np_predictor*dim;
     t_predictor += path_idx*np_predictor;
     x_new += path_idx*dim;*/

   if(idx < n_coef)
   {
      // workspace_coef[idx] = coef_orig[idx];
      // workspace_coef[idx] = (*t)*coef_orig[2*idx]
      //                     + (*one_minor_t)*coef_orig[2*idx+1];
      workspace_coef[idx] = coef_orig[idx];
   }
}
