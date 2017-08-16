#ifndef NEWTON_CU_
#define NEWTON_CU_

#include "eval.cu"
#include "mgs.cu"

__global__ void update_x_kernel ( GT* x, GT* sol, int dim )
{
   int BS = blockDim.x;
   int bidx = blockIdx.x*BS;
   int tidx = threadIdx.x;
   int idx = bidx + tidx;

   /* int path_idx = blockIdx.z;
      x_predictor += path_idx*np_predictor*dim;
      t_predictor += path_idx*np_predictor;
      x_new += path_idx*dim; */

   if(idx < dim)
   {
      x[idx] = x[idx] - sol[idx];
   }
}

__global__ void update_x_kernel
 ( GT* x, GT* sol, int dim, int workspace_size, int* x_t_idx )
{
   // int BS = blockDim.x;
   // int bidx = blockIdx.x*BS;
   int idx = threadIdx.x;
   // int idx = bidx + tidx;
   int path_idx = blockIdx.x;
   sol += path_idx*dim;
   x += path_idx*workspace_size + x_t_idx[path_idx]*dim;

   /* int path_idx = blockIdx.z;
      x_predictor += path_idx*np_predictor*dim;
      t_predictor += path_idx*np_predictor;
      x_new += path_idx*dim; */

   if(idx < dim)
   {
      x[idx] = x[idx] - sol[idx];
   }
}

__global__ void update_x_kernel_mult ( GT* x, GT* sol, int dim, int n_path )
{
   // int BS = blockDim.x;
   // int bidx = blockIdx.x*BS;
   // int idx = threadIdx.x;
   // int idx = bidx + tidx;
   // int path_idx = blockIdx.x;
   int path_idx = (gridDim.x*blockIdx.y+blockIdx.x)*blockDim.x + threadIdx.x;
   int var_idx = blockIdx.z;

   if(path_idx<n_path)
   {
      sol += path_idx;
      x += path_idx;
      x[var_idx*n_path] = x[var_idx*n_path] - sol[var_idx*n_path];
   }
}

__global__ void update_x_kernel
 ( GT* x, GT* sol, int dim, int workspace_size, int* x_t_idx,
   int* path_idx_mult )
{
   // int BS = blockDim.x;
   // int bidx = blockIdx.x*BS;
   int idx = threadIdx.x;
   // int idx = bidx + tidx;
   int path_idx = path_idx_mult[blockIdx.x];
   sol += blockIdx.x*dim;
   x += path_idx*workspace_size + x_t_idx[path_idx]*dim;

   /* int path_idx = blockIdx.z;
      x_predictor += path_idx*np_predictor*dim;
      t_predictor += path_idx*np_predictor;
      x_new += path_idx*dim; */

   if(idx < dim)
   {
      x[idx] = x[idx] - sol[idx];
   }
}

__global__ void update_x_kernel
 ( GT* x_mult, GT* sol, int dim, int workspace_size,
   int* x_t_idx, int* path_idx_mult, int n_predictor )
{
   // int BS = blockDim.x;
   // int bidx = blockIdx.x*BS;
   int idx = threadIdx.x;
   // int idx = bidx + tidx;
   int path_idx = path_idx_mult[blockIdx.x];
   sol += blockIdx.x*dim;
   // x += path_idx*workspace_size + x_t_idx[path_idx]*dim;
   x_mult += path_idx*dim*(n_predictor+1) + x_t_idx[path_idx]*dim;

   /* int path_idx = blockIdx.z;
      x_predictor += path_idx*np_predictor*dim;
      t_predictor += path_idx*np_predictor;
      x_new += path_idx*dim; */

   if(idx < dim) 
   {
      // x[idx] = x[idx] - sol[idx];
      x_mult[idx] = x_mult[idx] - sol[idx];
   }
}

__global__ void mult_x_init
 ( GT* x_array, GT* t_array, GT* alpha, GT* x_mult, GT* t_mult,
   GT* one_minor_t_mult, int* path_idx_mult, int* x_t_idx_mult,
   int n_path, int dim, int n_predictor )
{
   int t_idx = threadIdx.x;
   int BS = blockDim.x;
   int eval_idx = (gridDim.x*blockIdx.y+blockIdx.x)*BS+t_idx;
   if(eval_idx < n_path)
   {
      int path_idx = path_idx_mult[eval_idx];

      // GT* t = t_array + workspace_size*path_idx + x_t_idx_mult[path_idx];
      GT* tmp_t = t_array + path_idx*(n_predictor+1) + x_t_idx_mult[path_idx];
      // one_minor_t += workspace_size*path_idx;

      t_mult[eval_idx] = *tmp_t;
      one_minor_t_mult[eval_idx] = (*alpha)*(GT(1.0,0) - t_mult[eval_idx]);
      x_mult += eval_idx;
      // GT* x = x_array + workspace_size*path_idx + x_t_idx_mult[path_idx]*dim;
      GT* x_tmp = x_array + path_idx*dim*(n_predictor+1)
                + x_t_idx_mult[path_idx]*dim;

      for(int var_idx=0; var_idx<dim; var_idx++)
      {
         x_mult[var_idx*n_path] = x_tmp[var_idx];
      }
   }
}

__global__ void array_max_double_kernel
 ( GT* sol, int dim, int dimLog2, double* max_array, double* r_max_array,
   double* max_x, int workspace_size )
{
   __shared__ double x_norm[max_array_size];

   int j = threadIdx.x;
   // max for the norm
   int path_idx = blockIdx.x;
   sol += path_idx*workspace_size;
   // sol += path_idx + path_idx*workspace_size;

   x_norm[j] = sol[j].norm1_double();

   dimLog2 -= 1;
   int half_size = 1 << (dimLog2); // sum for the norm

   if(half_size > 16)
   {
      __syncthreads();
   }
   if(j + half_size < dim) 
   {
      if(x_norm[j] < x_norm[j+half_size]) 
      {
         x_norm[j] = x_norm[j+half_size];
      }
   }
   for(int k=0; k < dimLog2; k++) 
   {
      if(half_size > 16) 
      {
         __syncthreads();
      }
      half_size /= 2;
      if(j < half_size) 
      {
         if(x_norm[j] < x_norm[j+half_size]) 
         {
            x_norm[j] = x_norm[j+half_size];
         }
      }
   }
   if(j == 0) 
   {
      max_array[path_idx] = x_norm[0];
      r_max_array[path_idx] = x_norm[0]/max_x[path_idx];
   }
}

__global__ void max_relative_double_kernel
 ( GT* sol, int dim, int n_path_continuous,
   double* max_array, double* r_max_array, double* max_x,
   int workspace_size, int* path_idx_mult )
{
   // __shared__ double x_norm[max_array_size];
   int idx = blockIdx.x*blockDim.x + threadIdx.x;
   if(idx < n_path_continuous)
   {
      int path_idx = path_idx_mult[idx];
      sol += idx*dim;

      double x_norm = sol[0].norm1_double();

      for(int var_idx=1; var_idx<dim; var_idx++)
      {
         double tmp_x_norm = sol[var_idx].norm1_double();
         if(x_norm < tmp_x_norm)
         {
            x_norm = tmp_x_norm;
         }
      }
      max_array[path_idx] = x_norm;
      r_max_array[path_idx] = x_norm/max_x[path_idx];
   }
}

__global__ void max_relative_double_kernel2
 ( GT* sol, int dim, int n_path_continuous, double* max_array,
   double* max_array_last, double* r_max_array, double* max_x,
   int workspace_size, int* path_idx_mult )
{
   //__shared__ double x_norm[max_array_size];
   int idx = blockIdx.x*blockDim.x + threadIdx.x;
   if(idx < n_path_continuous)
   {
      int path_idx = path_idx_mult[idx];
      sol += path_idx*workspace_size;

      double x_norm = sol[0].norm1_double();

      for(int var_idx=1; var_idx<dim; var_idx++)
      {
         double tmp_norm_x = sol[var_idx].norm1_double();
         if(x_norm < tmp_norm_x)
         {
            x_norm = tmp_norm_x;
         }
      }
      max_array[path_idx] = x_norm;
      max_array_last[path_idx] = x_norm;
      r_max_array[path_idx] = x_norm/max_x[path_idx];
   }
}

__global__ void max_relative_double_kernel3
 ( GT* sol, int dim, int n_path_continuous, int n_path, int* path_idx_mult,
   double* max_array, double* max_array_last,
   double* r_max_array, double* max_x )
{
   // __shared__ double x_norm[max_array_size];
   int idx = blockIdx.x*blockDim.x + threadIdx.x;
   if(idx < n_path_continuous)
   {
      int path_idx = path_idx_mult[idx];
      sol += idx;

      double x_norm = sol[0].norm1_double();

      for(int var_idx=1; var_idx<dim; var_idx++)
      {
         double tmp_norm_x = sol[var_idx*n_path].norm1_double();
         if(x_norm < tmp_norm_x) 
         {
            x_norm = tmp_norm_x;
         }
      }
      max_array[path_idx] = x_norm;
      max_array_last[path_idx] = x_norm;
      r_max_array[path_idx] = x_norm/max_x[path_idx];
   }
}

// Not good for a lot path
__global__ void max_relative_double_kernel_tree
 ( GT* sol, int dim, int dimLog2, double* max_array, double* r_max_array,
   double* max_x, int workspace_size, int* path_idx_mult )
{
   __shared__ double x_norm[max_array_size];
   int j = threadIdx.x;
   // max for the norm
   int path_idx = path_idx_mult[blockIdx.x];
   sol += path_idx*workspace_size;
   // sol += path_idx + path_idx*workspace_size;

   x_norm[j] = sol[j].norm1_double();

   dimLog2 -= 1;
   int half_size = 1 << (dimLog2); // sum for the norm

   if(half_size > 16)
   {
      __syncthreads();
   }
   if(j + half_size < dim)
   {
      if(x_norm[j] < x_norm[j+half_size])
      {
         x_norm[j] = x_norm[j+half_size];
      }
   }
   for(int k=0; k < dimLog2; k++) 
   {
      if(half_size > 16) 
      {
         __syncthreads();
      }
      half_size /= 2;
      if(j < half_size) 
      {
         if(x_norm[j] < x_norm[j+half_size]) 
         {
            x_norm[j] = x_norm[j+half_size];
         }
      }
   }
   if(j == 0)
   {
      max_array[path_idx] = x_norm[0];
      r_max_array[path_idx] = x_norm[0]/max_x[path_idx];
   }
}

__global__ void array_max_double_kernel
 ( GT* sol, int dim, int dimLog2, double* max_array,
   int workspace_size, int* x_t_idx )
{
   __shared__ double x_norm[max_array_size];

   int j = threadIdx.x;
   // max for the norm
   int path_idx = blockIdx.x;
   sol += path_idx*workspace_size + dim*x_t_idx[path_idx];

   x_norm[j] = sol[j].norm1_double();

   dimLog2 -= 1;
   int half_size = 1 << (dimLog2);// sum for the norm

   if(half_size > 16) 
   {
      __syncthreads();
   }
   if(j + half_size < dim) 
   {
      if(x_norm[j] < x_norm[j+half_size]) 
      {
         x_norm[j] = x_norm[j+half_size];
      }
   }
   for(int k=0; k < dimLog2; k++) 
   {
      if(half_size > 16) 
      {
         __syncthreads();
      }
      half_size /= 2;
      if(j < half_size)
      {
         if(x_norm[j] < x_norm[j+half_size]) 
         {
            x_norm[j] = x_norm[j+half_size];
	 }
      }
   }
   if(j == 0)
   {
      max_array[path_idx] = x_norm[0];
   }
}

__global__ void max_x_double_kernel_tree
 ( GT* sol, int dim, int dimLog2, double* max_array, 
   int workspace_size, int* x_t_idx, int* path_idx_mult )
{
   __shared__ double x_norm[max_array_size];

   int j = threadIdx.x;
   // max for the norm
   int path_idx = path_idx_mult[blockIdx.x];
   sol += path_idx*workspace_size + dim*x_t_idx[path_idx];

   x_norm[j] = sol[j].norm1_double();

   dimLog2 -= 1;
   int half_size = 1 << (dimLog2); // sum for the norm

   if(half_size > 16) 
   {
      __syncthreads();
   }
   if(j + half_size < dim) 
   {
      if(x_norm[j] < x_norm[j+half_size]) 
      {
         x_norm[j] = x_norm[j+half_size];
      }
   }
   for(int k=0; k < dimLog2; k++) 
   {
      if(half_size > 16)
      {
         __syncthreads();
      }
      half_size /= 2;
      if(j < half_size)
      {
         if(x_norm[j] < x_norm[j+half_size]) 
         {
            x_norm[j] = x_norm[j+half_size];
         }
      }
   }
   if(j == 0) 
   {
      max_array[path_idx] = x_norm[0];
   }
}

__global__ void max_x_double_kernel
 ( GT* sol, int dim, int n_path_continuous, double* max_array,
   int workspace_size, int* x_t_idx, int* path_idx_mult )
{
   // __shared__ double x_norm[max_array_size];
   int idx = blockIdx.x*blockDim.x + threadIdx.x;
   if(idx < n_path_continuous)
   {
      int path_idx = path_idx_mult[idx];
      sol += path_idx*workspace_size + dim*x_t_idx[path_idx];

      double x_norm = sol[0].norm1_double();

      for(int var_idx=1; var_idx<dim; var_idx++)
      {
         double tmp_norm_x = sol[var_idx].norm1_double();
         if(x_norm < tmp_norm_x) 
         {
            x_norm = tmp_norm_x;
         }
      }
      max_array[path_idx] = x_norm;
   }
}

__global__ void max_x_double_kernel_align
 ( GT* sol, int dim, int n_path_continuous, double* max_array,
   int workspace_size, int* x_t_idx, int* path_idx_mult )
{
   // __shared__ double x_norm[max_array_size];
   int idx = blockIdx.x*blockDim.x + threadIdx.x;
   if(idx < n_path_continuous)
   {
      sol += idx;

      double x_norm = sol[0].norm1_double();

      for(int var_idx=1; var_idx<dim; var_idx++)
      {
         double tmp_norm_x = sol[var_idx*n_path_continuous].norm1_double();
         if(x_norm < tmp_norm_x) 
         {
            x_norm = tmp_norm_x;
         }
      }
      int path_idx = path_idx_mult[idx];
      max_array[path_idx] = x_norm;
   }
}

__global__ void zip_kernel
 ( GT* sol, int dim, int n_path_continuous, double* max_array,
   int workspace_size, int* x_t_idx, int* path_idx_mult,
   GT* newton_sol_mult, int* newton_sol_mult_idx,
   GT* one_minor_t, GT* one_minor_t_mult )
{
   // __shared__ double x_norm[max_array_size];
   int idx = blockIdx.x*blockDim.x + threadIdx.x;
   if(idx < n_path_continuous)
   {
      int path_idx = path_idx_mult[idx];
      sol += path_idx*workspace_size + dim*x_t_idx[path_idx];
      newton_sol_mult += idx;

      newton_sol_mult[0] = sol[0];

      for(int var_idx=1; var_idx<dim; var_idx++)
      {
         newton_sol_mult[var_idx*n_path_continuous] = sol[var_idx];
      }
      newton_sol_mult_idx[idx] = path_idx;
      path_idx_mult[idx] = idx;
      one_minor_t_mult[idx] = one_minor_t[path_idx*workspace_size];
   }
}

__global__ void max_x_double_zip_kernel
 ( GT* sol, int dim, int n_path_continuous, double* max_array ) 
{
   // __shared__ double x_norm[max_array_size];
   int idx = blockIdx.x*blockDim.x + threadIdx.x;
   if(idx < n_path_continuous)
   {
      sol += idx;
      double x_norm = sol[0].norm1_double();

      for(int var_idx=1; var_idx<dim; var_idx++)
      {
         double tmp_norm_x = sol[var_idx*n_path_continuous].norm1_double();
         if(x_norm < tmp_norm_x) 
         {
            x_norm = tmp_norm_x;
         }
      }
      max_array[idx] = x_norm;
   }
}

__global__ void max_relative_double_zip_kernel
 ( GT* sol, int dim, int n_path_continuous, double* max_array,
   double* r_max_array, double* max_x )
{
   // __shared__ double x_norm[max_array_size];
   int idx = blockIdx.x*blockDim.x + threadIdx.x;
   if(idx < n_path_continuous)
   {
      sol += idx;
      double x_norm = sol[0].norm1_double();

      for(int var_idx=1; var_idx<dim; var_idx++)
      {
         double tmp_norm_x = sol[var_idx*n_path_continuous].norm1_double();
         if(x_norm < tmp_norm_x) 
         {
            x_norm = tmp_norm_x;
         }
      }
      max_array[idx] = x_norm;
      r_max_array[idx] = x_norm/max_x[idx];
   }
}

__global__ void check_kernel
 ( double* max_f_val_gpu, double* r_max_f_val_gpu, double* max_delta_x_gpu,
   double* r_max_delta_x_gpu, int* path_success, int* success, 
   int* n_point_mult, int* x_t_idx_mult, int n_array, int workspace_size,
   int n_path, int* end_range,
   double err_min_round_off, double err_min_round_off_end )
{
   int path_idx = threadIdx.x+blockIdx.x*blockDim.x;
   if(path_idx<n_path)
   {
      double err_roundoff;
      if(end_range[path_idx]==1)
      {
         err_roundoff = err_min_round_off_end;
      }
      else
      {
         err_roundoff = err_min_round_off;
      }
      if(path_success[path_idx] == 0
         && ( max_f_val_gpu[path_idx] < err_roundoff
         || r_max_f_val_gpu[path_idx] < err_roundoff
         || max_delta_x_gpu[path_idx] < err_roundoff
         || r_max_delta_x_gpu[path_idx] < err_roundoff) )
      {
         success[path_idx] = 1;
         n_point_mult[path_idx]++;
         // remove %
         x_t_idx_mult[path_idx] = (x_t_idx_mult[path_idx]+1)%n_array;
      }
   }
}

__global__ void check_kernel
 ( double* max_delta_x_gpu, double* r_max_delta_x_gpu, int* success,
   int n_path, int* end_range, double err_min_round_off,
   double err_min_round_off_end )
{
   int path_idx = threadIdx.x+blockIdx.x*blockDim.x;
   if(path_idx<n_path)
   {
      double err_roundoff;
      if(end_range[path_idx]==1)
      {
         err_roundoff = err_min_round_off_end;
      }
      else
      {
         err_roundoff = err_min_round_off;
      }
      if(max_delta_x_gpu[path_idx] < err_roundoff
         || r_max_delta_x_gpu[path_idx] < err_roundoff)
      {
         success[path_idx] = 1;
      }
   }
}

__global__ void check_kernel
 ( double* max_f_val_gpu, double* r_max_f_val_gpu, double* max_f_val_last_gpu,
   int* success, int n_path, int* end_range, double err_min_round_off,
   double err_min_round_off_end )
{
   int path_idx = threadIdx.x+blockIdx.x*blockDim.x;
   if(path_idx<n_path)
   {
      if(max_f_val_gpu[path_idx]>max_f_val_last_gpu[path_idx])
      {
         success[path_idx] = -1;
      }
      double err_roundoff;
      if(end_range[path_idx]==1)
      {
         err_roundoff = err_min_round_off_end;
      }
      else
      {
         err_roundoff = err_min_round_off;
      }
      if(max_f_val_gpu[path_idx] < err_roundoff
         || r_max_f_val_gpu[path_idx] < err_roundoff)
      {
         success[path_idx] = 1;
      }
      else
      {
         max_f_val_last_gpu[path_idx] = max_f_val_gpu[path_idx];
      }
   }
}

__global__ void newton_init ( int* success, int n_path )
{
   int newton_path_idx = threadIdx.x+blockIdx.x*blockDim.x;
   if(newton_path_idx<n_path)
   {
      success[newton_path_idx] = 0;
   }
}

__global__ void check_kernel_mult
 ( double* max_delta_x_gpu, double* r_max_delta_x_gpu, int* success,
   int n_path, int* end_range, int* path_idx_mult,
   double err_min_round_off, double err_min_round_off_end )
{
   int newton_path_idx = threadIdx.x+blockIdx.x*blockDim.x;
   int path_idx = path_idx_mult[newton_path_idx];
   if(newton_path_idx<n_path)
   {
      double err_roundoff;
      if(end_range[path_idx]==1)
      {
         err_roundoff = err_min_round_off_end;
      }
      else
      {
         err_roundoff = err_min_round_off;
      }
      if(max_delta_x_gpu[path_idx] < err_roundoff
         || r_max_delta_x_gpu[path_idx] < err_roundoff)
      {
         success[path_idx] = 1;
      }
   }
}

__global__ void array_max_double_kernel
 ( GT* sol, int dim, int dimLog2, double* max_delta_x ) 
{
   __shared__ double delta_x[max_array_size];

   int j = threadIdx.x;
   // max for the norm
   delta_x[j] = sol[j].norm1_double();

   dimLog2 -= 1;
   int half_size = 1 << (dimLog2); // sum for the norm

   if(half_size > 16)
   {
      __syncthreads();
   }
   if(j + half_size < dim)
   {
      if(delta_x[j] < delta_x[j+half_size]) 
      {
         delta_x[j] = delta_x[j+half_size];
      }
   }
   for(int k=0; k < dimLog2; k++) 
   {
      if(half_size > 16)
      {
         __syncthreads();
      }
      half_size /= 2;
      if(j < half_size)
      {
         if(delta_x[j] < delta_x[j+half_size])
         {
            delta_x[j] = delta_x[j+half_size];
         }
      }
   }
   if(j == 0)
   {
      *max_delta_x = delta_x[0];
   }
}

bool newton_single
 ( GPUWorkspace& workspace, GPUInst& inst, Parameter path_parameter,
   bool end_range=false )
{
   bool debug = false; // debug = true;
   bool success = false;
   int rowsLog2 = log2ceil(inst.n_eq); // ceil for sum reduction
   int dimLog2 = log2ceil(inst.dim); // ceil for sum reduction

   double max_x;
   double max_f_val;
   double r_max_f_val;
   double max_delta_x;
   double r_max_delta_x;

   double* max_x_gpu;
   cudaMalloc((void **) &max_x_gpu, sizeof(double));

   double* max_f_val_gpu;
   cudaMalloc((void **) &max_f_val_gpu, sizeof(double));

   double* max_delta_x_gpu;
   cudaMalloc((void **) &max_delta_x_gpu, sizeof(double));

   double err_round_off;
   if(end_range==true)
   {
      err_round_off = path_parameter.err_min_round_off_refine;
   }
   else
   {
      err_round_off = path_parameter.err_min_round_off;
   }
   array_max_double_kernel<<<1, inst.dim>>>
      (workspace.x, inst.dim, dimLog2, max_x_gpu);

   cudaMemcpy(&max_x, max_x_gpu, sizeof(double), cudaMemcpyDeviceToHost);

   if(debug)
   {
      std::cout << "          max_x : " << max_x << std::endl;
   }
   eval(workspace, inst);
   inst.n_eval_GPU++;

   array_max_double_kernel<<<1, inst.n_eq>>>
      (workspace.f_val, inst.n_eq, rowsLog2, max_f_val_gpu);

   cudaMemcpy(&max_f_val, max_f_val_gpu,
              sizeof(double), cudaMemcpyDeviceToHost);

   r_max_f_val = max_f_val/max_x;

   if(debug)
   {
      std::cout << "   residual(a&r): " << max_f_val 
                << " " << r_max_f_val << std::endl;
   }
   if(max_f_val < err_round_off || r_max_f_val < err_round_off)
   {
      success = 1;
      return success;
   }
   double last_max_f_val = max_f_val;

   for(int i = 0; i < path_parameter.max_it; i++)
   {
      if(debug)
      {
         cout << "  Iteration " << i << endl;
      }
      if(inst.n_eq <= BS_QR)
      {
         mgs_small_with_delta
            (workspace.matrix, workspace.R, workspace.sol,
             inst.n_eq, inst.dim + 1, max_delta_x_gpu);
         /* CT* tmp_sol =  workspace.get_sol();
            for(int var_idx=0; var_idx<inst.dim; var_idx++)
            {
               std::cout << var_idx << " " << tmp_sol[var_idx];
            }*/
      }
      else
      {
         mgs_large_block
            (workspace.matrix, workspace.R, workspace.P,
             workspace.sol, inst.n_eq, inst.dim + 1);
         //mgs_large(workspace.V, workspace.R, workspace.sol, inst.n_eq,
         //          inst.dim+1);
         array_max_double_kernel<<<1,inst.dim>>>
            (workspace.sol, inst.dim, dimLog2, max_delta_x_gpu);
      }
      inst.n_mgs_GPU++;

      cudaMemcpy(&max_delta_x, max_delta_x_gpu,
                 sizeof(double), cudaMemcpyDeviceToHost);
      r_max_delta_x = max_delta_x/max_x;

      if(debug)
      {
         std::cout << " correction(a&r): " << max_delta_x
                   << " " << r_max_delta_x << std::endl;
      }
      update_x_kernel<<<inst.dim_grid, inst.dim_BS>>>
         (workspace.x, workspace.sol, inst.dim);

      if(max_delta_x < err_round_off || r_max_delta_x < err_round_off)
      {
         success = 1;
         break;
      }
      array_max_double_kernel<<<1, inst.dim>>>
         (workspace.x, inst.dim, dimLog2, max_x_gpu);

      cudaMemcpy(&max_x, max_x_gpu, sizeof(double), cudaMemcpyDeviceToHost);

      if(debug)
      {
         std::cout << "          max_x : " << max_x << std::endl;
      }
      eval(workspace, inst);
      inst.n_eval_GPU++;

      array_max_double_kernel<<<1, inst.n_eq>>>
         (workspace.f_val, inst.n_eq, rowsLog2, max_f_val_gpu);

      cudaMemcpy(&max_f_val, max_f_val_gpu,
                 sizeof(double), cudaMemcpyDeviceToHost);
      r_max_f_val = max_f_val/max_x;

      if(debug)
      {
         std::cout << "   residual(a&r): " << max_f_val
                   << " " << r_max_f_val << std::endl;
      }
      if (max_f_val > last_max_f_val)
      {
         success = 0;
         break;
      }
      if(max_f_val < err_round_off || r_max_f_val < err_round_off)
      {
         success = 1;
         break;
      }
      last_max_f_val = max_f_val;
   }
   return success;
}

// only difference is the function value check
bool newton_single2
 ( GPUWorkspace& workspace, GPUInst& inst, Parameter path_parameter,
   bool end_range=false )
{
   bool debug = false; // debug = true;
   bool success = false;
   int rowsLog2 = log2ceil(inst.n_eq); // ceil for sum reduction
   int dimLog2 = log2ceil(inst.dim);   // ceil for sum reduction

   double max_x;
   double max_f_val;
   double r_max_f_val;
   double max_delta_x;
   double r_max_delta_x;

   double* max_x_gpu;
   cudaMalloc((void **) &max_x_gpu, sizeof(double));

   double* max_f_val_gpu;
   cudaMalloc((void **) &max_f_val_gpu, sizeof(double));

   double* max_delta_x_gpu;
   cudaMalloc((void **) &max_delta_x_gpu, sizeof(double));

   double err_round_off;
   if(end_range==true)
   {
      err_round_off = path_parameter.err_min_round_off_refine;
   }
   else
   {
      err_round_off = path_parameter.err_min_round_off;
   }
   array_max_double_kernel<<<1, inst.dim>>>
      (workspace.x, inst.dim,dimLog2, max_x_gpu);

   cudaMemcpy(&max_x, max_x_gpu, sizeof(double), cudaMemcpyDeviceToHost);

   if(debug)
   {
      std::cout << "          max_x : " << max_x << std::endl;
   }
   eval(workspace, inst);
   inst.n_eval_GPU++;

   array_max_double_kernel<<<1, inst.n_eq>>>
      (workspace.f_val, inst.n_eq, rowsLog2, max_f_val_gpu);

   cudaMemcpy(&max_f_val, max_f_val_gpu,
              sizeof(double), cudaMemcpyDeviceToHost);

   r_max_f_val = max_f_val/max_x;

   if(debug)
   {
      std::cout << "   residual(a&r): " << max_f_val
                << " " << r_max_f_val << std::endl;
   }
   if(max_f_val < err_round_off || r_max_f_val < err_round_off)
   {
      success = 1;
      return success;
   }
   double last_max_f_val = max_f_val;

   for(int i = 0; i < path_parameter.max_it; i++) 
   {
      if(debug)
      {
         cout << "  Iteration " << i << endl;
      }
      if(inst.n_eq <= BS_QR) 
      {
         mgs_small_with_delta
            (workspace.matrix, workspace.R, workspace.sol,
             inst.n_eq, inst.dim + 1, max_delta_x_gpu);
      }
      else
      {
         mgs_large_block
            (workspace.matrix, workspace.R, workspace.P,
             workspace.sol, inst.n_eq, inst.dim + 1);
         // mgs_large(workspace.V, workspace.R, workspace.sol,
         // inst.n_eq, inst.dim+1);
         array_max_double_kernel<<<1,inst.dim>>>
            (workspace.sol, inst.dim, dimLog2, max_delta_x_gpu);
      }
      inst.n_mgs_GPU++;

      cudaMemcpy(&max_delta_x, max_delta_x_gpu,
                 sizeof(double), cudaMemcpyDeviceToHost);
      r_max_delta_x = max_delta_x/max_x;

      if(debug)
      {
         std::cout << " correction(a&r): " << max_delta_x
                   << " " << r_max_delta_x << std::endl;
      }
      update_x_kernel<<<inst.dim_grid, inst.dim_BS>>>
         (workspace.x, workspace.sol, inst.dim);

      if(max_delta_x < err_round_off || r_max_delta_x < err_round_off)
      {
         success = 1;
         break;
      }
      array_max_double_kernel<<<1, inst.dim>>>(workspace.x, inst.dim, \
				dimLog2, max_x_gpu);

      cudaMemcpy(&max_x, max_x_gpu, sizeof(double), cudaMemcpyDeviceToHost);

      if(debug)
      {
         std::cout << "          max_x : " << max_x << std::endl;
      }
      eval(workspace, inst);
      inst.n_eval_GPU++;

      array_max_double_kernel<<<1, inst.n_eq>>>
         (workspace.f_val, inst.n_eq, rowsLog2, max_f_val_gpu);

      cudaMemcpy(&max_f_val, max_f_val_gpu,
                 sizeof(double), cudaMemcpyDeviceToHost);
      r_max_f_val = max_f_val/max_x;

      if(debug) workspace.print_f_val();

      if(debug)
      {
         std::cout << "   residual(a&r): " << max_f_val
                   << " " << r_max_f_val << std::endl;
      }

      /* if (max_f_val > last_max_f_val)
         {
            success = 0;
            break;
         } */

      if(max_f_val < err_round_off || r_max_f_val < err_round_off)
      {
         success = 1;
         break;
      }
      last_max_f_val = max_f_val;
   }
   return success;
}

bool newton
 ( GPUWorkspace& workspace, GPUInst& inst, Parameter path_parameter,
   bool debug=false )
{
   debug = true;
   int path_idx_test = 0;

   int rowsLog2 = log2ceil(inst.n_eq); // ceil for sum reduction
   int dimLog2 = log2ceil(inst.dim); // ceil for sum reduction
   int n_path = workspace.n_path;

   for(int path_idx=0; path_idx<n_path; path_idx++)
   {
      workspace.newton_success_host[path_idx]
         = workspace.path_success_host[path_idx];
      // std::cout << "workspace.success_host[path_idx] = "
      // << workspace.success_host[path_idx] << std::endl;
   }
   cudaMemcpy(workspace.newton_success, workspace.newton_success_host,
              n_path*sizeof(int), cudaMemcpyHostToDevice);

   int n_path_continuous = 0;
   // std::cout << "newton_success" << std::endl;
   for(int path_idx=0; path_idx<n_path; path_idx++)
   {
      // std::cout << path_idx << " "
      // << workspace.newton_success_host[path_idx] << std::endl;
      if(workspace.newton_success_host[path_idx] == 0)
      {
         workspace.path_idx_host[n_path_continuous] = path_idx;
         n_path_continuous += 1;
      }
   }
   workspace.n_path_continuous = n_path_continuous;

   /*if(debug)
     {
        std::cout << "n_path_continuous" << std::endl;
        for(int path_idx=0; path_idx<n_path_continuous; path_idx++)
        {
           std::cout << path_idx << " "
                     << workspace.path_idx_host[path_idx] << std::endl;
        }
     }*/
   cudaMemcpy(workspace.path_idx, workspace.path_idx_host, n_path*sizeof(int),
              cudaMemcpyHostToDevice);

   dim3 max_grid = get_grid(n_path_continuous,inst.predict_BS,1);
   max_x_double_kernel<<<max_grid, inst.predict_BS>>>
      (workspace.x_array, inst.dim, n_path_continuous, workspace.max_x_gpu,
       workspace.workspace_size, workspace.x_t_idx_mult, workspace.path_idx);

   if(debug)
   {
      cudaMemcpy(workspace.max_x_host, workspace.max_x_gpu,
                 n_path*sizeof(double), cudaMemcpyDeviceToHost);

      for(int path_idx=0; path_idx<n_path; path_idx++)
      {
         // if(path_idx == path_idx_test){
         std::cout << "       max_x_gpu " << path_idx << "= "
                   << workspace.max_x_host[path_idx] << std::endl;
         //}
      }
   }
   eval(workspace, inst);
   inst.n_eval_GPU++;

   max_relative_double_kernel<<<max_grid, inst.predict_BS>>>
      (workspace.f_val, inst.n_eq, n_path_continuous,
       workspace.max_f_val_gpu, workspace.r_max_f_val_gpu,
       workspace.max_x_gpu, workspace.workspace_size, workspace.path_idx);

   if(debug)
   {
      cudaMemcpy(workspace.max_f_val_host, workspace.max_f_val_gpu,
                 n_path*sizeof(double), cudaMemcpyDeviceToHost);
      cudaMemcpy(workspace.r_max_f_val_host, workspace.r_max_f_val_gpu,
                 n_path*sizeof(double), cudaMemcpyDeviceToHost);

      for(int path_idx=0; path_idx<n_path; path_idx++)
      {
         std::cout << "       max_f_value " << path_idx << " = "
                   << workspace.max_f_val_host[path_idx]
                   << " " << workspace.r_max_f_val_host[path_idx]
                   << std::endl;
      }
   }
   dim3 check_grid = get_grid(n_path,inst.predict_BS,1);

   for(int it_idx = 0; it_idx < path_parameter.max_it; it_idx++)
   {
      if(debug)
      {
         std::cout << "  Iteration " << it_idx << std::endl;
      }
      if(inst.dim <= BS_QR)
      {
         mgs_small_dynamic_idx(workspace.V, workspace.R, workspace.sol,
            inst.n_eq, inst.dim+1, workspace.small_mgs_size,
            workspace.n_matrix, workspace.n_matrix_R, n_path_continuous,
            workspace.path_idx);
         // array_max_double_kernel<<<1,inst.dim>>>
         //    (workspace.sol, inst.dim, dimLog2, max_delta_x_gpu);
         /*
           CT** sol_gpu = new CT*[n_path];
           CT** matrix_gpu_q = new CT*[n_path];
           CT** matrix_gpu_r = new CT*[n_path];
           for(int path_idx=0; path_idx<n_path; path_idx++)
           {
              sol_gpu[path_idx] = workspace.get_sol(path_idx);
              for(int var_idx=0; var_idx<inst.dim; var_idx++)
              {
                 std::cout << path_idx << " " << var_idx
                           << " " << sol_gpu[path_idx][var_idx];
              }
              matrix_gpu_q[path_idx] = workspace.get_matrix(path_idx);
              matrix_gpu_r[path_idx] = workspace.get_matrix_r(path_idx);
           }*/
         max_relative_double_kernel<<<max_grid, inst.predict_BS>>>
            (workspace.sol, inst.dim, n_path_continuous,
             workspace.max_delta_x_gpu, workspace.r_max_delta_x_gpu,
             workspace.max_x_gpu, workspace.workspace_size,
             workspace.path_idx);
      }
      else
      {
         mgs_large_block(workspace.V, workspace.R, workspace.P,
                         workspace.sol, inst.n_eq, inst.dim+1);
         // array_max_double_kernel<<<1,inst.dim>>>
         //    (workspace.sol, inst.dim, dimLog2, max_delta_x_gpu);
      }
      if(debug)
      {
         cudaMemcpy(workspace.max_delta_x_host, workspace.max_delta_x_gpu,
                    n_path*sizeof(double), cudaMemcpyDeviceToHost);
         cudaMemcpy(workspace.r_max_delta_x_host,
                    workspace.r_max_delta_x_gpu, n_path*sizeof(double),
                    cudaMemcpyDeviceToHost);

         for(int path_idx=0; path_idx<n_path; path_idx++)
         {
            // if(path_idx == path_idx_test){
            std::cout << "       max_delta_x " << path_idx
                      << " = " << workspace.max_delta_x_host[path_idx]
                      << " " << workspace.r_max_delta_x_host[path_idx]
                      << std::endl;
            //}
         }
      }
      update_x_kernel<<<n_path_continuous, inst.dim>>>
         (workspace.x_array, workspace.sol, inst.dim, 
          workspace.workspace_size, workspace.x_t_idx_mult,
          workspace.path_idx);

      // workspace.print_x_mult();
      check_kernel<<<check_grid, inst.predict_BS>>>
         (workspace.max_delta_x_gpu, workspace.r_max_delta_x_gpu,
          workspace.newton_success, n_path, workspace.end_range,
          path_parameter.err_min_round_off,
          path_parameter.err_min_round_off_refine);

      cudaMemcpy(workspace.newton_success_host, workspace.newton_success,
                 n_path*sizeof(int), cudaMemcpyDeviceToHost);

      n_path_continuous = 0;
      // std::cout << "newton_success" << std::endl;
      for(int path_idx=0; path_idx<n_path; path_idx++)
      {
         // std::cout << path_idx << " "
         // << workspace.newton_success_host[path_idx] << std::endl;
         if(workspace.newton_success_host[path_idx] == 0)
         {
            workspace.path_idx_host[n_path_continuous] = path_idx;
            n_path_continuous += 1;
         }
      }
      workspace.n_path_continuous = n_path_continuous;
      max_grid = get_grid(n_path_continuous,inst.predict_BS,1);
      /*
        std::cout << "n_path_continuous" << std::endl;
        for(int path_idx=0; path_idx<n_path_continuous; path_idx++)
        {
           std::cout << path_idx << " "
                     << workspace.path_idx_host[path_idx] << std::endl;
        }
       */
      cudaMemcpy(workspace.path_idx, workspace.path_idx_host,
                 n_path*sizeof(int), cudaMemcpyHostToDevice);

      // std::cout << "Correct X:" << std::endl;
      // workspace.print_x_mult();

      max_x_double_kernel<<<max_grid, inst.predict_BS>>>
         (workspace.x_array, inst.dim, n_path_continuous,
          workspace.max_x_gpu, workspace.workspace_size,
          workspace.x_t_idx_mult, workspace.path_idx);

      if(debug)
      {
         cudaMemcpy(workspace.max_x_host, workspace.max_x_gpu,
                    n_path*sizeof(double), cudaMemcpyDeviceToHost);
         for(int path_idx=0; path_idx<n_path; path_idx++)
         {
            // if(path_idx == path_idx_test){
            std::cout << "       max_x_gpu " << path_idx << "= "
                      << workspace.max_x_host[path_idx] << std::endl;
            //}
         }
      }
      /*
        cudaMemcpy(max_x_host, max_x_gpu, n_path*sizeof(double),
                   cudaMemcpyDeviceToHost);

        for(int path_idx=0; path_idx<n_path; path_idx++)
        {
           std::cout << "       max_x_gpu " << path_idx << "= "
                     << max_x_host[path_idx] << std::endl;
        }
       */
      eval(workspace, inst);
      // inst.n_eval_GPU++;

      max_relative_double_kernel<<<max_grid, inst.predict_BS>>>
         (workspace.f_val, inst.n_eq, n_path_continuous,
          workspace.max_f_val_gpu, workspace.r_max_f_val_gpu,
          workspace.max_x_gpu, workspace.workspace_size, workspace.path_idx);

      check_kernel<<<check_grid, inst.predict_BS>>>
         (workspace.max_f_val_gpu, workspace.r_max_f_val_gpu,
          workspace.newton_success, n_path, workspace.end_range,
          path_parameter.err_min_round_off,
          path_parameter.err_min_round_off_refine);

      cudaMemcpy(workspace.newton_success_host, workspace.newton_success,
                 n_path*sizeof(int), cudaMemcpyDeviceToHost);

      n_path_continuous = 0;
      // std::cout << "newton_success" << std::endl;
      for(int path_idx=0; path_idx<n_path; path_idx++)
      {
         // std::cout << path_idx << " "
         // << workspace.newton_success_host[path_idx] << std::endl;
         if(workspace.newton_success_host[path_idx] == 0)
         {
            workspace.path_idx_host[n_path_continuous] = path_idx;
            n_path_continuous += 1;
         }
      }
      workspace.n_path_continuous = n_path_continuous;
      max_grid = get_grid(n_path_continuous,inst.predict_BS,1);
      std::cout << "n_path_continuous" << std::endl;
      for(int path_idx=0; path_idx<n_path_continuous; path_idx++)
      {
         std::cout << path_idx << " "
                   << workspace.path_idx_host[path_idx] << std::endl;
      }
      cudaMemcpy(workspace.path_idx, workspace.path_idx_host,
                 n_path*sizeof(int), cudaMemcpyHostToDevice);

      if(debug)
      {
         cudaMemcpy(workspace.max_f_val_host, workspace.max_f_val_gpu,
                    n_path*sizeof(double), cudaMemcpyDeviceToHost);
         cudaMemcpy(workspace.r_max_f_val_host, workspace.r_max_f_val_gpu,
                    n_path*sizeof(double), cudaMemcpyDeviceToHost);

         for(int path_idx=0; path_idx<n_path; path_idx++)
         {
            // if(path_idx == path_idx_test){
            std::cout << "       max_f_value " << path_idx
                      << " = " << workspace.max_f_val_host[path_idx]
                      << " " <<  workspace.r_max_f_val_host[path_idx]
                      << std::endl;
            //}
         }
      }
   }
   dim3 check_grid2 = get_grid(n_path,inst.predict_BS,1);
   check_kernel<<<check_grid, inst.predict_BS>>>
      (workspace.max_f_val_gpu, workspace.r_max_f_val_gpu,
       workspace.max_delta_x_gpu, workspace.r_max_delta_x_gpu,
       workspace.path_success, workspace.newton_success,
       workspace.n_point_mult, workspace.x_t_idx_mult, workspace.n_array,
       workspace.workspace_size, n_path, workspace.end_range,
       path_parameter.err_min_round_off,
       path_parameter.err_min_round_off_refine);

   cudaMemcpy(workspace.newton_success_host, workspace.newton_success,
              n_path*sizeof(int), cudaMemcpyDeviceToHost);

   if(debug)
   {
      std::cout << workspace.newton_success_host[path_idx_test] << std::endl;
   }
   return true;
}

bool newton_align
 ( GPUWorkspace& workspace, GPUInst& inst, Parameter path_parameter,
   bool debug=false )
{
   // std::cout << "Newton Align" << std::endl;
   // debug = true;
   int path_idx_test = 0;
   int debug_all = false;
   // debug_all = true;

   int rowsLog2 = log2ceil(inst.n_eq); // ceil for sum reduction
   int dimLog2 = log2ceil(inst.dim); // ceil for sum reduction
   int n_path = workspace.n_path;

   for(int path_idx=0; path_idx<n_path; path_idx++)
   {
      workspace.newton_success_host[path_idx]
         = workspace.path_success_host[path_idx];
      // std::cout << "workspace.success_host[path_idx] = "
      //           << workspace.success_host[path_idx] << std::endl;
   }
   cudaMemcpy(workspace.newton_success,workspace.newton_success_host,
              n_path*sizeof(int),cudaMemcpyHostToDevice);
   // std::cout << "newton_success" << std::endl;
   int n_path_continuous = 0;
   for(int path_idx=0; path_idx<n_path; path_idx++)
   {
      // std::cout << path_idx << " "
      //           << workspace.newton_success_host[path_idx] << std::endl;
      if(workspace.newton_success_host[path_idx] == 0)
      {
         workspace.path_idx_host[n_path_continuous] = path_idx;
         n_path_continuous += 1;
      }
   }
   cudaMemcpy(workspace.path_idx, workspace.path_idx_host,
              n_path*sizeof(int), cudaMemcpyHostToDevice);

   dim3 init_grid = get_grid(n_path_continuous, inst.coef_BS, 1);

   mult_x_init<<<init_grid, inst.coef_BS>>>
      (workspace.x_array, workspace.t_array, workspace.alpha_gpu,
       workspace.x_mult, workspace.newton_t_mult, workspace.one_minor_t,
       workspace.path_idx, workspace.x_t_idx_mult, n_path_continuous,
       inst.dim, workspace.n_predictor);

   dim3 max_grid = get_grid(n_path_continuous,inst.predict_BS,1);
   max_x_double_kernel_align<<<max_grid, inst.predict_BS>>>
      (workspace.x_mult, inst.dim, n_path_continuous, workspace.max_x_gpu,
       workspace.workspace_size, workspace.x_t_idx_mult, workspace.path_idx);

   if(debug)
   {
      cudaMemcpy(workspace.max_x_host, workspace.max_x_gpu,
                 n_path*sizeof(double), cudaMemcpyDeviceToHost);

      for(int path_idx=0; path_idx<n_path; path_idx++)
      {
         if(path_idx == path_idx_test || debug_all)
         {
            std::cout << "       max_x_gpu " << path_idx << "= "
                      << workspace.max_x_host[path_idx] << std::endl;
         }
      }
   }
   eval_mult(workspace, inst);
   inst.n_eval_GPU++;

   /*
     max_relative_double_kernel2<<<max_grid, inst.predict_BS>>>
        (workspace.f_val, inst.n_eq, n_path_continuous,
         workspace.max_f_val_gpu, workspace.max_f_val_last_gpu,
         workspace.r_max_f_val_gpu, workspace.max_x_gpu,
         workspace.workspace_size, workspace.path_idx);
    */

   max_relative_double_kernel3<<<max_grid, inst.predict_BS>>>
      (workspace.f_val, inst.n_eq, n_path_continuous, workspace.n_path,
       workspace.path_idx, workspace.max_f_val_gpu,
       workspace.max_f_val_last_gpu, workspace.r_max_f_val_gpu,
       workspace.max_x_gpu);

   if(debug)
   {
      cudaMemcpy(workspace.max_f_val_host, workspace.max_f_val_gpu,
                 n_path*sizeof(double), cudaMemcpyDeviceToHost);
      cudaMemcpy(workspace.r_max_f_val_host, workspace.r_max_f_val_gpu,
                 n_path*sizeof(double), cudaMemcpyDeviceToHost);

      for(int path_idx=0; path_idx<n_path; path_idx++)
      {
         if(path_idx == path_idx_test || debug_all)
         {
            std::cout << "       max_f_value " << path_idx << " = "
                      << workspace.max_f_val_host[path_idx] 
                      << " " << workspace.r_max_f_val_host[path_idx]
                      << std::endl;
         }
      }
   }
   dim3 check_grid = get_grid(n_path,inst.predict_BS,1);

   for(int it_idx = 0; it_idx < path_parameter.max_it; it_idx++) 
   {
      if(debug)
      {
         std::cout << "  Iteration " << it_idx << std::endl;
      }
      if(inst.dim <= BS_QR)
      {
         if(sizeof(GT) > 16)
         {
            mgs_small_idx(workspace.matrix_horizontal_mult, workspace.R,
                          workspace.sol, inst.n_eq, inst.dim+1,
                          workspace.n_matrix, workspace.n_matrix_R,
                          n_path_continuous, workspace.path_idx);
         }
         else
         {
            mgs_small_dynamic_idx
               (workspace.matrix_horizontal_mult, workspace.R, workspace.sol,
                inst.n_eq, inst.dim+1, workspace.small_mgs_size,
                workspace.n_matrix, workspace.n_matrix_R, n_path_continuous,
                workspace.path_idx);
         }
         max_relative_double_kernel<<<max_grid, inst.predict_BS>>>
            (workspace.sol, inst.dim, n_path_continuous,
             workspace.max_delta_x_gpu, workspace.r_max_delta_x_gpu,
             workspace.max_x_gpu, workspace.workspace_size,
             workspace.path_idx);
      }
      else
      {
         mgs_large_block(workspace.V, workspace.R, workspace.P,
                         workspace.sol, inst.n_eq, inst.dim+1);
         // array_max_double_kernel<<<1,inst.dim>>>(workspace.sol,
         // inst.dim, dimLog2, max_delta_x_gpu);
      }
      if(debug)
      {
         cudaMemcpy(workspace.max_delta_x_host, workspace.max_delta_x_gpu,
                    n_path*sizeof(double), cudaMemcpyDeviceToHost);
         cudaMemcpy(workspace.r_max_delta_x_host, workspace.r_max_delta_x_gpu,
                    n_path*sizeof(double), cudaMemcpyDeviceToHost);

         for(int path_idx=0; path_idx<n_path; path_idx++)
         {
            if(path_idx == path_idx_test|| debug_all)
            {
               std::cout << "       max_delta_x " << path_idx << " = "
                         << workspace.max_delta_x_host[path_idx]
                         << " " <<  workspace.r_max_delta_x_host[path_idx]
                         << std::endl;
            }
         }
      }
      update_x_kernel<<<n_path_continuous, inst.dim>>>
         (workspace.x_array, workspace.sol, inst.dim,
          workspace.workspace_size, workspace.x_t_idx_mult,
          workspace.path_idx, workspace.n_predictor);

      // workspace.print_x_mult();
      check_kernel<<<check_grid, inst.predict_BS>>>
         (workspace.max_delta_x_gpu, workspace.r_max_delta_x_gpu,
          workspace.newton_success, n_path, workspace.end_range,
          path_parameter.err_min_round_off,
          path_parameter.err_min_round_off_refine);

      cudaMemcpy(workspace.newton_success_host, workspace.newton_success,
                 n_path*sizeof(int), cudaMemcpyDeviceToHost);

      n_path_continuous = 0;
      // std::cout << "newton_success" << std::endl;
      for(int path_idx=0; path_idx<n_path; path_idx++)
      {
         // std::cout << path_idx << " "
         //  << workspace.newton_success_host[path_idx] << std::endl;
         if(workspace.newton_success_host[path_idx] == 0)
         {
            workspace.path_idx_host[n_path_continuous] = path_idx;
            n_path_continuous += 1;
         }
      }
      if(n_path_continuous==0)
      {
         break;
      }
      workspace.n_path_continuous = n_path_continuous;
      max_grid = get_grid(n_path_continuous,inst.predict_BS,1);
      /*
        std::cout << "n_path_continuous = " << n_path_continuous << " : ";
        for(int path_idx=0; path_idx<n_path_continuous; path_idx++)
        {
           std::cout << workspace.path_idx_host[path_idx] << ", ";
        }
        std::cout << std::endl;
       */
      cudaMemcpy(workspace.path_idx, workspace.path_idx_host,
                 n_path*sizeof(int), cudaMemcpyHostToDevice);

      // std::cout << "Correct X:" << std::endl;
      // workspace.print_x_mult();
      init_grid = get_grid(n_path_continuous, inst.coef_BS, 1);

      mult_x_init<<<init_grid, inst.coef_BS>>>
         (workspace.x_array, workspace.t_array, workspace.alpha_gpu,
          workspace.x_mult, workspace.newton_t_mult, workspace.one_minor_t,
          workspace.path_idx, workspace.x_t_idx_mult, n_path_continuous,
          inst.dim, workspace.n_predictor);

      max_x_double_kernel_align<<<max_grid, inst.predict_BS>>>
         (workspace.x_mult, inst.dim, n_path_continuous, workspace.max_x_gpu,
          workspace.workspace_size, workspace.x_t_idx_mult, 
          workspace.path_idx);

      if(debug)
      {
         cudaMemcpy(workspace.max_x_host, workspace.max_x_gpu,
                    n_path*sizeof(double), cudaMemcpyDeviceToHost);
         for(int path_idx=0; path_idx<n_path; path_idx++)
         {
            if(path_idx == path_idx_test || debug_all)
            {
               std::cout << "       max_x_gpu " << path_idx << "= "
                         << workspace.max_x_host[path_idx] << std::endl;
            }
         }
      }
      eval_mult(workspace, inst);
      // inst.n_eval_GPU++;

      /*
        GT* tmp_f_val = workspace.matrix_horizontal_mult + inst.n_eq*inst.dim;
        max_relative_double_kernel<<<max_grid, inst.predict_BS>>>
           (tmp_f_val, inst.n_eq, n_path_continuous, workspace.max_f_val_gpu,
            workspace.r_max_f_val_gpu, workspace.max_x_gpu,
            workspace.n_matrix, workspace.path_idx);
       */
      max_relative_double_kernel3<<<max_grid, inst.predict_BS>>>
         (workspace.f_val, inst.n_eq, n_path_continuous, workspace.n_path,
          workspace.path_idx, workspace.max_f_val_gpu,
          workspace.max_f_val_last_gpu, workspace.r_max_f_val_gpu, 
          workspace.max_x_gpu);

      if(debug)
      {
         cudaMemcpy(workspace.max_f_val_host, workspace.max_f_val_gpu,
                    n_path*sizeof(double), cudaMemcpyDeviceToHost);
         cudaMemcpy(workspace.r_max_f_val_host, workspace.r_max_f_val_gpu,
                    n_path*sizeof(double), cudaMemcpyDeviceToHost);

         for(int path_idx=0; path_idx<n_path; path_idx++)
         {
            if(path_idx == path_idx_test || debug_all)
            {
               std::cout << "       max_f_value " << path_idx << " = "
                         << workspace.max_f_val_host[path_idx]
                         << " " <<  workspace.r_max_f_val_host[path_idx]
                         << std::endl;
            }
         }
      }
      check_kernel<<<check_grid, inst.predict_BS>>>
         (workspace.max_f_val_gpu, workspace.r_max_f_val_gpu,
          workspace.max_f_val_last_gpu, workspace.newton_success, n_path,
          workspace.end_range, path_parameter.err_min_round_off, 
          path_parameter.err_min_round_off_refine);

      cudaMemcpy(workspace.newton_success_host, workspace.newton_success,
                 n_path*sizeof(int), cudaMemcpyDeviceToHost);

      n_path_continuous = 0;
      // std::cout << "newton_success" << std::endl;
      for(int path_idx=0; path_idx<n_path; path_idx++)
      {
         // std::cout << path_idx << " "
         //  << workspace.newton_success_host[path_idx] << std::endl;
         if(workspace.newton_success_host[path_idx] == 0)
         {
            workspace.path_idx_host[n_path_continuous] = path_idx;
            n_path_continuous += 1;
         }
      }
      if(n_path_continuous==0)
      {
         break;
      }
      workspace.n_path_continuous = n_path_continuous;
      max_grid = get_grid(n_path_continuous,inst.predict_BS,1);
      /*
        std::cout << "n_path_continuous = " << n_path_continuous << " : ";
        for(int path_idx=0; path_idx<n_path_continuous; path_idx++)
        {
           std::cout << workspace.path_idx_host[path_idx] << ", ";
        }
        std::cout << std::endl;
       */
      cudaMemcpy(workspace.path_idx, workspace.path_idx_host,
                 n_path*sizeof(int), cudaMemcpyHostToDevice);
   }
   dim3 check_grid2 = get_grid(n_path,inst.predict_BS,1);
   check_kernel<<<check_grid, inst.predict_BS>>>
      (workspace.max_f_val_gpu, workspace.r_max_f_val_gpu,
       workspace.max_delta_x_gpu, workspace.r_max_delta_x_gpu,
       workspace.path_success, workspace.newton_success,
       workspace.n_point_mult, workspace.x_t_idx_mult, workspace.n_array,
       workspace.workspace_size, n_path, workspace.end_range,
       path_parameter.err_min_round_off,
       path_parameter.err_min_round_off_refine);

   cudaMemcpy(workspace.newton_success_host, workspace.newton_success,
              n_path*sizeof(int), cudaMemcpyDeviceToHost);

   if(debug)
   {
      std::cout << workspace.newton_success_host[path_idx_test] << std::endl;
   }
   return true;
}

bool GPU_Newton
 ( CPUInstHom& hom, Parameter path_parameter, CT* cpu_sol0, CT cpu_t, 
   CT*& x_new, int n_path, int verbose )
{
   if(verbose > 0)
   {
      cout << "Newton ";
      cout << "max_it = " << path_parameter.max_it << endl;
      cout << "eps    = " << path_parameter.err_max_delta_x << endl;
   }
   // clock_t begin = clock();

   cuda_set();
   GPUInst inst(hom, n_path);
   GPUWorkspace workspace(inst.mon_pos_size,inst.n_coef,inst.n_constant,
      inst.n_eq,inst.dim,path_parameter.n_predictor,inst.alpha,
      inst.base_table_size);

   workspace.update_x_t_value(cpu_sol0, cpu_t);

   struct timeval start, end;
   long seconds, useconds;
   double timeMS_gpu;
   gettimeofday(&start, NULL);

   bool success = newton_single2(workspace, inst, path_parameter);

   x_new = workspace.get_x();

   gettimeofday(&end, NULL);
   seconds  = end.tv_sec  - start.tv_sec;
   useconds = end.tv_usec - start.tv_usec;
   timeMS_gpu = seconds*1000 + useconds/1000.0;

   if(verbose > 0)
   {
      std::cout.precision(8);
      cout << "Path GPU Newton    Time: "<< timeMS_gpu << endl;
   }
   // return success;
   return true;
}

#endif /* NEWTON_CU_ */
