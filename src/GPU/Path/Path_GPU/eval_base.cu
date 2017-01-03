__global__ void eval_base_table_kernel
 ( int dim, int* max_deg_base, int* base_table_start, GT* x, GT* base_table )
{
   int bidx = (gridDim.x*blockIdx.y+blockIdx.x)*blockDim.x;
   int tidx = threadIdx.x;
   int idx = bidx + tidx;

   if(idx < dim && max_deg_base[idx]>0)
   {
      GT* var_start = base_table + base_table_start[idx];
      GT var = x[idx];
      GT var_term = var;
      var_start[0] = var_term;
      for(int deg_idx=1; deg_idx<max_deg_base[idx]; deg_idx++)
      {
         var_term *= var;
         var_start[deg_idx] = var_term;
      }
   }
}

__global__ void eval_base_kernel
 ( int n_mon_base, unsigned short* pos, unsigned short* exp,
   int* mon_pos_start, int* base_table_start, GT* base_table, GT* coef )
{
   int bidx = (gridDim.x*blockIdx.y+blockIdx.x)*blockDim.x;
   int tidx = threadIdx.x;
   int idx = bidx + tidx;

   if(idx < n_mon_base) 
   {
      unsigned short* mon_pos = pos+mon_pos_start[idx];
      unsigned short* mon_exp = exp+mon_pos_start[idx];
      int base_start = *mon_exp++;
      int n_var = *mon_pos++;
      GT tmp_coef = coef[idx];
      for(int var_idx=base_start; var_idx<n_var; var_idx++)
      {
         tmp_coef *= base_table[base_table_start[mon_pos[var_idx]]
                   + mon_exp[var_idx]-2];
      }
      coef[idx] = tmp_coef;
   }
}
