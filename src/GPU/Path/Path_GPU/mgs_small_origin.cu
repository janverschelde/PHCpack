__global__ void mgs_small_normalize_kernel
 ( GT* v, GT* R, int dimR, int rows, int rowsLog2, int cols, 
   int pivot, int workspace_size=0 );

__global__ void mgs_small_normalize_kernel1
 ( GT* v, GT* R, int dimR, int rows, int rowsLog2, int cols,
   int pivot, int workspace_size=0 );

__global__ void mgs_small_normalize_kernel1_idx
 ( GT* v, GT* R, int dimR, int rows, int rowsLog2, int cols,
   int pivot, int workspace_size, int n_matrix, int n_matrix_R,
   int* path_idx );

__global__ void mgs_small_reduce_kernel
 ( GT* v, GT* R, int dimR, int rows, int rowsLog2, int cols,
   int pivot, int workspace_size=0 );

__global__ void mgs_small_reduce_kernel1
 ( GT* v, GT* R, int dimR, int rows, int rowsLog2, int cols,
   int pivot, int workspace_size=0 );

__global__ void mgs_small_reduce_kernel16
 ( GT* v, GT* R, int dimR, int rows, int rowsLog2, int cols,
   int pivot, int workspace_size=0 );

__global__ void mgs_small_reduce_kernel16_idx
 ( GT* v, GT* R, int dimR, int rows, int rowsLog2, int cols,
   int pivot, int workspace_size, int n_matrix, int n_matrix_R,
   int* path_idx );

__global__ void mgs_small_backsubstitution_kernel
 ( GT* R, GT*sol0, int dim, int workspace_size=0 );

__global__ void mgs_small_backsubstitution_kernel1
 ( GT* R, GT*sol0, int dim, int workspace_size=0 );

__global__ void mgs_small_backsubstitution_kernel1_idx
 ( GT* R, GT*sol0, int dim, int workspace_size, int n_matrix_R,
   int* path_idx );

__global__ void mgs_small_backsubstitution_max_kernel
 ( GT* R, GT* sol0, int dim, int dimLog2, double* max_delta_x );

__global__ void mgs_small_kernel
 ( GT* V, GT* R, GT* sol0, int rows, int rowsLog2, int half_size_init,
   int cols, int workspace_size );

__global__ void mgs_small_kernel_idx
 ( GT* V, GT* R, GT* sol0, int rows, int rowsLog2, int half_size_init,
   int cols, int workspace_size, int* path_idx );

__global__ void mgs_small_kernel_mult
 ( GT* V, GT* R, GT* sol, int rows, int rowsLog2, int half_size_init,
   int cols );

void mgs_small_mult
 ( GT* V, GT* R, GT* sol, int rows, int cols, int n_path )
{
   // int rows = dim;

   int rowsLog2 = log2ceil(rows);// ceil for sum reduction

   // int cols = dim + 1;
   // std::cout << "rows = " << rows << " cols = " << cols 
   // << " rowsLog2 = " << rowsLog2 << " n_path = " << n_path << std::endl;

   int half_size_init = 1 << (rowsLog2-1);
   mgs_small_kernel_mult<<<n_path,32>>>
      (V,R,sol,rows,rowsLog2-1,half_size_init,cols);
}

void mgs_small_idx
 ( GT* V, GT* R, GT* sol, int rows, int cols, int workspace_size,
   int n_path, int* path_idx )
{
   // int rows = dim;

   int rowsLog2 = log2ceil(rows);// ceil for sum reduction

   // int cols = dim + 1;

   // std::cout << "rows = " << rows << " cols = " << cols 
   //  << " rowsLog2 = " << rowsLog2 << "n_path = " << n_path << std::endl;

   int half_size_init = 1 << (rowsLog2-1);
   mgs_small_kernel_idx<<<n_path,32>>>
      (V,R,sol,rows,rowsLog2-1,half_size_init,cols,workspace_size, path_idx);
}

void mgs_small
 ( GT* V, GT* R, GT* sol, int rows, int cols, int workspace_size=0,
   int n_path=1 )
{
   // int rows = dim;

   int rowsLog2 = log2ceil(rows);// ceil for sum reduction

   // int cols = dim + 1;

   // std::cout << "rows = " << rows << " cols = " << cols 
   //  << " rowsLog2 = " << rowsLog2 << "n_path = " << n_path << std::endl;

   int half_size_init = 1 << (rowsLog2-1);
   mgs_small_kernel<<<n_path,32>>>
      (V,R,sol,rows,rowsLog2-1,half_size_init,cols,workspace_size);
}

void mgs_small11
 ( GT* V, GT* R, GT* sol, int rows, int cols, int workspace_size=0,
   int n_path=1 )
{
   int BS = rows; // XXX Temperary solution
   // int rows = dim;
   int rowsLog2 = log2ceil(rows);// ceil for sum reduction
   // int cols = dim + 1;

   std::cout << "rows = " << rows << " cols = " << cols
             << " rowsLog2 = " << rowsLog2 << std::endl;

   std::cout << "n_path = " << n_path << std::endl;
   int half_size_init = 1 << (rowsLog2-1);
   mgs_small_kernel<<<n_path,BS>>>
      (V,R,sol,rows,rowsLog2-1, half_size_init,cols,workspace_size);
}

void mgs_small1
 ( GT* V, GT* R, GT* sol, int rows, int cols, int workspace_size=0,
   int n_path=1 )
{
   int BS = rows; // XXX Temperary solution

   // int rows = dim;

   int rowsLog2 = log2ceil(rows);// ceil for sum reduction
   int dimR = cols*(cols+1)/2;
   // int cols = dim + 1;

   /*
     std::cout << "rows = " << rows << " cols = " << cols
               << " rowsLog2 = " << rowsLog2 << std::endl;*/
     int half_size_init = 1 << (rowsLog2-1);

     std::cout << "n_path = " << n_path << std::endl;

     for(int piv=0; piv<cols-1; piv++)
     {
        mgs_small_normalize_kernel1<<<n_path,BS>>>
        (V,R,dimR,rows,half_size_init,cols,piv, workspace_size);
        mgs_small_reduce_kernel16<<<n_path,32>>>
        (V,R,dimR,rows,half_size_init,cols,piv, workspace_size);
     }
     mgs_small_backsubstitution_kernel1<<<n_path,cols-1>>>
        (R,sol,cols-1, workspace_size);
}

void mgs_small1_idx
 ( GT* V, GT* R, GT* sol, int rows, int cols, int workspace_size,
   int n_matrix, int n_matrix_R, int n_path, int* path_idx )
{
   int BS = rows; // XXX Temperary solution
   // int rows = dim;
   int rowsLog2 = log2ceil(rows);// ceil for sum reduction
   int dimR = cols*(cols+1)/2;
   // int cols = dim + 1;

   /*
     std::cout << "rows = " << rows << " cols = " << cols
               << " rowsLog2 = " << rowsLog2 << std::endl;
    */
     int half_size_init = 1 << (rowsLog2-1);

     // std::cout << "n_path = " << n_path << std::endl;

     for(int piv=0; piv<cols-1; piv++) 
     {
        mgs_small_normalize_kernel1_idx<<<n_path,BS>>>
           (V,R,dimR,rows,half_size_init,cols,piv, workspace_size,
            n_matrix, n_matrix_R, path_idx);
        mgs_small_reduce_kernel16_idx<<<n_path,32>>>
           (V,R,dimR,rows,half_size_init,cols,piv, workspace_size,
            n_matrix, n_matrix_R, path_idx);
     }
     mgs_small_backsubstitution_kernel1_idx<<<n_path,cols-1>>>
        (R,sol,cols-1, workspace_size, n_matrix_R, path_idx);
}

__global__ void mgs_small_kernel_idx
 ( GT* V, GT* R, GT* sol, int rows, int rowsLog2, int half_size_init, 
   int cols, int workspace_size, int* path_idx_mult )
{
   /*
     __shared__ GT V_shared[272]; // contains pivot column
     __shared__ GT R_shared[153]; // contains pivot column
     __shared__ T prd[16];        // for norm of the pivot
     __shared__ GT tmp_col[32];   // for norm of the pivot
    */

   __shared__ GT V_shared[110]; // contains pivot column
   __shared__ GT R_shared[66];  // contains pivot column
   __shared__ T prd[10];        // for norm of the pivot
   __shared__ GT tmp_col[32];   // for norm of the pivot

   /*
     __shared__ GT V_shared[72]; // contains pivot column
     __shared__ GT R_shared[45]; // contains pivot column
     __shared__ T prd[8];        // for norm of the pivot
     __shared__ GT tmp_col[32];  // for norm of the pivot
    */

   int path_idx = path_idx_mult[blockIdx.x];

   V += path_idx*workspace_size;
   // R += blockIdx.x*workspace_size;

   // load matrix
   int t_idx = threadIdx.x;
   int v_idx=t_idx;

   int n_cols = 32/rows;
   int col_idx_rnd = t_idx/rows;
   int col_t_idx = t_idx - rows*col_idx_rnd;
   int n_rnd = (cols-1)/n_cols + 1;

   int col_idx = col_idx_rnd;
   for(int rnd_idx=0; rnd_idx<n_rnd; rnd_idx++)
   {
      if(col_idx_rnd < n_cols && col_idx < cols)
      {
         V_shared[v_idx] = V[v_idx];
      }
      v_idx += n_cols*rows;
      col_idx += n_cols;
   }

   // QR
   GT* piv = V_shared;
   v_idx=t_idx;
   for(int pivot=0; pivot<cols-1; pivot++)
   {
      // normalize
      if(t_idx<rows)
      {
         prd[t_idx] = piv[t_idx].real*piv[t_idx].real
            + piv[t_idx].imag*piv[t_idx].imag;
         if(t_idx + half_size_init < rows)
         {
            prd[t_idx] = prd[t_idx] + prd[t_idx+half_size_init];
         }
         if(t_idx < 4)
         {
            prd[t_idx] = prd[t_idx] + prd[t_idx+4];
         }
         if(t_idx < 2)
         {
            prd[t_idx] = prd[t_idx] + prd[t_idx+2];
         }
         if(t_idx == 0)
         {
            prd[0] = prd[0] + prd[1];
            prd[0] = sqrt(prd[0]);
            int indR = r_pos(pivot, pivot, cols);
            // (dimR-1) - (pivot*(pivot+1))/2 - (b*(b+1))/2 - b*(pivot+1);
            // R[indR].init_imag();
            // R[indR].real = prd[0];
            R_shared[indR].real = prd[0];
            R_shared[indR].init_imag();
         }
      }
      if(t_idx<rows)
      {
         piv[t_idx] /= prd[0];
      }
      // reduce
      int col_idx = col_idx_rnd + pivot + 1;
      int n_rnd = (cols -pivot-2)/n_cols + 1;
      int sum_pos = col_idx_rnd*rows;

      GT* shv = V_shared + col_idx*rows;
      for(int rnd_idx=0; rnd_idx<n_rnd; rnd_idx++)
      {
         if(col_idx < cols && col_idx_rnd<n_cols)
         {
            tmp_col[t_idx] = piv[col_t_idx].adj_multiple(shv[col_t_idx]);
            if(col_t_idx + half_size_init < rows)
            {
               tmp_col[t_idx] = tmp_col[t_idx] + tmp_col[t_idx+half_size_init];
            }
            if(col_t_idx < 4)
            {
               tmp_col[t_idx] = tmp_col[t_idx] + tmp_col[t_idx+4];
            }
            if(col_t_idx < 2)
            {
               tmp_col[t_idx] = tmp_col[t_idx] + tmp_col[t_idx+2];
            }
            if(col_t_idx < 1)
            {
               tmp_col[t_idx] = tmp_col[t_idx] + tmp_col[t_idx+1];
            }
            shv[col_t_idx] = shv[col_t_idx] - tmp_col[sum_pos]*piv[col_t_idx];
            // V[block*rows+t_idx] = shv[t_idx];
            if(col_t_idx == 0)
            {
               int indR = r_pos(pivot, col_idx, cols);
               // (dimR-1) - (pivot*(pivot+1))/2 - (b*(b+1))/2 - b*(pivot+1);
               // R[indR] = tmp_col[0];
               R_shared[indR] = tmp_col[sum_pos];
            }
            shv += rows*n_cols;
            col_idx += n_cols;
         }
      }
      piv += rows;
      v_idx += rows;
   }
   // back substitution
   if(t_idx<cols-1)
   {
      GT update = R_shared[t_idx];
      int dim = cols -1;
      GT* Rcl;
      for(int k=dim-1; k>=0; k--) // compute k-th component of solution
      {
         if(t_idx < k+1)
         {
            int ind = (dim - k)*(dim + 3 + k)/2;
            Rcl = R_shared+ind;
         }
         if(t_idx == k) tmp_col[t_idx] = update/Rcl[t_idx];
         // all other threads wait
         if(t_idx < k) update = update - tmp_col[k]*Rcl[t_idx]; // update
      }
      sol += path_idx*workspace_size;
      sol[t_idx] = tmp_col[t_idx];
   }
};

__global__ void mgs_small_kernel_mult
 ( GT* V, GT* R, GT* sol, int rows, int rowsLog2, int half_size_init, 
   int cols )
{
   __shared__ GT V_shared[110];// contains pivot column
   __shared__ GT R_shared[66];// contains pivot column
   __shared__ T prd[10];// for norm of the pivot
   __shared__ GT tmp_col[32];// for norm of the pivot

   int path_idx = blockIdx.x;
   int n_path = gridDim.x;

   V += path_idx;
   // R += blockIdx.x*workspace_size;

   // load matrix
   int t_idx = threadIdx.x;
   int v_idx = t_idx;

   int n_cols = 32/rows;
   int col_idx_rnd = t_idx/rows;
   int col_t_idx = t_idx - rows*col_idx_rnd;
   int n_rnd = (cols-1)/n_cols + 1;

   int col_idx = col_idx_rnd;
   for(int rnd_idx=0; rnd_idx<n_rnd; rnd_idx++)
   {
      if(col_idx_rnd < n_cols && col_idx < cols)
      {
         V_shared[v_idx] = V[v_idx*n_path];
      }
      v_idx += n_cols*rows;
      col_idx += n_cols;
   }
   // QR
   GT* piv = V_shared;
   v_idx=t_idx;
   for(int pivot=0; pivot<cols-1; pivot++)
   {
      // normalize
      if(t_idx<rows)
      {
         prd[t_idx] = piv[t_idx].real*piv[t_idx].real
            + piv[t_idx].imag*piv[t_idx].imag;
         if(t_idx + half_size_init < rows)
         {
            prd[t_idx] = prd[t_idx] + prd[t_idx+half_size_init];
         }
         if(t_idx < 4)
         {
            prd[t_idx] = prd[t_idx] + prd[t_idx+4];
         }
         if(t_idx < 2)
         {
            prd[t_idx] = prd[t_idx] + prd[t_idx+2];
         }
         if(t_idx == 0)
         {
            prd[0] = prd[0] + prd[1];
            prd[0] = sqrt(prd[0]);
            int indR = r_pos(pivot, pivot, cols);
            // (dimR-1) - (pivot*(pivot+1))/2 - (b*(b+1))/2 - b*(pivot+1);
            // R[indR].init_imag();
            // R[indR].real = prd[0];
            R_shared[indR].real = prd[0];
            R_shared[indR].init_imag();
         }
      }
      if(t_idx<rows)
      {
         piv[t_idx] /= prd[0];
      }
      // reduce
      int col_idx = col_idx_rnd + pivot + 1;
      int n_rnd = (cols -pivot-2)/n_cols + 1;
      int sum_pos = col_idx_rnd*rows;

      GT* shv = V_shared + col_idx*rows;
      for(int rnd_idx=0; rnd_idx<n_rnd; rnd_idx++)
      {
         if(col_idx < cols && col_idx_rnd<n_cols)
         {
            tmp_col[t_idx] = piv[col_t_idx].adj_multiple(shv[col_t_idx]);
            if(col_t_idx + half_size_init < rows)
            {
               tmp_col[t_idx] = tmp_col[t_idx] + tmp_col[t_idx+half_size_init];
            }
            if(col_t_idx < 4)
            {
               tmp_col[t_idx] = tmp_col[t_idx] + tmp_col[t_idx+4];
            }
            if(col_t_idx < 2)
            {
               tmp_col[t_idx] = tmp_col[t_idx] + tmp_col[t_idx+2];
            }
            if(col_t_idx < 1)
            {
               tmp_col[t_idx] = tmp_col[t_idx] + tmp_col[t_idx+1];
            }
            shv[col_t_idx] = shv[col_t_idx] - tmp_col[sum_pos]*piv[col_t_idx];
            // V[block*rows+t_idx] = shv[t_idx];
            if(col_t_idx == 0)
            {
               int indR = r_pos(pivot, col_idx, cols);
               // (dimR-1) - (pivot*(pivot+1))/2 - (b*(b+1))/2 - b*(pivot+1);
               // R[indR] = tmp_col[0];
               R_shared[indR] = tmp_col[sum_pos];
            }
            shv += rows*n_cols;
            col_idx += n_cols;
         }
      }
      piv += rows;
      v_idx += rows;
   }
   // back substitution
   if(t_idx<cols-1)
   {
      GT update = R_shared[t_idx];
      int dim = cols -1;
      GT* Rcl;
      for(int k=dim-1; k>=0; k--)  // compute k-th component of solution
      {
         if(t_idx < k+1)
         {
            int ind = (dim - k)*(dim + 3 + k)/2;
            Rcl = R_shared+ind;
         }
         if(t_idx == k) tmp_col[t_idx] = update/Rcl[t_idx];
         // all other threads wait
         if(t_idx < k) update = update - tmp_col[k]*Rcl[t_idx]; // update
      }
      sol += path_idx;
      sol[t_idx*n_path] = tmp_col[t_idx];
   }
};

__global__ void mgs_small_kernel
 ( GT* V, GT* R, GT* sol, int rows, int rowsLog2, int half_size_init,
   int cols, int workspace_size )
{
   __shared__ GT V_shared[110];// contains pivot column
   __shared__ GT R_shared[66];// contains pivot column
   __shared__ T prd[10];// for norm of the pivot
   __shared__ GT tmp_col[32];// for norm of the pivot

   V += blockIdx.x*workspace_size;
   // R += blockIdx.x*workspace_size;

   // load matrix
   int t_idx = threadIdx.x;
   int v_idx=t_idx;

   int n_cols = 32/rows;
   int col_idx_rnd = t_idx/rows;
   int col_t_idx = t_idx - rows*col_idx_rnd;
   int n_rnd = (cols-1)/n_cols + 1;

   int col_idx = col_idx_rnd;
   for(int rnd_idx=0; rnd_idx<n_rnd; rnd_idx++)
   {
      if(col_idx_rnd < n_cols && col_idx < cols)
      {
         V_shared[v_idx] = V[v_idx];
      }
      v_idx += n_cols*rows;
      col_idx += n_cols;
   }
   // QR
   GT* piv = V_shared;
   v_idx=t_idx;
   for(int pivot=0; pivot<cols-1; pivot++)
   {
      // normalize
      if(t_idx<rows)
      {
         prd[t_idx] = piv[t_idx].real*piv[t_idx].real
            + piv[t_idx].imag*piv[t_idx].imag;
         if(t_idx + half_size_init < rows) 
         {
            prd[t_idx] = prd[t_idx] + prd[t_idx+half_size_init];
         }
         if(t_idx < 4)
         {
            prd[t_idx] = prd[t_idx] + prd[t_idx+4];
         }
         if(t_idx < 2)
         {
            prd[t_idx] = prd[t_idx] + prd[t_idx+2];
         }
         if(t_idx == 0)
         {
            prd[0] = prd[0] + prd[1];
            prd[0] = sqrt(prd[0]);
            int indR = r_pos(pivot, pivot, cols);
            // (dimR-1) - (pivot*(pivot+1))/2 - (b*(b+1))/2 - b*(pivot+1);
            // R[indR].init_imag();
            // R[indR].real = prd[0];
            R_shared[indR].real = prd[0];
            R_shared[indR].init_imag();
         }
      }
      if(t_idx<rows)
      {
         piv[t_idx] /= prd[0];
      }
      // reduce
      int col_idx = col_idx_rnd + pivot + 1;
      int n_rnd = (cols -pivot-2)/n_cols + 1;
      int sum_pos = col_idx_rnd*rows;

      GT* shv = V_shared + col_idx*rows;
      for(int rnd_idx=0; rnd_idx<n_rnd; rnd_idx++)
      {
         if(col_idx < cols && col_idx_rnd<n_cols)
         {
            tmp_col[t_idx] = piv[col_t_idx].adj_multiple(shv[col_t_idx]);
            if(col_t_idx + half_size_init < rows)
            {
               tmp_col[t_idx] = tmp_col[t_idx] + tmp_col[t_idx+half_size_init];
            }
            if(col_t_idx < 4)
            {
               tmp_col[t_idx] = tmp_col[t_idx] + tmp_col[t_idx+4];
            }
            if(col_t_idx < 2)
            {
               tmp_col[t_idx] = tmp_col[t_idx] + tmp_col[t_idx+2];
            }
            if(col_t_idx < 1)
            {
               tmp_col[t_idx] = tmp_col[t_idx] + tmp_col[t_idx+1];
            }
            shv[col_t_idx] = shv[col_t_idx] - tmp_col[sum_pos]*piv[col_t_idx];
            // V[block*rows+t_idx] = shv[t_idx];
            if(col_t_idx == 0)
            {
               int indR = r_pos(pivot, col_idx, cols);
               // (dimR-1) - (pivot*(pivot+1))/2 - (b*(b+1))/2 - b*(pivot+1);
               // R[indR] = tmp_col[0];
               R_shared[indR] = tmp_col[sum_pos];
            }
            shv += rows*n_cols;
            col_idx += n_cols;
         }
      }
      piv += rows;
      v_idx += rows;
   }
   // back substitution
   if(t_idx<cols-1)
   {
      GT update = R_shared[t_idx];
      int dim = cols -1;
      GT* Rcl;
      for(int k=dim-1; k>=0; k--)  // compute k-th component of solution
      {
         if(t_idx < k+1)
         {
            int ind = (dim - k)*(dim + 3 + k)/2;
            Rcl = R_shared+ind;
         }
         if(t_idx == k) tmp_col[t_idx] = update/Rcl[t_idx];
         // all other threads wait
         if(t_idx < k) update = update - tmp_col[k]*Rcl[t_idx]; // update
      }
      sol += blockIdx.x*workspace_size;
      sol[t_idx] = tmp_col[t_idx];
   }
};

__global__ void mgs_small_kernel1
 ( GT* V, GT* R, GT* sol, int rows, int rowsLog2, int half_size_init,
   int cols, int workspace_size )
{
   __shared__ GT V_shared[110];// contains pivot column
   __shared__ GT R_shared[66];// contains pivot column
   __shared__ T prd[10];// for norm of the pivot
   __shared__ GT tmp_col[10];// for norm of the pivot

   V += blockIdx.x*workspace_size;
   // R += blockIdx.x*workspace_size;

   int t_idx = threadIdx.x;
   int v_idx=t_idx;
   for(int i=0; i<cols; i++)
   {
      V_shared[v_idx] = V[v_idx];
      v_idx += rows;
   }
   GT* piv = V_shared;
   v_idx=t_idx;
   for(int pivot=0; pivot<cols-1; pivot++)
   {
      // Normalization
      prd[t_idx] = piv[t_idx].real*piv[t_idx].real
         + piv[t_idx].imag*piv[t_idx].imag;
      if(t_idx + half_size_init < rows)
      {
         prd[t_idx] = prd[t_idx] + prd[t_idx+half_size_init];
      }
      if(t_idx < 4)
      {
         prd[t_idx] = prd[t_idx] + prd[t_idx+4];
      }
      if(t_idx < 2)
      {
         prd[t_idx] = prd[t_idx] + prd[t_idx+2];
      }
      if(t_idx == 0)
      {
         prd[0] = prd[0] + prd[1];
         prd[0] = sqrt(prd[0]);
         int indR = r_pos(pivot, pivot, cols);
         // (dimR-1) - (pivot*(pivot+1))/2 - (b*(b+1))/2 - b*(pivot+1);
         // R[indR].init_imag();
         // R[indR].real = prd[0];
         R_shared[indR].real = prd[0];
         R_shared[indR].init_imag();
      }
      __syncthreads();
      piv[t_idx] /= prd[0];

      // Reduce
      GT* shv = V_shared + (pivot+1)*rows;
      for(int block=pivot+1; block<cols; block++)
      {
         tmp_col[t_idx] = piv[t_idx].adj_multiple(shv[t_idx]);
         __syncthreads();
         if(t_idx + half_size_init < rows)
         {
            tmp_col[t_idx] = tmp_col[t_idx] + tmp_col[t_idx+half_size_init];
         }
         if(t_idx < 4)
         {
            tmp_col[t_idx] = tmp_col[t_idx] + tmp_col[t_idx+4];
         }
         if(t_idx < 2)
         {
            tmp_col[t_idx] = tmp_col[t_idx] + tmp_col[t_idx+2];
         }
         if(t_idx < 1)
         {
            tmp_col[t_idx] = tmp_col[t_idx] + tmp_col[t_idx+1];
         }
         shv[t_idx] = shv[t_idx] - tmp_col[0]*piv[t_idx];
         // V[block*rows+t_idx] = shv[t_idx];
         shv += rows;

         int indR = r_pos(pivot, block, cols);
         // (dimR-1) - (pivot*(pivot+1))/2 - (b*(b+1))/2 - b*(pivot+1);
         if(t_idx == 0)
         {
            // R[indR] = tmp_col[0];
            R_shared[indR] = tmp_col[0];
         }
      }
      piv += rows;
      v_idx += rows;
   }
   // Back substitution
   GT update = R_shared[t_idx];
   int dim = cols -1;
   GT* Rcl;
   for(int k=dim-1; k>=0; k--)  // compute k-th component of solution
   {
      if(t_idx < k+1)
      {
         int ind = (dim - k)*(dim + 3 + k)/2;
         Rcl = R_shared+ind;
      }
      if(t_idx == k) tmp_col[t_idx] = update/Rcl[t_idx];
      // all other threads wait
      if(t_idx < k) update = update - tmp_col[k]*Rcl[t_idx]; // update
   }
   sol += blockIdx.x*workspace_size;
   sol[t_idx] = tmp_col[t_idx];
};

void mgs_small2
 ( GT* V, GT* R, GT* sol, int rows, int cols, int workspace_size=0,
   int n_path=1 )
{
   int BS = rows; // XXX Temperary solution
   // int rows = dim;
   int rowsLog2 = log2ceil(rows);// ceil for sum reduction
   int dimR = cols*(cols+1)/2;
   // int cols = dim + 1;

   /*
     std::cout << "rows = " << rows << " cols = " << cols
               << " rowsLog2 = " << rowsLog2 << std::endl;
    */

   std::cout << "n_path = " << n_path << std::endl;

   for(int piv=0; piv<cols-1; piv++)
   {
      mgs_small_normalize_kernel<<<n_path,BS>>>
         (V,R,dimR,rows,rowsLog2,cols,piv, workspace_size);
      dim3 NB(cols-piv-1,1,n_path);
      mgs_small_reduce_kernel<<<NB,BS>>>
         (V,R,dimR,rows,rowsLog2,cols,piv, workspace_size);
   }
   mgs_small_backsubstitution_kernel<<<n_path,cols-1>>>
      (R,sol,cols-1, workspace_size);
}

void mgs_small_with_delta
 ( GT* V, GT* R, GT* sol, int rows, int cols, double* max_delta_x )
{
   int BS = rows; // XXX Temperary solution
   // int rows = dim;
   int rowsLog2 = log2ceil(rows);// ceil for sum reduction
   int dimR = cols*(cols+1)/2;
   int dimLog2 = log2ceil(cols-1);// ceil for sum reduction

   for(int piv=0; piv<cols-1; piv++)
   {
      mgs_small_normalize_kernel<<<1,BS>>>
        (V,R,dimR,rows,rowsLog2,cols,piv);
      mgs_small_reduce_kernel<<<cols-piv-1,BS>>>
        (V,R,dimR,rows,rowsLog2,cols,piv);
   }
   mgs_small_backsubstitution_max_kernel<<<1,cols-1>>>
      (R,sol,cols-1, dimLog2, max_delta_x);
}

/*
  int GPU_MGS
   ( const CPUInstHom& hom, CT*& sol_gpu, CT*& matrix_gpu_q,
     CT*& matrix_gpu_r, int n_predictor, CT* V, int n_path )
  {
     cout << "GPU Eval" << endl;

     // CUDA configuration
     cuda_set();

     GPUInst inst(hom, n_path);
     GPUWorkspace workspace
        (inst.n_workspace, inst.n_coef, inst.n_constant, inst.n_eq,
         inst.dim, n_predictor, inst.alpha);

     workspace.init_V_value(V);

     mgs_small(workspace.V, workspace.R, workspace.sol, hom.n_eq, hom.dim+1);

     sol_gpu = workspace.get_sol();
     matrix_gpu_q = workspace.get_matrix();
     matrix_gpu_r = workspace.get_matrix_r();

     cudaDeviceReset();
     return 0;
  }
 */

__global__ void mgs_small_normalize_kernel
 ( GT* v, GT* R, int dimR, int rows, int rowsLog2, int cols,
   int pivot, int workspace_size )
{
   // int b = blockIdx.x;
   int t_idx = threadIdx.x;
   // int block = b+pivot;    // column for reduction w.r.t. pivot
   // int i = block*rows + t_idx;// idx
   int L = pivot*rows + t_idx;

   __shared__ GT piv[shmemsize];// contains pivot column
   __shared__ T prd[shmemsize];// for norm of the pivot

   v += blockIdx.x*workspace_size;
   R += blockIdx.x*workspace_size;

   piv[t_idx] = v[L];
   prd[t_idx] = piv[t_idx].real*piv[t_idx].real
      + piv[t_idx].imag*piv[t_idx].imag;
   __syncthreads();

   rowsLog2 -= 1;
   int half_size = 1 << (rowsLog2); // sum for the norm
   if(t_idx + half_size < rows)
   {
      prd[t_idx] = prd[t_idx] + prd[t_idx+half_size];
   }
   for(int k=0; k < rowsLog2; k++)
   {
      if(half_size > 16)
      {
         __syncthreads();
      }
      half_size /= 2;
      if(t_idx < half_size)
      {
         prd[t_idx] = prd[t_idx] + prd[t_idx+half_size];
      }
   }
   if(t_idx == 0) prd[0] = sqrt(prd[0]);
   __syncthreads();

   piv[t_idx] /= prd[0];
   v[L] = piv[t_idx];
   if(t_idx == 0)
   {
      int indR = r_pos(pivot, pivot, cols);
      // (dimR-1) - (pivot*(pivot+1))/2 - (b*(b+1))/2 - b*(pivot+1);
      R[indR].init_imag();
      R[indR].real = prd[0];
   }
}

__global__ void mgs_small_normalize_kernel1
 ( GT* v, GT* R, int dimR, int rows, int half_size, int cols, int pivot,
   int workspace_size )
{
   // int b = blockIdx.x;
   int t_idx = threadIdx.x;
   // int block = b+pivot;        // column for reduction w.r.t. pivot
   // int i = block*rows + t_idx; // idx
   int L = pivot*rows + t_idx;

   __shared__ GT piv[32];// contains pivot column
   __shared__ T prd[32];// for norm of the pivot

   v += blockIdx.x*workspace_size;
   R += blockIdx.x*workspace_size;

   piv[t_idx] = v[L];
   prd[t_idx] = piv[t_idx].real*piv[t_idx].real
      + piv[t_idx].imag*piv[t_idx].imag;

   if(t_idx + half_size < rows)
   {
      prd[t_idx] = prd[t_idx] + prd[t_idx+half_size];
   }
   if(t_idx < 4)
   {
      prd[t_idx] = prd[t_idx] + prd[t_idx+4];
   }
   if(t_idx < 2)
   {
      prd[t_idx] = prd[t_idx] + prd[t_idx+2];
   }
   if(t_idx < 1)
   {
      prd[t_idx] = prd[t_idx] + prd[t_idx+1];
      prd[0] = sqrt(prd[0]);
   }
   piv[t_idx] /= prd[0];
   v[L] = piv[t_idx];
   if(t_idx == 0)
   {
      int indR = r_pos(pivot, pivot, cols);
      // (dimR-1) - (pivot*(pivot+1))/2 - (b*(b+1))/2 - b*(pivot+1);
      R[indR].init_imag();
      R[indR].real = prd[0];
   }
}

__global__ void mgs_small_normalize_kernel1_idx
 ( GT* v, GT* R, int dimR, int rows, int half_size, int cols, 
   int pivot, int workspace_size, int n_matrix, int n_matrix_R,
   int* path_idx )
{
   // int b = blockIdx.x;
   int t_idx = threadIdx.x;
   // int block = b+pivot;        // column for reduction w.r.t. pivot
   // int i = block*rows + t_idx; // idx
   int L = pivot*rows + t_idx;

   __shared__ GT piv[32];  // contains pivot column
   __shared__ T prd[32];   // for norm of the pivot

   v += path_idx[blockIdx.x]*n_matrix;
   R += blockIdx.x*n_matrix_R;

   piv[t_idx] = v[L];
   prd[t_idx] = piv[t_idx].real*piv[t_idx].real
      + piv[t_idx].imag*piv[t_idx].imag;

   if(t_idx + half_size < rows)
   {
      prd[t_idx] = prd[t_idx] + prd[t_idx+half_size];
   }
   if(t_idx < 4)
   {
      prd[t_idx] = prd[t_idx] + prd[t_idx+4];
   }
   if(t_idx < 2)
   {
      prd[t_idx] = prd[t_idx] + prd[t_idx+2];
   }
   if(t_idx < 1)
   {
      prd[t_idx] = prd[t_idx] + prd[t_idx+1];
      prd[0] = sqrt(prd[0]);
   }
   piv[t_idx] /= prd[0];
   v[L] = piv[t_idx];
   if(t_idx == 0)
   {
      int indR = r_pos(pivot, pivot, cols);
      // (dimR-1) - (pivot*(pivot+1))/2 - (b*(b+1))/2 - b*(pivot+1);
      R[indR].init_imag();
      R[indR].real = prd[0];
   }
}

__global__ void mgs_small_reduce_kernel1
 ( GT* v, GT* R, int dimR, int rows, int rowsLog2, int cols,
   int pivot, int workspace_size )
{
   int b = blockIdx.x;
   int t_idx = threadIdx.x;
   int L = pivot*rows + t_idx;

   __shared__ GT piv[32];// contains pivot column
   __shared__ GT shv[2][32];// for the reduction

   v += b*workspace_size;
   R += b*workspace_size;

   piv[t_idx] = v[L];
   rowsLog2 -= 1;

   for(int block=pivot+1; block<cols; block++)
   {
      int i = block*rows + t_idx; // idx

      shv[0][t_idx] = v[i];
      shv[1][t_idx] = piv[t_idx].adj_multiple(shv[0][t_idx]);

      int indR = r_pos(pivot, block, cols);
      // (dimR-1) - (pivot*(pivot+1))/2 - (b*(b+1))/2 - b*(pivot+1);
      int rowsLog2_tmp = rowsLog2;
      int half_size = 1 << (rowsLog2);// sum for the norm
      if(t_idx + half_size < rows)
      {
         shv[1][t_idx] = shv[1][t_idx] + shv[1][t_idx+half_size];
      }
      for(int k=0; k <rowsLog2_tmp; k++)
      {
         half_size /= 2;
         if(t_idx < half_size)
         {
            shv[1][t_idx] = shv[1][t_idx] + shv[1][t_idx+half_size];
         }
      }
      shv[0][t_idx] = shv[0][t_idx] - shv[1][0]*piv[t_idx];
      v[i] = shv[0][t_idx];

      if(t_idx == 0) R[indR] = shv[1][0];
   }
}

__global__ void mgs_small_reduce_kernel161
 ( GT* v, GT* R, int dimR, int rows, int half_size, int cols,
   int pivot, int workspace_size )
{
   int b = blockIdx.x;
   int t_idx = threadIdx.x;
   int L = pivot*rows + t_idx;

   __shared__ GT piv[32];// contains pivot column
   __shared__ GT shv[2][32];// for the reduction

   v += b*workspace_size;
   R += b*workspace_size;

   piv[t_idx] = v[L];

   for(int block=pivot+1; block<cols; block++)
   {
      int i = block*rows + t_idx; // idx

      shv[0][t_idx] = v[i];
      shv[1][t_idx] = piv[t_idx].adj_multiple(shv[0][t_idx]);

      // int rowsLog2_tmp = rowsLog2;
      if(t_idx + half_size < rows)
      {
         shv[1][t_idx] = shv[1][t_idx] + shv[1][t_idx+half_size];
      }
      if(t_idx < 4)
      {
         shv[1][t_idx] = shv[1][t_idx] + shv[1][t_idx+4];
      }
      if(t_idx < 2)
      {
         shv[1][t_idx] = shv[1][t_idx] + shv[1][t_idx+2];
      }
      if(t_idx < 1)
      {
         shv[1][t_idx] = shv[1][t_idx] + shv[1][t_idx+1];
      }
      shv[0][t_idx] = shv[0][t_idx] - shv[1][0]*piv[t_idx];
      v[i] = shv[0][t_idx];

      int indR = r_pos(pivot, block, cols);
      // (dimR-1) - (pivot*(pivot+1))/2 - (b*(b+1))/2 - b*(pivot+1);
      if(t_idx == 0) R[indR] = shv[1][0];
   }
}

__global__ void mgs_small_reduce_kernel16
 ( GT* v, GT* R, int dimR, int rows, int half_size, int cols,
   int pivot, int workspace_size )
{
   int b = blockIdx.x;
   int t_idx = threadIdx.x;
   int L = pivot*rows + t_idx;

   __shared__ GT piv[16];// contains pivot column
   __shared__ GT shv[2][32];// for the reduction

   v += b*workspace_size;
   R += b*workspace_size;

   // load pivot
   if(t_idx < rows)
   {
      piv[t_idx] = v[L];
   }
   int n_cols = 32/rows;
   int col_idx_rnd = t_idx/rows;
   int col_t_idx = t_idx - rows*col_idx_rnd;

   int col_idx = col_idx_rnd + pivot + 1;
   int n_rnd = (cols -pivot-2)/n_cols + 1;
   int sum_pos = col_idx_rnd*rows;

   GT* tmp_v = v + (pivot+1)*rows;

   for(int rnd_idx=0; rnd_idx<n_rnd; rnd_idx++)
   {
      if(col_idx < cols && col_idx_rnd<n_cols)
      {
         shv[0][t_idx] = tmp_v[t_idx];
         shv[1][t_idx] = piv[col_t_idx].adj_multiple(shv[0][t_idx]);
         if(col_t_idx + half_size < rows)
         {
            shv[1][t_idx] = shv[1][t_idx] + shv[1][t_idx+half_size];
         }
         if(col_t_idx < 4) 
         {
            shv[1][t_idx] = shv[1][t_idx] +shv[1][t_idx+4];
         }
         if(col_t_idx < 2) 
         {
            shv[1][t_idx] = shv[1][t_idx] + shv[1][t_idx+2];
         }
         if(col_t_idx < 1)
         {
            shv[1][t_idx] = shv[1][t_idx] + shv[1][t_idx+1];
         }
         tmp_v[t_idx] = shv[0][t_idx] - shv[1][sum_pos]*piv[col_t_idx];
         // V[block*rows+t_idx] = shv[t_idx];
         if(col_t_idx == 0)
         {
            int indR = r_pos(pivot, col_idx, cols);
            // (dimR-1) - (pivot*(pivot+1))/2 - (b*(b+1))/2 - b*(pivot+1);
            // R[indR] = tmp_col[0];
            R[indR] = shv[1][sum_pos];
         }
         col_idx += n_cols;
         tmp_v += n_cols*rows;
      }
   }
}

__global__ void mgs_small_reduce_kernel16_idx
 ( GT* v, GT* R, int dimR, int rows, int half_size, int cols, int pivot,
   int workspace_size, int n_matrix, int n_matrix_R, int* path_idx )
{
   __shared__ GT piv[16];    // contains pivot column
   __shared__ GT shv[2][32]; // for the reduction

   int t_idx = threadIdx.x;
   int L = pivot*rows + t_idx;

   v += path_idx[blockIdx.x]*n_matrix;
   R += blockIdx.x*n_matrix_R;

   // load pivot
   if(t_idx < rows)
   {
      piv[t_idx] = v[L];
   }

   int n_cols = 32/rows;
   int col_idx_rnd = t_idx/rows;
   int col_t_idx = t_idx - rows*col_idx_rnd;

   int col_idx = col_idx_rnd + pivot + 1;
   int n_rnd = (cols -pivot-2)/n_cols + 1;
   int sum_pos = col_idx_rnd*rows;

   GT* tmp_v = v + (pivot+1)*rows;

   for(int rnd_idx=0; rnd_idx<n_rnd; rnd_idx++)
   {
      if(col_idx < cols && col_idx_rnd<n_cols)
      {
         shv[0][t_idx] = tmp_v[t_idx];
         shv[1][t_idx] = piv[col_t_idx].adj_multiple(shv[0][t_idx]);
         if(col_t_idx + half_size < rows) 
         {
            shv[1][t_idx] = shv[1][t_idx] + shv[1][t_idx+half_size];
         }
         if(col_t_idx < 4) 
         {
            shv[1][t_idx] = shv[1][t_idx] +shv[1][t_idx+4];
         }
         if(col_t_idx < 2)
         {
            shv[1][t_idx] = shv[1][t_idx] + shv[1][t_idx+2];
         }
         if(col_t_idx < 1) 
         {
            shv[1][t_idx] = shv[1][t_idx] + shv[1][t_idx+1];
         }
         tmp_v[t_idx] = shv[0][t_idx] - shv[1][sum_pos]*piv[col_t_idx];
         // V[block*rows+t_idx] = shv[t_idx];
         if(col_t_idx == 0)
         {
            int indR = r_pos(pivot, col_idx, cols);
            // (dimR-1) - (pivot*(pivot+1))/2 - (b*(b+1))/2 - b*(pivot+1);
            // R[indR] = tmp_col[0];
            R[indR] = shv[1][sum_pos];
         }
         col_idx += n_cols;
         tmp_v += n_cols*rows;
      }
   }
}

__global__ void mgs_small_reduce_kernel
 ( GT* v, GT* R, int dimR, int rows, int rowsLog2, int cols,
   int pivot, int workspace_size )
{
   int b = blockIdx.x+1;
   int t_idx = threadIdx.x;
   int block = b+pivot;        // column for reduction w.r.t. pivot
   int i = block*rows + t_idx; // idx
   int L = pivot*rows + t_idx;

   __shared__ GT piv[shmemsize/2 + 15];    // contains pivot column
   __shared__ GT shv[2][shmemsize/2 + 15]; // for the reduction

   v += blockIdx.z*workspace_size;
   R += blockIdx.z*workspace_size;

   int indR = r_pos(pivot, pivot+b, cols);
   // (dimR-1) - (pivot*(pivot+1))/2 - (b*(b+1))/2 - b*(pivot+1);

   piv[t_idx] = v[L];

   shv[0][t_idx] = v[i];
   shv[1][t_idx] = piv[t_idx].adj()*shv[0][t_idx];

   __syncthreads();

   rowsLog2 -= 1;
   int half_size = 1 << (rowsLog2);// sum for the norm
   if(t_idx + half_size < rows) 
   {
      shv[1][t_idx] = shv[1][t_idx] + shv[1][t_idx+half_size];
   }
   for(int k=0; k < rowsLog2; k++)
   {
      if(half_size > 16)
      {
         __syncthreads();
      }
      half_size /= 2;
      if(t_idx < half_size)
      {
         shv[1][t_idx] = shv[1][t_idx] + shv[1][t_idx+half_size];
      }
   }
   __syncthreads();

   shv[0][t_idx] = shv[0][t_idx] - shv[1][0]*piv[t_idx];
   v[i] = shv[0][t_idx];
   if(t_idx == 0) R[indR] = shv[1][0];
}

__global__ void mgs_small_backsubstitution_kernel1
 ( GT* R, GT* sol0, int dim, int workspace_size )
{
   int t_idx = threadIdx.x;
   __shared__ GT sol[32];
   __shared__ GT Rcl[32];
   int ind;
   GT update;

   R += blockIdx.x*workspace_size;
   sol0 += blockIdx.x*workspace_size;

   update = R[t_idx];
   for(int k=dim-1; k>=0; k--)  // compute k-th component of solution
   {
      if(t_idx < k+1)
      {
         ind = (dim - k)*(dim + 3 + k)/2 + t_idx;
         Rcl[t_idx] = R[ind];
      }
      if(t_idx == k) sol[t_idx] = update/Rcl[t_idx]; // all other threads wait
      // __syncthreads();
      if(t_idx < k) update = update - sol[k]*Rcl[t_idx]; // update
   }
   sol0[t_idx] = sol[t_idx];
}

__global__ void mgs_small_backsubstitution_kernel1_idx
 ( GT* R, GT* sol0, int dim, int workspace_size, int n_matrix_R,
   int* path_idx )
{
   int t_idx = threadIdx.x;
   __shared__ GT sol[32];
   __shared__ GT Rcl[32];
   int ind;
   GT update;

   R += blockIdx.x*n_matrix_R;
   sol0 += blockIdx.x*dim;

   update = R[t_idx];
   for(int k=dim-1; k>=0; k--)  // compute k-th component of solution
   {
      if(t_idx < k+1)
      {
         ind = (dim - k)*(dim + 3 + k)/2 + t_idx;
         Rcl[t_idx] = R[ind];
      }
      if(t_idx == k) sol[t_idx] = update/Rcl[t_idx]; // all other threads wait
      // __syncthreads();
      if(t_idx < k) update = update - sol[k]*Rcl[t_idx];// update
   }
   sol0[t_idx] = sol[t_idx];
}

__global__ void mgs_small_backsubstitution_kernel
 ( GT* R, GT* sol0, int dim, int workspace_size )
{
   int t_idx = threadIdx.x;
   __shared__ GT sol[shmemsize/2];
   __shared__ GT Rcl[shmemsize/2];
   int ind;
   GT update;

   R += blockIdx.x*workspace_size;
   sol0 += blockIdx.x*workspace_size;

   update = R[t_idx];
   for(int k=dim-1; k>=0; k--)  // compute k-th component of solution
   {
      if(t_idx < k+1)
      {
         ind = (dim - k)*(dim + 3 + k)/2 + t_idx;
         Rcl[t_idx] = R[ind];
      }
      if(t_idx == k) sol[t_idx] = update/Rcl[t_idx]; // all other threads wait
      __syncthreads();
      if(t_idx < k) update = update - sol[k]*Rcl[t_idx];// update
   }
   sol0[t_idx] = sol[t_idx];
}

__global__ void mgs_small_backsubstitution_max_kernel
 ( GT* R, GT* sol0, int dim, int dimLog2, double* max_delta_x )
{
   int t_idx = threadIdx.x;
   __shared__ GT sol[shmemsize/2];
   __shared__ GT Rcl[shmemsize/2];
   __shared__ double delta_x[shmemsize/2];
   int ind;
   GT update;

   update = R[t_idx];
   for(int k=dim-1; k>=0; k--)  // compute k-th component of solution
   {
      if(t_idx < k+1)
      {
         ind = (dim - k)*(dim + 3 + k)/2 + t_idx;
         Rcl[t_idx] = R[ind];
      }
      if(t_idx == k) sol[t_idx] = update/Rcl[t_idx]; // all other threads wait
      __syncthreads();
      if(t_idx < k) update = update - sol[k]*Rcl[t_idx];// update
   }
   sol0[t_idx] = sol[t_idx];

   // max for the norm
   delta_x[t_idx] = sol[t_idx].norm1_double();
   __syncthreads();

   dimLog2 -= 1;
   int half_size = 1 << (dimLog2);// sum for the norm
   if(t_idx + half_size < dim)
   {
      if(delta_x[t_idx] < delta_x[t_idx+half_size]) 
      {
         delta_x[t_idx] = delta_x[t_idx+half_size];
      }
   }
   for(int k=0; k < dimLog2; k++)
   {
      if(half_size > 16)
      {
         __syncthreads();
      }
      half_size /= 2;
      if(t_idx < half_size)
      {
         if(delta_x[t_idx] < delta_x[t_idx+half_size])
         {
            delta_x[t_idx] = delta_x[t_idx+half_size];
         }
      }
   }
   if(t_idx == 0)
   {
      *max_delta_x = delta_x[0];
   }
}
