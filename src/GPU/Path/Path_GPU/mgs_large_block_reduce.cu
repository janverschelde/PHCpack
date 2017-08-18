#ifndef __PATH_GPU_MGS_LARGE_BLOCK_REDUCE_CU__
#define __PATH_GPU_MGS_LARGE_BLOCK_REDUCE_CU__

__global__ void matrix_init ( GT* V, int rows, int cols )
{
   int tidx = threadIdx.x;
   int bidx = (gridDim.x*blockIdx.y+blockIdx.x)*blockDim.x;
   int idx = tidx + bidx;
   if(idx<rows*cols)
   {
      V[idx] = GT(1,0);
      // V[idx] = GT(idx,0);
   }
}

__global__ void matrix_init_zero ( GT* V, int rows, int cols )
{
   int tidx = threadIdx.x;
   int bidx = (gridDim.x*blockIdx.y+blockIdx.x)*blockDim.x;
   int idx = tidx + bidx;
   if(idx<rows*cols)
   {
      V[idx] = GT(0,0);
   }
}

template <unsigned int row_block_size>
__global__ void mgs_large_reduce_block_kernel1_single_pivot
 ( GT* V, GT* P, int row_block, int col_block, int rows, int cols,
   int piv_col_block )
/*
 * Constant: matrix_block_row, matrix_block_pivot_col, matrix_block_reduce_col
 * row_block_size = matrix_block_row
 */
{
   int b = gridDim.x*blockIdx.y+blockIdx.x;
   int tidx = threadIdx.x;

   __shared__ GT pivot_block[matrix_block_row];    // contains pivot column
   __shared__ GT product_array[matrix_block_row];  // contains pivot column

   int row_block_idx = b-b/row_block*row_block;
   // b-(col_block_idx-1)*row_block
   int reduce_start_col
      = b/row_block*matrix_block_reduce_col + piv_col_block + 1;

   GT* pivot_block_start
      = V + piv_col_block*rows + row_block_idx*row_block_size;
   GT* reduce_block_start
      = V + reduce_start_col*rows + row_block_idx*row_block_size;
   int kernel_rows = rows - row_block_idx*row_block_size;

   // load pivot block and reduce block
   if(tidx < kernel_rows)
   {
      pivot_block[tidx] = pivot_block_start[tidx];
   }
   int reduce_block_size = min(cols-reduce_start_col, matrix_block_reduce_col);

   // Compute product
   for(int reduce_idx=0; reduce_idx<reduce_block_size; reduce_idx++)
   {
      if(tidx < kernel_rows)
      {
         product_array[tidx] = reduce_block_start[reduce_idx*rows+tidx];
         // product_array[tidx] = pivot_block[piv_idx*row_block_size
         //  +tidx].adj()*reduce_block[tidx];
         product_array[tidx]
            = pivot_block[tidx].adj_multiple(product_array[tidx]);
      }
      if(row_block_size > 64)
      {
         __syncthreads();
         if(tidx < 64 && tidx+64 < kernel_rows)
         {
            product_array[tidx] = product_array[tidx] + product_array[tidx+64];
         }
      }
      if(row_block_size > 32)
      {
         __syncthreads();
         if(tidx < 32 && tidx+32 < kernel_rows)
         {
            product_array[tidx] = product_array[tidx] + product_array[tidx+32];
         }
      }
      if(row_block_size > 16)
      {
         if(tidx < 16 && tidx+16 < kernel_rows)
         {
            product_array[tidx] = product_array[tidx] + product_array[tidx+16];
         }
      }
      if(row_block_size > 8)
      {
         if(tidx < 8 && tidx+8 < kernel_rows)
         {
            product_array[tidx] = product_array[tidx] + product_array[tidx+8];
         }
      }
      if(row_block_size > 4)
      {
         if(tidx < 4 && tidx+4 < kernel_rows)
         {
            product_array[tidx] = product_array[tidx] + product_array[tidx+4];
         }
      }
      if(row_block_size > 2)
      {
         if(tidx < 2 && tidx+2 < kernel_rows)
         {
            product_array[tidx] = product_array[tidx] + product_array[tidx+2];
         }
      }
      if(row_block_size > 1)
      {
         if(tidx < 1 && tidx+1 < kernel_rows)
         {
            product_array[0] = product_array[0] + product_array[1];
         }
      }
      if(tidx == 0)
      {
         int col_idx = reduce_start_col+reduce_idx;
         P[row_block_idx*cols+col_idx] = product_array[0];
      }
      if(row_block_size > 32)
      {
         __syncthreads();
      }
   }
}

template <unsigned int row_block_size>
__global__ void mgs_large_reduce_block_kernel1
 ( GT* V, GT* P, int row_block, int col_block, int rows, int cols,
   int piv_col_block )
/*
 * Constant: matrix_block_row, matrix_block_pivot_col, matrix_block_reduce_col
 * row_block_size = matrix_block_row
 */
{
   int b = gridDim.x*blockIdx.y+blockIdx.x;
   int tidx = threadIdx.x;

   __shared__ GT pivot_block[matrix_block_row*matrix_block_pivot_col];
   // contains pivot column
   __shared__ GT reduce_block[matrix_block_row];  // contains pivot column
   __shared__ GT product_array[matrix_block_row];  // contains pivot column

   int row_block_idx = b-b/row_block*row_block;
   // b-(col_block_idx-1)*row_block
   int reduce_start_col = b/row_block*matrix_block_reduce_col
      + (piv_col_block + 1)*matrix_block_pivot_col;

   GT* pivot_block_start = V + piv_col_block*matrix_block_pivot_col*rows
      + row_block_idx*row_block_size;
   GT* reduce_block_start = V + reduce_start_col*rows
      + row_block_idx*row_block_size;
   int kernel_rows = rows - row_block_idx*row_block_size;

   // load pivot block and reduce block
   for(int i=0; i<matrix_block_pivot_col; i++)
   {
      if(tidx < kernel_rows)
      {
         pivot_block[i*row_block_size+tidx] = pivot_block_start[i*rows+tidx];
      }
   }
   int reduce_block_size = min(cols-reduce_start_col, matrix_block_reduce_col);

   // Compute product
   for(int reduce_idx=0; reduce_idx<reduce_block_size; reduce_idx++)
   {
      if(tidx < kernel_rows)
      {
         reduce_block[tidx] = reduce_block_start[reduce_idx*rows+tidx];
      }
      for(int piv_idx=0; piv_idx<matrix_block_pivot_col; piv_idx++)
      {
         if(tidx<kernel_rows)
         {
            // product_array[tidx] = pivot_block[piv_idx*row_block_size
            //    +tidx].adj()*reduce_block[tidx];
            product_array[tidx] = pivot_block[piv_idx*row_block_size
               +tidx].adj_multiple(reduce_block[tidx]);
         }
         if(row_block_size > 128)
         {
            __syncthreads();
            if(tidx < 128 && tidx+128 < kernel_rows)
            {
               product_array[tidx] = product_array[tidx]
                  + product_array[tidx+128];
            }
         }
         if(row_block_size > 64)
         {
            __syncthreads();
            if(tidx < 64 && tidx+64 < kernel_rows)
            {
               product_array[tidx] = product_array[tidx]
                  + product_array[tidx+64];
            }
         }
         if(row_block_size > 32)
         {
            __syncthreads();
            if(tidx < 32 && tidx+32 < kernel_rows)
            {
               product_array[tidx] = product_array[tidx]
                  + product_array[tidx+32];
            }
         }
         if(row_block_size > 16)
         {
            if(tidx < 16 && tidx+16 < kernel_rows)
            {
               product_array[tidx] = product_array[tidx]
                  + product_array[tidx+16];
            }
         }
         if(row_block_size > 8)
         {
            if(tidx < 8 && tidx+8 < kernel_rows)
            {
               product_array[tidx] = product_array[tidx]
                  + product_array[tidx+8];
            }
         }
         if(row_block_size > 4)
         {
            if(tidx < 4 && tidx+4 < kernel_rows)
            {
               product_array[tidx] = product_array[tidx]
                  + product_array[tidx+4];
            }
         }
         if(row_block_size > 2)
         {
            if(tidx < 2 && tidx+2 < kernel_rows)
            {
               product_array[tidx] = product_array[tidx]
                  + product_array[tidx+2];
            }
         }
         if(row_block_size > 1)
         {
            if(tidx < 1 && tidx+1 < kernel_rows)
            {
               product_array[0] = product_array[0] + product_array[1];
            }
         }
         if(tidx == 0)
         {
            int col_idx = reduce_start_col+reduce_idx;
            P[row_block_idx*cols*matrix_block_pivot_col
               + col_idx*matrix_block_pivot_col+piv_idx] = product_array[0];
         }
         if(row_block_size > 32)
         {
            __syncthreads();
         }
      }
   }
}

//template <unsigned int n_th>
__global__ void mgs_large_reduce_block_kernel3_single_pivot
 ( GT* V, GT* R, int row_block, int col_block, int rows, int cols,
   int piv_col_block )
/*
 * Constant: matrix_block_row, matrix_block_pivot_col, matrix_block_reduce_col
 */
{
   int b = gridDim.x*blockIdx.y+blockIdx.x;
   int tidx = threadIdx.x;

   __shared__ GT pivot_block[matrix_block_row*matrix_block_pivot_col];
   // contains pivot column
   __shared__ GT reduce_block[matrix_block_row];  // contains reduce column
   // __shared__ GT norm_array[matrix_block_pivot_col];
   // contains pivot column

   // int col_block_idx = b/row_block + 1 + piv_col_block;
   int row_block_idx = b-b/row_block*row_block;
   // b-(col_block_idx-1)*row_block
   int reduce_start_col = b/row_block*matrix_block_reduce_col
      + (piv_col_block + 1);

   GT* pivot_block_start = V + piv_col_block*rows
      + row_block_idx*matrix_block_row;
   GT* reduce_block_start = V + reduce_start_col*rows
      + row_block_idx*matrix_block_row;
   int kernel_rows = rows - row_block_idx*matrix_block_row;

   int reduce_block_size = min(cols-reduce_start_col, matrix_block_reduce_col);

   if(tidx < kernel_rows)
   {
      // load pivot block and reduce block
      pivot_block[tidx] = pivot_block_start[tidx];

      for(int reduce_idx=0; reduce_idx<reduce_block_size; reduce_idx++)
      {
         reduce_block[tidx] = reduce_block_start[reduce_idx*rows+tidx];
         int start_r_pos = r_pos(piv_col_block*matrix_block_pivot_col,
                                 reduce_start_col+reduce_idx, cols);
         reduce_block_start[reduce_idx*rows+tidx]
            = reduce_block[tidx] - R[start_r_pos]*pivot_block[tidx];
      }
   }
}

//template <unsigned int n_th>
__global__ void mgs_large_reduce_block_kernel3
 ( GT* V, GT* R, int row_block, int col_block, int rows, int cols,
   int piv_col_block )
/*
 * Constant: matrix_block_row, matrix_block_pivot_col, matrix_block_reduce_col
 */
{
   int b = gridDim.x*blockIdx.y+blockIdx.x;
   int tidx = threadIdx.x;

   __shared__ GT pivot_block[matrix_block_row*matrix_block_pivot_col];
   // contains pivot column
   __shared__ GT reduce_block[matrix_block_row];  // contains reduce column
   __shared__ GT norm_array[matrix_block_pivot_col];  // contains pivot column

   // int col_block_idx = b/row_block + 1 + piv_col_block;
   int row_block_idx = b-b/row_block*row_block;
   // b-(col_block_idx-1)*row_block
   int reduce_start_col = b/row_block*matrix_block_reduce_col
      + (piv_col_block + 1)*matrix_block_pivot_col;

   GT* pivot_block_start = V + piv_col_block*matrix_block_pivot_col*rows
      + row_block_idx*matrix_block_row;
   GT* reduce_block_start = V + reduce_start_col*rows
      + row_block_idx*matrix_block_row;
   int kernel_rows = rows - row_block_idx*matrix_block_row;

   int reduce_block_size = min(cols-reduce_start_col, matrix_block_reduce_col);

   if(tidx < kernel_rows)
   {
      // load pivot block and reduce block
      for(int i=0; i<matrix_block_pivot_col; i++)
      {
         pivot_block[i*matrix_block_row+tidx] = pivot_block_start[i*rows+tidx];
      }
      for(int reduce_idx=0; reduce_idx<reduce_block_size; reduce_idx++)
      {
         reduce_block[tidx] = reduce_block_start[reduce_idx*rows+tidx];
         GT tmp = reduce_block[tidx];
         for(int piv_idx=0; piv_idx<matrix_block_pivot_col; piv_idx++)
         {
            //]tmp = tmp
            // - Norm[(b/row_block*matrix_block_reduce_col
            // +reduce_idx)*matrix_block_pivot_col
            // +piv_idx]*pivot_block[piv_idx*matrix_block_row+tidx];

            int start_r_pos = r_pos(piv_col_block*matrix_block_pivot_col
               +piv_idx,reduce_start_col+reduce_idx, cols);
            tmp = tmp
               - R[start_r_pos]*pivot_block[piv_idx*matrix_block_row+tidx];
         }
         reduce_block_start[reduce_idx*rows+tidx] = tmp;
      }
   }
}

__global__ void mgs_large_reduce_block_kernel2
 ( GT* P, GT* R, int row_block, int P_rows, int n_sum, 
   int col_block_size, int cols, int piv_block )
{
   int b = gridDim.x*blockIdx.y+blockIdx.x;
   int tidx = threadIdx.x;
   int idx = b*blockDim.x+tidx;

   /*
     int block = b+pivot; // column for reduction w.r.t. pivot
     int vBSind = 0;
     int powTwo;
     int indR = (dimR-1) - (pivot*(pivot+1))/2 - (b*(b+1))/2 - b*(pivot+1);
    */

   // Compute product
   if(idx<n_sum)
   {
      GT sum = P[idx];
      for(int i=1; i< row_block; i++)
      {
         sum = sum + P[idx+i*P_rows];
      }
      int col_idx = idx-idx/col_block_size*col_block_size;
      int r_x_idx = piv_block*col_block_size;
      int r_y_idx = r_x_idx + col_block_size + idx/col_block_size;
      // r_x_idx += col_idx;

      int start_r_pos = r_pos(r_x_idx,r_y_idx, cols);

      // Norm[idx] = sum;//GT(idx, 0);
      R[start_r_pos+col_idx] = sum;
   }
}

void mgs_large_reduce_block_rest
 ( GT* V, GT* R, GT* P, int rows, int cols, int pivot_block )
{
   int BS = matrix_block_pivot_col;

   // std::cout << "rows   = " << rows << std::endl;
   // std::cout << "cols   = " << cols << std::endl;
   // std::cout << "BS     = " << BS << std::endl;

   int row_block = (rows-1)/matrix_block_row+1;
   int col_block = (cols-1)/matrix_block_pivot_col+1;
   // last column doesn't count

   // std::cout << "row_block = " << row_block << std::endl;
   // std::cout << "col_block = " << col_block << std::endl;

   int NB = (rows*cols-1)/BS+1;
   // matrix_init<<<NB, BS>>>(V, rows, cols);

   // NB = (row_block*matrix_block_pivot_col*cols-1)/BS+1;
   // matrix_init_zero<<<NB, BS>>>(P, row_block, matrix_block_pivot_col*cols);

   // GT* P_host
   //    = (GT*)malloc(row_block*matrix_block_pivot_col*cols*sizeof(GT));
   // matrix_init<<<NB, BS>>>(P, col_block*row_block, cols-1);

   NB = row_block*((cols-(pivot_block+1)*matrix_block_pivot_col
                        -1)/matrix_block_reduce_col+1);
   // std::cout << "NB = " << NB << std::endl;
   dim3 NB3 = get_grid(NB,1);
   if(matrix_block_pivot_col==1)
   {
      mgs_large_reduce_block_kernel1_single_pivot
         <matrix_block_row><<<NB3,matrix_block_row>>>
         (V, P, row_block, col_block, rows, cols, pivot_block);
   }
   else
   {
      mgs_large_reduce_block_kernel1
         <matrix_block_row><<<NB3,matrix_block_row>>>
         (V, P, row_block, col_block, rows, cols, pivot_block);
   }

   /*
     cudaMemcpy(P_host, P, row_block*matrix_block_pivot_col*cols*sizeof(GT),
                cudaMemcpyDeviceToHost);

     std::cout << "------------- Matrix P -----------" << std::endl;
     for(int col_idx=0; col_idx<cols; col_idx++)
     {
        for(int block_col_idx=0; block_col_idx<matrix_block_pivot_col;
                block_col_idx++)
        {
           for(int row_block_idx=0; row_block_idx<row_block; row_block_idx++)
           {
              GT tmp_matrix = P_host[row_block_idx*cols*matrix_block_pivot_col
                 +col_idx*matrix_block_pivot_col+block_col_idx];
              std::cout << "col=" << col_idx << " col_block=" << block_col_idx
                        << " row_block=" << row_block_idx
                        << "   " << tmp_matrix.real << " + "
                        << tmp_matrix.imag << std::endl;
           }
        }
     }
    */

   NB = (matrix_block_pivot_col*cols-1)/BS+1;
   // GT* Norm;
   // cudaMalloc((void**)&Norm, matrix_block_pivot_col*cols*sizeof(GT));
   //matrix_init_zero<<<NB, BS>>>(Norm, matrix_block_pivot_col, cols);

   int BS_sum = 32;
   int P_rows = cols*matrix_block_pivot_col;
   int n_sum
      = (cols-(pivot_block+1)*matrix_block_pivot_col)*matrix_block_pivot_col;
   NB = (n_sum-1)/BS_sum+1;
   // std::cout << "n_sum = " << n_sum
   // << " start = " << (cols-1)*col_block << std::endl;
   NB3 = get_grid(NB,1);
   mgs_large_reduce_block_kernel2<<<NB3,BS_sum>>>
      (P+(pivot_block+1)*matrix_block_pivot_col*matrix_block_pivot_col, R,
       row_block, P_rows, n_sum, matrix_block_pivot_col, cols, pivot_block);

   /*
     std::cout << "------------- Norm -----------" << std::endl;
     GT* Norm_host = (GT*)malloc(matrix_block_pivot_col*cols*sizeof(GT));

     cudaMemcpy(Norm_host, Norm, matrix_block_pivot_col*cols*sizeof(GT),
                cudaMemcpyDeviceToHost);

     for(int i=0; i <cols*matrix_block_pivot_col; i++)
     {
        std::cout << i << " " << Norm_host[i].real << " + "
                  << Norm_host[i].imag << std::endl;
     }
    */

   NB = row_block*((cols-(pivot_block+1)*matrix_block_pivot_col
                        -1)/matrix_block_reduce_col+1);
   // std::cout << "NB = " << NB << std::endl;

   NB3 = get_grid(NB,1);
   if(matrix_block_pivot_col==1)
   {
      mgs_large_reduce_block_kernel3_single_pivot<<<NB3,matrix_block_row>>>
         (V, R, row_block, col_block, rows, cols, pivot_block);
   }
   else
   {
      mgs_large_reduce_block_kernel3<<<NB3,matrix_block_row>>>
         (V, R, row_block, col_block, rows, cols, pivot_block);
   }
}

#endif /*__PATH_GPU_MGS_LARGE_BLOCK_REDUCE_CU__*/
