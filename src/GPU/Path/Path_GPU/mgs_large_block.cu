#include "mgs_large_backsub.cu"
#include "mgs_large_norm.cu"
#include "mgs_large_row_reduce.cu"
#include "mgs_large_block_reduce.cu"

void mgs_large_reduce_block
 ( GT* V, GT* R, GT* P, GT* sol, int rows, int cols );

void mgs_large_block
 ( GT* V, GT* R, GT* P, GT* sol, int rows, int cols )
{
   mgs_large_reduce_block(V, R, P, sol, rows, cols);
   mgs_large_backsubstitution(R, sol, rows, cols);
}

void mgs_large_reduce_block
 ( GT* V, GT* R, GT* P, GT* sol, int rows, int cols)
{
   /*std::cout << "------------- Matrix P Seq-----------" << std::endl;
   GT* tmp_matrix = P_host;
   int tmp_idx = 0;
   for(int col_idx=0; col_idx<cols-1; col_idx++)
   {
      for(int col_block_idx=0; col_block_idx<col_block; col_block_idx++)
      {
         for(int row_block_idx=0; row_block_idx<row_block; row_block_idx++)
         {
            std::cout << tmp_idx++ << " "
                      << (*tmp_matrix).real << " + "
                      << (*tmp_matrix).imag << std::endl;
            tmp_matrix++;
         }
      }
   }*/

   int BS = min(BS_QR,rows);
   // int BS = 32;

   int rowsLog2 = log2ceil(rows);// ceil for sum reduction
   // int dimR = cols*(cols+1)/2;

   T* pivnrm; // norm of the pivot column
   cudaMalloc((void**)&pivnrm,sizeof(T));

   T* sums_global; // norm of the pivot column
   cudaMalloc((void**)&sums_global,maxrounds*sizeof(T));

   int rf = (rows-1)/BS+1;
   int rfLog2 = log2ceil(rf);
   int BSLog2 = log2ceil(BS);
   int lastBSLog2 = log2ceil(rows-BS*(rf-1));

   int col_block = (cols-1)/matrix_block_pivot_col + 1;

   // int NB = (rows*cols-1)/BS+1;
   // matrix_init<<<NB, BS>>>(V, rows, cols);

   // std::cout << "rf = " << rf << std::endl;
   // std::cout << "col_block = " << col_block << std::endl;

   for(int col_block_idx=0; col_block_idx<col_block; col_block_idx++)
   {
      int piv_start = col_block_idx*matrix_block_pivot_col;
      int piv_end;
      if(col_block_idx != col_block-1)
      {
         piv_end = (col_block_idx+1)*matrix_block_pivot_col;
      }
      else
      {
         piv_end = cols-1;
      }

      // std::cout << "piv_start = " << piv_start
      //           << ", piv_end = " << piv_end << std::endl;

      for(int piv=piv_start; piv<piv_end; piv++) 
      {
         mgs_large_normalize_kernel1<<<rf,BS>>>
            (V,R,rows,rowsLog2,cols,piv,rf,rfLog2,BS,BSLog2,
             pivnrm,lastBSLog2,sums_global);
         mgs_large_normalize_kernel2<<<rf,BS>>>
            (V,R,rows,rowsLog2,cols,piv,rf,rfLog2,BS,BSLog2,
             pivnrm,lastBSLog2,sums_global);
         // XXX BS should be greater than maxround
         mgs_large_row_reduce_kernel<<<piv_end-piv-1,BS>>>
            (V,R,cols,rows,rowsLog2, piv,rf,rfLog2,BS,BSLog2,
             pivnrm,lastBSLog2);
      }
      if(col_block_idx != col_block-1)
      {
         mgs_large_reduce_block_rest(V, R, P, rows, cols, col_block_idx);
      }
   }
}
