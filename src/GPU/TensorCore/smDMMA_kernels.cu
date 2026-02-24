/* The file smDMMA_kernels.cu contains the definitions of the functions with
 * prototypes in smDMMA_kernels.h. */

#include <iostream>
#include "smDMMA_dims.h"
#include "mma.h"

using namespace std;
using namespace nvcuda;

__global__ void simple_wmma_gemm
 ( double *a, double *b, double *c, double *d, int m_ld, int n_ld, int k_ld,
   double alpha, double beta )
/*
   Performs an MxNxK DGEMM (C=alpha*A*B + beta*C) assuming:
   1) Matrices are packed in memory.
   2) M, N and K are multiples of 8, 8 and 4 respectively. 
   3) A is row major, B is column major matrix.
   Note: This is a less performant version of the compute_dgemm kernel.
   It is designed for demonstration purposes only to show the CUDA WMMA API 
   use without relying on availability of the shared memory.  */
{
#if __CUDA_ARCH__ >= 800
    // Leading dimensions. Packed with no transpositions.
   int lda = k_ld;
   int ldb = k_ld;
   int ldc = n_ld;

   // tile using a 2D grid
   int warpM = (blockIdx.x * blockDim.x + threadIdx.x) / warpSize;
   int warpN = (blockIdx.y * blockDim.y + threadIdx.y);

   // declare the fragments
   wmma::fragment<wmma::matrix_a, M, N, K, double, wmma::row_major> a_frag;
   wmma::fragment<wmma::matrix_b, M, N, K, double, wmma::col_major> b_frag;
   wmma::fragment<wmma::accumulator, M, N, K, double> acc_frag;
   wmma::fragment<wmma::accumulator, M, N, K, double> c_frag;

   wmma::fill_fragment(acc_frag, 0.0f);

   for(int i = 0; i < k_ld; i += K) // loop over k
   {
      int aCol = i;
      int aRow = warpM * M;
      int bCol = warpN * N;
      int bRow = i;

      // bounds checking
      if(aRow < m_ld && aCol < k_ld && bRow < k_ld && bCol < n_ld)
      {
         // load the inputs
         wmma::load_matrix_sync(a_frag, a + aCol + aRow * lda, lda);
         wmma::load_matrix_sync(b_frag, b + bRow + bCol * ldb, ldb);
         // perform the matrix multiplication
         wmma::mma_sync(acc_frag, a_frag, b_frag, acc_frag);
      }
   }
   // load in the current value of c, scale it by beta,
   // and add this our result scaled by alpha
   int cCol = warpN * N;
   int cRow = warpM * M;

   if(cRow < m_ld && cCol < n_ld)
   {
      wmma::load_matrix_sync(c_frag,c+cCol+cRow*ldc,ldc,wmma::mem_row_major);

      for(int i=0; i < c_frag.num_elements; i++)
      {
         c_frag.x[i] = alpha * acc_frag.x[i] + beta * c_frag.x[i];
      }
      // Store the output
      wmma::store_matrix_sync(d+cCol+cRow*ldc,c_frag,ldc,wmma::mem_row_major);
   }
#endif
}

__host__ void GPU_DMMA
 ( double *A, double *B, double *C, double *D, double alpha, double beta )
{
   dim3 gridDim;
   dim3 blockDim;

   // blockDim.x must be a multple of warpSize 128x4 means we have
   // 16 warps and a block computes a 64x64 output tile
   blockDim.x = 128;
   blockDim.y = 4;

   gridDim.x = (M_GLOBAL + (M*blockDim.x/32-1))/(M*blockDim.x/32);
   gridDim.y = (N_GLOBAL + N*blockDim.y-1)/(N*blockDim.y);

   cout << "Computing ... using simple_wmma_gemm kernel ... ";
   simple_wmma_gemm<<<gridDim, blockDim>>>
      (A, B, C, D, M_GLOBAL, N_GLOBAL, K_GLOBAL, alpha, beta);
   cout << "done." << endl;
}
