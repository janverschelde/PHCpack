/* The file smDMMA_kernels.cu contains the definitions of the functions with
 * prototypes in smDMMA_kernels.h. */

#include <iostream>
#include <assert.h>
#include "smDMMA_dims.h"
#include "mma.h"

// copied from helper_cuda.h

#ifndef MAX
#define MAX(a, b) (a > b ? a : b)
#endif

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

__global__ void compute_dgemm
 ( const double *A, const double *B, const double *C, double *D,
   double alpha, double beta)
{
#if __CUDA_ARCH__ >= 800
   extern __shared__ double shmem[][CHUNK_K * K + SKEW_DOUBLE];

   // Warp and lane identification.
   const unsigned int warpId = threadIdx.x / WARP_SIZE;
   const unsigned int laneId = threadIdx.x % WARP_SIZE;

   // Offset in shared memory from which the B matrix is stored.
   const size_t shmem_idx_b_off = BLOCK_COL_TILES * M;

   // This pointer is used to access the C and D matrix tiles 
   // this warp computes.
   double *shmem_warp_tile_ptr = (double*)&shmem[0][0]
      + (warpId / BLOCK_ROW_WARPS) * SHMEM_STRIDE * N * BLOCK_ROW_WARPS
      + (warpId % BLOCK_ROW_WARPS) * SHMEM_OFFSET;

   // This pointer is used to stream the C and D matrices block-wide tile
   // to and from shared memory.
   double *shmem_warp_stream_ptr = (double*)&shmem[0][0]
      + warpId * SHMEM_STRIDE * N;

   // Adjust the beta scaler, as it'll be multiplied by alpha at the end of
   // each tile computation. Technically this is not generally correct
   // (may result in a loss of precision).
   // Zero still needs to be specially handled though.
   beta /= alpha;

   // Each CTA slides along the 64 x 64 tiles from the top left corner
   // of the matrix to the right and down, and selects the next tile 
   // to compute. Once there's no such tile, all warps in this CTA exit.
   for(unsigned int block_pos = blockIdx.x;; block_pos += gridDim.x)
   {
      const unsigned int block_tile_i
      = ((block_pos * BLOCK_ROW_TILES) / N_TILES) * (BLOCK_COL_TILES);
      const unsigned int block_tile_j
      = (block_pos * BLOCK_COL_TILES) % N_TILES;

      // Stop when there are no more D matrix tiles to compute in this CTA.
      if(block_tile_i >= M_TILES) break;

      // This warp's pointer to the C matrix data to copy memory 
      // from to shared memory.
      const size_t gmem_idx
      = (block_tile_i + warpId) * M * GLOBAL_MEM_STRIDE + block_tile_j * N;
      const double *src_gmem_warp_stream_ptr = &C[gmem_idx];

      // Stream multiple C tiles to shared memory.
#pragma unroll
      for(int i = 0; i < N; i++)
      {
         *((int4 *)(shmem_warp_stream_ptr + SHMEM_STRIDE * i) + laneId)
         = *((int4 *)(src_gmem_warp_stream_ptr + GLOBAL_MEM_STRIDE * i)
           + laneId);
      }
      __syncthreads();

      // These fragments will accumulate the result of A and B matrix
      // fragment multiplications along the K_GLOBAL dimension.
      wmma::fragment<wmma::accumulator, M, N, K, double>
         c[WARP_COL_TILES][WARP_ROW_TILES];

      // Load the C matrix tiles into fragments from shared memory.
#pragma unroll
      for(int i = 0; i < WARP_COL_TILES; i++)
      {
#pragma unroll
         for(int j = 0; j < WARP_ROW_TILES; j++)
         {
            const double *tile_ptr
            = shmem_warp_tile_ptr + i * SHMEM_STRIDE * N + j * N;

            wmma::load_matrix_sync(c[i][j],tile_ptr,SHMEM_STRIDE,C_LAYOUT);
         }
      }
      __syncthreads();

      // Scale the C matrix.
#pragma unroll
      for(int i = 0; i < WARP_COL_TILES; i++)
      {
#pragma unroll
         for(int j = 0; j < WARP_ROW_TILES; j++)
         {
#pragma unroll
            for(int t = 0; t < c[i][j].num_elements; t++)
            {
               c[i][j].x[t] *= beta;
            }
         }
      }
      // Select what warp copies what matrix to shared memory.
      // Warps 0-3 copy the A matrix, warps 4-7 copy the B matrix.
      const double *warp_ptr
         = (warpId < (WARPS_PER_BLOCK/2)) ? (&A[block_tile_i * M * K_GLOBAL]
           + M * K_GLOBAL * (warpId % (WARPS_PER_BLOCK/2)) * 2) :
           (&B[block_tile_j * N * K_GLOBAL]
             + N * K_GLOBAL * (warpId % (WARPS_PER_BLOCK/2)) * 2);

      // Go through the global K dimension by a fixed step at a time.
#pragma unroll
      for(int tile_k = 0; tile_k < K_TILES; tile_k += CHUNK_K)
      {
         // Copy slices of the A and B matrices to shared memory.
         // The first half of the warps in the CTA copy the A matrix,
         // the rest copy the B matrix.
         size_t shmem_idx = warpId < (WARPS_PER_BLOCK/2) ?
            (M * (warpId % (WARPS_PER_BLOCK/2)) * 2) :
            (N * (warpId % (WARPS_PER_BLOCK/2)) * 2 + shmem_idx_b_off);

         // First half of the warp copies the first row/column of the matrix,
         // the second half of the warp copies the next.
         const double *lane_ptr = warp_ptr + tile_k * K
            + (laneId / CHUNK_COPY_LINE_LANES) * K_GLOBAL;

         // Shift the second half of the warp to the next row/column
         // in the shared memory.
         shmem_idx += laneId / CHUNK_COPY_LINE_LANES;

#pragma unroll
         for(int i = 0; i < ((WARP_SIZE/2) / CHUNK_COPY_LINES_PER_WARP); i++)
         {
            // Copy 16 bytes at once in each lane.
            *((int4*)&shmem[shmem_idx][0] + (laneId % CHUNK_COPY_LINE_LANES))
            = *((int4*)lane_ptr +  (laneId % CHUNK_COPY_LINE_LANES));

            // Advance the global memory pointer and the shared memory index.
            lane_ptr = lane_ptr + K_GLOBAL * CHUNK_COPY_LINES_PER_WARP;
            shmem_idx += CHUNK_COPY_LINES_PER_WARP;
         }
         __syncthreads();

         // Compute a grid of C matrix tiles in each warp.
#pragma unroll
         for(int k_step = 0; k_step < CHUNK_K; k_step++)
         {
            wmma::fragment<wmma::matrix_a, M, N, K, double, wmma::row_major>
               a[WARP_COL_TILES];
            wmma::fragment<wmma::matrix_b, M, N, K, double, wmma::col_major>
               b[WARP_ROW_TILES];

#pragma unroll
            for(int i = 0; i < WARP_COL_TILES; i++)
            {
               size_t shmem_idx_a = (warpId/2) * M * 2 + (i * M);
               const double *tile_ptr = &shmem[shmem_idx_a][k_step * K];

               wmma::load_matrix_sync(a[i],tile_ptr,K*CHUNK_K + SKEW_DOUBLE);

#pragma unroll
               for(int j = 0; j < WARP_ROW_TILES; j++)
               {
                  if(i == 0)
                  {
                     // Load the B matrix fragment once,
                     // because it is going to be reused
                     // against the other A matrix fragments.
                     size_t shmem_idx_b = shmem_idx_b_off
                        + (WARP_ROW_TILES * N) * (warpId%2) + (j * N);
                     const double *tile_ptr = &shmem[shmem_idx_b][k_step * K];

                     wmma::load_matrix_sync(b[j], tile_ptr,
                        K * CHUNK_K + SKEW_DOUBLE);
                  }
                  wmma::mma_sync(c[i][j], a[i], b[j], c[i][j]);
               }
            }
         }
         __syncthreads();
      }
      // Store the D fragments to shared memory.
#pragma unroll
      for(int i = 0; i < WARP_COL_TILES; i++)
      {
#pragma unroll
         for(int j = 0; j < WARP_ROW_TILES; j++)
         {
            // Uniform, point-wise transformations of ALL fragment elements
            // by ALL threads in the warp are well-defined even though element
            // indices within fragment storage are not defined.
#pragma unroll
            for(int t = 0; t < c[i][j].num_elements; t++)
               c[i][j].x[t] *= alpha;

            double *tile_ptr = shmem_warp_tile_ptr + i*SHMEM_STRIDE*N + j*N;

            wmma::store_matrix_sync(tile_ptr, c[i][j], SHMEM_STRIDE, C_LAYOUT);
         }
      }
      __syncthreads();

      // Now that shared memory contains all the D tiles,
      // stream them to global memory.
      double *dst_gmem_warp_stream_ptr = &D[gmem_idx];

#pragma unroll
      for(int i = 0; i < N; i++)
      {
         *((int4*)(dst_gmem_warp_stream_ptr + GLOBAL_MEM_STRIDE * i) + laneId)
	 = *((int4*)(shmem_warp_stream_ptr + SHMEM_STRIDE * i) + laneId);
      }
      __syncthreads();
   }
#endif
}

__host__ void GPU_DMMA
 ( double *A, double *B, double *C, double *D, double alpha, double beta,
   cudaDeviceProp deviceProp, int kerneltype )
{
   dim3 gridDim;
   dim3 blockDim;

   // blockDim.x must be a multple of warpSize 128x4 means we have
   // 16 warps and a block computes a 64x64 output tile
   blockDim.x = 128;
   blockDim.y = 4;

   gridDim.x = (M_GLOBAL + (M*blockDim.x/32-1))/(M*blockDim.x/32);
   gridDim.y = (N_GLOBAL + N*blockDim.y-1)/(N*blockDim.y);

   if(kerneltype == 0)
   {
      cout << "Computing ... using simple_wmma_gemm kernel ... ";
      simple_wmma_gemm<<<gridDim, blockDim>>>
         (A, B, C, D, M_GLOBAL, N_GLOBAL, K_GLOBAL, alpha, beta);
      cout << "done." << endl;
   }
   else
   {
      enum
      {
         // Compute the right amount of shared memory to request.
         // We need shared memory to hold per-CTA C and D matrix tiles,
         // and to cache per-CTA chunks of the A and B matrices.
         // Therefore, the right amount to request is the maximum
         // of those two numbers.
         SHMEM_SZ = MAX(sizeof(double) * (BLOCK_COL_TILES * M)
                        * (CHUNK_K * K + SKEW_DOUBLE) * 2,
                        M * (BLOCK_ROW_WARPS * WARP_ROW_TILES)
                          * N * (BLOCK_COL_WARPS * WARP_COL_TILES)
                          * sizeof(double))
      };

      printf("Required shared memory size: %lu Kb\n", SHMEM_SZ / 1024UL);

      cout << "Computing ... using compute_dgemm kernel ... ";

      assert(cudaFuncSetAttribute(compute_dgemm,
             cudaFuncAttributeMaxDynamicSharedMemorySize, SHMEM_SZ) == 0);
      compute_dgemm<<<deviceProp.multiProcessorCount*2, THREADS_PER_BLOCK,
                      SHMEM_SZ>>>(A, B, C, D, alpha, beta);

      cout << "done." << endl;
   }
}
