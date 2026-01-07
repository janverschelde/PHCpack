/* Copyright (c) 2022, NVIDIA CORPORATION. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *  * Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *  * Neither the name of NVIDIA CORPORATION nor the names of its
 *    contributors may be used to endorse or promote products derived
 *    from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ``AS IS'' AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
 * OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

// Formatted version, mainly to fit 80 character wide terminal windows
// and replaced the checkCudaErrors with a plain assert().

// CUDA sample demonstrating a Double precision GEMM computation using the
// Warp Matrix Multiply and Accumulate API introduced in CUDA 11.0.

// In this program, the compute_dgemm kernel computes the result of 
// a matrix multiplication and addition: D = alpha * A * B + beta * C.
// The dimensions of both C and D matrices are M_GLOBAL x N_GLOBAL.
// The A matrix is M_GLOBAL x K_GLOBAL (row-major), the B matrix
// is K_GLOBAL x N_GLOBAL (column-major).
// In that kernel, each CTA computes one 64 x 64 tile of the resulting matrix
// per iteration. When the tile is computed, the CTA stores it to the global
// memory and begins a new iteration, selecting a new 64 x 64 tile to compute.
// Each CTA consists of eight warps. For the 64 x 64 tile,
// each warp computes eight 8 x 8 subtiles, organized in a 2 x 4 
// two-dimensional array.
// Warps compute the 8 x 8 subtiles using nvcuda::wmma::mma_sync operations
// by moving through the K_GLOBAL dimension of the A and B matrices and
// accumulating the intermediate result in the local thread state.

#include <assert.h>
#include <stdio.h>
#include <cuda.h>
#include <mma.h>

// GPU configuration.

#define WARP_SIZE 32

// MMA matrix tile dimensions.

#define M 8
#define N 8
#define K 4

// GEMM configuration.

#define M_TILES 1024
#define N_TILES 1024
#define K_TILES 1024

#define M_GLOBAL (M * M_TILES)
#define N_GLOBAL (N * N_TILES)
#define K_GLOBAL (K * K_TILES)

#define C_LAYOUT wmma::mem_row_major

// Implementation constants.

#define WARPS_PER_BLOCK 8
#define THREADS_PER_BLOCK (WARP_SIZE * WARPS_PER_BLOCK)

#define BLOCK_ROW_WARPS 2
#define BLOCK_COL_WARPS 4

#define WARP_ROW_TILES 4
#define WARP_COL_TILES 2

#define BLOCK_ROW_TILES (WARP_ROW_TILES * BLOCK_ROW_WARPS)
#define BLOCK_COL_TILES (WARP_COL_TILES * BLOCK_COL_WARPS)

// The only kernel remaining does DMMA non-shmem using simple kernel.

using namespace nvcuda;

__host__ void init_host_matrices
  ( double *a, double *b, double *c, int nbrange )
/*
 * Defines the matrices a, b, and c, as follows:
 *   a is upper triangular with 1 on diagonal and above the diagonal,
 *     row major with M_GLOBAL rows and K_GLOBAL columns,
 *   b is upper triangular with increasing sequences of nbrange numbers,
 *     each time doubling the previous number,
 *     column major with K_GLOBAL rows and N_GLOBAL columns,
 *   c is the zero matrix, M_GLOBAL rows, N_GLOBAL columns. */
{
   for(int i=0; i<K_GLOBAL; i++) // #rows: M_GLOBAL > K_GLOBAL
   {
      for(int j=0; j<i; j++) a[i*K_GLOBAL+j] = 0.0;
      for(int j=i; j<K_GLOBAL; j++) a[i*K_GLOBAL+j] = 1.0;
   }
   for(int i=K_GLOBAL; i<M_GLOBAL; i++) // row major with M_GLOBAL rows
      for(int j=0; j<K_GLOBAL; j++) a[i*K_GLOBAL+j] = 0.0;

   for(int i=0; i<K_GLOBAL; i++) // #columns: N_GLOBAL > K_GLOBAL
   {
      for(int j=0; j<=i; j++) b[i*K_GLOBAL+j] = 1.0;
      for(int j=i+1; j<K_GLOBAL; j++) b[i*K_GLOBAL+j] = 0.0;
   }
   for(int i=K_GLOBAL; i<N_GLOBAL; i++) // #columns: N_GLOBAL > K_GLOBAL
   {
      for(int j=0; j<K_GLOBAL; j++) 
      {
         if(j % nbrange == 0)
            b[i*K_GLOBAL+j] = 1.0;
         else
            b[i*K_GLOBAL+j] = 2.0*b[i*K_GLOBAL+j-1];
      }
   }

   for(int t=0; t<M_GLOBAL*N_GLOBAL; t++) c[t] = 0.0;
}

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

__host__ void matMultiplyOnHost
 ( double *A, double *B, double *C, float alpha, float beta,
   int numARows, int numAColumns, int numBRows, int numBColumns,
   int numCRows, int numCColumns )
/*
 * Multiplies the matrices A and B into C on the host. */
{
   for(int i = 0; i < numCRows; i++)
   {
      for(int j = 0; j < numCColumns; j++)
      {
         double temp = 0.0;

         for(int k = 0; k < numAColumns; k++)
         {
            // B matrix is column major. A matrix is row major.
            temp += A[i*numAColumns+k] * B[j*numBRows+k];
         }
         C[i*numCColumns+j] = temp*alpha + beta*C[i*numCColumns+j];
        }
    }
}

int main ( int argc, char **argv )
{
   int compare_host = 0;
   int number_range = 40;
   if(argc > 1) 
   {
      compare_host = 1;
      number_range = atoi(argv[1]);
      printf("number range : %d\n", number_range);
   }

   printf("Initializing ...\n");

   int dev = 0; // findCudaDevice(argc, (const char **)argv);

   cudaDeviceProp deviceProp;
   // replaced checkCudaErrors by assert()
   assert((cudaGetDeviceProperties(&deviceProp, dev)) == 0);

   // Double precision Tensor cores require 
   // a GPU of Ampere (SM8X) architecture or higher.
   if(deviceProp.major < 8)
   {
      printf("dmmaTensorCoreGemm requires SM 8.0 or higher.  Exiting ...\n");
      return -1;
   }
   else
      printf("... running on a SM 8.0 or higher ...\n");

   printf("M: %d (%d x %d)\n", M_GLOBAL, M, M_TILES);
   printf("N: %d (%d x %d)\n", N_GLOBAL, N, N_TILES);
   printf("K: %d (%d x %d)\n", K_GLOBAL, K, K_TILES);

   double *A_h = NULL;
   double *B_h = NULL;
   double *C_h = NULL;

   double *result_hD = NULL;
   double *result_host = NULL;


   A_h = (double*) malloc(sizeof(double) * M_GLOBAL * K_GLOBAL);
   B_h = (double*) malloc(sizeof(double) * K_GLOBAL * N_GLOBAL);
   C_h = (double*) malloc(sizeof(double) * M_GLOBAL * N_GLOBAL);

   if(compare_host > 0)
   {
      result_hD   = (double*) malloc(sizeof(double) * M_GLOBAL * N_GLOBAL);
      result_host = (double*) malloc(sizeof(double) * M_GLOBAL * N_GLOBAL);
   }
   double *A = NULL;
   double *B = NULL;
   double *C = NULL;
   double *D = NULL;
   // replaced checkCudaErrors by assert()
   assert((cudaMalloc((void**)&A,sizeof(double)*M_GLOBAL*K_GLOBAL)) == 0);
   assert((cudaMalloc((void**)&B,sizeof(double)*N_GLOBAL*K_GLOBAL)) == 0);
   assert((cudaMalloc((void**)&C,sizeof(double)*M_GLOBAL*N_GLOBAL)) == 0);
   assert((cudaMalloc((void**)&D,sizeof(double)*M_GLOBAL*N_GLOBAL)) == 0);

   assert(((unsigned long long)A) % 128 == 0);
   assert(((unsigned long long)B) % 128 == 0);
   assert(((unsigned long long)C) % 128 == 0);
   assert(((unsigned long long)D) % 128 == 0);

   printf("Initializing matrices on host ...\n");
   init_host_matrices(A_h, B_h, C_h, number_range);

   printf("Preparing data for GPU ...\n");

   assert((cudaMemcpy(A, A_h, sizeof(double) * M_GLOBAL * K_GLOBAL,
                      cudaMemcpyHostToDevice)) == 0);
   assert((cudaMemcpy(B, B_h, sizeof(double) * N_GLOBAL * K_GLOBAL,
                      cudaMemcpyHostToDevice)) == 0);
   assert((cudaMemcpy(C, C_h, sizeof(double) * M_GLOBAL * N_GLOBAL,
                      cudaMemcpyHostToDevice)) == 0);
   assert((cudaMemset(D, 0, sizeof(double) * M_GLOBAL * N_GLOBAL)) == 0);

   const double alpha = 1.0f; // 1.1f;
   const double beta = 1.0f; // 1.2f;

   cudaEvent_t start, stop;

   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   cudaEventRecord(start);

   dim3 gridDim;
   dim3 blockDim;

   // blockDim.x must be a multple of warpSize 128x4 means we have
   // 16 warps and a block computes a 64x64 output tile
   blockDim.x = 128;
   blockDim.y = 4;

   gridDim.x = (M_GLOBAL + (M*blockDim.x/32-1))/(M*blockDim.x/32);
   gridDim.y = (N_GLOBAL + N*blockDim.y-1)/(N*blockDim.y);

   printf("Computing ... using simple_wmma_gemm kernel ...\n");
   simple_wmma_gemm<<<gridDim, blockDim>>>
      (A, B, C, D, M_GLOBAL, N_GLOBAL, K_GLOBAL, alpha, beta);

   if(compare_host > 0)
   {
      assert(cudaMemcpy(result_hD,D,sizeof(double)*M_GLOBAL*N_GLOBAL,
                        cudaMemcpyDeviceToHost) == 0);
   }
   assert((cudaEventRecord(stop)) == 0);
   assert((cudaEventSynchronize(stop)) == 0);

   if(compare_host > 0)
   {
      printf("Verifying correctness of the computations ...\n");

      memcpy(result_host, C_h, sizeof(double) * M_GLOBAL * N_GLOBAL);

      printf("Computing matrix matrix multiplication on host ...\n");
      matMultiplyOnHost(A_h,B_h,result_host,alpha,beta,M_GLOBAL,K_GLOBAL,
                        K_GLOBAL,N_GLOBAL,M_GLOBAL,N_GLOBAL);

      printf("  host : %.15le\n", result_host[N_GLOBAL-1]);
      printf("device : %.15le\n", result_hD[N_GLOBAL-1]);

      size_t number_of_matches = 0;
      for(int i=0; i<N_GLOBAL*M_GLOBAL; i++)
      {
         if(fabs(result_hD[i] - result_host[i]) > 0.1f)
         {
            printf("mismatch i=%d result_hD=%f result_host=%f\n", i,
                   result_hD[i], result_host[i]);
            break;
         }
         else
         {
            number_of_matches++;
         }
      }
      printf("number_of_matches = %zu out of = %d \n",
             number_of_matches, N_GLOBAL*M_GLOBAL);
      free(result_hD);
      free(result_host);
   }
   float milliseconds = 0;

   cudaEventElapsedTime(&milliseconds, start, stop);

   printf("Time: %f ms\n", milliseconds);
   printf("FP64 TFLOPS: %.2f\n",
          (((double)M_GLOBAL*N_GLOBAL*K_GLOBAL*2)/(milliseconds/1000.))/1e12);

   free(A_h);
   free(B_h);
   free(C_h);
   assert((cudaFree((void*)A)) == 0);
   assert((cudaFree((void*)B)) == 0);
   assert((cudaFree((void*)C)) == 0);
   assert((cudaFree((void*)D)) == 0);

   return 0;
}
