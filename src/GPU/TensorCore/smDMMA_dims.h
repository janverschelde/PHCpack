/* Dimension for the matrix multiplication with a simple tensor core kernel. */

#ifndef __SMDMMA_DIMS_H__
#define __SMDMMA_DIMS_H__

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

// use of shared memory

#if SHARED_MEMORY_LIMIT_64K
// With only 64 Kb shared memory available, 
// we can fit 8x16-tile chunks of each the A and B matrix data,
// that are (M = 8) * (K = 4) * 8 * (CHUNK_K = 16) * sizeof(double)
// = 32 Kb each.  But we cannot account the 4 Kb total skew overhead,
// without which the performance would be severely impacted.
// So we choose to reduce the chunk size in half,
// i.e. the amount of A and B matrix data we cache in shared memory.
// Accordingly, this doubles the number of outer iterations across 
// the global K dimension, which only slightly impacts the performance.
#define CHUNK_K 8
#else
#define CHUNK_K 16
#endif

#define CHUNK_LINE_BYTES (CHUNK_K * K * sizeof(double))
#define WARP_COPY_BYTES (WARP_SIZE * sizeof(int4))
#define CHUNK_COPY_LINES_PER_WARP (WARP_COPY_BYTES / CHUNK_LINE_BYTES)
#define CHUNK_COPY_LINE_LANES (WARP_SIZE / CHUNK_COPY_LINES_PER_WARP)

#define GLOBAL_MEM_STRIDE N_GLOBAL

#define SHMEM_STRIDE (N * BLOCK_ROW_TILES)
#define SHMEM_OFFSET (N * WARP_ROW_TILES)

#define SKEW_DOUBLE 4

#endif
