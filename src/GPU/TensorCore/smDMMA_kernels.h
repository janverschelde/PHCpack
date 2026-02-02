/* Code to initialize and multiply the test matrices on the device. */

#ifndef __SMDMMA_KERNELS_H__
#define __SMDMMA_KERNELS_H__

__global__ void simple_wmma_gemm
 ( double *a, double *b, double *c, double *d, int m_ld, int n_ld, int k_ld,
   double alpha, double beta );
/*
   Performs an MxNxK DGEMM (C=alpha*A*B + beta*C) assuming:
   1) Matrices are packed in memory.
   2) M, N and K are multiples of 8, 8 and 4 respectively. 
   3) A is row major, B is column major matrix.
   Note: This is a less performant version of the compute_dgemm kernel.
   It is designed for demonstration purposes only to show the CUDA WMMA API 
   use without relying on availability of the shared memory.  */

__host__ void GPU_DMMA
 ( double *A, double *B, double *C, double *D, double alpha, double beta );
/*
 * Launches the kernels to multiply the matrices in A and B. */

#endif
